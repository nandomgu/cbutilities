#' Refresh derived cellinfo fields while preserving annotations
#'
#' Rebuilds flattened `cells`/`markers`, gaps, and annotations from `cell.list`/`markerlist`, carrying forward previous annotation columns where possible.
#'
#' @param cellinfo A cellinfo list object.
#' @param cell.group.label Fallback label for the primary cell-group annotation column.
#' @param marker.group.label Fallback label for the primary marker-group annotation column.
#'
#' @return Refreshed `cellinfo` list.
#'
#' @seealso [refresh.cellinfo.hard()]
#'
#' @export
refresh.cellinfo <- function(cellinfo, cell.group.label="cell.group", marker.group.label="marker.group"){
  fcat("readjusting")
  
  cellinfo$markers= cellinfo$markerlist %>% Reduce(c, .)
  cellinfo$cells= cellinfo$cell.list %>% Reduce(c, .)
  
  if(is.null(names(cellinfo$cell.list))){
   names(cellinfo$cell.list) = 1:length(cellinfo$cell.list)
  }
  
  fcat("making references")
  cell.reference=lapply(1:length(cellinfo$cell.list), function(nn) rep(names(cellinfo$cell.list)[nn], length(cellinfo$cell.list[[nn]])   )) %>% Reduce(c, .)  
  marker.reference=lapply(1:length(cellinfo$markerlist), function(nn) rep(names(cellinfo$markerlist)[nn], length(cellinfo$markerlist[[nn]])   )) %>% Reduce(c, .)
  
  ##making sure there arent deduped genes
    fcat("removing potential duplicates in markers")
  tff=!duplicated(cellinfo$markers)
  cellinfo$markers=cellinfo$markers[tff]
  marker.reference=marker.reference[tff]
  
  fcat("removing potential duplicates in cells")
  tfc=!duplicated(cellinfo$cells)
  cellinfo$cells=cellinfo$cells[tfc]
  cell.reference=cell.reference[tfc]
  
  ## reassemble deduplicated lists
   fcat("reassembling deduped lists")
  markernames=cellinfo$markerlist %>% names
  cellnames=cellinfo$cell.list %>% names
  
  cellinfo$markerlist=lapply(names(cellinfo$markerlist), function(x) cellinfo$markers[marker.reference==x]) %>% givename(., markernames)
  cellinfo$cell.list=lapply(names(cellinfo$cell.list), function(x) cellinfo$cells[cell.reference==x]) %>% givename(., cellnames)
  
    fcat("designing cell and marker gaps for heatmaps")
  cellinfo$gaps.markers=lapply(1:length(cellinfo$markerlist), function(x) rep(x, length(cellinfo$markerlist[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which
  cellinfo$gaps.cells= lapply(1:length(cellinfo$cell.list), function(x) rep(x, length(cellinfo$cell.list[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which
  
  fcat("making annotations")
  if(!is.null(cellinfo$clustering.variable)){
    cell.group.label=cellinfo$clustering.variable
    
    if(zsoverlap(names(cellinfo$markerlist), names(cellinfo$cell.list))==1 ){
      fcat("Labels of cells and markers are shared.\n Assuming", cellinfo$clustering.variable,"as common variable for markers and cells")
      marker.group.label=cellinfo$clustering.variable
    }
  }else{
    fcat("No clustering variable. Checking metadata labels")
    if(!is.null(cellinfo$cell.metadata)){
      fcat("Cell metadata found")
      cell.group.label=cellinfo$cell.metadata[1]
    }
    
    if(!is.null(cellinfo$marker.metadata)){
      fcat("Marker metadata found")
      marker.group.label=cellinfo$marker.metadata[1]
    }
  }
  
  #if(is.null(cellinfo$cell.annotation)){}
  previous.cell.annotation=cellinfo$cell.annotation
  fresh.cell.annotation=as.data.frame(cell.reference) %>% giverownames(., cellinfo$cells) %>% givecolnames(., 1, cell.group.label)
  
  
  
  if(is.null(previous.cell.annotation) || length(intersect(rownames(previous.cell.annotation), cellinfo$cells))==0){
    cellinfo$cell.annotation= fresh.cell.annotation
  }else{
    #bring in the previous cell information dat for the new cells.
    othercols=setdiff(colnames(previous.cell.annotation), cell.group.label) %>% unique
    
    added.cell.annotation=previous.cell.annotation %>% getrows(., rownames(fresh.cell.annotation)) %>% select(all_of(othercols))
    
    cellinfo$cell.annotation=cbind(fresh.cell.annotation, added.cell.annotation)
    
  }
  
  ## annotation of markers
  previous.marker.annotation=cellinfo$marker.annotation
  fresh.marker.annotation= as.data.frame(marker.reference) %>% giverownames(., cellinfo$markers) %>% givecolnames(., 1, marker.group.label)
  
  if(is.null(previous.marker.annotation)|| length(intersect(rownames(previous.marker.annotation), cellinfo$markers))==0){
    cellinfo$marker.annotation= fresh.marker.annotation
  }else{
    #bring in the previous marker information dat for the new cells, except of the first column
    othercols.markers=setdiff(colnames(previous.marker.annotation), marker.group.label) %>% unique
    
    added.marker.annotation=previous.marker.annotation %>% getrows(., rownames(fresh.marker.annotation)) %>% select(all_of(othercols.markers))
    cellinfo$marker.annotation=cbind(fresh.marker.annotation, added.marker.annotation)
    
  }
  
  
  cellinfo$cell.metadata= colnames(cellinfo$cell.annotation)
  cellinfo$marker.metadata= colnames(cellinfo$marker.annotation)
  
  cellinfo
}

#' Hard-refresh cellinfo annotations from group labels
#'
#' Like [refresh.cellinfo()], but rebuilds cell/marker annotation tables from scratch using the provided group labels and sets `cell.group.variable` / `marker.group.variable`.
#'
#' @param cellinfo A cellinfo list object.
#' @param cell.group.label Primary cell-group annotation column name.
#' @param marker.group.label Primary marker-group annotation column name.
#'
#' @return Hard-refreshed `cellinfo` list.
#'
#' @seealso [refresh.cellinfo()], [cellinfo.init()]
#'
#' @export
refresh.cellinfo.hard <- function(cellinfo, cell.group.label="cell.group", marker.group.label="marker.group"){
  

    fcat("readjusting")
  
  cellinfo$markers= cellinfo$markerlist %>% Reduce(c, .)
  cellinfo$cells= cellinfo$cell.list %>% Reduce(c, .)
  
  if(is.null(names(cellinfo$cell.list))){
   names(cellinfo$cell.list) = 1:length(cellinfo$cell.list)
  }
  
  fcat("making references")
  cell.reference=lapply(1:length(cellinfo$cell.list), function(nn) rep(names(cellinfo$cell.list)[nn], length(cellinfo$cell.list[[nn]])   )) %>% Reduce(c, .)  
  marker.reference=lapply(1:length(cellinfo$markerlist), function(nn) rep(names(cellinfo$markerlist)[nn], length(cellinfo$markerlist[[nn]])   )) %>% Reduce(c, .)
  
  ##making sure there arent deduped genes
    fcat("removing potential duplicates in markers")
  tff=!duplicated(cellinfo$markers)
  cellinfo$markers=cellinfo$markers[tff]
  marker.reference=marker.reference[tff]
  
  fcat("removing potential duplicates in cells")
  tfc=!duplicated(cellinfo$cells)
  cellinfo$cells=cellinfo$cells[tfc]
  cell.reference=cell.reference[tfc]
  
  ## reassemble deduplicated lists
   fcat("reassembling deduped lists")
  markernames=cellinfo$markerlist %>% names
  cellnames=cellinfo$cell.list %>% names
  
  cellinfo$markerlist=lapply(names(cellinfo$markerlist), function(x) cellinfo$markers[marker.reference==x]) %>% givename(., markernames)
  cellinfo$cell.list=lapply(names(cellinfo$cell.list), function(x) cellinfo$cells[cell.reference==x]) %>% givename(., cellnames)
  
    fcat("designing cell and marker gaps for heatmaps")
  cellinfo$gaps.markers=lapply(1:length(cellinfo$markerlist), function(x) rep(x, length(cellinfo$markerlist[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which
  cellinfo$gaps.cells= lapply(1:length(cellinfo$cell.list), function(x) rep(x, length(cellinfo$cell.list[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which
  
  
   #if(is.null(cellinfo$cell.annotation)){}

  cellinfo$cell.annotation=as.data.frame(cell.reference) %>% giverownames(., cellinfo$cells) %>% givecolnames(., 1, cell.group.label)
  
  cellinfo$marker.annotation= as.data.frame(marker.reference) %>% giverownames(., cellinfo$markers) %>% givecolnames(., 1, marker.group.label)
  
  cellinfo$cell.metadata= colnames(cellinfo$cell.annotation)
  cellinfo$marker.metadata= colnames(cellinfo$marker.annotation)
  cellinfo$cell.group.variable=ifelse(cell.group.label=="cell.group", NULL, cell.group.label)
  cellinfo$marker.group.variable=ifelse(cell.group.label=="marker.group", NULL, marker.group.label)
  cellinfo
  
}
