#' Keep the last n cells in each cellinfo group
#'
#' Trims each vector in `cellinfo$cell.list` to at most `n` cells (keeping the end of each vector), then refreshes derived fields.
#'
#' @param cellinfo A cellinfo list object.
#' @param n Integer; maximum number of cells to retain per group.
#'
#' @return Updated `cellinfo` list.
#'
#' @seealso [adjust.markers()], [refresh.cellinfo()]
#'
#' @export
adjust.cells <- function(cellinfo, n){
  cellinfo$cell.list=lapply(cellinfo$cell.list, function(x){
    
    if(n<=length(x)){
      removenas(x[(length(x)-n):length(x)])
    }else{
      x
    }
    
    
    
  }) %>% givename(., cellinfo$cell.list %>% names) 
  cellinfo=refresh.cellinfo(cellinfo)
  cellinfo
}

#' Keep the first n markers in each marker group
#'
#' Trims each vector in `cellinfo$markerlist` to at most `n` markers and refreshes derived fields.
#'
#' @param cellinfo A cellinfo list object.
#' @param n Integer; maximum markers per group.
#'
#' @return Updated `cellinfo` list.
#'
#' @seealso [adjust.cells()], [removemarkers()]
#'
#' @export
adjust.markers <- function(cellinfo, n){
  cellinfo$markerlist= adjust.markerlist(cellinfo$markerlist, n)
  cellinfo=refresh.cellinfo(cellinfo)
  cellinfo
}

#' Remove markers from all marker groups
#'
#' Deletes specified genes from every vector in `cellinfo$markerlist` and refreshes derived fields.
#'
#' @param cellinfo A cellinfo list object.
#' @param markers Character vector of genes to remove.
#'
#' @return Updated `cellinfo`.
#'
#' @export
removemarkers <- function(cellinfo, markers){
  nms=names(cellinfo$markerlist)
  cellinfo$markerlist=lapply(1:length(cellinfo$markerlist), function(x){
    gns=cellinfo$markerlist[[x]]
    
    gns[!(gns %in% markers)]
    
  }) %>% givename(., nms)
  
  refresh.cellinfo(cellinfo)
}

#' Subset cellinfo to selected group labels
#'
#' Keeps only marker/cell groups whose names are in `labels` (intersected with available markerlist names) and refreshes derived fields.
#'
#' @param cellinfo A cellinfo list object.
#' @param labels Character vector of group labels to retain.
#'
#' @return Subsetted and refreshed `cellinfo`.
#'
#' @name subset.cellinfo
#' @usage subset.cellinfo(cellinfo, labels)
#' @export subset.cellinfo
subset.cellinfo <- function(cellinfo, labels){
  fcat("restricting to labels only present in list")
  tf=lapply(labels, function(x) x %in% (cellinfo$markerlist %>% names)) %>% Reduce(c, .)
  labels2=labels[tf]
  
  cellinfo$markerlist=cellinfo$markerlist[labels2]
  cellinfo$cell.list=cellinfo$cell.list[labels2]
  cellinfo$input.cluster.ids=labels2
  refresh.cellinfo(cellinfo)
}

#' Resegregate cells by one or more annotation variables
#'
#' Reorders and splits cells according to `variables`, rebuilding `cell.list` and a reduced `cell.annotation`, then refreshes the object.
#'
#' @param cellinfo A cellinfo list object.
#' @param variables Character vector of annotation columns; first column defines split groups.
#'
#' @return Resegregated and refreshed `cellinfo`.
#'
#' @export
segregate.cells <- function(cellinfo, variables){
rearranged.cells.df.list=cellinfo$cell.annotation  %>% names2col(., "cellid") %>% arrange(!!!syms(variables) ) %>% group_split(., !!sym(variables[1]))

rear.names=lapply(rearranged.cells.df.list, function(x) x %>% pull(!!sym(variables[1])) %>% unique) %>% Reduce(c, .)

rearranged.cell.annot=cellinfo$cell.annotation  %>% names2col(., "cellid") %>% arrange(!!!syms(variables)) %>% group_split(., !!sym(variables[1])) %>% bind_rows %>% as.data.frame %>%  col2names(., "cellid") 


rearranged.cell.list=rearranged.cells.df.list %>% lapply(., function(x) x %>% as.data.frame %>% pull(cellid)) %>% givenames(., rear.names)
rearranged.cells=rearranged.cells.df.list %>% lapply(., function(x) x %>% as.data.frame %>% pull(cellid))
cellinfo$cell.list=rearranged.cell.list
cellinfo$cell.annotation=rearranged.cell.annot %>% select(-cellid) %>% dplyr::select(all_of(variables))
cellinfo$cells=rearranged.cells

refresh.cellinfo(cellinfo)
}

#' Keep cellinfo groups whose names match a pattern
#'
#' Subsets `cellinfo$cell.list` to names matching `pattern` and refreshes derived fields.
#'
#' @param cellinfo A cellinfo list object.
#' @param pattern Regular expression passed to `grepl()`.
#'
#' @return Filtered and refreshed `cellinfo`.
#'
#' @export
keep.cellgroups.pattern <- function(cellinfo, pattern){
 cellinfo$cell.list = cellinfo$cell.list[grepl(pattern, names(cellinfo$cell.list))] 
refresh.cellinfo(cellinfo)
 }
