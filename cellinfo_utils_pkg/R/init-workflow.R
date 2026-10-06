#' Initialize a cellinfo object from a Seurat object
#'
#' Creates a new cellinfo list with sampled cells and marker lists for a clustering variable, hard-refreshes annotations, and optionally seriates cells.
#'
#' @note Requires companion helpers such as `setmarkers()`, `cellinfo.getcells()`, `cellinfo.get.markerlist()`, and `harmonise.all.markers()` from the parent mlutils codebase.
#'
#' @param so A Seurat object.
#' @param seriate Logical; whether to seriate cells with module scores.
#' @param ncells Number of cells to sample per group.
#' @param clusvar Metadata column defining cell groups.
#' @param markerset Marker set key / label used for marker grouping.
#' @param ... Additional arguments reserved for compatibility.
#'
#' @return A new `cellinfo` list.
#'
#' @seealso [refresh.cellinfo.hard()], [cellinfo.seriatecells.seurat()]
#'
#' @export
cellinfo.init <- function(so, seriate=T, ncells=100,clusvar="seurat_clusters",markerset=clusvar, ...){
 cellinfo=list()
 so=setmarkers(so, markerset=markerset)
   cellinfo$cell.list=cellinfo.getcells(so, clusvar=clusvar, ncells=ncells)
   cellinfo$markerlist=cellinfo.get.markerlist(harmonise.all.markers(so), clusvar=clusvar)
   #if(is.numeric.like(names(cellinfo$markerlist))){
  # names(cellinfo$markerlist)=paste_("cluster", names(cellinfo$markerlist))
    # }
   #if(is.numeric.like(names(cellinfo$cell.list))){
    # names(cellinfo$cell.list)=paste_("cluster", names(cellinfo$cell.list))
   #}
   fcat("pt1")
   cellinfo<-refresh.cellinfo.hard(cellinfo, cell.group.label=clusvar, marker.group.label=markerset)
   fcat("cell group variable is", cellinfo$cell.group.variable)
   if(seriate){
     fcat("pt2")
cellinfo=cellinfo.seriatecells.seurat(cellinfo, so, new.cells=F)  
  
   }
   
cellinfo   
}

#' End-to-end cellinfo heatmap workflow for one dataset
#'
#' Seriates cells for a grouping variable and builds a labelled heatmap, optionally returning the Seurat object as well.
#'
#' @note Requires companion helper `seriatecells()` from mlutils.
#'
#' @param so A Seurat object.
#' @param grouping.variable Metadata column defining groups.
#' @param assay Assay used for visualisation.
#' @param markers Optional marker table; otherwise uses markers stored on the object.
#' @param recalculate.markers Logical; reserved for compatibility.
#' @param marker.params Optional marker-calculation parameters.
#' @param colorlist Annotation colors.
#' @param ncells Cells per group for seriation.
#' @param return.seurat If `TRUE`, include the Seurat object in the returned list.
#'
#' @return List from [cellinfo.heatmap()] with `return.everything=TRUE`, optionally including `seuratobject`.
#'
#' @export
cellinfo.workflow.singleds <- function(so, grouping.variable="seurat_clusters", assay="RNA", markers=NULL, recalculate.markers=F, marker.params=NULL, colorlist=NULL, ncells=300, return.seurat=F){
  
  if(assay=="RNA"){
    fcat("Warning: assay is set as RNA. consider using a SCT assay for better visualisation") 
  }
  # housekeeping.
  #check that there is an overlap between categories of clusvar and categories of group1 in markers table  
  ###############################################################################
  #Step A. locate or recalculate markers. 
  ###############################################################################
  
  if(is.null(so$RNA@misc$top_markers) & is.null(markers)){
    fcat("no markers found. Please calculate markers and pass them to argunment markers or to so$RNA@misc$top_markers") 
    
  }
  ##############
  # Step b.  generate cellinfo with marker signatures. 
  #############
  if(is.null(colorlist)){
    clusterss=NULL  
  }else{
    clusterss=names(colorlist[[grouping.variable]])
    
  }
  
  cellinfo= seriatecells(so, clusvar=grouping.variable, clusters= clusterss, meth="seurat", extended.output=T, deduped=T, ncells=ncells)
  
  outs= cellinfo.heatmap(so, cellinfo,return.everything=T, genes.to.label=adjust.markers(cellinfo, 2)$markers, colorlist=colorlist)
  if(return.seurat){
    outs[["seuratobject"]]=so 
  }
  
  outs
  
}
