#' Plot a variable along within-group cell ranks
#'
#' Creates a ggplot path plot of a stored `cellinfo$values` variable against each cell relative rank within its group.
#'
#' @param cellinfo A cellinfo list object with populated `values`.
#' @param variable Name of the value series under `cellinfo$values`.
#' @param scale.x Logical; currently used when building ranks.
#' @param clusters Group names to include.
#' @param colorlist Named color list keyed by `groupvar`.
#' @param groupvar Grouping variable name for colors.
#'
#' @return A ggplot object.
#'
#' @export
cellinfo.ordered.cell.plot <- function(cellinfo, variable, scale.x=T, clusters=names(cellinfo$cell.list), colorlist=NULL, groupvar=cellinfo$clustering.variable){
  
  if(is.null(colorlist)){
    
    colorlist=list()
    
    catnames=names(cellinfo$cell.list)
    colorlist[[groupvar]]= randomcolors(length(catnames)) %>% givename(., catnames)
  }
  
  # structure data in df for ggplot
  
  
  plotdf=lapply(clusters, function(gr){  
    numcells= length(cellinfo$cell.list[[gr]])
    
    if(scale.x){
      tot=numcells 
    }else{
      tot=1 
    }
    
    
    cellranks=(1:numcells)/numcells
    vals=cellinfo$values[[variable]][[gr]]
    
    mt1=cellranks %>% as.data.frame %>% givecolnames(., nms="relative.rank.in.group")
    mt1[, variable]=vals
    mt2=giverownames(mt1, cellinfo$cell.list[[gr]]) %>% dplyr::mutate(group=!!gr)
    
    mt2    
  }) %>% Reduce(rbind, .)
  
  ggplot(plotdf)+geom_path(aes(x=relative.rank.in.group, y=!!sym(variable), color=group))+theme_classic()+xlab("Cell's rank in group")+NoLegend()+scale_color_manual(values=colorlist[[groupvar]])#+ylab(!!sym(variable))
  
}

#' FeaturePlot selected markers from a cellinfo object
#'
#' Runs Seurat `FeaturePlot` on markers retained after [adjust.markers()].
#'
#' @note This legacy function currently references `so.mes` inside the body; callers may need to align object names.
#'
#' @param so A Seurat object (note: legacy body references `so.mes` and may need local adaptation).
#' @param cellinfo A cellinfo list object.
#' @param markers.per.group Markers retained per group for plotting.
#' @param umap.number Optional index into UMAP reductions; default uses the last UMAP.
#'
#' @return A FeaturePlot / patchwork-like plot object.
#'
#' @export
cellinfo.plotfeatures <- function(so, cellinfo, markers.per.group=3,umap.number=NULL){
  

   nms=names(so.mes@reductions)
   unms=nms[grepl("umap", nms)]
   
  if(is.null(umap.number)){
   rr=unms[length(unms)]# using the last umap generated
  }else{
   rr=unms[umap.number] 
  }
  
plt=FeaturePlot(so.mes, features=adjust.markers(cellinfo, markers.per.group)$markers, reduction=rr)  

plt}

#' Summarise marker expression for a bubble/point-style table
#'
#' Computes mean expression and percent-expressed summaries for cellinfo markers across groups. Returns a long data frame suitable for ggplot bubble plots.
#'
#' @note Requires companion helpers `cellinfo.getcells()` and `join_meta_exp2()` from mlutils.
#'
#' @param so A Seurat object.
#' @param cellinfo A cellinfo list object.
#' @param ncells Cells sampled per cluster when building the expression matrix.
#' @param groupby Grouping metadata column.
#' @param assay Assay used for expression extraction.
#'
#' @return A summarised tibble/data frame with mean expression and percent expressed.
#'
#' @export
cellinfo.pointplot <- function(so, cellinfo, ncells=1000000,groupby=cellinfo$cell.group.variable,  assay=DefaultAssay(so)){
  library(viridis)
  library(tidyr)
  
  fcat("Warning: pointplot by default is made with all the cells in the seurat object.\n
       if you would like to use the cellinfo cells instead set allcells=F") 
 
   cell.list=cellinfo.getcells(so,"seurat_clusters",  ncells = ncells)
   markers=cellinfo$markers
  so=GetResidual(so, markers)
   mat=join_meta_exp2(so, genes=markers, cells=Reduce(c, cell.list), assay=assay) 
   #fcat("the following genes have not been included:")
   #fcat(Reduce(pastec, setdiff(make.names(markers), colnames(mat))))
   
   mat2=mat %>% pivot_longer(., setdiff(colnames(mat), colnames(so@meta.data)), names_to="gene", values_to="expression") %>% group_by(!!sym(groupby), gene) %>% summarise(mn=mean(expression), pct_expressed=sum(expression>0)/n())
   
   #lst=list(plot=ggplot(mat2, aes(y=factor(gene, levels=intersect( rev(make.names(markers)), mat2$gene)), x=!!sym(groupby), color=mn, size=pct_expressed ))+geom_point()+scale_color_viridis(), 
            
          #  mat=mat2)
   #lst
   mat2
}

#' Scatter module score vs marker percentage by group
#'
#' Adds module scores and marker percentages for each marker set, then builds a patchwork of scatter plots coloured by cell group.
#'
#' @note Requires `AddModuleScore3()`, `AddMarkerPercentages()`, and optionally `install_and_load()` / patchwork from the parent project.
#'
#' @param so A Seurat object.
#' @param cellinfo A cellinfo list object.
#' @param colorlist Named annotation colors (defaults to global `allcolors`).
#'
#' @return A combined ggplot/patchwork object.
#'
#' @name plot.signature.profiles
#' @usage plot.signature.profiles(so, cellinfo, colorlist = allcolors)
#' @export plot.signature.profiles
plot.signature.profiles <- function(so, cellinfo, colorlist=allcolors){
install_and_load("patchwork")
fcat("calculating signatures")
so=AddModuleScore3(so, cellinfo$markerlist)  
fcat("calculating perentages")
so=AddMarkerPercentages(so, cellinfo$markerlist)  

if(is.null(colorlist[[cellinfo$cell.group.variable]])){
  fcat("creating colors for categories")
  colorlist[[cellinfo$cell.group.variable]]<-rainbow(length(cellinfo$cell.list)) %>% givename(., names(cellinfo$cell.list))
allcolors[[cellinfo$cell.group.variable]]<<-rainbow(length(cellinfo$cell.list)) %>% givename(., names(cellinfo$cell.list))
}

plt=lapply(1:length(cellinfo$markerlist), function(x){
  
ggplot(so@meta.data, aes(x=!!sym(paste0("mscore_group_",names(cellinfo$markerlist)[x])), y=!!sym(paste0("pct_group_",names(cellinfo$markerlist)[x])), color=!!sym(cellinfo$cell.group.variable)))+geom_point()+scale_color_manual(values=allcolors[[cellinfo$cell.group.variable]])

  }) %>% Reduce('+', .)
plt
}

#' Build sparse column labels for cell-group midpoints
#'
#' Returns a character vector aligned to concatenated cells, with group names placed at the middle index of each group and `NA` elsewhere. Useful as heatmap column labels.
#'
#' @param cellinfo A cellinfo list object.
#'
#' @return Character vector of labels.
#'
#' @export
prepare.cellgroup.labels <- function(cellinfo){
 
  label.list=lapply(names(cellinfo$cell.list), function(name) {
    vec_length <- length(cellinfo$cell.list[[name]])
    middle <- ceiling(vec_length / 2)
    out <- rep(NA_character_, vec_length)
    out[middle] <- name
    out
  })
  
  label.list %>% Reduce(c, .)
  
  }
