#' Draw a ComplexHeatmap from a cellinfo object
#'
#' Builds a row-scaled expression heatmap for `cellinfo$markers` x `cellinfo$cells`, using cell and marker annotations for column/row annotation tracks.
#'
#' @param so A Seurat object.
#' @param cellinfo A cellinfo list object.
#' @param assay Assay providing the data layer.
#' @param name Optional heatmap title.
#' @param pheatmap.params Currently unused; reserved for compatibility.
#' @param genes.to.label Optional gene names annotated on the right.
#' @param color.labels.by Optional marker-annotation column used to color gene labels.
#' @param colorlist Named list of annotation colors.
#' @param cmap Three-color vector defining the expression palette.
#' @param colorbar.max Absolute max for scaled color breaks.
#' @param colorbar.spacing Spacing between color breaks.
#' @param show.legend Logical; retained for API compatibility.
#' @param return.everything If `TRUE`, return list(heatmap, cellinfo, matrix).
#' @param return.cellinfo If `TRUE` (and not `return.everything`), return updated cellinfo.
#'
#' @return A ComplexHeatmap object, updated cellinfo, or a list depending on return flags.
#'
#' @seealso [heatmap.workflow()], [refresh.cellinfo()]
#'
#' @export
cellinfo.heatmap <- function(so,
                           cellinfo,
                           assay="RNA",
                           name=NULL,
                           pheatmap.params=NULL, ## currently not used
                           genes.to.label=NULL, 
                           color.labels.by=NULL,
                           colorlist=NULL,
                           cmap=NULL,
                           colorbar.max=2,
                           colorbar.spacing=0.1,
                           show.legend=TRUE,
                           return.everything=F,
                           return.cellinfo=F){
  
  #function refresh cell info to update ubject
  
  
  
  # heatmap colors
  
  if(is.null(cmap)){
    cmap=c("#0B9988", "#f7f7f7", "#F9AD03") 
  }
  
  colorpalette=colorRampPalette(cmap, space="rgb")
  
  ####Heatmap color schemes
  #petrol gold
  
  br=seq(-colorbar.max,colorbar.max,colorbar.spacing)  ## numeric breaks of the color bins
  colls=colorpalette(length(br))   
  
  #cellinfo= refresh.cellinfo(cellinfo)
  
  
  
  
  
  fcat("getting data matrix v2...")
  # upgraded for seurat >v4
mat<- GetAssayData(so, assay = assay, layer = "data")[cellinfo$markers, cellinfo$cells] %>% as.matrix
  fcat("finding invariant and absent genes...")
  absent.genes= setdiff(cellinfo$markers, rownames(mat))
  variant.rows=get.variant.rows(mat)
  invariant.genes=mat[!variant.rows, ] %>% rownames 
  fcat("Warning: the following invariant genes have been removed:\n", invariant.genes %>% dput)
  cellinfo2=removemarkers(cellinfo, markers=union(invariant.genes, absent.genes))
  
  if(!is.null(cellinfo$cell.metadata)){
    fcat("adding cell annotation from variables in metadata...")
    cellinfo2=cellinfo2 %>% add.cell.annotation(so, ., cellinfo2$cell.metadata, assay=assay)
    
     
      }
  
  
  
  if(!is.null(genes.to.label)){
    fcat("finding genes to label...")
    # original markers
    gns=cellinfo2$markers
    
    
    geneguide=  1:length(gns)
    names(geneguide)=gns
    #mrkrs are genes of interest
    
    
    
    if(!is.null(color.labels.by) & !is.null(colorlist)){
      
      gene.order.info=geneguide[genes.to.label]
      
      
      getcol=Vectorize(function(x) 
        if(x %in% (allcolors[[color.labels.by]] %>% names)){
          allcolors[[color.labels.by]][x]}else{
            NA
          }, USE.NAMES=F)
      
      getposition=Vectorize(function(x) geneguide[x], USE.NAMES=F)
      
      colref= (cellinfo2$marker.annotation %>% names2col(., "gene") %>% filter(gene %in% genes.to.label) %>% dplyr::mutate(color=getcol(!!sym(color.labels.by)), position= getposition(gene) )  )
      genecols=colref %>% arrange(position) %>% pull(color) %>% unname
      
      
      linesgp= grid::gpar(col=genecols)
      labelsgp= grid::gpar(col=genecols)
    }else{
      linesgp=NULL
      labelsgp=NULL
    }
    ha = ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at = geneguide[genes.to.label] %>% unname %>% as.numeric, 
                                                                       labels = genes.to.label,
                                                                       lines_gp= linesgp, 
                                                                       labels_gp= labelsgp))
    
    
    #geneguide[genes.to.label] %>% as.character %>% dput
    
  }else{ha=NULL}
  
  fcat("getting data matrix v1...")
  
  mat2<<-GetAssayData(so, assay = assay, layer = "data")[cellinfo2$markers, cellinfo2$cells] %>% as.matrix

  
  cats=cellinfo2$cell.metadata
  
  for(ct in cats){
    
    if(!is.null(allcolors[[ct]])){  
      cellinfo2$cell.annotation[[ct]]= factor(cellinfo2$cell.annotation[[ct]], levels=colorlist[[ct]] %>% names) 
    }
  }
  
  #setting factor for cell group variable
  cellinfo2$cell.annotation[[cellinfo2$cell.group.variable]]= factor(cellinfo2$cell.annotation[[cellinfo2$cell.group.variable]], levels=colorlist[[cellinfo2$cell.group.variable]] %>% names) 
  
    cats2=cellinfo2$marker.metadata
  
  for(ct in cats2){
    
    if(!is.null(allcolors[[ct]])){  
      cellinfo2$marker.annotation[[ct]]= factor(cellinfo2$marker.annotation[[ct]], levels=colorlist[[ct]] %>% names) 
    }
  }
  
  
  ph= ComplexHeatmap::pheatmap(mat2,
                               scale="row",
                               cluster_row=FALSE,
                               show_colnames=FALSE,
                               show_rownames=FALSE,
                               cluster_col=FALSE, 
                               breaks=br,
                               col=colls, 
                               annotation_col=cellinfo2$cell.annotation,
                               annotation_row=cellinfo2$marker.annotation,
                               gaps_col= cellinfo2$gaps.cells,
                               gaps_row= cellinfo2$gaps.markers,
                               annotation_colors=colorlist,
                               border_color=NA,
                               right_annotation=ha, 
                               fontsize=5,
                               main=name,
                               labels_col=F,#prepare.cellgroup.labels(cellinfo2)
                               #show_heatmap_legend=show.legend
                               #column_gap = unit(.2, "mm")
                               annotation_names_row = FALSE,
                               annotation_names_col = FALSE
  )
  
  
  if(return.everything){
    
    list(heatmap=ph, cellinfo=cellinfo2, matrix=mat2)
  }else{
    
    if(return.cellinfo){
      return(cellinfo2) 
    }else{
      
      return(ph)
    }
    
  }
  
}

#' Annotate, segregate, and save a cellinfo heatmap
#'
#' Convenience wrapper around annotation, segregation, natural ordering, and [cellinfo.heatmap()], writing a PDF via `tpdf()`.
#'
#' @note Uses global `allcolors`, `landmark.genes`, and helper `tpdf()` from the parent project when available.
#'
#' @param so A Seurat object.
#' @param cellinfo A cellinfo list object.
#' @param vars Variables used to segregate cells.
#' @param colorlist Named annotation colors (defaults to global `allcolors`).
#' @param extra.genes Extra gene labels.
#' @param extra.annotation Extra annotation variables to import.
#' @param assay Assay name.
#' @param nmarkers Markers per group to label.
#' @param return.heatmap If `TRUE`, return heatmap instead of cellinfo.
#' @param label Optional filename label prefix.
#'
#' @return Updated `cellinfo` or heatmap object.
#'
#' @seealso [cellinfo.group.signatures()], [cellinfo.heatmap()]
#'
#' @export
heatmap.workflow <- function(so, cellinfo,vars, colorlist=allcolors, extra.genes=NULL, extra.annotation=NULL, assay="SCT", nmarkers=3, return.heatmap=F, label=""){
cellinfo2=add.cell.annotation(so, cellinfo=cellinfo , vars = vars, assay=assay)
cellinfo3=cellinfo2 %>% segregate.cells(., variables = c(cellinfo$cell.group.variable, vars)) %>% cellinfo.order.groups.natural
if(is.null(colorlist[[cellinfo3$cell.group.variable]])){
  
  colorlist[[cellinfo3$cell.group.variable]]<-rainbow(length(cellinfo3$cell.list)) %>% givename(., names(cellinfo3$cell.list))
allcolors[[cellinfo3$cell.group.variable]]<<-rainbow(length(cellinfo3$cell.list)) %>% givename(., names(cellinfo3$cell.list))
}
fcat(dput(colorlist[[cellinfo3$cell.group.variable]]))

if(!is.null(extra.annotation)){
  cellinfo3=cellinfo3 %>% add.cell.annotation(so, ., extra.annotation, assay=assay, overwrite=T)
}
#allcolors[[cellinfo3$marker.group.variable]]=rainbow(length(cellinfo3$markerlist)) %>% givename(., names(cellinfo3$markerlist))

cellinfo3$cell.annotation[[cellinfo3$cell.group.variable]]=factor(cellinfo3$cell.annotation[[cellinfo3$cell.group.variable]], levels=ifelse(!is.null(colorlist[[cellinfo3$cell.group.variable]]), colorlist[[cellinfo3$cell.group.variable]] %>% names, gtools::mixedsort(unique(cellinfo3$cell.annotation[[cellinfo3$cell.group.variable]]))))


hm=cellinfo.heatmap(so=so, cellinfo=cellinfo3, assay=assay, genes.to.label = c((cellinfo3 %>% adjust.markers(., nmarkers))$markers, landmark.genes, extra.genes) , colorlist=colorlist)

tpdf(paste0("heatmap_",label,ifelse(label=="", "", "_"), Project(so),"_by_", cellinfo$cell.group.variable),  sca=6)
print(hm)
dev.off()

if(!return.heatmap){
cellinfo3}else{
  hm
}
}

#' Build segregated cellinfo and plot a signature heatmap
#'
#' Annotates cells by `vars`, segregates groups, optionally adds module scores, draws a ComplexHeatmap via [cellinfo.heatmap()], and returns either the updated cellinfo or the heatmap.
#'
#' @param so A Seurat object.
#' @param cellinfo A cellinfo list object.
#' @param vars Variables used to segregate cells.
#' @param colorlist Named list of annotation color vectors. Defaults to global `allcolors` when available.
#' @param extra.genes Additional gene labels for the heatmap.
#' @param extra.annotation Optional extra annotation variables to import.
#' @param assay Assay name used for expression.
#' @param nmarkers Number of markers per group to label.
#' @param return.heatmap If `TRUE`, return the heatmap object instead of cellinfo.
#'
#' @return Updated `cellinfo` or a heatmap object when `return.heatmap` is `TRUE`.
#'
#' @seealso [heatmap.workflow()], [cellinfo.heatmap()]
#'
#' @export
cellinfo.group.signatures <- function(so, cellinfo,vars, colorlist=allcolors, extra.genes=NULL, extra.annotation=NULL, assay="SCT", nmarkers=3, return.heatmap=F){
  cellinfo2=add.cell.annotation(so, cellinfo=cellinfo , vars = vars, assay=assay)
  cellinfo3=cellinfo2 %>% segregate.cells(., variables = c(cellinfo$cell.group.variable, vars)) %>% cellinfo.order.groups.natural

  #get module scores per cell

  so <- AddModuleScore3(clean.module.scores(so), cellinfo3$markerlist )

  cls=c(grepl(rownames(so), "mscore_group", value=T),cellinfo$cell.group.variable, vars)
 so@meta.data[, cls ] %>% pivot_longer(.,  grepl(rownames(so), "mscore_group", value=T), names_to="signature", values_to="value") %>%


  if(is.null(colorlist[[cellinfo3$cell.group.variable]])){

    colorlist[[cellinfo3$cell.group.variable]]<-rainbow(length(cellinfo3$cell.list)) %>% givename(., names(cellinfo3$cell.list))
    allcolors[[cellinfo3$cell.group.variable]]<<-rainbow(length(cellinfo3$cell.list)) %>% givename(., names(cellinfo3$cell.list))
  }
  fcat(dput(colorlist[[cellinfo3$cell.group.variable]]))

  if(!is.null(extra.annotation)){
    cellinfo3=cellinfo3 %>% add.cell.annotation(so, ., extra.annotation, assay=assay, overwrite=T)
  }
  #allcolors[[cellinfo3$marker.group.variable]]=rainbow(length(cellinfo3$markerlist)) %>% givename(., names(cellinfo3$markerlist))

  cellinfo3$cell.annotation[[cellinfo3$cell.group.variable]]=factor(cellinfo3$cell.annotation[[cellinfo3$cell.group.variable]], levels=ifelse(!is.null(colorlist[[cellinfo3$cell.group.variable]]), colorlist[[cellinfo3$cell.group.variable]] %>% names, gtools::mixedsort(unique(cellinfo3$cell.annotation[[cellinfo3$cell.group.variable]]))))


  hm=cellinfo.heatmap(so=so, cellinfo=cellinfo3, assay=assay, genes.to.label = c((cellinfo3 %>% adjust.markers(., nmarkers))$markers, landmark.genes, extra.genes) , colorlist=colorlist)

  tpdf(paste0("heatmap_",Project(so),"_by_", cellinfo$cell.group.variable),  sca=6)
  print(hm)
  dev.off()

  if(!return.heatmap){
    cellinfo3}else{
      hm
    }
}
