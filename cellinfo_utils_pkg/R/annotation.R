#' Add cell annotation columns to a cellinfo object
#'
#' Imports metadata and/or gene expression values from a Seurat object into `cellinfo$cell.annotation` for the cells already represented in `cellinfo`.
#'
#' @param so A Seurat object.
#' @param cellinfo A cellinfo list object.
#' @param vars Character vector of metadata columns and/or gene names to add.
#' @param assay Assay to use for gene values. Defaults to `DefaultAssay(so)`.
#' @param overwrite Logical; if `FALSE`, only new annotation columns are added.
#'
#' @return Updated `cellinfo` list.
#'
#' @seealso [remove.cell.annotation()], [refresh.cellinfo()]
#'
#' @export
add.cell.annotation <- function(so, cellinfo, vars, assay=NULL, overwrite=T){
if(is.null(assay)){
 assay=DefaultAssay(so) 
  
}
  
  gene.vars=vars[vars %in% rownames(so)]
  non.gene.vars=vars[!(vars %in% rownames(so))]
  meta.vars= non.gene.vars[non.gene.vars %in% (metadata(so) %>% colnames)]
  non.meta.vars= non.gene.vars[!(non.gene.vars %in% (metadata(so) %>% colnames))]
  
  absent.vars=intersect(non.gene.vars, non.meta.vars)
  if(length(absent.vars)>0){
    fcat("Warning: the following variables are absent from the dataset: ", paste(absent.vars, collapse=" "))
  }
  
  clls=cellinfo$cell.annotation %>% rownames
  
  
  new.annotation=metadata(so)[clls,] %>% select(all_of(meta.vars))
  #%>% givecolnames(., nms=meta.vars)
  
  #if overwrite is false, we only incorporate new columns
  if(overwrite==F){
    addition.vars=setdiff(meta.vars, cellinfo$metadata)
  }else{
    #if overwrite is true then we allow replacement of old variables
    addition.vars=meta.vars 
  }
  
  for(mv in addition.vars){
    cellinfo$cell.annotation[[mv]]= new.annotation[[mv]]
  }
  
  
  if(length(gene.vars)>0){
    fcat("alert0")
    genemat=join_meta_exp2(so, genes=gene.vars, cells=clls, assay=assay, layer="data")
    for(gn in gene.vars){
      genemat=fillmat(genemat, gn)
    }
    #%>%  givecolnames(., nms=gene.vars)
    cellinfo$cell.annotation=cellinfo$cell.annotation %>% cbind(., genemat[clls,] %>% select(all_of(gene.vars)) )
  }
  
  cellinfo$cell.metadata=colnames(cellinfo$cell.annotation)
  cellinfo
  
  
}

#' Remove columns from cellinfo cell annotations
#'
#' Drops the requested columns from `cellinfo$cell.annotation` and updates `cell.metadata`.
#'
#' @param cellinfo A cellinfo list object.
#' @param variables Character vector of annotation columns to remove.
#'
#' @return Updated `cellinfo`.
#'
#' @export
remove.cell.annotation <- function(cellinfo, variables){
  
  cellinfo$cell.annotation= cellinfo$cell.annotation %>% select(-all_of(variables)) 
  
  cellinfo$cell.metadata=  setdiff(colnames(cellinfo$cell.annotation), variables)
  cellinfo 
}

#' Remove columns from cellinfo marker annotations
#'
#' Drops the requested columns from `cellinfo$marker.annotation` and updates `marker.metadata`.
#'
#' @param cellinfo A cellinfo list object.
#' @param variables Character vector of annotation columns to remove.
#'
#' @return Updated `cellinfo`.
#'
#' @export
remove.marker.annotation <- function(cellinfo, variables){
  
  cellinfo$marker.annotation= cellinfo$marker.annotation %>% select(-all_of(variables)) 
  
  cellinfo$marker.metadata=  setdiff(colnames(cellinfo$marker.annotation), variables)
  cellinfo 
}

#' Rename primary grouping variables in a cellinfo object
#'
#' Updates clustering / cell / marker group labels and the corresponding first annotation column names.
#'
#' @param cellinfo A cellinfo list object.
#' @param global If provided, set both cell and marker primary labels to this value.
#' @param cells.group Optional new primary cell-group label.
#' @param markers.group Optional new primary marker-group label.
#'
#' @return Updated `cellinfo`.
#'
#' @name update.group.vars
#' @usage update.group.vars(cellinfo, global = NULL, cells.group = NULL, markers.group = NULL)
#' @export update.group.vars
update.group.vars <- function(cellinfo, global=NULL, cells.group=NULL, markers.group=NULL){
  if(!is.null(global)){
    cellinfo$clustering.variable=global
    cellinfo$cell.metadata[1]=global
    cellinfo$marker.metadata[1]=global
    subst=colnames(cellinfo$cell.annotation)
    subst[1]=global
    colnames(cellinfo$cell.annotation)=subst
    substm=colnames(cellinfo$marker.annotation)
    substm[1]=global
    colnames(cellinfo$marker.annotation)=substm
    
  }else{
    
    if(!is.null(cells.group)){
      cellinfo$clustering.variable=cells.group
      cellinfo$cell.metadata[1]=cells.group
      
      subst=colnames(cellinfo$cell.annotation)
      subst[1]=cells.group
      colnames(cellinfo$cell.annotation)=subst
      
    }
    
    if(!is.null(markers.group)){
      
      
      cellinfo$marker.metadata[1]=markers.group
      
      substm=colnames(cellinfo$marker.annotation)
      substm[1]=markers.group
      colnames(cellinfo$marker.annotation)=substm
      
    }
    
  }
  cellinfo
  
}
