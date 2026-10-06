#' Filter cellinfo cells by an annotation value set
#'
#' Retains only cells whose annotation column `varr` is in `vals`. If the column is missing and a Seurat object is provided, the column is imported first.
#'
#' @param cellinfo A cellinfo list object.
#' @param varr Annotation column name used for filtering.
#' @param vals Values of `varr` to keep.
#' @param so Optional Seurat object used to import missing annotations.
#' @param assay Assay forwarded to [add.cell.annotation()] when importing.
#'
#' @return Filtered and refreshed `cellinfo` list.
#'
#' @export
cellinfo.filter.cells <- function(cellinfo, varr, vals, so=NULL, assay=NULL){

  if(!(varr %in% colnames(cellinfo$cell.annotation))){
    
  
  if (!is.null(so)){
    fcat("seurat object detected. importing annotations")
    cellinfo= cellinfo %>% add.cell.annotation(so, ., vars = varr, assay=assay)
  }else{
   fcat(varr, "not detected in cellinfo and no seurat object provided. please provide seurat object. returning unfiltered") 
    return(cellinfo)
  }}
    
  
   
    filtered.cells= cellinfo$cell.annotation %>% dplyr::filter(!!sym(varr) %in% vals ) %>% rownames
    
    cellinfo$cell.list=lapply(cellinfo$cell.list, function(x){
      x[x %in% filtered.cells]
      
    })
    
  refresh.cellinfo(cellinfo)   
    
  }

#' Extract per-group metadata or expression values
#'
#' Fills `cellinfo$values[[variable]]` with either metadata or expression values for each group in `cellinfo$cell.list`.
#'
#' @param cellinfo A cellinfo list object.
#' @param variable Metadata column or gene name.
#' @param so Seurat object providing values.
#'
#' @return Updated `cellinfo` with `values` populated for `variable`.
#'
#' @export
cellinfo.get.values <- function(cellinfo, variable, so){
  
  metacols= so %>% metadata  %>% colnames
  genenames= rownames(so)
  
  if(variable %in% metacols){
    
    for(xx in names(cellinfo$cell.list)){
      
      cellinfo$values[[variable]][[xx]]= (so %>% metadata)[cellinfo$cell.list[[xx]], variable] 
    }
    
  }
  
  
  if(variable %in% genenames){
    
    
    metamat=join_meta_exp(so, genes=variable, assay=DefaultAssay(so))
    for(xx in names(cellinfo$cell.list)){ 
      cellinfo$values[[variable]][[xx]]=  metamat[cellinfo$cell.list[[xx]],variable ]
      
      
    }
    
  }
  cellinfo
}
