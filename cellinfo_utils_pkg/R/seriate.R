#' Reorder cells within groups using Seurat module scores
#'
#' Scores each group with its markers via `AddModuleScore3()`, then reorders (and optionally resamples) cells in `cellinfo$cell.list` by score.
#'
#' @note Requires `AddModuleScore3()` and `get.module.name()` from the parent mlutils codebase.
#'
#' @param cellinfo A cellinfo list object.
#' @param so A Seurat object.
#' @param new.cells If `TRUE`, take the top `ncells` by score; if `FALSE`, reorder existing listed cells.
#' @param ncells Number of cells to keep per group when `new.cells` is `TRUE`.
#' @param return.scores If `TRUE`, also return module-score metadata.
#'
#' @return Refreshed `cellinfo`, or a list with cellinfo and scores.
#'
#' @seealso [seriategenes()], [refresh.cellinfo.hard()]
#'
#' @export
cellinfo.seriatecells.seurat <- function(cellinfo, so , new.cells=T, ncells=100, return.scores=F){
    
          markerlist.original=cellinfo$markerlist
          fcat("cleaning up markers...")
          markerlist=cellinfo$markerlist[lapply(cellinfo$markerlist, function(x) length(x)>0) %>% Reduce(c, .)  ]
          
          sdi=setdiff(names(markerlist.original), names(markerlist))
          if(length(sdi)>0){
           fcat("no markers were found for the following clusters:", paste(sdi, collapse=",")) 
          }
            celltotals=list()
            

            fcat("Using Seurat Module Score method to rank cells...")
            
            temp.metadata=so@meta.data
            
            so <- AddModuleScore3(so, markerlist )
            
            
            ordcellist=lapply(1:length(markerlist), function(x) {
               fcat("rearranging cells of", names(markerlist)[x])
              modulename= get.module.name(names(markerlist)[x])
              fcat("modulescore name", modulename)
              #sumcells=apply(so[markerlist[[x]],so %>% metadata %>% pull(), 2, sum)
              fcat("available module scores are", so@meta.data %>% select(contains("mscore_group")) %>% colnames %>% paste(., collapse=","))
              fcat("available marker names are", paste(names(markerlist) , collapse=","))
              fcat("cellinfo$cell.group.variable is", cellinfo$cell.group.variable)
              val=names(markerlist)[x]
              all.ordered.cells<<-so@meta.data %>% dplyr::filter(!!sym(cellinfo$cell.group.variable)==!!val) 
              all.ordered.cells= all.ordered.cells %>% dplyr::select(!!sym(modulename)) %>% arrange(!!sym(modulename)) 
               
              fcat("made allorderedcells", modulename)
              all.ordered.cells=all.ordered.cells %>% rownames
              if(!new.cells){
               fcat("no new cells for", modulename)
                # here we retrieve the reordered cells on the original list
               
                outputcells=all.ordered.cells[all.ordered.cells %in% cellinfo$cell.list[[names(markerlist)[x]]]]
                
              }else{
                l=length(all.ordered.cells)
               # retrieve the strongest cells in each cluster 
                if(l<ncells){ncells=l}
                  
                 outputcells=all.ordered.cells[(l-ncells+1):l] 
                
                
              }
              
              #celltotals[[x]]<<- rep(x, length(outputcells))
             outputcells
            }
            ) %>% givename(., names(markerlist))
            
            ####################################################################
            # after reordering cells for the clusters fr which there are markers,
            #we modify the cells for those specific clusters, leaving the
            #possibility to keep some cell clusters untouched if there are no
            # no corresponding markers for them
            ####################################################################
          
            fcat("updating new cell information on cellinfo...")
            for(j in 1:length(ordcellist)){
            
           cellinfo$cell.list[[names(ordcellist)[j]]]=ordcellist[[j]] 
            
          }
          ## markerlist has been cleaned to remove all sorts of faulty marker sets  so we replace it
          fcat("updating revised marker sets on cellinfo...")
            cellinfo$markerlist=markerlist  
            
          if(return.scores){
           return(list(cellinfo=cellinfo, module.scores=so@meta.data %>% select(all_of(paste0(newnames, 1:length(newnames)))))) 
          }else{
           return(refresh.cellinfo.hard(cellinfo, cell.group.label=cellinfo$cell.metadata[1], marker.group.label=cellinfo$marker.metadata[1]))
          }
            
}

#' Seriate markers within each group using PCA seriation
#'
#' Orders genes inside each marker set based on scaled expression patterns over the linked cells (or all cells). Optionally reorders marker groups.
#'
#' @param so A Seurat object.
#' @param cellinfo A cellinfo list object.
#' @param assay Assay providing the data matrix.
#' @param seriate.groups If `TRUE`, also reorder marker groups by a summary score.
#' @param link.markers.cells If `TRUE`, seriate each marker set using its paired cell group.
#'
#' @return Updated `cellinfo` with seriated markers.
#'
#' @seealso [cellinfo.seriatecells.seurat()]
#'
#' @export
seriategenes <- function(so, cellinfo, assay="SCT", seriate.groups=F, link.markers.cells=T){
  
  #mat=tryCatch({so[cellinfo$markers, ][[assay]]@data %>% as.matrix}, function(e) {
  # warning("problem extracting assay data matrix. extracting counts instead")
  #so[st, dscells][["RNA"]]@counts %>% as.matrix
  #})
  
  mat=so[cellinfo$markers, ][[assay]]@data %>% as.matrix
  absent.genes= setdiff(cellinfo$markers, rownames(mat))
  vr= get.variant.rows(mat) 
  
  cellinfo= removemarkers(cellinfo, union(mat[!vr,] %>% rownames, absent.genes))
  mat=mat[vr, ] 
  clusmeans=c()
  ct=1
  fcat("Seriating genes")
  ordered.gene.list=lapply(1: length(cellinfo$markerlist), function(nn){
    
    st= cellinfo$markerlist[[nn]]
    if(!link.markers.cells){
      dscells=colnames(mat)
    }else{
      dscells=cellinfo$cell.list[[nn]]
    }
    if(length(st)>1){
      ss=seriation::get_order(
        #assumes an ordered correspondence between cell markers and cells in cellinfo
        seriation::seriate(scale((mat[st, dscells]) %>% t) %>%t)
        , method="PCA")
      clusmeans[ct]<<- mean(colMeans(scale((mat[st, dscells]) %>% t) %>%t)*(sapply(1:length(dscells), function(x) x**2, USE.NAMES=F)))
      ct<<-ct+1
      st[ss %>% unname]}else{
        
        clusmeans[ct]<<- mean(scale((mat[st, dscells]) %>% as.numeric)*(sapply(1:length(dscells), function(x) x**2, USE.NAMES=F)))
        
        st
      }
  }) %>% givename(., cellinfo$markerlist %>% names)
  fcat("clusmeans", clusmeans)
  

  if(seriate.groups){
    names(clusmeans)= 1:length(clusmeans)
    neworder=as.integer(clusmeans %>% sort %>% names)
    newnames=names(cellinfo$markerlist)[neworder]
    repmt<<-lapply(neworder, function(tt) ordered.gene.list[[tt]]) %>% givename(., newnames)
    cellinfo$markerlist= repmt
  }else{
    
    cellinfo$markerlist=ordered.gene.list 
  }
  
  cellinfo.final=refresh.cellinfo(cellinfo)
  
  cellinfo.final
}

#' Order cellinfo groups with natural/mixed sorting
#'
#' Reorders `cell.list` and `markerlist` using `gtools::mixedorder()` on their names, then refreshes derived fields.
#'
#' @param cellinfo A cellinfo list object.
#'
#' @return Reordered and refreshed `cellinfo`.
#'
#' @export
cellinfo.order.groups.natural <- function(cellinfo){
cellinfo$cell.list=cellinfo$cell.list[cellinfo$cell.list %>% names %>% mixedorder]
cellinfo$markerlist=cellinfo$markerlist[cellinfo$markerlist %>% names %>% mixedorder]
refresh.cellinfo(cellinfo)
}
