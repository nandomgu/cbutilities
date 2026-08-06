################################################################################
# cellinfo utility functions
# Extracted from mlutils_setup2.R: top-level functions that reference cellinfo.
# Ordered alphabetically by function name.
################################################################################

add.cell.annotation= function(so, cellinfo, vars, assay=NULL, overwrite=T){
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


adjust.cells= function(cellinfo, n){
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


adjust.markers= function(cellinfo, n){
  cellinfo$markerlist= adjust.markerlist(cellinfo$markerlist, n)
  cellinfo=refresh.cellinfo(cellinfo)
  cellinfo
}


cellinfo.filter.cells=function(cellinfo, varr, vals, so=NULL, assay=NULL){

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


cellinfo.get.values= function(cellinfo, variable, so){
  
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


cellinfo.group.signatures=function(so, cellinfo,vars, colorlist=allcolors, extra.genes=NULL, extra.annotation=NULL, assay="SCT", nmarkers=3, return.heatmap=F){
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


cellinfo.heatmap= function(so,
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


cellinfo.init= function(so, seriate=T, ncells=100,clusvar="seurat_clusters",markerset=clusvar, ...){
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


cellinfo.order.groups.natural=function(cellinfo){
cellinfo$cell.list=cellinfo$cell.list[cellinfo$cell.list %>% names %>% mixedorder]
cellinfo$markerlist=cellinfo$markerlist[cellinfo$markerlist %>% names %>% mixedorder]
refresh.cellinfo(cellinfo)
}


cellinfo.ordered.cell.plot= function(cellinfo, variable, scale.x=T, clusters=names(cellinfo$cell.list), colorlist=NULL, groupvar=cellinfo$clustering.variable){
  
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


cellinfo.plotfeatures=function(so, cellinfo, markers.per.group=3,umap.number=NULL){
  

   nms=names(so.mes@reductions)
   unms=nms[grepl("umap", nms)]
   
  if(is.null(umap.number)){
   rr=unms[length(unms)]# using the last umap generated
  }else{
   rr=unms[umap.number] 
  }
  
plt=FeaturePlot(so.mes, features=adjust.markers(cellinfo, markers.per.group)$markers, reduction=rr)  

plt}


cellinfo.pointplot=function(so, cellinfo, ncells=1000000,groupby=cellinfo$cell.group.variable,  assay=DefaultAssay(so)){
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


cellinfo.seriatecells.seurat=function(cellinfo, so , new.cells=T, ncells=100, return.scores=F){
    
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


cellinfo.workflow.singleds=function(so, grouping.variable="seurat_clusters", assay="RNA", markers=NULL, recalculate.markers=F, marker.params=NULL, colorlist=NULL, ncells=300, return.seurat=F){
  
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


heatmap.workflow=function(so, cellinfo,vars, colorlist=allcolors, extra.genes=NULL, extra.annotation=NULL, assay="SCT", nmarkers=3, return.heatmap=F, label=""){
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


keep.cellgroups.pattern=function(cellinfo, pattern){
 cellinfo$cell.list = cellinfo$cell.list[grepl(pattern, names(cellinfo$cell.list))] 
refresh.cellinfo(cellinfo)
 }


plot.signature.profiles=function(so, cellinfo, colorlist=allcolors){
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


prepare.cellgroup.labels=function(cellinfo){
 
  label.list=lapply(names(cellinfo$cell.list), function(name) {
    vec_length <- length(cellinfo$cell.list[[name]])
    middle <- ceiling(vec_length / 2)
    out <- rep(NA_character_, vec_length)
    out[middle] <- name
    out
  })
  
  label.list %>% Reduce(c, .)
  
  }


refresh.cellinfo=function(cellinfo, cell.group.label="cell.group", marker.group.label="marker.group"){
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


refresh.cellinfo.hard=function(cellinfo, cell.group.label="cell.group", marker.group.label="marker.group"){
  

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


remove.cell.annotation=function(cellinfo, variables){
  
  cellinfo$cell.annotation= cellinfo$cell.annotation %>% select(-all_of(variables)) 
  
  cellinfo$cell.metadata=  setdiff(colnames(cellinfo$cell.annotation), variables)
  cellinfo 
}


remove.marker.annotation=function(cellinfo, variables){
  
  cellinfo$marker.annotation= cellinfo$marker.annotation %>% select(-all_of(variables)) 
  
  cellinfo$marker.metadata=  setdiff(colnames(cellinfo$marker.annotation), variables)
  cellinfo 
}


removemarkers=function(cellinfo, markers){
  nms=names(cellinfo$markerlist)
  cellinfo$markerlist=lapply(1:length(cellinfo$markerlist), function(x){
    gns=cellinfo$markerlist[[x]]
    
    gns[!(gns %in% markers)]
    
  }) %>% givename(., nms)
  
  refresh.cellinfo(cellinfo)
}


segregate.cells=function(cellinfo, variables){
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


seriategenes=function(so, cellinfo, assay="SCT", seriate.groups=F, link.markers.cells=T){
  
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


subset.cellinfo=function(cellinfo, labels){
  fcat("restricting to labels only present in list")
  tf=lapply(labels, function(x) x %in% (cellinfo$markerlist %>% names)) %>% Reduce(c, .)
  labels2=labels[tf]
  
  cellinfo$markerlist=cellinfo$markerlist[labels2]
  cellinfo$cell.list=cellinfo$cell.list[labels2]
  cellinfo$input.cluster.ids=labels2
  refresh.cellinfo(cellinfo)
}


update.group.vars=function(cellinfo, global=NULL, cells.group=NULL, markers.group=NULL){
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
