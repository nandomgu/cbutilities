sca=2.5
fcat= function(...) cat(paste(...,sep=" "), "\n")
msg=fcat
################################################################################
# load libraries and packages
################################################################################
library(reticulate)
use_condaenv("scvi", required = TRUE)
py_config()


packages <- c("devtools", "dplyr","data.table","tidyr", "ggplot2","simpleCache","DESeq2","Seurat", 
               "ggfortify", "ggExtra", "gridExtra", "ggforce", "ggrastr", "alluvial", "mltools", "GSVA", "tinytex", "hdf5r", "glmGamPoi", "gtools", "patchwork", "STRINGdb", "pcaMethods",
              #"velocyto-team/velocyto.R",
               "AUCell", "RcisTarget","GENIE3", 
              "zoo", "mixtools", 
              #"rbokeh",
              #"DT", "NMF",
              "ComplexHeatmap",
              #"R2HTML", 
              "Rtsne",
              #"doMC", "doRNG", 
              #"aertslab/SCopeLoomR", 
              "sccore", 
              #"canceRbits/DElegate", 
              "harmony", "zellkonverter", "viridis", "concaveman"
              )

# Function to check, install if necessary, and load packages

dir.create( "~/mnt_out/containerlibrary/", recursive=T)
.libPaths(c(.libPaths(), "~/mnt_out/containerlibrary/"))

################################################################################
# installing and loading packages decision tree
################################################################################

install_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    # Try installing with BiocManager first
    tryCatch({
      BiocManager::install(
        pkg,
        dependencies = TRUE,
        update = FALSE,
        lib = "~/mnt_out/containerlibrary/"
      )
    }, error = function(e) {
      message(paste(
        "BiocManager could not install", pkg,
        "- trying install.packages instead."
      ))
      install.packages(
        pkg,
        dependencies = TRUE,
        lib = "~/mnt_out/containerlibrary/"
      )
    })
    
    # Try loading the package after installation
    tryCatch(
      {
        library(pkg, character.only = TRUE)
      },
      error = function(e) {
        warning(
          "package could not be loaded. check for a different name",
          call. = FALSE
        )
      }
    )
  }
}
# Apply the function to each package
devtools::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
sapply(packages, install_and_load)



################################################################################
#formatting functions
################################################################################

remove.pipes=function(x) gsub("\\|","-",x)

givenames=
  function(x, nam) {
    names(x)= nam
    x
  } 

givecolnames=function(x, ind, nms){
  colnames(x)[ind]= nms
  x
  
}


givename=
  function(x, nam) {
    names(x)= nam
    x
  } 
          
################################################################################
# assembling workspace variables
################################################################################

prep.dir.for.out.root= function(x) gsub("/home/rstudio/mnt_out|~/mnt_out", config$out_root_host, x)

prep.dir.for.in.root= function(x) gsub(config$out_root_host, config$out_root, x)

prep.dir.for.project.root= function(x) gsub("~", config$project_root_host, x)

make.bigwigpath=function(rootdir) paste0(rootdir, "bwa/merged_library/bigwig/")
make.profilepath= function(x) paste0(x, "/bwa/merged_library/deeptools/plotprofile/")

allcolors=list(stdcolors = c(`2` = "#bebada", `4` = "#fccde5", `19` = "#80b1d3", 
                             `13` = "#FEEDC3FF", `8` = "#ccebc5", `11` = "#FF7BBAFF", `7` = "#E6F288FF", 
                             `1` = "#B5C8FAFF", `20` = "#ffed6f", `14` = "#FDA6D8FF", `18` = "#999999", 
                             `15` = "#d9d9d9", `6` = "#b3de69", `12` = "#fdb462", `16` = "#bc80bd", 
                             `5` = "#B4F9FDFF", `9` = "#fb8072", `3` = "#EE8B6EFF"), 
               #celltype = c(HSC = "#ccebc5", pHSC = "#EE8B6EFF",`CD34 Cord Blood` = "#fb8072",`CD34 Bone Marrow` = "#d9d9d9",
                # MPP = "#FDA6D8FF", LMPP = "#bebada",  GMP = "#80b1d3",Mono = "#fccde5", CMP = "#FEEDC3FF", MEP = "#E6F288FF",
                # Ery = "#bc80bd",CLP = "#b3de69",Bcell = "#B4F9FDFF",CD4Tcell = "#999999",  CD8Tcell = "#B5C8FAFF",  NKcell = "#fdb462", 
  #Round1 = "red", Round2 = "blue"), 
  celltype=c(HSC="#CCEBC5", MPP="#008000", `CD34 Cord Blood` = "#32CD32", `CD34 Bone Marrow` = "#32CD32",  CLP="#1BCED0", LMPP="#446BDE", Bcell="#B4F9FD",CD4Tcell="#619EA0", CD8Tcell="#B5C8FA",
    NKcell="#0747A3",   GMP="#80B1D3", Mono="#F7F7C6",
    CMP="#FF8C00", MEP= "#FFD700" , Ery = "#bc80bd", Round1 = "red", Round2 = "blue"),
               condition_replicate = c(`1-r1` = "#a6cee3", 
                                       `1-r2` = "darkblue", `2-r1` = "#b2df8a", `2-r2` = "#33a02c", 
                                       `3-r1` = "#fb9a99", `3-r2` = "#e31a1c", `4-r1` = "cyan", `4-r2` = "cyan4", 
                                       `5-r1` = "#cab2d6", `5-r2` = "#6a3d9a", `6-r1` = "#ffff99", `6-r2` = "brown4", 
                                       `7-r1` = "darkgoldenrod1", `7-r2` = "darkgoldenrod", `8-r1` = "deeppink", 
                                       `8-r2` = "#bd2398"),
               is.target = c(test = "blue", Ery = "red", 
                             other = "grey"), condition=c(CTRL = "#FF0000", TSR3kd = "#00FFFF")
  )


                     


################################################################################
#function to import all variables from the config to make them arrays when they are separated by commas
################################################################################


config.to.vectors=function(config){
dataset.info=lapply(names(config), function(x){
  #fcat(x)
  if(!is.null(config[[x]]) && is.character(config[[x]])){
  strsplit(config[[x]], split=",")[[1]]
  }else{
   config[[x]]
  }
}) %>% givename(., names(config))
  
dataset.info$is_reference=as.logical(dataset.info$is_reference) 
dataset.info
}                               

################################################################################
#actually import the config file and create vectors when needed, so all other 
#functions have access to the config
################################################################################

config <- yaml::read_yaml(file = 'config.yaml')

# modify the line below if the cachedir must be changed on the go

#config$cachedir<- config$cachedir

setCacheDir(config$cachedir)
## transform config variables into vectors based on comma separation  
dataset.info=config.to.vectors(config)



################################################################################
#add current version to a character string
################################################################################

addversion= function(x) paste_(x, "version", analysis.version)


addplotpath= function(x) paste0(config$plotpath, x)

################################################################################
# Labeling arrays and matrices
################################################################################




################################################################################
#force a named vector into a matrix
################################################################################
coerce.vector.to.matrix=function(vec)  matrix(vec, ncol = 1, dimnames = list(names(vec), NULL))

################################################################################
# global image parameters
################################################################################

w=250
pw=2.5
rr=600

################################################################################
#timestamp png or pdf plots
################################################################################

timestamp= function(nm, date.use=NULL)  ifelse(is.null(date.use), paste0(Sys.Date(),"_", nm), paste0(date.use,"_", nm))



################################################################################
# add prefix and suffix to a name vector
################################################################################

addprefix= (function(cc, prefix="C")  paste0(prefix, cc)) %>% Vectorize


################################################################################
# modifying the plot labels in ggplots in  a simple fashion
################################################################################

#resizing overall text in the plot
resizetext=function(x) theme_classic(base_size = x)

rotatex= function(x) theme(axis.text.x = element_text(angle = x, hjust=1, vjust=0.5))


################################################################################
# project-specific colors and cell types
################################################################################


################################################################################
# Trimming column names from a feature counts table in which the names have a long bam file suffix
################################################################################

cleannames <- function(a, sep="_"){
  ncn<- a %>%  lapply(., function(x) {strsplit(x, split=sep)[[1]][1]})  %>% unlist
  ncn
}


cleancolnames <- function(a){
  ncn<- a %>% colnames %>% lapply(., function(x) {strsplit(x, split="_")[[1]][1]})  %>% unlist
  colnames(a)<- ncn
  a
}


################################################################################
#FUNCTIONS TO EDIT METADATA FIELDS
################################################################################


assignhealth<- Vectorize(function(x){
  out=0
  if(any(grepl('Donor|4983|6792|CB102|BM-CD34', x))) {return("Healthy"); out=1} 
  if(any(grepl('SU', x))) {return("AML"); out=1} 
  if(out==0){
    return(x) 
  }
}, USE.NAMES=F)

assignhealth2<- Vectorize(function(x){
  out=0
  if(any(grepl('leukem', x))) {return("AML"); out=1} 
  else {return("Healthy"); out=1} 
  if(out==0){
    return(x) 
  }
}, USE.NAMES=F)



assigncelltype<- Vectorize(function(x){
  out=0
  if(any(grepl('HSC', x))) return("HSC")
  if(any(grepl('Erythroblast', x))) {return("Erythroblast"); out=1} 
  if(any(grepl('Bcell', x))) {return("Bcell") ; out=1}
  if(any(grepl('LSC', x))) {return("LSC") ; out=1}
  if(any(grepl('CLP', x))) {return("CLP") ; out=1}
  if(any(grepl('CMP', x))) {return("CMP") ; out=1}
  if(any(grepl('Leuk', x))) {return("Leuk") ; out=1}
  if(any(grepl('CD8', x))) {return("CD8_Tcell") ; out=1}
  if(any(grepl('CD4', x))) {return("CD4_Tcell") ; out=1}
  if(any(grepl('Nkcell', x)) || any(grepl('NK', x))) {return("NKcell"); out=1}
  if(any(grepl('GMP', x))) {return("GMP") ; out=1}
  if(any(grepl('MPP', x))) {return("MPP") ; out=1}
  if(any(grepl('MEP', x))) {return("MEP") ; out=1}
  if(any(grepl('CD34', x))) {return("CD34") ; out=1}
  if(out==0){
    return("unlabeled") 
  }
}, USE.NAMES=F)

pcadf=function(sc){
  cbind(metadata(sc), Embeddings(sc, reduction="pca"))
}

getdate= function() format(Sys.time(), "%Y%b%d")


setlevels=function(so, varr, levls){
  Idents(so)=varr
  so@active.ident <- factor(x = so@active.ident, levels = levls)
  so}


col2names= function(dff, coll){
  rownames(dff)= dff[, coll]
  dff
}


names2col=function(dff, coln="gene") {
  dff[, coln]= rownames(dff) 
  dff
}



paste_=function(...) paste(..., sep="_")


showpalette=function(p, pname="generic"){
  sca=10
  rr=600
  w=250
  tpng(paste_("palette", pname), res=rr, he=sca*w, wi=sca*w*2)
  
  pl=ggplot(1:length(p) %>% as.data.frame %>% dplyr::mutate(h=1, colors=p, name=names(p)), aes(x=factor(.), y=factor(h), label=colors))+
    geom_col(aes(fill=factor(.)))+
    scale_fill_manual(values=namenums(p))+
    theme_classic()+
    coord_cartesian(ylim=c(0,1))+
    geom_text(y=0.5, color="black", angle=90)+
    #NoAxes()+
    NoLegend()+
    ggtitle(pname)
  
  if(is.null(names(p))){
    print(pl )
  }else{
    print(pl+geom_text(inherit.aes=T,aes(label=name), y=1.1, color="black", angle=90))  
  }
  dev.off()
  
  return(pl+geom_text(inherit.aes=T,aes(label=name), y=1.1, color="black", angle=90))
}




#getsourcename<- Vectorize(function(xp) metadata %>% as.data.frame %>% filter(Experiment==xp) %>% pull(source_name), USE.N=F)
#sranames=lapply(ddsexpts, getsourcename) %>% Reduce(c, .)

#getbioname=Vectorize(function(xp) samplesheet %>% filter(experiment_accession==xp) %>% pull(sample_description), USE.N=F)
#bionames=lapply(ddsexpts, getbioname) %>% Reduce(c, .)
#colnames(dds)= editname(colnames(dds))


###########################
#Below: imports from scutils_lean
###########################



##########################################################################################################################
#general functions
##########################################################################################################################
###########################################################################################################################
cat("adding general functions...\n")
#permutate an array
shuffle=function(x) sample(x, length(x), replace=F)

#prepare gene DE table by removing NAs, filtering by cutoff p-value
filterorderDE<-function(x, thresh=0.05, lfthresh=1, netchange=FALSE, top=NA){
  
  if(netchange==TRUE){
    x=x[is.na(x[,pvaluecol])==FALSE & x[,pvaluecol ]<thresh & abs(x[,lfccol])>lfthresh,]
    ord=order(abs(x[,lfccol]), decreasing=T)
  }else{
    x=x[is.na(x[,pvaluecol])==FALSE & x[,pvaluecol ]<thresh & x[,lfccol]>lfthresh,]
    ord=order(x[, lfccol], decreasing=TRUE)
  }
  if(!is.na(top)){
    x[ord[1:top], ]
  }else{
    x[ord, ]
  }
}


reddf=function(sc, reduc="pca"){
  cbind(metadata(sc), Embeddings(sc, reduction=reduc))
}

harmonydf=function(sc){
  cbind(metadata(sc), Embeddings(sc, reduction="harmony"))
}

umapdf=function(sc, reductions=c("umap"), change.names=NULL){
  cbind(metadata(sc), lapply(1:length(reductions), function(x) Embeddings(sc, reduction=reductions[x])) %>% Reduce(cbind, .)   )
  
  
}

pcadf=function(sc){
  cbind(metadata(sc), Embeddings(sc, reduction="pca"))
}

getdate= function() format(Sys.time(), "%Y%b%d")


setlevels=function(so, varr, levls){
  Idents(so)=varr
  so@active.ident <- factor(x = so@active.ident, levels = levls)
  so}


col2names= function(dff, coll){
  rownames(dff)= dff[, coll]
  dff
}

condrbind=function(x,y){
  
  if(is.null(x)){
    return(y)}else{
      if(is.null(y)){return(x)}
      else{
        return(rbind(x,y))
      }}
  
}

condcbind=function(x,y){
  
  if(is.null(x)){
    return(y)}else{
      if(is.null(y)){return(x)}
      else{
        return(cbind(x,y))
      }}
  
}



#recursive function to merge a bunch of clusters at the same time
# function takes several vectors c(a,b,c...), c(d,e,f)... so that clusters in each vector will be merged together into a single category. 
#output: a final vector defining new labels for the merged clusters. 

triagecells.multi= function(so, classvar="seurat_clusters", ...){
  grps=list(...)# assemble groups into a list
  
  allgroups=lapply(grps, function(x) combineclusters2(so, x))
  mergername=function(x, refnames){
    if(x==0){
      return(refnames)
    }else{
      refnames[allgroups[[x]]]=LETTERS[x]
      return(mergername(x-1, refnames))
    }
  }
  out=as.character(so@meta.data[, classvar])
  final=mergername(length(allgroups), out)
  final= factor(final)
  final=factor(final, labels=1:length(levels(final)))
  final
}



###############
#integrate names into the dataframe
###############
namestocol=function(df){
  df[, "id"]=rownames(df)
  df
}

###############
#shortcut paste with "_" in between
###############

paste_=function(...) paste(..., sep="_")

###############
#open png or pdf with a timestamp (not close)
###############
timestamp= function(nm) paste0(Sys.Date(),"_", nm)


tpng=function(nm,path=config$plotpath,  ...) png(paste0(path, timestamp(nm), ".png"),...)

tpdf=function(nm, path=config$plotpath,pw=2.5, sca=2.5, width=1, height=1, ...) pdf(paste0(path, timestamp(nm), ".pdf"), wi=width*pw*sca, he=height*pw*sca, ...)

#############################################################
#                   
#############################################################

replacepoint= function(num) paste(strsplit(as.character(num), split='\\.')[[1]], collapse="p")


#######################################################################################
#nearest neighbor operations.
#######################################################################################
#get a particular variable from nearest neighbors of a cell. 
#######################################################################################
getnnvar= function(so, cll, vr) metadata(so)[so@graphs$RNA_nn[cll, ]==1, vr]
getsnnvar= function(so, cll, vr) metadata(so)[so@graphs$RNA_snn[cll, ]==1, vr]

getnn= function(so, cll) metadata(so)[so@graphs$RNA_snn[cll, ]==1, ]



#######################################################################################
#get data slots from a seurat object
#######################################################################################

getcounts<- function(so) so@assays$RNA@counts



#######################################################################################
# Cluster merging operations
########################################################################################
####look for cells that belong to any of the cluster ids in the list and true them. specifically for scM
combineclusters= function(lst, outvar="seurat_clusters") Reduce('|', lapply(lst, function(el) metadata(scM)[,outvar]==el))
### more general to be applied to any seurat object
combineclusters2= function(so, lst, outvar="seurat_clusters") Reduce('|', lapply(lst, function(el) metadata(so)[,outvar]==el))
####################################
#merge clusters in list lst, in factor column
####################################
renamelabel=function(column, ..., renameto){
  grps=list(...)
  lapply(grps, function(lst){
    groupp=combineclusters(lst)
    sapply(1:length(scM@meta.data[,column]), function(x) ifelse(groupp[x],renameto, as.character(scM@meta.data[x, column])), USE.NAMES=F)
  })
}


triagecells.multi2=function(so, groupvar="seurat_clusters", clean.names=T, ...){
  #... are the clusters in the group one after the other, ordered arbitrarily and grouped using c() when they are meant to be merged
  # warning: if some cluster groups overlap in any extent, cell indices will appear in two or more occasions and will be sequentially replaced!
  
  grps=list(...)
  grps2=grps %>% lapply(., function(x) paste0(addprefix(x, "^"), "$"))
  
  original=so %>% metadata %>% pull(get(groupvar))
  final.arr= rep(NA, length(original))
  indlist=lapply(grps2, function(gr) grepl( gr %>% paste0(., collapse="|"), original) %>% which )
  qlist=lapply(grps2, function(gr) gr %>% paste0(., collapse="|"))
  vector.stamp=function( vector, ind){
    if (ind<= length(grps)){ 
      vector[indlist[[ind]]]=ind
      return(vector.stamp( vector, ind+1))
    }else{
      vector[is.na(vector)]= paste_("og", original[is.na(vector)])
      
 
      
      return(vector) 
    }
  }
  
  out=vector.stamp(final.arr, 1)
  
       if(clean.names){
       reflist= 1:length(unique(out)) %>% givenames(., unique(out))
            
       out=reflist[out] %>% unname %>% as.character
        
      }
  out
}


###########################################
#replace each value in an array for another one.
###########################################

generalreplace=function(value, reff, target){
  
  target[which(reff==value)]
  
}

##############################################
#making the row names into a gene column
##############################################

names2col=function(dff, coln="gene") {
  dff[, coln]= rownames(dff) 
  dff
}

##########################################################
#mohamed's function to accomodate the slingshot data into a dataframe
##########################################################

tidy_ss_output <- function(ss_output){
  library(tidyr)
  library(dplyr)
  library(tibble)
  #cell embeddings
  cell_embedd <- ss_output@reducedDim
  #coordinates of the curve
  curve_coord <- lapply(names(ss_output@curves), function(x){
    `[[`(ss_output@curves[[x]],1) %>% 
      as.data.frame() %>% 
      rename_with(.,~paste0(x,"_",.))
  }) %>%
    bind_cols()
  #order of cell over infered curves
  cell_order <- lapply(names(ss_output@curves), function(x){
    `[[`(ss_output@curves[[x]],3) %>% 
      as.data.frame() %>% 
      rename_with(.,~paste0(x,"_lambda"))
  }) %>%
    bind_cols()
  #bind outputs in a single data frame
  cbind(cell_embedd,curve_coord,cell_order) %>% 
    as.data.frame() %>% 
    rownames_to_column("cell_id")
}


#############################################################################################################################################
#      PLOTTING FUNCTIONS
#############################################################################################################################################


#### modify hex colors up or down. 


dehash=function(color) substr(color, 2,7)
color2dec=function(color) hex2dec(substr(color, 2,7))
dec2color=function(dec) paste0('#', dec2hex(dec))
color2=function(dec) paste0('#', dec2hex(dec))

colormore=function(color, pts){
  
  newdec=color2dec(color)+pts
  
  if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
  if(newdec<hex2dec("000000")){newdec="000000"}
  
  return(dec2color(newdec))
}

colorless=function(color, pts){ #virtually the same as above but automatically subtracts the number so you don't have to worry about signs. 
  
  newdec=color2dec(color)-pts
  
  if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
  if(newdec<hex2dec("000000")){newdec="000000"}
  return(dec2color(newdec))
}


colormix=function(color1, color2){ #mixing two colors!
  getsubcolors= function(co) lapply(c(1,3,5), function(x) substr(dehash(co), x, x+1))
  
  
  newdec=color2dec(color)+colordec(color2)
  
  if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
  if(newdec<hex2dec("000000")){newdec="000000"}
  return(dec2color(newdec))
}

colorsubtract=function(color1, color2){ #subtracting two colors!
  
  newdec=abs(color2dec(color)-colordec(color2))
  
  if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
  if(newdec<hex2dec("000000")){newdec="000000"}
  return(dec2color(newdec))
}



colorbleach=function(color, pts){ ###ADDS PTS TO EACH OF THE r g b SECTIONS
  
  
  newdec=color2dec(color)+pts
  
  if(newdec>hex2dec("FFFFFF")){newdec="FFFFFF"}
  if(newdec<hex2dec("000000")){newdec="000000"}
  
  return(dec2color(newdec))
}

#########################################################
#ACCESSORY FUNCTIONS
#########################################################

metadata= function(x) x@meta.data


seuratprocessing=function(o){
  o=NormalizeData(o)
  o=FindVariableFeatures(o, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  o=ScaleData(o)
  o=RunPCA(o, dims=cfig$PCA_DIMS)
  o=ProjectDim(o)
  #o=RunUMAP(o, reduction="pca", dims=cfg$PCA_DIMS)
  o=FindNeighbors(o)
  o=FindClusters(o)
}

SCTprocessing=function(o){
  o=NormalizeData(o)
  o=FindVariableFeatures(o, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  o=ScaleData(o)
  o=RunPCA(o, dims=cfig$PCA_DIMS)
  o=ProjectDim(o)
  #o=RunUMAP(o, reduction="pca", dims=cfg$PCA_DIMS)
  o=FindNeighbors(o)
  o=FindClusters(o)
  o
}

seuratprocessing_umap=function(o){
  o=NormalizeData(o)
  o=FindVariableFeatures(o, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  o=ScaleData(o)
  o=RunPCA(o, dims=cfig$PCA_DIMS)
  o=ProjectDim(o)
  o=RunUMAP(o, reduction="pca", dims=cfg$PCA_DIMS)
  o=FindNeighbors(o)
  o=FindClusters(o)
  o
}
getpname=function(x) unique(x@meta.data[, "orig.ident"])


see.clusterings=function(so){
  grep("seurat_clusters", colnames(so@meta.data), value=T)
}

set.clustering=function(so, num=NULL){
all.clusterings=see.clustrerings(so)
if(!is.null(num)){
so@meta.data$seurat_clusters_previous=so@meta.data$seurat_clusters
so@meta.data= so@meta.data %>% dplyr::mutate(seurat_clusters= !!sym(all.clusterings[num]), current_clustering=all.clusterings[num])
so
}else{
fcat("no clustering performed.")  
}
}


reset.clustering=function(so,  original=F){
all.clusterings=see.clustrerings(so)

if("seurat_clusters_previous" %in% all.clusterings){
so@meta.data$seurat_clusters=so@meta.data$seurat_clusters_previous
so@meta.data$current_clustering=NULL
}else{
  fcat("no previous clustering explicitly registered")
  if("seurat_clusters_original" %in% all.clusterings && original==T){
    fcat("resetting to original clustering...")
so@meta.data$seurat_clusters=so@meta.data$seurat_clusters_original
}else{
  fcat("No resetting was performed. You can try setting original=T")
  fcat("possible changes:")
  fcat(paste(all.clusterings, collapse=","))
}}
  so

}


seuratprocessing.sct=function(so=NULL,rcache=NULL, project=NULL, verbose=T, umapdims=1:10, vfeatures=5000, label="", reduction.label="", reduction.for.neighbors=NULL, resolution.clusters=0.8,refresh=T, reprocess=F, recluster=T, harmony.vars=NULL, vars.to.regress.sct=NULL, nn.cluster=20, nn.umap=30L, min.dist=0.3, repulsion.strength=1, cluster.algorithm=1, cluster.method="igraph"){
  library(reticulate)
#use_condaenv("py311", required = TRUE)

 if(label !=""){
  label=make.names(label)
  reduction.label=make.names(reduction.label)
 }
 if(!is.null(rcache)){
  simpleCache(rcache, assignToVar="so", reload=T) 
   so@meta.data$umapdims=paste(umapdims, collapse=";")
 }
  
  
 
   
  sctname=paste0("SCT", label)
  pcaname=paste0("pca", label)
  umapname=paste0("umap", reduction.label)
  clustername=paste0("seurat_clusters_", label, reduction.label)
  harmony.reduction=paste_("harmony", pcaname)
  harmony.umap=paste_(umapname, "harmony") 
 graphname=ifelse(label!="", paste0("SCT_", label), "SCT")
  
fcat("Defining names for this run:",
    "SCT assay: ", sctname, "\n",
     "PCA: ", pcaname, 
     "UMAP", umapname,
     "cluster variable", clustername, 
     "if harmony:\n",
     "reduction: ", harmony.reduction,"\n", 
     "harmony umap", harmony.umap, "\n",
     "graph:\n", graphname, "\n")

 if(label==""){
   label= "none"
 }
 
if(!is.null(so) & is.null(project)){
 if(Project(so)=="" ||is.null(Project(so)) ){
   Project(so)= paste0("so_", analysis.version, "_", substr(digest::digest(so)), 1,3)
 }


project=Project(so)  
  
}


# execution
  
  #simpleCache(paste0("processed_dataset_id",label,"_", Project(so), "_vfeatures_", vfeatures) %>% addversion, {1}, assignToVar="m", reload=T)
  
  #if(m!=1){ #if it is not one, the label is not free and must be changed
  #  fcat("Aborting: cache with label already exists. Stopping execution to avoid overwriting")
  #  return(m)
  #}else{
   
    fcat("loading so cache")
 sctcachename=paste0("processed_dataset_id",label,"_", project, "_vfeatures_", vfeatures) %>% addversion   
  
cache.exists=file.exists(file.path(config$cachedir, paste0(sctcachename, ".RData")))
#if reprocessing is requested while cache exists,or
# if the assay is not present or cache does not exist
  if((reprocess && cache.exists) || ( (!is.null(so) && (!any(grepl(sctname, names(so@assays))))) && !cache.exists)    ){
  
    rl=F
  }else{
    rl=T}
rc=!rl


fcat("cache", sctcachename, ifelse(cache.exists, "exists", "doesn't exist"), "recreate is" , rc, "reload is", rl)
if(refresh || is.null(so)   ){
    simpleCache(sctcachename, {
      
  fcat("running SCTransform...")
  so=SCTransform(so, verbose=verbose, variable.features.n=vfeatures, new.assay.name=sctname, vars.to.regress=vars.to.regress.sct)
  fcat("Calculating PCA...")
  so=RunPCA(so, assay=sctname, verbose=verbose, reduction.name=pcaname )
  
so
  }, assignToVar="so", reload=rl, recreate=rc)
}
    
    simpleCache(paste0("umap_dataset_id",label,"_", Project(so), "umapname", umapname), {
   fcat("Calculating UMAP...")
  so=RunUMAP(so, reduction=pcaname, dims=umapdims, n.neighbors=nn.umap, return.model=T, verbose=verbose, reduction.name=umapname, min.dist=min.dist, repulsion.strength=repulsion.strength )
  so@reductions[[umapname]]  
  },assignToVar="umapreduction", reload=rl, recreate=rc)
 
  so@reductions[[umapname]]=umapreduction 
  
      if(!is.null(harmony.vars)){
      
        simpleCache(paste0("umap_harmony_dataset_id",label,"_", Project(so), "reduction_name", harmony.reduction), {    
        
          
    so=RunHarmony(so, verbose=verbose, group.by.vars=harmony.vars, reduction=pcaname, reduction.save=harmony.reduction, assay.use=sctname)
    so=RunUMAP(so, reduction=harmony.reduction, dims=umapdims, return.model=T, verbose=verbose, reduction.name=harmony.umap)
    list(harmony=so@reductions[[harmony.reduction]], umapharmony=so@reductions[[harmony.umap]])  
        },assignToVar="harmonyreductions", reload=rl, recreate=rc)
        
        so@reductions[[harmony.reduction]]=harmonyreductions$harmony
        so@reductions[[harmony.umap]]=harmonyreductions$umapharmony
    
    fcat("model is only being returned for the Harmony umap. for a model on a different reduction please use seuratprocessing.sct")
    if(is.null(reduction.for.neighbors)){
      fcat("defaulting neighbors reduction to harmony PCA with id", harmony.reduction)
      reduction.for.neighbors=harmony.reduction
    }
    
   
       
  graphname=paste_(graphname, "harmony")
  
  
  
  }else{
      
    
          if(is.null(reduction.for.neighbors)){
      fcat("defaulting neighbors reduction to", pcaname)
      reduction.for.neighbors=pcaname
    }
  }
  
  
  

  #additional checks
  if(!reprocess){
    #if the assay already exists
       if(any(grepl(sctname, names(so@assays)))){
    
     fcat("reprocessing aborted: assay", sctname, "already present while reprocess=FALSE. Using existing assays for nn and clustering")
     
     
   }
    
  
       if(any(grepl(paste(pcaname, umapname, sep="|"), names(so@reductions)))){
    
     fcat("reprocessing aborted: reductions with label",grepl(sctname, names(so@reductions))  ,"already present while reprocess=FALSE. Using existing assays for nn and clustering")
     
       }
    
  }
  

    ## finish constructing the graph name(s)
    graphname=paste0(paste_(graphname,reduction.for.neighbors), c("nn","snn"))  
  #calculate how many dimensions should be used for neeighbors depending on whether the reduction is umap or not
    
    if(recluster || !(graphname %in% names(so@graphs))){
    
  if(grepl("umap", reduction.for.neighbors)){neighdims=c(1,2)}else{neighdims=umapdims}
  
  so=FindNeighbors(so, verbose=verbose, k.param=nn.cluster, dims=neighdims, reduction=reduction.for.neighbors, graph.name=graphname)
  
  if("seurat_clusters" %in% colnames(so@meta.data)){
   
   if("current_clustering" %in% colnames(so@meta.data)){
   so@meta.data=so@meta.data %>% dplyr::mutate(previous_clustering=current_clustering) #this is a string
   so@meta.data=so@meta.data %>% dplyr::mutate(seurat_clusters_previous=seurat_clusters) 
   
   }else{
     so@meta.data=so@meta.data %>% dplyr::mutate(seurat_clusters_original=seurat_clusters) 
     so@meta.data=so@meta.data %>% dplyr::mutate(previous_clustering="seurat_clusters_original")
   }
  }
  
  so=FindClusters(so, verbose=verbose, resolution=resolution.clusters, graph.name = graphname[2],algorithm=cluster.algorithm, method=cluster.method )
  so@meta.data=so@meta.data %>% dplyr::mutate(!!sym(paste0("seurat_clusters_on_", reduction.for.neighbors)):=seurat_clusters, current_clustering=paste0("seurat_clusters_on_", reduction.for.neighbors))
    }#recluster
 # so

#fcat("Processing finished. SCT stored under ", sctname, "\nPCA stored under ", pcaname, "\numap stored under ", umapname, "\ngraph stored under ", graphname[1], "\n", )

#if(!is.null(harmony.vars)){
#fcat("Harmony stored under ", harmony.reduction, harmony.umap, "\n")  
 
 so 
}


seuratprocessing.sct.harmony=function(so,rcache=NULL, verbose=T, umapdims=1:10, vfeatures=2000, label="", reduction.for.neighbors=NULL, group.by.vars=NULL, resolution.clusters=0.8){
  
    
 if(!is.null(rcache)){
  simpleCache(rcache, assignToVar="so", reload=T) 
 }
  
  so@meta.data$umapdims=paste(umapdims, collapse=";")
 simpleCache(paste0("processed_dataset_id",label, "_", Project(so), "_",digest::digest(so), "_vfeatures_", vfeatures) %>% addversion, {
   
  sctname=paste0("SCT", label)
  pcaname=paste0("pca", label)
  umapname=paste0("umap", label)
  graphname=paste0(ifelse(label!="", paste0("SCT_", label), "SCT"), c("_nn","_snn"))
  clustername=paste0("seurat_clusters_", label)
  
   so=SCTransform(so, verbose=verbose, variable.features.n=vfeatures, new.assay.name=sctname)

  so=RunPCA(so, assay=sctname, verbose=verbose, reduction.name=pcaname )
  
  if(!is.null(group.by.vars)){
    so=RunHarmony(so, verbose=verbose, group.by.vars=group.by.vars, reduction=pcaname, reduction.save=paste_("harmony", pcaname), assay.use=sctname )
    so=RunUMAP(so, reduction=paste_("harmony", pcaname), dims=umapdims, return.model=T, verbose=verbose, reduction.name=paste_(umapname, "harmony") )
    fcat("model is only being returned for the Harmony umap. for a model on a different reduction please use seuratprocessing.sct")
    if(is.null(reduction.for.neighbors)){
      fcat("defaulting neighbors reduction to harmony PCA with id", paste_("harmony", pcaname))
      reduction.for.neighbors=paste_("harmony", pcaname)
    }
    
  graphname=paste_(graphname, "harmony")
  }
  
  so=RunUMAP(so, reduction=pcaname, dims=umapdims, return.model=F, verbose=verbose, reduction.name=umapname )
  
  #calculate how many dimensions should be used for neeighbors depending on whether the reduction is umap or not
  if(grepl("umap", reduction.for.neighbors)){neighdims=c(1,2)}else{neighdims=umapdims}
  
  so=FindNeighbors(so, verbose=verbose,  dims=neighdims, reduction=reduction.for.neighbors, graph.name=graphname)
  
  if(any(grepl("seurat_clusters", colnames(so@meta.data)))){
    
    num=grepl("seurat_clusters", colnames(so@meta.data)) %>% sum
    
   so@meta.data=so@meta.data %>% dplyr::mutate(!!sym(paste_("seurat_clusters_original", num+1)):=seurat_clusters) 
  }
  
  so=FindClusters(so, verbose=verbose, resolution=resolution.clusters, graph.name = graphname[2])
  so
}, assignToVar="so", recreate=T)
 so 
}





seuratprocessing_markers1=function(so, verbose=F, umapdims=30, mincpm=3, minfc=1, fdr=0.05, selecttop=5, vfeatures=2000){
  so=SCTransform(so, verbose=verbose, variable.features.n=vfeatures)
  so=RunPCA(so, verbose=verbose)
  so=FindNeighbors(so, verbose=verbose)
  so=FindClusters(so, verbose=verbose)
  #so=RunUMAP(so, dims=1:di, verbose=verbose, return.model=T)
  
  
  cat(paste0("finding markers using the Riffle permutation test...\n"))  
  
  m=cb_pos_markers(counts=so@assays$RNA@counts, grouping=so %>% metadata %>% pull("seurat_clusters") )
  
  clusters=m %>% pull(group) %>% unique %>% as.numeric %>% sort
  ## arrange and filter best markers
  
  if(is.null(selecttop)){
    lt=lapply(clusters, function(x) m %>% filter(group==x,FDR<=fdr, logFC>=minfc, logCPM>=mincpm) %>% arrange(-logFC))
    topstr=""  
  }else{
    lt=lapply(clusters, function(x) m %>% filter(group==x,FDR<=fdr, logFC>=minfc, logCPM>=mincpm) %>% arrange(-logFC) %>% head(n=selecttop))
  }
  mp=lt %>% Reduce(rbind,.)
  
  cat("Number of filtered markers: ", nrow(mp), "(", (nrow(mp)/nrow(m))*100, "% markers post filtering)\n")
  so$RNA@misc$markers <- m
  so$RNA@misc$top_markers <- mp
  cat("Dataset", nm, ": calculating residuals for missing genes in scale data...\n")
  so <- GetResidual(so, features = so$RNA@misc$top_markers$feature, verbose = F)
  
  gc()
  so
}


seuratprocessing.markers2=function(so, verbose=T, method="deseq", group_column="seurat_clusters", replicate_column=NULL, padj=0.05, minfc=0.5,selecttop=NULL,minrate=0.5, fullreload=T, fullrecreate=F, getresidual=F, min.cells.per.group=3, vfeatures=2000){
  
  
  so=SCTransform(so, verbose=verbose, variable.features.n=vfeatures)
  so=RunPCA(so, verbose=verbose)
  so=FindNeighbors(so, verbose=verbose)
  so=FindClusters(so, verbose=verbose)
  
  fcat(paste0("finding markers using the DElegate package...\n") ) 
  #simpleCache(paste_("markers_DElegate_seuratproject",Project(so), "on_groups_of", group_column, "splitbyreps", replicate_column, "method", method), {
  
  #function (object, meta_data = NULL, group_column = NULL, replicate_column = NULL, 
  #method = "edger", min_rate = 0.05, min_fc = 1, lfc_shrinkage = NULL, 
  #verbosity = 1)   
  if(min.cells.per.group){
    
    countstab=so %>% metadata %>% group_by(!!sym(group_column)) %>% summarise(counts=n()) %>% as.data.frame
    
    goodcats=countstab[ countstab$counts>=min.cells.per.group, group_column]
    fcat("enough cells found for", paste(goodcats, collapse=","))
    
    so = so[, (metadata(so) %>% filter(!!sym(group_column) %in% goodcats) %>% rownames)]
    
  }
  
  
  m= DElegate::FindAllMarkers2(object=so, group_column=group_column, replicate_column=replicate_column, method=method, min_rate=minrate)
  m
  #m=cb_pos_markers(scwt@assays$RNA@counts, grouping=sma %>% metadata %>% pull() )
  #}, assignToVar="m", reload=fullreload, recreate=fullrecreate)
  
  #print(m)
  
  if(group_column %in% (allcolors %>% names)){
    clusters=allcolors[[group_column]][allcolors[[group_column]] %>% names %in% (so %>% metadata %>% pull(group_column) %>% unique)] %>% names
  }else{
    clusters=so %>% metadata %>% pull(group_column) %>% unique 
    
  }
  
  ## arrange and filter best markers.
  ## currently filtering such that the minimum rate for in group cells is at least minrate. 
  pa=padj
  mi=minfc
  if(is.null(selecttop)){
    lt=lapply(clusters, function(x) m %>% filter(group1==x,padj<=pa, log_fc>=mi, rate1>=minrate) %>% arrange(-log_fc))
    topstr=""  
  }else{
    lt=lapply(clusters, function(x) m %>% filter(group1==x,padj<=pa, log_fc>=mi, rate1>=minrate) %>% arrange(-log_fc) %>% head(n=selecttop))
  }
  mp=lt %>% Reduce(rbind,.)
  
  cat("Number of filtered markers: ", nrow(mp), "(", (nrow(mp)/nrow(m))*100, "% markers post filtering)\n")
  so[["RNA"]]@misc$markers <- m
  so[["RNA"]]@misc$top_markers <- mp
  if(getresidual){
    so <- GetResidual(so, features = so$RNA@misc$top_markers$feature, verbose = F)
  }
  so
}






seuratmarkers.delegate=function(so, rcache=NULL, method="deseq", group_column="seurat_clusters", slot=group_column, replicate_column=NULL, padj=0.05, minfc=0.5,selecttop=NULL,minrate=0.5, getresidual=T, min.cells.per.group=10, label="", verbosity=1, reprocess=F){
  fcat(paste0("finding markers using the DElegate package...\n")) 
 if(!is.null(rcache)){
  simpleCache(rcache, assignToVar="so", reload=T) 
 }
  
 misc.object=so$RNA@misc 
 so$RNA@misc=list()
 
  so.digest=digest::digest(so)
  
  fulldigest=digest::digest(package.vars(so.digest, method, replicate_column,group_column,padj, minfc,selecttop, minrate,getresidual, min.cells.per.group))
  
  markercache=paste_("markers_DElegate_so_id_",Project(so), label, 
                #so.digest,
                #"_analysis_", 
                "groupby", group_column, "reps", replicate_column, "minfc", minfc,"minrate",  minrate, "mincells", min.cells.per.group) %>% addversion
  
  fcat("markercache", markercache,ifelse(check.cache(markercache), "exists", "doesn't exist"))
  
  if(!reprocess){
   rl=T 
  }else{
    rl=F}
  rc=!rl
  
  simpleCache(markercache, {
  
  #function (object, meta_data = NULL, group_column = NULL, replicate_column = NULL, 
  #method = "edger", min_rate = 0.05, min_fc = 1, lfc_shrinkage = NULL, 
  #verbosity = 1)   
  if(min.cells.per.group){
    
    countstab=so %>% metadata %>% group_by(!!sym(group_column)) %>% summarise(counts=n()) %>% as.data.frame
    
    goodcats=countstab[ countstab$counts>=min.cells.per.group, group_column]
    fcat("enough cells found for", paste(goodcats, collapse=","))
    
    so = so[, (metadata(so) %>% filter(!!sym(group_column) %in% goodcats) %>% rownames)]
    
  }
  
  
  m<<- DElegate::FindAllMarkers2(object=so, group_column=group_column, replicate_column=replicate_column, method=method, min_rate=minrate, min_fc=0.1, verbosity=verbosity)
  m
  #m=cb_pos_markers(scwt@assays$RNA@counts, grouping=sma %>% metadata %>% pull() )
  }, assignToVar="m", reload=rl, recreate=rc)
  m<<-m
  #print(m)
  
  ## arrange and filter top markers.
  ## currently filtering such that the minimum rate for in group cells is at least minrate. 
  pa=padj
  mi=minfc
  
  mp=m %>% group_split(group1) %>% lapply(., function(x) x %>% filter(log_fc>=minfc, padj<=0.05, rate1>=minrate) %>% arrange(-log_fc)) %>% Reduce(rbind,. )
  mb=m %>% group_split(group1) %>% lapply(., function(x) x %>% filter(log_fc<=-minfc, padj<=0.05) %>% arrange(log_fc)) %>% Reduce(rbind,. )

  cat("Number of filtered markers: ", nrow(mp), "(", (nrow(mp)/nrow(m))*100, "% markers post filtering)\n")
  so[["RNA"]]@misc=misc.object
  so[["RNA"]]@misc[[slot]]$markers <- m
  so[["RNA"]]@misc[[slot]]$top_markers <- mp
  so[["RNA"]]@misc[[slot]]$negative_markers <- mb
  if(getresidual){
    fcat("getting residual")
    so <- GetResidual(so, features = so[["RNA"]]@misc[[slot]]$markers, verbose = F)
  }
  so
}


seuratmarkers.roc=function(so, rcache=NULL, group_column="seurat_clusters", slot=group_column, replicate_column=NULL, padj=0.05, minfc=0.5,selecttop=NULL,minrate=0.5, getresidual=F, min.cells.per.group=10, label="", column.names.delegate=T){
  
  method="roc"
  fcat(paste0("finding markers using the Seurat package (roc method)...\n")) 
 if(!is.null(rcache)){
  simpleCache(rcache, assignToVar="so", reload=T) 
 }
 
 misc.object=so$RNA@misc 
 so$RNA@misc=list()
  so.digest=digest::digest(so)
  
  fulldigest=digest::digest(package.vars(so.digest, method, replicate_column, group_column,padj, minfc,selecttop, minrate,getresidual, min.cells.per.group))
  
  simpleCache(paste0("markers_seuratroc_so_id_", label, "_",fulldigest, "_groupby_", group_column, "_method_", method) %>% addversion, {
  
  #function (object, meta_data = NULL, group_column = NULL, replicate_column = NULL, 
  #method = "edger", min_rate = 0.05, min_fc = 1, lfc_shrinkage = NULL, 
  #verbosity = 1)   
  if(min.cells.per.group){
    
    countstab=so %>% metadata %>% group_by(!!sym(group_column)) %>% summarise(counts=n()) %>% as.data.frame
    
    goodcats=countstab[ countstab$counts>=min.cells.per.group, group_column]
    fcat("enough cells found for", paste(goodcats, collapse=","))
    
    so = so[, (metadata(so) %>% filter(!!sym(group_column) %in% goodcats) %>% rownames)]
    
  }
  
  Idents(so)=group_column
  m<<- Seurat::FindAllMarkers(object=so, test.use="roc")
  m 

  #>                 p_val avg_log2FC pct.1 pct.2    p_val_adj cluster     gene
  }, assignToVar="m", reload=T)
  m<<-m


 ###adding variables from the delegate tables
   m= m %>% dplyr::mutate( log_fc=avg_log2FC, rate1=pct.1, rate2=pct.2, group1=cluster,  feature=gene) 
    

  m=m %>% dplyr::mutate(pct.ratio=rate1/rate2) %>% arrange(group1, -pct.ratio)
  mp=m %>% group_split(group1) %>% lapply(., function(x) x %>% filter(log_fc>=minfc, padj<=0.05) %>% arrange(-pct.ratio)) %>% Reduce(rbind,. ) %>% as.data.frame
  mb=m %>% group_split(group1) %>% lapply(., function(x) x %>% filter(log_fc<=-minfc, padj<=0.05) %>% arrange(pct.ratio)) %>% Reduce(rbind,. ) %>% as.data.frame

  
  
  cat("Number of filtered markers: ", nrow(mp), "(", (nrow(mp)/nrow(m))*100, "% markers post filtering)\n")
  so[["RNA"]]@misc=misc.object
  so[["RNA"]]@misc[[slot]]$markers_roc <- m
  so[["RNA"]]@misc[[slot]]$top_markers_roc <- mp
  so[["RNA"]]@misc[[slot]]$negative_markers_roc <- mb
  if(getresidual){
    so <- GetResidual(so, features = so$RNA@misc[[slot]]$top_markers$feature, verbose = F)
  }
  so
}




################################################################################
# calculate markers using DE, each versus each other
################################################################################


seuratmarkers.delegateDE=function(so, method="deseq", group_column="seurat_clusters", compare="all_vs_all", replicate_column=NULL, padj=0.05, minfc=0.5,selecttop=NULL,minrate=0.5, fullreload=T, fullrecreate=F, getresidual=F, min.cells.per.group=3, lfc_shrinkage=NULL){
  #the other option for "compare" is each vs rest"
  fcat(paste0("finding markers using the DElegate package...\n") ) 
  simpleCache(paste_("markers_DElegate_seuratproject",Project(so), "on_groups_of", group_column, "splitbyreps", replicate_column, "method", method), {
    
    #function (object, meta_data = NULL, group_column = NULL, replicate_column = NULL, 
    #method = "edger", min_rate = 0.05, min_fc = 1, lfc_shrinkage = NULL, 
    #verbosity = 1)   
    if(min.cells.per.group){
      
      countstab=so %>% metadata %>% group_by(!!sym(group_column)) %>% summarise(counts=n()) %>% as.data.frame
      
      goodcats=countstab[ countstab$counts>=min.cells.per.group, group_column]
      fcat("enough cells found for", paste(goodcats, collapse=","))
      
      so = so[, (metadata(so) %>% filter(!!sym(group_column) %in% goodcats) %>% rownames)]
      
    }
    
    
    m= DElegate::findDE(object=so, group_column=group_column, compare=compare, replicate_column=replicate_column, method=method,lfc_shrinkage=NULL)
    m
    #m=cb_pos_markers(scwt@assays$RNA@counts, grouping=sma %>% metadata %>% pull() )
  }, assignToVar="m", reload=fullreload, recreate=fullrecreate)
  
  #print(m)
  
  if(group_column %in% (allcolors %>% names)){
    clusters=allcolors[[group_column]][allcolors[[group_column]] %>% names %in% (so %>% metadata %>% pull(group_column) %>% unique)] %>% names
  }else{
    clusters=so %>% metadata %>% pull(group_column) %>% unique 
    
  }
  
  ## arrange and filter best markers.
  ## currently filtering such that the minimum rate for in group cells is at least minrate. 
  pa=padj
  mi=minfc
  if(is.null(selecttop)){
    lt=lapply(clusters, function(x) m %>% filter(group1==x,padj<=pa, log_fc>=mi, rate1>=minrate) %>% arrange(-log_fc))
    topstr=""  
  }else{
    lt=lapply(clusters, function(x) m %>% filter(group1==x,padj<=pa, log_fc>=mi, rate1>=minrate) %>% arrange(-log_fc) %>% head(n=selecttop))
  }
  mp=lt %>% Reduce(rbind,.)
  
  cat("Number of filtered markers: ", nrow(mp), "(", (nrow(mp)/nrow(m))*100, "% markers post filtering)\n")
  so[["RNA"]]@misc$markers <- m
  so[["RNA"]]@misc$top_markers <- mp
  if(getresidual){
    so <- GetResidual(so, features = so$RNA@misc$top_markers$feature, verbose = F)
  }
  so
}



####### get only induction from barcode


recreateAll <- F

#dA <- fread(metapath)
#dA <- merge(dA[bsf_name!="",], fread(barcodepath)[,.(multiseq_id, multiseq_sequence)], by="multiseq_id", all.x=T)
#setkey(dA, sample_name)

#################make roc curve
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}


############Quality Ccontrol QC plots function

add.mito= function(so){ 
  # for old seurat
#percent.mt<- PercentageFeatureSet(so, pattern = "^MT-") %>%givecolnames(., nms="percent.mt")  
# for new seurat  
so<- PercentageFeatureSet(so, pattern = "^MT-", col.name='percent.mt')# %>% coerce.vector.to.matrix %>% givecolnames(., nms="percent.mt")
so
}  

logvarname= function(x) paste0("log10_", x)



qcplots=function(so, 
                 plotred=FALSE, 
                 plotsize=10, 
                 plotname="qcplots",
                 plot.title=NULL,
                 demux= NULL, 
                 topx=10000,
                 topf=20000, 
                 ww=1, 
                 hh=1,
                 colorby="seurat_clusters",
                 format="png",
                 nott=scales::comma,
                 path="./", 
                 cutoffs=NULL, 
                 cellnumbers=NULL, 
                 seriate.method=NULL,
                 seriate.pars=list(ncells=500,nmarkers=NULL),
                 plot.empty=NULL, 
                 groupvar="group1",
                 lfcvar="log_fc", 
                 nmarkers.heatmap=20, 
                 sz=0.1, 
                 colorlist=NULL){
  library(grid)
fcat("qcplots_checkpoint_1")
so <- add.mito(so)
so <- bin.counts(so, bin_size=1000, max_bin=15000, colname="bin.counts")
Idents(so)=colorby  

   markers= so@assays$RNA@misc$seurat_clusters$top_markers
  
   if(is.null(markers)){
    fcat("There are no markers! please add markers under @assays$RNA@misc$top_markers" )
     
     
     
   cells=NULL
     }else{
if(is.null(seriate.method)){            
cells <- cellinfo.init(so, seriate=F)
}else{
cells<<-seriatecells(so, clusvar=colorby, meth=seriate.method, groupvarname=groupvar, lfcvarname = lfcvar, ncells=seriate.pars$ncells, nmarkers=seriate.pars$nmarkers)


}
     }   
  
  
if(is.null(plot.title)){
  plot.title=Project(so)
  }

  
  require(scales)
  tsz=7
  
  
  so$log10_nCount_RNA <- log10(so$nCount_RNA)
  so$log10_nFeature_RNA <- log10(so$nFeature_RNA)
fcat("qcplots_checkpoint_3")
qc1=ggplot(so@meta.data, aes(x=nCount_RNA))+
  geom_histogram(binwidth=20, alpha=0.5)+
  labs(x = "Read depth",y = "Frequency")+
  guides(fill=guide_legend())+
  coord_cartesian(xlim=c(0, topf))+
  scale_y_continuous(labels = nott)+
  scale_x_continuous(labels = nott)+
  theme_classic()+
  labs(title="Reads per barcode")+theme_classic()+theme(axis.text=element_text(size=tsz), axis.title=element_text(size=tsz));
fcat("qcplots_checkpoint_4")
#distribution of counts per barcode
qc2=ggplot(so@meta.data, aes(x=nFeature_RNA))+
  labs(x = "Num. Genes",y = "Frequency")+
  geom_histogram(binwidth=20, alpha=0.5)+
  labs(title="Genes per barcode")+
  scale_y_continuous(labels = nott)+
  scale_x_continuous(labels = nott)+
theme_classic()+
coord_cartesian(xlim=c(0, topx))+theme_classic()+theme(axis.text=element_text(size=tsz), axis.title=element_text(size=tsz));

fcat("qcplots_checkpoint_5")

if(is.null(colorlist)){
 colorlist=list() 
}
if(!(colorby %in% names(colorlist)) || length(colorlist[[colorby]]) != so %>% metadata %>% pull(!!sym(colorby)) %>% unique %>% length){
  library(gtools)
  cats=so %>% metadata %>% pull(!!sym(colorby)) %>% unique %>% mixedsort
  colorlist[[colorby]]=  rainbow(cats %>% length)
  names(colorlist[[colorby]])= cats
}

qc3 <- FeatureScatter(so, group.by=colorby, feature1 = "nFeature_RNA", feature2 = "percent.mt", raster=T, pt.size=sz)+
  scale_color_manual(values=colorlist[[colorby]])+
  labs(x = "Num. Genes",y = "% Mitochondrial")+
  ggtitle("")+theme_classic()+
  scale_x_continuous(labels = scientific)+
  theme(axis.text.x=element_text(size=tsz), axis.text.y=element_text(size=tsz),axis.title=element_text(size=tsz))+NoLegend()

fcat("qcplots_checkpoint_6")
qc4 <- FeatureScatter(so, group.by=colorby,feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=T, pt.size=sz)+
  scale_color_manual(values=colorlist[[colorby]])+
  labs(x = "Read depth",y = "Num.Genes")+
  ggtitle("")+
  scale_y_continuous(labels = nott)+
  scale_x_continuous(labels = nott)+
  theme_classic()+
  theme(axis.text.x=element_text(size=tsz), axis.text.y=element_text(size=tsz), axis.title=element_text(size=tsz),legend.position="top",legend.title=element_blank())

fcat("qcplots_checkpoint_7")
qc5 <- FeatureScatter(so, group.by=colorby,feature1 = "nCount_RNA", feature2 = "percent.mt", raster=T, pt.size=sz)+
  scale_color_manual(values=colorlist[[colorby]])+
  labs(x = "Read depth",y ="% Mitochondrial")+
  ggtitle("")+theme_classic()+
  scale_x_continuous(labels = scientific)+
  theme(axis.text.x=element_text(size=tsz), axis.text.y=element_text(size=tsz), axis.title=element_text(size=tsz))+NoLegend()  


if(plotred!=FALSE){

  dms=10 
     umap1=DimPlot(so, group.by = colorby, reduction=plotred, label=T, label.size=6, label.color="black", cols= colorlist[[colorby]], raster=T, pt.size=sz)+
     ggtitle(paste0(plot.title, ": ", plotred, "\n", ncol(so), " cells"))+NoAxes()+NoLegend()#+NoLegend()#+scale_color_manual( values=qccols);tsne1 

   umapmt=FeaturePlot(so, feature="percent.mt", reduction=plotred, raster=T, pt.size=sz)+
     ggtitle("% Mitochondrial")+NoAxes()#+NoLegend()#+scale_color_manual( values=qccols);tsne1 
   
      umapreads=FeaturePlot(so,  feature="nCount_RNA",reduction=plotred, raster=T, pt.size=sz)+
     ggtitle("Reads/BC")+NoAxes()#+NoLegend()#+scale_color_manual( values=qccols);tsne1 
      
            umapgenes=FeaturePlot(so,feature="nFeature_RNA",  reduction=plotred, raster=T, pt.size=sz)+
     ggtitle("Genes/BC")+NoAxes()#+NoLegend()#+scale_color_manual( values=qccols);tsne1 
            
            fcat("qcplots_checkpoint_8")
            
            vlnreads= ggplot(so %>% metadata, aes(x=factor(!!sym(colorby)), fill=factor(!!sym(colorby)), y=nCount_RNA))+geom_violin()+scale_y_log10(labels = label_number())+scale_fill_manual(values=colorlist[[colorby]])+theme_classic()+NoLegend()
            vlngenes= ggplot(so %>% metadata, aes(x=factor(!!sym(colorby)), fill=factor(!!sym(colorby)), y=nFeature_RNA))+geom_violin()+scale_fill_manual(values=colorlist[[colorby]])+theme_classic()+NoLegend()
            vlnmt= ggplot(so %>% metadata, aes(x=factor(!!sym(colorby)), fill=factor(!!sym(colorby)), y=percent.mt))+scale_fill_manual(values=colorlist[[colorby]])+geom_violin()+theme_classic()+NoLegend()
            
            fcat("qcplots_checkpoint_9")
            if( !is.null(cutoffs)){
              lcolor="red"
             vlnreads= vlnreads+geom_hline(yintercept= cutoffs$counts, color=lcolor)+ggtitle(paste("count cutoff:", cutoffs$counts))
             vlngenes= vlngenes+geom_hline(yintercept= cutoffs$features, color=lcolor)+ggtitle(paste("feature cutoff:", cutoffs$features))
             vlnmt= vlnmt+geom_hline(yintercept= cutoffs$mito, color=lcolor)+ggtitle(paste("% mito cutoff:", cutoffs$mito))
            }
             fcat("qcplots_checkpoint_10")
             if( !is.null(cellnumbers)){
            ncolor="black"
               sumdata= so %>% metadata %>% group_by(seurat_clusters) %>% summarise(counts=n(), nCount_RNA=median(nCount_RNA), nFeature_RNA=median(nFeature_RNA), percent.mt=median(percent.mt))
                  
             vlnreads= vlnreads+geom_text(data=sumdata, aes( label=paste0("n=",counts)), angle=90, col=ncolor)
             vlngenes= vlngenes+geom_text(data=sumdata, aes( label=paste0("n=",counts)), angle=90, col=ncolor)
             vlnmt= vlnmt+geom_text(data=sumdata, aes( label=paste0("n=",counts)), angle=90, col=ncolor)
            }
            fcat("qcplots_checkpoint_11")
      if(!is.null(demux)){ #&& 
        
      demuxplot=ggplot(so %>% metadata, aes(x=seurat_clusters))+geom_bar(aes(fill=!!sym(demux)), position="fill")+theme_classic()
      fcat("qcplots_checkpoint_12")
            
      }else{demuxplot=ggplot()+ggtitle("No demux/ categories to display")}
            

      

   if(!is.null(so$RNA@misc$seurat_clusters$top_markers)){
fcat("Preparing heatmap...")
             
 cat("checkpoint5...\n")                       
ci=cellinfo.init(so)
hm=heatmap.workflow(so, ci, assay=DefaultAssay(so), vars=c("seurat_clusters"), extra.annotation = c("bin.counts"), nmarkers = 1, return.heatmap=T, colorlist=colorlist)
#ht_grob <- grid.grabExpr(draw(hm))
ht_grob <- grid.grabExpr({draw(hm, show_annotation_legend = FALSE)}, wrap.grobs = TRUE)
ph <- wrap_elements(full = ht_grob)

#ph <- DoHeatmap(so, features = adjust.markers(cells, nmarkers.heatmap)$markers, group.colors=colorlist[[colorby]], slot = "scale.data", cells = cells$cells) + NoLegend()
}else{ph=ggplot()+ggtitle("No marker heatmap to display")}

                 layout="AAAAMMMMM
             AAAALLLLL
                    AAAALLLLL
                    AAAALLLLL
                    AAAALLLLL
                    EEBBIIIII
                    EEBBIIIII
                    FFCCJJJJJ
                    FFCCJJJJJ
                    GGDDKKKKK
                    GGDDKKKKK
                    "
    
    qc_assembly=wrap_plots(A = umap1, B = qc3, C = qc4,
                  D = qc5, E = umapmt, F=umapreads, G=umapgenes, I=vlnmt, J=vlnreads, K=vlngenes, L=ph,M=demuxplot,
                   design = layout)
  #}else{

    
  #qc_assembly=((umap1)/(qc3+qc4+qc5))|(umapmt+umapreads+umapgenes)|(vlnmt/(vlnreads/vlngenes))#/(qc6+qc7+qc8)
#}
###################
#plotting with tsne  
###################
tpng(plotname, path=path, width=w*plotsize*ww, height=w*plotsize*hh, res=600)
print(qc_assembly)
dev.off()

}else{

#ppc1=DimPlot(object = subset(so, subset = percent.mt<10), cols= c("grey", "blue"), reduction = "pca");
###################
#plotting without tsne  
###################
  
layout <- "
####NN
AABBNN
AABBNN
CCDDEE
CCDDEE
"
  qc_assembly=wrap_plots(A = qc1, qc2 = qc2, C = qc3, D=qc4+NoLegend(), E=qc5, N=qc4,  design = layout)
  

png.pdf.prop=12
if(format=="png"){
tpng(plotname, path=path, width=w*plotsize*ww/png.pdf.prop, height=w*plotsize*hh/png.pdf.prop, res=600)
print(qc_assembly)
dev.off()
}
if(format=="pdf"){
  tpdf(plotname, path=path, width=pw*plotsize*ww/png.pdf.prop, height=pw*plotsize*hh/png.pdf.prop)
  print(qc_assembly)
  dev.off()
}

}
list(plot=qc_assembly, cells=cells)
}

################################################################################
# make qc analysis for a list of seurat objects
################################################################################


qc.analysis=function(sos, qcpars, demux="hto_demux", label="" ,w=0.7, h=0.7){
qcs=lapply(names(sos), function(x){
  if(is.null(sos[[x]])){return(NULL)}else{
  qc=qcplots(sos[[x]], plotred="umapX", demux =demux, cutoffs = qcpars[[x]],sz=1)$plot
    tpdf(paste0("qcplots_",label, "_",x), width=pw*w, height=pw*h)
    print(qc)
    dev.off()

  }
  qc
})  %>% givenames(., names(sos) )
qcs
}





fcat("0.1")
os <- function() {
  require(dplyr)
  gc();
  message("Objects in MB:")
  objects = ls(envir=.GlobalEnv);
  classes = sapply(objects, function(object) { class(get(object))[1] });
  sizes = sapply(objects, function (object) { object.size(get(object)) } )
  a = data.frame(MB=sizes/1e6, class=classes)
  ord = order(sizes, decreasing=TRUE)
  a2 = a[ord,];
  a2 = a2[! a2$class == "function",]
  print(head(a2, 30))
  message("Sum: ", signif(sum(sizes/1e6),4), " MB (", signif(sum(sizes/1e9),4), " GB)")
  # fcat("total: ", a2 %>% summarise(total=sum(MB)/1000) %>% pull(total) %>% as.numeric , "GB" 
  a2
}
#os()



###########################################################################
# FUNCTION TO KEEP THE BARE BONES MINIMUM SEURAT DATA, REMOVING ALL UNNNECESSARY ASSAYS ETC.
###########################################################################

slimseurat=function(so){
  DietSeurat( so,counts = TRUE,data = TRUE,scale.data = TRUE,features = NULL,
              assays = c("RNA","refAssay"))#,dimreducs = c("pca", "umap"))
}


###data=FALSE curently not supported
slimerseurat=function(so){
  DietSeurat( so,counts = TRUE,scale.data = FALSE,features = NULL,
              assays = "RNA", dimreducs=NULL)#,dimreducs = c("pca", "umap"))
}


############################################################################
#  function to get a plotting table with some genes and some cells together with the metadata
############################################################################
join_meta_exp=function(so, genes, cells=NULL, normrow=FALSE, lognormrow=FALSE, metacols=NULL, assay="RNA"){
  
  if(is.null(cells)){cells= colnames(so)}
  
  m=as.data.frame(t(as.matrix(so[genes,cells ]@assays[[assay]]@data)))
  if(normrow){
    m=apply(m, 1, normvec)
  }else{
    if(lognormrow){
      m=apply(m, 1, lognormvec)
    }
    
  }
  
  if (is.null(metacols)){metacols=colnames(metadata(so))}  
  cbind( metadata(so)[cells, metacols], m)
  
}

#############################################################################
#function to attach metadata and gene reads.
#############################################################################


genesdf=function(so, genes) as.matrix((so@assays$RNA@data))[genes,  ]    


####### function to change the legend size inside guide
legendsize=function(x) guide_legend(override.aes = list(size = x) )

##############################################################################
#resizing overall text in the plot
resizetext=function(x) theme_classic(base_size = x)

rotatex= function(x) theme(axis.text.x = element_text(angle = x))


################################################################################
#Heatmap accessory functions
################################################################################

###function to find the names of the genes in each cluster
genesinclusters=function(phe, hvarmat,  top=NULL){
  
  lapply(row_order(phe), function(x) {
    
    if(!is.null(top)){
      rownames(hvarmat)[x][1:top]
    }else{
      rownames(hvarmat)[x]
    }
  }
  )
}



topnclustmembers=function(phe, top=3){
  
  Reduce(c, lapply(row_order(phe), function(x) {
    x[1:top]
    
  }
  ))
}
###############################################################

getheatmapcolors=function(ph, var, orientation='top_annotation'){
  base::eval(glue::glue("ph@{orientation}@anno_list${var}@color_mapping@colors"))
  #ph@{orientation}@anno_list${var}@color_mapping@colors
  }

#####################################################################

namenums=function(x, prefix=NULL, zero=F){
  if(is.null(prefix)){
    
    if(zero){n=1}else{n=0}
    names(x)= (1:length(x))-n; return(x)}else{
      names(x)= paste0(prefix, 1:length(x))
    }
  
  x
}

########################################################################

showpalette=function(p, pname="generic"){
  sca=10
  rr=600
  w=250
  tpng(paste_("palette", pname), res=rr, he=sca*w, wi=sca*w*2)
  
  pl=ggplot(1:length(p) %>% as.data.frame %>% dplyr::mutate(h=1, colors=p, name=names(p)), aes(x=factor(.), y=factor(h), label=colors))+
    geom_col(aes(fill=factor(.)))+
    scale_fill_manual(values=namenums(p))+
    theme_classic()+
    coord_cartesian(ylim=c(0,1))+
    geom_text(y=0.5, color="black", angle=90)+
    #NoAxes()+
    NoLegend()+
    ggtitle(pname)
  
  if(is.null(names(p))){
    print(pl )
  }else{
    print(pl+geom_text(inherit.aes=T,aes(label=name), y=1.1, color="black", angle=90))  
  }
  dev.off()
  
  return(pl+geom_text(inherit.aes=T,aes(label=name), y=1.1, color="black", angle=90))
}
#showpalette(p)

##################################################################################
#GGUMAP wrapper  and derivatives
##################################################################################

adjust.seurat.colors=function(so){
 if("seurat_clusters" %in% colnames(so@meta.data)){
   num=length(so@meta.data$seurat_clusters %>%  unique)
   #spectral
   #return(colorRampPalette(RColorBrewer::brewer.pal(20, "Spectral"), space="rgb")(num) %>% givename(., 0:(length(num)-1)))
   #rainbow
   return(rainbow(num) %>% givename(., 0:(num-1)))
 }else{
   
   return(colorRampPalette(RColorBrewer::brewer.pal(20, "Spectral"), space="rgb")(25))
 }
  
}




ggumap=function(so, umap.df=NULL, colorby="seurat_clusters", labelby="seurat_clusters", sz=0.02, ssca=300, colorlist=NULL, reduction=DefaultDimReduc(so), legend.size=3, color.labels.by="defaultcolor", glasswork.by=NULL, glasswork.params=NULL, glasswork.df=NULL){
  
  require(ggrepel)
  require(ggrastr)
  
reduction.keys=Embeddings(so, reduction) %>% colnames
dim1=reduction.keys[1]
dim2=reduction.keys[2]

  #dim1=paste0(reduction.key, 1)
  #dim2=paste0(reduction.key, 2)
  
  ##############################################################################
  # sourcing data
  ##############################################################################
  fcat(2)
  if(!is.null(umap.df)){
    cats=umap.df %>% pull(!!sym(colorby)) %>% unique
    dat= umap.df %>% dplyr::mutate(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% dplyr::mutate(none="")
  }else{
    cats=so %>% metadata %>% pull(!!sym(colorby)) %>% unique
    dat=umapdf(so, reductions=reduction) %>% dplyr::mutate(UMAP_1=!!sym(dim1), UMAP_2=!!sym(dim2)) %>% dplyr::mutate(none="")
  }
  ##############################################################################
  # preparing colorlist
  ##############################################################################
  fcat(3)
  # colorlist doesn't exist.
  if(is.null(colorlist)){
  colorlist<-list()  
  }
  
  if (is.null(colorlist[[colorby]])){
    colorlist[[colorby]]= rainbow(length(cats)) %>% givename(., gtools::mixedsort(cats))
    fcat("colors used:")
    
    colorlist[[colorby]] %>% dput
  }else{
  
  ### actual categories are smaller than total number of colors
  if (length(colorlist[[colorby]] %>% names) > length(cats) && any(cats %in% colorlist[[colorby]] %>% names) ){
    fcat("adjusting colorlist to present labels")
    
    colorlist[[colorby]]= allcolors[[colorby]][cats]
  }
  ### there are more categories than colors present
  if (length(colorlist[[colorby]] %>% names) < length(cats) && any(cats %in% colorlist[[colorby]] %>% names) ){
    fcat("adjusting colorlist to present labels")
    missinglabels= setdiff(cats, colorlist[[colorby]] %>% names)
    
    colorlist[[colorby]]= c(allcolors[[colorby]][cats], randomcolors(missinglabels) %>% givename(., missinglabels))
  }
  colorlist[["defaultcolor"]]=c(one="black")
  colorlist[[color.labels.by]]=c(colorlist[[color.labels.by]] )
  colorlist[[colorby]]=c(colorlist[[colorby]] )
  if(colorby=="seurat_clusters"){
  colorlist[["seurat_clusters"]]=adjust.seurat.colors(so)
  }
  }
  ##############################################################################
  # constructing plot
  ##############################################################################
  
  fcat(4)
  
  dat2=dat %>% group_by(!!sym(labelby)) %>% summarise(UMAP_1=jitter(mean(UMAP_1)), UMAP_2=jitter(mean(UMAP_2)), defaultcolor="one")
  #catch some pesky history dependent bug with geom_point_rast
  
  gpuc<- tryCatch(
  {
    fcat(5)
    # First try with geom_point_rast
   print(ggplot(dat, aes(x=UMAP_1, y=UMAP_2, label=!!sym(labelby)))+
    geom_point(size=sz, aes(color=!!sym(colorby)))+
    scale_color_manual(values=colorlist[[colorby]], guide = guide_legend(override.aes = list(size = legend.size, shape=15) ))+
    theme_classic()+
    geom_text(data= dat2,  color="black", size=sz*ssca))+
     ggtitle(reduction)
    
    #ggplot(dat, aes(x=UMAP_1, y=UMAP_2, label=!!sym(labelby)))+
    #geom_point(size=sz, aes(color=!!sym(colorby)))+
    #scale_color_manual(values=colorlist[[colorby]], guide = guide_legend(override.aes = list(size = legend.size, shape=15) ))+
    #theme_classic()+
    #geom_text(data= dat2,  color="black", size=sz*ssca)
  },
  error = function(e) {
    if (grepl("Empty raster", e$message)) {
      message("geom_point_rast failed (Empty raster). Falling back to geom_point().")
        ggplot(dat, aes(x=UMAP_1, y=UMAP_2, label=!!sym(labelby)))+
    geom_point(size=sz, aes(color=!!sym(colorby)))+
    scale_color_manual(values=c(colorlist[[colorby]]), guide = guide_legend(override.aes = list(size = legend.size, shape=15) ))+
    theme_classic()+
    geom_text(data= dat2,  color="black", size=sz*ssca)+
          ggtitle(reduction)
    } else {
      stop(e) # rethrow unexpected errors
    }
  }
)

  
  if(!is.null(glasswork.by)){
    fcat(6)
    if(is.null(glasswork.params)){
      fcat("making glasswork hulls")
      glasswork.params=get.default.glasswork.params()
      glasswork.params$clus=dat %>% pull(!!sym(glasswork.by)) %>% unique %>% as.character
      fcat("categories to glasswork:", glasswork.params$clus)
    }
    
    #make glasswork hull df
    if(is.null(glasswork.df)){
       fcat("making glasswork df")
    glasswork.df=clusterhull4(so, clusvar = glasswork.by, fillvar = "null", reduction=reduction)
    
    }
    
    glasswork.df= glasswork.df %>% dplyr::mutate(glasswork.classes=!!sym(glasswork.by))
    print(glasswork.df)
    #make sure that we dont add twice the same colors
    if(colorby != glasswork.by){
      fcat("colorby and glasswork by are not the same")
    colorlist[["glasswork.classes"]]=colorlist[[glasswork.by]]
    }
    
    gpuc=gpuc+geom_polygon(data=glasswork.df %>% dplyr::filter(!!sym(glasswork.by) %in% glasswork.params$clus), inherit.aes=F, aes(x=UMAP_1, y=UMAP_2, color=glasswork.classes), fill=NA)#+
    scale_color_manual(values=c(colorlist[[colorby]], colorlist[["glasswork.classes"]]))
  }
  
    
  gpuc
  
  
  
  }



#################################

addprefix= (function(cc, prefix="C")  paste0(prefix, cc)) %>% Vectorize(., USE.NAMES=F)



givecolnames=function(x, ind=NULL, nms){
  if (is.null(ind)){
    ind= 1:length(nms) 
  }
  colnames(x)[ind]= nms
  x
  
}

giverownames=function(x, nms){
  rownames(x)= nms
  x
  
}

#####

findgeneindex=function(so, x)   (so %>% rownames) %>% grepl(paste0("^", x,"$"), .) %>% which


########################

na2zero=function(x) ifelse(is.na(x),0, x) 



maptoref= function(query, refdataset, mapvar, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NULL, mappinglabel=paste0(querylabel,"_to_",reflabel)){
  ###if query or ref do not have an SCT refAssay, compute them
  nameflag=0
  if ((intersect(query %>% colnames, refdataset %>% colnames) %>% length)>0){
    colnames(query)= paste0("query_", colnames(query))
    
    colnames(refdataset)=paste0("ref_", colnames(refdataset))
    nameflag=1
  }             
  
  cat("Making reference...\n\n\n")
  DefaultAssay(query)="refAssay"
  DefaultAssay(refdataset)="refAssay"
  refdataset <- RunPCA(refdataset, features = VariableFeatures(refdataset), verbose=F)
  refdataset <- RunUMAP(refdataset, reduction = "pca", dims = pdms, , assay="refAssay", return.model=TRUE, verbose=F)
  ref <- AzimuthReference(
    object = refdataset,
    refUMAP = "umap",
    refDR = "pca",
    refAssay = "refAssay",
    metadata = mapvar,
    dims = pdms,
    k.param = 31,
    reference.version = "1.0.0"
  )
  cat("Finding transfer anchors...\n\n\n")
  anchors <- FindTransferAnchors(
    reference = ref,
    query = query,
    k.filter = NA,
    reference.neighbors = "refdr.annoy.neighbors",
    reference.assay = "refAssay",
    query.assay = "refAssay",
    reference.reduction = "refDR",
    normalization.method = "SCT",
    features = intersect(rownames(x = ref), VariableFeatures(object = query)),
    dims = pdms,
    n.trees = 20,
    mapping.score.k = 100
  )
  cat("Transfering data...\n\n\n")
  query <- TransferData(
    reference = ref,
    query = query,
    dims = pdms,
    anchorset = anchors,
    refdata = mapvar,
    n.trees = 20,
    store.weights = TRUE
  )
  
  
  cat("Integrating embeddings...\n\n\n")
  # Calculate the embeddings of the query data on the reference SPCA
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = ref,
    query = query,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE,
    dims.to.integrate=pdms
  )
  
  cat("Finding query's neighbors in ref...\n\n\n")
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(ref[["refDR"]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  cat("Azimuth nntransform...\n\n\n")
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = ref[[]]
  )
  
  # Project the query to the reference UMAP.
  cat("projecting query into the reference umap...\n\n\n")
  query[["mapfun.umap"]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = ref[["refUMAP"]],
    reduction.key = 'mapfun.umap_',
    dims=pdms
  )
  
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )
  
  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  refdataset[["origin"]]=reflabel
  
  refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_1")], col.name="refumap_1")
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_2")], col.name="refumap_2")
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
  if(!is.na(customcol)){
    getcols=function(x) x[, c("origin", "refumap_1", "refumap_2",nm1, nm2, nm3, customcol)]
  }else{
    getcols=function(x) x[, c("origin", "refumap_1", "refumap_2",nm1, nm2, nm3)  ]
  }
  if(return.merge.metadata){
    
    
    #B %>% left_join(A )
    rbind(query %>% metadata %>% getcols , refdataset %>% metadata %>% expandcols %>% getcols)
    
  }else{
    if (return.query.metadata){ 
      if (nameflag==1){
        colnames(query)= sapply(strsplit(colnames(query), split="query_"), function(x) x[2])
      }
      
      
      if(return.lean){
        
        return(
          
          query %>% metadata  %>% getcols
          
        )}else{ return(query %>% metadata)}
      
      
      
    }
  }
}



maptoref2= function(query, refdataset, ref.assay="refAssay", query.assay="refAssay",  mapvar, normalization.method="SCT", reference.pca="pca", reference.umap="umap", project.query=FALSE,n.trees=20, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NA, mappinglabel=paste0(querylabel,"_to_",reflabel), process.reference=F){
  ###if query or ref do not have an SCT refAssay, compute them
  nameflag=0
  if ((intersect(query %>% colnames, refdataset %>% colnames) %>% length)>0){
    colnames(query)= paste0("query_", colnames(query))
    
    colnames(refdataset)=paste0("ref_", colnames(refdataset))
    nameflag=1
  }             
  
  
  cat("Making reference...\n\n\n")
  DefaultAssay(query)=query.assay
  
  if(process.reference==T){ # this assumes that the assay exists
    
    if(!(ref.assay %in% refdataset@assays %>% names)){
      error(paste0("Reference assay", ref.assay, " does not exist. Make sure to include this assay before running."))
    }else{
      
      DefaultAssay(refdataset)=ref.assay  
      refdataset <- RunPCA(refdataset, features = VariableFeatures(refdataset), verbose=F, reduction.name="mappingpca")
      refdataset <- RunUMAP(refdataset, reduction = "mappingpca", dims = pdms, , assay=ref.assay, return.model=TRUE, verbose=F)
    }
  }
  
  cat("Finding transfer anchors...\n\n\n")
  
  
  anchors <- FindTransferAnchors(
    reference = refdataset,
    query = query,
    k.filter = NA,
    reduction="pcaproject", #
    reference.neighbors = NULL,
    reference.assay = ref.assay,
    query.assay = query.assay,
    reference.reduction = reference.reduction,
    project.query=project.query,
    normalization.method = normalization.method,
    features = intersect(rownames(x = refdataset), VariableFeatures(object = query)),
    dims = pdms,
    n.trees = n.trees,
    mapping.score.k = 100
  )
  cat("Transfering data...\n\n\n")
  query <- TransferData(
    reference = refdataset,
    query = query,
    dims = pdms,
    anchorset = anchors,
    refdata = mapvar,
    n.trees = n.trees,
    store.weights = TRUE
  )
  
  
  cat("Integrating embeddings...\n\n\n")
  # Calculate the embeddings of the query data on the reference SPCA
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = refdataset,
    query = query,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE,
    dims.to.integrate=pdms,
    new.reduction.name="integrated_dr"
  )
  
  cat("Finding query's neighbors in ref...\n\n\n")
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(refdataset[["mappingpca"]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  cat("Azimuth nntransform...\n\n\n")
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = refdataset[[]]
  )
  
  # Project the query to the reference UMAP.
  cat("projecting query into the reference umap...\n\n\n")
  query[["mapfun.umap"]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = refdataset[["refUMAP"]],
    reduction.key = 'mapfun.umap_',
    dims=pdms
  )
  
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )
  
  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  refdataset[["origin"]]=reflabel
  
  refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_1")], col.name="refumap_1")
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_2")], col.name="refumap_2")
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
  if(!is.na(customcol)){
    getcols=function(x) x[, c("origin", "refumap_1", "refumap_2",nm1, nm2, nm3, customcol)]
  }else{
    getcols=function(x) x[, c("origin", "refumap_1", "refumap_2",nm1, nm2, nm3)  ]
  }
  if(return.merge.metadata){
    
    
    #B %>% left_join(A )
    rbind(query %>% metadata %>% getcols , refdataset %>% metadata %>% expandcols %>% getcols)
    
  }else{
    if (return.query.metadata){ 
      if (nameflag==1){
        colnames(query)= sapply(strsplit(colnames(query), split="query_"), function(x) x[2])
      }
      
      
      if(return.lean){
        
        return(
          
          query %>% metadata  %>% getcols
          
        )}else{ return(query %>% metadata)}
      
      
      
    }
  }
}



maptoref.pt1= function(query, 
                       refdataset, 
                       ref.assay="refAssay",
                       query.assay="refAssay",
                       mapvar, 
                       normalization.method="SCT",
                       reference.pca="pca", 
                       reference.umap="umap",
                       reference.reduction=reference.pca,
                       reference.neighbors=NULL,
                       project.query=FALSE,
                       n.trees=20,
                       pdms=1:50,
                       return.merge.metadata=T, 
                       return.query.metadata=F,
                       return.lean=T,
                       querylabel="query", 
                       reflabel="reference", 
                       customcol=NULL, 
                       mappinglabel=paste0(querylabel,"_to_",reflabel),
                       process.reference=F){
  ###if query or ref do not have an SCT refAssay, compute them
  
  nameflag=0
  if ((intersect(colnames(query), colnames(refdataset)) %>% length )>0){
    msg("labelling cell ids according to dataset of origin...")
    query=RenameCells(query, new.name= paste0("query_", colnames(query)))
    
    refdataset=RenameCells(refdataset, new.name= paste0("ref_", colnames(refdataset)))
    nameflag=1
  }             
  fcat("finished labeling cell ids")
  
  cat("Making reference...\n\n\n")
  DefaultAssay(query)=query.assay
  
  if(process.reference==T){ # this assumes that the assay exists
    
    if(!(ref.assay %in% refdataset@assays %>% names)){
      error(paste0("Reference assay", ref.assay, " does not exist. Make sure to include this assay before reprocessing."))
    }else{
      
      DefaultAssay(refdataset)=ref.assay  
      refdataset <- RunPCA(refdataset, features = VariableFeatures(refdataset), verbose=F, reduction.name="mappingpca")
      refdataset <- RunUMAP(refdataset, reduction = "mappingpca", dims = pdms, , assay=ref.assay, return.model=TRUE, verbose=F)
    }
  }else{
    fcat("setting", reference.pca, "and", reference.umap, "as mapping reductions")
    refdataset@reductions$mappingpca= refdataset@reductions[[reference.pca]]
    refdataset@reductions$mappingumap= refdataset@reductions[[reference.umap]]
  }
  
  # if(!is.null(reference.neighbors){
  #  nnname=paste0(umapred, "_snn")
  #  fcat("Attempting to use", nnname, "as reference neighbors...")
  # reference.neighbors=nnname
  #}
  
  if(is.null(reference.reduction)){
    fcat("Setting", reference.pca, "as the reference reduction")
    reference.reduction=reference.pca
  }
  
  
  cat("Finding transfer anchors...\n\n\n")
  
  
  anchors <- FindTransferAnchors(
    reference = refdataset,
    query = query,
    k.filter = NA,
    reduction="pcaproject", # method of projection
    reference.neighbors = reference.neighbors,
    reference.assay = ref.assay,
    query.assay = query.assay,
    reference.reduction = reference.reduction,
    project.query=project.query,
    normalization.method = normalization.method,
    features = intersect(rownames(x = refdataset), VariableFeatures(object = query)),
    dims = pdms,
    n.trees = n.trees,
    mapping.score.k = 100
  )
  cat("Transfering data...\n\n\n")
  query <- TransferData(
    reference = refdataset,
    query = query,
    dims = pdms,
    anchorset = anchors,
    refdata = mapvar,
    n.trees = n.trees,
    store.weights = TRUE
  )
  list(query=query, anchors=anchors, ref=refdataset)
  
  
  
}





maptoref.pt2= function(query, refdataset, anchors, ref.assay="refAssay", query.assay="refAssay",  mapvar, normalization.method="SCT", reference.pca="pca", reference.umap="umap",reference.reduction=reference.pca, reference.neighbors=NULL, project.query=FALSE,n.trees=20, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NULL, mappinglabel=paste0(querylabel,"_to_",reflabel), process.reference=F){
  
  
  query=Seurat::MapQuery(
    anchor=anchors,
    query=query,
    reference=refdataset,
    refdata = mapvar,
    new.reduction.name = "mapfun.pca",
    reference.reduction = reference.reduction,
    reference.dims = pdms,
    query.dims = pdms,
    reduction.model = reference.umap,
    transferdata.args = list(),
    integrateembeddings.args = list(),
    projectumap.args = list(
      reduction.name=reference.umap,
      reduction.key= prepare.rk(reference.umap)
    ),
    verbose = TRUE
  )   
  
  
  
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )    
  
  
  fcat("preparing the merge step 1...")
  query[["origin"]]=querylabel
  refumap1=prepare.rk(reference.umap) %>% paste0(., "1")
  refumap2=prepare.rk(reference.umap) %>% paste0(., "2")
  fcat("preparing the merge step 2...")
  #refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  fcat("preparing the merge step 3...")
  fcat("adding", refumap1)
  query=AddMetaData(query, reddf(query, reference.umap)[, refumap1], col.name=refumap1)
  fcat("adding", refumap2)
  query=AddMetaData(query, reddf(query, reference.umap)[, refumap2], col.name=refumap2)
  fcat("preparing the merge step 4...")
  fcat("adding", paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  fcat("adding", paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  fcat("adding", paste0("mapfun_",mapvar, "_mapping.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  fcat("preparing the merge step 5...")
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  fcat("preparing the merge step 6...")
  query

  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
  if(!is.null(customcol)){
    getcols=function(x) x[, c("origin", refumap1, refumap2,nm1, nm2, nm3, customcol)]
  }else{
    getcols=function(x) x[, c("origin", refumap1, refumap2,nm1, nm2, nm3)  ]
  }
  if(return.merge.metadata){
    
    fcat("preparing the merge step 7...")
    
    missing_cols <- setdiff(c(refumap1, refumap2,nm1, nm2, nm3, customcol), colnames(df))
    fcat("missing cols:", missing_cols)
    #B %>% left_join(A )
    fcat("refumap1",refumap1,  "refumap2", refumap2, "nm1", nm1,"nm2", nm2, "nm3",nm3,"customcol", customcol)
    
    rbind(query %>% metadata %>% getcols , refdataset %>% metadata %>% expandcols %>% getcols)
    
    
    
    
  }else{
    if (return.query.metadata){ 
      #if (nameflag==1){
      # colnames(query)= sapply(strsplit(colnames(query), split="query_"), function(x) x[2])
      #}
      
      
      if(return.lean){
        fcat("preparing the merge step 8...")
        return(
          query %>% metadata  %>% getcols
          
        )}else{
          fcat("preparing the merge step 9...")
          return(query %>% metadata)}
      
      
      
    }
  }   
  
}
  
  






maptoref.workflow=function(query, refdataset, ...){
  #function(query, refdataset, ref.assay="refAssay", query.assay="refAssay",  mapvar, normalization.method="SCT", reference.pca="pca", reference.umap="umap", project.query=FALSE,n.trees=20, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NA, mappinglabel=paste0(querylabel,"_to_",reflabel), process.reference=F)
  mappedcells=maptoref.pt1(query, refdataset, ...)
  mapmetadata=maptoref.pt2(query=mappedcells$query, anchors=mappedcells$anchors, refdataset=mappedcells$ref, ...)
  mapmetadata
  
}

maptoref.workflow1=function(query, refdataset, ...){
  #function(query, refdataset, ref.assay="refAssay", query.assay="refAssay",  mapvar, normalization.method="SCT", reference.pca="pca", reference.umap="umap", project.query=FALSE,n.trees=20, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NA, mappinglabel=paste0(querylabel,"_to_",reflabel), process.reference=F)
  mappedcells=maptoref.pt1(query, refdataset, ...)
  mappedcells
}  

maptoref.workflow2=function(mappedcells, ...){
  mapmetadata=maptoref.pt2(query=mappedcells$query, anchors=mappedcells$anchors, refdataset=mappedcells$ref, ...)
  mapmetadata
  
}


labeltransfer= function(query, 
                        refdataset, 
                        ref.assay="refAssay",
                        query.assay="refAssay",
                        mapvar, 
                        normalization.method="SCT",
                        reference.pca="pca", 
                        reference.umap="umap",
                        reference.reduction=reference.pca,
                        reference.neighbors=NULL,
                        project.query=FALSE,
                        n.trees=20,
                        pdms=1:50,
                        return.lean=T,
                        querylabel="query", 
                        reflabel="reference", 
                        customcol=NULL, 
                        mappinglabel=paste0(querylabel,"_to_",reflabel),
                        process.reference=F){
  ###if query or ref do not have an SCT refAssay, compute them
  
  nameflag=0
  if ((intersect(colnames(query), colnames(refdataset)) %>% length )>0){
    msg("labelling cell ids according to dataset of origin...")
    query=RenameCells(query, new.name= paste0("query_", colnames(query)))
    
    refdataset=RenameCells(refdataset, new.name= paste0("ref_", colnames(refdataset)))
    nameflag=1
  }             
  fcat("finished labeling cell ids")
  
  cat("Making reference...\n\n\n")
  DefaultAssay(query)=query.assay
  
  if(process.reference==T){ # this assumes that the assay exists
    
    if(!(ref.assay %in% refdataset@assays %>% names)){
      error(paste0("Reference assay", ref.assay, " does not exist. Make sure to include this assay before reprocessing."))
    }else{
      
      DefaultAssay(refdataset)=ref.assay  
      refdataset <- RunPCA(refdataset, features = VariableFeatures(refdataset), verbose=F, reduction.name="mappingpca")
      refdataset <- RunUMAP(refdataset, reduction = "mappingpca", dims = pdms, , assay=ref.assay, return.model=TRUE, verbose=F)
    }
  }else{
    fcat("setting", reference.pca, "and", reference.umap, "as mapping reductions")
    refdataset@reductions$mappingpca= refdataset@reductions[[reference.pca]]
    refdataset@reductions$mappingumap= refdataset@reductions[[reference.umap]]
  }
  
  # if(!is.null(reference.neighbors){
  #  nnname=paste0(umapred, "_snn")
  #  fcat("Attempting to use", nnname, "as reference neighbors...")
  # reference.neighbors=nnname
  #}
  
  if(is.null(reference.reduction)){
    fcat("Setting", reference.pca, "as the reference reduction")
    reference.reduction=reference.pca
  }
  
  
  cat("Finding transfer anchors...\n\n\n")
  
  
  anchors <- FindTransferAnchors(
    reference = refdataset,
    query = query,
    k.filter = NA,
    reduction="pcaproject", # method of projection
    reference.neighbors = reference.neighbors,
    reference.assay = ref.assay,
    query.assay = query.assay,
    reference.reduction = reference.reduction,
    project.query=project.query,
    normalization.method = normalization.method,
    features = intersect(rownames(x = refdataset), VariableFeatures(object = query)),
    dims = pdms,
    n.trees = n.trees,
    mapping.score.k = 100
  )
  cat("Transfering data...\n\n\n")
  query <- TransferData(
    reference = refdataset,
    query = query,
    dims = pdms,
    anchorset = anchors,
    refdata = mapvar,
    n.trees = n.trees,
    store.weights = TRUE
  )
  #list(query=query, anchors=anchors, ref=refdataset)
  
  
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )    
  
  
  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  
  #refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
  if(!is.null(customcol)){
    getcols=function(x) x[, c("origin", nm1, nm2, nm3, customcol)]
  }else{
    getcols=function(x) x[, c("origin",nm1, nm2, nm3)  ]
  }
  
  
  if(return.lean){
    
    return(
      
      query %>% metadata  %>% getcols
      
    )}else{ return(query %>% metadata)}
  
  
  
}   












maptoref.pt3= function(query, refdataset, anchors, ref.assay="refAssay", query.assay="refAssay",  mapvar, normalization.method="SCT", reference.pca="pca", reference.umap="umap", project.query=FALSE,n.trees=20, pdms=1:50, return.merge.metadata=T, return.query.metadata=F, return.lean=T, querylabel="query", reflabel="reference", customcol=NA, mappinglabel=paste0(querylabel,"_to_",reflabel), process.reference=F){
  
  cat("Finding query's neighbors in ref...\n\n\n")
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(refdataset[[reference.pca]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  
  
  
  
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  cat("Azimuth nntransform...\n\n\n")
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = refdataset[[]]
  )
  
  # Project the query to the reference UMAP.
  cat("projecting query into the reference umap...\n\n\n")
  query[["mapfun.umap"]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = refdataset[["refUMAP"]],
    reduction.key = 'mapfun.umap_',
    dims=pdms
  )
  
  cat("adding mapping score...\n\n\n")
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )
  
  cat("preparing the merge...\n\n\n")
  query[["origin"]]=querylabel
  refdataset[["origin"]]=reflabel
  
  refdataset=AddMetaData(refdataset, reddf(ref, "refUMAP")[, c("refumap_1","refumap_2")])
  refdataset=AddMetaData(refdataset, metadata(refdataset)[, c(mapvar)], col.name=paste0("mapfun_", mapvar))
  
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_1")], col.name="refumap_1")
  query=AddMetaData(query, reddf(query, "mapfun.umap")[, c("mapfunumap_2")], col.name="refumap_2")
  query=AddMetaData(query, metadata(query)[, c("predicted.id")], col.name= paste0("mapfun_",mapvar))
  query=AddMetaData(query, metadata(query)[, c("predicted.id.score")], col.name= paste0("mapfun_",mapvar, "_predicted.id.score"))
  query=AddMetaData(query, metadata(query)[, c("mapping.score")], col.name= paste0("mapfun_",mapvar, "_mapping.score"))
  
  nm1=paste0("mapfun_",mapvar)
  nm2=paste0("mapfun_",mapvar, "_predicted.id.score")
  nm3=paste0("mapfun_",mapvar, "_mapping.score")
  
  
  expandcols=function(x){
    x[,nm2]=NA
    x[,nm3]=NA
    
    x}
  
  
  if(!is.null(customcol)){
    getcols=function(x) x[, c("origin",nm1, nm2, nm3, customcol)]
  }else{
    getcols=function(x) x[, c("origin",nm1, nm2, nm3)  ]
  }
  if(return.merge.metadata){
    
    
    #B %>% left_join(A )
    rbind(query %>% metadata %>% getcols , refdataset %>% metadata %>% expandcols %>% getcols)
    
  }else{
    if (return.query.metadata){ 
      if (nameflag==1){
        colnames(query)= sapply(strsplit(colnames(query), split="query_"), function(x) x[2])
      }
      
      
      if(return.lean){
        
        return(
          
          query %>% metadata  %>% getcols
          
        )}else{ return(query %>% metadata)}
      
      
      
    }
  }
}

#######################################################################################################
#Function to generate hulls from umap clusters that light up according to some metric
#######################################################################################################
clusterhullfill=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar=NULL, fillcols=NULL, size=0.05){
  # if squares...
  #maxx=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_1) %>% max
  #maxy=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>%pull(UMAP_2) %>% max
  #minx=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_1) %>% min
  #miny=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_2) %>% min
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf
    
  }
  
  if(is.null(clus)){
    clus=umap.df %>% pull(clusvar) %>% unique %>% sort
    
  }
  
  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% filter(!!sym(clusvar)==clus[cc]) %>% select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    clhull=cl[chull(cl), ]
    
    if(is.null(pols)){
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      pols<-append(pols, geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }
    
    if(cc==length(clus)){
      return(pols)
    }else{
      return(catpols(cc+1, pols))}
    
    
  }
  
  catpols(1)
}



################################################################################
#
################################################################################

clusterhullfill2=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar=NULL, fillcols=NULL, size=0.05, concavity=2){
  # if squares...
  #maxx=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_1) %>% max
  #maxy=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>%pull(UMAP_2) %>% max
  #minx=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_1) %>% min
  #miny=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% pull(UMAP_2) %>% min
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf
    
  }
  
  if(is.null(clus)){
    clus=umap.df %>% pull(clusvar) %>% unique %>% sort
    
  }
  
  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% filter(!!sym(clusvar)==clus[cc]) %>% select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    #clhull=cl[chull(cl), ]
    clhull=concaveman(cl, concavity=concavity)
    
    if(is.null(pols)){
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      pols<-append(pols, geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }
    
    if(cc==length(clus)){
      return(pols)
    }else{
      return(catpols(cc+1, pols))}
    
    
  }
  
  catpols(1)
}

#############################################################################
#Code for seurat utilities which are useful
#############################################################################
NoLegend=function (...)
{
  no.legend.theme <- theme(legend.position = "none", validate = TRUE,
                           ...)
  return(no.legend.theme)
}


NoAxes=function (..., keep.text = FALSE, keep.ticks = FALSE)
{
  blank <- element_blank()
  no.axes.theme <- theme(axis.line.x = blank, axis.line.y = blank,
                         validate = TRUE, ...)
  if (!keep.text) {
    no.axes.theme <- no.axes.theme + theme(axis.text.x = blank,
                                           axis.text.y = blank, axis.title.x = blank, axis.title.y = blank,
                                           validate = TRUE, ...)
  }
  if (!keep.ticks) {
    no.axes.theme <- no.axes.theme + theme(axis.ticks.x = blank,
                                           axis.ticks.y = blank, validate = TRUE, ...)
  }
  return(no.axes.theme)
}



#######################importing manually functions to run hull sample from mvgps

hull_sample=function (X, num_grid_pts = 500, grid_type = "regular", 
                      trim_hull = FALSE, trim_quantile = NULL) 
{
  X_rslt <- X_check(X)
  assign("X", X_rslt$X)
  assign("m", X_rslt$m)
  grid_type <- match.arg(grid_type, choices = c("regular", 
                                                "random", "hexagonal"))
  if (trim_hull == TRUE) {
    if (is.null(trim_quantile)) 
      stop("trim_hull set to TRUE but trim_quantile not specified.", 
           call. = FALSE)
    if (trim_quantile < 0.5 | trim_quantile > 1) 
      stop("trim_quantile must be between [0.5, 1]", 
           call. = FALSE)
    trim_upper <- apply(X, 2, quantile, trim_quantile)
    trim_lower <- apply(X, 2, quantile, 1 - trim_quantile)
    X_trim <- sapply(seq_len(m), function(x) {
      ifelse(X[, x] > trim_upper[x], NA, ifelse(X[, x] < 
                                                  trim_lower[x], NA, X[, x]))
    })
    colnames(X_trim) <- colnames(X)
    X <- na.omit(X_trim)
  }
  if (m == 2) {
    hpts <- chull(X)
    hpts <- c(hpts, hpts[1])
    hpts_vs <- as.matrix(X[hpts, ])
    m <- Polygon(hpts_vs)
    ps <- Polygons(list(m), 1)
    sps <- SpatialPolygons(list(ps))
    sp_grid_pts <- spsample(sps, n = num_grid_pts, type = grid_type)
    grid_pts <- coordinates(sp_grid_pts)
    colnames(grid_pts) <- colnames(X)
  }
  else {
    hpts <- geometry::convhulln(X)
    hpts_ind <- unique(c(hpts))
    hpts_vs <- X[hpts_ind, ]
    grid_pts <- NULL
  }
  return(list(hpts_vs = hpts_vs, grid_pts = grid_pts, X = X))
}



X_check=function (X) 
{
  X <- as.matrix(X)
  m <- ncol(X)
  if (!any(apply(X, 2, is.numeric))) 
    stop("X must be numeric", call. = FALSE)
  if (m < 2) 
    stop("Exposure is not multivariate", call. = FALSE)
  return(list(X = X, m = m))
}


####################################################################
#function to make cluster hulls filled with points that signal a correspondence of cells, with points jittered to minimise overplotting effects
####################################################################


clusterhullpoint=function(so, umap.df= NULL, query.df, clusvar, clus, linecolor="white" ,fillvar=NULL, colorpoints="stage", fillcols=NULL, size=0.05, pointsize=0.05, reduction.key="UMAP_", concavity=2, downsample=1, rescale.cells=NULL){
  require(sf)
  
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf %>% dplyr::mutate(UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))
    
  }else{
    umap.df= umap.df  %>% dplyr::mutate(UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))
  }
  
  
  destinyclusters=query.df %>% group_by( !!sym(mapfun(clusvar))) %>% summarise(counts=n()) %>% as.data.frame %>% arrange(-counts)
  
  ##collect all the destinies of each cell in each  together
  
  if(is.null(rescale.cells)){
    des=destinyclusters %>% dplyr::mutate(temptot= sum(counts), frac=counts/temptot, newcounts= ceiling(frac*temptot) ) %>% col2names(., mapfun(clusvar))
  }else{
    des=destinyclusters %>% dplyr::mutate(temptot= sum(counts), newtot=!!rescale.cells, frac=counts/temptot, newcounts= ceiling(frac*newtot) ) %>% col2names(., mapfun(clusvar))
    
  }
  
  
  fcat("printing normal des:")
  print(destinyclusters)
  if(downsample!=1){
    fcat("printing downsampled des:")
    des= des %>% dplyr::mutate(counts=ceiling(counts*downsample))
    print(des)
  }
  getclusterfrac=function(x) des[x, "frac"]
  umap.df[["clusterfrac"]]=NA
  umap.df[["clusterfrac"]]=sapply(umap.df %>% pull(clusvar), getclusterfrac) %>% unname
  
  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% dplyr::mutate(null="one") %>% filter(!!sym(clusvar)==clus[cc]) %>% dplyr::select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    #clhull=cl[chull(cl), ] #obsolete but keep in case we need to revert to convex hulls
    
    
    
    clhull=concaveman(as.matrix(cl[, c("UMAP_1", "UMAP_2")]), concavity=concavity)
    
    rownames(clhull)= 1:nrow(clhull) %>% addprefix(., "hullpoint_")
    
    clhull=as.data.frame(clhull) %>% givecolnames(., nms=c("UMAP_1", "UMAP_2")) %>% dplyr::mutate(null="one")
    
    
    dists=dist(rbind(clhull[, c("UMAP_1", "UMAP_2")], cl[,c("UMAP_1", "UMAP_2")])) %>% as.matrix
    ld= nrow(dists)
    lc= nrow(clhull)
    inds=lapply(rownames(clhull), function(r){ names(dists[r, (lc+1):ld])[findmin(dists[r,(lc+1):ld])]}) %>% unlist
    clhull[[fillvar]]= cl[inds, ][[fillvar]]
    clhull[[mapfun(varr)]]= cl[inds, ][[mapfun(varr)]]
    clhull[[varr]]= cl[inds, ][[mapfun(varr)]]
    
    
    #print(clhull)
    ##getting number of points and sampling the hull space uniformly for the number of points found there
    fcat("number of cells in", clus[cc], ":", query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow())
    #at this point downsampling should already have been taken care of
    numcells=des %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% pull(newcounts)
    originalcells=des %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% pull(counts)
    if(length(numcells)==0){
      numcells=0
    }
    #numcells=query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow()
    subsample=(numcells*downsample) %>% round
    #fcat("numcells is" , numcells)
    if(!is.null(colorpoints) & !is.na(numcells) & numcells>2){
      #hsamples=hull_sample(X=clhull[, c(1,2)] %>% as.matrix, query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow() , grid_type="random")$grid_pts
      
      polygon = sf::st_polygon(list(clhull[, c("UMAP_1","UMAP_2")] %>% data.matrix))
      
      # Sample 50 random points within the polygon
      hsamples = sf::st_sample(polygon, size=numcells)
      hsamples= sf::st_coordinates(hsamples) %>% givecolnames(., nms=c("UMAP_1","UMAP_2"))
      
      
      fcat("samples computed...")
      fcat("number of cells in hull sample for", clus[cc], ":", hsamples %>% nrow())
      #print(hsamples)
      ###making matrix of mock coordinates for each cell
      
      if(subsample==originalcells){
        rplce=F}else{rplce=T}
      pointsdf=query.df %>% dplyr::mutate(null=NA) %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% sample_n(subsample, replace=rplce) %>% dplyr::mutate(UMAP_1=hsamples[,1], UMAP_2=hsamples[,2]) 
      
      pointsdf[[varr]] = pointsdf[[mapfun(varr)]]
      
      if(is.null(pols)){
        fcat("creating poygon list...")
        pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size),
                   geom_point_rast(data=pointsdf, inherit.aes=T, size=pointsize,aes(x=UMAP_1, y=UMAP_2, color=!!sym(colorpoints))))
      }else{
        fcat("appending poygon to list...")
        pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size),
                                geom_point_rast(data=pointsdf,inherit.aes=T,size=pointsize, aes(x=UMAP_1, y=UMAP_2, color=!!sym(colorpoints)))))
      }
      
      if(cc==length(clus)){
        fcat("finishing up. Attaching color scheme...")
        return(append(pols, scale_color_manual(values=allcolors[[colorpoints]])))
        #return(pols)
      }else{
        return(catpols(cc+1, pols))}
      
    }else{###if there are no query cells in cluster
      
      
      if(is.null(pols)){
        fcat("creating poygon list...")
        pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
      }else{
        fcat("appending poygon to list...")
        pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size)))
      }
      
      if(cc==length(clus)){
        #msg("finishing up. Attaching color scheme...")
        fcat("finishing up. Attaching color scheme...")
        return(append(pols, scale_color_manual(values=allcolors[[colorpoints]])))
      }else{
        return(catpols(cc+1, pols))}
      
      
      
      
      
      
    }
    
    
    
  }
  
  catpols(1)
}



####just make cluster hulls
clusterhull=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar="null", colorpoints="stage", fillcols=NULL, size=0.05, pointsize=0.05, reduction.key="UMAP_"){
  require(sp)
  
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf     
  }
  
  umap.df= umap.df %>% dplyr::mutate(null=NA, UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))
  
  
  if(is.null(clus)){
    clus= umap.df %>% pull(clusvar) %>% unique %>% sort
    
  }
  
  if( is.null(allcolors[[fillvar]])){
    
    allcolors[[fillvar]]= randomcolors((umap.df %>% pull(fillvar) %>% unique %>% length) )
  }
  
  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% dplyr::mutate(null="one") %>% filter(!!sym(clusvar)==clus[cc]) %>% select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    clhull=cl[chull(cl), ]
    
    #print(clhull)
    ##getting number of points and sampling the hull space uniformly for the number of points found there
    #fcat("number of cells in", clus[cc], ":", query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow())
    #numcells=0
    #numcells=query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow()
    
    
    #if(!is.null(colorpoints) && numcells>2){
    #hsamples=hull_sample(X=clhull[, c(1,2)] %>% as.matrix, query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow() , grid_type="random")$grid_pts
    
    #fcat("samples computed...")
    #fcat("number of cells in hull sample for", clus[cc], ":", hsamples %>% nrow())
    #print(hsamples)
    ###making matrix of mock coordinates for each cell
    #pointsdf=query.df %>% dplyr::mutate(null=NA) %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% dplyr::mutate(UMAP_1=hsamples[,1], UMAP_2=hsamples[,2])
    
    
    if(is.null(pols)){
      fcat("creating poygon list...")
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      fcat("appending poygon to list...")
      pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size)))
    }
    
    if(cc==length(clus)){
      fcat("finishing up. Attaching color scheme...")
      return(list(pols
                  #, scale_color_manual(values=allcolors[[colorpoints]])
                  , scale_fill_manual(values=allcolors[[fillvar]]) 
      ) %>% Reduce(append, .)
      )
      #return(pols)
    }else{
      return(catpols(cc+1, pols))}
  }
  catpols(1)
}

################################################################################
#uses concaveman for cluster hulls
################################################################################

clusterhull2=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar="null", colorpoints="stage", fillcols=NULL, size=0.05, pointsize=0.05, reduction.key="UMAP_", concavity=2){
  require(sp)
  
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf     
  }
  
  umap.df= umap.df %>% dplyr::mutate(null=NA, UMAP_1= !!sym(paste0(reduction.key, "1")), UMAP_2=!!sym(paste0(reduction.key, "2")))
  
  
  if(is.null(clus)){
    clus= umap.df %>% pull(clusvar) %>% unique %>% sort
    
  }
  
  if( is.null(allcolors[[fillvar]])){
    
    allcolors[[fillvar]]= randomcolors((umap.df %>% pull(fillvar) %>% unique %>% length) )
  }
  
  catpols=function(cc, pols=NULL){ 
    #if(is.null(fillvar)){
    cl=umap.df %>% dplyr::mutate(null="one") %>% filter(!!sym(clusvar)==clus[cc]) %>% select(UMAP_1, UMAP_2, !!sym(fillvar))
    #fillvar=NA
    #}else{
    #cl=so %>% umapdf %>% filter(!!sym(clusvar)==clus) %>% select(UMAP_1, UMAP_2, fillvar)
    #}
    
    #clhull=cl[chull(cl), ]
    
    clhull=concaveman(as.matrix(cl[, c("UMAP_1", "UMAP_2")]), concavity=concavity)
    
    rownames(clhull)= 1:nrow(clhull) %>% addprefix(., "hullpoint_")
    
    clhull=as.data.frame(clhull) %>% givecolnames(., nms=c("UMAP_1", "UMAP_2")) %>% dplyr::mutate(null="one")
    # fillvar on the concave matrix gets assiged the color of the closest point
    
    
    
    
    if(fillvar != "null"){
      dists=dist(rbind(clhull[, c("UMAP_1", "UMAP_2")], cl[,c("UMAP_1", "UMAP_2")])) %>% as.matrix
      ld= nrow(dists)
      lc= nrow(clhull)
      clhull[[fillvar]]= cl[lapply(rownames(clhull), function(r){
        names(dists[r, (lc+1):ld])[findmin(dists[r,(lc+1):ld])]  
        
      }) %>% unlist, ][[fillvar]]
      
    }
    
    
    
    #print(clhull)
    ##getting number of points and sampling the hull space uniformly for the number of points found there
    #fcat("number of cells in", clus[cc], ":", query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow())
    #numcells=0
    #numcells=query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow()
    
    
    #if(!is.null(colorpoints) && numcells>2){
    #hsamples=hull_sample(X=clhull[, c(1,2)] %>% as.matrix, query.df %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% nrow() , grid_type="random")$grid_pts
    
    #fcat("samples computed...")
    #fcat("number of cells in hull sample for", clus[cc], ":", hsamples %>% nrow())
    #print(hsamples)
    ###making matrix of mock coordinates for each cell
    #pointsdf=query.df %>% dplyr::mutate(null=NA) %>% filter(!!sym(mapfun(clusvar))==clus[cc]) %>% dplyr::mutate(UMAP_1=hsamples[,1], UMAP_2=hsamples[,2])
    
    
    if(is.null(pols)){
      fcat("creating poygon list...")
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      fcat("appending poygon to list...")
      pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size)))
    }
    
    if(cc==length(clus)){
      fcat("finishing up. Attaching color scheme...")
      return(list(pols
                  #, scale_color_manual(values=allcolors[[colorpoints]])
                  , scale_fill_manual(values=allcolors[[fillvar]]) 
      ) %>% Reduce(append, .)
      )
      #return(pols)
    }else{
      return(catpols(cc+1, pols))}
  }
  catpols(1)
}


################################################################################
# non recursive version of clusterhull2 to avoid weird R issues
################################################################################


clusterhull3=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar="null", colorpoints="stage", fillcols=NULL, size=0.05, pointsize=0.05, reduction="umap", reduction.key="UMAP_", concavity=2){
  require(sp)
  
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf(., reduction=reduction)     
  }
  
  redkeys=Embeddings(so, reduction) %>% colnames
  dim1=redkeys[1]
  dim2=redkeys[2]
  
  umap.df= umap.df %>% dplyr::mutate(null=NA, UMAP_1= !!sym(dim1), UMAP_2=!!sym(dim2))
  
  
  if(is.null(clus)){
    clus= umap.df %>% pull(clusvar) %>% unique %>% sort
    
  }
  
  if( is.null(allcolors[[fillvar]])){
    
    allcolors[[fillvar]]= randomcolors((umap.df %>% pull(fillvar) %>% unique %>% length) )
  }
  
  pols=NULL
  polygon=function (x, y = NULL, density = NULL, angle = 45, border = NULL, 
                    col = NA, lty = par("lty"), ..., fillOddEven = FALSE) 
  {
    ..debug.hatch <- FALSE
    xy <- xy.coords(x, y, setLab = FALSE)
    if (is.numeric(density) && all(is.na(density) | density < 
                                   0)) 
      density <- NULL
    if (!is.null(angle) && !is.null(density)) {
      polygon.onehatch <- function(x, y, x0, y0, xd, yd, ..debug.hatch = FALSE, 
                                   ...) {
        if (..debug.hatch) {
          points(x0, y0)
          arrows(x0, y0, x0 + xd, y0 + yd)
        }
        halfplane <- as.integer(xd * (y - y0) - yd * (x - 
                                                        x0) <= 0)
        cross <- halfplane[-1L] - halfplane[-length(halfplane)]
        does.cross <- cross != 0
        if (!any(does.cross)) 
          return()
        x1 <- x[-length(x)][does.cross]
        y1 <- y[-length(y)][does.cross]
        x2 <- x[-1L][does.cross]
        y2 <- y[-1L][does.cross]
        t <- (((x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - 
                                                      x1))/(xd * (y2 - y1) - yd * (x2 - x1)))
        o <- order(t)
        tsort <- t[o]
        crossings <- cumsum(cross[does.cross][o])
        if (fillOddEven) 
          crossings <- crossings%%2
        drawline <- crossings != 0
        lx <- x0 + xd * tsort
        ly <- y0 + yd * tsort
        lx1 <- lx[-length(lx)][drawline]
        ly1 <- ly[-length(ly)][drawline]
        lx2 <- lx[-1L][drawline]
        ly2 <- ly[-1L][drawline]
        segments(lx1, ly1, lx2, ly2, ...)
      }
      polygon.fullhatch <- function(x, y, density, angle, ..debug.hatch = FALSE, 
                                    ...) {
        x <- c(x, x[1L])
        y <- c(y, y[1L])
        angle <- angle%%180
        if (par("xlog") || par("ylog")) {
          warning("cannot hatch with logarithmic scale active")
          return()
        }
        usr <- par("usr")
        pin <- par("pin")
        upi <- c(usr[2L] - usr[1L], usr[4L] - usr[3L])/pin
        if (upi[1L] < 0) 
          angle <- 180 - angle
        if (upi[2L] < 0) 
          angle <- 180 - angle
        upi <- abs(upi)
        xd <- cos(angle/180 * pi) * upi[1L]
        yd <- sin(angle/180 * pi) * upi[2L]
        if (angle < 45 || angle > 135) {
          if (angle < 45) {
            first.x <- max(x)
            last.x <- min(x)
          }
          else {
            first.x <- min(x)
            last.x <- max(x)
          }
          y.shift <- upi[2L]/density/abs(cos(angle/180 * 
                                               pi))
          x0 <- 0
          y0 <- floor((min(y) - first.x * yd/xd)/y.shift) * 
            y.shift
          y.end <- max(y) - last.x * yd/xd
          while (y0 < y.end) {
            polygon.onehatch(x, y, x0, y0, xd, yd, ..debug.hatch = ..debug.hatch, 
                             ...)
            y0 <- y0 + y.shift
          }
        }
        else {
          if (angle < 90) {
            first.y <- max(y)
            last.y <- min(y)
          }
          else {
            first.y <- min(y)
            last.y <- max(y)
          }
          x.shift <- upi[1L]/density/abs(sin(angle/180 * 
                                               pi))
          x0 <- floor((min(x) - first.y * xd/yd)/x.shift) * 
            x.shift
          y0 <- 0
          x.end <- max(x) - last.y * xd/yd
          while (x0 < x.end) {
            polygon.onehatch(x, y, x0, y0, xd, yd, ..debug.hatch = ..debug.hatch, 
                             ...)
            x0 <- x0 + x.shift
          }
        }
      }
      if (missing(col) || is.null(col)) {
        col <- par("fg")
      }
      else if (any(is.na(col))) {
        col[is.na(col)] <- par("fg")
      }
      if (is.null(border)) 
        border <- col
      if (is.logical(border)) {
        if (!is.na(border) && border) 
          border <- col
        else border <- NA
      }
      start <- 1
      ends <- c(seq_along(xy$x)[is.na(xy$x) | is.na(xy$y)], 
                length(xy$x) + 1)
      num.polygons <- length(ends)
      col <- rep_len(col, num.polygons)
      if (length(border)) 
        border <- rep_len(border, num.polygons)
      if (length(lty)) 
        lty <- rep_len(lty, num.polygons)
      if (length(density)) 
        density <- rep_len(density, num.polygons)
      angle <- rep_len(angle, num.polygons)
      i <- 1L
      for (end in ends) {
        if (end > start) {
          if (is.null(density) || is.na(density[i]) || 
              density[i] < 0) 
            .External.graphics(C_polygon, xy$x[start:(end - 
                                                        1)], xy$y[start:(end - 1)], col[i], NA, lty[i], 
                               ...)
          else if (density[i] > 0) {
            polygon.fullhatch(xy$x[start:(end - 1)], xy$y[start:(end - 
                                                                   1)], col = col[i], lty = lty[i], density = density[i], 
                              angle = angle[i], ..debug.hatch = ..debug.hatch, 
                              ...)
          }
          i <- i + 1
        }
        start <- end + 1
      }
      .External.graphics(C_polygon, xy$x, xy$y, NA, border, 
                         lty, ...)
    }
    else {
      if (is.logical(border)) {
        if (!is.na(border) && border) 
          border <- par("fg")
        else border <- NA
      }
      .External.graphics(C_polygon, xy$x, xy$y, col, border, 
                         lty, ...)
    }
    invisible()
  }
  for(cc in 1:length(clus)){
    
    cl=umap.df %>% dplyr::mutate(null="one") %>% filter(!!sym(clusvar)==clus[cc]) %>% dplyr::select(UMAP_1, UMAP_2, !!sym(fillvar))
    
    
    clhull=concaveman(as.matrix(cl[, c("UMAP_1", "UMAP_2")]), concavity=concavity)
    
    rownames(clhull)= 1:nrow(clhull) %>% addprefix(., "hullpoint_")
    
    clhull=as.data.frame(clhull) %>% givecolnames(., nms=c("UMAP_1", "UMAP_2")) %>% dplyr::mutate(null="one")
    # fillvar on the concave matrix gets assiged the color of the closest point
    
    if(fillvar != "null"){
      dists=dist(rbind(clhull[, c("UMAP_1", "UMAP_2")], cl[,c("UMAP_1", "UMAP_2")])) %>% as.matrix
      ld= nrow(dists)
      lc= nrow(clhull)
      clhull[[fillvar]]= cl[lapply(rownames(clhull), function(r){
        names(dists[r, (lc+1):ld])[findmin(dists[r,(lc+1):ld])]  
        
      }) %>% unlist, ][[fillvar]]
      
    }
    
    if(is.null(pols)){
      fcat("creating poygon list...")
      pols<-list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size))
    }else{
      fcat("appending poygon to list...")
      pols<-append(pols, list(geom_polygon(data=clhull, inherit.aes=FALSE,aes(x=UMAP_1, y=UMAP_2, fill=!!sym(fillvar)), color=linecolor, size=size)))
    }
    
    if(cc==length(clus)){
      fcat("finishing up. Attaching color scheme...")
      return(list(pols
                  #, scale_color_manual(values=allcolors[[colorpoints]])
                  , scale_fill_manual(values=allcolors[[fillvar]]) 
      ) %>% Reduce(append, .)
      )
      #return(pols)
    }
  } 
  return(pols) 
}



clusterhull4=function(so, umap.df= NULL, clusvar, clus=NULL, linecolor="white" ,fillvar="null", colorpoints="stage", fillcols=NULL, size=0.05, pointsize=0.05, reduction="umap", reduction.key="UMAP_", concavity=2){
  require(sp)
  
  
  if(is.null(umap.df)){
    umap.df= so %>% umapdf(., reduction=reduction)     
  }
  
  redkeys=Embeddings(so, reduction) %>% colnames
  dim1=redkeys[1]
  dim2=redkeys[2]
  
  umap.df= umap.df %>% dplyr::mutate(null=NA, UMAP_1= !!sym(dim1), UMAP_2=!!sym(dim2))
  
  
  if(is.null(clus)){
    clus= umap.df %>% pull(clusvar) %>% unique %>% sort
    
  }
  
  if( is.null(allcolors[[fillvar]])){
    
    allcolors[[fillvar]]= randomcolors((umap.df %>% pull(fillvar) %>% unique %>% length) )
  }
  
  pols=NULL
  polygon=function (x, y = NULL, density = NULL, angle = 45, border = NULL, 
                    col = NA, lty = par("lty"), ..., fillOddEven = FALSE) 
  {
    ..debug.hatch <- FALSE
    xy <- xy.coords(x, y, setLab = FALSE)
    if (is.numeric(density) && all(is.na(density) | density < 
                                   0)) 
      density <- NULL
    if (!is.null(angle) && !is.null(density)) {
      polygon.onehatch <- function(x, y, x0, y0, xd, yd, ..debug.hatch = FALSE, 
                                   ...) {
        if (..debug.hatch) {
          points(x0, y0)
          arrows(x0, y0, x0 + xd, y0 + yd)
        }
        halfplane <- as.integer(xd * (y - y0) - yd * (x - 
                                                        x0) <= 0)
        cross <- halfplane[-1L] - halfplane[-length(halfplane)]
        does.cross <- cross != 0
        if (!any(does.cross)) 
          return()
        x1 <- x[-length(x)][does.cross]
        y1 <- y[-length(y)][does.cross]
        x2 <- x[-1L][does.cross]
        y2 <- y[-1L][does.cross]
        t <- (((x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - 
                                                      x1))/(xd * (y2 - y1) - yd * (x2 - x1)))
        o <- order(t)
        tsort <- t[o]
        crossings <- cumsum(cross[does.cross][o])
        if (fillOddEven) 
          crossings <- crossings%%2
        drawline <- crossings != 0
        lx <- x0 + xd * tsort
        ly <- y0 + yd * tsort
        lx1 <- lx[-length(lx)][drawline]
        ly1 <- ly[-length(ly)][drawline]
        lx2 <- lx[-1L][drawline]
        ly2 <- ly[-1L][drawline]
        segments(lx1, ly1, lx2, ly2, ...)
      }
      polygon.fullhatch <- function(x, y, density, angle, ..debug.hatch = FALSE, 
                                    ...) {
        x <- c(x, x[1L])
        y <- c(y, y[1L])
        angle <- angle%%180
        if (par("xlog") || par("ylog")) {
          warning("cannot hatch with logarithmic scale active")
          return()
        }
        usr <- par("usr")
        pin <- par("pin")
        upi <- c(usr[2L] - usr[1L], usr[4L] - usr[3L])/pin
        if (upi[1L] < 0) 
          angle <- 180 - angle
        if (upi[2L] < 0) 
          angle <- 180 - angle
        upi <- abs(upi)
        xd <- cos(angle/180 * pi) * upi[1L]
        yd <- sin(angle/180 * pi) * upi[2L]
        if (angle < 45 || angle > 135) {
          if (angle < 45) {
            first.x <- max(x)
            last.x <- min(x)
          }
          else {
            first.x <- min(x)
            last.x <- max(x)
          }
          y.shift <- upi[2L]/density/abs(cos(angle/180 * 
                                               pi))
          x0 <- 0
          y0 <- floor((min(y) - first.x * yd/xd)/y.shift) * 
            y.shift
          y.end <- max(y) - last.x * yd/xd
          while (y0 < y.end) {
            polygon.onehatch(x, y, x0, y0, xd, yd, ..debug.hatch = ..debug.hatch, 
                             ...)
            y0 <- y0 + y.shift
          }
        }
        else {
          if (angle < 90) {
            first.y <- max(y)
            last.y <- min(y)
          }
          else {
            first.y <- min(y)
            last.y <- max(y)
          }
          x.shift <- upi[1L]/density/abs(sin(angle/180 * 
                                               pi))
          x0 <- floor((min(x) - first.y * xd/yd)/x.shift) * 
            x.shift
          y0 <- 0
          x.end <- max(x) - last.y * xd/yd
          while (x0 < x.end) {
            polygon.onehatch(x, y, x0, y0, xd, yd, ..debug.hatch = ..debug.hatch, 
                             ...)
            x0 <- x0 + x.shift
          }
        }
      }
      if (missing(col) || is.null(col)) {
        col <- par("fg")
      }
      else if (any(is.na(col))) {
        col[is.na(col)] <- par("fg")
      }
      if (is.null(border)) 
        border <- col
      if (is.logical(border)) {
        if (!is.na(border) && border) 
          border <- col
        else border <- NA
      }
      start <- 1
      ends <- c(seq_along(xy$x)[is.na(xy$x) | is.na(xy$y)], 
                length(xy$x) + 1)
      num.polygons <- length(ends)
      col <- rep_len(col, num.polygons)
      if (length(border)) 
        border <- rep_len(border, num.polygons)
      if (length(lty)) 
        lty <- rep_len(lty, num.polygons)
      if (length(density)) 
        density <- rep_len(density, num.polygons)
      angle <- rep_len(angle, num.polygons)
      i <- 1L
      for (end in ends) {
        if (end > start) {
          if (is.null(density) || is.na(density[i]) || 
              density[i] < 0) 
            .External.graphics(C_polygon, xy$x[start:(end - 
                                                        1)], xy$y[start:(end - 1)], col[i], NA, lty[i], 
                               ...)
          else if (density[i] > 0) {
            polygon.fullhatch(xy$x[start:(end - 1)], xy$y[start:(end - 
                                                                   1)], col = col[i], lty = lty[i], density = density[i], 
                              angle = angle[i], ..debug.hatch = ..debug.hatch, 
                              ...)
          }
          i <- i + 1
        }
        start <- end + 1
      }
      .External.graphics(C_polygon, xy$x, xy$y, NA, border, 
                         lty, ...)
    }
    else {
      if (is.logical(border)) {
        if (!is.na(border) && border) 
          border <- par("fg")
        else border <- NA
      }
      .External.graphics(C_polygon, xy$x, xy$y, col, border, 
                         lty, ...)
    }
    invisible()
  }
  
  clhulls=lapply(1:length(clus), function(cc){
  #for(cc in 1:length(clus)){
    
    cl=umap.df %>% dplyr::mutate(null="one") %>% filter(!!sym(clusvar)==clus[cc]) %>% dplyr::select(UMAP_1, UMAP_2, !!sym(fillvar), !!sym(clusvar))
    
    
    clhull=concaveman(as.matrix(cl[, c("UMAP_1", "UMAP_2")]), concavity=concavity)
    
    rownames(clhull)= 1:nrow(clhull) %>% addprefix(., "hullpoint_")
    
    clhull=as.data.frame(clhull) %>% givecolnames(., nms=c("UMAP_1", "UMAP_2")) %>% dplyr::mutate(null="one", !!sym(clusvar):=clus[cc])
    # fillvar on the concave matrix gets assiged the color of the closest point
    
    if(fillvar != "null"){
      dists=dist(rbind(clhull[, c("UMAP_1", "UMAP_2")], cl[,c("UMAP_1", "UMAP_2")])) %>% as.matrix
      ld= nrow(dists)
      lc= nrow(clhull)
      clhull[[fillvar]]= cl[lapply(rownames(clhull), function(r){
        names(dists[r, (lc+1):ld])[findmin(dists[r,(lc+1):ld])]  
        
      }) %>% unlist, ][[fillvar]]
      
    }
    clhull
  })
  clhulls %>% Reduce(rbind, .)
}
  




randomcolors= function(n)
  
  lapply(1:n, function(x) paste0("#",sample(c(0:9,LETTERS[1:6]), 1),
                                 sample(c(0:9,LETTERS[1:6]), 1),
                                 sample(c(0:9,LETTERS[1:6]), 1),
                                 sample(c(0:9,LETTERS[1:6]), 1),
                                 sample(c(0:9,LETTERS[1:6]), 1),
                                 sample(c(0:9,LETTERS[1:6]), 1))) %>% Reduce(c, .)




demapfun=function(x) strsplit(x, split="mapfun_")[[1]][2]
mapfun=function(x) x %>% addprefix(., "mapfun_")


#####################################

prefixids.smart= function(so){
  colnames(so@assays$RNA@counts)= paste_(Project(so),  colnames(so@assays$RNA@counts))
  so
}
prefixids= function(so, pref=""){
  if ( !("patientID" %in% (so %>% metadata %>% colnames)) ){
    if ( ("orig.ident" %in% (so %>% metadata %>% colnames)) ){
      RenameCells(so, new.names= so %>% metadata %>% dplyr::mutate(newcellid=paste0(pref,orig.ident,"_", so %>% colnames)) %>% pull(newcellid) )
    }else{
      RenameCells(so, new.names= so %>% metadata %>% dplyr::mutate(newcellid=paste0(pref, so %>% colnames)) %>% pull(newcellid) )
    }}else{
      RenameCells(so, new.names= so %>% metadata %>% dplyr::mutate(newcellid=paste0(pref, patientID, "_", so %>% colnames)) %>% pull(newcellid) )
    }
  
}

##assuming there is a patientID
prefixids= function(so, pref=""){
  RenameCells(so, new.names= so %>% metadata %>% dplyr::mutate(newcellid=paste0(pref, so %>% colnames)) %>% pull(newcellid) )
  
}


prefixids.matrix= function(mat, prefix){
  colnames(mat)= paste_(prefix,  colnames(mat))
  mat
}

prefixids.meta= function(mat, prefix){
  rownames(mat)= paste_(prefix,  rownames(mat))
  mat
}

alph=77
applyalpha=function(color) color %>% paste0(., alph) 


################################################################################
# SPECIFIC NCNB2 functions
################################################################################

# split the hto_demux name into different bits. 


extractname = Vectorize(function(x) {
  #cat(x, "\n")
  x=gsub("-", "_", x)
  if (!grepl("_", x)) {
    return(c(NA, NA, NA, NA) %>% givename(., c("cell.line", "condition", "stage", "replicate")))
  } else{
    #HTO-H7_17q_day3_rep1_CMO
    
    #small replacement on the condition section of the string
    #cat("1")
    if (!grepl("^H7", x)) {
      x0=paste0("H7_", x)
    }else{x0=x}
    
    x1 = gsub("H7WT", "H7_WT", x0)
    x2 = gsub("17q1qMYCN_", "17qx_", x1) #transform all the potential17q1qs into 17qx
    x3 = gsub("17qMYCN_", "17q1qMYCNdox_", x2) #transform all the potential mycns into mycndox
    x4 = gsub("MYCNdox_", "MYCN_", x3) #transform all dox into mycn
    x5 = gsub("17qx", "17q1q", x4) #tranform all qx into 17q1q
    
    
    x6 = gsub("_CMO", "", x5)
    x7= gsub("rep","R", x6)
    x8= gsub("day", "D", x7)
    
    xf=x8
    
    #make a list with all the parts thereafter
    xs1 = strsplit(xf, split = "_")[[1]]
    #print(x2)
    names(xs1) = c("cell.line", "condition", "stage", "replicate")
    xs2 = as.list(xs1)
    xs2$condition = paste0("c", xs2$condition)
    
    xs2 %>% unlist
  }
  
}, USE.NAMES=F)

################################################################################
# screen print shortcuts
################################################################################


#concatenate a message for the screen, add screen format

fcat= function(...) cat(paste(...,sep=" "), "\n")
msg=fcat
showdataset=function(so)
  fcat("Dataset", so$dsname %>% unique %>% paste(., collapse="-")) 


################################################################################
# dataset filtering on the go
################################################################################

filterds=function(so, dbvar="scdblscore", dbcutoff=NULL, filter.empty.drops=NULL, demuxvar="sample_name",qc.pars=NULL, keepvar=NULL, keep.classes=NULL){
  if(!is.null(qc.pars)){
    qcpars=qc.pars 
  }
  
  nm=make.names(Project(so))
  if(!(nm %in% names(qcpars)) && !(nm %in% make.names(names(qcpars)))){
    st=paste(nm, "does not have an assigned QC filter")
    stop(st)
  }else{
    
    fcat("Updating filter metadata...")  
    so=AddMetaData(so, so %>% metadata %>% dplyr::mutate(counts.cutoff=qcpars[[nm]]$counts, features.cutoff= qcpars[[nm]]$features, mito.cutoff=qcpars[[nm]]$mito, qcc.pass=((nFeature_RNA>=features.cutoff & nCount_RNA >= counts.cutoff & percent.mt<= mito.cutoff )))  %>% dplyr::select(counts.cutoff, features.cutoff, mito.cutoff, qcc.pass))
  }
  if(demuxvar %in% (so %>% metadata %>% colnames)){
    so=so[, so %>% metadata %>% filter(!(!!sym(demuxvar) %in% c("Negative", "Doublet", "Multiplet", "Unassigned", "Blank", NA)), nCount_RNA>=qcpars[[nm]]$counts, nFeature_RNA>=qcpars[[nm]]$features, percent.mt<=qcpars[[nm]]$mito ) %>% rownames] 
    
  }else{
    so=so[, so %>% metadata %>% filter(nCount_RNA>=qcpars[[nm]]$counts, nFeature_RNA>=qcpars[[nm]]$features, percent.mt<=qcpars[[nm]]$mito ) %>% rownames] 
  }  
  
  if( !is.null(dbcutoff) & (dbvar %in% (so %>% metadata %>% colnames))){
    so=so[, so %>% metadata %>% dplyr::filter(!!sym(dbvar)<=dbcutoff) %>% rownames]
    so=AddMetaData(so, metadata= so %>% metadata %>% dplyr::mutate(scDblFinder.cutoff=dbcutoff) %>% select(scDblFinder.cutoff))
  }
  if( !is.null(filter.empty.drops) & ("isemptydroplet" %in% (so %>% metadata %>% colnames))){
    so=so[, so %>% metadata %>% filter(isemptydroplet==FALSE) %>% rownames]
  }
  
  
  so   
}


filteremptydrops=function(so, filter.empty.drops=NULL){
  nm=Project(so)
  
  if( !is.null(filter.empty.drops) & ("isemptydroplet" %in% (so %>% metadata %>% colnames))){
    so=so[, so %>% metadata %>% filter(isemptydroplet==FALSE) %>% rownames]
  }
  
  
  so   
}



stagename= Vectorize(function(nm) stages[[nm]], USE.NAMES=F)


################################################################################
#reorder cells in clusters according to different methods (e.g., seriation or signature, i.e. how trong the markers are)
################################################################################

cellinfo.getcells=function(so, clusvar="seurat_clusters",ncells=100, clusters=NULL, fullrecreate=T, replicatevar=NULL, fullreload=F, ...){ #replicatevar=NULL,  clusters=NULL){
  
  Idents(so)=clusvar
  if(is.null(clusters)){
    clusids= so[[clusvar]] %>% pull(clusvar) %>% as.character %>% unique %>% sort
  }else{
    clusids=clusters 
  }
  
  fcat("Getting cells for each cluster...")
  cellist=lapply(as.character(clusids), function(x) {
    out=WhichCells(so, idents=x, downsample=ncells)
    #fcat("total number of cells:", length(out))
    out
  }
  
  ) %>% givename(., as.character(clusids))
  
  cellist
}


cellinfo.getmarkers=function(so, clusvar="seurat_clusters", selecttop=NULL){
  
  if(is.null(so$RNA@misc$top_markers)){  
    so=seuratmarkers.delegate(so, minfc=1, group_column=clusvar, replicate_column=replicatevar, fullrecreate=fullrecreate, fullreload=fullreload, ...)
  }else{
    
    
  }
  
}



################################################################################
# deduplicate marker list forcing the marker to be assigned to the highest expressing cluster
################################################################################

      deduplicate.markers=function(markers, gene.field="feature", logfc.field="log_fc", group.field="group1"){
      deduplicated_markers <- markers %>%
     group_by(!!sym(gene.field)) %>%
     slice_max(order_by = !!sym(logfc.field), n = 1, with_ties = FALSE) %>%
      ungroup() %>% as.data.frame %>% arrange(!!sym(group.field), -!!sym(logfc.field))
      deduplicated_markers
      }

harmonise.marker.table=function(m){
if("cluster" %in% grep("group1|cluster", colnames(m), value=T)){
 m = m %>% dplyr::mutate(group1=cluster, feature=gene, log_fc=avg_log2FC, rate1=pct.1, rate2=pct.2) 
}  
m  
}
harmonise.all.markers=function(so){
  
for(j in names(so$RNA@misc)){
 
  if(j %in% c("markers", "top_markers")){
    so$RNA@misc[[j]]=harmonise.marker.table(so$RNA@misc[[j]])
  }else{
    so$RNA@misc[[j]]$markers=harmonise.marker.table(so$RNA@misc[[j]]$markers)
    so$RNA@misc[[j]]$top_markers=harmonise.marker.table(so$RNA@misc[[j]]$top_markers)
    }

}
so}  
  



get.cellinfo.markerlist=function(so, clusvar="seurat_clusters", clusters=NULL, markers=NULL, ass="SCT", groupvarname="group1",lfcvarname="log_fc", extended.output=T, deduped=T, nmarkers=NULL){
  
  Idents(so)=clusvar
  if(is.null(clusters)){
    clusids= so[[clusvar]] %>% pull(clusvar) %>% as.character %>% unique %>% sort
  }else{
    clusids=clusters 
  }
  
  
  markers= so@assays$RNA@misc$top_markers
  
  if(is.null(markers)){
    stop("There are no markers! please add markers under @assays$RNA@misc$top_markers\n" ) 
  }else{
    fcat("Splitting groups...")
    markerlist0=lapply(clusids, function(x) markers %>% filter(!!sym(groupvarname)==x)) %>% givename(., clusids)
    
    hasgenes= Vectorize(function(tabb) as.logical(nrow(tabb)>0), USE.NAMES=F)
    
    
    fcat("current length of markerlist0", length(markerlist0))
    which.have.genes=hasgenes(markerlist0)
    markerlist0=markerlist0[which.have.genes]
    fcat("current length of markerlist0 after filtering", length(markerlist0))
    fcat("these clusters have genes:", paste(markerlist0 %>% names, collapse=" ")) 
    
    fcat("current length of clusids", length(clusids))
    clusids=as.character(clusids)[which.have.genes]
    fcat("current length of clusids after filter", length(clusids))
    names(markerlist0)= as.character(clusids) 
    
    
    if(deduped==T){

      ###when duplicates, keep the gene as marker in the cluster where it shows highest log_fold change
      fcat("Assigning genes in more than one cluster to the highest expressing cluster...")    
      #NEW CODE!!! hopefully more efficient
      
       markers.revised=deduplicate.markers(markerlist0 %>% Reduce(rbind, .), logfc.field=lfcvarname, group.field = groupvarname)
    }else{
      markers.revised=markers # forwarding the original table
    }
    
    ############################################################################
    # assemble final markerlist
    ############################################################################
    
        if(is.null(nmarkers)){    
      markerlist=lapply(as.character(clusids), function(x) markers.revised %>% filter(!!sym(groupvarname)==x) %>% pull(feature))
    }else{
      markerlist=lapply(as.character(clusids), function(x) markers.revised %>% filter(!!sym(groupvarname)==x) %>% pull(feature) %>% head(nmarkers))
    }
    
    names(markerlist)=as.character(clusids)  
    
    
    names(markerlist)=as.character(clusids)  
}
markerlist  

}

cellinfo.get.markerlist=get.cellinfo.markerlist



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


seriatecells=function(so, clusvar="seurat_clusters", clusters=NULL, meth="seurat", markers=NULL, ncells=500, ass="SCT", groupvarname="group1",lfcvarname="log_fc", extended.output=T, deduped=T, nmarkers=NULL){
  
markerlist=get.cellinfo.markerlist(so, clusvar=clusvar, clusters=NULL, markers=markers,  groupvarname="group1",lfcvarname="log_fc", extended.output=T, deduped=T, nmarkers=nmarkers)
  
  markergaps=lapply(1:length(markerlist), function(x) rep(x, length(markerlist[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which
  
  fcat("Getting cells for each cluster...")
  cellist=lapply(as.character(clusids), function(x) {
    out=WhichCells(so, idents=x, downsample=ncells)
    #fcat("total number of cells:", length(out))
    out
  }
  
  ) %>% givename(., as.character(clusids))
  
  
  ordcellist=tryCatch({
    
    if (meth=="rankcounts"){
      celltotals=list()
      ordcellist=lapply(1:length(cellist), function(x) {
        
        newcells= metadata(so)[cellist[[x]], ] %>% select(!!sym(clusvar), nCount_RNA) %>% arrange(-nCount_RNA) %>% slice_head(n=ncells.show) %>% rownames
        celltotals[[x]]<<- rep(x, length(newcells))
        newcells
      }
      )
      
      
    }else{
      if (meth=="signature"){
        fcat("Using signature method to rank cells...")
        celltotals=list()
        ordcellist=lapply(1:length(cellist), function(x) {
          
          print(markerlist[[x]])
          
          sumcells=apply(so[markerlist[[x]],cellist[[x]]]@assays[[ass]]@counts, 2, sum)
          
          
          #arranged from weakest to strongest  
          newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
          celltotals[[x]]<<- rep(x, length(newcells))
          newcells
        }
        ) %>% givename(., clusids)
        
      }else{
        if (meth=="seriate"){
          library(seriation)
          celltotals=list()
          ordcellist=lapply(1:length(cellist), function(x) {
            
            
            ss=seriation::get_order(seriate(as.matrix(so[markerlist[[x]],cellist[[x]]]@assays$RNA@data), margin=2, method="BEA_TSP"))
            celltotals[[x]]<<- rep(x, length(ss))
            cellist[[x]][ss]
          }
          )  %>% givename(., clusids)
          
          
        }else{
          
          if(meth=="seurat"){
            fcat("Using Seurat Module Score method to rank cells...")
            celltotals=list()
            
          
            newnames=names(markerlist) %>% make.names
            so <- AddModuleScore(so, markerlist, name=newnames )
            
            
            ordcellist=lapply(1:length(clusids), function(x) {
              
              print(markerlist[[x]])
              
              modulename= paste0(newnames[x], x)
              
              #sumcells=apply(so[markerlist[[x]],so %>% metadata %>% pull(), 2, sum)
              
              newcells=so@meta.data[cellist[[x]], ] %>% filter(!!sym(clusvar)==clusids[x]) %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
              
              #arranged from weakest to strongest  
              #newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
              celltotals[[x]]<<- rep(x, length(newcells))
              newcells
            }
            ) %>% givename(., clusids)
            
            
          }else{
            
            if(meth=="seurat2"){
              fcat("Using Seurat Module Score method to rank cells and retrieveng strongest", ncells, "cells in clusters...")
              celltotals=list()
              
              so <- AddModuleScore(so, markerlist, name=names(markerlist) )
              
              
              ordcellist=lapply(1:length(clusids), function(x) {
                
                print(markerlist[[x]])
                
                modulename= paste0(clusids[x], x)
                
                #sumcells=apply(so[markerlist[[x]],so %>% metadata %>% pull(), 2, sum)
                newcells=so@meta.data %>% filter(!!sym(clusvar)==clusids[x]) %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
                newcells=newcells[1:ncells]
                
                #arranged from weakest to strongest  
                #newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
                celltotals[[x]]<<- rep(x, length(newcells))
                newcells
              }
              )
              
              
              
            }
            
          }
          
          
        }}}
    ordcellist} , error=function(e){fcat("There was trouble with seriating the cells. returning unseriated cells");
      print(e); 
      cellist
    })
  
  celltotals=lapply(1:length(ordcellist), function(x) rep(x, length(ordcellist[[x]])))
  
  fcat("making cell vector...")
  cellvector=Reduce(c, ordcellist)
  fcat("getting the gap positions...")
  cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which  
  
  cell.reference=lapply(1:length(ordcellist), function(nn) rep(names(ordcellist)[nn], length(ordcellist[[nn]])   )) %>% Reduce(c, .)  
  marker.reference=lapply(1:length(markerlist), function(nn) rep(names(markerlist)[nn], length(markerlist[[nn]])   )) %>% Reduce(c, .)
  
  cell.annotation= as.data.frame(cell.reference) %>% giverownames(cellvector) %>% givecolnames(., 1, clusvar)
  marker.annotation= as.data.frame(marker.reference) %>% giverownames(markervector) %>% givecolnames(., 1, clusvar)
  
  
  
  
  if(extended.output==T){
    list(cells=cellvector,
         cell.list=ordcellist,
         gaps.cells=cellgaps,
         gaps.markers=markergaps,
         markers=markervector,
         markerlist=markerlist,
         cell.annotation=cell.annotation,
         marker.annotation=marker.annotation,
         clustering.variable=clusvar,
         input.cluster.ids=clusids,
         method=meth,
         cell.metadata=clusvar,
         marker.metadata=clusvar
    )}else{
      cellvector
    }
}



seriatecells2=function(so, clusvar="seurat_clusters", clusters=NULL, meth="signature", markers=NULL, ncells=100, ass="SCT", groupvarname="group",lfcvarname="logFC", extended.output=F, deduped=T, groups.match.markers=T){
  
  Idents(so)=clusvar
  if(is.null(clusters)){
    clusids= so[[clusvar]] %>% pull(clusvar) %>% as.character %>% unique %>% sort
  }else{
    clusids=clusters 
  }
  markers.are.list=F
  markers.are.table=F
  ## check if markers have been provided externally
  if(!is.null(markers) & class(markers)=="list"){
    markers.are.list=T
    markers.are.table=F
    fcat("External marker list has been provided.")
  }else{
    markers.are.list=F
    if(!is.null(markers) & (class(markers) %in% c("data.table", "matrix","data.frame"))){
      fcat("External marker table has been provided.")
      markers.are.table=T
    }
  }
  ### CASE WHEN EXTERNAL LIST OF MARKERS IS PROVIDED
  if(markers.are.list){    
    markerlist0=markers
    
    if(!is.null(names(markerlist0))){# if markerlist has names
      clusids=names(markerlist0)
      
      #all(clusters %in% (ref.cellinfo$markerlist %>% names))
      
    }else{ ### if markerlist has no names just make some
      fcat("external markers have no names. Automatically generating names...")
      clusids=(1:length(markerlist0)) %>% paste0("C", .)
      names(markerlist0)= clusids
      
    }}
  
  # assemble the markers into a simple dataframe to perform duplication operations
  if(deduped==T){
    fcat("deduplicating marker list...")
    fullgns=lapply(1:length(markerlist0), function(x) markerlist0[[x]]) %>% Reduce(c, .)
    fullgns.groups=lapply(1:length(markerlist0), function(x) rep(names(markerlist0)[x], length(markerlist0[[x]])) ) %>% Reduce(c, .)
    
    # in the absence of GE information, remove duplicates in the order they appear first time marker appears is considered the 
    #most relevant. 
    dupd= duplicated(fullgns)
    markers.df=as.data.frame(fullgns.groups[!dupd]) %>% giverownames(fullgns[!dupd]) %>% givecolnames(., 1, nms=clusvarname) 
    
    markerlist= lapply(as.character(clusids), function(x)  markers.df %>% filter(!!sym(clusvarname)==!!x) %>% rownames)
    names(markerlist)=as.character(clusids)
  }
  
  ### NO EXTERNAL MARKERS PROVIDED: MARKER TABLE EXPECTED WITHIN OBJECT
  if(!markers.are.list & !markers.are.table & is.null(markers)){
    #if external markers are not provided at all, start from the top markers table
    markers= so@assays$RNA@misc$top_markers
    if(is.null(markers)){
      markers= so@assays$RNA@misc$markers
      
      if(is.null(markers)){
        #if top_markers is not there, then throw an error
        stop("There are no markers! please add markers under @assays$RNA@misc$top_markers\n" ) 
      }else{
        markers.are.table=T 
      }
      
    }else{
      markers.are.table=T
    }
    
    if(markers.are.table){
      # if there is a table, collect the DE markers for each group. 
      fcat("Splitting groups and retrieving top markers from table...")
      markerlist0=lapply(clusids, function(x) markers %>% filter(!!sym(groupvarname)==x))
      
      hasgenes= Vectorize(function(tabb) as.logical(nrow(tabb)>0), USE.NAMES=F)
      fcat("current length of markerlist0", length(markerlist0))
      which.have.genes=hasgenes(markerlist0)
      markerlist0=markerlist0[which.have.genes]
      fcat("current length of markerlist0 after filtering", length(markerlist0))
      fcat("current length of clusids", length(clusids))
      clusids=as.character(clusids)[which.have.genes]
      fcat("current length of clusids after filter", length(clusids))
      names(markerlist0)= clusids 
      
      
      fcat("Markers found for the following groups:", paste(clusids, collapse=","))
      
    } 
    
    #removing duplicates when a marker DE table is provided (not a list)
    if(deduped==T & markers.are.table){
      fullgns= lapply(markerlist0, function(x) x %>% pull(feature)) %>% Reduce(c, .)
      
      ###when duplicates, keep the gene as marker in the cluster where it shows highest log_fold change
      fcat("Assigning genes in more than one cluster to the highest expressing cluster...")    
      deduped=lapply(unique(fullgns), function(x){
        lapply(markerlist0, function(y) y %>% filter(feature==x)) %>% Reduce(condrbind, .) %>% arrange(-!!sym(lfcvarname)) %>% head(n=1)  %>% select(feature,!!sym(lfcvarname), !!sym(groupvarname)) 
        
      }) %>% Reduce(rbind, .)
      
      deduped2= deduped %>% dplyr::mutate(ncluster=sapply(deduped[[groupvarname]], function(x) generalreplace(x, clusids,1:length(clusids)), USE.NAMES=F )) %>% arrange(ncluster)
      
      #arrange genes from highest log fold change to lowest per cluster
      markers.revised= deduped2 %>% group_by(ncluster) %>% arrange(-!!sym(lfcvarname), .by_group=T) %>% ungroup()
      
      
    }
    
    
    markergaps= markers.revised %>% pull(ncluster) %>% as.integer %>% diff %>% as.logical %>% which
    #markervector=markers.revised %>% pull(feature) 
    
    markerlist=lapply(as.character(clusids), function(x) markers.revised %>% filter(!!sym(groupvarname)==x) %>% pull(feature))
    
    names(markerlist)=as.character(clusids)  
    
  }
  
  ### COMMON DOWNSTREAM PROCESSING AFTER DEFINING THE MARKER LIST 
  
  markervector= markerlist %>% Reduce(c, .)
  markergaps=lapply(1:length(markerlist), function(x) rep(x, length(markerlist[[x]]))) %>% Reduce(c, .) %>% diff %>% as.logical %>% which
  
  
  fcat("Getting cells for each cluster...")
  cellclasses=so %>% metadata %>% pull(clusvar) %>% unique
  fcat("cellclasses are", cellclasses)
  cellist=lapply(cellclasses, function(x) {
    
    #  if(!(x %in% ( cellclasses ))  ){
    #  NULL}else{
    out=WhichCells(so, idents=x, downsample=ncells)
    #fcat("total number of cells:", length(out))
    out
    #}
  }
  
  ) %>% givename(., cellclasses)
  
  #fcat("removing clusters for which no cells were found")
  #nonnulls=lapply(cellist, function(x) is.null(x)) %>% unlist %>% unname
  #cellist= cellist[!nonnulls] %>% givename(., as.character(clusids[!nonnulls]))
  
  
  if (meth=="rankcounts"){
    celltotals=list()
    ordcellist=lapply(1:length(cellist), function(x) {
      
      newcells= metadata(so)[cellist[[x]], ] %>% select(!!sym(clusvar), nCount_RNA) %>% arrange(-nCount_RNA) %>% slice_head(n=ncells.show) %>% rownames
      celltotals[[x]]<<- rep(x, length(newcells))
      newcells
    }
    )
    
    
  }else{
    if (meth=="signature"){
      fcat("Using signature method to rank cells...")
      celltotals=list()
      ordcellist=lapply(1:length(cellist), function(x) {
        
        print(markerlist[[x]])
        
        sumcells=apply(so[markerlist[[x]],cellist[[x]]]@assays[[ass]]@counts, 2, sum)
        
        
        #arranged from weakest to strongest  
        newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
        celltotals[[x]]<<- rep(x, length(newcells))
        newcells
      }
      )
      
      cellvector=Reduce(c, ordcellist)
      fcat("getting the gap positions...")
      cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which
    }else{
      if (meth=="seriate"){
        library(seriation)
        celltotals=list()
        ordcellist=lapply(1:length(cellist), function(x) {
          
          
          ss=seriation::get_order(seriate(as.matrix(so[markerlist[[x]],cellist[[x]]]@assays$RNA@data), margin=2, method="BEA_TSP"))
          celltotals[[x]]<<- rep(x, length(ss))
          cellist[[x]][ss]
        }
        ) 
        
        
      }else{
        
        if(meth=="seurat"){
          fcat("Using Seurat Module Score method to rank cells...")
          celltotals=list()
          
          newnames=names(markerlist) %>% make.names
          #names(newnames)=names(markerlist)
          fcat("newnames is ", newnames)
          so=tryCatch({ AddModuleScore(so, markerlist, name=newnames)}, error=function(e){
            fcat("trying 12 bins because of numerical issues...")
            AddModuleScore(so, markerlist, nbin=12, name=newnames)})
          
          
          if(groups.match.markers){
            ordcellist=lapply(1:length(markerlist), function(x) {
              
              print(names(markerlist)[x])
              
              modulename= paste0(newnames[x], x)
              markername=names(markerlist)[x]
              #sumcells=apply(so[markerlist[[x]],so %>% metadata %>% pull(), 2, sum)
              newcells=so@meta.data[cellist[[markername]], ] %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
              if(!is.null(ncells) & ncells<= length(newcells)){
                #get the ncells strongest cells
                newcells=newcells[(length(newcells)-ncells):length(newcells)]
              }
              #arranged from weakest to strongest  
              #newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
              celltotals[[x]]<<- rep(x, length(newcells))
              newcells
            }
            )
            
            cellvector=Reduce(c, ordcellist)
            fcat("getting the gap positions...")
            cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which
          }else{ # if cell groups do not match markers, we just retrieve the ncells with the strongest signature
            fcat("cells and markers do not have corresponding ids. Getting cells with strongest signature")
            ordcellist=lapply(1:length(clusids), function(x) { #we evaluate the strength of all signatures regardless of the cell label
              
              print(names(markerlist)[x])
              
              modulename= paste0(newnames[clusids[x]], x)
              
              #arrange from biggest to smallest
              newcells=so@meta.data %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
              
              if(!is.null(ncells) & ncells<= length(newcells)){
                #get the ncells strongest cells
                newcells=newcells[(length(newcells)-ncells):length(newcells)]
              }
              #arranged from weakest to strongest  
              #newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
              celltotals[[x]]<<- rep(x, length(newcells))
              newcells
            }
            )
            
            cellvector=Reduce(c, ordcellist)
            fcat("getting the gap positions...")
            cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which
            
          } 
        }else{ #another method different from above
          
          if(meth=="seurat2"){
            fcat("Using Seurat Module Score method to rank cells and retrieveng strongest", ncells, "cells in clusters...")
            celltotals=list()
            
            so <- AddModuleScore(so, markerlist, name=names(markerlist) )
            
            
            ordcellist=lapply(1:length(clusids), function(x) {
              
              print(markerlist[[x]])
              
              modulename= paste0(clusids[x], x)
              
              #sumcells=apply(so[markerlist[[x]],so %>% metadata %>% pull(), 2, sum)
              newcells=so@meta.data %>% filter(!!sym(clusvar)==clusids[x]) %>% select(modulename) %>% arrange(!!sym(modulename)) %>% rownames
              newcells=newcells[1:ncells]
              
              #arranged from weakest to strongest  
              #newcells= metadata(so)[cellist[[x]], ] %>% dplyr::mutate(markersum= sumcells, signature=markersum/nCount_RNA) %>% select(!!sym(clusvar), signature) %>% arrange(signature) %>% slice_head(n=ncells) %>% rownames
              celltotals[[x]]<<- rep(x, length(newcells))
              newcells
            }
            )
            
            cellvector=Reduce(c, ordcellist)
            fcat("getting the gap positions...")
            cellgaps= celltotals %>% Reduce(c,.) %>% diff %>% as.logical %>% which
            
            
          }
          
        }
        
        
      }
      
    }
  }
  if(extended.output==T){
    list(cells=cellvector,cell.list=ordcellist , gaps.cells=cellgaps, gaps.markers=markergaps, markers=markervector, markerlist=markerlist)}else{
      cellvector
    }
}

################################################################################
# classify cells into singlet, doublet or negative based on the hto_demux
################################################################################

classifydemux= Vectorize(function(x){
  if(x!="Negative" & x!= "Doublet" & x!= "Unassigned"& x!= "Blank" & x!= "Multiplet"){
    return("Singlet")
  }else{
    return(x)
  }
}, USE.NAMES=F)


################################################################################
# retrieve a doublet scoring filter for a dataset
################################################################################

getdbcutoff= Vectorize(function(nm){
  dbcutoffs[[nm]]
  
}, USE.NAMES=F)


################################################################################
# paste two strings together with a comma in between (great with Reduce())
################################################################################


pastec=function(a,b) paste0(c(a,b), collapse=",")


################################################################################
# a more comprehensive condition taking into account the demux state, negative or doublet. 
################################################################################

get_condition_demux=Vectorize(function(condition, demux){
  
  if(is.na(condition)){
    
    if(!(is.na(demux))){
      return(demux)
    }else{
      return("Blank")
      
    }
    
  }else{
    return(condition)
  }
  
}, USE.NAMES=F)


################################################################################
# remove the last character from a string
################################################################################
chopstring=function(st) substr(st, 1, nchar(st)-1)


################################################################################
#fix bad names in data frame columns
################################################################################

fixnames= function(x){nms= x %>% colnames %>% gsub("#", "", .) %>% make.names; colnames(x)=nms; fcat("new names are", nms); x}


################################################################################
#get the replicate
################################################################################
getrep= Vectorize(function(x) repmap[[x]], USE.NAMES=F)


################################################################################
# trim final one in the cell names
###############################################################################

trimone= Vectorize(function(x) strsplit(x, split="-")[[1]][1]) 


################################################################################
#
################################################################################

is.loaded= function(pkgname) !requireNamespace(pkgname, quietly = TRUE)

################################################################################
# make square legend
################################################################################

squarelegend = function(varr, size=1, shape=0, colorlist=allcolors){
  
  arglist=list()
  arglist[[varr]]= guide_legend(
    override.aes = list(
      shape = shape,  # 0: square, 1: circle
      size = size,   # Customize point sizes
      color = allcolors[[varr]],  # Customize point colors
      label = allcolors[[varr]] %>% names # Customize legend labels
    )
  )
  
  do.call(guides, arglist)
  
}

################################################################################
# compare category changes for cells across several category variables
################################################################################

transferplot= function(so=NULL, df=NULL,  labelvars, colorlist=allcolors, facetvar="all_cells", facet.ncol=NULL, facet.nrow=NULL, column.width=0.5, alpha=0.5, scale=2, width=3, height=2, plotlabel=""){
  
  
  sca=scale  
  freqvar=facetvar
  library(ggalluvial)
  library(ggrepel)
  
  #prepare colors if unavailable

  
  #lapply(labelvars, function(x) allcolors[[k]])    
  if(!is.null(df)){
    met=df %>% dplyr::mutate(all_cells="All cells")
  }else{
    if(!is.null(so)){  
      met= so %>% metadata %>% dplyr::mutate(all_cells="All cells")  
    }else{
      error("Please add either a Seurat object or a data frame as input!")
    }}
  
  
  fcat("preparing colors")
  for(k in 1:length(labelvars)){
    
    if(is.null(colorlist) || is.null(colorlist[[labelvars[k]]])){
      fcat(as.character(labelvars[k]))
      cats=met %>% pull(!!sym(as.character(labelvars[k]))) %>% unique
      fcat("unique categories:", paste(cats, collapse=" "))
      colorlist[[as.character(labelvars[k])]]<- randomcolors( cats %>% length) %>% givename(., cats)
      #allcolors[[labelvars[k]]]<<- colorlist[[labelvars[k]]]
      
    }
  }
  fcat("Colors used for plotting:")
  cat(code.for.list.of.vectors(colorlist[labelvars]))
  
  
   fcat("creating factors")
  
  # factors suspended for now
  for(k in 1:length(labelvars)){
    
    met = met %>% dplyr::mutate(!!labelvars[k]:=factor(!!sym(labelvars[k]), levels=colorlist[[labelvars[k]]] %>% names ) )
    
  }
 
    fcat("calculating frequencies t1")

  t1=met %>% dplyr::select(c(freqvar,labelvars)) %>% group_by(!!!syms(c(freqvar, labelvars))) %>% summarise( counts=n())  %>% group_by(!!sym(freqvar)) %>% dplyr::mutate(Freq=counts/sum(counts)) %>% ungroup() 
   fcat("calculating frequencies t2")
  t2= t1 %>% dplyr::mutate(ind= 1:nrow(t1)) %>% pivot_longer(., cols=labelvars, names_to="originalcolumn", values_to="cell.label")
  
   fcat("making plot")
  tpdf(path=config$plotpath, paste0("alluvialplot_labeltransfer_",plotlabel, "_", paste(labelvars, collapse="-")), width=pw*5/2*scale*width, height=pw*scale*height)
  ga=ggplot(t2,
            aes(y = Freq, x=factor(originalcolumn, levels=labelvars), stratum=cell.label, alluvium=ind, fill=cell.label))+ 
    geom_flow() +
    geom_stratum(width = column.width, alpha=alpha) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)))+#, min.segment.length=0.1) +
    scale_x_discrete(expand= c(.1, .1), labels = labelvars)+
    facet_wrap(~factor(get(freqvar), levels=(colorlist[[freqvar]] %>% names)), ncol=facet.ncol, nrow=facet.nrow)+
    scale_fill_manual(values=colorlist[labelvars] %>% Reduce(c,.))+
    theme_classic()
  print(ga+NoLegend())
  dev.off()
  
  ga+NoLegend()
}

################################################################################
# get highest loadings from pca
################################################################################

gethighloadings=function(so, comps=1:5, nfeatures=10, reduction="pca", reduction.key="PC_"){
  rk=reduction.key
  compnms= paste0(rk, comps)
  flist= lapply(compnms, function(x){
    
    list(positive=
           so@reductions[[reduction]]@feature.loadings %>% as.data.frame %>% select(!!sym(x)) %>% arrange(-!!sym(x)) %>% head(nfeatures) %>% rownames
         , negative=
           so@reductions[[reduction]]@feature.loadings %>% as.data.frame %>% select(!!sym(x)) %>% arrange(!!sym(x)) %>% head(nfeatures) %>% rownames
    ) 
    
  }) %>% givename(compnms)
  flist
}


################################################################################

################################################################################
code.for.list.of.vectors <- function(list_of_vectors) {
  code0 <- "list(\n"
  codelist=list()
  nms=names(list_of_vectors)
  # Loop through each vector in the list
  for (i in seq_along(list_of_vectors)) {
    
    if(is.null(nms)){
      vector_name <- paste0("vector", i)
    }else{
      vector_name <- nms[i] 
    }
    code <-  paste0("  ", vector_name, " = c(\n")
    
    # Loop through each element in the vector
    for (element_name in names(list_of_vectors[[i]])) {
      element_value <- toString(list_of_vectors[[i]][element_name])
      code <- paste0(code, paste0("   ", element_name, " = '", element_value, "',\n"))
    }
    
    code <- paste0(code, ")\n")
    codelist[[i]]=code
  }
  code2= codelist %>% Reduce(c, .) %>% paste(., collapse=",")
  
  code3 <- paste0(code0, code2, ")\n") 
  
  return(code3 %>% gsub(",\n)", ")", .))
}

# Example usage:
#vector1 <- c("name1" = "value1", "name2" = "value2")
#vector2 <- c("name3" = "value3", "name4" = "value4")
#list_of_vectors <- list(vector1, vector2)

#code_string <- code.for.list.of.vectors(allcolors[labelvars])
#code_string
#cat(code_string)

################################################################################
# remove a field from a list
################################################################################

dropfield= function(lisst, x){ lisst[[x]]=NULL; lisst} 


################################################################################
# create an arbitrary colormap
################################################################################

breakColors = function(breaks, colors, center=0, tol=0.001)
{
  ## In case of explicit color definitions
  nbreaks = length(breaks)
  nclass  = nbreaks - 1
  if (!is.function(colors)) {
    ncolors = length(colors)
    if (ncolors > nclass) {
      warning("more colors than classes: ignoring ", ncolors-nclass, " last colors")
      colors = colors[1:nclass]
    } else if (nclass > ncolors) {
      stop(nclass-ncolors, " more classes than colors defined")
    }
  } else {
    ## Are the classes symmetric and of same lengths?
    clens = diff(breaks)
    aclen = mean(clens)
    if (aclen==0) stop("Dude, your breaks are seriously fucked up!")
    relerr = max((clens-aclen)/aclen)
    if ( (center %in% breaks) & (relerr < tol) ) { ## yes, symmetric
      ndxcen = which(breaks==center)
      kneg = ndxcen -1 
      kpos = nbreaks - ndxcen
      kmax = max(kneg, kpos)
      colors = colors(2*kmax)
      if (kneg < kpos) {
        colors = colors[ (kpos-kneg+1) : (2*kmax) ]
      } else if (kneg > kpos) {
        colors = colors[ 1 : (2*kmax - (kneg-kpos)) ]
      }
    } else {                                      ## no, not symmetric
      colors = colors(nclass)
    }
  }
  colors
}

################################################################################

removegene=function(genes,x) genes[genes!=x]

################################################################################
# format recuciton key based on an arbitrary id.
################################################################################


prepare.rk=function(x)   gsub("\\.", "", x) %>% gsub("_", "", .) %>% paste0(., "_") 


################################################################################
#
################################################################################


sling_cell_df <- function(sling_ctr){
  raw_pseudotime <- slingPseudotime(sling_ctr)
  order_df <- apply(raw_pseudotime, 2, function(x){
    x[!is.na(x)] <- rank(x[!is.na(x)])
    x
  }) %>%
    as.data.frame()%>% 
    rename_with(.,~paste0("order","_",.)) %>% 
    dplyr::mutate(mean_order = apply(., 1, mean, na.rm = TRUE))
  
  
  pseudotime_df <- apply(raw_pseudotime, 2, function(x){
    x1 <- x[!is.na(x)]
    x[!is.na(x)] <- (x1-min(x1))/(max(x1)-min(x1))*100
    x
  })%>%
    as.data.frame()%>% 
    rename_with(.,~paste0("pseudotime","_",.))%>% 
    dplyr::mutate(mean_pseudotime = apply(., 1, mean, na.rm = TRUE))
  
  df <- bind_cols(raw_pseudotime,pseudotime_df, order_df) %>% 
    as.data.frame()
  
  return(df)
}


################################################################################
# ENTROPY FUNCTION
################################################################################

#https://rpubs.com/philjet/shannonentropy
#entropy function
entropy <- function(target) {
  freq <- table(target)/length(target)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  -sum(vec * log2(vec))
}



################################################################################
# findmin
################################################################################

findmin = function(arr){
  
  m= min(abs(arr))
  return(which(arr==m))
  
}


################################################################################
# get stage number from a "DX" format
################################################################################

get.stage.number= Vectorize(
  function(x) substr(x, 2, 10) %>% as.numeric
)



################################################################################
# dplyr - ready get rows
################################################################################


getrows=function(x, rows) x %>% names2col(., "nm") %>% filter(nm %in% rows) %>% select(-nm)


################################################################################
# generalreplace, grep version
################################################################################

#value may be longer than elements in reff and so may contain reff
generalreplace.grepl=function(vec, refs, targets, return.sums=T){
  subss=list()
  finalvector=vec
  for(n in 1:length(refs)){
    reff=refs[n]
    which.match=grepl(reff, vec)
    subss[[n]]=rep(FALSE, length(vec))
    subss[[n]][which.match]=TRUE
    finalvector=gsub(reff, targets[n], finalvector )
    
  }
  allsums=subss %>% Reduce('+', .)  
  if(any(allsums)>1){
    warning("Same value substituted multiple times")
  }
  if(return.sums==T){
    list(final.vector=finalvector, number.replacements=allsums)
  }else{
    finalvector 
  }
  
}
################################################################################
#
################################################################################

reorderdf= function(df, varr, catorder, return.list=T){
  
  l=lapply(catorder, function(x){
    df %>% filter(!!sym(varr)==x)
  })
  if(return.list){
  }else{
    l=l %>% Reduce(rbind,.)
  }
  l
}





################################################################################
# wrapper to take 2 seurat objects with a mapfun variable and show the points of one in another's clusters (groups)
################################################################################

#project.to.glasswork=function(ref, query, )
################################################################################
# functions that emulate matlab zeros and ones
################################################################################
zeros <- function(rows, cols) {
  return(matrix(0, nrow = rows, ncol = cols))
}

ones <- function(rows, cols) {
  return(matrix(1, nrow = rows, ncol = cols))
}


remove.zerorows= function(matt){

  if(length(dim(matt))<=1){
    print("Argument not a matrix. returning as is")
    matt
    }else{
    matt[apply(matt, 1, sum)!=0,]}

}
  
  remove.zerocols=function(matt){
      if(length(dim(matt))<=1){
    print("Argument not a matrix. returning as is")
    matt
        }else{
    matt[,apply(matt, 2, sum)!=0]}}

get.variant.rows= function(matt){
  
  #non.zero.rows=apply(matt, 1, sum)!=0
  variant.rows=apply(matt, 1, sd)!=0
  #print(cbind(non.zero.rows, variant.rows))
  variant.rows
  
}

get.invariant.rows= function(matt){
  
  #non.zero.rows=apply(matt, 1, sum)!=0
  invariant.rows=apply(matt, 1, sd)==0
  #print(cbind(non.zero.rows, variant.rows))
  invariant.rows
  
}

get.nonzerorows= function(matt) (apply(matt, 1, sum)!=0) %>% unname  


fill.missing.markers.matrix=function(matt, markernames, fill.value=0){
  missing.ones=NULL 
  
  newlines=lapply(1:length(markernames), function(k){ 
    
    nl=NULL
    if(!(markernames[k] %in% rownames(matt))){
      
      
      
      #missing.ones=c(missing.ones, markernames[k])
      
      nc=length(colnames(matt))
      nl=matrix(data=0, nrow=1, ncol= nc, dimnames=list(markernames[k], colnames(matt)))
    }
    nl
  }) %>% Reduce(rbind, .)
  #fcat("number of columns of matt", length(colnames(matt)), "number of columns of newlines", ncol(newlines))
  matt=rbind(matt, newlines)
  #if(!is.null(missing.ones)){
  
  #fcat("missing genes filled:", paste(missing.ones, collapse=",")) 
  #}
  matt
}



################################################################################
#for plotting purposes, slightly jitter one position of a row
################################################################################


random.glitch.row=function(mat, rownumber=x, val=0.00001){
  
  s=sample(1: ncol(mat), 1)
  mat[  rownumber,  s]=mat[  rownumber,  s]+val
  
  mat
  
}

################################################################################
# randomly glitch all invariant rows for plotting purposes
################################################################################

glitch.invariants=function(matt, val=0.00001){
  
  invariant.rows= !get.variant.rows(matt)   
  
  
  for(k in 1:length(invariant.rows)){
    set.seed(invariant.rows[k])
    matt=random.glitch.row(matt, rownumber=invariant.rows[k], val=val)
    
    
  }
  set.seed(42)
  matt
}

################################################################################
#create a heatmap based directly from a cellinfo object
################################################################################
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



################################################################################
#refresh cellinfo object
################################################################################

  ##############################################################################
  # in this version, no assumptions are made and all information is defaulted.
  # all previous annotations are thrown away
  ##############################################################################
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




## this version attempts to preserve previous annotations 
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

################################################################################
#update clustering variable
################################################################################
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

################################################################################
# seriate genes given a cellinfo object
################################################################################

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



################################################################################
# cellinfo related functions
################################################################################
adjust.markerlist= function(markerlist, n) lapply(markerlist, function(x) removenas(x[1:n])) %>% givename(., markerlist %>% names)


adjust.markerlist.tail= function(markerlist, n) lapply(markerlist, function(x) removenas(x[(length(x)-n):length(x)])) %>% givename(., markerlist %>% names)



removenas=function(x) x[!is.na(x)]
adjust.markers= function(cellinfo, n){
  cellinfo$markerlist= adjust.markerlist(cellinfo$markerlist, n)
  cellinfo=refresh.cellinfo(cellinfo)
  cellinfo
}


################################################################################
# adjust cells (get top)
################################################################################

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


################################################################################
# cellinfo collect values for a variable inside the cellinfo
################################################################################

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


################################################################################
# cellinfo single dataset workflow
################################################################################


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

################################################################################
# cellinfo- create ranked cell plot... single variable
################################################################################

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




################################################################################
#get markers
################################################################################

get.markers=function(so){
  
  so$RNA@misc$markers 
}

get.top.markers=function(so){
  
  so$RNA@misc$top_markers 
}

get.significant.markers=function(so){
  
  so$RNA@misc$significant_markers
}


################################################################################
# format gene vector for enrichr
################################################################################

format.enrichr=function(x) x %>% paste(., collapse="\n") %>% cat


################################################################################
# function to subset a cellinfo to a few subgroups.
################################################################################
subset.cellinfo=function(cellinfo, labels){
  fcat("restricting to labels only present in list")
  tf=lapply(labels, function(x) x %in% (cellinfo$markerlist %>% names)) %>% Reduce(c, .)
  labels2=labels[tf]
  
  cellinfo$markerlist=cellinfo$markerlist[labels2]
  cellinfo$cell.list=cellinfo$cell.list[labels2]
  cellinfo$input.cluster.ids=labels2
  refresh.cellinfo(cellinfo)
}

fillmat=function(genemat, xx){ 
  if(is.null(genemat[[xx]]))
  {genemat[[xx]]=0}
  
  return(genemat)}

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

################################################################################
# remove genes from a cell info object
################################################################################
removemarkers=function(cellinfo, markers){
  nms=names(cellinfo$markerlist)
  cellinfo$markerlist=lapply(1:length(cellinfo$markerlist), function(x){
    gns=cellinfo$markerlist[[x]]
    
    gns[!(gns %in% markers)]
    
  }) %>% givename(., nms)
  
  refresh.cellinfo(cellinfo)
}


################################################################################

#szymkiewicz_simpson_coefficient to calculate size of overlaps between vectors
zsoverlap <- function(vector1, vector2) {
  
  # Calculate the size of the intersection
  intersection_size <- length(intersect(vector1 , vector2))
  #fcat(intersection_size)
  # Calculate the minimum size of the two sets
  if(length(vector1)==0 || length(vector1)==0){
    fcat("Warning: some of the sets are zero length")
  }
  min_size <- min(length(vector1), length(vector2))
  #fcat(min_size)
  # Calculate the Szymkiewicz–Simpson coefficient
  szymkiewicz_simpson <- intersection_size / min_size
  
  return(szymkiewicz_simpson)
}


################################################################################
# filter cells based on a prediciton score for a mapping.
################################################################################

filterps= function(so, varbl="seurat_clusters", th) so[, so %>% metadata %>% filter(!!sym(paste0(mapfun(varbl), "_predicted.id.score"))>=!!th ) %>% rownames]



################################################################################
# selectively remove heatmap legends from a ComlpexHeatmap object
################################################################################


removelegends= function(hm, legends, where="top"){
  
  annoslot=paste0(where, "_annotation")
  #for(ann in c("top_annotation", "right_annotation", "left_annotation", "bottom_annotation")){
  for(leg in legends){
    slot(hm, annoslot)@anno_list[[leg]]@show_legend=F
  }
  
  hm 
  #}
}


################################################################################
#   classify cells into tf, ligands and receptors
################################################################################


simpleCache(paste_("genes", "db", "human_tfs"), {
# m <- fread("http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv")
hocomoco.tfs= c("AHR", "AIRE", "ALX1", "ALX3", "ALX4", "ANDR", "ANDR", "ANDR", 
                "AP2A", "AP2B", "AP2C", "AP2D", "ARI3A", "ARI5B", "ARNT2", "ARNT", 
                "ARX", "ASCL1", "ASCL2", "ATF1", "ATF2", "ATF2", "ATF2", "ATF3", 
                "ATF4", "ATF6A", "ATF7", "ATOH1", "BACH1", "BACH2", "BARH1", 
                "BARH2", "BARX1", "BARX2", "BATF3", "BATF", "BATF", "BC11A", 
                "BCL6B", "BCL6", "BHA15", "BHE22", "BHE23", "BHE40", "BHE41", 
                "BPTF", "BRAC", "BRAC", "BRCA1", "BSH", "CDC5L", "CDX1", "CDX2", 
                "CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG", "CEBPZ", "CENPB", 
                "CLOCK", "COE1", "COT1", "COT1", "COT2", "COT2", "CPEB1", "CR3L1", 
                "CR3L2", "CREB1", "CREB3", "CREB5", "CREM", "CRX", "CTCFL", "CTCF", 
                "CUX1", "CUX2", "CXXC1", "DBP", "DDIT3", "DLX1", "DLX2", "DLX3", 
                "DLX4", "DLX5", "DLX6", "DMBX1", "DMRT1", "DPRX", "DRGX", "DUX4", 
                "DUXA", "E2F1", "E2F2", "E2F3", "E2F4", "E2F4", "E2F5", "E2F6", 
                "E2F7", "E2F8", "E4F1", "EGR1", "EGR2", "EGR2", "EGR3", "EGR4", 
                "EHF", "ELF1", "ELF2", "ELF3", "ELF5", "ELK1", "ELK3", "ELK4", 
                "EMX1", "EMX2", "EOMES", "EPAS1", "ERG", "ERR1", "ERR2", "ERR3", 
                "ESR1", "ESR1", "ESR2", "ESR2", "ESX1", "ETS1", "ETS2", "ETV1", 
                "ETV2", "ETV3", "ETV4", "ETV5", "ETV6", "ETV7", "EVI1", "FEV", 
                "FEZF1", "FIGLA", "FLI1", "FLI1", "FOSB", "FOSL1", "FOSL2", "FOS", 
                "FOXA1", "FOXA2", "FOXA3", "FOXB1", "FOXC1", "FOXC2", "FOXD1", 
                "FOXD2", "FOXD3", "FOXF1", "FOXF2", "FOXG1", "FOXH1", "FOXI1", 
                "FOXJ2", "FOXJ3", "FOXJ3", "FOXK1", "FOXL1", "FOXM1", "FOXO1", 
                "FOXO3", "FOXO4", "FOXO6", "FOXP1", "FOXP2", "FOXP3", "FOXQ1", 
                "FUBP1", "GABPA", "GATA1", "GATA1", "GATA2", "GATA2", "GATA3", 
                "GATA4", "GATA5", "GATA6", "GBX1", "GBX2", "GCM1", "GCM2", "GCR", 
                "GCR", "GFI1B", "GFI1", "GLI1", "GLI2", "GLI3", "GLIS1", "GLIS2", 
                "GLIS3", "GMEB2", "GRHL1", "GRHL2", "GSC2", "GSC", "GSX1", "GSX2", 
                "HAND1", "HAND1", "HBP1", "HEN1", "HES1", "HES5", "HES7", "HESX1", 
                "HEY1", "HEY2", "HIC1", "HIC2", "HIF1A", "HINFP", "HLF", "HLTF", 
                "HMBX1", "HME1", "HME2", "HMGA1", "HMGA2", "HMX1", "HMX2", "HMX3", 
                "HNF1A", "HNF1B", "HNF1B", "HNF4A", "HNF4G", "HNF6", "HOMEZ", 
                "HSF1", "HSF1", "HSF2", "HSF4", "HTF4", "HXA10", "HXA11", "HXA13", 
                "HXA1", "HXA2", "HXA5", "HXA7", "HXA9", "HXB13", "HXB1", "HXB2", 
                "HXB3", "HXB4", "HXB6", "HXB7", "HXB8", "HXC10", "HXC11", "HXC12", 
                "HXC13", "HXC6", "HXC8", "HXC9", "HXD10", "HXD11", "HXD12", "HXD13", 
                "HXD3", "HXD4", "HXD8", "HXD9", "ID4", "IKZF1", "INSM1", "IRF1", 
                "IRF2", "IRF3", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "IRX2", 
                "IRX3", "ISL1", "ISL2", "ISX", "ITF2", "JDP2", "JUNB", "JUND", 
                "JUN", "KAISO", "KAISO", "KAISO", "KLF12", "KLF13", "KLF14", 
                "KLF15", "KLF16", "KLF1", "KLF3", "KLF4", "KLF5", "KLF6", "KLF8", 
                "KLF9", "LBX2", "LEF1", "LHX2", "LHX3", "LHX4", "LHX6", "LHX8", 
                "LHX9", "LMX1A", "LMX1B", "LYL1", "MAFA", "MAFB", "MAFF", "MAFF", 
                "MAFG", "MAFG", "MAFK", "MAFK", "MAF", "MAF", "MAX", "MAZ", "MAZ", 
                "MBD2", "MCR", "MECP2", "MEF2A", "MEF2B", "MEF2C", "MEF2D", "MEIS1", 
                "MEIS1", "MEIS2", "MEIS3", "MEOX1", "MEOX2", "MESP1", "MGAP", 
                "MITF", "MIXL1", "MLXPL", "MLX", "MNX1", "MSX1", "MSX2", "MTF1", 
                "MXI1", "MXI1", "MYBA", "MYBB", "MYB", "MYCN", "MYC", "MYF6", 
                "MYNN", "MYOD1", "MYOD1", "MYOG", "MZF1", "NANOG", "NANOG", "NDF1", 
                "NDF2", "NF2L1", "NF2L2", "NFAC1", "NFAC1", "NFAC2", "NFAC3", 
                "NFAC4", "NFAT5", "NFE2", "NFIA", "NFIA", "NFIB", "NFIC", "NFIC", 
                "NFIL3", "NFKB1", "NFKB2", "NFYA", "NFYB", "NFYC", "NGN2", "NKX21", 
                "NKX22", "NKX23", "NKX25", "NKX28", "NKX31", "NKX32", "NKX61", 
                "NKX61", "NKX62", "NOBOX", "NOTO", "NR0B1", "NR1D1", "NR1D1", 
                "NR1H2", "NR1H3", "NR1H3", "NR1H4", "NR1H4", "NR1I2", "NR1I2", 
                "NR1I3", "NR1I3", "NR2C1", "NR2C2", "NR2E1", "NR2E3", "NR2F6", 
                "NR4A1", "NR4A2", "NR4A3", "NR5A2", "NR6A1", "NRF1", "NRL", "OLIG1", 
                "OLIG2", "OLIG2", "OLIG3", "ONEC2", "ONEC3", "OSR2", "OTX1", 
                "OTX2", "OVOL1", "OVOL2", "OZF", "P53", "P53", "P5F1B", "P63", 
                "P63", "P73", "P73", "PATZ1", "PATZ1", "PAX1", "PAX2", "PAX3", 
                "PAX4", "PAX5", "PAX6", "PAX7", "PAX8", "PBX1", "PBX1", "PBX2", 
                "PBX3", "PBX3", "PDX1", "PDX1", "PEBB", "PHX2A", "PHX2B", "PIT1", 
                "PITX1", "PITX2", "PITX3", "PKNX1", "PLAG1", "PLAL1", "PO2F1", 
                "PO2F2", "PO2F3", "PO3F1", "PO3F2", "PO3F3", "PO3F4", "PO4F1", 
                "PO4F2", "PO4F3", "PO5F1", "PO5F1", "PO6F1", "PO6F2", "PPARA", 
                "PPARA", "PPARD", "PPARG", "PPARG", "PRD14", "PRDM1", "PRDM4", 
                "PRDM6", "PRGR", "PRGR", "PROP1", "PROX1", "PRRX1", "PRRX2", 
                "PTF1A", "PTF1A", "PURA", "RARA", "RARA", "RARA", "RARB", "RARG", 
                "RARG", "RARG", "RAX2", "RELB", "REL", "REST", "RFX1", "RFX1", 
                "RFX2", "RFX2", "RFX3", "RFX4", "RFX5", "RFX5", "RHXF1", "RORA", 
                "RORG", "RREB1", "RUNX1", "RUNX2", "RUNX3", "RXRA", "RXRA", "RXRB", 
                "RXRG", "RX", "SALL4", "SCRT1", "SCRT2", "SHOX2", "SHOX", "SIX1", 
                "SIX2", "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMCA1", "SMCA5", 
                "SNAI1", "SNAI2", "SOX10", "SOX10", "SOX11", "SOX13", "SOX15", 
                "SOX17", "SOX18", "SOX1", "SOX21", "SOX2", "SOX2", "SOX3", "SOX4", 
                "SOX5", "SOX7", "SOX8", "SOX9", "SOX9", "SP1", "SP1", "SP2", 
                "SP2", "SP3", "SP4", "SP4", "SPDEF", "SPI1", "SPIB", "SPIC", 
                "SPZ1", "SRBP1", "SRBP2", "SRF", "SRY", "STA5A", "STA5B", "STAT1", 
                "STAT1", "STAT2", "STAT3", "STAT4", "STAT6", "STF1", "SUH", "TAF1", 
                "TAL1", "TAL1", "TBP", "TBR1", "TBX15", "TBX19", "TBX1", "TBX20", 
                "TBX21", "TBX2", "TBX3", "TBX4", "TBX5", "TCF7", "TEAD1", "TEAD2", 
                "TEAD3", "TEAD4", "TEF", "TF2LX", "TF65", "TF7L1", "TF7L2", "TFAP4", 
                "TFCP2", "TFDP1", "TFE2", "TFE3", "TFEB", "TGIF1", "TGIF2", "THA11", 
                "THAP1", "THA", "THA", "THB", "THB", "TLX1", "TWST1", "TWST1", 
                "TYY1", "TYY2", "UBIP1", "UNC4", "USF1", "USF2", "VAX1", "VAX2", 
                "VDR", "VDR", "VENTX", "VEZF1", "VEZF1", "VSX1", "VSX2", "WT1", 
                "WT1", "XBP1", "Z324A", "Z354A", "ZBED1", "ZBT14", "ZBT17", "ZBT18", 
                "ZBT48", "ZBT49", "ZBT7A", "ZBT7B", "ZBTB4", "ZBTB4", "ZBTB6", 
                "ZEB1", "ZEP1", "ZEP2", "ZF64A", "ZFHX3", "ZFP28", "ZFP42", "ZFP82", 
                "ZFX", "ZFX", "ZIC1", "ZIC2", "ZIC3", "ZIC4", "ZIM3", "ZKSC1", 
                "ZKSC3", "ZN121", "ZN134", "ZN134", "ZN136", "ZN140", "ZN143", 
                "ZN148", "ZN214", "ZN219", "ZN232", "ZN250", "ZN257", "ZN260", 
                "ZN263", "ZN263", "ZN264", "ZN274", "ZN281", "ZN282", "ZN317", 
                "ZN320", "ZN322", "ZN329", "ZN331", "ZN333", "ZN335", "ZN335", 
                "ZN341", "ZN341", "ZN350", "ZN350", "ZN382", "ZN384", "ZN394", 
                "ZN394", "ZN410", "ZN418", "ZN418", "ZN423", "ZN436", "ZN449", 
                "ZN467", "ZN490", "ZN502", "ZN524", "ZN528", "ZN547", "ZN549", 
                "ZN554", "ZN554", "ZN563", "ZN563", "ZN582", "ZN586", "ZN589", 
                "ZN652", "ZN667", "ZN680", "ZN708", "ZN708", "ZN713", "ZN740", 
                "ZN768", "ZN770", "ZN770", "ZN784", "ZN816", "ZN816", "ZNF18", 
                "ZNF41", "ZNF41", "ZNF76", "ZNF85", "ZNF85", "ZNF8", "ZSC16", 
                "ZSC22", "ZSC31", "ZSCA4")
#union(ccbr.tfs, hocomoco.tfs)  
hocomoco.tfs
}, assignToVar="tfList",  recreate=F)

fcat("1")
istf= function(x) x %in% tfList 
#lrtable=readRDS(file.path(params$resourcepath, "cell_cell_interactions/CellTalkDB/human_lr_pair.rds"))

#isligand=function(x) x %in%  (lrtable %>% pull(ligand_gene_symbol))
#isreceptor=function(x) x %in%  (lrtable %>% pull(receptor_gene_symbol))
#islrpair=function(x) x %in%  (lrtable %>% pull(lr_pair))


colorgene=function(x) {
  
  if(istf(x)){return("black")}else{
    if(isligand(x)){return("blue")}else{
      if(isreceptor(x) || (x %in% c("CRABP1", "RARRES1", "RARRES2", "CRABP2"))){return("red")}else{return("grey")}
    }
  }}

colorgene2=function(x) {
  
  if(istf(x)){return("black")}else{return("grey")}}





#alosteric protein database
#https://mdl.shsmu.edu.cn/ASD2023Common/static_file/archive_2023/ASD_Release_202309_AS.tar.gz

################################################################################
#
################################################################################

divide_legend_columns <- function(plot, num_columns, legend.position="bottom") {
  # Get the current legend
  current_legend <- plot + theme(legend.position = legend.position)  # Adjust the position if needed
  
  # Modify the legend to have x columns
  modified_legend <- current_legend +
    guides(color = guide_legend(nrow = 1, ncol = num_columns))
  
  return(modified_legend)
}



################################################################################
# No Legend title
################################################################################
NoLegendTitle <- function(plot) {
  # Modify the legend to remove the title
  guides(guide_legend(title = NULL))
  
}
################################################################################

remove_x_axis <- function() {
  # Modify the theme to remove x-axis
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
}

################################################################################
# scale text size
################################################################################

scale.text=function(scaling_factor){
  theme(
    text = element_text(size =  scaling_factor),
    axis.title = element_text(size =  scaling_factor),
    axis.text = element_text(size = scaling_factor),
    axis.title.x = element_text(size = scaling_factor),
    axis.text.x = element_text(size =  scaling_factor),
    axis.title.y = element_text(size =  scaling_factor),
    axis.text.y = element_text(size =  scaling_factor),
    legend.text = element_text(size =  scaling_factor),
    legend.title = element_text(size =  scaling_factor),
    plot.title = element_text(size =  scaling_factor),
    plot.subtitle = element_text(size =  scaling_factor),
    plot.caption = element_text(size =  scaling_factor)
  )
}



scale.text.general=function(plt, sca){
  plt$theme$text$size= plt$theme$text$size*sca
  plt
}


################################################################################
# masking types in kameneva et al cell types
################################################################################
masktype=Vectorize(function(x) if(x =="SCP" || x=="sympathoblasts" || x=="mesenchyme"){return(x)}else{return("other")}, USE.NAMES=F)

masktype3=Vectorize(function(x) if(x =="SCP" || x=="sympathoblasts" || x=="mesenchyme" || x=="chromaffin"){return(x)}else{return("other")}, USE.NAMES=F)

################################################################################
#replace markers within a seurat object
################################################################################

replacemarkers=function(so, newmarkers){
  if(!is.null(newmarkers)){
    so$RNA@misc$top_markers=newmarkers 
    return(so)
  }else{
    return(so)
  }}


rmquery=Vectorize(function(x){if(grepl("query",x)){
  strsplit(x, split="query_")[[1]][2]
}else{x} }, USE.NAMES=F)



rmquery.matrix= function(x){
  x %>% names2col(., "cellid") %>% dplyr::mutate(cellid=rmquery(cellid)) %>% col2names(., "cellid")
}




################################################################################
#  remove duplicates from a vector
################################################################################

deduplicate= function(x) x[!duplicated(x)]


################################################################################
# quick print a matrix's contents
################################################################################

quickprint=function(mat, wide=F, limit=5) if(wide){ mat[1:limit, ]}else{mat[1:limit, 1:limit]}


################################################################################
# fit hill function 
################################################################################

hillfunction <- function(x, k, n, beta) beta*((x**n)/((k**n)+(x**n)))

#model1 <- nls(fluorI ~ eDecay(t,myA,myT), data=ExpData, start=list(myA=10,myT=5))



################################################################################
# gglabel (auxiliary side bar plot to make a legend
################################################################################

manual.legend=function(colorlist, pname, size=7, trailing.spaces=15, bar.width=.6){
  
  
  p=colorlist[[pname]] %>% rev
  
  pl=ggplot(1:length(p) %>% as.data.frame %>% dplyr::mutate(h=.1, colors=p, name=names(p)), aes(x=factor(.), y=factor(h)))+
    geom_col(aes(fill=factor(.)))+
    scale_fill_manual(values=namenums(p))+
    coord_cartesian( ylim=c(0,10))+
    theme_classic()+
    #geom_text(y=0.5, color="black", angle=90)+
    NoAxes()+
    NoLegend()+
    ggtitle(paste0(paste(rep(" ", trailing.spaces), collapse=""), pname))
  
  return(pl+geom_text(inherit.aes=T,aes(label=name), y=1.05, size=size, color="black", hjust=0)+coord_flip(ylim=c(0,1/bar.width)))
}

################################################################################
# print an extract of a matrix without getting into huge columns or huge rows overprinting trouble
################################################################################

quickview=function(mat, wide=F, limit=5) if(wide){ mat[1:limit, ]}else{mat[1:limit, 1:limit]}

################################################################################
# roc analysis to quantify the signal of targets within selected peaks
################################################################################
make.roc.curves=function(targets, targetpeaks, pca.data, classvar="celltype"){


####make a function that takes the values for a peak for all samples and calculates a roc curve


is.target= function(x) ifelse(x %in% targets, TRUE, FALSE)

fcat("collecting embeddings...")
#ae= get.all.embeddings(selected.pc)
ae= pca.data %>% dplyr::mutate(class=is.target(!!sym(classvar))) %>% select(class, !!sym(classvar))
sp= get.sample.peakvals(peakid=targetpeaks[, "Geneid"], samples= rownames(ae)) %>% cbind(ae, .)


################################################################################
# generate roc.curves based on the individual promoter-embedded peaks to dplyr::select those that on its own do a great job of dplyr::selecting target cells. 
################################################################################
fcat("starting roc curve loop...")
#selected.peaks=peak.gene.results[[selected.pc]][[region.selection]]$Geneid[1:total.intervals] # dplyr::select peaks with overlapping promoters
#selected.names=peak.gene.results[[selected.pc]][[region.selection]]$gene_symbol[1:total.intervals] # dplyr::select peaks with overlapping promoters
fcat("Applying ROC to intervals")
nl=0
allrocs=lapply(1:nrow(targetpeaks), function(xx){
pp=targetpeaks[xx, "Geneid"]

#trying new function to make roc calculation more efficient
#fcat("starting to calculate roc matrix")
rocmat=get.roc.matrix(sp, signal.col=pp, class.col="class")
#fcat("done calculating roc")
#rocmat=roc.curve(sp[, c(pp, "class")],step=0.001, true.value=TRUE, false.value=FALSE, color="black", line.width=1, plot.symbol=20,make.plot=F )
#dev.off()
au=trapezium.auc2(rocmat)
#fcat(pp, nn, au)
if(nl==50){
cat("\n")
  nl<<-0
}else{
  nl<<-nl+1
cat(".")
}

list(
  #plot=ggplot(rocmat %>% as.data.frame)+geom_path(aes(x=FPs, y=TPs))+ggtitle(nn)+geom_abline(slope=1, intercept=0, color="red", linetype="dashed")+theme_classic()+NoAxes(),
     auc=au, 
     Geneid=pp,
     gene_symbol=nn)



}) 
allrocs
}

################################################################################
# roc analysis to quantify the signal of a collective list of targets within selected peaks
################################################################################

make.roc.curves.all=function(targetlist, targetpeaks, pca.data){


####make a function that takes the values for a peak for all samples and calculates a roc curve


is.target= function(x, targets) fifelse(x %in% targets, TRUE, FALSE)

fcat("collecting embeddings...")
#ae= get.all.embeddings(selected.pc)

sp= get.sample.peakvals(peakid=targetpeaks[, "Geneid"], samples= rownames(pca.data)) %>% cbind(., pca.data %>% select(!contains("PC")))



################################################################################
# generate roc.curves based on the individual promoter-embedded peaks to dplyr::select those that on its own do a great job of dplyr::selecting target cells. 
################################################################################
fcat("starting roc curve loop...")
#selected.peaks=peak.gene.results[[selected.pc]][[region.selection]]$Geneid[1:total.intervals] # dplyr::select peaks with overlapping promoters
#selected.names=peak.gene.results[[selected.pc]][[region.selection]]$gene_symbol[1:total.intervals] # dplyr::select peaks with overlapping promoters
fcat("Applying ROC to intervals")
nl=0
allrocs=lapply(1:nrow(targetpeaks), function(xx){
  if(nl==0){
    cat(xx, ".", targetpeaks[xx, "Geneid"])
  }
alltargets=lapply(1:length(targetlist), function(xy){

classvar=targetlist[[xy]]$classvar
spt =sp %>% dplyr::mutate(class=is.target(!!sym(classvar), targetlist[[xy]]$targets))   


pp=targetpeaks[xx, "Geneid"]
#fcat(pp)
#trying new function to make roc calculation more efficient
#fcat("starting to calculate roc matrix")
rocmat=get.roc.matrix(spt, signal.col=pp, class.col="class")
#fcat("done calculating roc")
#rocmat=roc.curve(sp[, c(pp, "class")],step=0.001, true.value=TRUE, false.value=FALSE, color="black", line.width=1, plot.symbol=20,make.plot=F )
#dev.off()
au=trapezium.auc2(rocmat)
#fcat(pp, nn, au)


list(
  #plot=ggplot(rocmat %>% as.data.frame)+geom_path(aes(x=FPs, y=TPs))+ggtitle(nn)+geom_abline(slope=1, intercept=0, color="red", linetype="dashed")+theme_classic()+NoAxes(),
     auc=au, 
     Geneid=pp,
     gene_symbol=nn, 
     target=names(targetlist)[xy]
     )

}) %>% bind_rows
  
if(nl==50){
cat("\n")
  nl<<-0
}else{
  nl<<-nl+1
cat(".")
}
  alltargets
}) %>% bind_rows
allrocs
}

################################################################################
# load GSEA gene signatures
################################################################################

load.gene.signatures=function(sigs=NULL, recreate=F){
  
library(hypeR)
library(stringr)
#library(signatureSearch)
    #updated chemical perturbations signatures
    simpleCache::simpleCache("all_gene_signatures" %>% addversion, {
        #cp.sigs= read_gmt(getExternalFile("l1000_cp.gmt", "https://lincs-dcic.s3.amazonaws.com/LINCS-sigs-2021/gmt/l1000_cp.gmt"))
  
 # lig.sigs= read_gmt(getExternalFile("l1000_lig.gmt","https://lincs-dcic.s3.amazonaws.com/LINCS-sigs-2021/gmt/l1000_lig.gmt"))

      
    gsetpars=list(
        GOBP=list(set=msigdb_gsets(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"), name="GOBP"),
    CELLTYPES=list(set=msigdb_gsets(species="Homo sapiens", category="C8"), name="CELLTYPES"),
    PANGLAO=list(set=hypeR::enrichr_gsets(genesets = "PanglaoDB_Augmented_2021"), name="PANGLAODB"),
    DESCARTES=list(set=hypeR::enrichr_gsets(genesets = "Descartes_Cell_Types_and_Tissue_2021"), name="DESCARTES"),
    WIKIPATHWAYS=list(set=hypeR::enrichr_gsets(genesets ="WikiPathway_2021_Human"), name="WIKIPATHWAYS"),
    GOMF=list(set=msigdb_gsets(species="Homo sapiens", category="C5", subcategory="GO:MF")), 
    #GPDOWN=list(set=hypeR::enrichr_gsets(genesets = "Gene_Perturbations_from_GEO_down"), name="GPDOWN"), 
    #GPUP=list(set=hypeR::enrichr_gsets(genesets = "Gene_Perturbations_from_GEO_up"), name="GPUP"),
      TFPERT=list(set=hypeR::enrichr_gsets(genesets = "TF_Perturbations_Followed_by_Expression"), name="TFPERT"),
      #=list(set=hypeR::enrichr_gsets(genesets = "Kinase_Perturbations_from_GEO_down"), name="KINPERTDOWN"),
      #KINPERTUP=list(set=hypeR::enrichr_gsets(genesets = "Kinase_Perturbations_from_GEO_up"), name="KINPERTUP"),  
      #CHEMPERTDOWN=list(set=hypeR::gsets$new(cp.sigs[grepl("down$", names(cp.sigs))], name="LINCS1000_chempert_down_downloaded", version="v2021"), name="CHEMPERTDOWN"),
      #CHEMPERTUP=list(set=hypeR::gsets$new(cp.sigs[grepl("up$", names(cp.sigs))], name="LINCS1000_chempert_up_downloaded", version="v2021"), name="CHEMPERTUP"),
      #CHEMPERTDOWN=list(set=hypeR::enrichr_gsets(genesets = "LINCS_L1000_Chem_Pert_down"), name="CHEMPERTDOWN"),
      #CHEMPERTUP=list(set=hypeR::enrichr_gsets(genesets = "LINCS_L1000_Chem_Pert_up"), name="CHEMPERTUP"),
      LIGPERTDOWN=list(set=hypeR::enrichr_gsets(genesets = "Ligand_Perturbations_from_GEO_down"), name="LIGPERTDOWN"),
      LIGPERTUP=list(set=hypeR::enrichr_gsets(genesets = "Ligand_Perturbations_from_GEO_up"), name="LIGPERTUP")
      #LIGPERTLINCSDN=list(set=hypeR::gsets$new(lig.sigs[grepl("down$", names(lig.sigs))], name="LINCS1000_ligandpert_down_downloaded", version="v2021"), name="LIGPERTLINCSDN"),
      #LIGPERTLINCSUP=list(set=hypeR::gsets$new(lig.sigs[grepl("up$", names(lig.sigs))], name="LINCS1000_ligandpert_up_downloaded", version="v2021"), name="LIGPERTLINCSUP"),
      #CGP=list(set=msigdb_gsets(species="Homo sapiens",category="C2", subcategory = "CGP" ), name="CGP")
    #ARCHS4CELL_LINES=list(set=ARCHCELLLINES, name="ARCHS4_Cell-Lines")
    
  )
    }, assignToVar="gsetpars", recreate=recreate)
  
    fcat("Existent gene sets are:", paste(names(gsetpars), collapse=","))
if(is.null(sigs)){
  gsetpars
  }else{
  gsetpars[sigs]
  } 

}

################################################################################
# perform hypeR gene set enrichments and plots
################################################################################

get.hyper.enrichments= function(markerlist, gsetpars=load.gene.signatures(), test="hypergeometric", statt="JACCARD", markersid="",plotting.format="spaced5", colorlist=NULL, numbars=5, title="", make.barplot=F, make.heatmap=F, background=23467){
  library(hypeR)
    library(stringr)
    if(is.null(colorlist)){
  colorlist=randomcolors(markerlist %>% names %>% length) %>% givename(.,markerlist %>% names )
  fcat(dput(colorlist))
    }
allhyperdfs=lapply(gsetpars, function(gp){
    hyperwtdf=NULL
    gset=gp$set
    gsetname=gp$name
    
    get.gsgenes= Vectorize(function(genesetname) gset$genesets[[genesetname]], USE.NAMES=F)
    
    
    
    #tryCatch({
    hyperwtdf=lapply(1:length(markerlist), function(y){
                       cat(gsetname, "gene set",names(markerlist)[y], "\n")
                       
                       
                       x=markerlist[[y]] 
                       if(length(x)==0){
                         return()}
                       hyp=hypeR::hypeR(signature=x, genesets=gset, test=test, plotting=TRUE, background=background)
                       
                       if (statt=="JACCARD"){
                         df=hyp$data %>% dplyr::mutate(jaccard_index= overlap/(signature+geneset-overlap)) %>% dplyr::arrange(-jaccard_index) 
                       }else{
                         if(statt=="PVAL"){
                           df=hyp$data %>% dplyr::mutate(jaccard_index= overlap/(signature+geneset-overlap)) %>% dplyr::arrange(pval) 
                           
                         }    else {
                           
                           if(statt=="FDR"){
                             df=hyp$data %>% dplyr::mutate(jaccard_index= overlap/(signature+geneset-overlap)) %>% dplyr::arrange(fdr) 
                             
                           }
                         
                         }
                         
                       }
                       df =df %>% dplyr::mutate(geneset.size=geneset) 
                       df =df %>% dplyr::mutate(geneset=gsetname) 
                       df$clusterlabel=names(markerlist)[y]
                       #df$cluster=y+seq(0, 1-(1/numbars), 1/numbars)
                       df$guide=1
                       df
                     }) %>% Reduce(rbind,.)
    
    if(!is.null(hyperwtdf)){
      hyperwtdf[, "label"]= gsub("_", " ", hyperwtdf[, "label"]) %>% str_to_sentence
      hyperwtdf= hyperwtdf %>% names2col(., "original.geneset")
      pparsets=list(
        tight=list(barw=0.98, tsz=200, scf=0.3, sca=1.3, tallf=8,fatf=5, parstring="tight"),
        spaced=list(barw=0.95, tsz=180, scf=0.3, sca=1.3, tallf=12,fatf=5, parstring="spaced"),
        spaced2=list(barw=0.95, tsz=180, scf=0.3, sca=1.3, tallf=12,fatf=15, parstring="superwide"),
        spaced3=list(barw=0.90, tsz=160, scf=0.3, sca=1.3, tallf=18,fatf=15, parstring="superwide"), 
        spaced4=list(barw=0.90, tsz=80, scf=0.1, sca=8, tallf=25,fatf=15, parstring="fit"),
        spaced5=list(barw=0.90, tsz=20, scf=0.1, sca=8, tallf=25,fatf=15, parstring="fit")
      )
      parsel=pparsets[[plotting.format]]
      fatf=parsel$fatf
      tallf=parsel$tallf
      tsz=parsel$tsz
      gtsz=parsel$tsz*parsel$scf
      barw=parsel$barw
      sca=parsel$sca
      parstring=parsel$parstring
      allstops=hyperwtdf$cluster
      
      ### add some formatting to the labels
      
      hyperwtdf=hyperwtdf %>% dplyr::mutate(original.label=label, short.label=summarisename(label, 30), label= str_replace(label, paste0(str_to_sentence(gsetname)," "),""))
      
      ##we use pdf because sometimes the text in the categories is very long and goes beyond the edge of the artboards which we unclip in illustrator. 
      # all genes in marker list
    }
    
  hyperwtdf
})


  allhyperdfs=allhyperdfs %>% bind_rows %>% dplyr::filter(overlap>0)
  allhyperdfs=allhyperdfs %>% dplyr::mutate(pval.direct=dhyper(overlap,geneset.size, background-geneset.size, signature)) %>% group_by(clusterlabel, geneset) %>% dplyr::mutate(padj=p.adjust(pval, method="BH"), tlog.pval=ifelse(pval.direct>=0.05, 0, -log(pval.direct)), tlog.padj=ifelse(padj<=0.05,-log(padj), 0)) %>% ungroup() %>% as.data.frame %>% dplyr::select(short.label, geneset, clusterlabel, overlap, jaccard_index, background, geneset.size, signature, hits,pval,pval.direct, padj, fdr, tlog.pval, tlog.padj,original.label)
}


  

hyper.analysis=function(so, cellinf, sigs=c("DESCARTES", "CELLTYPES", "PANGLAO", "TFPERT"), hyper.res=NULL, npergroup=5, label=""){
gsetpars=load.gene.signatures(sigs = sigs)
if(is.null(hyper.res)){
simpleCache(paste0("hyper_enrichments_id_",label, Project(so), "_groupvar_",cellinf$marker.group.variable, "sigs", paste(sigs,collapse="_" ))  %>% addversion, {
hyper.res=get.hyper.enrichments(cellinf$markerlist, gsetpars = gsetpars)
hyper.res
}, assignToVar="hyper.res", reload=T)
}
  
clusterlevs=allcolors[[cellinf$marker.group.variable]] %>% names
allgenesets=hyper.res$geneset %>% unique %>% removenas
hplots=lapply(allgenesets, function(gs){
  #first select all the significant short labels
sgs= hyper.res %>% filter(geneset==gs, padj<=0.05) %>% pull(short.label) %>% unique
#then filter the entire table based on those labels to catch those gene sets that are a bit enriched in other groups despite not being significant

hyper.res.gs=hyper.res %>% dplyr::filter(short.label %in% sgs) %>% group_split(factor(clusterlabel, levels=clusterlevs)) %>% lapply(., function(x) head(x, npergroup)) %>% Reduce(rbind, .) %>% select(clusterlabel, short.label,original.label, jaccard_index, hits, tlog.padj) %>% as.data.frame %>% dplyr::mutate(short.label=summarisename(original.label, total.print.length=50, first.words=6))

#hyper.res.gs=hyper.res %>% filter(geneset==gs, padj<=0.05) %>% group_split(factor(clusterlabel, levels=clusterlevs)) %>% lapply(., function(x) head(x, npergroup)) %>% Reduce(rbind, .) %>% select(clusterlabel, short.label,original.label, jaccard_index, hits, tlog.padj) %>% as.data.frame %>% dplyr::mutate(short.label=summarisename(original.label, total.print.length=50, first.words=6))

cats=hyper.res.gs$short.label %>% unique %>% rev
ppl<-ggplot(hyper.res.gs, aes(x=factor(short.label, levels=cats), y=factor(clusterlabel, levels=clusterlevs ), size=jaccard_index, color=tlog.padj))+geom_point()+theme_minimal()+rotatex(90)+scale_color_viridis()+theme(axis.text.x = element_text(hjust = 1, vjust=0.5))+ggtitle(paste(Project(so), gs, "enrichments per" ,cellinf$cell.group.variable))+labs(y=cellinf$cell.group.variable, x="gene set")+coord_flip()

tpdf(paste0("pointplot_cluster_geneset_enrichments_",label,"_", Project(so),"_", cellinf$marker.group.variable,"_", gs), 
#wi=pw*sca*.4,he=pw*sca*.2)
wi=pw*sca*.2,he=pw*sca*.05*npergroup)
print(ppl)
dev.off()
ppl
})
hplots

}







################################################################################
# extract colors from a ggplot, or any other parameter for that matter
################################################################################
getcolors=function(plt, type){
plot_data <- ggplot_build(plt)
  
  # Extract the colors from the data
  colors <- unique(unlist(lapply(plot_data$data, function(layer) {
    layer[[type]]
  })))
  colors
}

################################################################################
# find if an array is empty
################################################################################
isempty= Vectorize(function(x) length(x)==0, USE.NAMES=F)

################################################################################
# coerce a one dimensional named vector into a data frame
################################################################################

vector.as.df=function(mt, nm="value"){
 mt2=list(peakid=unname(mt))
 mt2=mt2[[1]] %>% data.frame %>% giverownames(., names(mt)) %>% givecolnames(., nms=nm)
 mt2
}


################################################################################
# add a random set of colors to colorlist if nonexistent
################################################################################


 create.colors.if.empty=function(colorlist, colorvar, dff, seed=42){
   set.seed(seed)
if(is.null(colorlist[[colorvar]])){
  caats= dff %>% pull(colorvar) %>% unique
 colorlist[[colorvar]]= randomcolors(length(caats)) %>% givename(., caats)
}
  set.seed(42)
    return(colorlist)
 }

################################################################################
# import metadata from dataset number j as described in the config
################################################################################

get.metadata.from.dataset=function(j, ids=NULL, exclude= dataset.info$exclude_cell_type, celltype.var=dataset.info$cell_type_variable, dataset.info=dataset.info){  
  
metadata.main<- fread(dataset.info$dataset_paths_metadata[j]) %>% as.data.frame %>% fixcolnames %>% dplyr::mutate(dsname=dataset.info$dataset_name[j], dscategory=ifelse(dataset.info$is_reference[j], "reference", "test"), dsdatatype="ATAC", Experiment=!!sym(dataset.info$metadata_sample_name_col[j])) %>% col2names(., dataset.info$metadata_sample_name_col[j])

if(!dataset.info$is_reference[j]){
  
 metadata.main= metadata.main %>% dplyr::mutate(!!sym(dataset.info$cell_type_variable):=paste0("test_",!!sym(dataset.info$cell_group_vars[j]) ) )
}
  
rownames(metadata.main)=splitsamplename(rownames(metadata.main))

if(is.null(ids)){
  ids=rownames(metadata.main) 
}
#retrieve quality control files

if(any(grepl("multiqc", list.files(dataset.info$dataset_paths_nfcore[j])))){
  
metadata.qc= fread(file=paste0( dataset.info$dataset_paths_nfcore[j], "/multiqc/", config$peak_width, "/multiqc_data/multiqc_samtools_stats_1.txt")) %>% as.data.frame 
rownames(metadata.qc) <- metadata.qc %>% pull(Sample) %>% splitsamplename


metadata.frip=import.multiqc.x.x.info(j, filename="multiqc_mlib_frip_score-plot.txt", varname="frip.score")
metadata.peak.count=import.multiqc.x.x.info(j, filename="multiqc_mlib_peak_count-plot.txt", varname="peak.count")

#bind all matrices together

metadata.complete <- bind_cols(metadata.main[ids,],metadata.qc[ids,], metadata.frip[ids, ,drop=F],metadata.peak.count[ids, ,drop=F] ) %>% giverownames(., nms=ids)


}else{
 warning(paste0("dataset.", dataset.info$dataset_name[j], ": Please verify the path containing the multiqc output. \nProceeding without global quality control measures"))
  metadata.frip <- NULL
  metadata.qc <- NULL
  metadata.peak.count<-NULL
  metadata.complete=metadata.main
}


metalist=list(
  metadata.main=metadata.main,
  metadata.frip=metadata.frip,
  metadata.qc=metadata.qc,
  metadata.complete=metadata.complete,
  metadata.peak.count=metadata.peak.count)

if(!is.null(exclude)){
metalist=metadata.exclude.groups(metalist, is.reference=dataset.info$is_reference[j], celltype.var=celltype.var, exclude= exclude)  
  
}
metalist

}


################################################################################
# exclude dataset groups from metadata
################################################################################
metadata.exclude.groups=function(metalist, is.reference, celltype.var=dataset.info$cell_type_variable, exclude= dataset.info$exclude_cell_type){
  
excluded.ids=NULL
if(is.reference){
    excluded.ids<-metalist$metadata.complete %>% filter(!!sym(celltype.var) %in% exclude) %>% rownames
 metadata.filtered<- metalist$metadata.complete %>% filter(!(!!sym(celltype.var) %in% exclude))
 included.ids=setdiff(rownames(metalist$metadata.complete), excluded.ids)

}else{
 metadata.filtered<-metalist$metadata.complete
 included.ids=rownames(metalist$metadata.complete)
  
   
}

list(metadata.complete=metadata.filtered,
  metadata.main=metalist$metadata.main[included.ids,],
  metadata.qc=metalist$metadata.qc[included.ids,], 
  metadata.frip=metalist$metadata.frip[included.ids,,drop=F],
  metadata.peak.count=metalist$metadata.peak.count[included.ids,,drop=F],
  excluded.ids=excluded.ids, 
  included.ids=included.ids)
}

################################################################################
# exclude datasets below a specific frip score
################################################################################
metadata.filter.frip=function(metalist, frip.threshold=0.05){
  
excluded.ids=NULL

    excluded.ids<-metalist$metadata.complete %>% filter(frip.score<frip.threshold) %>% rownames
 metadata.filtered<- metalist$metadata.complete %>% filter(frip.score>=frip.threshold)
 included.ids=setdiff(rownames(metalist$metadata.complete), excluded.ids)



list(metadata.complete=metadata.filtered,
  metadata.main=metalist$metadata.main[included.ids,],
  metadata.qc=metalist$metadata.qc[included.ids,], 
  metadata.frip=metalist$metadata.frip[included.ids,,drop=F],
  metadata.peak.count=metalist$metadata.peak.count[included.ids, , drop=F],
  excluded.ids=metalist$excluded.ids,
  excluded.ids.frip=excluded.ids,
  included.ids=included.ids)
}

################################################################################
# expand colors based on an updated dataframe that presumably describes more categories than those present in the colorlist
################################################################################

expand.colors=function(colorlist, df, varr){ 
  if(varr %in% colnames(df) && length(unique(df[,varr])>0)){
  extracats=setdiff(as.character(unique(df[, varr])), names(colorlist[[varr]])) %>% removenas() ;
  if(length(extracats)>0){
  colorlist[[varr]]=c(colorlist[[varr]], randomcolors(length(extracats)) %>% givename(., extracats)); 
  }
  }
  colorlist 
  
  }

update.all.colors=function(colorlist, df){
for(j in names(colorlist)){
colorlist=expand.colors(colorlist, df, j)
}
  colorlist
}


################################################################################
# import multiqc x,x matrix information
################################################################################


import.multiqc.x.x.info= function(j, filename="multiqc_mlib_peak_count-plot.txt", varname="peak.count"){
  
  
  ff= fread(file=paste0( dataset.info$dataset_paths_nfcore[j], "/multiqc/", config$peak_width, "/multiqc_data/", filename)) %>% as.matrix 
ff <- ff %>% giverownames(., ff[, "Sample"])
rownames(ff) <- rownames(ff) %>% splitsamplename
colnames(ff) <- colnames(ff) %>% splitsamplename
ff <- lapply(rownames(ff), function(x){l=list(); l[[x]]=as.numeric(ff[x,x]); l}) %>% unlist %>% as.data.frame %>% givecolnames(., nms=varname) 
ff

}
  
  
  
################################################################################
# import counts
################################################################################
  
   
  import.counts=function(j, dataset.info){ 
    simpleCache(paste0("count_list_dataset", dataset.info$dataset_name[j], "_path_", digest::digest(c(dataset.info$dataset_paths_counts_self[j],dataset.info$dataset_paths_counts_ref[j]) ) ) %>% addversion, {
if(dataset.info$is_reference[j] & sum(dataset.info$is_reference)==1){
  dataset.is.ref=1
counts.to.ref= fread(file= dataset.info$dataset_paths_counts_self[j], skip=1) %>% as.data.frame %>% col2names(., "Geneid")  
}else{
#when there is more than one reference, it will be assumed that each subreference has its peaks
#quantified in the global reference via the nfcore pipeline and the path is to the corresponding allmergedcounts_withreference 
counts.to.ref= fread(file=paste0( dataset.info$dataset_paths_counts_ref[j], "/allmergedcounts_withreference.txt"), skip=1) %>% as.data.frame %>% col2names(., "Geneid")
}

counts.to.self= fread(file=dataset.info$dataset_paths_counts_self[j], skip=1) %>% as.data.frame %>% col2names(., "Geneid")


colnames(counts.to.self)= counts.to.self %>% colnames %>% splitsamplename
colnames(counts.to.ref)= counts.to.ref %>% colnames %>% splitsamplename

list(counts.to.self=counts.to.self, counts.to.ref=counts.to.ref)

    }, assignToVar="counts.list", reload=T)
 counts.list     
      
  }



################################################################################
# gene peak and tfbs annotation routines
################################################################################

prepare.motif.annotation.materials=function(){
  fcat("importing gtf...")
geneAnnot=rtracklayer::import(dataset.info$gtfpath)

#get the annotation of only protein coding genes. perhaps a bit restrictive for this pipeline
selGeneAnnot <- geneAnnot %>% as.data.frame %>% dplyr::filter(type=="gene", gene_type=="protein_coding") %>% GRanges

#imported from ncnb2_code/atac/04,  doesn't work 
#selGeneAnnot <- geneAnnot[(as.character(seqnames(geneAnnot)) %in% (peaksDt[,unique(chrom)]) & geneAnnot$type=="gene" & geneAnnot$gene_type=="protein_coding"),]

## motif annotation files
  fcat("importing motifs frome HOCOMOCO v11...")
motifFile <- getExternalFile("HOCOMOCOv11_full_HUMAN_mono_jaspar_format.txt", "https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_jaspar_format.txt") 
motifFileAnnot <- getExternalFile("HOCOMOCOv11_full_annotation_HUMAN_mono.tsv", "https://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv") 
#### from 04 annotate_genes



#### from annotate motifs
  fcat("preparing motif scanning auziliar variables...")
simpleCache::simpleCache("hocomoco_v11", {
	hocomoco <- universalmotif::read_jaspar(motifFile)
	hocomoco <- XMatrixList(lapply(hocomoco, function(m) universalmotif::convert_motifs(m, class = "TFBSTools-PFMatrix")), use.names = F, type = "PFMatrixList",    matrixClass = "PFMatrix")
	names(hocomoco) <- lget(hocomoco, "name", slot=TRUE)
	return(hocomoco)
}, assignToVar="mtfSet")

hocoAnnot <- fread(motifFileAnnot)
setkey(hocoAnnot, Model)

print(mtfSet)
print(head(hocoAnnot))

# T should be TBXT, the others are to be dplyr::filtered out (both on Y chromosome)
hocoAnnot[`Transcription factor`=="T", `Transcription factor`:="TBXT"]

mtfSet <- mtfSet[hocoAnnot[names(mtfSet), `Transcription factor`] %in% geneAnnot$gene_name]

mtfNames <- hocoAnnot[names(mtfSet), `Transcription factor`] 
names(mtfNames) <- names(mtfSet)
mtfLabels <- sprintf("%s (%s)", mtfNames, names(mtfNames))
names(mtfLabels) <- names(mtfSet)

relevant.chromosomes= c(paste0("chr", 1:22), "chrX")


fcat("Variables ready.")
list(mtfSet=mtfSet, 
  mtfNames=mtfNames, 
  mtfLabels=mtfLabels, 
  selGeneAnnot=selGeneAnnot,
  relevant.chromosomes=relevant.chromosomes)
  
  
  
}

################################################################################
#
################################################################################

get.global.peak.annotations=function(){
  
  library(TFBSTools)
library(motifmatchr)
  
  mtf.prep.list=prepare.motif.annotation.materials()

  list2env(mtf.prep.list, .GlobalEnv)
  ################################################################################
# arrange the peaks with the highest PC loadings into a GRanges object
################################################################################
  
  fcat("formatting peak annotation as GRanges...")
simpleCache("peaks_dt_global_annotations" %>% addversion, {
peaks_dt_global=norm.peak.vals %>% as.data.frame %>% dplyr::mutate(seqnames=Chr,start=Start, end=End) %>% dplyr::select(Geneid,seqnames, start, end) %>% filter(seqnames %in% !!relevant.chromosomes) %>% GRanges
#peaks_dt_global=norm.peak.vals %>% as.data.frame %>% dplyr::mutate(seqnames=Chr,start=Start, end=End) %>% dplyr::select(Geneid,seqnames, start, end) %>% GRanges



fcat("Annotating TSS and overlapping...")
peak.gene.annotations=list(
peaks_dt_global=peaks_dt_global ,  
closest_tss_global=geneAssignmentStrategies$closest_tss(peaks_dt_global, gene_annot=selGeneAnnot, window_size=MAX_DIST_GENE),
    
overlapping_promoters_global=geneAssignmentStrategies$promo(peaks_dt_global, gene_annot=selGeneAnnot)
  )
peak.gene.annotations
}, assignToVar="peak.gene.annotations", reload=T)

 list2env(peak.gene.annotations, .GlobalEnv)
 
simpleCache(paste0("matched_motifs_allpeaks_version", analysis.version), {
#mtfMat <- versionedCache(paste0("motif_counts_", dplyr::selected_pc), instruction={
    fcat("matching motifs...")
	if(!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=T)){
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update=FALSE, ask=FALSE)
	}
  
	  mtfs <- motifmatchr::matchMotifs(mtfSet, peaks_dt_global, genome=config$genome_build, out="scores") 
	
	    fcat("counting motifs...")
	  mtf_mat <- as.matrix(motifmatchr::motifCounts(mtfs))
	rownames(mtf_mat) <- peaks_dt_global$Geneid
	 
	 fcat("formatting output...")
	trimtfname= Vectorize(function(x) strsplit(x, split="_HUMAN")[[1]][1], USE.NAMES=F)
	
matrix.tfnames= trimtfname(colnames(mtf_mat))	
mtf_mat=mtf_mat %>% givecolnames(., nms=matrix.tfnames)
	
mtf_mat}, assignToVar="mtf_mat", reload=T)
	#return(mtfMat)
#}, buildEnvir=list(
 # mtf_set=mtfSet, 
  #peaks_dt=peaksDt[,.(chrom,start,end,rid)],
  #genome_build=config$genome_build
#))

simpleCache("mtf_mat_gathered" %>% addversion, {
  
  unique.tfs=unique(colnames(mtf_mat))

mtf.mat.summed.list=lapply(unique.tfs, function(tf){
   fcat(tf)
  mtchs=colnames(mtf_mat)==tf
  if(sum(mtchs)==1){
    mtf_mat[,colnames(mtf_mat)==tf] %>% coerce.to.matrix
  }else{
     apply(mtf_mat[,colnames(mtf_mat)==tf], 1, sum) %>% coerce.to.matrix %>% giverownames(., nms=rownames(mtf_mat)) %>% givecolnames(.,nms= tf)
  }

}) 

mtf_mat_gathered=mtf.mat.summed.list %>% Reduce(cbind, .) %>% givecolnames(., nms=unique.tfs)
mtf_mat_gathered

}, assignToVar="mtf_mat_gathered", reload=T)



simpleCache(paste_("motif_matrix_calculations_matrixid", digest::digest(mtf_mat_gathered)) %>% addversion, {
  
    list(mtf_props=apply(mtf_mat_gathered>0,2, function(x) sum(x)/length(x)),
mtf_sums=apply(mtf_mat_gathered>0,2, function(x) sum(x)),
nonmtf_sums=apply(mtf_mat_gathered==0,2, function(x) sum(x)))
  }, assignToVar="mtf.props.list", reload=T)


  fcat("Done.")
list(peak.gene.annotations=peak.gene.annotations, 
  mtf_mat=mtf_mat, 
  mtf_mat_gathered=mtf_mat_gathered,
  mtf.props.list=mtf.props.list
  )


}



################################################################################
# functions to  obtain embeddings
################################################################################

get.all.embeddings= function(pc) pca.data %>% dplyr::arrange(!!sym(pc)) %>% dplyr::select(!!sym(pc), !!sym(classvar))
get.sample.peakvals= function(peakid, samples=NULL){
  
  if(!is.null(samples)){ mt=norm.peak.vals.t[ samples, peakid]}else{ mt=norm.peak.vals.t[ !colnames(norm.peak.vals) %in% c("Geneid", "Chr", "Start", "End", "Strand", "Length"),peakid ]}

if(length(peakid)==1){
  mt2=vector.as.df(mt, peakid)
}else{
rn=rownames(mt)  
mt2=apply(mt, 2, as.numeric)
rownames(mt2)=rn
}
rm(mt)
gc()
mt2
}


################################################################################
# functions to make roc curves
################################################################################

make.roc.curves=function(targets, targetpeaks, pca.data, classvar="celltype"){


####make a function that takes the values for a peak for all samples and calculates a roc curve


is.target= function(x) ifelse(x %in% targets, TRUE, FALSE)

fcat("collecting embeddings...")
#ae= get.all.embeddings(selected.pc)
ae= pca.data %>% dplyr::mutate(class=is.target(!!sym(classvar))) %>% select(class, !!sym(classvar))
sp= get.sample.peakvals(peakid=targetpeaks[, "Geneid"], samples= rownames(ae)) %>% cbind(ae, .)


################################################################################
# generate roc.curves based on the individual promoter-embedded peaks to dplyr::select those that on its own do a great job of dplyr::selecting target cells. 
################################################################################
fcat("starting roc curve loop...")
#selected.peaks=peak.gene.results[[selected.pc]][[region.selection]]$Geneid[1:total.intervals] # dplyr::select peaks with overlapping promoters
#selected.names=peak.gene.results[[selected.pc]][[region.selection]]$gene_symbol[1:total.intervals] # dplyr::select peaks with overlapping promoters
fcat("Applying ROC to intervals")
nl=0
allrocs=lapply(1:nrow(targetpeaks), function(xx){
pp=targetpeaks[xx, "Geneid"]

#trying new function to make roc calculation more efficient
#fcat("starting to calculate roc matrix")
rocmat=get.roc.matrix(sp, signal.col=pp, class.col="class")
#fcat("done calculating roc")
#rocmat=roc.curve(sp[, c(pp, "class")],step=0.001, true.value=TRUE, false.value=FALSE, color="black", line.width=1, plot.symbol=20,make.plot=F )
#dev.off()
au=trapezium.auc2(rocmat)
#fcat(pp, nn, au)
if(nl==50){
cat("\n")
  nl<<-0
}else{
  nl<<-nl+1
cat(".")
}

list(
  #plot=ggplot(rocmat %>% as.data.frame)+geom_path(aes(x=FPs, y=TPs))+ggtitle(nn)+geom_abline(slope=1, intercept=0, color="red", linetype="dashed")+theme_classic()+NoAxes(),
     auc=au, 
     Geneid=pp,
     gene_symbol=nn)



}) 
allrocs
}

################################################################################
# functions to make roc curves for all targets
################################################################################

make.roc.curves.all=function(targetlist, targetpeaks, pca.data){


####make a function that takes the values for a peak for all samples and calculates a roc curve

is.target= function(x, targets) fifelse(x %in% targets, TRUE, FALSE)

fcat("collecting embeddings...")
#ae= get.all.embeddings(selected.pc)

sp= get.sample.peakvals(peakid=targetpeaks[, "Geneid"], samples= rownames(pca.data)) %>% cbind(., pca.data %>% select(!contains("PC")))



################################################################################
# generate roc.curves based on the individual promoter-embedded peaks to dplyr::select those that on its own do a great job of dplyr::selecting target cells. 
################################################################################
fcat("starting roc curve loop...")
#selected.peaks=peak.gene.results[[selected.pc]][[region.selection]]$Geneid[1:total.intervals] # dplyr::select peaks with overlapping promoters
#selected.names=peak.gene.results[[selected.pc]][[region.selection]]$gene_symbol[1:total.intervals] # dplyr::select peaks with overlapping promoters
fcat("Applying ROC to intervals")
nl=0
allrocs=lapply(1:nrow(targetpeaks), function(xx){
  #if(nl==0){
  #  cat(xx, ".", targetpeaks[xx, "Geneid"])
  #}
alltargets=lapply(1:length(targetlist), function(xy){

classvar=targetlist[[xy]]$classvar
spt =sp %>% dplyr::mutate(class=is.target(!!sym(classvar), targetlist[[xy]]$targets))


pp=targetpeaks[xx, "Geneid"]
#fcat(pp)
#trying new function to make roc calculation more efficient
#fcat("starting to calculate roc matrix")
rocmat=get.roc.matrix(spt, signal.col=pp, class.col="class")
#fcat("done calculating roc")
#rocmat=roc.curve(sp[, c(pp, "class")],step=0.001, true.value=TRUE, false.value=FALSE, color="black", line.width=1, plot.symbol=20,make.plot=F )
#dev.off()
au=trapezium.auc2(rocmat)
#fcat(pp, nn, au)


list(
  #plot=ggplot(rocmat %>% as.data.frame)+geom_path(aes(x=FPs, y=TPs))+ggtitle(nn)+geom_abline(slope=1, intercept=0, color="red", linetype="dashed")+theme_classic()+NoAxes(),
     auc=au, 
     Geneid=pp,
     #gene_symbol=nn, 
     target=names(targetlist)[xy]
     )

}) %>% bind_rows
  
#if(nl==50){
#cat("\n")
#  nl<<-0
#}else{
#  nl<<-nl+1
#cat(".")
#}
  alltargets
}) %>% bind_rows
allrocs
}




make.roc.curves.fast=function(sample.peak.values, pca.data){


####make a function that takes the values for a peak for all samples and calculates a roc curve

is.target= function(x, targets) fifelse(x %in% targets, TRUE, FALSE)

fcat("collecting embeddings...")
#ae= get.all.embeddings(selected.pc)

sp= get.sample.peakvals(peakid=targetpeaks[, "Geneid"], samples= rownames(pca.data)) %>% cbind(., pca.data %>% select(!contains("PC")))



################################################################################
# generate roc.curves based on the individual promoter-embedded peaks to dplyr::select those that on its own do a great job of dplyr::selecting target cells. 
################################################################################
fcat("starting roc curve loop...")
#selected.peaks=peak.gene.results[[selected.pc]][[region.selection]]$Geneid[1:total.intervals] # dplyr::select peaks with overlapping promoters
#selected.names=peak.gene.results[[selected.pc]][[region.selection]]$gene_symbol[1:total.intervals] # dplyr::select peaks with overlapping promoters
fcat("Applying ROC to intervals")
nl=0
allrocs=lapply(1:nrow(targetpeaks), function(xx){
  #if(nl==0){
  #  cat(xx, ".", targetpeaks[xx, "Geneid"])
  #}
alltargets=lapply(1:length(targetlist), function(xy){

classvar=targetlist[[xy]]$classvar
spt =sp %>% dplyr::mutate(class=is.target(!!sym(classvar), targetlist[[xy]]$targets))


pp=targetpeaks[xx, "Geneid"]
#fcat(pp)
#trying new function to make roc calculation more efficient
#fcat("starting to calculate roc matrix")
rocmat=get.roc.matrix(spt, signal.col=pp, class.col="class")
#fcat("done calculating roc")
#rocmat=roc.curve(sp[, c(pp, "class")],step=0.001, true.value=TRUE, false.value=FALSE, color="black", line.width=1, plot.symbol=20,make.plot=F )
#dev.off()
au=trapezium.auc2(rocmat)
#fcat(pp, nn, au)


list(
  #plot=ggplot(rocmat %>% as.data.frame)+geom_path(aes(x=FPs, y=TPs))+ggtitle(nn)+geom_abline(slope=1, intercept=0, color="red", linetype="dashed")+theme_classic()+NoAxes(),
     auc=au, 
     Geneid=pp,
     #gene_symbol=nn, 
     target=names(targetlist)[xy]
     )

}) %>% bind_rows
  
#if(nl==50){
#cat("\n")
#  nl<<-0
#}else{
#  nl<<-nl+1
#cat(".")
#}
  alltargets
}) %>% bind_rows
allrocs
}


################################################################################
# coerce to matrix (similar to drop)
################################################################################


coerce.to.matrix=function(vec){
 if(is.vector(vec)){ 
  matrix(vec, ncol = 1, dimnames = list(names(vec), NULL))
 }else{
  vec 
 }
}

################################################################################

################################################################################
isempty.list= function(ll) lapply(ll , function(x) ifelse(length(x)<=1, T, F)) %>% unlist


################################################################################
# get gene set genes
################################################################################



################################################################################
# making enrichment heatmap
################################################################################
make.enrichment.heatmap=function(x, hyperdf, gsetpars=load.gene.signatures(), setid=NULL, cluster=F){

  
     if(is.null(setid)){
       
       setid=hyperdf %>% pull(clusterlabel) %>% unique %>% paste(., collapse="-")
     }
  
      fct=hyperdf$original.geneset %>% unique
      
      get.gsgenes= Vectorize(function(genesetname, gsetgroup) gsetpars[[gsetgroup]]$set$genesets[[genesetname]], USE.NAMES=F)
      
      get.genesetgroup= function(gsname, hyperdf=hyperdf) as.data.frame(hyperdf) %>% dplyr::filter(original.geneset==!!gsname) %>% pull(geneset) %>% unique
      
        if(length(x)>0 && length(fct)>0 ){
          #fcat("x has more than ozero elements")
        gsgenemat=matrix(nrow=length(x), ncol=length(fct), data=0)
      rownames(gsgenemat)=x
      colnames(gsgenemat)=fct
      
      for(fc in 1:nrow(hyperdf)){
        
        fcgenes=get.gsgenes(hyperdf[fc, "original.geneset"],hyperdf[fc, "geneset"] ) %>% as.vector
        
         if(!is.null(fcgenes)){
           #fcat("fcgenes is not null. making dataframe")
           odf=list(gene=x, geneset=fc, is.in.markerlist=(x %in% fcgenes)) %>% data.frame
           podf= odf %>% dplyr::filter(is.in.markerlist)
           #fcat("rows of podf", nrow(podf), "...")
           if(nrow(podf)>0){
           for(j in 1:nrow(podf)){
           gsgenemat[podf[j, "gene"], podf[j, "geneset"]]<-1
             
             
             }
           }
         }
      }
        
      
        }
           
            if(length(dim( remove.zerorows(gsgenemat))) <=1  || length(dim( remove.zerocols(gsgenemat))) <=1 ){
            warning("problem with final matrix dimensions (either too many 0's or wrong dimensions) Heatmap cannot be created") 
         
         hm=NULL 
         gsgenemat=NULL
         original.gset.names=NULL
       }else{
         
         
         
         
         gsgenemat=remove.zerocols(remove.zerorows(gsgenemat))
        
       }
        
      if(length(dim(gsgenemat)) >1){
        
      fcat("getting gset names")
      original.gset.names=colnames(gsgenemat)     
      colnames(gsgenemat)= lapply(colnames(gsgenemat), function(x) paste0(get.genesetgroup(x, hyperdf), "-", x)) %>% Reduce(c, .)
      
     # fcat("drawing heatmap...")
      colnames(gsgenemat)=summarisename(colnames(gsgenemat), dataset.info$gsea_total_print_length)
      #quickprint(gsgenemat) 
      hm=ComplexHeatmap::pheatmap(gsgenemat, cluster_row=cluster, cluster_col=cluster, border_color = NA, main=paste0("gene set ", setid, ": hit incidence in all gene sets"))
      #hm=NULL  #temporary patch
       
      }else{
      
         warning("either the markerlist or the gene set list is zero. Heatmap cannot be created.") 
        
      
         hm=NULL 
         gsgenemat=NULL
         original.gset.names=NULL
        }
      list(heatmap=hm, matrix=gsgenemat, original.gset.names=original.gset.names)
      
}


################################################################################
# scrutinise the gene sets for a keyword
################################################################################
get.relevant.genesets=function( keyword, gsetpars, exclude="placeholder"){
  
allgsetnames=lapply(gsetpars, function(x) names(x$set$genesets)) %>% Reduce(c, .) 
allgsets=lapply(gsetpars, function(x) x$set$genesets)  %>% Reduce(append, .) %>% givename(., allgsetnames)

relevant.gsetnames=allgsetnames[grepl(keyword, allgsetnames, ignore.case=T)]
relevant.gsetnames=relevant.gsetnames[!grepl(exclude, relevant.gsetnames)]

list(gene.sets=allgsets[relevant.gsetnames], relevant.geneset.names=relevant.gsetnames)
}



get.gene.gsea.relevance=function(genes,keyword,gsetpars,exclude="placeholder", label="",simple.output=T){

  geneset.search=get.relevant.genesets( keyword, gsetpars, exclude=exclude)


all.incidences=lapply(genes, function(gene) (lapply(geneset.search$gene.sets, function(x) gene %in% x) %>% unlist %>% sum)/ length(geneset.search$gene.sets)) %>% unlist() %>% givename(., genes)

if(simple.output){
all.incidences   
}else{
  list(incidence.percentage=all.incidences, relevant.genesets=geneset.search, df=coerce.vector.to.matrix(all.incidences) %>% givecolnames(., nms=paste0("incidence.in.", label, ".sets")) %>% as.data.frame )
}
  
}

################################################################################
#
################################################################################

################################################################################
#
################################################################################


coerce.vector.to.matrix=function(vec)  matrix(vec, ncol = 1, dimnames = list(names(vec), NULL))



################################################################################
# analysis of whether some peaks are in the range of the reference
################################################################################

range.analysis.workflow=function(target.peaks,target_cell_type, test.ids,vsd, pca.data , mtf_mat_gathered, sample.arrangement, keyword=config$gsea_celltype_searchterm, label=config$target_cell_type, exclude=config$gsea_celltype_exclude,enrichments=NULL,below.percentile=.1, pct.always.below=1){


target.ids=pca.data %>% filter(celltype==!!target_cell_type) %>% pull(Experiment)

test.vals=vsd[target.peaks , test.ids ]
ref.vals=vsd[target.peaks , target.ids]

################################################################################
# direct method
################################################################################
#rangetests=lapply( 1:length(target.peaks), function(j) {
#fcat(j, ":", peak.collection$roc[j])
 # rng=range(ref.vals[j, ])
  
  #isinrange= function(x, rng) ifelse(x>rng[1] && x<=rng[2], 0, isabove(x, rng))
 # isabove= function(x, rng) ifelse(x>rng[2], 1, isbelow(x, rng))
  #isbelow= function(x, rng) ifelse(x<=rng[1], -1, NA)
  #isexpressed=function(x) ifelse(x>=min.expression, -1, -2)
  
 # result=lapply(test.vals[target.peaks[j],], function(y) isinrange(y, rng)) %>% unlist()
 # result
  
#}) %>% Reduce(rbind, .) %>% givecolnames(.,nms= colnames(test.vals)) %>% giverownames(., target.peaks)

################################################################################
#percentile method
################################################################################
rangetests=lapply( 1:length(target.peaks), function(j) {

  ec=ecdf(ref.vals[j, ])
  result=ec(test.vals[target.peaks[j],test.ids])
  result=ifelse(result>=below.percentile, -1,0)
  
}) %>% Reduce(rbind,.) %>% givecolnames(.,nms= colnames(test.vals)) %>% giverownames(., target.peaks)

################################################################################
#  prepare the enrichments
################################################################################

howmany.peaks.in.range=apply(rangetests[target.peaks, test.ids], 2, function(x) sum(x==0))
howmany.peaks.below.range=apply(rangetests[target.peaks, test.ids], 2, function(x) sum(x==-1)) 
howmany.peaks.above.range=apply(rangetests[target.peaks, test.ids], 2, function(x) sum(x==1))

peak.range.count.df=list(in.range=howmany.peaks.in.range, below=howmany.peaks.below.range, above=howmany.peaks.above.range) %>% as.data.frame %>% cbind(., pca.data[names(howmany.peaks.in.range), ]) %>% pivot_longer(., cols=c("in.range", "below", "above"), names_to="peak.class", values_to="peak.class.count")

#ggplot(peak.range.count.df, aes(x=tmnt.full, y=peak.class.count, fill=peak.class))+geom_col(position="fill")

peak.freq.below.range=apply(rangetests[target.peaks, test.ids], 1, function(x) sum(x==-1)/length(x)) %>% sort(., decreasing=T)
peaks.always.below.range=peak.freq.below.range[peak.freq.below.range>= pct.always.below] %>% names

peak.freq.in.range=apply(rangetests[target.peaks,test.ids], 1, function(x) sum(x==0)/length(x)) %>% sort(., decreasing=T)
peaks.always.in.range=peak.freq.in.range[peak.freq.in.range==1] %>% names

peak.freq.above=apply(rangetests[target.peaks,test.ids], 1, function(x) sum(x==1)/length(x)) %>% sort(., decreasing=T)
peaks.always.above=peak.freq.in.range[peak.freq.above==1] %>% names

peak.freq.in.or.above=apply(rangetests[target.peaks,test.ids], 1, function(x) (sum(x==0)+sum(x==1))/length(x)) %>% sort(., decreasing=T)
peaks.in.or.above=peak.freq.in.or.above[peak.freq.in.or.above==1] %>% names

#best.sample.peaks.in.range=names(howmany.peaks.in.range)[1]
#which.peaks.in.range=rangetests[all.peaks.celltypes.auc[[target_cell_type]], best.sample.peaks.in.range]==0
#which.peaks.below.range=rangetests[all.peaks.celltypes.auc[[target_cell_type]], best.sample.peaks.in.range]==-1


#peaks.in.range= rownames(rangetests[which.peaks.in.range,])


################################################################################
# make heatmaps for the different classes of peaks.
################################################################################
#categorise peaks with many sorts of active and inactive peaks
fcat("organising peaks according to sensitivity")
above.partitions=lapply(seq(1, 0, -0.05), function(x)  peak.freq.above[peak.freq.above>x] %>% names) %>% Reduce(c, .)
in.range.partitions=lapply(seq(1, 0, -0.05), function(x)  peak.freq.in.range[peak.freq.in.range>x] %>% names) %>% Reduce(c, .)
below.partitions=lapply(seq(0, 1, 0.05), function(x)  peak.freq.below.range[peak.freq.below.range<x] %>% names) %>% Reduce(c, .)
pks=c(above.partitions, in.range.partitions, below.partitions) %>% unique

#define boundary between fractional peaks and eternally low peaks. 
gaps.peaks= length(pks)

pks.complete=c(pks, peaks.always.below.range)

#order arranged by number of peaks in range
smpls=names(howmany.peaks.in.range) %>% rev

# matrix subset to the relevant peaks being shown
range.mat=rangetests[pks.complete, ]

#contingency for lack of a pre-arrangement of samples
if(!is.null(sample.arrangement)){
  sample.arrangement=colnames(range.mat)
}

#seriating the matrix
range.mat.row.ord=seriation::get_order(seriation::seriate(range.mat, margin=1, method="PCA"))
range.mat.col.ord=seriation::get_order(seriation::seriate(range.mat, margin=2, method="PCA"))

fcat("plotting sensitive/insensitive peaks in different ways...")

br=seq(-1,1,1)  ## numeric breaks of the color bins
#colls=c("aquamarine4", "aquamarine", "darkgoldenrod1")  
colls=c("cyan4", "cyan3", "white")  

sca=3
tpdf(path=config$plotpath, paste0("heatmap_inrange_to_", target_cell_type, "rowseriation_pct.threshold_",below.percentile ) %>% addversion, wi=pw*sca*.6, he=pw*sca)
print(ComplexHeatmap::pheatmap(range.mat[range.mat.row.ord, sample.arrangement],
  cluster_row=F,
  cluster_col=F,
  annotation_col=pca.data[sample.arrangement %>% rev,dataset.info$annotvars],
  show_rownames=F,
  show_colnames=F,
  annotation_colors=allcolors,
  use_raster=T,
  breaks=br, 
  col=colls, 
  gaps_row=gaps.peaks 
  ) )
dev.off()

fcat("printing peaks to screen")
peaklist=list(inrange=pks, below.range=peaks.always.below.range) 
#fcat(peaklist %>% dput)

lapply(1:length(peaklist), function(pk){
  
  fwrite(norm.peak.vals[peaklist[[pk]], 2:4] %>% unname %>% giverownames(., NULL), file=paste0("~/ml2cell_code/metadata/", names(peaklist)[pk], ".bed"), sep="\t", row.names=F)
})


################################################################################
#  calculate the enrichments
################################################################################

fcat("calculating enrichments of sensitive/insensitive peaks...")
tfbs.enrichment.list=list(
#highly.matched.peaks=test.tfbs.enrichment(peak.freq.in.or.above[peak.freq.in.or.above>dataset.info$peak_sensitivity_threshold] %>% names, mtf_mat_gathered, target.name="highly.matched.peaks"),
below.range=test.tfbs.enrichment(peak.freq.below.range[peak.freq.below.range==1] %>% names, mtf_mat_gathered, target.name="always.below.range"),
#often.below.range=test.tfbs.enrichment(peak.freq.below.range[peak.freq.below.range>=dataset.info$peak_sensitivity_threshold] %>% names, mtf_mat_gathered, target.name="often.below.range"),
#high.treatment.sensitivity=test.tfbs.enrichment(peak.freq.below.range[peak.freq.below.range>=1-dataset.info$peak_sensitivity_threshold & peak.freq.below.range<=dataset.info$peak_sensitivity_threshold] %>% names, mtf_mat_gathered, target.name="high.treatment.sensitivity") 
  )




fcat("making data frame...")
tfbs.enrichment.df=tfbs.enrichment.list %>% bind_rows()
rownames(tfbs.enrichment.df)=NULL


################################################################################
# prepare contrasts for venn diagram
################################################################################

tfbs.enrichment.list2=tfbs.enrichment.list[lapply(tfbs.enrichment.list, function(x)  !is.null(x)) %>% Reduce(c, .)]
fcat("Making Venn Diagrams of TFs associated with each sensitivity group...")

mtfs.contrast=lapply(tfbs.enrichment.list2, function(x) x %>% pull(gene_symbol))  %>% givename(., names(tfbs.enrichment.list2))


gven=ggVennDiagram(mtfs.contrast, label_alpha=0 )
#ggplot2::scale_fill_gradient(low="blue",high = "pink")


fcat("extracting exclusive tfs for each sensitivity group...")


list.tfbs.exclusive=lapply(1:length(mtfs.contrast), function(x){
  setdiff(mtfs.contrast[[x]],    mtfs.contrast[1:length(mtfs.contrast)!=x] %>% Reduce(c, .)  ) 
}) %>% givename(., names(mtfs.contrast))


#list.tfbs.shared=list(shared.tfbs=intersect(mtfs.contrast$highly.matched.peaks, mtfs.contrast$below.range), 
 # all.tfbs=union(mtfs.contrast$highly.matched.peaks, mtfs.contrast$below.range))



fcat("calculating gene set enrichments for each sensitivity group's TFs over all genesets...")
if(is.null(enrichments)){
  
simpleCache(paste0("enrichments_targets_target_",target_cell_type, "_", digest::digest(target.peaks),"_list_", digest::digest(list.tfbs.exclusive), "_gsetpars_", digest::digest(gsetpars)) %>% addversion, {   
enrichments=lapply(1:length(gsetpars), function(x) get.hyper.enrichments(list.tfbs.exclusive, gsetpars = gsetpars[x], make.barplot=F, make.heatmap = F))  %>% givename(., names(gsetpars))

enrichments},assignToVar="enrichments", reload=T) 
}

fcat("starting to prepare perturbation heatmaps...")

  perts=c("LIGPERTUP|LIGPERTLINCSUP|CHEMPERTUP|KINPERTUP", "LIGPERTDOWN|LIGPERTLINCSDN|CHEMPERTDOWN|KINPERTDOWN", "CGP")
  
  
perturbation.heatmaps=lapply(names(list.tfbs.exclusive), function(cluslab){

fcat("sensitivity group", cluslab, "...")

lapply(perts, function(searchterms){
  
fcat("Looking for groups enrichments in ", searchterms, "...")
endf=lapply(enrichments[grepl(searchterms, names(enrichments))],  function(x) x$df) %>% Reduce(rbind, .) 

#we make the enrichment heatmap matrix whereby we see what gene sets share the same hits for biological insight.
#we keep the matrix to make a more tailored heatmap later. 

fcat("preparing the enrichment matrix...")

eh=make.enrichment.heatmap(x=list.tfbs.exclusive[[cluslab]], gsetpars=gsetpars, endf %>% dplyr::filter(clusterlabel==!!cluslab))

if( !is.null(eh$matrix) ){
mat=remove.zerocols(remove.zerorows(eh$matrix %>% t))


gset.annot=endf %>% dplyr::filter(clusterlabel==!!cluslab) %>% filter(original.geneset %in% eh$original.gset.names) %>% select(jaccard_index, pval)


fcat("finding the relative relevance of each gene in the genesets for cell type ", target_cell_type, " ...")

relev=get.gene.gsea.relevance(colnames(mat), gsetpars, simple.output=F, keyword=keyword, label=target_cell_type, exclude=exclude)
tf.annot=relev$df %>% givecolnames(., nms=paste0("relevance.for.", target_cell_type)) 

#
fcat("ranking the perturbations to put those hitting the most relevant TFs on top")
#find the profiles that most closely involve target-relevant genes
hammings=lapply(1:nrow(mat), function(j){
sum(!(tf.annot[,1]>0 & as.logical(mat[j,])))
}) %>% unlist
horder=order(hammings, decreasing=F)
col_fun <- colorRampPalette(c(0, 1), c( "white", "green"))

fcat("preparing annotations for the perturbation heatmap...")
column_ha <- HeatmapAnnotation(
  annotation = anno_simple(tf.annot[,1]),
  show_annotation_name = TRUE
)


row_ha <- ComplexHeatmap::rowAnnotation( jaccard.index = anno_barplot(gset.annot))
fcat("executing the perturbation heatmap...")
hm.perturbations=ComplexHeatmap::pheatmap(mat[horder, ], col=c("purple", "yellow"), cluster_row=F, main=paste0("potential perturbations to ", ifelse(grepl("UP", searchterms), "activate", "deactivate")," TFs at ",cluslab , " peaks"), annotation_col=tf.annot, right_annotation=row_ha)
tpdf(path=config$plotpath, paste0("heatmap_tfbs_",  "_target_", target_cell_type, cluslab,"_perturbation_signatures", remove.pipes(searchterms)) %>% addversion, wi=pw*sca*3, he=pw*sca*1)
print(hm.perturbations)
dev.off()

}else{
  mat=NULL
  hm.perturbations=NULL
  }
  
  list(matrix=mat, heatmap=hm.perturbations)

}) %>% givename(.,perts %>% remove.pipes())  
  
}) %>% givename(., names(list.tfbs.exclusive))

list(
  range.mat=range.mat,
  tfbs.enrichment.list=tfbs.enrichment.list, 
  tfbs.per.region.all=mtfs.contrast,
  tfbs.per.region.exclusive=list.tfbs.exclusive,
  perturbation.heatmaps=perturbation.heatmaps, 
  sensitive.peaks=pks, 
  insensitive.peaks=peaks.always.below.range, 
  peak.range.count.df,
  venn.diagram=gven, 
  enrichments=enrichments, 
  
  )

}


################################################################################
# function to get enrichments of a list and store it as a cache
################################################################################

get.enrichments.cache=function(list.tfbs, gsetpars, cores=10){
    
  fcat("calculating gene set enrichments for each sensitivity group's TFs over all genesets...")

  cachename=paste0("enrichments_list_", digest::digest(list.tfbs), "_gsetpars_", digest::digest(gsetpars)) %>% addversion
simpleCache(cachename, { 
  
enrichments=parallel::mclapply(1:length(gsetpars), function(x, ll, gp){

simpleCache(paste0("enrichment_", digest::digest(ll), "_gset_", digest::digest(gp[x])) %>% addversion, {   
enr= get.hyper.enrichments(ll, gsetpars = gp[x], make.barplot=F, make.heatmap = F)
enr
}, assignToVar="enr", reload=T)
enr  
},ll=list.tfbs,gp=gsetpars, mc.cores=cores, mc.cleanup=T)

lapply(enrichments, function(x) x$df) %>% Reduce(rbind, .)

},assignToVar="enrichments", recreate=T) 

list(df=enrichments, gsets= names(gsetpars), cachename=cachename)  
}
################################################################################
# module to make perturbation heatmaps based on a tf list
################################################################################

################################################################################
# prepare deeptools script programatically by selecting sample groups from the workspace
# and automatically adding a reference sampe to cross compare
################################################################################


prepare.deeptools.script= function(
      test.class.variable,
  test.class,
    peak.list,
  peaks.id="unnamed-peaks",
  norm.peak.vals,
  meta,
  dataset.info,
  test.filter=NULL,
  invert=F,
  max.test.samples=NULL,
    reftype=config$target_cell_type,
  ref.sample.number=1,
  ref.sample.id=NULL,
  outdir=  NULL, 
  outdir.internal=NULL,
  outfile.matrix=NULL,
outfile.script= NULL, 
  outfile.sortedregions=NULL, 
  num.clusters=4, 
  append.to.file=F, 
  cluster.using.samples=1,
  sort.using.samples=1,
  heatmap.height=28, 
  bin.size=25, 
  specific.samples=NULL,
  yMin=0, 
  yMax=1, 
  zMin=0,
  zMax=2, 
  num.processors=1,
  selected.bigwigs=NULL, #provide the exact path of the bigwigs and ignore everything else
  date.use=NULL, 
  replace=F){
  
  fcat("preparing paths")
  test.class.formatted=gsub(" ", "-", test.class)
  
  newdir=timestamp(paste(c("deeptools_ref",reftype, test.class.variable, gsub("\\|","-", test.class.formatted), peaks.id ), collapse="_"), date.use=date.use)

  
  if(is.null(outdir.internal)){
   coredir=paste0(dataset.info$out_root, "/deeptoolsrerun/")
   
      if(!dir.exists(coredir)){
  dir.create(coredir)
  }
        outdir.internal=paste0(coredir, newdir)
  }
  
  outdir.external=outdir.internal %>% prep.dir.for.out.root
  
  if(!dir.exists(outdir.internal)){
  dir.create(outdir.internal)
  }
  
  if(is.null(outfile.matrix)){
    outfile.matrix=paste0(outdir.external, "/integrated_matrix.gz")
  }else{
  outifle.matrix=outfile.matrix %>% prep.dir.for.out.root
  }
  
    if(is.null(outfile.sortedregions)){
    outfile.sortedregions=paste0(outdir.external, "/sortedregions.bed")
    }else{
    outifle.sortedregions=outfile.sortedregions %>% prep.dir.for.out.root
    }
  
    if(is.null(outfile.script)){
    outfile.script="~/ml2cell_code/bash/deeptools_to_refsample.sh" 
    }
  
    outfile.script.archive=paste0(outdir.internal, "/run_deeptools.sh")
    outfile.scripts.list=list(bash=outfile.script, archive=outfile.script.archive)
    
  
  fcat("preparing peak region files")
  
  peakfilenames=lapply(1:length(peak.list), function(y){
  peak.list.name.formatted=gsub(" ", "-", names(peak.list)[y])  
  peakfilename.internal=paste0(outdir.internal,"/", peak.list.name.formatted, ".bed")
  peakfilename.external=paste0(outdir.external,"/", peak.list.name.formatted, ".bed")
  fcat(y, "extracting bed data")
  matt=norm.peak.vals[peak.list[[y]], 2:4]
  rownames(matt)=NULL
  fcat(y, "writing bed file")
  fwrite(matt, file=peakfilename.internal, row.names=F, col.names=F, sep="\t")
    fcat("making list")
  list(internal=peakfilename.internal, external=peakfilename.external, name= paste(names(peak.list)[y], collapse="_") )
}) %>% givename(., names(peak.list))


  allpeakfilenames.external= lapply(peakfilenames, function(x) x$external)  
  
  fcat("finding reference sample") # pending: merge test samples into one file, but this depends on having mergebigwig
  
  refbigwigpath=dataset.info$dataset_paths_nfcore[as.logical(dataset.info$is_reference)] %>% make.bigwigpath
  pprofilepath=dataset.info$dataset_paths_nfcore[as.logical(dataset.info$is_reference)] %>% make.profilepath 

  all.bigwig.paths.reference=list.files(path=refbigwigpath, pattern=".bigWig") %>% paste0(refbigwigpath,.) %>% prep.dir.for.out.root
  all.dmatrix.paths.reference=list.files(path=pprofilepath, pattern=".scale_regions.computeMatrix.mat.gz") %>% paste0(pprofilepath, .)

  if(is.null(ref.sample.id)){
  refsample=(meta %>% filter(celltype==!!reftype) %>% pull(Experiment))[ref.sample.number]
}else{
   refsample=meta %>% filter(celltype==!!reftype) %>%  grepl(ref.sample.id, Experiment) %>% pull(Experiment)
}
  
bigwig.refsample=grep(refsample, all.bigwig.paths.reference, value=T)
dmatrix.refsample=grep(refsample, all.dmatrix.paths.reference, value=T)
  
  fcat("finding all bigwig files and assempling their paths assembling their paths")
  
  
  if(is.null(selected.bigwigs)){
 # assemble the source paths where to find the bigwigs 
    sourcepaths.bigwig.root=dataset.info$dataset_paths_nfcore %>% make.bigwigpath 
 
  # look for the collection of bigwig files and make their paths
  all.bigwigpaths.parts=lapply(1:length(sourcepaths.bigwig.root), function(x){
  corefiles=list.files(path=sourcepaths.bigwig.root[x], pattern=".bigWig")
  
 list(filenames=corefiles, filepaths=corefiles %>% paste0(sourcepaths.bigwig.root[x], .) )
})
  
  
bigwigpaths.assembled= lapply(all.bigwigpaths.parts, function(x) x$filepaths) %>% Reduce(c, .) %>% prep.dir.for.out.root

#with all the possible paths assemble, now retrieve the relevant ones based on filters

if(is.null(specific.samples)){
selected.bigwigs= grep( (meta %>% filter(grepl(!!test.class, !!sym(test.class.variable))) %>% pull(Experiment) %>% paste(., collapse="|")), bigwigpaths.assembled, value=T )
}else{

  
selected.bigwigs= lapply(specific.samples, function(x){
  
  grep( (meta %>% filter(Experiment==x) %>% pull(Experiment) %>% paste(., collapse="|")), bigwigpaths.assembled, value=T )
  

}) %>% Reduce(c,.)

}


if(!is.null(test.filter)){
  
    selected.bigwigs= grep(test.filter, selected.bigwigs, invert=invert, value=T)
     
  }

  fcat("checking to adjust the total number of test samples displayed")
  if(!is.null(max.test.samples)){
    if(max.test.samples<=length(selected.bigwigs)){
   selected.bigwigs=selected.bigwigs[1:max.test.samples] 
    }
  }
}else{
  
 bigwig.refsample="" 
}
################################################################################
# WRITING SRCIPT  
################################################################################
  fcat("writing script")
  
    outfile.template="~/mnt_out/matrix_B.gz"
  
################################################################################
# prepare computeMAtrix command
################################################################################

computematrix.command=c(
  "computeMatrix reference-point -a 3000 -b 3000 --sortRegions descend -S ", bigwig.refsample, 
  selected.bigwigs,  
  " -bs ", bin.size,
    " --missingDataAsZero",
  " --referencePoint center",
  " -p", num.processors,
  " -R ",
  allpeakfilenames.external,"-o",
   outfile.matrix,
  " --skipZeros --outFileSortedRegions",
  outfile.sortedregions
  )  %>% paste(., collapse=" ")  

################################################################################
# prepare plotHeatmap command
################################################################################


heatmap.command= paste0(
  "plotHeatmap -m ", 
  outfile.matrix, 
  " --sortUsingSamples ",sort.using.samples,
  " --zMin ", zMin, " --zMax ", zMax, " ", 
  " --yMin ", yMin," --yMax ", yMax, " ",
    "  --outFileSortedRegions ", outdir.external, "/sorted_regions_heatmap.bed ",
  "  --outFileNameMatrix ", outdir.external, "/heatmap_matrix.gz ",
  " --heatmapHeight ", heatmap.height,
  " -o ", 
  prep.dir.for.out.root(dataset.info$plotpath), 
  paste_(newdir,"integrated_heatmap.pdf")
  )

heatmap.command.cluster= paste0(
  "plotHeatmap -m ", 
  outfile.matrix, 
  " --sortUsingSamples ",sort.using.samples,
  " --zMin ", zMin, " --zMax ", zMax, " ", 
  " --yMin ", yMin," --yMax ", yMax, " ",
    " --hclust ", num.clusters, " --clusterUsingSamples ", cluster.using.samples,
  "  --outFileSortedRegions ", outdir.external, "/sorted_regions_heatmap_clustered.bed ",
  "  --outFileNameMatrix ", outdir.external, "/heatmap_matrix_clustered.gz ",
  " --heatmapHeight ", heatmap.height,
  " -o ", 
  prep.dir.for.out.root(dataset.info$plotpath), 
  paste_(newdir,"integrated_heatmap_clustered.pdf")
  )



fcat("creating output scripts...")
#there is a working script and an archive script going to bash and a script going to the deeptools output so it can be permantently reproduced later on

#conditions to create the script
# if there is no date specified,  if there is a date specified but replace is true, 
if((is.null(date.use)) || (!is.null(date.use) && replace) ){
lapply(outfile.scripts.list, function(of){
if(!append.to.file){
  fwrite(list("#!/usr/bin/bash"), file=of, append=F)  
}

fwrite(list(computematrix.command), file=of, append=append.to.file)  
fwrite(list(heatmap.command), file=of, append=T)    
fwrite(list(heatmap.command.cluster), file=of, append=T)    
 
})
}

list(computematrix.command=computematrix.command,
  heatmap.command=heatmap.command, 
  heatmap.command.cluster=heatmap.command.cluster,
  path.matrix=paste0(outdir.internal, "/integrated_matrix.gz"),
  path.heatmap.matrix=paste0(outdir.internal, "/heatmap_matrix.gz"),
  path.heatmap.matrix.clustered=paste0(outdir.internal, "/heatmap_matrix_clustered.gz"),
  path.sortedregions=paste0(outdir.internal, "/sorted_regions_heatmap.bed"), 
  path.sortedregions.clustered=paste0(outdir.internal, "/sorted_regions_heatmap_clustered.bed"), 
  selected.bigwigs=selected.bigwigs,
  all.bigwigs=c(bigwig.refsample, selected.bigwigs),
  outdir.internal=outdir.internal, 
  outdir.external=outdir.external,
  outfile.scripts.list=outfile.scripts.list
  )

}

################################################################################
#prepare bigwigmerge command programatically
################################################################################

prepare.bigwigmerge= function(
      test.class.variable,
  test.class,
    peak.list,
  peaks.id="unnamed-peaks",
  norm.peak.vals,
  meta,
  dataset.info,
  test.filter=NULL,
  invert=F,
  max.test.samples=NULL,
    reftype=config$target_cell_type,
  ref.sample.number=1,
  ref.sample.id=NULL,
  outdir=  NULL, 
  outdir.internal=NULL,
  outfile.matrix=NULL,
outfile.script= NULL, 
  outfile.sortedregions=NULL, 
  num.clusters=4, 
  append.to.file=F, 
  cluster.using.samples=1,
  sort.using.samples=1,
  heatmap.height=28, 
  bin.size=25, 
  specific.samples=NULL,
  yMin=0, 
  yMax=1, 
  zMin=0,
  zMax=2
  ){
  
  fcat("preparing paths")
  
  test.class.formatted=gsub(" ", "-", test.class)
  newdir=timestamp(paste(c("bigwigmerge", test.class.variable, gsub("\\|","-", test.class.formatted), peaks.id ), collapse="_"))

  
  if(is.null(outdir.internal)){
   coredir=paste0(dataset.info$out_root, "/bigwigmerge/")
   
      if(!dir.exists(coredir)){
  dir.create(coredir)
  }
        outdir.internal=paste0(coredir, newdir)
  }
  
  outdir.external=outdir.internal %>% prep.dir.for.out.root
  
  if(!dir.exists(outdir.internal)){
  dir.create(outdir.internal)
  }
  
  if(is.null(outfile.matrix)){
    outfile.matrix=paste0(outdir.external, "/integrated_matrix.gz")
  }else{
  outifle.matrix=outfile.matrix %>% prep.dir.for.out.root
  }
  
    if(is.null(outfile.sortedregions)){
    outfile.sortedregions=paste0(outdir.external, "/sortedregions.bed")
    }else{
    outifle.sortedregions=outfile.sortedregions %>% prep.dir.for.out.root
    }
  
    if(is.null(outfile.script)){
    outfile.script="~/ml2cell_code/bash/mergebigwig.sh"
    }
  
  chr.sizes.file=prep.dir.for.project.root("~/ml2cell_code/metadata/GRCh38_EBV.chrom.sizes.tsv")
  

  fcat("finding reference sample") # pending: merge test samples into one file, but this depends on having mergebigwig
  
  refbigwigpath=dataset.info$dataset_paths_nfcore[as.logical(dataset.info$is_reference)] %>% make.bigwigpath
  pprofilepath=dataset.info$dataset_paths_nfcore[as.logical(dataset.info$is_reference)] %>% make.profilepath 

  all.bigwig.paths.reference=list.files(path=refbigwigpath, pattern=".bigWig") %>% paste0(refbigwigpath,.) %>% prep.dir.for.out.root
  all.dmatrix.paths.reference=list.files(path=pprofilepath, pattern=".scale_regions.computeMatrix.mat.gz") %>% paste0(pprofilepath, .)

  if(is.null(ref.sample.id)){
  refsample=(meta %>% filter(celltype==!!reftype) %>% pull(Experiment))[ref.sample.number]
}else{
   refsample=meta %>% filter(celltype==!!reftype) %>%  grepl(ref.sample.id, Experiment) %>% pull(Experiment)
}
  
bigwig.refsample=grep(refsample, all.bigwig.paths.reference, value=T)
dmatrix.refsample=grep(refsample, all.dmatrix.paths.reference, value=T)
  
  fcat("finding all bigwig files and assempling their paths assembling their paths")
  
 # assemble the source paths where to find the bigwigs 
    sourcepaths.bigwig.root=dataset.info$dataset_paths_nfcore %>% make.bigwigpath 
 
  # look for the collection of bigwig files and make their paths
  all.bigwigpaths.parts=lapply(1:length(sourcepaths.bigwig.root), function(x){
  corefiles=list.files(path=sourcepaths.bigwig.root[x], pattern=".bigWig")
  
 list(filenames=corefiles, filepaths=corefiles %>% paste0(sourcepaths.bigwig.root[x], .) )
})
  
  
bigwigpaths.assembled= lapply(all.bigwigpaths.parts, function(x) x$filepaths) %>% Reduce(c, .) %>% prep.dir.for.out.root

#with all the possible paths assemble, now retrieve the relevant ones based on filters

if(is.null(specific.samples)){
selected.bigwigs= grep( (meta %>% filter(grepl(!!test.class, !!sym(test.class.variable))) %>% pull(Experiment) %>% paste(., collapse="|")), bigwigpaths.assembled, value=T )
}else{

  
selected.bigwigs= lapply(specific.samples, function(x){
  
  grep( (meta %>% filter(Experiment==x) %>% pull(Experiment) %>% paste(., collapse="|")), bigwigpaths.assembled, value=T )
  

}) %>% Reduce(c,.)

}


if(!is.null(test.filter)){
  
    selected.bigwigs= grep(test.filter, selected.bigwigs, invert=invert, value=T)
     
  }

  fcat("checking to adjust the total number of test samples displayed")
  if(!is.null(max.test.samples)){
    if(max.test.samples<=length(selected.bigwigs)){
   selected.bigwigs=selected.bigwigs[1:max.test.samples] 
    }
  }
  
  
  heatmap.command= paste0(
  "plotHeatmap -m ", 
  outfile.matrix, 
  " --sortUsingSamples ",sort.using.samples,
  " --zMin ", zMin, " --zMax ", zMax, " ", 
  " --yMin ", yMin," --yMax ", yMax, " ",
    "  --outFileSortedRegions ", outdir.external, "/merged_bigwigs_bedgraph.bg ",
  "  --outFileNameMatrix ", outdir.external, "/heatmap_matrix_clustered.gz ",
  " --heatmapHeight ", heatmap.height,
  " -o ", 
  prep.dir.for.out.root(dataset.info$plotpath), 
  paste_(newdir,"integrated_heatmap.pdf")
  )
  
  ##############################################################################
  # prepare escript
  ##############################################################################
  
  out.bedgraph.nonsorted=paste0(outdir.external, "/merged_bigwigs_bedgraph.bedGraph ")
  out.bedgraph.sorted=paste0(outdir.external, "/merged_bigwigs_bedgraph_sorted.bedGraph ")
  
  mbw.command= paste0("bigWigMerge ",
  paste(selected.bigwigs, collapse=" "),  " ", out.bedgraph.nonsorted )
  
sort.command=paste0("sort -k1,1 -k2,2n ", out.bedgraph.nonsorted, ">", out.bedgraph.sorted)
    
 final.bigwig=paste0(outdir.external, "/merged_bigwigs_", test.class.formatted, ".bw" )
btbb.command=paste0("bedGraphToBigWig ", out.bedgraph.sorted ," ", chr.sizes.file, " ", final.bigwig)




#fwrite(list( "#!/usr/bin/bash"), file=outfile.script, append=append.to.file)
fwrite(list( mbw.command), file=outfile.script, append=append.to.file)
fwrite(list( sort.command), file=outfile.script, append=T, quote=F)
fwrite(list( btbb.command), file=outfile.script, append=T, quote=F)


list(merge.bigwig.command=mbw.command, bedgraphtobigwig.command=btbb.command, sort.command=sort.command, final.bigwig=final.bigwig)



}





umap.louvain.clustering <- function(umap_result, resolution=1) {
  
  # Step 1: Extract the nearest neighbors matrix from the UMAP result
  nn_matrix <- umap_result$nn$euclidean$idx
  
  # Step 2: Create an edge list from the nearest neighbors matrix
  # Create pairs of observations based on their nearest neighbors
  edges <- do.call(rbind, lapply(1:nrow(nn_matrix), function(i) {
    cbind(i, nn_matrix[i, ])
  }))
  
  # Step 3: Create a graph using the edge list
  g <- graph_from_edgelist(edges, directed = FALSE)
  
  # Step 4: Apply the Louvain clustering algorithm to the graph
  louvain_result <- cluster_louvain(g, resolution=resolution)
  
  # Step 5: Extract cluster membership
  cluster_membership <- membership(louvain_result)
  
  # Return the cluster membership vector
  return(cluster_membership)
}


################################################################################
# functions to process information from the deeptools output
################################################################################
get.sample.coords= function(k, extend.from.center=NULL, cols.per.sample){
 startpoint=cols.per.sample*(k-1)+1
 endpoint=cols.per.sample*k
  
 if(!is.null(extend.from.center)){
   ctr=round((endpoint-startpoint)/2)
   
   (ctr-extend.from.center):(ctr+extend.from.center)
   
 }else{
 startpoint:endpoint
 }
}


deeptools.sample.umap=function(smpl, extend.from.center=30, md=.3, nn=20, ncomp=60, rs=1, matth=matth, matts=matts){ 
  submat=matth[,get.sample.coords(1, extend.from.center=extend.from.center)]
mattpca=prcomp(submat, scale=T)
umap.deeptools.object <- umap(mattpca$x[, 1:ncomp], init = "random",  min_dist = md,n_neighbors = nn, repulsion_strength=rs, ret_model=T)
umap.df= umap.deeptools.object$embedding %>% as.data.frame %>% givecolnames(., ind=c(1,2), nms=c("UMAP1", "UMAP2")) %>% cbind(.,matts) %>% dplyr::mutate(value.at.peak.center=submat[, round(extend.from.center/2)])
(ggplot(umap.df, aes(x=UMAP1, y=UMAP2, color=deepTools_group))+geom_point())+ggplot(umap.df, aes(x=UMAP1, y=UMAP2, color=value.at.peak.center))+geom_point()+scale_color_viridis(limits = c(0, 3), oob = scales::squish)
}


deeptools.cross.samples.umap=function(samples, extend.from.center=30, md=.3, nn=20, ncomp=60, rs=1,coll="median", cluster.resolution=1, matth=matth, matts=matts, cols.per.sample){ 
  submat.list=lapply(1:length(samples), function(xx){
  matth[,get.sample.coords(samples[xx], extend.from.center=extend.from.center, cols.per.sample=cols.per.sample)]
  }) 
  

  submat=submat.list %>% Reduce(cbind, .)
  
  submat.medians= apply(submat, 1, median)
  
  submat.centeredmedians=lapply(submat.list, function(mm) mm[, round(extend.from.center/2), drop=T]) %>% Reduce(cbind, .) %>% apply(., 1, median)
      
      peak.center.medians=lapply(1:length(samples), function(xx){
  matth[,get.sample.coords(samples[xx], extend.from.center=extend.from.center, cols.per.sample=cols.per.sample)] 
  }) %>% Reduce(cbind, .)
  
mattpca=prcomp(submat, scale=T)
umap.deeptools.object <- umap(mattpca$x[, 1:ncomp], init = "random",  min_dist = md,n_neighbors = nn, repulsion_strength=rs, ret_model=T, ret_nn=T)

louvain.result=umap.louvain.clustering(umap.deeptools.object, resolution=cluster.resolution)

hclust.result=hclust(dist(mattpca$x[, 1:ncomp]))$labels


umap.df= umap.deeptools.object$embedding %>% as.data.frame %>% givecolnames(., ind=c(1,2), nms=c("UMAP1", "UMAP2")) %>% cbind(.,matts) %>% dplyr::mutate( median= submat.medians, median.centers=submat.centeredmedians, louvain.cluster=factor(as.vector(louvain.result)), hcluster=hclust.result) 

clusternames=unique(umap.df$louvain.cluster)

peak.list=lapply(clusternames, function(ll){

  umap.df %>% filter(louvain.cluster==!!ll) %>% pull(Geneid)  
  
}) %>% givename(., clusternames)

plt=(ggplot(umap.df, aes(x=UMAP1, y=UMAP2, color=deepTools_group))+geom_point())+ggplot(umap.df, aes(x=UMAP1, y=UMAP2, color=!!sym(coll)))+geom_point()

list(umap.df=umap.df, plot=plt, peak.list=peak.list, submat=submat
  )

}


add.peak.annotation.to.mat=function(matt, start.ind=1){
  
  lapply(1:nrow(matt), function(y){
  
  annotation.peaks %>% dplyr::filter(Chr==matt[y, start.ind], Start==matt[y,start.ind+1], End==!!matt[y, start.ind+2])
  
  
}) %>% Reduce(rbind, .) %>% bind_cols(matt,.)

}


peak.clustering.workflow=function(matfile, matfile.heatmap, matfile.sortedregions, nsamples,which.samples.umap, umap.pars, annotation.peaks){

  
list2env(umap.pars)  
library(uwot)
library(igraph)
#library(leiden)
# import the matrix
fcat("importing matrices")
matfile=matfile %>% prep.dir.for.in.root
matfile.heatmap= matfile.heatmap %>% prep.dir.for.in.root
matfile.sortedregions= matfile.sortedregions %>% prep.dir.for.in.root

matt=fread(file=matfile, header=F)
matth=fread(file=matfile.heatmap, header=F) %>% as.data.frame
matts=fread(file=matfile.sortedregions, header=T) %>% as.data.frame


relevant.cols= 7:ncol(matth)
matth.annotation=matth[, 1:6]
matth=matth[rownames(matts), relevant.cols]
cols.per.sample=length(relevant.cols)/nsamples


fcat("annotating reordered peaks")

simpleCache(paste0("reordered_peaks_annotated_matrixid", digest::digest(matts)) %>% addversion, {
matts.annotated=add.peak.annotation.to.mat(matts)
}, assignToVar="matts.annotated", reload=T)
gaps.row=factor(matts[, "deepTools_group"]) %>% as.integer %>% diff %>% as.logical %>% which

 
#ComplexHeatmap::pheatmap(matth[,get.sample.coords(1, extend.from.center=100)], scale="none", cluster_row=F, cluster_col=F, breaks=seq(0,3,.02), gaps_row=gaps.row, annotation_row=matts[, "deepTools_group", drop=F])



#deeptools.sample.umap(1, 30, ncomp=10, md=1)

################################################################################
#make a umap now with each peak containing the information of multiple samples at once
################################################################################



fcat("processing cross sample umap")
peak.clustering.result=deeptools.cross.samples.umap(which.samples.umap, extend.from.center=umap.pars$extend.from.center, ncomp=umap.pars$ncomp, md=umap.pars$md, rs = umap.pars$rs,nn=umap.pars$nn, coll=umap.pars$coll, cluster.resolution=umap.pars$cluster.resolution, matth=matth,  matts=matts.annotated, cols.per.sample=cols.per.sample)
fcat("Done. retrieving result")

fcat("preparing heatmap")

hm=ComplexHeatmap::pheatmap(peak.clustering.result$submat %>% as.matrix, scale="row", cluster_row=F, cluster_col=F, breaks=seq(0,3,.02), gaps_row=gaps.row, annotation_row=peak.clustering.result$umap.df[, c("deepTools_group", "louvain.cluster"), drop=F])
gc()

list(peak.clustering.result=peak.clustering.result, 
  cross.sample.umap=deeptools.cross.samples.umap, 
  matth=matth, 
  matth.annotation=matth.annotation,
  matts=matts.annotated, 
  umap.pars=umap.pars, 
  heatmap=hm
  
)}


get.section=function(matt, section, num.sections=4){
 blocksize=round(nrow(matt)/num.sections) 
cc=1 #record of current matrix index
bc=1 #record of current section
blocks=list()
while(cc< nrow(matt)){ # go on while section is not reached
 if((cc+blocksize)>=nrow(matt)){
   blocks[[bc]]=cc:nrow(matt)
 }else{
   blocks[[bc]]=cc:(cc+blocksize)
 }
 cc=cc+blocksize+1
 bc=bc+1
}  
    sectionmatt= matt[blocks[section] %>% Reduce(c, .), ]
    sectionmatt
 }
  


################################################################################
# parse deeptools script to find all paths inside
################################################################################

### function in progress.
# Function to parse the script and retrieve paths after specific flags
parse_script_paths <- function(script_path) {
  # Read the script from the provided path
  script_lines <- readLines(script_path)
  
  # Initialize the list to store the paths
  final.list <- list()
  
  # Define the flags you are looking for
  flag_patterns <- list(
    "-S" = "-S\\s+([^\\s+]+)",  # Single path after -S
    "--outFileSortedRegions" = "--outFileSortedRegions\\s+([^\\s]+)",  # Single path after --outFileSortedRegions
    "-R" = "-R\\s+([^\\s]+(\\s+[^-\\s][^\\s]+)*)"  # Multiple paths after -R (before next flag)
  )
  
  # Loop through each line in the script
  for (line in script_lines) {
    # For each flag, find matches using regular expressions
    for (flag in names(flag_patterns)) {
      pattern <- flag_patterns[[flag]]
      
      # Match the pattern to extract the paths
      match <- regmatches(line, regexec(pattern, line))
      
      # If a match is found
      if (length(match[[1]]) > 1) {
        # If the flag is -R (multiple paths), split the result by space
        if (flag == "-R") {
          paths <- unlist(strsplit(match[[1]][2], "\\s+"))
        } else {
          paths <- match[[1]][2]  # Single path for other flags
        }
        
        # Store the paths in the list using the flag as the name
        final.list[[gsub("-", "", flag)]] <- paths
      }
    }
  }
  
  # Return the final list of paths
  return(final.list)
}

################################################################################
#two level clustering workflow
################################################################################

two.level.clustering.workflow= function(matt, cor.threshold, clusk, clusk2, label="", allcolors=NULL){  
corr.tfs=cor(matt)
high.cor.tfs=colnames(corr.tfs)[apply(abs(corr.tfs)>=cor.threshold, 2, function(x) sum(x)>10) ] 
corr.tfs2=corr.tfs[high.cor.tfs, high.cor.tfs]

hm.corr=ComplexHeatmap::pheatmap(corr.tfs2, show_rownames=F, show_colnames=F) %>% draw
################################################################################
# define cluster numbers that will define the tf groupings
################################################################################

tfannot.corr=as.data.frame(cutree(as.hclust(column_dend(hm.corr)), k=clusk)) %>% givecolnames(., nms="tfcorr.cluster") %>% dplyr::mutate(tfcorr.cluster=factor(tfcorr.cluster))

clusters.use.corr=tfannot.corr %>% group_by(tfcorr.cluster) %>% summarise(counts=n()) %>% dplyr::filter(counts>2) %>% pull(tfcorr.cluster) %>% as.vector

hm.corr2=ComplexHeatmap::pheatmap(corr.tfs2, annotation_col=tfannot.corr[colnames(corr.tfs2),,drop=F],
  annotation_row=tfannot.corr[colnames(corr.tfs2),,drop=F], show_rownames=F, show_colnames=F) %>% draw()

################################################################################
# define (and fix) cluster colors
################################################################################


simpleCache(paste_("tf_corr_cluster_colors_numclusters", clusk,"clustersid", digest::digest(clusters.use.corr)) %>% addversion, {
 tfcols=hm.corr2@ht_list[[1]]@top_annotation@anno_list$tfcorr.cluster@color_mapping@colors 
 tfcols
}, assignToVar="tfcols",reload=T)
allcolors$TF.cluster=tfcols
names(allcolors$TF.cluster)=paste0("TF.cluster.", names(tfcols))
allcolors$tfcorr.cluster=tfcols


hm.corr3=ComplexHeatmap::pheatmap(corr.tfs2, annotation_col=tfannot.corr[colnames(corr.tfs2),,drop=F],
  annotation_row=tfannot.corr[colnames(corr.tfs2),,drop=F], show_rownames=F, show_colnames=F, annotation_colors=allcolors, use_raster=T) %>% draw()
################################################################################
# hide one of the legends
################################################################################
hm.corr3@ht_list[[1]]@left_annotation@anno_list$tfcorr.cluster@show_legend=F

################################################################################
# plot clusters at the first level 
################################################################################

################################################################################
#further divide samples into two more clusters
################################################################################

tfannot.corr.2=lapply(tfannot.corr %>% pull(tfcorr.cluster) %>% unique, function(xx){
  subsection=tfannot.corr %>% dplyr::filter(tfcorr.cluster==xx) 
  tfs=subsection %>% rownames
  
  # clustering via the correlation
  tr=hclust(dist(corr.tfs2[tfs, tfs] %>% t))
  # clustering via euclidean distance of the tfbs incidences across blocks
  tr2=hclust(dist(matt[, tfs] %>% t))
  
  tfannot.level2.corr=cutree(tr, k=clusk2) %>% coerce.to.matrix %>% givecolnames(., nms="tfcluster.level2.corr") 
  tfannot.level2.euclidean=cutree(tr2, k=clusk2) %>% coerce.to.matrix %>% givecolnames(., nms="tfcluster.level2.euclidean")
  
  tf.annot.level3= bind_cols(subsection, tfannot.level2.corr, tfannot.level2.euclidean) %>% arrange(tfcorr.cluster, tfcluster.level2.corr, tfcluster.level2.euclidean) %>% dplyr::mutate(TF.cluster=paste_(tfcorr.cluster, tfcluster.level2.corr, tfcluster.level2.euclidean))
}) %>% Reduce(rbind,.)


hm.corr4=ComplexHeatmap::pheatmap(corr.tfs2, annotation_col=tfannot.corr.2[colnames(corr.tfs2),,drop=F],
  annotation_row=tfannot.corr.2[colnames(corr.tfs2),,drop=F], show_rownames=F, show_colnames=F) %>% draw()

################################################################################
# hide one of the legends
################################################################################
hm.corr4@ht_list[[1]]@left_annotation@anno_list$tfcorr.cluster@show_legend=F



################################################################################
# prepare heatmap resources
################################################################################

#tfannotation correlation rearranged by correlation cluster
# to use the original correlation clusters, use tfcorr.cluster as a variable
tfacr=tfannot.corr.2 %>% arrange(TF.cluster) %>% select(TF.cluster) %>% dplyr::mutate(TF.cluster=paste0("TF.cluster.", TF.cluster))
tfacr.names=tfannot.corr.2 %>% arrange(TF.cluster) %>% rownames

clusters.use.corr=tfacr %>% group_by(TF.cluster) %>% summarise(counts=n()) %>% dplyr::filter(counts>2) %>% pull(TF.cluster) %>% as.vector





simpleCache(paste_("tf_corr_cluster_colors_numclusters_twolevels", clusk,clusk2, "clustersid", digest::digest(clusters.use.corr)) %>% addversion, {
 tfcols=hm.corr4@ht_list[[1]]@top_annotation@anno_list$TF.cluster@color_mapping@colors 
 tfcols
}, assignToVar="tfcols",reload=T)
allcolors$TF.cluster=tfcols
names(allcolors$TF.cluster)=paste0("TF.cluster.", names(tfcols))
allcolors$tfcorr.cluster=tfcols



hm.corr5=ComplexHeatmap::pheatmap(corr.tfs2, annotation_col=tfacr[colnames(corr.tfs2),,drop=F],
  annotation_row=tfacr[colnames(corr.tfs2),,drop=F], annotation_colors=allcolors, show_rownames=F, show_colnames=F) %>% draw()


sca=2
tpdf(path=config$plotpath, paste0("heatmap_high_correlation_tfs_numclusters_twolevels_",label, clusk, clusk2, "mincor", cor.threshold) %>% addversion, wi=pw*sca, he=pw*sca)
print(hm.corr5)
dev.off()


################################################################################
#define tf groups with more than 2 tfs to find perturbations for
################################################################################

tfbs.assoc.list.corr=lapply(clusters.use.corr, function(x){
  
tff=tfacr %>% filter(TF.cluster==!!x) %>% rownames

#seriate

tff[seriation::get_order(seriation::seriate(matt[, tff], margin=2, method="PCA"))]

}) %>% givename(., clusters.use.corr)

newtff=tfbs.assoc.list.corr %>% Reduce(c, .)

hmtf=ComplexHeatmap::pheatmap(matt[,newtff ]  , cluster_row=F, cluster_col=F, scale="none", annotation_col=tfacr[newtff,,drop=F], show_rownames=F, show_colnames=F, annotation_row=1:nrow(matt) %>% coerce.to.matrix %>% as.data.frame %>% givecolnames(., nms="block") %>% dplyr::mutate(megablock=factor(cut(block,breaks=seq(0,200, 30), labels=F))), annotation_colors=allcolors, use_raster=T) %>% draw()

sca=2
tpdf(path=config$plotpath, paste_("heatmap_chromatinblocks_vs_tfclusters_two_levels",label, clusk,clusk2, "mincor", cor.threshold) %>% addversion, wi=pw*sca, he=pw*sca)
print(hmtf)
dev.off()


list(heatmap.correlations=hm.corr5,
  heatmap.counts=hmtf,
  tfbs.assoc.list.corr=tfbs.assoc.list.corr
  

  )

}


################################################################################
#
################################################################################
prepare.name.annotation=function(matt, genes.to.label, row.where="right", column.where="bottom", linecolor=NULL, labelcolor=NULL){
  
     linesgp= grid::gpar(col=linecolor)
      labelsgp= grid::gpar(col=labelcolor)
  
  fcat("preparing guide vector list")
  guidelist=lapply(c(1,2), function(namedim){
  
    if(namedim==1){
      geneguide=1:dim(matt)[1]
      }else{
        
      geneguide= 1:dim(matt)[2]}
    
    fcat("geneguide is of length", length(geneguide))
    fcat("dimnames(matt)[[",namedim, "]] is of length", length(dimnames(matt)[[namedim]]))
    
    if(is.null(dimnames(matt)[[namedim]])){
      fcat("dimnames not found. filling in")
      dimnames(matt)[[namedim]]=1:dim(matt)[namedim]
    }
    
    names(geneguide)=dimnames(matt)[[namedim]]
    
    geneguide
  })
    #mrkrs are genes of interest
 #dput(guidelist)
  
      getposition.row=Vectorize(function(x) guidelist[[1]][x], USE.NAMES=F)
      getposition.col=Vectorize(function(x) guidelist[[2]][x], USE.NAMES=F)
      
      fcat("preparing annotation objects")
      selected.rows=intersect(names(guidelist[[1]]), genes.to.label)
      selected.cols=intersect(names(guidelist[[2]]), genes.to.label)
      
      fcat("selected rows")
      print(selected.rows)
      
      
       fcat("selected columns")
      print(selected.cols)
      
      
      
    if(length(selected.rows)>0){
      
    row.annomark=  ComplexHeatmap::anno_mark(at = guidelist[[1]][selected.rows] %>% unname %>% as.numeric, 
                                                                       labels = selected.rows,
                                                                       lines_gp= linesgp, 
                                                                       labels_gp= labelsgp)
      
      
    rowannot = ComplexHeatmap::HeatmapAnnotation(which="row", foo = row.annomark)
    }else{
     rowannot=NULL 
     row.annomark=NULL
    }
    
    
     if(length(selected.cols)>0){
       column.annomark=ComplexHeatmap::anno_mark(at = guidelist[[2]][selected.cols] %>% unname %>% as.numeric, 
                                                                       labels = selected.cols,
                                                                       lines_gp= linesgp,
                                                                       labels_gp= labelsgp, 
                                                                                  side="bottom")
      columnannot=ComplexHeatmap::HeatmapAnnotation(which="column", collabels = column.annomark ) 
      }else{
        columnannot=NULL
        column.annomark=NULL
        }

        
list(row=rowannot, row.annomark=row.annomark, column=columnannot, column.annomark=column.annomark)
}
    




getlast= function(x, n) x[(length(x)-(n-1)): length(x)]



zero2na= function(x) ifelse(x==0, NA, x)


################################################################################
# transform a factor into a logical matrix of true or false
################################################################################

create_logical_matrix <- function(factor_vector) {
  # Get the unique levels (categories) in the factor
  categories <- levels(factor_vector)
  
  # Create a matrix where each entry is TRUE if it matches the category, FALSE otherwise
  logical_matrix <- outer(factor_vector, categories, FUN = "==")
  
  # Set column names to be the category names
  colnames(logical_matrix) <- categories
  
  return(logical_matrix)
}



create_colored_list <- function(factor_vector) {
  # Get the unique levels (categories) of the factor
  categories <- levels(factor_vector)
  
  # Create a named list where each element has the vector c("TRUE"="firebrick", "FALSE"="white")
  colored_list <- setNames(
    replicate(length(categories), c("TRUE" = "firebrick", "FALSE" = "white"), simplify = FALSE),
    categories
  )
  
  return(colored_list)
}

################################################################################
#package variables in environment into a list
################################################################################
package.vars <- function(...) {
  # Capture the input arguments as a list
  args <- list(...)
  # Use names of arguments from their deparse expressions
  names(args) <- sapply(substitute(list(...))[-1], deparse)
  return(args)
}



################################################################################
# pca random forest workflow
################################################################################

random.forest.workflow= function(envlist, dataset.info, recreate=F){

  
config=dataset.info  
test.datasets=dataset.info$dataset_name[as.logical(dataset.info$is_protocol)]
allcolors=envlist$allcolors
pca.data=envlist$pca.data
refids=envlist$refids

  
gini.threshold=dataset.info$randomforest_gini_threshold

set.seed(dataset.info$randomseed)

generic.pc.list=paste0("PC", 1:length(envlist$refids))


 
# Train the Random Forest classifier
traininglevels=intersect(names(allcolors[[dataset.info$cell_type_variable]]), unique(pca.data[refids,dataset.info$cell_type_variable ]))   
fcat(traininglevels)

simpleCache(paste_("randomforest_model_seed",config$randomseed, "initialcomponents", config$randomforest_initial_components, "target",config$cell_type_variable, "testds_",paste(test.datasets, collapse="-"),"version", analysis.version), { 

rf_model <- randomForest(x = pca.data[refids, generic.pc.list], y = factor(pca.data[refids, config$cell_type_variable], levels=traininglevels ), ntree = config$randomforest_ntrees)

rf_model
}, assignToVar="rf_model", reload=T)

# Print the trained model
print(rf_model)

rf.err.rate.generic=rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]

confusion.training=rf_model$confusion %>% as.data.frame %>% dplyr::select(-class.error)
confusion.training2=lapply(1:nrow(confusion.training), function(x) confusion.training[x, ]/sum(confusion.training[x,])) %>% Reduce(rbind,. )


sca=1.5
rr=600
hmpgen=ComplexHeatmap::pheatmap(confusion.training2[traininglevels, traininglevels], cluster_col=F, cluster_row=F, main=paste0("Training set. All PCs. Error:", rf.err.rate.generic))
tpdf(path=config$plotpath, paste0("heatmap_randomForest_trainingset_confusionmatrix_analysisversion", analysis.version), wi=pw*sca, he=pw*sca)
print(hmpgen)
dev.off()

# Get variable importance
var_importance <- importance(rf_model)

# Print variable importance
print(var_importance)

# Plot variable importance
sca=1.5
var_imp_plot <- varImpPlot(rf_model)
arr.varimp <- var_imp_plot %>% as.data.frame %>% names2col(., "PC")%>% dplyr::mutate(rank=1:n()) %>% dplyr::filter(MeanDecreaseGini>=1, rank<=20) %>% arrange(MeanDecreaseGini)
ginilevels.generic=arr.varimp %>% pull(PC) 
tpdf(path=config$plotpath, paste_("var.importance.plot.generic", analysis.version), wi=pw*sca, he=pw*sca*1)

print(
  
 ggplot(arr.varimp, aes(x=MeanDecreaseGini, y=factor(PC, levels=ginilevels.generic)))+geom_point(color="black")+cowplot::theme_cowplot()+cowplot::background_grid()+
    geom_vline(xintercept=gini.threshold, color="red")
)  
dev.off()
#var_imp_plot
####prediction of test samples
test_data=pca.data %>% filter(dsname %in% test.datasets) %>% dplyr::select( all_of(generic.pc.list))

testids=pca.data %>% filter(dsname %in% test.datasets)  %>% arrange(tmnt.full) %>% rownames

# Make predictions on the test set
predictions.prob <- predict(rf_model, test_data, type="prob")[testids,]
predictions <- predict(rf_model, test_data)
# Create the confusion matrix
conf_matrix <- table(Actual = rownames(test_data), Predicted = predictions)

color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Define the range for the color bar
# For example, we want the colors to range from -2 to 2
#breaks <- seq(0, 1, length.out=101)
sca=1.5
rr=600

hmpprobs.generic=ComplexHeatmap::pheatmap(predictions.prob, cluster_col=F, cluster_row=F, main=paste0("test prediction probabilities. error:",rf.err.rate.generic))
tpdf(path=config$plotpath, paste0("heatmap_randomForest_test_generic_predictionprobs_analysisversion", analysis.version), wi=pw*sca*1.1, he=pw*sca)
print(hmpprobs.generic)
dev.off()


hmpprobs.generic

################################################################################
# selecting pcs based on their Mean decrease Gini
################################################################################
var.importance.arranged=var_importance %>% as.data.frame %>% arrange(-MeanDecreaseGini)

simpleCache(paste0("informative.pcs_rfmodel", digest::digest(rf_model), "analysisversion_", analysis.version), {
informative.pcs=rownames(var.importance.arranged)[var.importance.arranged>=gini.threshold]
informative.pcs}, assignToVar="informative.pcs",reload=T) 


################################################################################
# use the informative  pcs togenerate a second cleaner model with sleected features
################################################################################


simpleCache(paste_("randomforest_model_informativepcs_seed",config$randomseed, "components", config$randomforest_initial_components, paste(informative.pcs, collapse="-"), "target",config$cell_type_variable, "testds_",paste(test.datasets, collapse="-"),"version", analysis.version), { 

rf_model_informative <- randomForest(x = pca.data[refids, informative.pcs], y = factor(pca.data[refids, config$cell_type_variable], levels=traininglevels), ntree = config$randomforest_ntrees)
rf_model_informative
}, assignToVar="rf_model_informative", reload=T)



rf.err.rate.informative=rf_model_informative$err.rate[nrow(rf_model_informative$err.rate), "OOB"]

print(rf_model_informative)

confusion.trainingi=rf_model_informative$confusion %>% as.data.frame %>% dplyr::select(-class.error)
confusion.training2i=lapply(1:nrow(confusion.trainingi), function(x) confusion.trainingi[x, ]/sum(confusion.trainingi[x,])) %>% Reduce(rbind,. )


sca=1.5
rr=600
hm.training.informative=ComplexHeatmap::pheatmap(confusion.training2i[traininglevels, traininglevels], cluster_col=F, cluster_row=F, main=paste0("Training set (Informative PCs). Error:", rf.err.rate.informative))
tpdf(path=config$plotpath, "heatmap_randomForest_trainingset_confusionmatrix_informativepcs_analysisversion" %>% addversion, wi=pw*sca, he=pw*sca)
print(hm.training.informative)
dev.off()
print(hm.training.informative)

# Get variable importance
var_importancei <- importance(rf_model_informative)

# Print variable importance
print(var_importancei)

# Plot variable importance
var_imp_ploti <- varImpPlot(rf_model_informative)
print(var_imp_ploti)



# Make predictions on the test set
predictions.probi <- predict(rf_model_informative, test_data, type="prob")[testids,]
predictionsi <- predict(rf_model_informative, test_data)
# Create the confusion matrix
conf_matrixi <- table(Actual = rownames(test_data), Predicted = predictionsi)


# Define the range for the color bar
# For example, we want the colors to range from -2 to 2
#breaks <- seq(0, 1, length.out=101)
edm2.samples=pca.data %>% dplyr::filter(dsname %in% test.datasets, tmnt2!="NANANA") %>% rownames
edm1.samples=pca.data %>% dplyr::filter(dsname %in% test.datasets,tmnt2=="NANANA") %>% rownames


ery.order.edm1=predictions.probi[edm1.samples, "Ery"] %>% sort(., decreasing=T) %>% names
ery.order.edm2=predictions.probi[edm2.samples, "Ery"] %>% sort(., decreasing=T) %>% names

sca=1.5
rr=600
hm.inf1=ComplexHeatmap::pheatmap(predictions.probi[ery.order.edm1, traininglevels], cluster_col=F, cluster_row=F, main="test prediction probabilities", annotation_row=pca.data[ery.order.edm1, dataset.info$annotvars ]  %>% dplyr::mutate(rbc.round1=factor(rbc.round1)),annotation_colors=allcolors, show_rownames=F, breaks=seq(0,.3,0.02))
hm.inf2=ComplexHeatmap::pheatmap(predictions.probi[ery.order.edm2,traininglevels ], cluster_col=F, cluster_row=F, main="test prediction probabilities", annotation_row=pca.data[ery.order.edm2, dataset.info$annotvars ] %>% dplyr::mutate(rbc.round1=factor(rbc.round1)), annotation_colors = allcolors, show_rownames = F, breaks=seq(0,.3,.02))

#round(edm1.samples/nrow(predictions.probi))
tpdf(path=config$plotpath, paste0("heatmap_randomForest_test_EDM1_predictionprobs_informativepcs_analysisversion", analysis.version), wi=pw*sca*1.5, he=pw*sca)
print(hm.inf1)
dev.off()

tpdf(path=config$plotpath, paste0("heatmap_randomForest_test_EDM2_predictionprobs_informativepcs_analysisversion", analysis.version), wi=pw*sca*1.5, he=pw*sca)
print(hm.inf2)
dev.off()


  result_list=list(
    var_importance= var_importance,
    var_importancei=    var_importancei,
    rf_model=    rf_model,
    rf_model_informative=    rf_model_informative,
    rf.err.rate.generic=    rf.err.rate.generic,
confusion.training2=confusion.training2,
predictions.prob=predictions.prob,
hmpprobs.generic=hmpprobs.generic,
var.importance.arranged=var.importance.arranged,
var.importance.arranged=var.importance.arranged,
confusion.training2i=confusion.training2i,
rf.err.rate.informative=rf.err.rate.informative,
predictions.probi=predictions.probi[, traininglevels],
    conf_matrix=conf_matrix,
conf_matrixi=conf_matrixi, 
    informative.pcs=informative.pcs
    )


  
  
  
}


################################################################################
# peak roc analysis workflow
################################################################################


compute.aucs.workflow <- function(input_list) {
  
library(parallel)
library(digest)
library(pROC)
library(simpleCache)
library(data.table)
  
  
  # Extract variables from the input list
  positives.celltypes <- input_list$positives.celltypes
  peaks <- input_list$peaks
  meta <- input_list$meta
  refids <- input_list$refids
  analysis.version <- input_list$analysis.version
  max.peaks.pull <- input_list$max.peaks.pull
  norm.peak.vals <- input_list$norm.peak.vals
  mc.cores<- input_list$mc.cores
  
  
  get.sample.peakvals=function(peakid, samples=NULL){
  
  if(!is.null(samples)){ mt=norm.peak.vals[ peakid, samples]}else{ mt=norm.peak.vals[peakid, !colnames(norm.peak.vals) %in% c("Geneid", "Chr", "Start", "End", "Strand", "Length")]}
  
  mt %>% t
  
}
  
  
  # Helper functions
  is.target <- function(x, targets) fifelse(x %in% targets, 1, 0)
  
  make.target.vector <- function(x) is.target(meta[refids, x$classvar], x$targets)
  
  # Precompute target vectors
  target.vector.list <- lapply(positives.celltypes, make.target.vector)
  
  # Precompute peak matrix
  sp <- get.sample.peakvals(peakid = peaks, samples = refids) %>% 
    cbind(., meta[refids, ])

  
  # Compute collective AUCs
  simpleCache(paste0("collective_aucs_celltypes_all_peaksid_", digest::digest(peaks)) %>% addversion, {
    
    parallel.auc.df <- mclapply(names(positives.celltypes), function(x, spv, pks, tvl, av) {
      
      simpleCache(paste0("collective_aucs_celltypes_top", max.peaks.pull, "perPC_target", x, 
                         "_peaks_id", digest::digest(pks), "analysis", av), {
        aucs <- lapply(pks, function(peak) {
          roc(response = tvl[[x]], predictor = spv[, peak], direction = "<", quiet = TRUE)$auc
        }) %>% unlist()
        
        data.frame(list(Geneid = pks, auc = aucs, target = x))
        
      }, assignToVar = "chunk.aucs", reload = TRUE)
      chunk.aucs
    }, spv = sp, pks = peaks, tvl = target.vector.list, mc.cores = mc.cores, av = analysis.version, mc.cleanup = TRUE)
    
    collective.aucs.celltypes <- Reduce(rbind, parallel.auc.df)
    
  }, assignToVar = "collective.aucs.celltypes", reload = TRUE)
  
  # Return the final object
  return(collective.aucs.celltypes)
}






################################################################################
# function to ensure that some packages are loaded, else install them, else fail
################################################################################

ensure.packages <- function(packages) {
  # Ensure BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Function to load or install a package
  load_or_install <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Package '%s' is not installed. Attempting to install...", pkg))
      tryCatch(
        {
          BiocManager::install(pkg, ask = FALSE, update = FALSE)
          if (!requireNamespace(pkg, quietly = TRUE)) {
            install.packages(pkg)
          }
        },
        error = function(e) {
          stop(sprintf("Failed to install package '%s': %s", pkg, e$message), call. = FALSE)
        }
      )
    }
    # Try loading the package
    if (!require(pkg, character.only = TRUE)) {
      stop(sprintf("Failed to load package '%s' after installation attempts.", pkg), call. = FALSE)
    }
  }
  
  # Apply the function to each package
  invisible(lapply(packages, load_or_install))
}




get.perturbation.heatmaps=function(list.tfbs, gsetpars=load.gene.signatures(), enrichments=NULL, keyword=config$gsea_celltype_searchterm, label=config$target_cell_type, exclude=config$gsea_celltype_exclude, tf.annot.external=NULL, colorlist=NULL, cluster=F){
if(is.null(enrichments)){
enrichments=get.enrichments.cache(list.tfbs, gsetpars)
}
perts=c("LIGPERTUP|LIGPERTLINCSUP", "LIGPERTDOWN|LIGPERTLINCSDN", "CGP")
perturbation.heatmaps=lapply(names(list.tfbs), function(cluslab){
fcat("sensitivity group", cluslab, "...")
lapply(perts, function(searchterms){
fcat("Looking for groups enrichments in ", searchterms, "...")
endf= enrichments$df %>% filter(grepl(!!searchterms, geneset))
#we make the enrichment heatmap matrix whereby we see what gene sets share the same hits for biological insight.
#we keep the matrix to make a more tailored heatmap later.
fcat("preparing the enrichment matrix...")
eh=make.enrichment.heatmap(x=list.tfbs[[cluslab]], gsetpars=gsetpars, endf %>% dplyr::filter(clusterlabel==!!cluslab))
if( !is.null(eh$matrix) ){
mat=remove.zerocols(remove.zerorows(eh$matrix %>% t))
gset.annot=endf %>% dplyr::filter(clusterlabel==!!cluslab) %>% filter(original.geneset %in% eh$original.gset.names) %>% select(jaccard_index, pval)
fcat("finding the relative relevance of each gene in the genesets for cell type ", target_cell_type, " ...")
relev=get.gene.gsea.relevance(colnames(mat), gsetpars, simple.output=F, keyword=keyword, label=target_cell_type, exclude=exclude)
tf.annot=relev$df %>% givecolnames(., nms=paste0("relevance.for.", target_cell_type))
if(!is.null(tf.annot.external)){
tf.annot=tf.annot %>% bind_cols(., tf.annot.external[rownames(tf.annot),,drop=F ])
}
#
fcat("ranking the perturbations to put those hitting the most relevant TFs on top")
#find the profiles that most closely involve target-relevant genes
hammings=lapply(1:nrow(mat), function(j){
sum(!(tf.annot[,1]>0 & as.logical(mat[j,])))
}) %>% unlist
horder=order(hammings, decreasing=F)
#col_fun <- colorRamp2(c(0, 1), c( "white", "green"))
fcat("preparing annotations for the perturbation heatmap...")
column_ha <- HeatmapAnnotation(
annotation = anno_simple(tf.annot[,1]),
show_annotation_name = TRUE
)
row_ha <- ComplexHeatmap::rowAnnotation( jaccard.index = anno_barplot(gset.annot))
fcat("executing the perturbation heatmap...")
hm.perturbations=ComplexHeatmap::pheatmap(mat[horder, ], col=c("darkmagenta", "yellow"), cluster_row=F, cluster_col=cluster, main=paste0("potential perturbations to ", ifelse(grepl("UP", searchterms), "activate", "deactivate")," TFs at ",cluslab , " peaks"), annotation_col=tf.annot, right_annotation=row_ha, annotation_colors=colorlist)
tpdf(path=config$plotpath, paste0("heatmap_tfbs_",  "_target_", target_cell_type, label, cluslab,"_perturbation_signatures", remove.pipes(searchterms)) %>% addversion, wi=pw*sca*3, he=pw*sca*1)
print(hm.perturbations)
dev.off()
}else{
mat=NULL
hm.perturbations=NULL
}
list(matrix=mat, heatmap=hm.perturbations)
}) %>% givename(.,perts %>% remove.pipes())
}) %>% givename(., names(list.tfbs))
perturbation.heatmaps}





###############################################################################
# get the ath submatrix of a larger matrix equally divided into n submatrices
###############################################################################


get.submatrix <- function(mat, a, n) {
  # Check if inputs are valid
  if (!is.matrix(mat)) stop("Input must be a matrix.")
  if (a <= 0 || a > n) stop("Argument 'a' must be between 1 and 'n'.")
  if (n <= 0) stop("Argument 'n' must be a positive integer.")
  
  # Determine column ranges for splitting
  total_cols <- ncol(mat)
  col_ranges <- split(1:total_cols, cut(1:total_cols, n, labels = FALSE))
  
  # Extract the requested submatrix
  submatrix <- mat[, col_ranges[[a]], drop = FALSE]
  return(submatrix)
}

# get an interval of -span:span of the middle of a matrix
get.half.interval=function(matt, span=10){ md=round(ncol(matt)/2); matt[,(md-span):(md+span)]}

get.center.sums=function(matt, nsubmats=3, span=10){
  lapply(1:nsubmats, function(x)
get.half.interval(get.submatrix(matt %>% as.matrix, x, nsubmats), span) %>% rowSums
)
}

simple.heatmap=function(matt, ...){ ComplexHeatmap::pheatmap(matt, cluster_row=F, cluster_col=F, scale="none", ...)}
row.heatmap=function(matt, ...){ ComplexHeatmap::pheatmap(matt, cluster_row=F, cluster_col=F, scale="row", ...)}
cluster.heatmap=function(matt, ...){ ComplexHeatmap::pheatmap(matt, cluster_row=T, cluster_col=T, scale="none", ...)}





ratio.to.first <- function(input_matrix) {
  # Ensure the input is a matrix
  if (!is.matrix(input_matrix)) {
    stop("Input must be a matrix.")
  }
  
  # Get the number of columns in the matrix
  num_cols <- ncol(input_matrix)
  
  # Ensure the matrix has at least two columns
  if (num_cols < 2) {
    stop("Matrix must have at least two columns.")
  }
  
  # Extract the first column
  first_col <- input_matrix[, 1]
  
  # Compute the average of the second-to-last column
  second_to_last_col <- input_matrix[, num_cols - 1]
  avg_second_to_last <- mean(second_to_last_col)
  
  # Divide the average of the second-to-last column by the first column
  result <- avg_second_to_last / first_col
  
  return(result)
}



dump.table=function(x, ...){
  variable_name <- deparse(substitute(x))
  fwrite(x, file=paste0(config$out_root, "/metadata/", variable_name, ".csv"), ...)

}

dump.table2=function(x, ...){
  variable_name <- deparse(substitute(x))
  fwrite(x, file=paste0("~/mnt_out/metadata/", variable_name, ".csv"), ...)

}

###############################################################################
# seriate columns of a matrix
##############################################################################

seriatecols=function(x) x[,seriation::get_order(seriation::seriate(x, margin=2, method="PCA"))]

seriaterows=function(x) x[seriation::get_order(seriation::seriate(x, margin=1, method="PCA")),]

cumulative.lengths <- function(lst) {
  # Calculate the lengths of each vector
  lengths <- sapply(lst, length)
  
  # Compute the cumulative sum of lengths
  cumulative <- cumsum(lengths) %>% unname
  
  # Return the result as a named list

  
  
  return(cumulative)
}

################################################################################
# function to generate treatment columns
################################################################################

generate_treatment_columns <- function(mat) {
  # Convert matrix to a data frame for easier handling
  df <- as.data.frame(mat)
  
  # Get the suffixes from column names
  suffixes <- unique(sub(".*\\.", "", colnames(df)))
  
  # Process each suffix
  for (suf in suffixes) {
    # Find columns with the current suffix
    cols_with_suffix <- grep(paste0("\\.", suf, "$"), colnames(df), value = TRUE)
    
    # Sort columns alphabetically and paste row values together
    df[[paste0("tmnt.", suf)]] <- apply(df[cols_with_suffix], 1, function(row) {
      # Replace NAs with "NA", then paste
      paste0(ifelse(is.na(row), "NA", row), collapse = "")
    })
  }
  
  # Create the tmnt.full column by combining all "tmnt.*" columns
  treatment_columns <- grep("^tmnt\\.\\d+$", colnames(df), value = TRUE)
  df$tmnt.full <- apply(df[treatment_columns], 1, function(row) {
    paste(row, collapse = "_")
  })
  
  # Return the updated matrix
  return(df)
}



#############################################
# second version generates round column
#############################################
generate_treatment_columns2 <- function(data) {
  # Identify the unique suffixes in the column names
  suffixes <- unique(gsub(".*\\.", "", grep("\\.", colnames(data), value = TRUE)))
  
  # Loop through the suffixes to create the treatment columns
  for (suffix in suffixes) {
    # Identify the columns with the current suffix
    cols_with_suffix <- grep(paste0("\\.", suffix, "$"), colnames(data), value = TRUE)
    # Sort columns alphabetically and concatenate their values
    data[[paste0("tmnt.", suffix)]] <- apply(data[cols_with_suffix], 1, function(row) {
      paste0(ifelse(is.na(row), "NA", row), collapse = "")
    })
  }
  
  # Add the "tmnt.full" column
  tmnt_cols <- grep("^tmnt\\.\\d+$", colnames(data), value = TRUE)
  data$tmnt.full <- apply(data[tmnt_cols], 1, function(row) {
    paste(row, collapse = "_")
  })
  
  # Add the "round" column
  data$round <- apply(data[tmnt_cols], 1, function(row) {
    # Determine the last non-"NANANA" column
    last_non_na_round <- max(which(row != "NANANA"), na.rm = TRUE)
    if (is.finite(last_non_na_round)) {
      paste0("Round", last_non_na_round)
    } else {
      NA # Case where all columns are "NANANA"
    }
  })
  
  return(data)
}



get.last.element=function(x) x[length(x)]

get.element= function(a,rr, cc ) x[rr, cc]
get.row= function(a, rr) x[,cc]
get.column= function(a, cc) x[,cc ]

get.nonzero= function(x) x[x!=0]



################################################################################
# sum rows with the same id
################################################################################
sum.rows.by.id <- function(mat, ids) {
  fcat("enter sum.rows.by.id")
  # Convert to data frame if rownames are provided as IDs
  if (!missing(ids)) {
    mat <- cbind(ids=ids, as.data.frame(mat))
  }
  
  # Group by the 'ids' column and sum all numeric columns
  aggregated.matrix <- aggregate(. ~ ids, data =mat, FUN = sum)
  
  # Convert back to a matrix if necessary
  rownames(aggregated.matrix) <- aggregated.matrix$ids
  aggregated.matrix <- as.matrix(aggregated.matrix %>% select(-ids))# Remove the 'ids' column
  fcat("exit sum.rows.by.id")
  return(aggregated.matrix)
}


################################################################################
# search files matching strings in a vector inside a file
################################################################################


search.files.in.folder <- function(strings, folder_path) {
  # Ensure the folder path exists
  if (!dir.exists(folder_path)) {
    stop("The specified folder does not exist.")
  }
  
  # Get the list of files in the folder
  all_files <- list.files(folder_path, full.names = FALSE)
  
  # Create a function to find matching files for each string
  find_matches <- function(string) {
    matches <- grep(string, all_files, value = TRUE)
    if (length(matches) == 0) {
      return(NA)  # If no match, return NA
    } else {
      return(paste(matches, collapse = ","))  # If multiple matches, concatenate with commas
    }
  }
  
  # Apply the function to each string
  result <- sapply(strings, find_matches)
  
  # Return the vector of results
  return(result)
}



################################################################################
#  get cis peaks
################################################################################


get.cis.peaks=function(args){
 #find the closest peaks
  
  present.peaks=intersect(rownames(args$mtf_mat_gathered), args$peaks)
  
      if(length(present.peaks)< length(args$peaks)){
   fcat("warning: not all peaks found in motif matrix") 
  }
  closest.genes=args$closest_tss_global %>% filter(Geneid %in% present.peaks) %>% arrange(factor(Geneid, levels=present.peaks)) %>% select(Geneid, gene_symbol) 

  
  relevant.motif.freqs= apply(args$mtf_mat_gathered[present.peaks, , drop=F]>0 , 2, sum) %>% get.nonzero %>% sort(., decreasing=T)
  
  relevant.motif.names= apply(args$mtf_mat_gathered[present.peaks, , drop=F]>0 , 2, sum) %>% get.nonzero %>% sort(., decreasing=T) %>% names
  
  
  ## sum motifs from regions that are nearest to the same tf  
sum.mat=sum.rows.by.id(mtf_mat_gathered[closest.genes$Geneid, relevant.motif.names], closest.genes$gene_symbol) 

sum.mat=sum.mat[(rownames(sum.mat) %in% relevant.motif.names),]  
  ### which tfs regulate their own regions
  
  which.are.cis=intersect(closest.genes %>% pull(gene_symbol), relevant.motif.names)
  list(
      gene.motif.mat=sum.mat,
    closest.genes=closest.genes,
    relevant.motif.names=relevant.motif.names,
    which.are.cis=which.are.cis,
  percentage.cis.regions= sum(closest.genes$gene_symbol %in% which.are.cis)/nrow(closest.genes),
  percentage.cis.tfs= sum(relevant.motif.names %in% which.are.cis)/length(relevant.motif.names)
  
  )
  

  
}


get.tfs=function(args){
 #find the closest peaks
  
  closest.genes=args$closest_tss_global %>% filter(Geneid %in% args$peaks) %>% arrange(factor(Geneid, levels=args$peaks)) %>% select(Geneid, gene_symbol) %>% dplyr::mutate(is.tf=istf(gene_symbol)) 

list(
  peak.info=closest.genes,
 tf.info= closest.genes %>% filter(is.tf),
 counts= closest.genes %>% filter(is.tf) %>% group_by(gene_symbol) %>% summarise(counts=n()) %>% arrange(-counts),
 tfs=closest.genes %>% filter(is.tf) %>% group_by(gene_symbol) %>% summarise(counts=n()) %>% arrange(-counts) %>% pull(gene_symbol) %>% unique
  )
  
}
  
 
stdplot=function(x, title, w=1, h=1, sca=1.5){
 tpdf(path=config$plotpath, title, wi=pw*sca*w, he=pw*sca*h)
  print(x)
  dev.off()
}

# Same as stdplot, but writes PNG. Uses global `w` (pixel base) instead of `pw`
# (inch base) for sizing. Default res is high quality (600 dpi).
stdpng=function(x, title, w=1, h=1, sca=1.5, res=600){
  # Function args `w`/`h` are multipliers; global `w` is the pixel base scale.
  base <- get("w", envir = .GlobalEnv)
  tpng(
    title,
    path = config$plotpath,
    width = base * sca * w,
    height = base * sca * h,
    res = res
  )
  print(x)
  dev.off()
}
  

  
  
capvalues=function(x, lim=5) ifelse(x>=lim, lim, x)



get.targets= function(tf, pks){
  
submat=mtf_mat_gathered[intersect(rownames(mtf_mat_gathered), pks), tf, drop=F]  

  targets=rownames(submat[submat[, 1]>0,, drop=F])
}


get.targets.mtfs= function(tf, pks){
  
submat=mtf_mat_gathered[intersect(rownames(mtf_mat_gathered), pks), tf, drop=F]  

  mtfs=submat[submat[, 1]>0,1]
  list(
  numbers=mtfs, 
  mean=mean(mtfs),
  std=sd(mtfs),
  median=median(mtfs),
  p25=quantile(mtfs, .25),
   p75=quantile(mtfs, .75), 
  )
  
}

shared.regulation= function(tf1, tf2, pks, mtf_mat_gathered){
  
submat=mtf_mat_gathered[intersect(rownames(mtf_mat_gathered), pks), c(tf1, tf2 ), drop=F] %>% as.data.frame %>% dplyr::mutate(both=(!!sym(tf1)>0 & !!sym(tf2)>0), either=(!!sym(tf1)>0 | !!sym(tf2)>0))

pcts=submat %>% summarise( both.pct= sum(both)/n(), either.pct=sum(either)/n())

both.peaks=submat %>% dplyr::filter(both) %>% rownames
either.peaks=submat %>% dplyr::filter(either) %>% rownames
both= length(both.peaks) 
either=length(either.peaks)

      either.pct=pcts$either.pct
    both.pct=pcts$both.pct
  package.vars(
    tf1,
    tf2,
    both.peaks,
    either.peaks,
    both,
    both.pct,
    either,
either.pct
  )
  
}


filter.matrix <- function(cor_matrix, threshold) {
  # Check if the input is a matrix
  if (!is.matrix(cor_matrix)) {
    stop("The input must be a correlation matrix.")
  }
  
  # Ensure the diagonal is excluded
  diag(cor_matrix) <- 0
  
  # Identify the rows and columns to keep based on the threshold
  keep_indices <- apply(cor_matrix, 1, function(row) any(abs(row) >= threshold))
  
  # Filter the matrix by the identified indices
  filtered_matrix <- cor_matrix[keep_indices, keep_indices]
  diag(filtered_matrix) <- 1
  return(filtered_matrix)
}



list.vals.in.cols= function(x){
 
  vals=lapply(colnames(x), function(xx){
  x[, xx] %>% unname  
  }) %>% Reduce(c, .) %>% as.vector() %>% unique 
  
  vals
   
}



make.cache.code <- function(var, tt, varstring)  {
  # Construct the string for simpleCache code
  code <- sprintf(
    "simpleCache('%s', { %s <- %s }, assignToVariable = TRUE, reload = TRUE)", 
    tt, varstring, var
  )
  
  # Return the generated code as a string
  return(code)
}




make.hyper.barplot=function(hdf, genesets=c("LIGPERTDOWN", "LIGPERTLINCSDN"), numbars=10, color="violet",  x="jaccard_index", colorvar="clusterlabel"){
prep.df=hdf %>% dplyr::filter(tlog.pval>0, geneset %in% !!genesets) %>% group_split(geneset, clusterlabel) %>% lapply(., function(x) x %>% arrange(-jaccard_index) %>% group_by(geneset) %>%  dplyr::mutate(rank=1:n(), cpval=categorize.p.value(pval)) %>% head(numbars)) %>% Reduce(rbind, .) %>% as.data.frame

 bp= ggplot(prep.df, aes(x=rank,  y=!!sym(x)) )+geom_col( alpha=0.7, aes(fill=!!sym(colorvar), label=original.label))+geom_text(inherit.aes=F, aes(label=original.label, x=rank),hjust=0, y=0 )+
#geom_text(inherit.aes=F, color="white",aes(label=hits, x=rank,hjust=0, y=0 ,vjust=-.5))+
   rotatex(90)+
  facet_grid(rows=vars(factor(clusterlabel)),cols=vars(factor(geneset, levels=genesets)), scales="free_x")+
  coord_flip()+theme_minimal()+scale_x_reverse()
 bp
  
}

make.hyper.barplot.padj=function(hdf, genesets=c("LIGPERTDOWN", "LIGPERTLINCSDN"), numbars=10, color="violet"){
prep.df=hdf %>% dplyr::filter(tlog.pval>0, geneset %in% !!genesets) %>% group_split(geneset, clusterlabel) %>% lapply(., function(x) x %>% arrange(-jaccard_index) %>% group_by(geneset) %>%  dplyr::mutate(rank=1:n(), cpval=categorize.p.value(pval)) %>% head(numbars)) %>% Reduce(rbind, .) %>% as.data.frame

 bp= ggplot(prep.df, aes(x=rank,  y=jaccard_index) )+geom_col(aes(fill=tlog.padj), alpha=0.7)+geom_text(inherit.aes=F, aes(label=short.label, x=rank,hjust=0, y=0 ))+
   #geom_text(inherit.aes=F, color="white",aes(label=hits, x=rank,hjust=0, y=0 ,vjust=-.5))+
  rotatex(90)+
   #scale_fill_viridis(option="C")+
   scale_fill_gradient(low="#DDDDDD", high=color)+
  #scale_fill_manual(values=allcolors$motif.clusters)+
  #scale_color_manual(values=allcolors$motif.clusters)+
  facet_grid(rows=vars(factor(clusterlabel)),cols=vars(factor(geneset, levels=genesets)), scales="free_x")+
  coord_flip()+theme_minimal()+scale_x_reverse()
 bp
  
}

make.hyper.barplot.pval=function(hdf, genesets=c("LIGPERTDOWN", "LIGPERTLINCSDN"), numbars=10, color="violet"){
prep.df=hdf %>% dplyr::filter(tlog.pval>0, geneset %in% !!genesets) %>% group_split(geneset, clusterlabel) %>% lapply(., function(x) x %>% arrange(-jaccard_index) %>% group_by(geneset) %>%  dplyr::mutate(rank=1:n(), cpval=categorize.p.value(pval)) %>% head(numbars)) %>% Reduce(rbind, .) %>% as.data.frame

 bp= ggplot(prep.df, aes(x=rank,  y=jaccard_index) )+geom_col(aes(fill=tlog.pval), alpha=0.7)+geom_text(inherit.aes=F, aes(label=short.label, x=rank,hjust=0, y=0 ))+
   #geom_text(inherit.aes=F, color="white",aes(label=hits, x=rank,hjust=0, y=0 ,vjust=-.5))+
  rotatex(90)+
   #scale_fill_viridis(option="C")+
   scale_fill_gradient(low="#DDDDDD", high=color)+
  #scale_fill_manual(values=allcolors$motif.clusters)+
  #scale_color_manual(values=allcolors$motif.clusters)+
  facet_grid(rows=vars(factor(clusterlabel)),cols=vars(factor(geneset, levels=genesets)), scales="free_x")+
  coord_flip()+theme_minimal()+scale_x_reverse()
 bp
  
}

make.hyper.barplot.cluster=function(hdf, genesets=c("LIGPERTDOWN", "LIGPERTLINCSDN"), numbars=10){
prep.df=hdf %>% dplyr::filter(tlog.pval>0, geneset %in% !!genesets) %>% group_split(geneset, clusterlabel) %>% lapply(., function(x) x %>% arrange(-jaccard_index) %>% group_by(geneset) %>%  dplyr::mutate(rank=1:n(), cpval=categorize.p.value(pval)) %>% head(numbars)) %>% Reduce(rbind, .) %>% as.data.frame

 bp= ggplot(prep.df, aes(x=rank,  y=jaccard_index) )+geom_col(aes(fill=clusterlabel), alpha=0.7)+geom_text(inherit.aes=F, aes(label=short.label, x=rank,hjust=0, y=0 ))+
   #geom_text(inherit.aes=F, color="white",aes(label=hits, x=rank,hjust=0, y=0 ,vjust=-.5))+
  rotatex(90)+
   #scale_fill_viridis(option="C")+
  # scale_fill_gradient(low="#DDDDDD", high="violet")+
  scale_fill_manual(values=allcolors$motif.clusters)+
  #scale_color_manual(values=allcolors$motif.clusters)+
  facet_grid(rows=vars(factor(clusterlabel)),cols=vars(factor(geneset, levels=genesets)), scales="free_x")+
  coord_flip()+theme_minimal()+scale_x_reverse()
 bp
  
}


################################################################################
# prepare cache command
################################################################################

make.cache.name<- function(title, analysis.version=config$analysis.version, varlist) {
  # Capture additional named arguments
  argnames<- names(varlist)
  # Format each variable as "nameValue"
  args_formatted <- lapply(1:length(argnames), function(nm) paste_(argnames[nm], varlist[[nm]]))
  
  # Concatenate all components with "_"
  filename <- paste(c(title, args_formatted, 'analysisversion', analysis.version), collapse = "_")
  
  return(filename)
}



################################################################################
# recreation/ reloading controls
################################################################################  
fullrecreate.list=list(
  cors=T,
  rgeneposint=T,
  fulltable=T
)

fullreload.list=list(
  cors=F,
  rgeneposint=F,
  fulltable=F)

fullrecreate=function(x){
  
  if(!(x %in% names(fullrecreate.list))){
    return(F)
  }else{return(fullrecreate.list[[x]])}
}

fullreload=function(x){
  
  if(!(x %in% names(fullreload.list))){
    return(T)
  }else{return(fullreload.list[[x]])}
}  


################################################################################
# get var name
################################################################################

getvarname=function(var) deparse(substitute(var))

#########################################
# clinical data formatting functions
#########################################



pctile= function(varname, val){
 vals=gpcdf %>% pull(!!sym(varname))
  dd= ecdf(vals) 
  v=dd(val)
  
  v
}

classify= Vectorize(function(v, lowthresh, topthresh, label.high="high", label.med="med", label.low="low"){
  if (v>topthresh){return(label.high)}else{
  if (v<lowthresh){return(label.low)}else{
  return(label.med) }}
}, USE.NAMES=F)

classify.percentiles=function(vec, lowthresh, topthresh,  label.high="high", label.med="med", label.low="low"){

    classify(ecdf(vec)(vec),lowthresh, topthresh, label.high=label.high, label.med=label.med, label.low=label.low)
 
       
}

# classify <- function(x, low_thresh, high_thresh) {
#   if (!is.numeric(x) || !is.numeric(low_thresh) || !is.numeric(high_thresh)) {
#     stop("All inputs must be numeric.")
#   }
#   
#   if (low_thresh > high_thresh) {
#     stop("low_thresh should be less than or equal to high_thresh.")
#   }
#   
#   # Use vectorized conditions
#   classification <- ifelse(x < low_thresh, "low",
#                     ifelse(x <= high_thresh, "med", "high"))
#   
#   return(classification)
# }
# 


format.censoring=Vectorize(function(x) if(x=="Censored"){return(0)}else{
  if(x %in% c("Event", "Progression", "Relapse", "Death")){return(1)}
  }, USE.NAMES=F)




add.vital.status=function(x){
 if(!("Vital_Status" %in% colnames(x))){
x= x %>% mutate(Vital_Status=format.status(death.from.disease))
 }
   x
}



harmonise.colname= function(df, name.options, target.name){
   df %>% dplyr::mutate(!!sym(target.name) := df[, name.options[data.table::first(which(name.options %in% colnames(df))) ]])
    
    }


format.r2.meta<- function(met){ met2=met %>% t; data.frame(data=met2[3:nrow(met2),]) %>% givecolnames(., nms=make.names(met2[1,]))}




harmonise.values <- function(df, search, renameto, lst=NULL, default.to=NA, replace=F) {

#check if the output column already exists
  duplicate_positions <- which(names(df) == renameto)
  
  if(length(duplicate_positions) > 0 && replace==F){
      warning(sprintf("Column %s already found. returning as is.",renameto) )
      return(df)
  }
  
  # Check if there are duplicates
  if (length(duplicate_positions) > 1) {
    warning(
      sprintf("Duplicate column names found for '%s' at positions: %s", 
              renameto, paste(duplicate_positions, collapse = " ")))
        df2 <- df[, duplicate_positions, drop = FALSE]
          df[,duplicate_positions[1:length(dupicate_positions)]]=NULL
    
  }
names(df)=make.unique(names(df))  
  

# Step 1: Find matching column names
matched_cols <- names(df)[grepl(paste(search, collapse="|"), names(df))]
fcat("found matching columns to ", search,":", matched_cols)
# Step 2: Rename the first matching column (if any) to `renameto`
if (length(matched_cols) > 0) {
df=df %>% dplyr::mutate(!!sym(renameto):= !!sym(matched_cols[1]))
} else {
fcat("No matching column found in `df` for", search," terms. Defaulting column to", default.to)
    
df <- df %>% dplyr::mutate(!!sym(renameto):=NA)
}
# Step 3: Replace values in the renamed column using the `lst` lookup
if (!is.null(lst)) {
df[[renameto]] <- unlist(lapply(df[[renameto]], function(x) {
if (x %in% names(lst)) {
return(lst[[x]])  # Replace if found in `lst`
} else {
return(x)  # Keep original if not in `lst`
}
}))



    # if(!is.null(df2) && (renameto %in% colnames(df2))){
    #     fcat("Warning: duplicate output columns already found in dataframe") 
    #   df=cbind(df, df2)
    #   names(df)=make.unique(names(df))
    # }
}
return(df)
}





#############################
#  find duplicate row name in a matrix, fix them, warn about them and print them to screen
#############################

fix.duplicate.rownames = function(mat) {
  # Ensure input is a matrix
  if (!is.matrix(mat)) stop("Input must be a matrix.")
  
  row_names <- rownames(mat)
  
  # Find duplicate row names
  dup_names <- row_names[duplicated(row_names) | duplicated(row_names, fromLast = TRUE)]
  
  if (length(dup_names) > 0) {
    warning("Duplicates found in row names. Modifying them with suffixes...")
    
    name_counts <- table(row_names)  # Count occurrences
    name_index <- setNames(rep(0, length(name_counts)), names(name_counts))  # Track count
    
    for (i in seq_along(row_names)) {
      row_name <- row_names[i]
      if (name_counts[row_name] > 1) {  # If it's a duplicate
        name_index[row_name] <- name_index[row_name] + 1
        new_name <- paste0(row_name, "_dup", name_index[row_name])
        row_names[i] <- new_name  # Assign new name
        
        # Print the modified duplicate row
        cat(sprintf("Row %d: %s -> %s\n", i, row_name, new_name))
      }
    }
    
    rownames(mat) <- row_names  # Update row names
  }
  
  return(mat)
}




make.unique.names = function(vec) {
  # Count occurrences of each name
  name_counts <- table(vec)
  name_index <- setNames(rep(0, length(name_counts)), names(name_counts))  # Track count
  unique_vec <- vec  # Copy original vector
  
  for (i in seq_along(vec)) {
    name <- vec[i]
    if (name_counts[name] > 1) {  # If name is duplicated
      name_index[name] <- name_index[name] + 1
      unique_vec[i] <- paste0(name, ".dup", name_index[name])
    }
  }
  
  return(unique_vec)
}

################################################################################
# import datasets obtained via r2
################################################################################

import.r2= function(filename, path=file.path(config$data_root,"bulkrna_patients"), id.col=NULL, dataset.id=paste_("untitled",timestamp(""))){
    
lines <- readLines(file.path(path, filename))   
 
# Split into two arrays
fcat("reading file...extracting metadata and gene expression...")
metadata_lines <- lines[grepl("^#", lines)]    # Lines that start with '#'
non_comment_lines <- lines[!grepl("^#", lines)] 

fcat("formatting metadata...")
    
metadata.df <- as.data.frame(do.call(rbind, strsplit(sub('#','',metadata_lines), "\\t")), stringsAsFactors = FALSE)

metadata.df=format.r2.meta(metadata.df) 
gc()

fcat("formatting counts...")
counts<- do.call(rbind, strsplit(non_comment_lines, "\\t"))

counts<- as.data.frame(counts, stringsAsFactors = FALSE)
counts=counts[2:nrow(counts), , drop=FALSE]  %>% dplyr::select(-2) 


rownames(counts)=make.unique.names(counts$V1)
counts = counts %>% select(-1)

counts=apply(counts, c(1,2), as.numeric)

if(!is.null(id.col)){
fcat("bringing names from metadata")
metadata.df<-metadata.df %>% col2names(., id.col) 
colnames(counts)=rownames(metadata.df)
}


rm(non_comment_lines)
rm(metadata_lines)
gc()
fcat("Done.")
outlist=list(counts=counts, metadata=metadata.df, dsname=dataset.id)
}


has.survival.data=function(analysis.list){
    lapply(c('Event_Free_Survival', 'Event_Free_Survival_Time_in_Days', 'Overall_Survival', 'Overall_Survival_Time_in_Days'), function(x)
 ifelse( length(analysis.list$gpcdf[,x ] %>% is.na %>% table %>% names)<2, 0,1)
) %>% Reduce(c,.)
}


################################################################################
#signature surivival analysis after calculating singatures
################################################################################
format.survival=function(analysis.list, sig, survivaltype, lowthresh, highthresh){
    

survmat=analysis.list$gpcdf


survmat=survmat %>% as.data.frame %>% dplyr::mutate(
                         Overall_Survival_Time_in_Days=as.integer(Overall_Survival_Time_in_Days),
                         Event_Free_Survival_Time_in_Days=as.integer(Event_Free_Survival_Time_in_Days),
                         Overall_Survival_Time_in_Years=Overall_Survival_Time_in_Days/365,
                         Event_Free_Survival_Time_in_Years=Event_Free_Survival_Time_in_Days/365,
                         Survival=!!sym(survivaltype), 
                         Time=!!sym(paste0(survivaltype, "_Time_in_", "Years")),
                         COG.Risk.Group=ifelse(COG.Risk.Group=="High Risk","High Risk", "Low-intermediate risk")
                         )

survmat[,"percentile"]=  classify.percentiles(analysis.list$gpcdf[, sig], lowthresh, highthresh)
percentile.column.name=paste0("percentile.", sig)
survmat[, percentile.column.name]=survmat[, "percentile"]

survmat}

signature.survival.analysis=function(analysis.list,sigs, lowthresh, highthresh, survivaltype="Overall_Survival", model=T, class.col="classify", adjustvar=NULL, exclude.med=T){ 

    plotlab=list("Overall_Survival"="OS", "Event_Free_Survival"="EFS")
    sca=2.5  
    allsurvs=lapply(sigs, function(sig){
  
    fcat("calculating survival dataset", analysis.list$dsname, "signature", sig)
  
      
        
survmat<-format.survival(analysis.list, sig, survivaltype, lowthresh, highthresh)
 
 
 

#for ease, track back the name of the percentile variable
percentile.column.name=paste0("percentile.", sig)


################################################################################
#generating percentiles for the signature
################################################################################

################################################################################
# generate combined variable stratifier
################################################################################
#variable specifying what specific stratum it is
fcat("0")
survmat[, "def"]=1
varname=ifelse(is.null(adjustvar), "def", adjustvar )
varstring=    paste0("stratum.", sig,'.',varname)
survmat=survmat %>% dplyr::mutate(temp=paste( !!sym(percentile.column.name),!!sym(varname) , sep="_"))
survmat=survmat %>% dplyr::mutate(!!sym(varstring):=temp)


#variable to be called generically called stratum within the code for practicality
survmat=survmat %>% dplyr::mutate(stratum=temp)

fcat("1")
if(class.col=="classify"){
survmat[ "classify"]= survmat[,percentile.column.name]

survmat[, "classification.variable"]=percentile.column.name
}else{
survmat[,"classify"]=  survmat[, class.col] 
survmat[, "classification.variable"]=class.col
}

################################################################################
# filtering
################################################################################

 fcat("size before filtering", nrow(survmat), ncol(survmat))
 

 cols.to.check=c(class.col, "gender", adjustvar, "Survival","Time")
 
 survmat=survmat[complete.cases(survmat[, cols.to.check]), ]
 
 fcat("size after filtering", nrow(survmat), ncol(survmat))



if(exclude.med){
    survmat=survmat %>% dplyr::filter(classify!="med")
}

obs=survmat %>% group_by(classify, inss_stage, MYCN_status, COG.Risk.Group, age.processed) %>% summarise(counts=n())
fcat("observation counts:")
print(obs)
fcat("3")


survmat$classify= factor(survmat$classify, levels=c("low", "high"))
survmat$classify=relevel(survmat$classify, ref="low")
#diagnose.cox.covariates(survmat, predictor="classify", covariates=cols.to.check)

fcat("4")


survtest=survival::survdiff(formula = survival::Surv(time=Time, event=as.numeric(Survival)) ~ classify,
data = survmat)



 
 
 
pval=p_value <- 1 - pchisq(survtest$chisq, df = length(survtest$n) - 1)
pstring= paste("\n p=", pval %>% sprintf("%.3f", .))




fcat("6")



if(class.col=="stratum"){
    fcat("7. plotting survival by stratum")

    survplot=ggsurvplot(fit = survival::Surv(time=Time, as.numeric(Survival)) ~ classify, 
    data = survmat, color="classify", risk.table=T, break.time.by=2 ,title=paste(plotlab[[survivaltype]],"strata", pstring))
}else{
    fcat("7.5")
    

 fcat("survmat dims", nrow(survmat), ncol(survmat))

survplot=explicit.survplot.classify(survmat, sig, pstring, survivaltype)
 }
fcat("8")
compbars=ggplot(survmat, aes(x=factor(classify, levels=c("low","high"))))+geom_bar(aes(fill=MYCN_status), position="fill")+theme_classic()+coord_flip()+NoLegend()

plot.title.string=paste_(survivaltype, analysis.list$dsname, "signature",sig, "by",class.col, "high", highthresh, "low", lowthresh)

if(class.col=="stratum"){
    fcat("8.5")
height.factor=1+ 0.03*(length(survmat[[class.col]] %>% table)-2)   
width.factor=0.7+ 0.01*(max(nchar(survmat[[class.col]]) %>% table)-4) 
}else{
 height.factor=1 
 width.factor=0.7
}


tpdf(path=paste0(create.storage.path(timestamp(paste_('survival_curves', analysis.version)), analysis.list$dsname,survivaltype,sig, paste0(lowthresh,"_", highthresh)), '/'), paste_("survival_plots_", plot.title.string), he=pw*sca*height.factor,wi=pw*sca*.7)
print(survplot)
dev.off()

coxforest=NULL
invariant.adjustvar=NULL
invariant.gender=NULL
complete.formula=NULL

covariates=vector()
flags=vector()
        

for(j in 1:length(adjustvar)){
if(!is.null(adjustvar[j])){
if(length(unique(analysis.list$gpcdf[[adjustvar[j]]]))==1){
 flags[j]=F
 
}else{flags[j]=T
covariates=c(covariates,paste0("+ ", adjustvar[j]))

}
}
}    
covariates=paste(covariates[flags] %>% removenas, collapse="")

if(model){
    
    fcat("9. performing model...")
library(survminer)    
    
    
        fcat("10. preparing formula...")

        coreformula=cox.formula <- paste("Surv(time=Time, event=as.numeric(Survival)) ~", class.col)

        

invariant.adjustvars=!flags %>% givename(., adjustvar)
        
        
if(!any(flags)){
    

fcat("9.5. Warning. with covariates. Proceeding to model only predictor coefficients")
    
complete.formula=coreformula


}else{
    
    ## perform this code when variables to adjust are clear of problems
    fcat("10. ")
    complete.formula=paste(coreformula, covariates)
    
    
}

fcat(complete.formula)
    cox.formula <- as.formula(complete.formula)

 
    fcat("11")
coxmodel <- coxph(cox.formula, data = survmat)


    fcat("12. plotting adjusted survival curve")


coxadjusted=make.adjustedplot(analysis.list, coxmodel, survmat, survivaltype, lowthresh, topthresh, sig, pstring, plot.title.string, width.factor) 

fcat("13 making coxforest")

coxforest=manual.ggforest(coxmodel, title=analysis.list$dsname, analysis.list=analysis.list, analysis.version=analysis.version, survivaltype=survivaltype, sig=sig, lowthresh=lowthresh, highthresh=topthresh, width.factor=width.factor, sca=sca, plot.title.string=plot.title.string)

}else{
 coxforest=NULL 
 coxmodel=NULL}


l=list(survival.data=survmat, survival.curve=survplot, test=survtest, counts=obs, cox.model=ifelse(model, list(model=coxmodel, forest.plot= coxforest, formula=complete.formula), NA), invariant.adjustvar=invariant.adjustvars)
l
}) %>% givename(., sigs)

pathh=create.storage.path(timestamp(paste_('survival_curves', analysis.version)), analysis.list$dsname,survivaltype)    
saveRDS(allsurvs, file= paste0(pathh,"/survival_analysis_output.rds"))    
allsurvs
}



explicit.survplot.classify= function(survmat, sig, pstring, survivaltype) {
    
  survplot = ggsurvplot(
    fit = survminer::surv_fit(Surv(time = Time, event = as.numeric(Survival)) ~ classify, 
                            data = survmat),
    color = "classify",
    palette = allcolors[["signature.group"]],
    risk.table = TRUE,
    break.time.by = 2,
    title = paste(plotlab[[survivaltype]], sig, pstring)
  )
  return(survplot)
}


make.adjustedplot<- function(analysis.list, coxmodel, 
  survmat, survivaltype, lowthresh, highthresh, 
  sig, pstring, plot.title.string, width.factor, sca
) {
  # Run ggadjustedcurves safely
  coxadjusted <- tryCatch({
    ggadjustedcurves(
      coxmodel,
      data = survmat %>% dplyr::filter(classify != "med"),
      variable = "classify",
      method = "marginal",
      color = "classify",
      palette = allcolors[["signature.group"]],
      risk.table = TRUE,
      break.time.by = 2,
      title = paste(plotlab[[survivaltype]], sig, pstring)
    )
  }, error = function(e) {
    message("Error in ggadjustedcurves: ", e$message)
   e$message
  })

  # Plot only if successful
  if (!is.character(coxadjusted)) {
    tpdf(
      path = paste0(
        create.storage.path(
          timestamp(paste_("survival_curves", analysis.version)),
          analysis.list$dsname,
          survivaltype,
          sig,
          paste0(lowthresh, "_", highthresh)
        ), '/'
      ),
      paste_("adjusted_KMplot", analysis.list$dsname,survivaltype, sig, analysis.version ),
      he = pw * sca,
      wi = pw * sca * width.factor
    )
    print(coxadjusted)
    dev.off()
  }
  
  coxadjusted
}


calculate.signatures=function(analysis.list, signatures){   
datasetid=digest::digest(analysis.list)
ds0=analysis.list$dsname
simpleCache(make.cache.name("gpcdf_with_signatures", analysis.version, package.vars(ds0, datasetid)),{    
#################################################################################
# calculating signatures
################################################################################

gsva_output=GSVA::gsva(analysis.list$vsd %>% apply(., c(1,2), as.numeric),
                         signatures,
                         method = "ssgsea")

signature.scores=gsva_output %>% t 

sigs=signature.scores %>% colnames
for(mv in sigs){
analysis.list$gpcdf[, mv]= signature.scores[, mv] %>% unname # incorporating 
}

analysis.list$signatures=sigs
gc()
analysis.list
}, assignToVar="analysis.list", reload=T)
analysis.list
}










make.cox.forest <- function(coxmodel, title="Hazard Ratio", plot.title.string) {
  
  coxforest <- tryCatch({
    plt=ggforest(coxmodel, data = survmat)
  }, error = function(e) {
    message("ggforest failed, using manual.ggforest instead.")
    plt=manual.ggforest(coxmodel)
  })
  
  tpdf(
    path = paste0(
      create.storage.path(timestamp(paste_('survival_curves', analysis.version)), 
                           analysis.list$dsname, survivaltype, sig, 
                           paste0(lowthresh, "_", highthresh)), 
      '/'
    ),
    paste_("coxforestplot_", plot.title.string),
    he = pw * sca,
    wi = pw * sca * width.factor
  )
  
  print(plt)
  dev.off()
  plt
}







manual.ggforest <- function(coxmodel, title = "Hazard Ratios", analysis.list, analysis.version, survivaltype, sig, lowthresh, highthresh, width.factor, sca, plot.title.string) {
    library(broom)    # for tidy()
  
  # 1. Extract tidy summary of model
  cox_summary <- broom::tidy(coxmodel, exponentiate = TRUE, conf.int = TRUE)
  
  # 2. Clean it up
  cox_summary <- cox_summary %>%
    filter(term != "(Intercept)") %>%  # if intercept exists, drop it
    mutate(term = factor(term, levels = rev(term)))  # reverse order for plotting
  
  # 3. Basic plot
  plt=ggplot(cox_summary, aes(x = estimate, y = term)) +
    geom_point(size = 3, color = "blue") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "black") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    theme_minimal(base_size = 14) +
    xlab("Hazard Ratio (log scale)") +
    ylab("") +
    scale_x_log10() +  # hazard ratios are usually on a log scale
    ggtitle(title) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(hjust = 1)
    )
  
    tpdf(
    path = paste0(
      create.storage.path(timestamp(paste_('survival_curves', analysis.version)), 
                           analysis.list$dsname, survivaltype, sig, 
                           paste0(lowthresh, "_", highthresh)), 
      '/'
    ),
    paste_("coxforestplot_", plot.title.string),
    he = pw * sca,
    wi = pw * sca * width.factor
  )
  
  print(plt)
  dev.off()
  plt
}




diagnose.cox.covariates <- function(data, predictor, covariates, time = "Time", event = "Survival") {
  library(survival)

  for (covar in covariates) {
    cat("\nChecking variable:", covar, "\n")
    
    # Check if variable exists
    if (!covar %in% colnames(data)) {
      cat("  ❌ Variable does not exist in the dataset.\n")
      next
    }

    var_data <- data[[covar]]
    
    # Check for only one unique value
    unique_vals <- unique(na.omit(var_data))
    if (length(unique_vals) == 1) {
      cat("  ⚠️ Variable has only one unique value:", unique_vals, "\n")
    }
    
    # Check for NAs
    na_count <- sum(is.na(var_data))
    if (na_count > 0) {
      cat("  ⚠️ Variable has", na_count, "missing values (NAs).\n")
    }
    
    # Check association with predictor
    if (is.factor(var_data) || is.character(var_data)) {
      tab <- table(data[[predictor]], var_data)
      if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
        cat("  ⚠️ Perfect association (zero cells in contingency table with predictor).\n")
      }
    } else {
      # For numeric variables, check if variance exists in each level of predictor
      group_var <- tapply(var_data, data[[predictor]], function(x) var(na.omit(x)))
      if (any(sapply(group_var, function(x) is.na(x) || x == 0))) {
        cat("  ⚠️ No variance in covariate within one or more levels of the predictor.\n")
      }
    }
    
    # Check for perfect separation in survival model
    formula <- as.formula(paste0("Surv(", time, ", ", event, ") ~ classify +", covar))
    fit <- try(coxph(formula, data = data), silent = TRUE)
    if (inherits(fit, "try-error")) {
      cat("  ❌ Cox model failed to fit. Possible perfect separation or other modeling issue.\n")
    } else if (any(is.na(coef(fit)))) {
      cat("  ⚠️ Model coefficients contain NA — might be due to singularities or collinearity.\n")
    }
  }
}



update.analysis<- function (analysis.list) {  
    dsname = analysis.list$dsname  
    package.id = digest::digest(analysis.list)  
    analysis.list$package.id = package.id   
    simpleCache(make.cache.name("full_analysis_components_r2",         analysis.version, package.vars(dsname, package.id)),         {            analysis.list        }, assignToVar = "analysis.list", reload = T)  
    
    analysis.list}



subset.analysis<- function(analysis.list, varname, alternatives, invert=F) {
  # Extract dataframe and matrix from the list
  df <- analysis.list$gpcdf
  mat <- analysis.list$vsd
  
  # Filter dataframe based on 'varname' column
  if(!invert){
      invertstring=""
  filtered_df <- df %>% dplyr::filter(!!sym(varname) %in% alternatives)
  }else{
invertstring="not_"      
   filtered_df <- df %>% dplyr::filter(!(!!sym(varname) %in% alternatives))   
  }
  # Subset matrix by selecting columns corresponding to the row names of filtered dataframe
  filtered_mat <- mat[, (colnames(mat) %in% rownames(filtered_df)), drop = FALSE]
  
  # Return modified list
  analysis.list$vsd=filtered_mat
  analysis.list$gpcdf=filtered_df
 
  newname=paste0(c(varname,"_", invertstring, alternatives), collapse="")
  ##editing the dataset name but making sure the editing happens only once
analysis.list$dsname=paste0(gsub(newname, "", analysis.list$dsname), newname)
analysis.list
}



################################################################################
# create storage folder chain, fully customiseable
################################################################################
create.storage.path <- function(..., base_dir="~/mnt_out/plots/") {
  # Define base directory
  
  
  # Construct the full path
  full_path <- file.path(base_dir, ...)
  
  # Create the directories (recursive = TRUE ensures all subdirectories are created)
  dir.create(full_path, recursive = TRUE, showWarnings = FALSE)
  
  # Return the full path
  return(normalizePath(full_path, mustWork = FALSE))
}


is.column.all.na <- function(df, colname) {
  # Check if the column exists
  if (!colname %in% names(df)) {
    stop("Column not found in dataframe.")
  }
  
  # Check if all values in the column are NA
  return(all(is.na(df[[colname]])))
}

check.variable.categories=function(df, var){
 df[[var]]  %>% table %>% length
}





fcat("V")
# Function to get the next available snapshot number
get.next.snapshot_number <- function(cache_prefix, cache_dir = tempdir()) {
    # List all existing cache files in the directory
    existing_caches <- list.files(cache_dir, pattern = paste0("^", cache_prefix, "_[0-9]+\\.RData$"))
    
    # Extract numbers from filenames
    snapshot_numbers <- as.numeric(gsub(paste0(cache_prefix, "_([0-9]+)\\.RData"), "\\1", existing_caches))
    
    # Find the next available number
    next_number <- ifelse(length(snapshot_numbers) == 0, 1, max(snapshot_numbers, na.rm = TRUE) + 1)
    
    return(next_number)
}

# Function to create a numbered snapshot of an object
save.snapshot <- function(object, cache_prefix, cache_dir = tempdir()) {
    # Get the next available snapshot number
    snapshot_number <- get_next_snapshot_number(cache_prefix, cache_dir)
    
    # Generate the cache name
    cache_name <- paste0(cache_prefix, "_", snapshot_number)
    
    # Save snapshot using simpleCache
    simpleCache(cache_name, { object }, cacheDir = cache_dir, saveOnExit = TRUE)
    
    return(cache_name)
}


################################################################################
# collect all self made functions and store them into a script
################################################################################


collect.functions <- function(filename = "loaded_functions.R") {
  # Create the R directory if it doesn't exist
  dir.create("./R", showWarnings = FALSE, recursive = TRUE)

  # Full path for saving
  filepath <- file.path("./R", filename)

  # Identify user-defined functions in the global environment
  user_functions <- ls(envir = globalenv())[sapply(ls(envir = globalenv()), function(x) is.function(get(x, envir = globalenv())))]

  # Save functions to a human-readable R script
  sink(filepath)
  for (f in user_functions) {
    cat("\n", f, "<- ", deparse(get(f)), "\n", sep = "")
  }
  sink()

  message("Functions saved to: ", filepath)
}

fcat("W")


cache.object <- function(object, id=NULL, analysis.version=analysis.version, recreate = FALSE) {
  # Get cache directory from global options
  
    if(is.null(id)){
        id=getvarname(object)
    }
    cache_dir <- getOption("RCACHE.DIR", default = tempdir())  # Default to tempdir if not set
  
  # Ensure the directory exists
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  # Generate a unique identifier using object's digest
  obj_hash <- digest(object, algo = "md5")  # MD5 hash of the object
  
  # Construct filename and full path
  file_identifier <- paste0(id, "_", obj_hash, "_", analysis.version)
  file_name <- paste0(file_identifier, ".RData")
  file_path <- file.path(cache_dir, file_name)
  
  # Check if file already exists
  if (file.exists(file_path) && !recreate) {
    message("File already exists and 'recreate' is FALSE. Skipping save.")
  } else {
    # Save the object to disk
    save(object, file = file_path)
    message("Object cached: ", file_path)
  }
  
  # Return file details
  return(list(
    path = file_path,
    identifier = file_identifier
  ))
}



hardpaste=function(x) paste(x, collapse="")


################################################################################
# print funciton args
################################################################################
print.function.args <- function(func_name) {
       # Check if the function exists in the current environment
       if (!exists(func_name, mode = "function")) {
         stop(paste("Function", func_name, "does not exist."))
       }
       
       # Get the function object
       func <- get(func_name, mode = "function")
       
       # Extract the formals (argument list)
       args <- formals(func)
       
       # Convert to a string representation
       arg_str <- paste0(names(args), 
                         ifelse(sapply(args, is.symbol), "", 
                                paste0(" = ", sapply(args, deparse))), 
                         collapse = ", ")
       
       # Print the function declaration line
       cat(sprintf("%s <- function(%s)\n", func_name, arg_str))
} 
     
     
summarisename= Vectorize(function(nm, total.print.length=30, first.words=4){
  #this function takes the first "first" words of the name, then if this exceeds total.print.length, 
  #we slice the middle section. 
  nm2=strsplit(nm, split=" ")[[1]]
  nm=nm2[1:(ifelse(first.words>length(nm2), length(nm2), first.words))] %>% paste(., collapse=" ")
  
  if(nchar(nm)>total.print.length){
    
    extend=round(total.print.length/2)
    a=substr(nm, 1, extend)
    b=substr(nm, nchar(nm)-extend, nchar(nm))
    out=paste0(a,"...", b)
    
  }else{
    out=nm
  }
  return(out)
}, USE.NAMES=F)


################################################################################
# cellinfo segregate cells within groups
################################################################################


################################################################################
# generate color variations around a base color
################################################################################
color.variations <- function(x, z = 10, m = 10) {
  # Helper to clamp values between 0 and 255
  clamp <- function(v) pmin(pmax(v, 0), 255)
  
  # Remove leading '#' if present
  x <- gsub("^#", "", x)
  
  # Parse RGB components
  r <- strtoi(substr(x, 1, 2), 16L)
  g <- strtoi(substr(x, 3, 4), 16L)
  b <- strtoi(substr(x, 5, 6), 16L)
  
  variations <- replicate(m, {
    # Sample deviations
    delta <- sample(-z:z, 3, replace = TRUE)
    
    # Apply and clamp
    r_new <- clamp(r + delta[1])
    g_new <- clamp(g + delta[2])
    b_new <- clamp(b + delta[3])
    
    # Convert back to hex
    sprintf("#%02X%02X%02X", r_new, g_new, b_new)
  })
  
  return(variations)
}

################################################################################
# sort colors by brightness
################################################################################

sort.by.sum <- function(colors) {
  # Remove "#" if present
  hex_clean <- gsub("^#", "", colors)
  
  # Extract R, G, B components and convert to decimal
  rgb_values <- t(sapply(hex_clean, function(hex) {
    r <- strtoi(substr(hex, 1, 2), base = 16)
    g <- strtoi(substr(hex, 3, 4), base = 16)
    b <- strtoi(substr(hex, 5, 6), base = 16)
    c(r, g, b)
  }))
  
  # Compute sum of RGB for each color
  rgb_sums <- rowSums(rgb_values)
  
  # Order colors by RGB sum
  colors_sorted <- colors[order(rgb_sums)]
  
  return(colors_sorted)
}


sort.by.brightness <- function(colors, descending = FALSE) {
  # Remove "#" if present
  hex_clean <- gsub("^#", "", colors)
  
  # Extract R, G, B components
  rgb_matrix <- t(sapply(hex_clean, function(hex) {
    r <- strtoi(substr(hex, 1, 2), 16L)
    g <- strtoi(substr(hex, 3, 4), 16L)
    b <- strtoi(substr(hex, 5, 6), 16L)
    c(r, g, b)
  }))
  
  # Calculate brightness
  brightness <- 0.2126 * rgb_matrix[, 1] + 
                0.7152 * rgb_matrix[, 2] + 
                0.0722 * rgb_matrix[, 3]
  
  # Sort by brightness
  sorted_indices <- order(brightness, decreasing = descending)
  return(colors[sorted_indices])
}

################################################################################
# cellinfo filter cells
################################################################################


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
  
    


################################################################################
# rebalance groups in a data frame
################################################################################



rebalance.by.group <- function(df, group_var, balance_var , return.ids.list=F) {
df=df %>% names2col(., "cellid")
  mins=df  %>%
    group_by(!!sym(group_var), !!sym(balance_var)) %>%
    mutate(nn = n()) %>%
    ungroup() %>%
    group_by(!!sym(group_var)) %>%
    summarise(min_n = min(nn)) %>% as.data.frame
  
  balanced=lapply(1:nrow(mins), function(x){
    
    localmin=mins[x, "min_n"]
    localclus=mins[x, group_var]
    df %>% dplyr::filter(!!sym(group_var)== localclus) %>% group_by(!!sym(balance_var)) %>%
    slice_sample(n = localmin, replace=F) %>%
    ungroup() 
}) %>% Reduce(rbind, .) 
  if(return.ids.list){
   
    idlist=balanced %>% group_split(!!sym(group_var)) %>% lapply(., function(x) x %>% pull(cellid )) 
    idlist.names=balanced %>% group_split(!!sym(group_var)) %>% lapply(., function(x) x %>% pull(!!sym(group_var)) %>% unique) %>% Reduce(c, .)
    idlist=idlist %>% givenames(., idlist.names)
    return(idlist)
    
  }else{
  return(balanced %>% as.data.frame %>% col2names(., "cellid") %>% select(-cellid))
    }
    
}
     
     
################################################################################
#segregate cells in a cellinfo object
################################################################################

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


cellinfo.segregate.cells=segregate.cells
################################################################################
# seriate cells with seurat
################################################################################

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



################################################################################
# save cellbygene h5ad
################################################################################

save.cellxgene=function(so, name="singlecellobject", path=paste0(config$out_root), reduction.use="umap", overwrite=F){
  
  
  ##cell by gene does not work well with numerical seurat clusters
  so@meta.data = so@meta.data %>% dplyr::mutate(seurat_clusters=paste0("cluster_", seurat_clusters))
  ##cellxgene can only use the default umap projection as the one shown on the platform so we rename the desired one
  if(reduction.use!="umap"){
    so@reductions$umap=so@reductions[[reduction.use]]
    colnames(so@reductions$umap@cell.embeddings)=c("UMAP_1", "UMAP_2")
  }
  
  savename=file.path(path, paste0(name, ".h5Seurat"))
SaveH5Seurat(so, filename = savename, overwrite=overwrite )
Convert(savename, dest = "h5ad", overwrite = overwrite )
}








################################################################################
# save cellbygene h5ad
################################################################################

# SeuratDisk::Convert only understands legacy Assay (v3/v4) objects.
# Assay5 (Seurat v5) must be joined and cast, or h5ad exports lack layers.
.assay5.to.assay <- function(assay_obj, assay_name = "assay") {
  if (!inherits(assay_obj, "Assay5")) {
    return(assay_obj)
  }
  assay_obj <- tryCatch(
    SeuratObject::JoinLayers(assay_obj),
    error = function(e) assay_obj
  )
  counts_mat <- tryCatch(
    SeuratObject::LayerData(assay_obj, layer = "counts"),
    error = function(e) NULL
  )
  data_mat <- tryCatch(
    SeuratObject::LayerData(assay_obj, layer = "data"),
    error = function(e) NULL
  )
  if (!is.null(counts_mat)) {
    new_assay <- CreateAssayObject(counts = counts_mat)
    if (!is.null(data_mat) &&
        identical(dim(data_mat), dim(counts_mat)) &&
        identical(rownames(data_mat), rownames(counts_mat)) &&
        identical(colnames(data_mat), colnames(counts_mat))) {
      new_assay <- SetAssayData(new_assay, slot = "data", new.data = data_mat)
    }
    return(new_assay)
  }
  # Fallback: coerce class if counts layer was unavailable
  tryCatch(
    as(assay_obj, Class = "Assay"),
    error = function(e) {
      stop(
        "Could not convert Assay5 '", assay_name,
        "' to legacy Assay for SeuratDisk export: ", conditionMessage(e)
      )
    }
  )
}


# Subset an assay to a feature set so SeuratDisk can place it in h5ad layers
# (Convert only copies other assays to layers when dims match X).
.subset.assay.features <- function(assay_obj, features, layer_slot = "data") {
  features <- features[features %in% rownames(assay_obj)]
  if (length(features) == 0) {
    return(NULL)
  }
  counts_mat <- tryCatch(
    GetAssayData(assay_obj, slot = "counts")[features, , drop = FALSE],
    error = function(e) NULL
  )
  data_mat <- tryCatch(
    GetAssayData(assay_obj, slot = "data")[features, , drop = FALSE],
    error = function(e) NULL
  )
  if (is.null(counts_mat) && is.null(data_mat)) {
    return(NULL)
  }
  if (is.null(counts_mat)) {
    new_assay <- CreateAssayObject(counts = data_mat)
  } else {
    new_assay <- CreateAssayObject(counts = counts_mat)
    if (!is.null(data_mat)) {
      new_assay <- SetAssayData(new_assay, slot = "data", new.data = data_mat)
    }
  }
  # SeuratDisk writes other assays to layers from scale.data (if present) else
  # data. For scvi layer="RNA" we usually want raw counts in that layer slot.
  if (identical(layer_slot, "counts") && !is.null(counts_mat)) {
    new_assay <- SetAssayData(new_assay, slot = "data", new.data = counts_mat)
  }
  new_assay
}

# Features SeuratDisk will use for AnnData X when converting `assay_name`.
# When scale.data exists it becomes X (typical SCT export: HVGs only), not all
# assay rownames — layers must match this feature set exactly.
.h5ad.convert.x.features <- function(so, assay_name) {
  assay_obj <- so[[assay_name]]
  scale_mat <- tryCatch(
    GetAssayData(assay_obj, assay = assay_name, slot = "scale.data"),
    error = function(e) NULL
  )
  if (!is.null(scale_mat) && length(scale_mat) > 0 && nrow(scale_mat) > 0) {
    return(rownames(scale_mat))
  }
  data_mat <- tryCatch(
    GetAssayData(assay_obj, assay = assay_name, slot = "data"),
    error = function(e) NULL
  )
  if (!is.null(data_mat) && length(data_mat) > 0 && nrow(data_mat) > 0) {
    return(rownames(data_mat))
  }
  rownames(assay_obj)
}

.h5ad.convert.x.source <- function(so, assay_name) {
  assay_obj <- so[[assay_name]]
  scale_mat <- tryCatch(
    GetAssayData(assay_obj, assay = assay_name, slot = "scale.data"),
    error = function(e) NULL
  )
  if (!is.null(scale_mat) && length(scale_mat) > 0 && nrow(scale_mat) > 0) {
    return("scale.data")
  }
  "data"
}

#' Export a Seurat object to h5ad via SeuratDisk for cellxgene / scvi
#'
#' SeuratDisk's h5ad Convert maps the Convert `assay` to AnnData `X` and maps
#' *other* assays to `adata.layers` **only if** they share the same feature
#' dimensions as `X`. When `scale.data` is present on the Convert assay (typical
#' for SCT), **X uses only those scaled features** (e.g. 5000 HVGs), not every
#' row in the assay. Other assays must be subset to that same feature set.
#'
#' @param so Seurat object
#' @param name Output basename (writes name.h5Seurat and name.h5ad)
#' @param path Output directory
#' @param reduction.use Reduction copied to `umap` for cellxgene display
#' @param overwrite Overwrite existing files
#' @param assay Assay written to AnnData `X` (default: DefaultAssay). Other
#'   assays become layers when feature dims match. For scvi with layer="RNA",
#'   keep this as an SCT/processed assay so `RNA` is exported as a layer.
#' @param ensure.rna.layer If TRUE and an RNA assay exists, subset RNA to the
#'   features that SeuratDisk will write to `X` so it can export as layers["RNA"].
#' @param rna.layer.slot Slot copied into RNA's `data` for layer export: `"counts"`
#'   (default, for scvi) or `"data"`.
save.cellxgene2=function(so, name="singlecellobject", path=paste0(config$out_root),
                        reduction.use="umap", overwrite=F, assay=NULL,
                        ensure.rna.layer=TRUE, rna.layer.slot=c("counts", "data")){
  rna.layer.slot <- match.arg(rna.layer.slot)

  ##cell by gene does not work well with numerical seurat clusters
  so@meta.data = so@meta.data %>% dplyr::mutate(seurat_clusters=paste0("cluster_", seurat_clusters))
  ##cellxgene can only use the default umap projection as the one shown on the platform so we rename the desired one
  if(reduction.use!="umap"){
    so@reductions$umap=so@reductions[[reduction.use]]
    colnames(so@reductions$umap@cell.embeddings)=c("UMAP_1", "UMAP_2")
  }

  # --- Seurat v5 / SeuratDisk compatibility ---------------------------------
  # Convert every Assay5 to legacy Assay so SaveH5Seurat/Convert can see
  # counts/data slots and emit adata.layers.
  for (assay_name in Assays(so)) {
    so[[assay_name]] <- .assay5.to.assay(so[[assay_name]], assay_name)
  }

  convert_assay <- assay %||% DefaultAssay(so)
  if (!convert_assay %in% Assays(so)) {
    stop("assay='", convert_assay, "' not found in Seurat object")
  }
  DefaultAssay(so) <- convert_assay
  x_source <- .h5ad.convert.x.source(so, convert_assay)
  x_features <- .h5ad.convert.x.features(so, convert_assay)
  n_assay_features <- nrow(so[[convert_assay]])

  if (length(x_features) < n_assay_features) {
    fcat(
      "AnnData X will use ", length(x_features), " ", x_source,
      " features from '", convert_assay, "' (assay has ", n_assay_features,
      " total features). Layers must match this set."
    )
  }

  # SeuratDisk only adds *other* assays to layers when nrow matches X exactly.
  if (isTRUE(ensure.rna.layer) && "RNA" %in% Assays(so) && convert_assay != "RNA") {
    rna_aligned <- .subset.assay.features(
      so[["RNA"]],
      x_features,
      layer_slot = rna.layer.slot
    )
    if (is.null(rna_aligned)) {
      warning(
        "Could not align RNA to AnnData X features from '", convert_assay,
        "'; layers['RNA'] may be missing in the h5ad."
      )
    } else if (nrow(rna_aligned) != length(x_features)) {
      warning(
        "RNA shares only ", nrow(rna_aligned), "/", length(x_features),
        " AnnData X features from '", convert_assay,
        "'. SeuratDisk requires identical dims for layers; RNA layer may be skipped."
      )
      so[["RNA"]] <- rna_aligned
    } else {
      so[["RNA"]] <- rna_aligned
      fcat(
        "Aligned RNA to", length(x_features),
        "features matching AnnData X (", convert_assay, " ", x_source,
        ") so it can export as layers['RNA']"
      )
    }
  }

  if (convert_assay == "RNA") {
    warning(
      "Convert assay is 'RNA', so RNA will be AnnData X, not layers['RNA']. ",
      "scvi.setup_anndata(layer='RNA') expects a layer; prefer assay=<SCT assay> ",
      "with ensure.rna.layer=TRUE (default)."
    )
  }

  savename=file.path(path, paste0(name, ".h5Seurat"))
  SaveH5Seurat(so, filename = savename, overwrite=overwrite)
  Convert(savename, dest = "h5ad", overwrite = overwrite, assay = convert_assay)
  invisible(file.path(path, paste0(name, ".h5ad")))
}





filter.marker.sets <- function(seurat_obj, marker_list) {
  # Get all genes in the Seurat object
  genes_in_obj <- rownames(seurat_obj)
  
  # Keep only marker sets that overlap with Seurat genes
  filtered_list <- lapply(marker_list, function(gene_set) {
    intersect(gene_set, genes_in_obj)
  })
  
  # Remove empty gene sets
  filtered_list <- filtered_list[lengths(filtered_list) > 0]
  
  return(filtered_list)
}
 
################################################################################
# improved module score wrapper around seurat's AddModuleScore function that gets rid of small lists and surgically removes pesky numbers
################################################################################

AddModuleScore3=function(object, features,rcache=NULL,label="", ...){
  if(is.null(rcache)){
    rcache=paste0("modulescore_metadata_id_", Project(object),label, "_", length(features), "groups")  %>% addversion
  }
  
  features=filter.marker.sets(object, features)
  fcat("total number of viable sets to test:", length(as.list(features)))
  
  simpleCache(rcache, {
  features=features[lapply(features, function(x) if(length(x)<1){F}else{T}) %>% Reduce(c, .)]
  names(features)=paste0("mscore_group_",names(features), "_s") #add an artificial flag to be later removed, thereby patching AddModuleScore's bizarre name modification
  for(k in names(features)){
  if(k %in% colnames(object@meta.data)){ 
    warning(cat(k, "is already present in object. Previous will be replaced"))
   object@meta.data= object@meta.data %>% dplyr::select(-!!sym(k))
    }
  }
  object=AddModuleScore(object, features, name=names(features))

  #AddModuleScore adds numbers to each signature's name so we find those new names
modnames=paste0(names(features), 1:length(features))
fcat("cleaning up ")

allcols=colnames(object@meta.data)

which.to.change=(allcols %in% modnames)
allcols[which.to.change]=gsub("_s[0-9]+$", "", allcols[which.to.change] )

#clean.names=gsub("mscore_group_", "", names(allcols[which.to.change]))

#if the cluster names are numbers after trimming, then leave the ms_cluster prefix in.
#if(!is.number.like(clean.names)){
#  allcols[which.to.change]=clean.names
#}

colnames(object@meta.data)=allcols
object@meta.data[, allcols[which.to.change],drop=F]
}, assignToVar="metamat", reload=T)
object=AddMetaData(object, metamat)
object
}    


################################################################################
# functions for data import
################################################################################


add.row.suffix=function(df, suf="-1"){
  
  rn=rownames(df)
  rn=paste0(rn, suf)
  
  rownames(df)=rn
  
}

import.dataset=function(folder){
CreateSeuratObject(Read10X_h5(file.path(config$data_root, folder, "filtered_feature_bc_matrix.h5")), project=strsplit(folder, split="/")[[1]][2])
}

import.clones=function(folder){
fread(file.path(config$data_root, folder, "clones.txt")) %>% as.data.frame %>% dplyr::mutate(cell_id=paste0(cell_id, "-1"), clone_nr_ds=paste_(strsplit(folder, split="/")[[1]][2],clone_nr )) %>% col2names(., "cell_id") 
}


################################################################################
# set markers from a collection of options
################################################################################

setmarkers<- function(so, markerset=1){
  so$RNA@misc$markers=so$RNA@misc[[markerset]]$markers
  so$RNA@misc$top_markers=so$RNA@misc[[markerset]]$top_markers
  so
}



plot.cluster.sizes <- function(seurat_obj, cluster_var = "seurat_clusters") {
  # Check if the clustering variable exists
  if (!cluster_var %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Clustering variable", cluster_var, "not found in metadata."))
  }
  
  # Count cells per cluster
  cluster_counts <- seurat_obj@meta.data %>%
    group_by(!!sym(cluster_var)) %>% 
    summarise(n=n()) %>%
    arrange(desc(n)) %>%
    mutate(cluster_var := factor(!!sym(cluster_var), levels = !!sym(cluster_var)))  # preserve order
  
  # Plot
  ggplot(cluster_counts, aes(x = !!sym(cluster_var), y = n)) +
    geom_bar(stat = "identity", fill = "#4C72B0") +
    labs(x = "Cluster", y = "Number of cells", title = "Cell count per cluster") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}



is.numeric.like <- function(x) {
  all(!is.na(suppressWarnings(as.numeric(x))))
}

get.module.name=function(x) paste0("mscore_group_", x)




cellinfo.order.groups.natural=function(cellinfo){
cellinfo$cell.list=cellinfo$cell.list[cellinfo$cell.list %>% names %>% mixedorder]
cellinfo$markerlist=cellinfo$markerlist[cellinfo$markerlist %>% names %>% mixedorder]
refresh.cellinfo(cellinfo)
}



keep.cellgroups.pattern=function(cellinfo, pattern){
 cellinfo$cell.list = cellinfo$cell.list[grepl(pattern, names(cellinfo$cell.list))] 
refresh.cellinfo(cellinfo)
 }


################################################################################
# subset specific clusters
################################################################################

subset.clusters=function(so, clusvar, clusters){
  
 so[ , so@meta.data %>% filter(!!sym(clusvar) %in% clusters) %>% rownames]  
  
}


################################################################################
# get function defaults for quick debugging
################################################################################

get.function.defaults <- function(fun) {
  # If input is a string, get the function by name
  if (is.character(fun)) {
    fun <- get(fun, mode = "function", envir = parent.frame())
  }

  # Get formals (arguments and defaults)
  formals_list <- formals(fun)

  # Replace missing default values with NULL
  defaults <- lapply(formals_list, function(x) {
    if (is.symbol(x)) return(NULL)
    eval(x)
  })

  return(defaults)
}

################################################################################
# function to retrieve the loaded genes
################################################################################

get.loaded.genes=function(so, npcs=5, ngenes=5){
lapply(1:npcs, function(x) Loadings(so) %>% as.data.frame %>% arrange(-!!sym(paste0("PC_",x))) %>% head(ngenes) %>% rownames) %>% Reduce(c, .)
}
################################################################################
# workflow t make a tidy heatmap out of a seurat object and its respective cellinfo
################################################################################
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




################################################################################
# calculate signature averages per group. 
################################################################################

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





join_meta_exp2 <- function(
    seurat_obj,
    genes,
    assay = "SCT",
    layer = "data", 
    cells=NULL
) {
  
  if(is.null(cells)){
   cells=colnames(seurat_obj) 
  }
  # ensure assay exists
  if (!assay %in% Assays(seurat_obj)) {
    stop(paste("Assay", assay, "not found in Seurat object"))
  }
  # ensure genes exist
  genes_present <- genes[genes %in% rownames(seurat_obj[[assay]])]
  genes_missing <- setdiff(genes, genes_present)
  if (length(genes_missing) > 0) {
    warning("These genes were not found and will be skipped: ",
            paste(genes_missing, collapse = ", "))
  }
  if (length(genes_present) == 0) {
    stop("None of the requested genes were found in the assay")
  }
  # extract expression matrix
  expr <- GetAssayData(seurat_obj, assay = assay, layer = layer)[genes_present, , drop = FALSE]
  # transpose to cells x genes
  expr_df <- as.data.frame(t(as.matrix(expr)))
  # add to metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = expr_df)
  return(seurat_obj@meta.data[cells, ])
}





################################################################################
# For each cell, calculate what percentage of the markers of each signatutre are expressed.
################################################################################

AddMarkerPercentages <- function(seurat_obj, marker_list, assay = NULL, slot = "data") {
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }
  
  # Extract expression matrix
  expr_mat <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  for (marker_name in names(marker_list)) {
    markers <- marker_list[[marker_name]]
    
    # Keep only markers present in data
    markers_present <- intersect(markers, rownames(expr_mat))
    
    if (length(markers_present) == 0) {
      warning(paste("No markers found in data for", marker_name))
      seurat_obj[[marker_name]] <- NA
      next
    }
    
    # Expression logical matrix: TRUE if expressed (>0)
    expr_logical <- expr_mat[markers_present, , drop = FALSE] > 0
    
    # Calculate % expressed per cell
    percent_expr <- Matrix::colSums(expr_logical) / length(markers_present) * 100
    
    # Add to metadata
    seurat_obj[[paste0("pct_group_",marker_name)]] <- percent_expr
  }
  
  return(seurat_obj)
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


################################################################################
#
################################################################################


marker.umap<- function(seurat_obj, 
                            markers, 
                            colorby, 
                            neighbors = 30, 
                            min_dist = 0.3, 
                            repulsion_strength = 1, mscorevar=NULL) {
  
  # Ensure markers are present in the SCT assay
  markers_present <- intersect(markers, rownames(seurat_obj[["RNA"]]))
  if (length(markers_present) == 0) {
    stop("None of the markers are present in the RNA assay.")
  }
  
  # Create new assay with only these markers
  new_assay <- seurat_obj[["SCT"]][markers_present, ]
  seurat_obj[["marker_assay"]] <- CreateAssayObject(
                                                    data = new_assay)
  
  DefaultAssay(seurat_obj) <- "marker_assay"
  
  # Run PCA & UMAP
  seurat_obj <- ScaleData(seurat_obj, features = markers_present, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, features = markers_present, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:min(length(markers_present), 10),
                        n.neighbors = neighbors,
                        min.dist = min_dist,
                        repulsion.strength = repulsion_strength,
                        verbose = FALSE)
  
  # Make UMAP plot
  if(colorby=="seurat_clusters"){cols= adjust.seurat.colors(seurat_obj)}else{ allcolors[[colorby]]}
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = colorby) + 
    ggtitle(paste("UMAP on markers:", paste(markers_present, collapse = ", ")))+scale_color_manual(values=cols)
  if(!is.null(mscorevar)){
  p2 <- FeaturePlot(seurat_object, reduction="umap", features=mscorevar)
  }else{p2<-FeaturePlot(seurat_object, reduction="umap", features=markers_present[1])}
  # Cleanup to free memory
  gc()
  
  return(p)
}


signature_expression_by_percentile <- function(seurat_obj, signature_genes, var, nbreaks = 5) {
  
  # Keep only genes present in object
  signature_genes <- intersect(signature_genes, rownames(seurat_obj))
  if (length(signature_genes) == 0) stop("None of the signature genes are present in the object.")
  
  # Extract metadata variable
  meta_var <- seurat_obj[[var]][, 1]
  if (!is.numeric(meta_var)) stop("The variable must be numeric to calculate percentiles.")
  
  # Create percentile bins
  breaks <- quantile(meta_var, probs = seq(0, 1, length.out = nbreaks + 1), na.rm = TRUE)
  bin_labels <- paste0("P", 1:nbreaks)
  seurat_obj$percentile_bin <- cut(meta_var, breaks = breaks, labels = bin_labels, include.lowest = TRUE)
  
  # Get expression data (logical: > 0 = expressed)
  expr_data <- GetAssayData(seurat_obj, slot = "data")[signature_genes, ]
  expr_binary <- expr_data > 0
  
  # Calculate percentage of expressing cells per bin & gene
  results <- data.frame()
  for (gene in signature_genes) {
    tmp <- data.frame(
      gene = gene,
      bin = bin_labels,
      pct_expr = sapply(bin_labels, function(b) {
        cells_in_bin <- WhichCells(seurat_obj, expression = percentile_bin == b)
        if (length(cells_in_bin) == 0) return(0)
        mean(expr_binary[gene, cells_in_bin]) * 100
      })
    )
    results <- rbind(results, tmp)
  }
  
  # Plot stacked barplot
  p <- ggplot(results, aes(x = gene, y = pct_expr, fill = bin)) +
    geom_bar(stat = "identity", position = "stack") +
    ylab("% Cells Expressing") + xlab("Gene") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}


library(Seurat)

add_zero_cells <- function(seurat_obj, n_cells = 100, cell_prefix = "zeroCell",rerun=F ) {
  
  # Get gene names and create matrix of zeros
  gene_names <- rownames(seurat_obj)
  zero_mat <- Matrix::Matrix(0, nrow = length(gene_names), ncol = n_cells, sparse = TRUE)
  rownames(zero_mat) <- gene_names
  colnames(zero_mat) <- paste0(cell_prefix, "_", seq_len(n_cells))
  
  # Add new cells to the raw counts
  counts_combined <- cbind(GetAssayData(seurat_obj, slot = "counts"), zero_mat)
  seurat_obj <- SetAssayData(seurat_obj, slot = "counts", new.data = counts_combined)
  
  # Add dummy metadata
  zero_meta <- data.frame(row.names = colnames(zero_mat), dummy = TRUE, seurat_clusters="mock_zero")
  meta_combined <- rbind(seurat_obj@meta.data, zero_meta)
  meta_combined$dummy[is.na(meta_combined$dummy)] <- FALSE
  seurat_obj@meta.data <- meta_combined
  
  
  # Re-run transformations so new cells are included
  
  if(rerun){
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
  }
  return(seurat_obj)
}



plot_marker_expression_fraction <- function(seurat_obj, markers, group_var, colorlist=allcolors) {
  # Ensure grouping variable exists in metadata
  if (!(group_var %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Grouping variable", group_var, "not found in Seurat metadata."))
  }
  
  # Extract expression matrix for markers
  expr <- FetchData(seurat_obj, vars = c(markers, group_var))
  
  # Reshape to long format
  expr_long <- expr %>%
    pivot_longer(cols = all_of(markers), names_to = "gene", values_to = "expr") %>%
    mutate(expr_binary = expr > 0) # TRUE/FALSE expression indicator
  
  # Compute % of cells expressing per group per gene
  summary_df <- expr_long %>%
    group_by(!!sym(group_var), gene) %>%
    summarise(
      pct_expressing = 100 * sum(expr_binary) / n(),
      .groups = "drop"
    )
  
  # Plot
  p <- ggplot(summary_df, aes(x = gene, y = pct_expressing, fill = !!sym(group_var))) +
    geom_col(position = "dodge") +
    theme_minimal(base_size = 14) +
    ylab("% cells expressing") +
    xlab("Marker gene") +
    #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values=colorlist[[group_var]])+facet_wrap(~gene, scales="free_x", nrow=5)
  
  return(p)
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



################################################################################
#n relabel catrgories
################################################################################

reset.categories<- function(seurat_obj, var_name, new_var_name = NULL) {
  # Extract variable from metadata
  var <- seurat_obj@meta.data[[var_name]]
  
  # Ensure it's a factor
  var <- as.factor(var)
  
  # Sort the levels and reassign them to 1...n
  new_labels <- setNames(seq_along(sort(levels(var))), levels(var))
  
  # Apply mapping
  var_numeric <- as.character(new_labels[as.character(var)])
  
  # If new_var_name not provided, overwrite
  if (is.null(new_var_name)) {
    seurat_obj@meta.data[[var_name]] <- var_numeric
  } else {
    seurat_obj@meta.data[[new_var_name]] <- var_numeric
  }
  
  return(seurat_obj)
}


################################################################################
# smooth based on neighbors so w don't have random cells flying in between.
################################################################################

relabel_by_neighbors <- function(seurat_obj, group_var, graph_name = "RNA_snn") {
  # Extract neighbor graph
  g <- seurat_obj@graphs[[graph_name]]
  if (is.null(g)) stop(paste("Graph", graph_name, "not found in Seurat object."))
  
  # Ensure grouping variable exists
  if (!(group_var %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Grouping variable", group_var, "not found in metadata."))
  }
  
  labels <- seurat_obj@meta.data[[group_var]]
  new_labels <- labels
  
  # Go through each cell
  for (i in seq_len(ncol(g))) {
    neighbors <- which(g[, i] > 0)  # neighbor indices
    if (length(neighbors) == 0) next
    
    neighbor_labels <- labels[neighbors]
    
    # Count frequencies
    tab <- table(neighbor_labels)
    max_count <- max(tab)
    top_labels <- names(tab[tab == max_count])
    
    # Only reassign if a unique majority exists
    if (length(top_labels) == 1 && top_labels != labels[i]) {
      new_labels[i] <- top_labels
    }
  }
  
  # Add new labels to metadata
  seurat_obj@meta.data[[paste0(group_var, "_smoothed")]] <- new_labels
  
  return(seurat_obj)
}


################################################################################
# mutual information between two genes
################################################################################


mutual_info_genes <- function(seurat_obj, gene_x, gene_y, assay = "RNA", slot = "data") {
  library(infotheo)
  # Extract expression matrix
  fcat("getting assay data")
  mat <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  fcat("checking if genes are present")
  # Check if genes exist
  if (!(gene_x %in% rownames(mat))) stop(paste("Gene", gene_x, "not found in Seurat object."))
  
  present.genes= gene_y[gene_y %in% rownames(mat)]
  if (length(present.genes) ==0) stop(paste("Genes", setdiff(gene_y, present.genes), "not found in Seurat object."))
  
  # Binarize expression (1 = expressed, 0 = not expressed)
   fcat("binarising matrix")
  mat2=mat[c(gene_x,present.genes), ]>0
  rownames(mat2)=c(gene_x,present.genes)
  mis=lapply(1:length(present.genes), function(gy){
    cat(".")
    if(gy %% 100 ==0){
     cat("\n") 
    }
  # Compute mutual information
  mi <- mutinformation(mat2[gene_x, ], mat2[present.genes[gy],])
  }) %>% Reduce(c,. ) %>% givename(., present.genes)
  
  gc()
  return(mis)
}

#########################################
pasteall=function(x) paste(x, collapse="_")


################################################################################
# cellinfo split clusters by condition
################################################################################


geomtype <- function(plot) {
  if (length(plot$layers) == 0) {
    return("no geoms")
  }
  # Extract the class of the geom object from the first layer
  gclass <- class(plot$layers[[1]]$geom)[1]
  # ggplot geoms are like "GeomPoint", "GeomBar", etc.
  # Strip the "Geom" prefix and lowercase
  gsub("^Geom", "", gclass) |>
    tolower()
}


################################################################################
# write cellranger-style outs folder from a seurat object
################################################################################


write.cellranger.outs <- function(seurat_obj, outdir) {
  # Load needed packages
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
  
  # Ensure output directory exists
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # Extract raw counts (sparse matrix preferred)
  counts <- Seurat::GetAssayData(seurat_obj, slot = "counts")
  
  # Barcodes
  barcodes <- colnames(counts)
  barcodes_file <- file.path(outdir, "barcodes.tsv.gz")
  writeLines(barcodes, gzfile(barcodes_file))
  
  # Features (10x format expects: gene_id, gene_name, feature_type)
  features <- rownames(counts)
  features_df <- data.frame(
    gene_id = features,
    gene_name = features,
    feature_type = "Gene Expression"
  )
  features_file <- file.path(outdir, "features.tsv.gz")
  write.table(
    features_df,
    file = gzfile(features_file),
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  
  # Matrix (sparse market matrix format)
  matrix_file_con <- gzfile(file.path(outdir, "matrix.mtx.gz"), "w")
  Matrix::writeMM(counts, file = matrix_file_con)
  close(matrix_file_con)
  message("CellRanger-like files written to: ", outdir)
  invisible(TRUE)
}

################################################################################
# from cluster object  
################################################################################
cluster.to.list <- function(clusters) {
  split(names(clusters), clusters)
}



################################################################################
#
################################################################################
filter_df <- function(df, vars, values, ops) {
  stopifnot(
    length(vars) == length(values),
    length(vars) == length(ops)
  )
  
  for (i in seq_along(vars)) {
    var <- vars[i]
    val <- values[i]
    op  <- ops[i]
    
    if (!op %in% c(">", ">=", "<", "<=")) {
      stop(paste("Invalid operator:", op))
    }
    
    df <- df[do.call(op, list(df[[var]], val)), ]
  }
  
  return(df)
}  





remove.pattern.columns <- function(df, pattern) {
  df[, !startsWith(names(df), pattern), drop = FALSE]
}




cut2 <- function(x, breaks, ..., include.lowest = TRUE, right = TRUE) {
  # Standard cut inside range
  f <- cut(x, breaks = breaks, include.lowest = include.lowest, right = right, ...)
  # Add outside categories
  f <- as.character(f)
  f[x < min(breaks)] <- "smaller"
  f[x > max(breaks)] <- "larger"
  # Make factor with ordered levels
  inside_levels <- levels(cut(breaks, breaks = breaks, include.lowest = include.lowest, right = right))
  levels_all <- c("smaller", inside_levels, "larger")
  factor(f, levels = levels_all, ordered = TRUE)
}




#load.markers=function(so, id=NULL, rcache=NULL, label=""){
 # m=
  
  
#}



################################################################################
#
################################################################################


group.proportion.barplot <- function(
    so,
    cluster_col,
    condition_col,
    rank_condition=NULL, colorlist=list(), height="n", #can also be proportion
    cluster_order=NULL
) {
  # Check inputs
  meta=so@meta.data
  stopifnot(
    cluster_col %in% colnames(meta),
    condition_col %in% colnames(meta)
  )
  
  # Summarize counts and proportions
  df_summary <- meta %>%
    dplyr::mutate(
      cluster = !!sym(cluster_col),
      condition = !!sym(condition_col)
    ) %>%
    dplyr::count(cluster, condition) %>%
    complete(cluster, condition, fill = list(n=0)) %>%
    group_by(cluster) %>% dplyr::mutate(proportion = n/sum(n)) %>% ungroup() %>% as.data.frame()
    
    print(df_summary)
    #group_by(cluster, condition) %>%
    #summarise(
     # total_cells = sum(n()),
    #) %>% group_by(cluster) %>% dplyr::mutate(proportion = total_cells/sum(total_cells)) %>% ungroup()
  #fcat("df_summary")
  #print(df_summary %>% as.data.frame)
if(is.null(rank_condition)){
  rank_condition=df_summary[[condition_col]][1]
}
  
  if(is.null(colorlist)){
    colorlist=make.colorlist(so@meta.data, col_names=c(condition_col, cluster_col))
  }  
  
  
  # Get ordering of clusters by proportion of chosen condition
  if(is.null(cluster_order)){
    cluster_order <- df_summary %>%
    dplyr::filter(condition == !!rank_condition) %>%
    arrange(desc(proportion)) %>%
    pull(cluster)
  }
  fcat(cluster_order)
  
  # Apply ordering
  df_summary <- df_summary %>%
    mutate(cluster = factor(cluster, levels = cluster_order))
  
  # Plot
  list(plot=ggplot(df_summary, aes(x = cluster, y = !!sym(height) , fill = condition))+geom_col()+scale_fill_manual(values=colorlist[[condition_col]]), 
       cluster_order=cluster_order, 
       df=df_summary
       )
}



pastemarkers=Vectorize(function(x){
  print(as.vector(x))
  result <- tryCatch(
    {
      if (as.vector(x) %in% names(ci$markerlist)) {
        ci$markerlist[[x]][
          1:ifelse(length(ci$markerlist[[x]]) >= 5, 5, length(ci$markerlist[[x]]))
        ] %>% 
          Reduce(pastec, .)
      } else {
        NA
      }
    },
    error = function(e) {
      NA
    }
  )
  
  result}, USE.NAMES=F)



plot.module.scores=function(so, modules=NULL , red=DefaultDimReduc(so), ...){
  
  if(is.null(modules)){
    modules=grep("mscore_group", colnames(so@meta.data), value=T) %>% mixedsort
  }
  
  ggfeature(so, features=modules, ...)
  }

################################################################################
#. AAA search for the value of a matrix in querycolumn when the value in refcol is lookup.value
################################################################################
lookup.synonym <- function(df, lookup.value, refcolumn, querycolumn){
  df[[querycolumn]][df[[refcolumn]] %in% lookup.value] %>% unique
  
}

lookup.synonym.human <- Vectorize(function(lookup.value){
  x=mouse.orthologs[["Gene.name"]][mouse.orthologs[["Mouse.gene.name"]] == lookup.value] %>% unique
  if(is.null(x)) return(NA) else return(x)
}, USE.NAMES=F)


lookup.synonym.mouse <- Vectorize(function(lookup.value){
  x=mouse.orthologs[["Gene.name"]][mouse.orthologs[["Mouse.gene.name"]] == lookup.value] %>% unique
  if(is.null(x)) return(NA) else return(x)
}, USE.NAMES=F)


coalesce_name <- function(x, y) {
  x <- ifelse(is.na(x) | x == "", NA, x)
  y <- ifelse(is.na(y) | y == "", NA, y)
  
  out <- ifelse(
    !is.na(x) & is.na(y), toupper(x),
    ifelse(
      is.na(x) & !is.na(y), toupper(y),
      ifelse(
        !is.na(x) & !is.na(y) & toupper(x) == toupper(y), toupper(x),
        toupper(y)
      )
    )
  )
  
  out <- ifelse(is.na(out), "noname", out)
  out
}





library(Seurat)
library(ggplot2)
library(patchwork)

run_umap_grid <- function(
    so,
    param_x, #. VALID PARAMETERS: c("neighbors", "min.dist", "dims", "repulsion.strength")
    param_y,
    nsteps = 5,
    param_ranges = list(
      neighbors = c(20, 100),
      min.dist = c( 0.30, 1),
      dims = c(10, 50),
      repulsion.strength = c( 1,5)
    ),
    reduction = "pca",
    colorby="seurat_clusters",
    defaults=list(
      neighbors = 50,
      min.dist =0.3,
      dims = 1:50,
      repulsion.strength = 1
    )
) {
  
  valid_params <- c("neighbors", "min.dist", "dims", "repulsion.strength")
  
  param.alias <- c(neighbors="nn", min.dist="md", dims=reduction, repulsion.strength="rs")
  stopifnot(param_x %in% valid_params)
  stopifnot(param_y %in% valid_params)
  stopifnot(param_x != param_y)
  
  # Generate parameter sequences
  x_vals <- seq(param_ranges[[param_x]][1],
                param_ranges[[param_x]][2],
                length.out = nsteps)
  
  y_vals <- seq(param_ranges[[param_y]][1],
                param_ranges[[param_y]][2],
                length.out = nsteps)
  
  plots <- list()
  
  for (i in seq_along(x_vals)) {
    for (j in seq_along(y_vals)) {
      
      params <- defaults
      
      params[[param_x]] <- x_vals[i]
      params[[param_y]] <- y_vals[j]
      
      reduction_name <- paste0(
        "umap_",
        param.alias[param_x], "_", round(x_vals[i], 2), "_",
        param.alias[param_y], "_", round(y_vals[j], 2)
      )
      
      so <- RunUMAP(
        so,
        reduction = reduction,
        dims = if (param_x == "dims" || param_y == "dims") 1:params$dims else params$dims,
        n.neighbors = params$neighbors,
        min.dist = params$min.dist,
        repulsion.strength = params$repulsion.strength,
        reduction.name = reduction_name,
        reduction.key = prepare.rk(reduction_name), 
        return.model=F
      )
      
      p <- ggumap(so, colorby=colorby, reduction=reduction_name, reduction.key=prepare.rk(reduction_name))+
        ggtitle(
          paste0(
            param_x, "=", round(x_vals[i], 2),
            ", ",
            param_y, "=", round(y_vals[j], 2)
          )
        ) +
        theme(plot.title = element_text(size = 8))
      
      plots[[reduction_name]] <- p
    }
  }
  
  grid_plot <- wrap_plots(plots, ncol = nsteps)
  print(grid_plot)
  
  return(list(reductions=so@reductions, plotgrid=grid_plot, plots=plots, nsteps=nsteps))
}



plot_umap_grid <- function(
    so, 
    umaplist, 
    colorby="seurat_clusters", colorlist=NULL){
  pholder=so@reductions
  so@reductions=umaplist$reductions[grepl("^umap", names(umaplist$reductions) )]

plots=lapply(names(so@reductions), function(reduction_name){
p <- ggumap(so, colorby=colorby, reduction=reduction_name, reduction.key=prepare.rk(reduction_name), colorlist=colorlist)+
  ggtitle(
    reduction_name 
  )+
  theme(plot.title = element_text(size = 8))

 p

}) %>% givename(names(so@reductions))

umaplist$plotgrid <- wrap_plots(plots, ncol = umaplist$nsteps)
so@reductions=pholder
umaplist$plots=plots
gc()

return(umaplist)
}



clean.umaplist=function(umaplist){
umaplist$plots=NULL
umaplist$plotgrid=NULL
umaplist   
}




plot.module.scores=function(so, modules=NULL , red=DefaultDimReduc(so)){
  
  if(is.null(modules)){
    modules=grep("mscore_group", colnames(so@meta.data), value=T)
  }
  lapply(FeaturePlot(so, features=modules, combine=F, reduction=red, pt.size=2, raster=T), function(x)  x+scale_color_viridis()+NoLegend()+NoAxes()) %>% Reduce('+', .)
}




################################################################################
# neighobrs and cluster processing
################################################################################

seuratprocess.clustering <- function(
    object,
    process.params=NULL,
    reduction=NULL,
    resolution=NULL,
    n.neighbors=NULL,
    clusterlabel=NULL,
    reprocess = FALSE
) {
  
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  
  if (!reduction %in% names(object@reductions)) {
    stop(paste0("Reduction '", reduction, "' not found in object"))
  }
  
  # Programmatic base name
  base_name <- paste0(
    reduction,
    "_nn", n.neighbors,
    "_res", gsub("\\.", "p", resolution)
  )
  
  graph_name <- paste0(base_name, "_graph")
  cluster_colname <- paste0(base_name, "_clusters")
  labeled_copy_name <- paste0("seuratclusters_", clusterlabel)
  
  # Detect existing content
  graph_exists <- graph_name %in% names(object@graphs)
  cluster_exists <- cluster_colname %in% colnames(object@meta.data)
  label_exists <- labeled_copy_name %in% colnames(object@meta.data)
  
  # Abort logic
  if (!reprocess && (graph_exists || cluster_exists || label_exists)) {
    message("Existing graph/cluster/label detected. Skipping processing.")
    return(object)
  }
  
  # If forced reprocess, remove existing artifacts
  if (reprocess) {
    
    # Remove graphs with matching base_name or label
    graphs_to_remove <- grep(
      pattern = paste0(base_name, "|", clusterlabel),
      x = names(object@graphs),
      value = TRUE
    )
    
    if (length(graphs_to_remove) > 0) {
      object@graphs[graphs_to_remove] <- NULL
    }
    
    # Remove metadata columns with matching base_name or label
    cols_to_remove <- grep(
      pattern = paste0(base_name, "|", clusterlabel),
      x = colnames(object@meta.data),
      value = TRUE
    )
    
    if (length(cols_to_remove) > 0) {
      object@meta.data[, cols_to_remove] <- NULL
    }
  }
  
  # Run neighbor graph
  object <- FindNeighbors(
    object = object,
    reduction = reduction,
    dims = 1:ncol(Embeddings(object, reduction)),
    k.param = n.neighbors,
    graph.name = graph_name,
    verbose = FALSE
  )
  
  # Run Leiden clustering
  object <- FindClusters(
    object = object,
    graph.name = graph_name,
    resolution = resolution,
    algorithm = 4,  # Leiden
    cluster.name = cluster_colname,
    verbose = FALSE
  )
  
  # Create labeled copy of seurat_clusters
  object[[labeled_copy_name]] <- object[[cluster_colname]][, 1]
  
  message("Completed clustering: ", cluster_colname)
  
  return(object)
}



clean.module.scores=function(so){
  so@meta.data=so@meta.data %>% dplyr::select(-contains("mscore_group"))
   so}

################################################################################
# add embeding to a seurat object
################################################################################

add.embedding <- function(
    seurat_obj,
    file,
    reduction_name = "importedembedding",
    assay = NULL,
    key = NULL
) {
  # Load required package
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required")
  }
  
  # Read CSV (assumes first column = cell names / rownames)
  df <- read.csv(file, row.names = 1, check.names = FALSE)
  
  # Ensure cells match
  common_cells <- intersect(rownames(df), colnames(seurat_obj))
  
  if (length(common_cells) == 0) {
    stop("No matching cell names between CSV and Seurat object")
  }
  
  # Subset and reorder to match Seurat object
  df <- df[common_cells, , drop = FALSE]
  df <- df[colnames(seurat_obj), , drop = FALSE]
  
  # Convert to matrix
  emb <- as.matrix(df)
  
  # Default key (prefix for dimensions like custom_1, custom_2)
  if (is.null(key)) {
    key <- paste0(reduction_name, "_")
  }
  
  # Create DimReduc object
  dr <- Seurat::CreateDimReducObject(
    embeddings = emb,
    key = key,
    assay = if (is.null(assay)) Seurat::DefaultAssay(seurat_obj) else assay
  )
  
  # Add to Seurat object
  seurat_obj[[reduction_name]] <- dr
  
  return(seurat_obj)
}



################################################################################
#
################################################################################

load.markers=function(so, rcache, label="imported_markers", minfc=0.5, minrate=0.5){
  minfc=0.5
  minrate=0.5
  simpleCache(rcache, assignToVar="m", reload=T )
  mp=m %>% group_split(group1) %>% lapply(., function(x) x %>% filter(log_fc>=minfc, padj<=0.05, rate1>=minrate) %>% arrange(-log_fc)) %>% Reduce(rbind,. )
  mb=m %>% group_split(group1) %>% lapply(., function(x) x %>% filter(log_fc<=-minfc, padj<=0.05) %>% arrange(log_fc)) %>% Reduce(rbind,. )
  so$RNA@misc[[label]]$markers=m
  so$RNA@misc[[label]]$top_markers=mp
  so$RNA@misc[[label]]$negative_markers=mb
  so
}


################################################################################
# function specific for all micropatterns
################################################################################
annotate.treatment <- function(seurat_obj, column) {
  vals <- seurat_obj@meta.data[[column]]
  
  result <- sapply(vals, function(x) {
    x <- as.character(x)
    
    # Rule 1: FGF2
    if (grepl("FGF2", x, ignore.case = TRUE)) {
      return("FGF2")
    }
    
    # Rule 2: IWP conditions
    if (grepl("noIWP", x, ignore.case = TRUE)) {
      return("noIWP")
    }
    if (grepl("plusIWP", x, ignore.case = TRUE)) {
      return("plusIWP")
    }
    
    # Rule 3: RA / DAPT combinations
    has_RA <- grepl("_RA", x, ignore.case = F)
    has_DAPT <- grepl("DAPT", x, ignore.case = F)
    
    if (has_RA && has_DAPT) {
      return("RA_DAPT")
    } else if (has_RA) {
      return("RA")
    } else if (has_DAPT) {
      return("DAPT")
    }
    
    # Rule 4: fate labels
    if (grepl("neural_crest", x, ignore.case = TRUE)) {
      return("neural_crest")
    }
    if (grepl("neural_plate", x, ignore.case = TRUE)) {
      return("neural_plate")
    }
    if (grepl("placode", x, ignore.case = TRUE)) {
      return("placode")
    }
    if (grepl("all_fates", x, ignore.case = TRUE)) {
      return("all_fates")
    }
    
    # Default if nothing matches
    return("BMP4")
  })
  
  return(unname(result))
}





add.pattern.flag <- function(seurat_obj, colname, pattern, new_col=NULL) {
  # Extract metadata
  meta <- seurat_obj@meta.data
  
  # Create dynamic column name
  if(is.null(new_col)){
  new_col <- paste0("has.", pattern)
  }
  # Add column
  meta <- meta %>%
    mutate(!!new_col := stringr::str_detect(.data[[colname]], pattern ))
  
  # Assign back
  seurat_obj@meta.data <- meta
  
  return(seurat_obj)
}
################################################################################
#
################################################################################

FindSubclusters <- function(
    seurat_obj,
    groupvar,
    groupid,
    resolution.clusters = 0.5,
    new_colname = NULL,
    reduction.for.neighbors="pca",
    ...
) {
  # Check inputs
  if (!groupvar %in% colnames(seurat_obj@meta.data)) {
    stop("groupvar not found in Seurat metadata.")
  }
  
  # Default new column name
  if (is.null(new_colname)) {
    new_colname <- paste0(groupvar, "_subclustered")
  }
  
  # Identify cells to subset
  cells.use <- colnames(seurat_obj)[seurat_obj@meta.data[[groupvar]] == groupid]
  
  if (length(cells.use) == 0) {
    stop("No cells found for the specified groupid.")
  }
  
  # Subset object
  sub_obj <- subset(seurat_obj, cells = cells.use)
  fcat("finding neighbors")
  # Run clustering on subset
  
  sub_obj <- Seurat::FindNeighbors(sub_obj, dims=1:length(colnames(seurat_obj@reductions[[reduction.for.neighbors]])), reduction=reduction.for.neighbors, graph.name="sgraph")
  fcat("finding clusters. graph name:")
  fcat(sub_obj@graphs %>% names %>% Reduce(pastec, .))
  sub_obj <- Seurat::FindClusters(sub_obj, resolution = resolution.clusters, graph.name="sgraph")
  
  # Extract new cluster identities
  sub_clusters <- sub_obj@meta.data %>% dplyr::mutate(!!sym(new_colname):=paste0(groupid, "_sc_", seurat_clusters) )
  
  # Create new labels: e.g. "3_subcluster_1"
  #new_labels <- paste0(groupid, "_sc_", sub_clusters)
  #names(new_labels) <- colnames(sub_obj)
  #print(new_labels)
  # Initialize new metadata column as a copy of original groupvar
  seurat_obj[[new_colname]] <- as.character(seurat_obj@meta.data[[groupvar]])
  
  # Replace labels only for subset cells
  seurat_obj@meta.data[rownames(sub_clusters), new_colname] <- sub_clusters[,new_colname]
  
  return(seurat_obj)
}



jitter.colors <- function(colors, sca = 0.1, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Validate scale
  if (sca < 0 || sca > 1) {
    stop("sca must be between 0 and 1")
  }
  
  # Convert scale to RGB integer range
  max_jitter <- round(255 * sca)
  
  # Function to jitter one color
  jitter_one <- function(col) {
    
    # Remove leading #
    hex <- gsub("^#", "", col)
    
    # Extract RGB pairs
    r_hex <- substr(hex, 1, 2)
    g_hex <- substr(hex, 3, 4)
    b_hex <- substr(hex, 5, 6)
    
    # Hex -> decimal
    r <- strtoi(r_hex, base = 16)
    g <- strtoi(g_hex, base = 16)
    b <- strtoi(b_hex, base = 16)
    
    # Random jitter
    r <- r + sample(-max_jitter:max_jitter, 1)
    g <- g + sample(-max_jitter:max_jitter, 1)
    b <- b + sample(-max_jitter:max_jitter, 1)
    
    # Clamp to valid RGB range
    r <- max(0, min(255, r))
    g <- max(0, min(255, g))
    b <- max(0, min(255, b))
    
    # Decimal -> hex
    rgb(
      r = r,
      g = g,
      b = b,
      maxColorValue = 255
    )
  }
  
  # Apply to all colors
  vapply(colors, jitter_one, character(1))
}



check.cache= function(x) file.exists(file.path(config$cachedir, paste0(x, ".RData")))  
  





shared.marker.plot=function(so, gene1, gene2, assay=DefaultAssay(so), groupvars, groupvars2=NULL){
  library(gtools)
  
  
  genedat.wide=join_meta_exp2(so, genes=c(gene1, gene2), assay=assay) %>% as.data.frame
  
  bintilecounts=genedat.wide %>% dplyr::mutate(!!sym(gene1):=ifelse(!!sym(gene1)>0, 1, 0), !!sym(gene2):=ifelse(!!sym(gene2)>0,1,0 )) %>% group_by(across(all_of(c(groupvars,groupvars2, gene1, gene2)))) %>% summarise(counts=n(), mngene1=mean(!!sym(gene1)),mngene2=mean(!!sym(gene2)) ) %>% group_by(across(all_of(c(!!groupvars, !!groupvars2)))) %>% dplyr::mutate(counts.normalised= counts/sum(counts)) %>% ungroup %>% data.frame
  
  ordr=bintilecounts[[groupvars]] %>% unique %>% mixedsort
  
  
  if(!is.null(groupvars2)){
    plt=ggplot(bintilecounts, aes(x=!!sym(gene1), y=!!sym(gene2)))+geom_tile(aes(fill=counts.normalised))+facet_wrap(vars(factor(.data[[groupvars]], levels=ordr), !!sym(groupvars2)) )
  }else{
    plt=ggplot(bintilecounts, aes(x=!!sym(gene1), y=!!sym(gene2)))+geom_tile(aes(fill=counts.normalised))+facet_wrap(vars(factor(.data[[groupvars]], levels=ordr)) )
  }
  
  
  plt2=plt+
    scale_fill_viridis(name = paste0(groupvars, "\nfraction in group"))+
    theme_classic()+
    geom_text(aes(label = counts), color = "white", size = 3)+
    scale_x_continuous(breaks = c(0, 1), labels = c("0", ">0"))+
    scale_y_continuous(breaks = c(0, 1), labels = c("0", ">0"))
  
  tpdf(path=config$plotpath, paste0("tiles_counts2d_", paste(gene1,gene2, collapse="-")), wi=pw*sca*1.2, he=pw*sca*1.3)
  print(plt2)
  dev.off()
  
  list(df=bintilecounts, plot=plt2)
}

relevel.colors=function(allcolors){
  nms=names(allcolors)
  allcolors=lapply(nms, function(x){
    
    newnames=gtools::mixedsort(names(allcolors[[x]])) 
    allcolors[[x]][newnames]   
  }) %>% givename(., nms)
  
  allcolors
}


################################################################################
# generate a geom_point plot of markers on the basis of a cellinfo object
################################################################################

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


################################################################################
#. count the number of  cells within a range of rna counts. 
################################################################################

bin.counts <- function(seurat_obj, bin_size = 1000, max_bin=NULL, colname="bin.counts") {
  fcat("processing binned counts...")
  ncount <- seurat_obj@meta.data$nCount_RNA
  
  # Determine the max bin needed to cover the highest nCount value
  if(is.null(max_bin)){
    max_bin <- ceiling(max(ncount, na.rm = TRUE) / bin_size) * bin_size
  }
  # Create breaks: 0, bin_size, 2*bin_size, ..., max_bin
  breaks <- c(seq(0, max_bin, by = bin_size), Inf)
  #fcat(breaks)
  # Create labels like "<1000", "<2000", etc.
  labels <- c(paste0("<", breaks[-1]))
  #fcat(labels)
  # Bin nCount_RNA using cut()
  bin_var <- cut(
    ncount,
    breaks = breaks,
    labels = labels,
    right = TRUE,
    include.lowest = TRUE, 
    
  )
  
  # Add to metadata as an ordered factor (preserves bin order for plotting)
  seurat_obj@meta.data[[colname]] <- factor(bin_var, levels = labels, ordered = TRUE)
  
  return(seurat_obj)
}




get.library=Vectorize(function(x) expinfo$libraries[[x]], USE.NAMES=F)


################################################################################
# hox genes workflow
################################################################################


hox.workflow=function(so){  
  allhox= rownames(so)[so %>% rownames  %>% grepl("^HOX",.)]
  
  
  prephox= function(x) if(grepl('-', x)) return(substr(x, 4,8)) else( paste0(substr(x, 1,3), substr(x, 5,5), substr(x, 4,4)))
  alphahox= sapply(allhox, prephox) %>% sort %>% names
  boxhox=allhox %>% mixedsort
  
  boxes=c("A","B","C","D")
  purehox= boxhox[!grepl("-AS", boxhox)]
  purehoxletters= substr(purehox, 4,4)
  purehoxnumbers= substr(purehox, 5,7)
  purehox.orderednumbers=lapply(boxes, function(box){
    purehoxnumbers[purehoxletters==box] %>% as.numeric %>% sort
    
  }) %>% Reduce(c, .)
  
  
  purehox.clusteredbynumbers=lapply(1:13, function(b){
    purehox[purehoxnumbers==b]
    
  }) %>% Reduce(c, .)
  
  get.letter=Vectorize(function(x) substr(x, 4,4), USE.NAMES=F)
  get.number=Vectorize(function(x) substr(x, 5,6) %>% as.numeric, USE.NAMES=F)
  
  finalhox=paste("HOX", purehoxletters, purehox.orderednumbers, sep="")
  
  list(
    hoxtab= join_meta_exp2(so, genes=finalhox, assay=DefaultAssay(so)),
    hoxnames=finalhox)
}




#' Create rainbow color mappings for categorical columns in a data frame
#'
#' For each specified column name, finds all unique values (categories) and
#' assigns each a distinct color from `rainbow()`. Returns a list of named
#' color vectors, one per column, where the vector names are the category
#' values and the vector elements are hex color codes.
#'
#' @param df A data frame.
#' @param col_names Character vector of column names to process. If NULL,
#'   all columns in `df` are used.
#' @param sort_levels Logical. Whether to sort unique values before assigning
#'   colors (for reproducible/consistent ordering). Default TRUE.
#' @param na_color If not NULL, NA values found in a column will be assigned
#'   this color and included in the output vector (named "NA"). If NULL,
#'   NAs are dropped before computing uniques. Default NULL.
#' @param ... Additional arguments passed to `rainbow()` (e.g. `alpha`, `s`, `v`).
#'
#' @return A named list, one entry per column in `col_names`, each containing
#'   a named character vector of colors keyed by category.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   group = c("A", "B", "A", "C"),
#'   batch = c("x", "y", "x", "y")
#' )
#' color_list <- make_rainbow_palettes(df, col_names = c("group", "batch"))
#' color_list <- make_rainbow_palettes(df)  # all columns
#' }
make.colorlist=function(df,
                                  col_names = NULL,
                                  sort_levels = TRUE,
                                  na_color = NULL,
                                  ...) {
  
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame.")
  }
  
  # If col_names is NULL, use all columns in the data frame
  if (is.null(col_names)) {
    col_names <- colnames(df)
  }
  
  # Check for columns that don't exist in df
  missing_cols <- setdiff(col_names, colnames(df))
  if (length(missing_cols) > 0) {
    warning("The following columns were not found in 'df' and will be skipped: ",
            paste(missing_cols, collapse = ", "))
    col_names <- intersect(col_names, colnames(df))
  }
  
  if (length(col_names) == 0) {
    stop("No valid columns to process.")
  }
  
  result_list <- list()
  
  for (col in col_names) {
    values <- df[[col]]
    
    has_na <- any(is.na(values))
    unique_vals <- unique(values[!is.na(values)])
    
    # Coerce to character for consistent naming (handles factors, numerics, etc.)
    unique_vals <- as.character(unique_vals)
    
    if (sort_levels) {
      unique_vals <- mixedsort(unique_vals)
    }
    
    n <- length(unique_vals)
    
    if (n == 0) {
      warning("Column '", col, "' has no non-NA unique values; skipping.")
      next
    }
    
    colors <- rainbow(n, ...)
    names(colors) <- unique_vals
    
    # Optionally add a color for NA category
    if (has_na && !is.null(na_color)) {
      colors <- c(colors, setNames(na_color, "NA"))
    }
    
    result_list[[col]] <- colors
  }
  
  return(result_list)
}



find.cache<- function(folder="~/mnt_out/rcache/",
                              id,
                              type = NULL,
                              feeling.lucky = FALSE,
                              recursive = FALSE, 
                      extension=F) {

  if (!dir.exists(folder)) {
    stop("Folder does not exist: ", folder)
  }

  # Get all .RData files (case-insensitive extension match)
  all_files <- list.files(
    path = folder,
    pattern = "\\.RData$",
    ignore.case = TRUE,
    full.names = TRUE,
    recursive = recursive
  )

  if (length(all_files) == 0) {
    stop("No .RData files found in folder: ", folder)
  }

  # Filter by id (match against basename, not full path)
  matches <- all_files[grepl(id, basename(all_files), fixed = TRUE)]

  # Optionally filter further by type
  if (!is.null(type)) {
    matches <- matches[grepl(type, basename(matches), fixed = TRUE)]
  }

  if (length(matches) == 0) {
    stop("No .RData files found matching id = '", id,
         if (!is.null(type)) paste0("' and type = '", type) else "",
         "' in folder: ", folder)
  }

  if (length(matches) > 1 && !feeling.lucky) {
    fcat(
      "Multiple .RData files matched (", length(matches), " found):\n",
      paste("  -", matches, collapse = "\n"),
      "\nSet feeling.lucky = TRUE to automatically load the first match, ",
      "or refine 'id'/'type' to narrow the search."
    )
  }
  

   
ifelse(!extension, tools::file_path_sans_ext(basename(matches)), basename(matches))  
}
#   file_to_load <- matches[1]
# 
#   if (length(matches) > 1 && feeling.lucky) {
#     message("Multiple matches found, loading the first one (feeling.lucky = TRUE): ",
#             basename(file_to_load))
#   }
# 
#   # Load into a fresh environment to safely capture the object(s)
#   env <- new.env()
#   loaded_names <- load(file_to_load, envir = env)
# 
#   if (length(loaded_names) == 1) {
#     return(get(loaded_names, envir = env))
#   } else {
#     # If the RData file contains multiple objects, return them as a named list
#     result <- mget(loaded_names, envir = env)
#     message("Note: '", basename(file_to_load), "' contains ", length(loaded_names),
#             " objects (", paste(loaded_names, collapse = ", "),
#             "). Returning them as a named list.")
#     return(result)
#   }
# }


get.default.glasswork.params=function(){
  
      glasswork.params=list(
        size=1,
        fillvar="null",
        linecolor="#000000",
        concavity=2, 
        fillcols="#FFFFFF00"
      )
glasswork.params
      }

## calculate catagory purity on nearest neighborhoods

nnpurity <- function(seurat_obj, graph.name, col_name, categ) {
library(Seurat)
library(Matrix)
  # --- checks ---
  if (!graph.name %in% names(seurat_obj@graphs)) {
    stop(paste0("Graph '", graph.name, "' not found. Available graphs: ",
                paste(names(seurat_obj@graphs), collapse = ", ")))
  }
  if (!col_name %in% colnames(seurat_obj@meta.data)) {
    stop(paste0("Column '", col_name, "' not found in meta.data"))
  }

  # --- pull graph as sparse matrix (cells x cells) ---
  graph <- as(seurat_obj@graphs[[graph.name]], "CsparseMatrix")

  # --- align metadata to graph's cell order ---
  meta_vals <- seurat_obj@meta.data[[col_name]]
  names(meta_vals) <- rownames(seurat_obj@meta.data)
  meta_vals <- meta_vals[rownames(graph)]

  is_categ <- as.numeric(meta_vals == categ)

  # --- binarize (handles weighted SNN graphs too) and drop self-edges ---
  graph_bin <- graph
  graph_bin@x[graph_bin@x > 0] <- 1
  diag(graph_bin) <- 0

  # --- vectorized neighbor counts ---
  n_neighbors        <- Matrix::rowSums(graph_bin)
  n_categ_neighbors  <- as.numeric(graph_bin %*% is_categ)

  purity <- ifelse(n_neighbors > 0, 100 * n_categ_neighbors / n_neighbors, NA)

  # --- write back to metadata ---
  col_out <- paste0("nnpurity_", categ)
  seurat_obj[[col_out]] <- purity[colnames(seurat_obj)]  # re-align just in case

  return(seurat_obj)
}


################################################################################
#
################################################################################


#' Resample cells within groups of a Seurat object (with replacement)
#'
#' For each level of `group_var` in the object's metadata, cells are resampled
#' WITH replacement so that every group ends up with exactly `n` cells. Every
#' property of a sampled cell (all assays, all layers/slots: counts, data,
#' scale.data, and all metadata columns) is copied exactly -- only the cell
#' name is changed, since Seurat requires unique cell names and a cell drawn
#' more than once can't keep its original name.
#'
#' Note: dimensional reductions, graphs (e.g. neighbor graphs), and images are
#' NOT copied over. These are computed from the original cell set/relationships
#' and become invalid once cells are duplicated/resampled -- you should
#' recompute them (e.g. NormalizeData/FindVariableFeatures/RunPCA/RunUMAP,
#' FindNeighbors) on the returned object if you need them.
#'
#' @param seurat_obj A Seurat object.
#' @param group_var  Character. Name of a metadata column to group cells by.
#' @param n          Integer. Target number of cells per group after resampling.
#'                   Must be >= the size of the largest group.
#' @param copy.reductions Logical (default TRUE). If TRUE, for every existing
#'                   dimensional reduction (PCA, UMAP, etc.) the cell embedding
#'                   row of each sampled cell is duplicated and renamed to match,
#'                   so the returned object keeps those reductions populated for
#'                   the new (resampled) set of cells. Feature loadings (if any)
#'                   are copied through unchanged, since they don't depend on
#'                   individual cells. Note that these reductions are NOT
#'                   recomputed -- a duplicated cell simply gets the same
#'                   embedding coordinates it originally had; this is meant as a
#'                   convenience for visualization/lookup, not a substitute for
#'                   rerunning RunPCA/RunUMAP if you need reductions computed on
#'                   the resampled data itself.
#'
#' @return A new Seurat object with length(unique(groups)) * n cells.
resample.seurat <- function(seurat_obj, group_var, n, copy.reductions = TRUE) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The 'Seurat' package is required.")
  }
  library(Seurat)

  meta <- seurat_obj@meta.data

  if (!group_var %in% colnames(meta)) {
    stop(sprintf("'%s' not found in the object's metadata.", group_var))
  }

  # Collapse Seurat v5 multi-layer assays (e.g. split by sample) into single
  # layers so the rest of the function can treat v4 and v5 objects uniformly.
  seurat_obj <- tryCatch(
    SeuratObject::JoinLayers(seurat_obj),
    error = function(e) seurat_obj
  )

  meta <- seurat_obj@meta.data
  groups <- split(rownames(meta), meta[[group_var]], drop = TRUE)
  max_size <- max(lengths(groups))

  if (n < max_size) {
    stop(sprintf(
      "n = %d is smaller than the largest group ('%s', %d cells). n must be >= %d.",
      n, names(groups)[which.max(lengths(groups))], max_size, max_size
    ))
  }

  # --- Draw the resampled cell names (with replacement) per group ---
  sampled_old_names <- unlist(
    lapply(groups, function(cells) sample(cells, size = n, replace = TRUE)),
    use.names = FALSE
  )

  new_cell_names <- paste0("cell_", seq_along(sampled_old_names))

  # --- Duplicate metadata rows to match ---
  new_meta <- meta[sampled_old_names, , drop = FALSE]
  rownames(new_meta) <- new_cell_names

  # --- Helper: pull a layer/slot's matrix in a v4/v5-safe way ---
  get_layer_safe <- function(assay_obj, layer_name) {
    out <- tryCatch(
      SeuratObject::GetAssayData(assay_obj, layer = layer_name),
      error = function(e) NULL
    )
    if (is.null(out)) {
      out <- tryCatch(
        SeuratObject::GetAssayData(assay_obj, slot = layer_name),
        error = function(e) NULL
      )
    }
    out
  }

  assay_names <- SeuratObject::Assays(seurat_obj)
  default_assay <- DefaultAssay(seurat_obj)
  # Process the default assay first so CreateSeuratObject has counts to start from
  assay_names <- c(default_assay, setdiff(assay_names, default_assay))

  new_obj <- NULL

  for (assay_name in assay_names) {

    assay_obj <- seurat_obj[[assay_name]]

    layer_names <- tryCatch(
      SeuratObject::Layers(assay_obj),
      error = function(e) c("counts", "data", "scale.data")
    )
    if (is.null(layer_names)) layer_names <- c("counts", "data", "scale.data")

    layer_mats <- list()

    for (ln in layer_names) {
      mat <- get_layer_safe(assay_obj, ln)
      if (is.null(mat)) next
      if (is.null(dim(mat)) || ncol(mat) == 0) next

      # scale.data (and sometimes other layers) may only contain a subset of
      # cells/features (e.g. variable features only) -- only replicate the
      # sampled cells that are actually present in this particular layer.
      keep <- sampled_old_names %in% colnames(mat)
      if (!any(keep)) next

      mat_new <- mat[, sampled_old_names[keep], drop = FALSE]
      colnames(mat_new) <- new_cell_names[keep]
      layer_mats[[ln]] <- mat_new
    }

    if (length(layer_mats) == 0) next

    if (is.null(new_obj)) {
      # Prefer "counts" to seed the object; fall back to whatever is available
      counts_mat <- layer_mats[["counts"]]
      if (is.null(counts_mat)) counts_mat <- layer_mats[[1]]

      new_obj <- CreateSeuratObject(
        counts   = counts_mat,
        meta.data = new_meta,
        assay    = assay_name
      )

      for (ln in setdiff(names(layer_mats), "counts")) {
        new_obj <- tryCatch(
          SeuratObject::SetAssayData(new_obj, assay = assay_name, layer = ln, new.data = layer_mats[[ln]]),
          error = function(e) SeuratObject::SetAssayData(new_obj, assay = assay_name, slot = ln, new.data = layer_mats[[ln]])
        )
      }

    } else {
      new_assay <- CreateAssayObject(counts = layer_mats[["counts"]] %||% layer_mats[[1]])

      for (ln in setdiff(names(layer_mats), "counts")) {
        new_assay <- tryCatch(
          SeuratObject::SetAssayData(new_assay, layer = ln, new.data = layer_mats[[ln]]),
          error = function(e) SeuratObject::SetAssayData(new_assay, slot = ln, new.data = layer_mats[[ln]])
        )
      }

      new_obj[[assay_name]] <- new_assay
    }
  }

  DefaultAssay(new_obj) <- default_assay

  # --- Optionally replicate cell embeddings for every existing reduction ---
  if (isTRUE(copy.reductions)) {

    reduction_names <- SeuratObject::Reductions(seurat_obj)

    for (red_name in reduction_names) {

      red_obj <- seurat_obj[[red_name]]
      embeddings <- SeuratObject::Embeddings(red_obj)

      keep <- sampled_old_names %in% rownames(embeddings)
      if (!any(keep)) next

      emb_new <- embeddings[sampled_old_names[keep], , drop = FALSE]
      rownames(emb_new) <- new_cell_names[keep]

      # Feature loadings don't depend on individual cells, so carry them through as-is
      loadings <- tryCatch(SeuratObject::Loadings(red_obj), error = function(e) new(Class = "matrix"))
      proj_loadings <- tryCatch(SeuratObject::Loadings(red_obj, projected = TRUE), error = function(e) new(Class = "matrix"))

      assay_used <- tryCatch(SeuratObject::DefaultAssay(red_obj), error = function(e) default_assay)
      red_key <- tryCatch(SeuratObject::Key(red_obj), error = function(e) paste0(red_name, "_"))

      new_red <- SeuratObject::CreateDimReducObject(
        embeddings = emb_new,
        loadings   = loadings,
        projected  = proj_loadings,
        assay      = assay_used,
        key        = red_key
      )

      new_obj[[red_name]] <- new_red
    }
  }

  new_obj
}


# Small helper used above (avoids depending on rlang for `%||%`)
`%||%` <- function(a, b) if (is.null(a)) b else a




###############################################################################


plot.module.correlations=function(so, colorby="seurat_clusters", sca=5){


library(GGally)
library(dplyr)

mscore_cols <- grep("^mscore_group_", names(so@meta.data), value = TRUE) %>% mixedsort()

ggp=ggpairs(so@meta.data, columns = mscore_cols,
        # color points by a grouping variable if you have one
        mapping = aes(color = !!sym(colorby), alpha = 0.5),
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", size = 0.5)))

ggp %>% stdplot(.,title=paste_("mscore_corrplot", Project(so)) %>% addversion, sca=sca)

ggp
}



plot.module.densities=function(so, colorby="seurat_clusters", sca=1, colorlist=NULL, label=NULL){

if(is.null(colorlist)){
colorlist=make.colorlist(so@meta.data, col_names = colorby)
}
  
ms <- grep("^mscore_group_", names(so@meta.data), value = TRUE) %>% mixedsort()
    
ggp=ggplot(so@meta.data %>% pivot_longer(., cols=contains("mscore_group"), names_to="signature", values_to="score"), aes(x=score))+
  geom_density(aes(fill=factor(!!sym(colorby), levels=names(colorlist[[colorby]]))), alpha=0.3)+
  labs(fill = colorby)+
  facet_wrap(~factor(signature, levels=ms), scales="free")+
  scale_fill_manual(values=colorlist[[colorby]])

ggp %>% stdplot(.,title=paste_("densityplot_mscores_", paste0(ifelse(is.null(label), "", paste0(label, "_")), Project(so)), "coloredby", colorby) %>% addversion, sca=sca)

ggp
}


ggfeature=function(so, features, red=DefaultDimReduc(so), noaxes=T, nolegend=T, list=F, colorscale=NULL){
  

if(noaxes & nolegend){    
plts=FeaturePlot(so, features=features, combine = F, red=red) %>% lapply(., function(x) x+NoAxes()+NoLegend()+colorscale) 
}else{
  if(noaxes & !nolegend){    
plts=FeaturePlot(so, features=features, combine = F, red=red) %>% lapply(., function(x) x+NoAxes()+colorscale) 
}else{
   if(!noaxes & nolegend){    
plts=FeaturePlot(so, features=features, combine = F, red=red) %>% lapply(., function(x) x+NoLegend()+colorscale) 
   }else{
  
     if(!noaxes & !nolegend){    
plts=FeaturePlot(so, features=features, combine = F, red=red) %>% lapply(., function(x) x+colorscale) 
   }
     }

}
  }
  
  
  
if(list){
  
  return(plts)}else{
    return(plts %>% Reduce('+', .))
}}


################################################################################
# collect multi paths, tuned for this project
################################################################################

collect.multi.paths <- function(main_folder) {
  
  # Get all library folders within the main folder
  library_folders <- list.files(main_folder, full.names = TRUE, recursive = FALSE)
  h5_paths <- list()
  libs<-list()
  
  for (lib_folder in library_folders) {
    
    per_sample_outs <- file.path(lib_folder, "per_sample_outs")
  
    if (!dir.exists(per_sample_outs)) {
      warning(paste("No per_sample_outs folder found in:", lib_folder))
      next
    }
    
    # Get all sample folders within per_sample_outs
    sample_folders <- list.files(per_sample_outs, full.names = TRUE, recursive = FALSE)
    
    for (sample_folder in sample_folders) {
      
      h5_file <- file.path(sample_folder, "sample_filtered_feature_bc_matrix.h5")
      
      if (file.exists(h5_file)) {
        sample_name <- basename(sample_folder)
        h5_paths[[sample_name]] <- h5_file
        libs[[sample_name]]<-strsplit(lib_folder, split="/")[[1]] %>% get.last.element
      } else {
        warning(paste("h5 file not found in:", sample_folder))
      }
    }
  }
  
  return(list(h5=h5_paths, libraries=libs))
}





################################################################################
# Report builder: initiate / load / add / collect PDF figures for an Rmd report
#
# Workspace object: `reports` (list of labeled reports)
# reports[[label]] = list(
#   cache  = "report_<label>_<analysis.version>",
#   label  = <label>,
#   sections = list(<section>, ...)   # ordered figure/pdf entries
# )
################################################################################

# Resolve the analysis version used in report cache and folder names.
# Prefer the workspace analysis.version (set in Rmds after sourcing); else config.
.get.report.analysis.version <- function() {
  if (exists("analysis.version", envir = .GlobalEnv, inherits = FALSE)) {
    av <- get("analysis.version", envir = .GlobalEnv)
    if (!is.null(av) && !is.na(av) && nzchar(as.character(av)[1])) {
      return(as.character(av)[1])
    }
  }
  av <- config[["analysis.version"]]
  if (is.null(av) || is.na(av) || !nzchar(as.character(av)[1])) {
    stop("analysis.version is not set in the workspace or in config")
  }
  as.character(av)[1]
}

# NULL / empty -> "default"; must already be a valid make.names() identifier
.normalize.report.label <- function(label = NULL) {
  if (is.null(label) || length(label) == 0L || !nzchar(as.character(label)[1])) {
    return("default")
  }
  label <- as.character(label)[1]
  fixed <- make.names(label)
  if (!identical(fixed, label)) {
    stop(
      "label '", label, "' is not a valid make.names() identifier; ",
      "use '", fixed, "' instead"
    )
  }
  label
}

.get.report.plots.dir <- function() {
  plots_dir <- if (!is.null(config$plotpath) && nzchar(config$plotpath)) {
    config$plotpath
  } else {
    "~/mnt_out/plots/"
  }
  path.expand(plots_dir)
}

.get.report.cache.dir <- function() {
  cache_dir <- if (!is.null(config$cachedir) && nzchar(config$cachedir)) {
    config$cachedir
  } else {
    "~/mnt_out/rcache"
  }
  path.expand(cache_dir)
}

# Ensure global `reports` list exists; return it
.ensure.reports <- function() {
  if (!exists("reports", envir = .GlobalEnv, inherits = FALSE)) {
    assign("reports", list(), envir = .GlobalEnv)
  }
  get("reports", envir = .GlobalEnv)
}

# Empty report section template (one figure/pdf entry)
.empty.report.element <- function() {
  list(
    figure = NULL,
    pdf = NULL,
    png = NULL,
    description = "",
    h = 1,
    w = 1,
    sca = 1.5
  )
}

# Empty labeled report object
.empty.report <- function(label = "default") {
  label <- .normalize.report.label(label)
  list(
    cache = .get.report.cache.name(label),
    label = label,
    sections = list()
  )
}

# Locate a PDF by name under ~/mnt_out/plots (exact path, exact name, or best fuzzy match)
.find.report.pdf <- function(pdf_name, plots_dir = .get.report.plots.dir()) {
  if (is.null(pdf_name) || !nzchar(as.character(pdf_name)[1])) {
    stop("pdf name is empty")
  }
  pdf_name <- as.character(pdf_name)[1]

  if (file.exists(pdf_name)) {
    return(normalizePath(pdf_name))
  }

  direct <- file.path(plots_dir, pdf_name)
  if (file.exists(direct)) {
    return(normalizePath(direct))
  }

  if (!grepl("\\.pdf$", pdf_name, ignore.case = TRUE)) {
    with_ext <- file.path(plots_dir, paste0(pdf_name, ".pdf"))
    if (file.exists(with_ext)) {
      return(normalizePath(with_ext))
    }
  }

  stem <- sub("\\.pdf$", "", basename(pdf_name), ignore.case = TRUE)
  # Escape regex metacharacters so stems like "fig_1.0" match literally
  stem_re <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", stem)
  hits <- list.files(plots_dir, pattern = stem_re, full.names = TRUE, ignore.case = TRUE)
  hits <- hits[grepl("\\.pdf$", hits, ignore.case = TRUE)]
  if (length(hits) == 0) {
    stop("Could not find PDF '", pdf_name, "' under ", plots_dir)
  }
  hits[which.max(file.info(hits)$mtime)]
}

# Write a native plot object to an exact PDF path using stdplot sizing
.write.report.figure.pdf <- function(figure, dest_path, w = 1, h = 1, sca = 1.5) {
  w <- w %||% 1
  h <- h %||% 1
  sca <- sca %||% 1.5
  grDevices::pdf(dest_path, width = pw * sca * w, height = pw * sca * h)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  print(figure)
  invisible(dest_path)
}

# Write a native plot object to an exact PNG path using stdpng sizing
.write.report.figure.png <- function(figure, dest_path, w = 1, h = 1, sca = 1.5, res = 600) {
  w <- w %||% 1
  h <- h %||% 1
  sca <- sca %||% 1.5
  res <- res %||% 600
  base <- get("w", envir = .GlobalEnv)
  grDevices::png(
    dest_path,
    width = base * sca * w,
    height = base * sca * h,
    res = res
  )
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  print(figure)
  invisible(dest_path)
}

# Locate a PNG companion for a PDF name under the plots dir (optional)
.find.report.png <- function(pdf_name, plots_dir = .get.report.plots.dir()) {
  stem <- sub("\\.pdf$", "", as.character(pdf_name)[1], ignore.case = TRUE)
  candidates <- c(
    stem,
    paste0(stem, ".png"),
    file.path(plots_dir, paste0(basename(stem), ".png"))
  )
  for (cand in candidates) {
    if (file.exists(cand) && grepl("\\.png$", cand, ignore.case = TRUE)) {
      return(normalizePath(cand))
    }
  }
  tryCatch(
    {
      p <- .find.report.pdf(paste0(basename(stem), ".png"), plots_dir = plots_dir)
      if (grepl("\\.png$", p, ignore.case = TRUE)) return(p)
      NULL
    },
    error = function(e) NULL
  )
}

# Persist one labeled report into report_<label>_<analysis.version> under ~/mnt_out/rcache
.cache.report <- function(report_obj) {
  if (is.null(report_obj$label) || is.null(report_obj$cache)) {
    stop("report object must contain $label and $cache")
  }
  label <- .normalize.report.label(report_obj$label)
  cache_name <- report_obj$cache %||% .get.report.cache.name(label)
  report_obj$label <- label
  report_obj$cache <- cache_name

  reports <- .ensure.reports()
  reports[[label]] <- report_obj
  assign("reports", reports, envir = .GlobalEnv)

  cache_dir <- .get.report.cache.dir()
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # Stash under a stable temp name for simpleCache instruction lookup
  assign(".report.to.cache", report_obj, envir = .GlobalEnv)
  on.exit(rm(".report.to.cache", envir = .GlobalEnv), add = TRUE)

  simpleCache::simpleCache(
    cache_name,
    instruction = {
      get(".report.to.cache", envir = .GlobalEnv)
    },
    assignToVar = "report",
    recreate = TRUE,
    cacheDir = cache_dir
  )

  # Keep workspace `reports` as the multi-report container
  assign("reports", reports, envir = .GlobalEnv)
  invisible(report_obj)
}

#' Initiate an empty labeled report inside `reports`, and cache it
#'
#' @param label Report name. NULL becomes "default". Must be a valid make.names() id.
initiate_report <- function(label = NULL) {
  label <- .normalize.report.label(label)
  report_obj <- .empty.report(label)
  .cache.report(report_obj)
  message(
    "Initiated empty report '", label, "' and cached as '",
    report_obj$cache, ".RData' in ",
    .get.report.cache.dir()
  )
  invisible(get("reports", envir = .GlobalEnv))
}

#' Load a labeled report cache into `reports[[label]]`
#'
#' @param label Report name. NULL becomes "default".
#' @param initiate If TRUE (default) and no cache exists, run
#'   `initiate_report(label)`. If FALSE, error when the cache is missing.
load_report <- function(label = NULL, initiate = TRUE) {
  label <- .normalize.report.label(label)
  cache_name <- .get.report.cache.name(label)
  cache_dir <- .get.report.cache.dir()
  cache_file <- file.path(cache_dir, paste0(cache_name, ".RData"))
  if (!file.exists(cache_file)) {
    if (isTRUE(initiate)) {
      message(
        "No cached report for '", label, "'; initiating a new one at ",
        cache_file
      )
      return(initiate_report(label = label))
    }
    stop(
      "No cached report found at ", cache_file,
      ". Run initiate_report(label='", label, "') first."
    )
  }

  simpleCache::simpleCache(
    cache_name,
    assignToVar = "report",
    reload = TRUE,
    cacheDir = cache_dir
  )

  report_obj <- NULL
  if (exists("report", inherits = TRUE)) {
    report_obj <- get("report", inherits = TRUE)
  } else {
    e <- new.env(parent = emptyenv())
    load(cache_file, envir = e)
    if (exists("report", envir = e, inherits = FALSE)) {
      report_obj <- get("report", envir = e)
    } else {
      objs <- ls(envir = e)
      if (length(objs) == 1L) {
        report_obj <- get(objs[[1]], envir = e)
      } else {
        stop("Failed to load report from cache '", cache_name, "'")
      }
    }
  }

  # Normalize older flat caches (plain section list) into the new structure
  if (is.null(report_obj$sections)) {
    if (is.list(report_obj) && (is.null(names(report_obj)) || all(names(report_obj) %in% c("", NA)))) {
      sections <- report_obj
    } else if (is.list(report_obj) && (!is.null(report_obj$figure) || !is.null(report_obj$pdf))) {
      sections <- list(report_obj)
    } else {
      sections <- list()
      # If it looks like a named list of section-like elements, keep values
      if (is.list(report_obj) && length(report_obj) > 0 &&
          all(vapply(report_obj, function(x) {
            is.list(x) && (!is.null(x$figure) || !is.null(x$pdf) || !is.null(x$description))
          }, logical(1)))) {
        sections <- unname(report_obj)
      }
    }
    report_obj <- list(
      cache = cache_name,
      label = label,
      sections = sections
    )
  } else {
    report_obj$cache <- cache_name
    report_obj$label <- label
    if (is.null(report_obj$sections)) report_obj$sections <- list()
  }

  reports <- .ensure.reports()
  reports[[label]] <- report_obj
  assign("reports", reports, envir = .GlobalEnv)
  # Drop the temporary single-report binding if present
  if (exists("report", envir = .GlobalEnv, inherits = FALSE)) {
    rm("report", envir = .GlobalEnv)
  }

  message(
    "Loaded report '", label, "' (", cache_name, ") with ",
    length(report_obj$sections), " section(s)"
  )
  invisible(reports)
}

#' Add a PDF name (string) or native plot object to a labeled report
#'
#' @param a Character PDF name/path, or a printable R plot object
#' @param n Index within $sections; NULL appends
#' @param label Report name. NULL becomes "default".
#' @param description Optional Rmd text shown before the figure
#' @param w,h,sca Optional plot size parameters (defaults match stdplot)
#' @param cache If TRUE, write the updated report to disk cache after adding.
#'   Default FALSE: only update the in-memory `reports[[label]]`.
#' @param insert If TRUE (default) and `n` already exists, insert the new
#'   section at `n` and shift the old `n`, `n+1`, ... one slot downstream.
#'   If FALSE, overwrite section `n`.
add_to_report <- function(a, n = NULL, label = NULL, description = NULL, w = NULL, h = NULL, sca = NULL, cache = FALSE, insert = TRUE) {
  label <- .normalize.report.label(label)
  reports <- .ensure.reports()
  if (is.null(reports[[label]])) {
    cache_name <- .get.report.cache.name(label)
    cache_file <- file.path(.get.report.cache.dir(), paste0(cache_name, ".RData"))
    if (file.exists(cache_file)) {
      load_report(label = label, initiate = FALSE)
      reports <- get("reports", envir = .GlobalEnv)
    } else {
      stop(
        "No report '", label, "' in workspace and no cache at ", cache_file, ". ",
        "Run initiate_report(label='", label, "') first."
      )
    }
  }
  report_obj <- reports[[label]]
  if (is.null(report_obj$sections)) report_obj$sections <- list()
  sections <- report_obj$sections

  if (is.null(n)) {
    n <- length(sections) + 1L
  } else {
    n <- as.integer(n)
    if (is.na(n) || n < 1L) {
      stop("n must be a positive integer or NULL")
    }
    if (n > length(sections) + 1L) {
      stop(
        "n=", n, " is beyond the next append slot (length=",
        length(sections), "). Append with n=NULL or n=", length(sections) + 1L, "."
      )
    }
  }

  if (n <= length(sections) && isTRUE(insert)) {
    sections <- append(sections, list(.empty.report.element()), after = n - 1L)
  } else if (n > length(sections)) {
    sections[[n]] <- .empty.report.element()
  } else if (is.null(sections[[n]]) || !is.list(sections[[n]])) {
    sections[[n]] <- .empty.report.element()
  } else {
    defaults <- .empty.report.element()
    for (nm in names(defaults)) {
      if (is.null(sections[[n]][[nm]]) && nm %in% c("description", "h", "w", "sca")) {
        sections[[n]][[nm]] <- defaults[[nm]]
      }
    }
  }

  if (is.character(a)) {
    if (length(a) != 1L || !nzchar(a)) {
      stop("When a is a string, it must be a single non-empty PDF name or path")
    }
    sections[[n]]$pdf <- a
  } else {
    sections[[n]]$figure <- a
  }

  if (!is.null(description)) sections[[n]]$description <- description
  if (!is.null(w)) sections[[n]]$w <- w
  if (!is.null(h)) sections[[n]]$h <- h
  if (!is.null(sca)) sections[[n]]$sca <- sca

  has_pdf <- !is.null(sections[[n]]$pdf) && nzchar(as.character(sections[[n]]$pdf)[1])
  has_figure <- !is.null(sections[[n]]$figure)
  if (!has_pdf && !has_figure) {
    stop("Each report section needs at least one of $pdf or $figure")
  }

  report_obj$sections <- sections
  report_obj$label <- label
  report_obj$cache <- .get.report.cache.name(label)

  if (isTRUE(cache)) {
    .cache.report(report_obj)
  } else {
    reports[[label]] <- report_obj
    assign("reports", reports, envir = .GlobalEnv)
  }
  invisible(get("reports", envir = .GlobalEnv))
}

#' Write the current in-memory report to its disk cache
#'
#' @param label Report name. NULL becomes "default".
cache_report <- function(label = NULL) {
  label <- .normalize.report.label(label)
  reports <- .ensure.reports()
  if (is.null(reports[[label]])) {
    stop(
      "No report '", label, "' in workspace. ",
      "Run initiate_report(label='", label, "') or load_report(label='", label, "') first."
    )
  }
  .cache.report(reports[[label]])
  message("Cached report '", label, "' as '", reports[[label]]$cache, ".RData'")
  invisible(get("reports", envir = .GlobalEnv))
}

# Resolve a section $pdf reference to an existing file path.
# Prefers paths already inside report_dir (portable relative figures/...),
# then any prior report dirs, then absolute/cwd paths, then plots lookup.
.resolve.report.pdf.source <- function(pdf_ref, report_dirs = NULL) {
  if (is.null(pdf_ref) || !nzchar(as.character(pdf_ref)[1])) {
    stop("pdf name is empty")
  }
  pdf_ref <- as.character(pdf_ref)[1]
  pdf_ref_norm <- sub("^\\./", "", pdf_ref)

  for (rd in unique(as.character(report_dirs %||% character()))) {
    if (!nzchar(rd)) next
    for (ref in unique(c(pdf_ref_norm, pdf_ref))) {
      in_report <- file.path(rd, ref)
      if (file.exists(in_report)) {
        return(normalizePath(in_report))
      }
    }
  }

  if (file.exists(pdf_ref)) {
    return(normalizePath(pdf_ref))
  }
  if (file.exists(pdf_ref_norm)) {
    return(normalizePath(pdf_ref_norm))
  }

  .find.report.pdf(pdf_ref_norm)
}

# Resolve a PNG companion for a pdf_ref, preferring files inside report dirs.
.resolve.report.png.source <- function(pdf_ref, report_dirs = NULL) {
  pdf_ref <- as.character(pdf_ref)[1]
  pdf_ref_norm <- sub("^\\./", "", pdf_ref)
  stem <- sub("\\.pdf$", "", basename(pdf_ref_norm), ignore.case = TRUE)
  for (rd in unique(as.character(report_dirs %||% character()))) {
    if (!nzchar(rd)) next
    candidates <- c(
      file.path(rd, pdf_ref_norm),
      file.path(rd, pdf_ref),
      file.path(rd, "figures", paste0(stem, ".png")),
      file.path(rd, sub("\\.pdf$", ".png", pdf_ref_norm, ignore.case = TRUE))
    )
    for (cand in candidates) {
      if (file.exists(cand) && grepl("\\.png$", cand, ignore.case = TRUE)) {
        return(normalizePath(cand))
      }
    }
  }
  .find.report.png(pdf_ref_norm)
}

#' Collect report figures into ~/mnt_out/reports/YYMMDD_report_LABEL_ANALYSISVERSION
#'
#' For each section: copy $pdf from ~/mnt_out/plots as report_fig_N.pdf when
#' present (unless rewrite=TRUE and $figure exists, in which case regenerate);
#' otherwise render $figure with stdplot/stdpng into report_fig_N.pdf and
#' report_fig_N.png. PNGs are high-resolution (res=600) for HTML knitting.
#'
#' After writing, each section's `$pdf` / `$png` are set according to
#' `path.type`: relative (`./figures/report_fig_N.pdf`) or absolute (full path).
#' The same `path.type` is stored on the report and used by
#' `create_report_markdown()` for Rmd `root.dir` and `include_graphics` paths.
#'
#' @param rewrite If TRUE and both $pdf and $figure exist, regenerate from $figure
#' @param label Report name. NULL becomes "default".
#' @param res PNG resolution passed to stdpng / .write.report.figure.png
#' @param path.type `"relative"` (default) or `"absolute"`: how figure paths are
#'   stored on sections and later written into the Rmd.
collect.report <- function(rewrite = FALSE, label = NULL, res = 600,
                           path.type = "relative") {
  path.type <- match.arg(path.type, c("relative", "absolute"))
  label <- .normalize.report.label(label)
  reports <- .ensure.reports()
  if (is.null(reports[[label]])) {
    stop(
      "No report '", label, "' in workspace. ",
      "Run initiate_report(label='", label, "') or load_report(label='", label, "') first."
    )
  }
  report_obj <- reports[[label]]
  sections <- report_obj$sections %||% list()
  if (length(sections) == 0) {
    stop("Report '", label, "' has no sections; add elements with add_to_report() first.")
  }

  analysis_version <- .get.report.analysis.version()
  report_root <- path.expand("~/mnt_out/reports")
  folder_name <- paste0(format(Sys.Date(), "%y%m%d"), "_report_", label, "_", analysis_version)
  report_dir <- file.path(report_root, folder_name)
  figures_dir <- file.path(report_dir, "figures")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  # Prefer the new folder, then a previous collected folder (for relative $pdf)
  search_dirs <- unique(c(report_dir, report_obj$report_dir %||% character()))

  for (n in seq_along(sections)) {
    el <- sections[[n]]
    if (is.null(el) || !is.list(el)) {
      stop("reports[['", label, "']]$sections[[", n, "]] is not a valid report section")
    }

    has_pdf <- !is.null(el$pdf) && nzchar(as.character(el$pdf)[1])
    has_figure <- !is.null(el$figure)
    if (!has_pdf && !has_figure) {
      stop("reports[['", label, "']]$sections[[", n, "]] has neither $pdf nor $figure")
    }

    rel_pdf <- file.path("figures", paste0("report_fig_", n, ".pdf"))
    rel_png <- file.path("figures", paste0("report_fig_", n, ".png"))
    # Portable form used in Rmds when path.type = "relative"
    rel_pdf_dot <- paste0("./", gsub("\\\\", "/", rel_pdf))
    rel_png_dot <- paste0("./", gsub("\\\\", "/", rel_png))
    dest_pdf <- file.path(report_dir, rel_pdf)
    dest_png <- file.path(report_dir, rel_png)
    w <- el$w %||% 1
    h <- el$h %||% 1
    sca <- el$sca %||% 1.5
    title <- paste0("report_fig_", label, "_", n)

    write_from_figure <- function() {
      stdplot(el$figure, title, w = w, h = h, sca = sca)
      stdpng(el$figure, title, w = w, h = h, sca = sca, res = res)
      .write.report.figure.pdf(el$figure, dest_pdf, w = w, h = h, sca = sca)
      .write.report.figure.png(el$figure, dest_png, w = w, h = h, sca = sca, res = res)
    }

    copy_pdf_and_maybe_png <- function() {
      src <- .resolve.report.pdf.source(el$pdf, report_dirs = search_dirs)
      # Skip self-copy when $pdf already points at this report's figure
      if (!identical(normalizePath(src, mustWork = FALSE), normalizePath(dest_pdf, mustWork = FALSE)) ||
          !file.exists(dest_pdf)) {
        file.copy(src, dest_pdf, overwrite = TRUE)
      }
      png_src <- .resolve.report.png.source(el$pdf, report_dirs = search_dirs)
      if (!is.null(png_src) && file.exists(png_src)) {
        if (!identical(normalizePath(png_src, mustWork = FALSE), normalizePath(dest_png, mustWork = FALSE)) ||
            !file.exists(dest_png)) {
          file.copy(png_src, dest_png, overwrite = TRUE)
        }
      } else if (has_figure) {
        stdpng(el$figure, title, w = w, h = h, sca = sca, res = res)
        .write.report.figure.png(el$figure, dest_png, w = w, h = h, sca = sca, res = res)
      } else {
        warning(
          "No PNG available for section ", n,
          " (pdf='", el$pdf, "'). HTML markdown will need report_fig_", n, ".png."
        )
      }
    }

    if (has_pdf && has_figure) {
      if (isTRUE(rewrite)) {
        write_from_figure()
      } else {
        copy_pdf_and_maybe_png()
      }
    } else if (has_pdf) {
      copy_pdf_and_maybe_png()
    } else {
      write_from_figure()
    }

    if (identical(path.type, "relative")) {
      sections[[n]]$pdf <- rel_pdf_dot
      if (file.exists(dest_png)) {
        sections[[n]]$png <- rel_png_dot
      } else {
        sections[[n]]$png <- NULL
      }
      message("Wrote ", rel_pdf_dot)
      if (file.exists(dest_png)) message("Wrote ", rel_png_dot)
    } else {
      sections[[n]]$pdf <- normalizePath(dest_pdf)
      if (file.exists(dest_png)) {
        sections[[n]]$png <- normalizePath(dest_png)
      } else {
        sections[[n]]$png <- NULL
      }
      message("Wrote ", sections[[n]]$pdf)
      if (!is.null(sections[[n]]$png)) message("Wrote ", sections[[n]]$png)
    }
  }

  report_obj$sections <- sections
  report_obj$report_dir <- if (identical(path.type, "absolute")) {
    normalizePath(report_dir)
  } else {
    report_dir
  }
  report_obj$path.type <- path.type
  reports[[label]] <- report_obj
  assign("reports", reports, envir = .GlobalEnv)

  message("Report figures collected in ", report_dir, " (path.type=", path.type, ")")
  invisible(report_dir)
}

# List YYMMDD_report_LABEL_* folders; date="latest" or a YYMMDD / YYYYMMDD prefix
.find.report.folder.by.date <- function(label, date = "latest") {
  label <- .normalize.report.label(label)
  report_root <- path.expand("~/mnt_out/reports")
  if (!dir.exists(report_root)) {
    stop("Report root does not exist: ", report_root)
  }
  dirs <- list.dirs(report_root, recursive = FALSE, full.names = TRUE)
  # YYMMDD_report_<label>_...
  prefix_re <- paste0("^[0-9]{6}_report_", gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", label), "(_|$)")
  hits <- dirs[grepl(prefix_re, basename(dirs))]
  if (length(hits) == 0) {
    stop(
      "No report folder for label '", label, "' under ", report_root,
      ". Run save_report(label='", label, "') first."
    )
  }

  dates <- substr(basename(hits), 1, 6)
  if (is.null(date) || identical(as.character(date)[1], "latest")) {
    max_date <- max(dates)
    hits <- hits[dates == max_date]
    return(hits[which.max(file.info(hits)$mtime)])
  }

  date <- as.character(date)[1]
  date6 <- gsub("[^0-9]", "", date)
  if (nchar(date6) == 8) date6 <- substr(date6, 3, 8)
  if (nchar(date6) != 6) {
    stop("date must be 'latest' or a YYMMDD / YYYYMMDD string; got '", date, "'")
  }
  hits <- hits[dates == date6]
  if (length(hits) == 0) {
    stop("No report folder for label '", label, "' with date ", date6, " under ", report_root)
  }
  hits[which.max(file.info(hits)$mtime)]
}

.load.rdata.object <- function(path, preferred = NULL) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  objs <- ls(envir = e)
  if (length(objs) == 0) {
    stop("No objects in ", path)
  }
  if (!is.null(preferred) && preferred %in% objs) {
    return(get(preferred, envir = e))
  }
  get(objs[[1]], envir = e)
}

.list.report.fig.indices <- function(figures_dir) {
  if (!dir.exists(figures_dir)) {
    return(integer(0))
  }
  files <- list.files(
    figures_dir,
    pattern = "^report_fig_[0-9]+\\.(RData|rdata|pdf|png)$",
    ignore.case = TRUE
  )
  if (length(files) == 0) {
    return(integer(0))
  }
  idx <- as.integer(sub("^report_fig_([0-9]+)\\..*$", "\\1", files, ignore.case = TRUE))
  sort(unique(idx[!is.na(idx)]))
}

#' Save a labeled report into its collected folder (PDFs, PNGs, per-section RData)
#'
#' Runs `collect.report(label)` first, then writes each section as
#' `figures/report_fig_N.RData` (same stem as the PDF/PNG) and writes
#' `report_structure.RData` in the report folder with `$sections = list()`.
#'
#' @param label Report name. NULL becomes "default".
#' @param rewrite,res,path.type Passed to `collect.report()`.
save_report <- function(label = NULL, rewrite = FALSE, res = 600,
                        path.type = "relative") {
  label <- .normalize.report.label(label)
  report_dir <- collect.report(
    rewrite = rewrite, label = label, res = res, path.type = path.type
  )
  reports <- .ensure.reports()
  report_obj <- reports[[label]]
  figures_dir <- file.path(report_dir, "figures")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  sections <- report_obj$sections %||% list()
  for (n in seq_along(sections)) {
    section <- sections[[n]]
    rdata_path <- file.path(figures_dir, paste0("report_fig_", n, ".RData"))
    save(section, file = rdata_path)
    message("Wrote ", rdata_path)
  }

  report_structure <- report_obj
  report_structure$sections <- list()
  report_structure$report_dir <- report_dir
  structure_path <- file.path(report_dir, "report_structure.RData")
  save(report_structure, file = structure_path)
  message("Wrote ", structure_path)

  report_obj$report_dir <- report_dir
  reports[[label]] <- report_obj
  assign("reports", reports, envir = .GlobalEnv)
  invisible(report_dir)
}

#' Load a labeled report from its collected folder into `reports[[label]]`
#'
#' Reads `report_structure.RData`, then each `figures/report_fig_N.RData`.
#' `$pdf` / `$png` are set to portable paths relative to the report folder
#' (`figures/report_fig_N.pdf` / `.png`) when those files exist.
#' Missing section RData becomes a blank section (`$h/$w/$sca` defaults) and
#' still picks up PDF/PNG files if they exist.
#'
#' If no matching report folder exists and `initiate = TRUE` (default), runs
#' `initiate_report(label)` instead of erroring.
#'
#' @param label Report name. NULL becomes "default".
#' @param date `"latest"` (default) or a YYMMDD / YYYYMMDD folder date.
#' @param initiate If TRUE (default) and no folder exists, run
#'   `initiate_report(label)`. If FALSE, error when the folder is missing.
load_report_folder <- function(label = NULL, date = "latest", initiate = TRUE) {
  label <- .normalize.report.label(label)
  report_dir <- tryCatch(
    .find.report.folder.by.date(label, date = date),
    error = function(e) NULL
  )
  if (is.null(report_dir) || !dir.exists(report_dir)) {
    if (isTRUE(initiate)) {
      message(
        "No report folder for '", label, "' (date='", date,
        "'); initiating a new empty report"
      )
      return(initiate_report(label = label))
    }
    stop(
      "No report folder for label '", label, "' (date='", date, "'). ",
      "Run save_report(label='", label, "') or initiate_report(label='", label, "') first."
    )
  }

  figures_dir <- file.path(report_dir, "figures")
  structure_path <- file.path(report_dir, "report_structure.RData")

  if (file.exists(structure_path)) {
    report_obj <- .load.rdata.object(structure_path, preferred = "report_structure")
    if (!is.list(report_obj)) {
      stop("report_structure.RData did not contain a list")
    }
  } else {
    warning("No report_structure.RData in ", report_dir, "; starting from an empty report")
    report_obj <- .empty.report(label)
  }

  report_obj$label <- label
  report_obj$report_dir <- normalizePath(report_dir)
  report_obj$sections <- list()

  indices <- .list.report.fig.indices(figures_dir)
  if (length(indices) > 0) {
    for (n in seq_len(max(indices))) {
      rdata_path <- file.path(figures_dir, paste0("report_fig_", n, ".RData"))
      rel_pdf <- file.path("figures", paste0("report_fig_", n, ".pdf"))
      rel_png <- file.path("figures", paste0("report_fig_", n, ".png"))
      pdf_path <- file.path(report_dir, rel_pdf)
      png_path <- file.path(report_dir, rel_png)

      if (file.exists(rdata_path)) {
        el <- .load.rdata.object(rdata_path, preferred = "section")
        if (!is.list(el)) {
          warning("report_fig_", n, ".RData is not a list; using a blank section")
          el <- .empty.report.element()
          el$pdf <- ""
          el$png <- ""
        }
      } else {
        el <- .empty.report.element()
        el$pdf <- ""
        el$png <- ""
        el$h <- el$h %||% 1
        el$w <- el$w %||% 1
        el$sca <- el$sca %||% 1.5
      }

      if (file.exists(pdf_path)) {
        el$pdf <- rel_pdf
      } else if (is.null(el$pdf)) {
        el$pdf <- ""
      }
      if (file.exists(png_path)) {
        el$png <- rel_png
      } else if (is.null(el$png)) {
        el$png <- ""
      }

      report_obj$sections[[n]] <- el
    }
  }

  reports <- .ensure.reports()
  reports[[label]] <- report_obj
  assign("reports", reports, envir = .GlobalEnv)
  message(
    "Loaded report folder '", basename(report_dir), "' into reports[['", label,
    "']] with ", length(report_obj$sections), " section(s)"
  )
  invisible(get("reports", envir = .GlobalEnv))
}


.is.absolute.path <- function(path) {
  if (is.null(path) || length(path) == 0L) return(FALSE)
  path <- as.character(path)[1]
  if (!nzchar(path) || is.na(path)) return(FALSE)
  grepl("^(~|/|[A-Za-z]:)", path)
}

# Convert one path string relative to base_dir (report folder) or report_root.
.change.one.report.path <- function(path, type, base_dir = NULL, report_root = NULL) {
  if (is.null(path) || length(path) == 0L) return(path)
  path <- as.character(path)[1]
  if (is.na(path) || !nzchar(path)) return(path)

  type <- match.arg(type, c("relative", "absolute"))

  if (identical(type, "relative")) {
    if (!.is.absolute.path(path)) {
      return(gsub("\\\\", "/", path))
    }
    abs_path <- normalizePath(path.expand(path), mustWork = FALSE)
    # Prefer paths under the report folder
    if (!is.null(base_dir) && nzchar(as.character(base_dir)[1])) {
      base_exp <- path.expand(as.character(base_dir)[1])
      if (dir.exists(base_exp) || file.exists(base_exp)) {
        base_n <- normalizePath(base_exp)
        if (identical(abs_path, base_n)) return(".")
        prefix <- paste0(base_n, .Platform$file.sep)
        if (startsWith(abs_path, prefix)) {
          return(gsub("\\\\", "/", substring(abs_path, nchar(prefix) + 1L)))
        }
      }
    }
    # Else prefer paths under the reports root (e.g. report_dir itself)
    if (!is.null(report_root) && nzchar(as.character(report_root)[1])) {
      root_exp <- path.expand(as.character(report_root)[1])
      if (dir.exists(root_exp)) {
        root_n <- normalizePath(root_exp)
        if (identical(abs_path, root_n)) return(".")
        prefix <- paste0(root_n, .Platform$file.sep)
        if (startsWith(abs_path, prefix)) {
          return(gsub("\\\\", "/", substring(abs_path, nchar(prefix) + 1L)))
        }
      }
    }
    return(gsub("\\\\", "/", abs_path))
  }

  # type == "absolute"
  if (.is.absolute.path(path)) {
    exp <- path.expand(path)
    if (file.exists(exp) || dir.exists(exp)) {
      return(normalizePath(exp))
    }
    return(gsub("\\\\", "/", normalizePath(exp, mustWork = FALSE)))
  }

  # Relative: resolve against report folder first, then reports root, then plots
  candidates <- character()
  if (!is.null(base_dir) && nzchar(as.character(base_dir)[1])) {
    candidates <- c(candidates, file.path(path.expand(as.character(base_dir)[1]), path))
  }
  if (!is.null(report_root) && nzchar(as.character(report_root)[1])) {
    candidates <- c(candidates, file.path(path.expand(as.character(report_root)[1]), path))
  }
  candidates <- c(candidates, path)
  for (cand in unique(candidates)) {
    if (file.exists(cand) || dir.exists(cand)) {
      return(normalizePath(cand))
    }
  }

  if (grepl("\\.(pdf|png)$", path, ignore.case = TRUE)) {
    found <- tryCatch(
      {
        if (grepl("\\.png$", path, ignore.case = TRUE)) {
          .find.report.png(path)
        } else {
          .find.report.pdf(path)
        }
      },
      error = function(e) NULL
    )
    if (!is.null(found) && file.exists(found)) return(found)
  }

  if (length(candidates) > 0) {
    return(gsub("\\\\", "/", normalizePath(candidates[[1]], mustWork = FALSE)))
  }
  path
}

.change.report.object.paths <- function(report_obj, type, report_root = path.expand("~/mnt_out/reports")) {
  if (is.null(report_obj) || !is.list(report_obj)) return(report_obj)
  type <- match.arg(type, c("relative", "absolute"))

  # Resolve / normalize report_dir first so other fields can use it as base
  base_dir <- report_obj$report_dir
  if (!is.null(base_dir) && nzchar(as.character(base_dir)[1])) {
    base_raw <- as.character(base_dir)[1]
    if (!.is.absolute.path(base_raw) && !is.null(report_root)) {
      cand <- file.path(report_root, base_raw)
      if (dir.exists(cand)) base_dir <- cand
    } else if (dir.exists(path.expand(base_raw))) {
      base_dir <- path.expand(base_raw)
    }
  } else {
    base_dir <- NULL
  }

  # report_dir: relative to reports root; absolute = full path
  if (!is.null(report_obj$report_dir) && nzchar(as.character(report_obj$report_dir)[1])) {
    report_obj$report_dir <- .change.one.report.path(
      report_obj$report_dir,
      type = type,
      base_dir = NULL,
      report_root = report_root
    )
    if (identical(type, "absolute") && dir.exists(report_obj$report_dir)) {
      base_dir <- report_obj$report_dir
    } else if (identical(type, "relative") && !is.null(base_dir) && dir.exists(base_dir)) {
      # keep absolute base_dir for resolving sibling fields even if stored value is relative
      base_dir <- normalizePath(base_dir)
    }
  }

  report_path_fields <- c("rmd", "pdf", "html", "rmd_copy")
  for (field in report_path_fields) {
    if (!is.null(report_obj[[field]])) {
      report_obj[[field]] <- .change.one.report.path(
        report_obj[[field]],
        type = type,
        base_dir = base_dir,
        report_root = report_root
      )
    }
  }

  sections <- report_obj$sections %||% list()
  for (n in seq_along(sections)) {
    el <- sections[[n]]
    if (is.null(el) || !is.list(el)) next
    for (field in c("pdf", "png")) {
      if (!is.null(el[[field]])) {
        el[[field]] <- .change.one.report.path(
          el[[field]],
          type = type,
          base_dir = base_dir,
          report_root = report_root
        )
      }
    }
    sections[[n]] <- el
  }
  report_obj$sections <- sections
  report_obj
}

#' Convert path fields in `reports` to relative or absolute form
#'
#' Walks every labeled report and rewrites path-like fields:
#' - report-level: `$report_dir`, `$rmd`, `$pdf`, `$html`, `$rmd_copy`
#' - section-level: `$pdf`, `$png`
#'
#' Relative paths are anchored to the report folder when possible
#' (`figures/report_fig_N.pdf`, Rmd basename, etc.). `$report_dir` itself is
#' stored relative to `~/mnt_out/reports` when `type = "relative"`.
#'
#' @param type `"relative"` (default) or `"absolute"`.
#' @param label Optional single report name; NULL converts all reports.
#' @return The updated global `reports` list (invisibly).
change_report_paths <- function(type = "relative", label = NULL) {
  type <- match.arg(type, c("relative", "absolute"))
  reports <- .ensure.reports()
  report_root <- path.expand("~/mnt_out/reports")

  if (is.null(label)) {
    labels <- names(reports)
    labels <- labels[!is.na(labels) & nzchar(labels)]
  } else {
    labels <- .normalize.report.label(label)
    if (is.null(reports[[labels]])) {
      stop("No report '", labels, "' in workspace reports")
    }
  }

  if (length(labels) == 0) {
    warning("No labeled reports in workspace")
    return(invisible(reports))
  }

  for (lab in labels) {
    reports[[lab]] <- .change.report.object.paths(
      reports[[lab]],
      type = type,
      report_root = report_root
    )
    message("Updated paths in reports[['", lab, "']] to ", type)
  }

  assign("reports", reports, envir = .GlobalEnv)
  invisible(get("reports", envir = .GlobalEnv))
}

.find.report.dir <- function(label, report_dir = NULL) {
  label <- .normalize.report.label(label)
  if (!is.null(report_dir) && nzchar(as.character(report_dir)[1])) {
    report_dir <- path.expand(as.character(report_dir)[1])
    # Allow report_dir stored relative to ~/mnt_out/reports
    if (!dir.exists(report_dir)) {
      root <- path.expand("~/mnt_out/reports")
      cand <- file.path(root, report_dir)
      if (dir.exists(cand)) {
        report_dir <- cand
      } else {
        stop("report_dir does not exist: ", report_dir)
      }
    }
    return(normalizePath(report_dir))
  }

  reports <- .ensure.reports()
  stored <- reports[[label]]$report_dir
  if (!is.null(stored) && nzchar(as.character(stored)[1])) {
    stored <- as.character(stored)[1]
    if (dir.exists(path.expand(stored))) {
      return(normalizePath(path.expand(stored)))
    }
    root <- path.expand("~/mnt_out/reports")
    cand <- file.path(root, stored)
    if (dir.exists(cand)) {
      return(normalizePath(cand))
    }
  }

  analysis_version <- .get.report.analysis.version()
  report_root <- path.expand("~/mnt_out/reports")
  suffix <- paste0("_report_", label, "_", analysis_version)
  dirs <- list.dirs(report_root, recursive = FALSE, full.names = TRUE)
  hits <- dirs[endsWith(basename(dirs), suffix)]
  if (length(hits) == 0) {
    stop(
      "Could not find a collected report folder matching '*", suffix,
      "' under ", report_root, ". Run collect.report(label='", label,
      "') first, or pass report_dir= explicitly."
    )
  }
  hits[which.max(file.info(hits)$mtime)]
}

#' Write an Rmd from a labeled report and its collected figures folder
#'
#' For each section, writes the `$description` followed by a knitr chunk that
#' displays `./figures/report_fig_N.png` (or the absolute PNG path when
#' `path.type = "absolute"`). Figures are sized to the full HTML page width
#' (`out.width = "100%"`); `$w * $sca` and `$h * $sca` set `fig.width`/`fig.height`.
#'
#' The Rmd is written inside the collected report folder (next to `figures/`)
#' so an HTML knit can resolve those images. A copy is also saved under
#' `~/Rmd/report_markdowns/` (knit from the report folder, not the copy).
#'
#' @param label Report name. NULL becomes "default".
#' @param report_dir Folder returned by `collect.report()`. If NULL, uses
#'   `reports[[label]]$report_dir` or the newest matching folder under
#'   `~/mnt_out/reports`.
#' @param path.type `"relative"` (default) or `"absolute"`. Controls knitr
#'   `root.dir` (`.` vs the full report folder path) and figure paths in the Rmd
#'   (`./figures/...` vs absolute paths). Defaults to `reports[[label]]$path.type`
#'   when set by `collect.report()`.
create_report_markdown <- function(label = NULL, report_dir = NULL,
                                   path.type = NULL) {
  label <- .normalize.report.label(label)
  reports <- .ensure.reports()
  if (is.null(reports[[label]])) {
    stop(
      "No report '", label, "' in workspace. ",
      "Run initiate_report(label='", label, "') or load_report(label='", label, "') first."
    )
  }
  report_obj <- reports[[label]]
  sections <- report_obj$sections %||% list()
  if (length(sections) == 0) {
    stop("Report '", label, "' has no sections.")
  }

  if (is.null(path.type)) {
    path.type <- report_obj$path.type %||% "relative"
  }
  path.type <- match.arg(as.character(path.type)[1], c("relative", "absolute"))

  report_dir <- .find.report.dir(label, report_dir)
  report_dir_abs <- normalizePath(report_dir)
  report_dir_posix <- gsub("\\\\", "/", report_dir_abs)
  figures_dir <- file.path(report_dir_abs, "figures")
  if (!dir.exists(figures_dir)) {
    stop("No figures/ folder in ", report_dir_abs, ". Run collect.report(label='", label, "') first.")
  }

  analysis_version <- .get.report.analysis.version()
  folder_name <- basename(report_dir_abs)
  title_text <- paste("Report", label, analysis_version)
  title_text <- gsub('"', "'", title_text)

  if (identical(path.type, "relative")) {
    root_dir_line <- 'knitr::opts_knit$set(root.dir = ".")'
  } else {
    root_dir_line <- paste0('knitr::opts_knit$set(root.dir = "', report_dir_posix, '")')
  }

  lines <- c(
    "---",
    paste0('title: "', title_text, '"'),
    paste0("date: \"", format(Sys.Date(), "%Y-%m-%d"), "\""),
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float: true",
    "    self_contained: true",
    "  pdf_document:",
    "    toc: true",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(",
    "  echo = FALSE, warning = FALSE, message = FALSE,",
    "  out.width = if (knitr::is_latex_output()) '\\\\textwidth' else '100%'",
    ")",
    root_dir_line,
    "```",
    "",
    "```{css, echo=FALSE, eval=knitr::is_html_output()}",
    "img { max-width: 100%; height: auto; display: block; margin: 0.5em auto; }",
    ".figure { max-width: 100%; }",
    "```",
    ""
  )

  for (n in seq_along(sections)) {
    el <- sections[[n]]
    if (is.null(el) || !is.list(el)) {
      stop("reports[['", label, "']]$sections[[", n, "]] is not a valid report section")
    }

    desc <- el$description %||% ""
    if (is.character(desc) && length(desc) > 1L) {
      desc <- paste(desc, collapse = "\n")
    }
    desc <- as.character(desc)[1]
    if (is.na(desc)) desc <- ""

    w <- el$w %||% 1
    h <- el$h %||% 1
    sca <- el$sca %||% 1.5
    fig_w <- as.numeric(w) * as.numeric(sca)
    fig_h <- as.numeric(h) * as.numeric(sca)
    if (!is.finite(fig_w) || fig_w <= 0) fig_w <- 1.5
    if (!is.finite(fig_h) || fig_h <= 0) fig_h <- 1.5

    png_on_disk <- file.path(figures_dir, paste0("report_fig_", n, ".png"))
    if (!file.exists(png_on_disk)) {
      warning(
        "Missing PNG for section ", n, ": ", png_on_disk,
        ". Re-run collect.report(label='", label, "') to generate it via stdpng()."
      )
      next
    }

    if (identical(path.type, "relative")) {
      fig_path <- paste0("./figures/report_fig_", n, ".png")
    } else {
      fig_path <- gsub("\\\\", "/", normalizePath(png_on_disk))
    }

    if (nzchar(desc)) {
      lines <- c(lines, desc, "")
    }

    chunk_name <- paste0("fig-", n)
    # Fit figures to page/HTML width; preserve aspect via fig.width/fig.height.
    chunk_header <- paste0(
      "```{r ", chunk_name,
      ", echo=FALSE, warning=FALSE, message=FALSE",
      ", fig.width=", fig_w,
      ", fig.height=", fig_h,
      "}"
    )
    lines <- c(
      lines,
      chunk_header,
      paste0('knitr::include_graphics("', fig_path, '")'),
      "```",
      ""
    )
  }

  rmd_name <- paste0(folder_name, ".Rmd")
  rmd_path <- file.path(report_dir_abs, rmd_name)
  writeLines(lines, rmd_path, useBytes = TRUE)

  copy_dir <- path.expand("~/Rmd/report_markdowns")
  dir.create(copy_dir, recursive = TRUE, showWarnings = FALSE)
  copy_path <- file.path(copy_dir, rmd_name)
  file.copy(rmd_path, copy_path, overwrite = TRUE)

  report_obj$report_dir <- report_dir_abs
  report_obj$path.type <- path.type
  report_obj$rmd <- if (identical(path.type, "relative")) {
    rmd_name
  } else {
    rmd_path
  }
  reports[[label]] <- report_obj
  assign("reports", reports, envir = .GlobalEnv)

  message("Wrote report markdown: ", rmd_path, " (path.type=", path.type, ")")
  message("Copied report markdown: ", copy_path)
  invisible(list(rmd = rmd_path, copy = copy_path, report_dir = report_dir_abs, path.type = path.type))
}

# Multiply each section's $sca by `factor` on a copy of the report object
.scale.report.sections <- function(report_obj, factor) {
  factor <- as.numeric(factor)[1]
  if (!is.finite(factor) || factor <= 0) {
    stop("sca multiplier must be a positive number")
  }
  sections <- report_obj$sections %||% list()
  for (n in seq_along(sections)) {
    if (is.null(sections[[n]]) || !is.list(sections[[n]])) next
    sections[[n]]$sca <- (sections[[n]]$sca %||% 1.5) * factor
  }
  report_obj$sections <- sections
  report_obj
}

#' Collect and write markdown for every report in `reports`
#'
#' @param sca Scale control:
#'   - NULL: use each section's stored `$sca`
#'   - a single numeric scalar: multiply every section `$sca` by this value
#'   - a named list (e.g. `list(report1=3, report2=2)`): multiply `$sca` only
#'     for the named reports; others stay unchanged
#' @param rewrite Passed to `collect.report()` (regenerate PDFs from `$figure`)
#' @param res PNG resolution passed to `collect.report()`
#' @param path.type `"relative"` (default) or `"absolute"`: passed to
#'   `collect.report()` and `create_report_markdown()` to control Rmd
#'   `root.dir` and figure paths (`./figures/...` vs absolute).
make_all_reports <- function(sca = NULL, rewrite = FALSE, res = 600,
                             path.type = "relative") {
  path.type <- match.arg(path.type, c("relative", "absolute"))
  reports <- .ensure.reports()
  labels <- names(reports)
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels) == 0) {
    stop("No labeled reports in workspace. Run initiate_report() first.")
  }

  # Resolve per-label multipliers without mutating stored $sca permanently
  named_scales <- NULL
  global_scale <- NULL
  if (is.null(sca)) {
    global_scale <- 1
  } else if (is.list(sca) || (is.numeric(sca) && !is.null(names(sca)))) {
    if (is.null(names(sca)) || any(!nzchar(names(sca)))) {
      stop("When sca is a list, all elements must be named (report labels)")
    }
    named_scales <- lapply(as.list(sca), function(x) as.numeric(x)[1])
    unknown <- setdiff(names(named_scales), labels)
    if (length(unknown) > 0) {
      warning("sca names not found in reports: ", paste(unknown, collapse = ", "))
    }
  } else if (is.numeric(sca) && length(sca) == 1L) {
    global_scale <- as.numeric(sca)[1]
    if (!is.finite(global_scale) || global_scale <= 0) {
      stop("sca scalar must be a positive number")
    }
  } else {
    stop("sca must be NULL, a single numeric scalar, or a named list of scalars")
  }

  multiplier_for <- function(label) {
    if (!is.null(named_scales)) {
      if (label %in% names(named_scales)) return(named_scales[[label]])
      return(1)
    }
    global_scale
  }

  results <- list()
  for (label in labels) {
    orig <- reports[[label]]
    if (is.null(orig$sections) || length(orig$sections) == 0) {
      warning("Skipping empty report '", label, "'")
      next
    }

    mult <- multiplier_for(label)
    working <- orig
    if (!identical(as.numeric(mult), 1)) {
      working <- .scale.report.sections(orig, mult)
      message("Report '", label, "': multiplying section $sca by ", mult)
    }

    current <- get("reports", envir = .GlobalEnv)
    current[[label]] <- working
    assign("reports", current, envir = .GlobalEnv)

    report_dir <- collect.report(
      rewrite = rewrite, label = label, res = res, path.type = path.type
    )
    md <- create_report_markdown(
      label = label, report_dir = report_dir, path.type = path.type
    )

    # Restore original $sca only; keep path.type paths and report_dir/rmd
    current <- get("reports", envir = .GlobalEnv)
    collected <- current[[label]]
    if (!identical(as.numeric(mult), 1) &&
        length(orig$sections) > 0 &&
        length(collected$sections) > 0) {
      for (n in seq_along(orig$sections)) {
        if (n > length(collected$sections)) break
        if (is.null(collected$sections[[n]]) || !is.list(collected$sections[[n]])) next
        if (is.null(orig$sections[[n]]) || !is.list(orig$sections[[n]])) next
        collected$sections[[n]]$sca <- orig$sections[[n]]$sca
      }
    }
    current[[label]] <- collected
    assign("reports", current, envir = .GlobalEnv)

    results[[label]] <- md
  }

  invisible(results)
}

#' Knit every report in `reports` to PDF, and optionally HTML
#'
#' Uses the Rmd written by `create_report_markdown()` / `make_all_reports()`.
#' Outputs are written next to that Rmd in the collected report folder, and
#' copied to `~/Rmd/report_markdowns/`.
#'
#' @param pdf Knit PDF (default TRUE)
#' @param html Knit HTML (default FALSE)
knit_all_reports <- function(pdf = TRUE, html = FALSE) {
  if (!isTRUE(pdf) && !isTRUE(html)) {
    stop("Set pdf=TRUE and/or html=TRUE")
  }
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required. install.packages('rmarkdown')")
  }

  reports <- .ensure.reports()
  labels <- names(reports)
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels) == 0) {
    stop("No labeled reports in workspace. Run initiate_report() first.")
  }

  copy_dir <- path.expand("~/Rmd/report_markdowns")
  dir.create(copy_dir, recursive = TRUE, showWarnings = FALSE)
  results <- list()

  for (label in labels) {
    report_obj <- reports[[label]]
    if (is.null(report_obj$sections) || length(report_obj$sections) == 0) {
      warning("Skipping empty report '", label, "'")
      next
    }

    report_dir <- tryCatch(
      .find.report.dir(label, report_obj$report_dir),
      error = function(e) NULL
    )

    rmd <- report_obj$rmd
    if (!is.null(rmd) && nzchar(as.character(rmd)[1])) {
      rmd <- as.character(rmd)[1]
      if (!file.exists(rmd) && !is.null(report_dir)) {
        cand <- file.path(report_dir, rmd)
        if (file.exists(cand)) rmd <- cand
      }
    } else {
      rmd <- NULL
    }

    if (is.null(rmd) || !file.exists(rmd)) {
      if (is.null(report_dir)) {
        warning(
          "Skipping '", label,
          "': no Rmd or collected folder. Run make_all_reports() first."
        )
        next
      }
      md <- create_report_markdown(
        label = label,
        report_dir = report_dir,
        path.type = report_obj$path.type %||% "relative"
      )
      rmd <- md$rmd
    }

    out_dir <- dirname(rmd)
    result <- list(rmd = rmd, pdf = NULL, html = NULL)

    if (isTRUE(pdf)) {
      message("Knitting PDF for report '", label, "' ...")
      result$pdf <- rmarkdown::render(
        input = rmd,
        output_format = rmarkdown::pdf_document(toc = TRUE),
        output_dir = out_dir,
        quiet = TRUE,
        envir = new.env(parent = globalenv())
      )
      file.copy(result$pdf, file.path(copy_dir, basename(result$pdf)), overwrite = TRUE)
      message("Wrote PDF: ", result$pdf)
    }

    if (isTRUE(html)) {
      message("Knitting HTML for report '", label, "' ...")
      result$html <- rmarkdown::render(
        input = rmd,
        output_format = rmarkdown::html_document(
          toc = TRUE,
          toc_float = TRUE,
          self_contained = TRUE
        ),
        output_dir = out_dir,
        quiet = TRUE,
        envir = new.env(parent = globalenv())
      )
      file.copy(result$html, file.path(copy_dir, basename(result$html)), overwrite = TRUE)
      message("Wrote HTML: ", result$html)
    }

    current <- get("reports", envir = .GlobalEnv)
    current[[label]]$pdf <- result$pdf
    current[[label]]$html <- result$html
    assign("reports", current, envir = .GlobalEnv)

    results[[label]] <- result
  }

  invisible(results)
}



.get.report.cache.name <- function(label = NULL) {
  label <- .normalize.report.label(label)
  paste0("report_", label) %>% addversion
}


#' Save every labeled report in `reports` to disk
#'
#' Loops over all names in the global `reports` list and calls
#' `save_report()` for each non-empty report, writing figures plus
#' per-section `report_fig_N.RData` and `report_structure.RData`.
#'
#' @param rewrite,res,path.type Passed to `save_report()` / `collect.report()`.
#' @return Named list of report folder paths (invisibly).
save_all_reports <- function(rewrite = FALSE, res = 600, path.type = "relative") {
  path.type <- match.arg(path.type, c("relative", "absolute"))
  reports <- .ensure.reports()
  labels <- names(reports)
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels) == 0) {
    stop("No labeled reports in workspace. Run initiate_report() first.")
  }

  results <- list()
  for (label in labels) {
    report_obj <- reports[[label]]
    if (is.null(report_obj$sections) || length(report_obj$sections) == 0) {
      warning("Skipping empty report '", label, "'")
      next
    }
    message("Saving report '", label, "' ...")
    results[[label]] <- save_report(
      label = label,
      rewrite = rewrite,
      res = res,
      path.type = path.type
    )
  }

  if (length(results) == 0) {
    warning("No non-empty reports were saved")
  }
  invisible(results)
}


#' Run UMAP on module-score metadata columns and store as a Seurat reduction
#'
#' Finds all `so@meta.data` columns whose names start with `prefix` (default
#' `"mscore_group"`), uses their values as input features, and runs
#' `Seurat::RunUMAP()` with the supplied parameters. The resulting embedding is
#' stored in `so@reductions[[label]]` as a standard `DimReduc` object.
#'
#' @param so Seurat object.
#' @param label Name for the new reduction (default `"mscore_umap"`).
#' @param prefix Metadata column prefix to search for (default `"mscore_group"`).
#' @param assay Assay recorded on the `DimReduc` object; default is
#'   `DefaultAssay(so)`.
#' @param reduction.key Key prefix for UMAP dimensions (e.g. `"mscoreumap_"`).
#'   Default uses `prepare.rk(label)`.
#' @param n.neighbors,min.dist,spread,learning.rate,repulsion.strength,
#'   negative.sample.rate,n.epochs,seed.use,metric,n.components,return.model,
#'   verbose Standard `RunUMAP()` arguments.
#' @param ... Additional arguments passed to `Seurat::RunUMAP()`.
#' @return The Seurat object with `so@reductions[[label]]` set (invisibly).
run_mscore_umap <- function(
    so,
    label = "mscore_umap",
    prefix = "mscore_group",
    assay = NULL,
    reduction.key = NULL,
    n.neighbors = 30,
    min.dist = 0.3,
    spread = 1,
    learning.rate = 1,
    repulsion.strength = 1,
    negative.sample.rate = 5,
    n.epochs = NULL,
    seed.use = 42,
    metric = "cosine",
    n.components = 2,
    return.model = FALSE,
    verbose = TRUE,
    ...
) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required")
  }

  mscore_cols <- grep(paste0("^", prefix), colnames(so@meta.data), value = TRUE)
  if (length(mscore_cols) == 0) {
    stop(
      "No metadata columns starting with '", prefix,
      "' found in so@meta.data"
    )
  }
  mscore_cols <- mixedsort(mscore_cols)

  cells <- colnames(so)
  score_mat <- as.matrix(so@meta.data[cells, mscore_cols, drop = FALSE])
  mode(score_mat) <- "numeric"
  if (anyNA(score_mat)) {
    warning(
      "NA values in ", sum(is.na(score_mat)), " module-score entries; ",
      "imputing column medians for UMAP input"
    )
    for (j in seq_len(ncol(score_mat))) {
      na_idx <- is.na(score_mat[, j])
      if (any(na_idx)) {
        score_mat[na_idx, j] <- stats::median(score_mat[!na_idx, j], na.rm = TRUE)
      }
    }
  }

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(so)
  }
  if (is.null(reduction.key)) {
    reduction.key <- prepare.rk(label)
  }

  tmp_red <- paste0(".", label, "_features")
  so[[tmp_red]] <- Seurat::CreateDimReducObject(
    embeddings = score_mat,
    key = paste0(tmp_red, "_"),
    assay = assay
  )

  umap_args <- c(
    list(
      object = so,
      reduction = tmp_red,
      dims = seq_len(ncol(score_mat)),
      n.neighbors = n.neighbors,
      min.dist = min.dist,
      spread = spread,
      learning.rate = learning.rate,
      repulsion.strength = repulsion.strength,
      negative.sample.rate = negative.sample.rate,
      n.epochs = n.epochs,
      seed.use = seed.use,
      metric = metric,
      n.components = n.components,
      reduction.name = label,
      reduction.key = reduction.key,
      return.model = return.model,
      verbose = verbose
    ),
    list(...)
  )

  so <- do.call(Seurat::RunUMAP, umap_args)
  so@reductions[[tmp_red]] <- NULL

  if (verbose) {
    message(
      "Stored UMAP on ", length(mscore_cols), " module-score column(s) as '",
      label, "' (", paste(head(mscore_cols, 3), collapse = ", "),
      if (length(mscore_cols) > 3) ", ..." else "", ")"
    )
  }

  invisible(so)
}

# Clean a list of marker vectors: drop empty sets, dedupe within sets, name unnamed entries.
.normalize.marker.list <- function(marker.list) {
  if (!is.list(marker.list)) {
    stop("marker.list must be a list of character vectors")
  }
  if (is.null(names(marker.list)) || any(!nzchar(names(marker.list)))) {
    if (is.null(names(marker.list))) {
      names(marker.list) <- paste0("set_", seq_along(marker.list))
    } else {
      empty_names <- which(!nzchar(names(marker.list)))
      names(marker.list)[empty_names] <- paste0("set_", empty_names)
    }
  }
  marker.list <- lapply(marker.list, function(x) {
    x <- unique(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    x
  })
  empty <- lengths(marker.list) == 0
  if (any(empty)) {
    warning("Dropping empty marker set(s): ", paste(names(marker.list)[empty], collapse = ", "))
    marker.list <- marker.list[!empty]
  }
  if (length(marker.list) == 0) {
    stop("marker.list has no non-empty sets")
  }
  marker.list
}

#' Set overlap analysis for a list of marker gene vectors
#'
#' Summarizes pairwise and higher-order intersections among named marker lists,
#' and counts genes unique to each set.
#'
#' @param marker.list Named list of character vectors (gene symbols).
#' @param include.genes If TRUE, attach gene vectors to combination outputs.
#' @param min.combination.size Minimum gene count to keep in `combinations`.
#' @return A list with:
#'   \item{sizes}{Named integer vector of set sizes.}
#'   \item{unique_counts}{Named integer vector of genes unique to each set.}
#'   \item{unique_genes}{Named list of unique gene vectors per set.}
#'   \item{intersection_all}{Genes present in every set.}
#'   \item{intersection_all_count}{Length of `intersection_all`.}
#'   \item{pairwise}{Data frame of pairwise intersections.}
#'   \item{combinations}{Data frame of all non-empty membership patterns.}
#'   \item{combination_genes}{Named list of genes per pattern (if `include.genes`).}
#'   \item{membership}{Logical matrix: genes x sets.}
analyze.marker.sets <- function(marker.list, include.genes = TRUE,
                                min.combination.size = 1) {
  marker.list <- .normalize.marker.list(marker.list)
  set_names <- names(marker.list)
  n_sets <- length(set_names)

  all_genes <- unique(unlist(marker.list, use.names = FALSE))
  membership <- matrix(
    FALSE,
    nrow = length(all_genes),
    ncol = n_sets,
    dimnames = list(all_genes, set_names)
  )
  for (sn in set_names) {
    hits <- intersect(all_genes, marker.list[[sn]])
    if (length(hits) > 0) {
      membership[hits, sn] <- TRUE
    }
  }

  # Genes unique to exactly one set
  unique_genes <- lapply(set_names, function(sn) {
    all_genes[membership[, sn] & rowSums(membership) == 1L]
  })
  names(unique_genes) <- set_names

  intersection_all <- if (n_sets >= 2) {
    Reduce(intersect, marker.list)
  } else {
    marker.list[[1]]
  }

  # Membership pattern per gene (e.g. "A&B", "A&C", "B")
  patterns <- apply(membership, 1, function(row) {
    paste(set_names[row], collapse = "&")
  })
  patterns <- as.character(patterns)
  patterns[!nzchar(patterns)] <- "(none)"
  names(patterns) <- rownames(membership)

  # n_sets per pattern from membership (robust if set names contain "&")
  n_sets_by_pattern <- vapply(
    split(rowSums(membership), patterns),
    function(x) as.integer(x[[1]]),
    integer(1)
  )

  pat_table <- sort(table(patterns, useNA = "no"), decreasing = TRUE)
  if (length(pat_table) == 0) {
    combinations <- data.frame(
      pattern = character(),
      n_sets = integer(),
      n_genes = integer(),
      stringsAsFactors = FALSE
    )
    combination_genes <- list()
  } else {
    pat <- names(pat_table)
    combinations <- data.frame(
      pattern = pat,
      n_sets = unname(n_sets_by_pattern[pat]),
      n_genes = as.integer(unname(pat_table)),
      stringsAsFactors = FALSE
    )
    combinations <- combinations[combinations$n_genes >= min.combination.size, , drop = FALSE]
    rownames(combinations) <- NULL

    combination_genes <- split(names(patterns), patterns, simplify = FALSE)
    combination_genes <- lapply(combination_genes, unlist, use.names = FALSE)
    combination_genes <- combination_genes[combinations$pattern]
  }

  # Pairwise intersections
  pairwise <- if (n_sets >= 2) {
    pairs <- utils::combn(set_names, 2, simplify = FALSE)
    do.call(
      rbind,
      lapply(pairs, function(pair) {
        s1 <- pair[[1]]
        s2 <- pair[[2]]
        inter <- intersect(marker.list[[s1]], marker.list[[s2]])
        uni <- union(marker.list[[s1]], marker.list[[s2]])
        data.frame(
          set1 = s1,
          set2 = s2,
          n_intersection = length(inter),
          n_union = length(uni),
          jaccard = if (length(uni) > 0) length(inter) / length(uni) else NA_real_,
          stringsAsFactors = FALSE
        )
      })
    )
  } else {
    data.frame(
      set1 = character(),
      set2 = character(),
      n_intersection = integer(),
      n_union = integer(),
      jaccard = numeric(),
      stringsAsFactors = FALSE
    )
  }
  if (include.genes && nrow(pairwise) > 0) {
    pairwise$intersection_genes <- lapply(
      seq_len(nrow(pairwise)),
      function(i) intersect(marker.list[[pairwise$set1[i]]], marker.list[[pairwise$set2[i]]])
    )
  }

  sizes <- vapply(marker.list, length, integer(1))
  unique_counts <- vapply(unique_genes, length, integer(1))

  out <- list(
    sizes = sizes,
    unique_counts = unique_counts,
    unique_genes = unique_genes,
    intersection_all = intersection_all,
    intersection_all_count = length(intersection_all),
    pairwise = pairwise,
    combinations = combinations,
    membership = membership
  )
  if (include.genes) {
    out$combination_genes <- combination_genes
  }

  class(out) <- c("marker_set_analysis", "list")
  out
}

#' Print summary for `analyze.marker.sets()` output
print.marker_set_analysis <- function(x, ...) {
  cat("Marker set analysis (", length(x$sizes), " sets, ",
      sum(x$sizes), " total entries, ",
      length(unique(unlist(x$unique_genes, use.names = FALSE))), " distinct genes)\n",
      sep = "")
  cat("\nSet sizes:\n")
  print(x$sizes)
  cat("\nUnique genes per set:\n")
  print(x$unique_counts)
  cat("\nGenes in intersection of all sets:", x$intersection_all_count, "\n")
  if (nrow(x$combinations) > 0) {
    cat("\nTop membership patterns:\n")
    print(utils::head(x$combinations, 10))
  }
  if (nrow(x$pairwise) > 0) {
    cat("\nPairwise intersections:\n")
    print(x$pairwise[, setdiff(colnames(x$pairwise), "intersection_genes"), drop = FALSE])
  }
  invisible(x)
}


.scvi.script.path <- function(python_script = NULL) {
  if (!is.null(python_script)) {
    python_script <- path.expand(python_script)
    if (!file.exists(python_script)) {
      stop("Python script not found: ", python_script)
    }
    return(normalizePath(python_script))
  }

  candidates <- c(
    file.path(config$project_root, "python", "run_scvi.py"),
    file.path(getwd(), "python", "run_scvi.py"),
    file.path(getwd(), "..", "python", "run_scvi.py")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop(
      "Could not find python/run_scvi.py. Pass python_script= explicitly or ",
      "set config$project_root to the micropatterns repo root."
    )
  }
  normalizePath(found[1])
}


#' Run scVI via `python/run_scvi.py` using the reticulate conda env
#'
#' Trains scVI on an h5ad (e.g. from `save.cellxgene()` / `save.cellxgene2()`)
#' and writes latent embeddings to CSV and/or an updated h5ad. Uses
#' `reticulate::use_condaenv()` so the script runs in the same Python as other
#' reticulate chunks in R Markdown.
#'
#' @param input_h5ad Path to input `.h5ad`.
#' @param output_csv Optional path for scVI latent embedding CSV.
#' @param output_h5ad Optional path for output `.h5ad` with `X_scVI` in `.obsm`.
#' @param layer Expression layer for scVI (default `"RNA"`).
#' @param n_top_genes Number of highly variable genes (default 5000).
#' @param batch_key `batch_key` for HVG selection (default `"orig_ident"`).
#' @param categorical_covariates Categorical obs columns for scVI.
#' @param continuous_covariates Continuous obs columns for scVI.
#' @param max_epochs Optional training epoch limit (`NULL` = scVI default).
#' @param conda_env Conda environment name (default `"scvi"`).
#' @param python_script Path to `run_scvi.py` (auto-detected if `NULL`).
#' @param load_csv If `TRUE` and `output_csv` is set, return the embedding
#'   `data.frame` invisibly (still returns `output_csv` path as visible value).
#'
#' @return Invisibly, the path(s) to output file(s) as a character vector.
#'   If `load_csv = TRUE`, the embedding table is returned as the visible value.
#'
#' @examples
#' \dontrun{
#' # In R Markdown after source("R/mlutils_setup.R"):
#' library(reticulate)
#' use_condaenv("scvi", required = TRUE)
#'
#' h5ad <- save.cellxgene2(so, name = "mydata", path = config$out_root)
#' emb <- run_scvi(
#'   input_h5ad = h5ad,
#'   output_csv = file.path(config$out_root, "mydata_scvi.csv"),
#'   max_epochs = 50L,
#'   load_csv = TRUE
#' )
#' so <- add_scvi_to_seurat(so, emb)
#' }
run_scvi <- function(
  so,
  output_csv = "~/mnt_out/scvi_output.csv",
  output_h5ad = NULL,
  layer = "RNA",
  n_top_genes = 5000L,
  batch_key = "orig_ident",
  categorical_covariates = c("experiment", "day"),
  continuous_covariates = c("percent_mt", "S_Score", "G2M_Score"),
  max_epochs = NULL,
  conda_env = "scvi",
  python_script = NULL,
  load_csv = FALSE, 
  skip_write=FALSE
) {
  
  
  input_h5ad=path.expand("~/mnt_out/temp.h5ad")
  #input_h5ad <- normalizePath(path.expand(input_h5ad), mustWork = TRUE)
  
if(!skip_write){  
fcat("exporting h5ad object...")  
save.cellxgene2(so, name = "temp", path=config$out_root, overwrite = T, assay=DefaultAssay(so))
}
  
  if (!is.null(output_csv)) {
    output_csv <- path.expand(output_csv)
    dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  }
  if (!is.null(output_h5ad)) {
    output_h5ad <- path.expand(output_h5ad)
    dir.create(dirname(output_h5ad), recursive = TRUE, showWarnings = FALSE)
  }

  #reticulate::use_condaenv(conda_env, required = TRUE)
  python <- reticulate::py_config()$python
  script <- .scvi.script.path(python_script)

  args <- c(
    script,
    "-i", input_h5ad,
    #"--layer", layer,
    "--n-top-genes", as.character(as.integer(n_top_genes)),
    "--batch-key", batch_key
  )

  if (!is.null(output_csv)) {
    args <- c(args, "--output-csv", output_csv)
  }
  if (!is.null(output_h5ad)) {
    args <- c(args, "--output-h5ad", output_h5ad)
  }
  if (length(categorical_covariates)) {
    args <- c(args, "--categorical-covariates")
    args <- c(args, categorical_covariates)
  }
  if (length(continuous_covariates)) {
    args <- c(args, "--continuous-covariates")
    args <- c(args, continuous_covariates)
  }
  if (!is.null(max_epochs)) {
    args <- c(args, "--max-epochs", as.character(as.integer(max_epochs)))
  }

  fcat("Running scVI:", python, paste(shQuote(args), collapse = " "))
  status <- system2(python, args, stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop("run_scvi.py failed with exit code ", status)
  }

  outputs <- c(
    if (!is.null(output_csv)) output_csv,
    if (!is.null(output_h5ad)) output_h5ad
  )

  if (isTRUE(load_csv) && !is.null(output_csv)) {
    emb <- utils::read.csv(output_csv, row.names = 1, check.names = FALSE)
    return(emb)
  }

  #invisible(outputs)
  invisible(add.embedding(so, file=output_csv, reduction_name = "scvi"))
}


#' Add scVI latent embeddings to a Seurat object
#'
#' @param so Seurat object.
#' @param embedding Embedding matrix or data.frame (cells x scvi_*), or path to
#'   CSV written by `run_scvi()`.
#' @param reduction_name Name for the new DimReduc (default `"scvi"`).
#' @param key Prefix for embedding column names (default `"scvi_"`).
#' @param assay Assay for the DimReduc object (default: `DefaultAssay(so)`).
#' @param reorder If `TRUE`, reorder embedding rows to match `colnames(so)`.
add_scvi_to_seurat <- function(
  so,
  embedding,
  reduction_name = "scvi",
  key = "scvi_",
  assay = NULL,
  reorder = TRUE
) {
  if (is.character(embedding) && length(embedding) == 1L) {
    embedding <- utils::read.csv(embedding, row.names = 1, check.names = FALSE)
  }
  lat <- as.matrix(embedding)

  if (isTRUE(reorder)) {
    common <- intersect(colnames(so), rownames(lat))
    if (length(common) == 0) {
      stop("No overlapping cell names between Seurat object and scVI embedding.")
    }
    if (length(common) < ncol(so)) {
      warning(
        length(common), " of ", ncol(so),
        " Seurat cells matched scVI embeddings; subsetting Seurat object.",
        call. = FALSE
      )
      so <- so[, common]
    }
    lat <- lat[common, , drop = FALSE]
  }

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(so)
  }

  so[[reduction_name]] <- SeuratObject::CreateDimReducObject(
    embeddings = lat,
    key = key,
    assay = assay
  )
  so
}

make_scvi_bash_script <- function(
  label,
  path_script = file.path(config$project_root_host, "python"),
  path_out = config$out_root_host
) {
  label <- as.character(label)[1]
  script_dir <- "/home/rstudio/bash"
  dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
  script_file <- file.path(script_dir, paste0("run_scvi_", label, ".sh"))
  output_csv <- file.path(path_out, paste0("scvi_embedding_", label, ".csv"))
  cmd <- paste(
    "python",
    file.path(path_script, "run_scvi.py"),
    "-i",
    file.path(path_script, "temp.h5ad"),
    "--output-csv",
    output_csv
  )

  writeLines(c("#!/usr/bin/env bash", cmd), script_file)
  Sys.chmod(script_file, mode = "0755")
  fcat("Wrote:", script_file)
  output_csv
}

