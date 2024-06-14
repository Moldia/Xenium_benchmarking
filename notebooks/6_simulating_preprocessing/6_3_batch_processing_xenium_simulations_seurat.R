library("Seurat")

path='/media/sergio/xenium_b_and_heart/Xenium_benchmarking-main/data/scRNAseq_for_simulations_Xenium_like'
fils=list.files(path)

######## FUNCTIONS#####################
preprocessing<-function(brain,pcs=30,neighs=15,norm='SCT',scale='True',hvg='True',logt='True',ts=100,hvg_method='vst'){
if (norm=='SCT'){
  brain <- SCTransform(brain,method='glmGamPoi', verbose = FALSE,assay='Spatial')
}else{
  if (logt=='True'){
    brain <- NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = ts)
  }
  brain <- FindVariableFeatures(brain, selection.method = hvg_method, nfeatures = round(dim(brain)[1]*0.7))
  
  if (hvg=='True'){
    all.genes <- VariableFeatures(object=brain)
  } else{
    all.genes <- rownames(brain)
  }
  if (scale=='True'){
    brain <- ScaleData(brain, features = all.genes) 
  }
}
###WORK ON CLUSTERING#######
brain <- RunPCA(brain ,npcs = pcs, verbose = FALSE)
brain<- RunUMAP(brain, dims = 1:pcs, verbose = FALSE)
brain <- FindNeighbors(brain, dims = 1:pcs, verbose = FALSE)
#####LOUVAIN CLUSTERING######################
res=0.8
brain <- FindClusters(brain, verbose = FALSE,resolution=res,algorithm=1)
DimPlot(brain,reduction='umap',group.by='seurat_clusters')
###COMPARISON BETWEEN CLUSTERS
ncts=length(unique(brain@meta.data$cell_type))
current=length(unique(brain@meta.data$seurat_clusters))
# WE NOW ITERATE ON THE CLUSTERING RESOLUTION TO MATCH THE CELL TYPES
gap=0.1
while (abs(current-ncts)>3)
{
  if (current>ncts){
    res=res-gap
    brain <- FindClusters(brain, verbose = FALSE,resolution=res)
    current=length(unique(brain@meta.data$seurat_clusters))
  }
  else{
    res=res+gap
    brain <- FindClusters(brain, verbose = FALSE,resolution=res)
    current=length(unique(brain@meta.data$seurat_clusters))
  }
}
ctname<-paste0('norm',toString(norm),'_lg',toString(logt),'_ng',toString(pcs),'_pc',toString(pcs),'_ts',toString(ts),'_hvg_',toString(hvg),'_scale_',toString(scale),'_louvain')
result<-data.frame(brain@meta.data$seurat_clusters)
colnames(result)<-c(ctname)
allres<-result
#####SLM CLUSTERING######################
res=0.8
brain <- FindClusters(brain, verbose = FALSE,resolution=res,algorithm=3)
DimPlot(brain,reduction='umap',group.by='seurat_clusters')
###COMPARISON BETWEEN CLUSTERS
ncts=length(unique(brain@meta.data$cell_type))
current=length(unique(brain@meta.data$seurat_clusters))
# WE NOW ITERATE ON THE CLUSTERING RESOLUTION TO MATCH THE CELL TYPES
gap=0.1
while (abs(current-ncts)>3)
{
  if (current>ncts){
    res=res-gap
    brain <- FindClusters(brain, verbose = FALSE,resolution=res)
    current=length(unique(brain@meta.data$seurat_clusters))
  }
  else{
    res=res+gap
    brain <- FindClusters(brain, verbose = FALSE,resolution=res)
    current=length(unique(brain@meta.data$seurat_clusters))
  }
}
ctname<-paste0('norm',toString(norm),'_lg',toString(logt),'_ng',toString(pcs),'_pc',toString(pcs),'_ts',toString(ts),'_hvg_',toString(hvg),'_scale_',toString(scale),'_SLM')
result<-data.frame(brain@meta.data$seurat_clusters)
colnames(result)<-c(ctname)
allres=cbind(allres,result)
return(allres)
}


############ Start here###############
for (i in c(24:40)){
tag=fils[i]
print(tag)

#  
rowData <- read.csv(paste0(path,'/',tag,'/var.csv'), stringsAsFactors=FALSE)
colData <- read.csv(paste0(path,'/',tag,'/obs.csv'), stringsAsFactors=FALSE, row.names=1)
counts <- read.csv(paste0(path,'/',tag,'/exp.csv'),row.names=1, check.names=F, stringsAsFactors=FALSE)
rownames(rowData)<-rowData$gene_id
rownames(counts)<- rownames(colData)
colnames(counts)<-rownames(rowData)
counts[is.na(counts)]=0
counts[counts==-1]=0

pcs=30#dim(counts)[2]
neighs=15
norm='SCT'
scale='True'
hvg='True'
logt='True'
ts=100
hvg_method='vst'#'mvp','disp'


pcs=c(10,30,50)
tss=c(10,100,1000)
nopts=c(8,12,20)
scaleopt=c('True','False')
logopt=c('True','False')
hvgopt=c('vst','mvp','disp')
normt=c('SCT','True')
hvgTF=c('True','False')

iterative=0
for (pcopt in pcs){
  for (ts in tss){
    for (hvgtfs in hvgTF){
    for (lg in logopt){
      for (nro in normt){
      for (sca in scaleopt){
        for (hvgo in hvgopt){
          for (ngo in nopts){
  iterative=iterative+1
  print(toString(iterative))
  brain<-CreateSeuratObject(t(counts), project = "SeuratProject", assay = "Spatial",
                            min.cells = 0, min.features = 0, names.field = 0,names.delim = "_",
                            meta.data =as.data.frame(colData))#,
  if (nro=='SCT'){
    if (sca=='False' & lg=='False' & hvgo=='vst' & hvgtfs=='True' & ts==100){
      try(allres<-preprocessing(brain,pcs=pcopt,neighs=ngo,norm=nro,scale=sca,hvg=hvgtfs,logt=lg,ts=ts,hvg_method=hvgo))
      if (exists("allrescol")==FALSE){
        allrescol=allres
      } else{
        allrescol=cbind(allrescol,allres)
      }
    }
  }else{
    if (hvgtfs=='False'){
      if (hvgo=='vst'){
        try(allres<-preprocessing(brain,pcs=pcopt,neighs=ngo,norm=nro,scale=sca,hvg=hvgtfs,logt=lg,ts=ts,hvg_method=hvgo))
        if (exists("allrescol")==FALSE){
          allrescol=allres
        } else{
          allrescol=cbind(allrescol,allres)
        }
      }
    }else{
      try(allres<-preprocessing(brain,pcs=pcopt,neighs=ngo,norm=nro,scale=sca,hvg=hvgtfs,logt=lg,ts=ts,hvg_method=hvgo))
      if (exists("allrescol")==FALSE){
        allrescol=allres
      } else{
        allrescol=cbind(allrescol,allres)
      }
    }}}
    if (exists("allrescol")==TRUE){
    print(dim(allrescol))  
          }}}}}}}}

write.csv(allrescol,paste0(path,'/',tag,'/results_clustering.csv'))
rm('allrescol')
}





  
DimPlot(brain,reduction='umap',group.by='cell_type')
DimPlot(brain, reduction = "umap")


all.genes <- rownames(brain)
brain <- ScaleData(brain, features = all.genes)
