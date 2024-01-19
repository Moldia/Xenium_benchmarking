library("Seurat")


tag='human_brain'
#  


files=list.files('/media/sergio/xenium_b_and_heart/Xenium_benchmarking-main/data/formatted_for_R/unprocessed_adata_nuclei/')

tag<-"ms_brain_multisection1"
cells<-c(100,500,1000,5000)

resui<-c()
for (cell in cells){
print(cell)
start<-Sys.time()


print(tag)
rowData <- read.csv(paste0('/media/sergio/xenium_b_and_heart/Xenium_benchmarking-main/data/formatted_for_R/unprocessed_adata_nuclei/',tag,'/var.csv'), stringsAsFactors=FALSE)
colData <- read.csv(paste0('/media/sergio/xenium_b_and_heart/Xenium_benchmarking-main/data/formatted_for_R/unprocessed_adata_nuclei/',tag,'/obs.csv'), stringsAsFactors=FALSE, row.names=1)
counts <- read.csv(paste0('/media/sergio/xenium_b_and_heart/Xenium_benchmarking-main/data/formatted_for_R/unprocessed_adata_nuclei/',tag,'/exp.csv'),row.names=1, check.names=F, stringsAsFactors=FALSE)
#  

smpl=sample(rownames(colData),cell)
counts=counts[smpl,]
colData=colData[smpl,]

rownames(rowData)<-rowData$gene_id
rownames(counts)<- rownames(colData)
colnames(counts)<-rownames(rowData)

coord.df = data.frame(x=colData$x_centroid, y=colData$y_centroid)
rownames(coord.df) = rownames(colData)
brain<-CreateSeuratObject(t(counts), project = "SeuratProject", assay = "Spatial",
                          min.cells = 0, min.features = 0, names.field = 0,names.delim = "_",
                           meta.data =as.data.frame(colData))#,

brain@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coord.df
)

brain <- NormalizeData(brain)
all.genes <- rownames(brain)
brain <- ScaleData(brain, features = all.genes)

brain <- FindSpatiallyVariableFeatures(brain, r.metric = 5,
                                       x.cuts = NULL,
                                       y.cuts = NULL,
                                       verbose = TRUE,assay = "Spatial",
                                       selection.method = "moransi")

spat_variable=brain@assays$Spatial@meta.features
spat_variable=spat_variable[,c("MoransI_p.value","moransi.spatially.variable","moransi.spatially.variable.rank")]
colnames(spat_variable)<-c('pvalue','SVG','rank')





finish<-Sys.time()

timespent<-finish-start
timespent<-as.numeric(timespent)
print(timespent)
resui<-append(resui,timespent)
}


resu2<-c()
for (cell in cells){
  print(cell)
  start<-Sys.time()
  
  
  print(tag)
  rowData <- read.csv(paste0('/media/sergio/xenium_b_and_heart/Xenium_benchmarking-main/data/formatted_for_R/unprocessed_adata_nuclei/',tag,'/var.csv'), stringsAsFactors=FALSE)
  colData <- read.csv(paste0('/media/sergio/xenium_b_and_heart/Xenium_benchmarking-main/data/formatted_for_R/unprocessed_adata_nuclei/',tag,'/obs.csv'), stringsAsFactors=FALSE, row.names=1)
  counts <- read.csv(paste0('/media/sergio/xenium_b_and_heart/Xenium_benchmarking-main/data/formatted_for_R/unprocessed_adata_nuclei/',tag,'/exp.csv'),row.names=1, check.names=F, stringsAsFactors=FALSE)
  #  
  
  smpl=sample(rownames(colData),cell)
  counts=counts[smpl,]
  colData=colData[smpl,]
  
  rownames(rowData)<-rowData$gene_id
  rownames(counts)<- rownames(colData)
  colnames(counts)<-rownames(rowData)
  
  coord.df = data.frame(x=colData$x_centroid, y=colData$y_centroid)
  rownames(coord.df) = rownames(colData)
  brain<-CreateSeuratObject(t(counts), project = "SeuratProject", assay = "Spatial",
                            min.cells = 0, min.features = 0, names.field = 0,names.delim = "_",
                            meta.data =as.data.frame(colData))#,
  
  brain@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = coord.df
  )
  
  brain <- NormalizeData(brain)
  all.genes <- rownames(brain)
  brain <- ScaleData(brain, features = all.genes)
  
  
  brain <- FindSpatiallyVariableFeatures(brain, r.metric = 5,
                                         x.cuts = NULL,
                                         y.cuts = NULL,
                                         verbose = TRUE,assay = "Spatial",
                                         selection.method = "markvariogram")
  
  spat_variable=brain@assays$Spatial@meta.features
  colnames(spat_variable)<-c('rstat','SVG','rank')
  
  
  finish<-Sys.time()
  
  timespent<-finish-start
  timespent<-as.numeric(timespent)
  print(timespent)
  resu2<-append(resu2,timespent)
}

resui[4]=resui[4]*60
resu2[3]=resu2[3]*60
resu2[4]=resu2[4]*60

allresults<-data.frame(cells,resui,resu2)
colnames(allresults)<-c('cells','time_seurat_i','time_seurat_markvariogram')




write.csv(allresults,'/media/sergio/xenium_b_and_heart/Xenium_benchmarking-main/figures/times_svf/seurat_times.csv')


#SpatiallyVariableFeatures(brain,selection.method = 'markvariogram')