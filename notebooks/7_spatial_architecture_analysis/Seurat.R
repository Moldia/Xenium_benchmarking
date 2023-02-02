library("Seurat")


InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")



#  
rowData <- read.csv("D:/Xenium_benchmarking-main/Xenium_benchmarking-main/data/formatted_for_R/var_msbrain_1.csv", stringsAsFactors=FALSE)
colData <- read.csv("D:/Xenium_benchmarking-main/Xenium_benchmarking-main/data/formatted_for_R/obs_msbrain_1.csv", stringsAsFactors=FALSE, row.names=1)
counts <- read.csv("D:/Xenium_benchmarking-main/Xenium_benchmarking-main/data/formatted_for_R/exp_msbrain_1.csv",row.names=1, check.names=F, stringsAsFactors=FALSE)
#  
rownames(counts)<- rownames(colData)
rownames(rowData)<-rowData$index
colnames(counts)<-rownames(rowData)


coord.df = data.frame(x=colData$x_centroid, y=colData$y_centroid)
rownames(coord.df) = rownames(colData)


brain<-CreateSeuratObject(counts, project = "SeuratProject", assay = "Spatial",
                   min.cells = 0, min.features = 0, names.field = 1,
                   names.delim = "_", meta.data =colData)

brain@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = coord.df
)


brain <- FindSpatiallyVariableFeatures(brain, assay = "Spatial",
                                       selection.method = "markvariogram")