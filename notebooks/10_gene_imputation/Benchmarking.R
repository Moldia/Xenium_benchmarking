# Set working directory
setwd('D:/MolDiaFiles/GenePrediction/SpatialBenchmarking/')

# Load libraries
library(Seurat)
library(ggplot2)
library(reticulate)
np <- import("numpy")

# File paths
PATHDIR <- "../DataPreparation/newData/"
RNAfile <- "../DataPreparation/newData/scRNA_count.txt"
Spatialfile <- "../DataPreparation/newData/Insitu_count.txt"

# Read scRNA data
RNA <- read.table(RNAfile, sep = '\t', header = TRUE, row.names = 1, quote = "")
RNA <- CreateSeuratObject(counts = RNA, project = 'RNA', min.cells = 0, min.features = 0)

# Read spatial data
Spatial_orig <- t(read.table(Spatialfile, sep = '\t', header = TRUE, quote = ""))
Genes <- scan(Spatialfile, what = 'character', nlines = 1)
rownames(Spatial_orig) <- Genes
colnames(Spatial_orig) <- paste0(colnames(Spatial_orig), 1:ncol(Spatial_orig))

# Initialize result dataframe
Result <- as.data.frame(array(, dim = c(dim(Spatial_orig)[2], dim(Spatial_orig)[1])))
colnames(Result) <- Genes
rownames(Result) <- 1:ncol(Spatial_orig)
Result <- t(Result)

# Load train and test lists
train_list <- t(np$load(paste0(PATHDIR, 'train_list.npy'), allow_pickle = TRUE))
test_list <- t(np$load(paste0(PATHDIR, 'test_list.npy'), allow_pickle = TRUE))

# Define imputation function
run_imputation <- function(i) {
  genes.leaveout <- unlist(train_list[, i])
  feature.remove <- unlist(test_list[, i])
  print('We Used Test Genes : ')
  print(feature.remove)
  features <- unlist(train_list[, i])
  print(length(features))
  Spatial <- Spatial_orig[features, ]
  print(dim(Spatial))
  Spatial <- CreateSeuratObject(counts = Spatial, project = 'Spatial', min.cells = 0, min.features = 0)
  DN <- 30
  if ((length(features) - 1) < 30) {
    DN <- (length(features) - 1)
  }
  anchors <- FindTransferAnchors(reference = RNA, query = Spatial, features = features,
                                 reduction = 'cca', reference.assay = 'RNA', query.assay = 'RNA',
                                 k.filter = NA, dims = 1:DN)
  refdata <- GetAssayData(object = RNA, assay = 'RNA', slot = 'data')
  print('run Transfer')
  imputation <- TransferData(anchorset = anchors, refdata = refdata, weight.reduction = 'pca', dims = 1:DN)
  options(warn = -1)
  Imp_New_genes <- as.data.frame(imputation@data)[feature.remove, ]
  return(Imp_New_genes)
}

# Create directory for results
dir.create("XeniumResults/")

# Run imputation for each test case
for (i in 1:10) {
  res <- Result[unlist(test_list[, i]), ]
  res_emp <- apply(res, 1, function(row) any(is.na(row)))
  Result[unlist(test_list[, i]), ][res_emp, ] <- as.matrix(run_imputation(i))[res_emp, ]
}

# Write results to file
write.table(t(Result), paste0("XeniumResults/", 'Seurat_impute.csv'), sep = ',', quote = FALSE)

# Turn off warnings
warnings('off')
