

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width="100%",
  fig.width=7, fig.height=5,
  dpi=300, fig.path="figures/BayesSpace-",
  message=FALSE, warning=FALSE, error=FALSE
)

## ----setup--------------------------------------------------------------------
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)

## ----readVisium, eval=FALSE---------------------------------------------------
#  sce <- readVisium("path/to/spaceranger/outs/")

## ----download-----------------------------------------------------------------
melanoma <- getRDS(dataset="2018_thrane_melanoma", sample="ST_mel1_rep2")

## ----manual.sce, eval=FALSE---------------------------------------------------
#  library(Matrix)
#  
rowData <- read.csv("D:/Xenium_benchmarking-main/Xenium_benchmarking-main/data/formatted_for_R/var_msbrain_1.csv", stringsAsFactors=FALSE)
colData <- read.csv("D:/Xenium_benchmarking-main/Xenium_benchmarking-main/data/formatted_for_R/obs_msbrain_1.csv", stringsAsFactors=FALSE, row.names=1)
counts <- read.csv("D:/Xenium_benchmarking-main/Xenium_benchmarking-main/data/formatted_for_R/exp_msbrain_1.csv",row.names=1, check.names=F, stringsAsFactors=FALSE)
#  
rownames(counts)<- rownames(colData)
colnames(counts)<-rownames(rowData)


sce <- SingleCellExperiment(assays=list(counts=as(as.matrix(counts), "dgCMatrix")),
                              rowData=colData,
                              colData=rowData)

melanoma<-sc
  
## ----preprocess---------------------------------------------------------------
set.seed(102)
melanoma <- spatialPreprocess(melanoma, platform="ST", 
                              n.PCs=7, n.HVGs=2000, log.normalize=FALSE)

## ----tuning_q-----------------------------------------------------------------
melanoma <- qTune(melanoma, qs=seq(2, 10), platform="ST", d=7)
qPlot(melanoma)

## ----cluster------------------------------------------------------------------
set.seed(149)
melanoma <- spatialCluster(melanoma, q=4, platform="ST", d=7,
                           init.method="mclust", model="t", gamma=2,
                           nrep=1000, burn.in=100,
                           save.chain=TRUE)

## ----cluster.results----------------------------------------------------------
head(colData(melanoma))

## ----cluster.plot, fig.width=7, fig.height=5----------------------------------
clusterPlot(melanoma)

## ----cluster.plot.customize, fig.width=7, fig.height=5------------------------
clusterPlot(melanoma, palette=c("purple", "red", "blue", "yellow"), color="black") +
  theme_bw() +
  xlab("Column") +
  ylab("Row") +
  labs(fill="BayesSpace\ncluster", title="Spatial clustering of ST_mel1_rep2")

## ----enhance, eval=TRUE-------------------------------------------------------
melanoma.enhanced <- spatialEnhance(melanoma, q=4, platform="ST", d=7,
                                    model="t", gamma=2,
                                    jitter_prior=0.3, jitter_scale=3.5,
                                    nrep=1000, burn.in=100,
                                    save.chain=TRUE)

## ----enhance.results----------------------------------------------------------
head(colData(melanoma.enhanced))

## ----enhance.plot, eval=TRUE, fig.width=7, fig.height=5-----------------------
clusterPlot(melanoma.enhanced)

## ----enhanceFeatures----------------------------------------------------------
markers <- c("PMEL", "CD2", "CD19", "COL1A1")
melanoma.enhanced <- enhanceFeatures(melanoma.enhanced, melanoma,
                                     feature_names=markers,
                                     nrounds=0)

## ----enhanced.logcount--------------------------------------------------------
logcounts(melanoma.enhanced)[markers, 1:5]

## ----enhanced.rmse------------------------------------------------------------
rowData(melanoma.enhanced)[markers, ]

## ----enhanced.featurePlot-----------------------------------------------------
featurePlot(melanoma.enhanced, "PMEL")

## ----enhanced.markers, fig.width=12, fig.height=8-----------------------------
enhanced.plots <- purrr::map(markers, function(x) featurePlot(melanoma.enhanced, x))
patchwork::wrap_plots(enhanced.plots, ncol=2)

## ----compare.resolution, fig.width=16, fig.height=8---------------------------
spot.plots <- purrr::map(markers, function(x) featurePlot(melanoma, x))
patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=4)

## ----mcmcChain, eval=TRUE-----------------------------------------------------
chain <- mcmcChain(melanoma)
chain[1:5, 1:5]
