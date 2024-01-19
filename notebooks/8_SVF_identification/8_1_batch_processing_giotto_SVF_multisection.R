
library(Giotto)
library(smfishHmrf)
library(Matrix)

direct='D:/Xenium_benchmarking-main/data/formatted_for_R/unprocessed_adata_nuclei/'

files=list.files(direct)

tag<-"ms_brain_multisection1"
cells<-c(100,500,1000,5000)

cell=5000

resui<-c()
for (tag in files){
print(tag)
start<-Sys.time()

print(tag)
rowData <- read.csv(paste0(direct,tag,'/var.csv'), stringsAsFactors=FALSE)
colData <- read.csv(paste0(direct,tag,'/obs.csv'), stringsAsFactors=FALSE, row.names=1)
counts <- read.csv(paste0(direct,tag,'/exp.csv'),row.names=1, check.names=F, stringsAsFactors=FALSE)
#  

smpl=sample(rownames(colData),cell)
counts=counts[smpl,]
colData=colData[smpl,]
rownames(rowData)<-rowData$gene_id


rownames(counts)<- rownames(colData)
colnames(counts)<-rownames(rowData)
sp_locs<-colData[,c(2,3)]
row.names(colData)=colData$cell_id
#t(counts)
osm_test<-createGiottoObject(
  t(as.matrix(counts)),
  spatial_locs =sp_locs,
  cell_metadata = colData,
  gene_metadata = rowData
)

## filter
osm_test <- filterGiotto(gobject = osm_test,
                         expression_threshold = 1,
                         gene_det_in_min_cells = 10,
                         min_det_genes_per_cell = 10,
                         expression_values = c('raw'),
                         verbose = T)

## normalize
# 1. standard z-score way
osm_test <- normalizeGiotto(gobject = osm_test)
# 2. osmFISH way
raw_expr_matrix = osm_test@raw_exprs

norm_genes = (raw_expr_matrix/rowSums(raw_expr_matrix)) * nrow(raw_expr_matrix)
norm_genes_cells = t((t(norm_genes)/colSums(norm_genes)) * ncol(raw_expr_matrix))
osm_test@custom_expr = norm_genes_cells

## add gene & cell statistics
osm_test <- addStatistics(gobject = osm_test)
## add gene & cell statistics
osm_test <- addStatistics(gobject = osm_test)

# look at svg genes
osm_test <- createSpatialNetwork(gobject = osm_test)
kmtest = binSpect(osm_test, bin_method = 'kmeans')

kmtestsub<-kmtest[,c('adj.p.value')]
kmtestsub$SVG<-kmtest$adj.p.value<0.05
rownames(kmtestsub)<-kmtest$genes
kmtestsub$rank<-rank(kmtest$score)
colnames(kmtestsub)<-c('pvalue','SVG','rank')

write.csv(kmtestsub,paste0('D:/Xenium_benchmarking-main/figures/SVF/',tag,'__giotto_kmeans.csv'))




kmtest = binSpect(osm_test, bin_method = 'rank')

kmtestsub<-kmtest[,c('adj.p.value')]
kmtestsub$SVG<-kmtest$adj.p.value<0.05
rownames(kmtestsub)<-kmtest$genes
kmtestsub$rank<-rank(kmtest$score)
colnames(kmtestsub)<-c('pvalue','SVG','rank')

write.csv(kmtestsub,paste0('D:/Xenium_benchmarking-main/figures/SVF/',tag,'__giotto_rank.csv'))

}



tag<-"ms_brain_multisection1"
cells<-c(100,500,1000,5000,10000,50000)
times<-c()
#########################time calculation
for (cell in cells){
  print(tag)
  print(cell)
  start<-Sys.time()
  
  print(tag)
  rowData <- read.csv(paste0(direct,tag,'/var.csv'), stringsAsFactors=FALSE)
  colData <- read.csv(paste0(direct,tag,'/obs.csv'), stringsAsFactors=FALSE, row.names=1)
  counts <- read.csv(paste0(direct,tag,'/exp.csv'),row.names=1, check.names=F, stringsAsFactors=FALSE)
  #  
  
  smpl=sample(rownames(colData),cell)
  counts=counts[smpl,]
  colData=colData[smpl,]
  rownames(rowData)<-rowData$gene_id
  
  
  rownames(counts)<- rownames(colData)
  colnames(counts)<-rownames(rowData)
  sp_locs<-colData[,c(2,3)]
  row.names(colData)=colData$cell_id
  #t(counts)
  osm_test<-createGiottoObject(
    t(as.matrix(counts)),
    spatial_locs =sp_locs,
    cell_metadata = colData,
    gene_metadata = rowData
  )
  ## filter
  osm_test <- filterGiotto(gobject = osm_test,
                           expression_threshold = 1,
                           gene_det_in_min_cells = 10,
                           min_det_genes_per_cell = 10,
                           expression_values = c('raw'),
                           verbose = T)
  ## normalize
  # 1. standard z-score way
  osm_test <- normalizeGiotto(gobject = osm_test)
  # 2. osmFISH way
  raw_expr_matrix = osm_test@raw_exprs
  norm_genes = (raw_expr_matrix/rowSums(raw_expr_matrix)) * nrow(raw_expr_matrix)
  norm_genes_cells = t((t(norm_genes)/colSums(norm_genes)) * ncol(raw_expr_matrix))
  osm_test@custom_expr = norm_genes_cells
  
  ## add gene & cell statistics
  osm_test <- addStatistics(gobject = osm_test)
  ## add gene & cell statistics
  osm_test <- addStatistics(gobject = osm_test)
  
  # look at svg genes
  osm_test <- createSpatialNetwork(gobject = osm_test)
  kmtest = binSpect(osm_test, bin_method = 'kmeans')
  
  kmtestsub<-kmtest[,c('adj.p.value')]
  kmtestsub$SVG<-kmtest$adj.p.value<0.05
  rownames(kmtestsub)<-kmtest$genes
  kmtestsub$rank<-rank(kmtest$score)
  colnames(kmtestsub)<-c('pvalue','SVG','rank')
  
 # write.csv(kmtestsub,paste0('D:/Xenium_benchmarking-main/figures/SVF/',tag,'__giotto_kmeans.csv'))
  end<-Sys.time()
  
  time<-as.numeric(end-start)
  times<-append(times,time)
}




tag<-"ms_brain_multisection1"
cells<-c(100,500,1000,5000,10000,50000)
times2<-c()
#########################time calculation
for (cell in cells){
  print(tag)
  print(cell)
  start<-Sys.time()
  
  print(tag)
  rowData <- read.csv(paste0(direct,tag,'/var.csv'), stringsAsFactors=FALSE)
  colData <- read.csv(paste0(direct,tag,'/obs.csv'), stringsAsFactors=FALSE, row.names=1)
  counts <- read.csv(paste0(direct,tag,'/exp.csv'),row.names=1, check.names=F, stringsAsFactors=FALSE)
  #  
  
  smpl=sample(rownames(colData),cell)
  counts=counts[smpl,]
  colData=colData[smpl,]
  rownames(rowData)<-rowData$gene_id
  
  
  rownames(counts)<- rownames(colData)
  colnames(counts)<-rownames(rowData)
  sp_locs<-colData[,c(2,3)]
  row.names(colData)=colData$cell_id
  #t(counts)
  osm_test<-createGiottoObject(
    t(as.matrix(counts)),
    spatial_locs =sp_locs,
    cell_metadata = colData,
    gene_metadata = rowData
  )
  ## filter
  osm_test <- filterGiotto(gobject = osm_test,
                           expression_threshold = 1,
                           gene_det_in_min_cells = 10,
                           min_det_genes_per_cell = 10,
                           expression_values = c('raw'),
                           verbose = T)
  ## normalize
  # 1. standard z-score way
  osm_test <- normalizeGiotto(gobject = osm_test)
  # 2. osmFISH way
  raw_expr_matrix = osm_test@raw_exprs
  norm_genes = (raw_expr_matrix/rowSums(raw_expr_matrix)) * nrow(raw_expr_matrix)
  norm_genes_cells = t((t(norm_genes)/colSums(norm_genes)) * ncol(raw_expr_matrix))
  osm_test@custom_expr = norm_genes_cells
  
  ## add gene & cell statistics
  osm_test <- addStatistics(gobject = osm_test)
  ## add gene & cell statistics
  osm_test <- addStatistics(gobject = osm_test)
  
  # look at svg genes
  
  osm_test <- createSpatialNetwork(gobject = osm_test)
  kmtest = binSpect(osm_test, bin_method = 'rank')
  
  kmtestsub<-kmtest[,c('adj.p.value')]
  kmtestsub$SVG<-kmtest$adj.p.value<0.05
  rownames(kmtestsub)<-kmtest$genes
  kmtestsub$rank<-rank(kmtest$score)
  colnames(kmtestsub)<-c('pvalue','SVG','rank')
  
#  write.csv(kmtestsub,paste0('D:/Xenium_benchmarking-main/figures/SVF/',tag,'__giotto_kmeans.csv'))
  end<-Sys.time()
  
  time<-str(end-start)
  times2<-as.numeric(times2,time)
}





results<-data.frame(cells,times,times2)

colnames(results)<-c('cells','time_giotto_kmeans','time_giotto_rank')


write.csv(results,'D:/Xenium_benchmarking-main/figures/times_svf/giotto_times.csv')
















# write.csv(kmtestsub,paste0('D:/Xenium_benchmarking-main/figures/SVF/',tag,'__giotto_rank.csv'))



















## calculate frequently seen proximities
cell_proximities = cellProximityEnrichment(gobject = osm_test,
                                           'initial_annotation',
                                           number_of_simulations = 1000,spatial_network_name = "Delaunay_network")

spatial_network_name


## barplot
cellProximityBarplot(gobject = osm_test, CPscore = cell_proximities, min_orig_ints = 25, min_sim_ints = 25,
                     save_param = c(save_name = '12_a_barplot_cell_cell_enrichment_multisection'))
cellProximityHeatmap(gobject = osm_test, CPscore = cell_proximities, order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'),
                     save_param = c(save_name = '12_b_heatmap_cell_cell_enrichment_multisection', unit = 'in'))

## network
cellProximityNetwork(gobject = osm_test, CPscore = cell_proximities, only_show_enrichment_edges = T)#, remove_self_edges = T, only_show_enrichment_edges = T,
save_param = c(save_name = '12_c_network_cell_cell_enrichment_multisection')


# HMRF2

my_spatial_genes = kmtest[1:100]$genes

HMRF_spatial_genes = doHMRF(gobject = osm_test,k=9,spatial_genes=kmtest[1:200]$genes,spatial_dimensions = c("sdimx", "sdimy"),betas = c(20, 2, 50))


for(i in seq(20, 50, by = 2)) {
  viewHMRFresults2D(gobject = osm_test,
                    HMRFoutput = HMRF_spatial_genes,
                    k = 9, betas_to_view = i,
                    point_size = 2)
}

viewHMRFresults2D(gobject = osm_test,
                  HMRFoutput = HMRF_spatial_genes,
                  k = 9, betas_to_view = 20,
                  point_size = 2)




my_giotto_object = addHMRF(gobject = osm_test,
                           HMRFoutput = HMRF_spatial_genes,
                           k = 9, betas_to_add = c(28),
                           hmrf_name = 'HMRF')

# visualize selected hmrf result
giotto_colors = Giotto:::getDistinctColors(25)
names(giotto_colors) = 1:25
spatPlot(gobject = osm_test,
         point_size = 3, coord_fix_ratio = 1, cell_color_code = giotto_colors)








































##########################
my_giotto_object<-osm_test

# create network (required for binSpect methods)
my_giotto_object = createSpatialNetwork(gobject = my_giotto_object, minimum_k = 2)

# identify genes with a spatial coherent expression profile
km_spatialgenes = binSpect(my_giotto_object, bin_method = 'kmeans')



# create a directory to save your HMRF results to
hmrf_folder = paste0(getwd(),'/','11_HMRF/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)


# perform hmrf
my_spatial_genes = km_spatialgenes$genes
HMRF_spatial_genes = doHMRF(gobject = my_giotto_object,
                            expression_values = 'scaled',
                            spatial_genes = my_spatial_genes,
                            k = 9,
                            betas = c(10,10,10),
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_top100_k9_scaled_multisection'))

# check and visualize hmrf results
for(i in seq(10, 10, by = 10)) {
  viewHMRFresults2D(gobject = my_giotto_object,
                    HMRFoutput = HMRF_spatial_genes,
                    k = 9, betas_to_view = i,
                    point_size = 2)
}

my_giotto_object = addHMRF(gobject = my_giotto_object,
                           HMRFoutput = HMRF_spatial_genes,
                           k = 9, betas_to_add = c(10),
                           hmrf_name = 'HMRF')

# visualize selected hmrf result
giotto_colors = Giotto:::getDistinctColors(9)
names(giotto_colors) = 1:9
spatPlot(gobject = my_giotto_object, cell_color = 'HMRF_k9_b.28',
         point_size = 3, coord_fix_ratio = 1, cell_color_code = giotto_colors)
