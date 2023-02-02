
library(Giotto)
library(smfishHmrf)


#  
rowData <- read.csv("D:/Xenium_benchmarking-main/Xenium_benchmarking-main/data/formatted_for_R/var_msbrain_1.csv", stringsAsFactors=FALSE)
colData <- read.csv("D:/Xenium_benchmarking-main/Xenium_benchmarking-main/data/formatted_for_R/obs_msbrain_1.csv", stringsAsFactors=FALSE, row.names=1)
counts <- read.csv("D:/Xenium_benchmarking-main/Xenium_benchmarking-main/data/formatted_for_R/exp_msbrain_1.csv",row.names=1, check.names=F, stringsAsFactors=FALSE)
#  
rownames(counts)<- rownames(colData)
colnames(counts)<-rownames(rowData)

sp_locs<-colData[,c(3,4)]

t(counts)
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
norm_genes = (raw_expr_matrix/rowSums_giotto(raw_expr_matrix)) * nrow(raw_expr_matrix)
norm_genes_cells = t_giotto((t_giotto(norm_genes)/colSums_giotto(norm_genes)) * ncol(raw_expr_matrix))
osm_test@custom_expr = norm_genes_cells

## add gene & cell statistics
osm_test <- addStatistics(gobject = osm_test)

## add gene & cell statistics
osm_test <- addStatistics(gobject = osm_test)

# save according to giotto instructions
spatPlot(gobject = osm_test, cell_color = 'initial_annotation', point_size = 1.5,
         save_param = list(save_name = '2_a_original_clusters'))





showClusterHeatmap(gobject = osm_test, expression_values = 'custom', cluster_column = 'initial_annotation',
                   save_param = list(save_name = '4_e_heatmap', units = 'cm'),
                   row_names_gp = grid::gpar(fontsize = 6), column_names_gp = grid::gpar(fontsize = 6))



osm_test <- createSpatialNetwork(gobject = osm_test)

# look at svg genes

kmtest = binSpect(osm_test, bin_method = 'kmeans')

## calculate frequently seen proximities

cell_proximities = cellProximityEnrichment(gobject = osm_test,
                                           cluster_column = 'initial_annotation',
                                           number_of_simulations = 1000)
## barplot
cellProximityBarplot(gobject = osm_test, CPscore = cell_proximities, min_orig_ints = 25, min_sim_ints = 25,
                     save_param = c(save_name = '12_a_barplot_cell_cell_enrichment'))
cellProximityHeatmap(gobject = osm_test, CPscore = cell_proximities, order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'),
                     save_param = c(save_name = '12_b_heatmap_cell_cell_enrichment', unit = 'in'))

## network
cellProximityNetwork(gobject = osm_test, CPscore = cell_proximities, only_show_enrichment_edges = T)#, remove_self_edges = T, only_show_enrichment_edges = T,
                     save_param = c(save_name = '12_c_network_cell_cell_enrichment')


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
my_spatial_genes = km_spatialgenes[1:100]$genes
HMRF_spatial_genes = doHMRF(gobject = my_giotto_object,
                            expression_values = 'scaled',
                            spatial_genes = my_spatial_genes,
                            k = 9,
                            betas = c(10,10,10),
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_top100_k9_scaled'))

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
