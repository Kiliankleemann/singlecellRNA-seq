if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

list.of.packages <- c('tidyverse', 'Seurat', 'Signac', 'enrichR', 'openxlsx', 'patchwork', 'data.table', 'dittoSeq', 'ggplot2','scCATCH', 'clustree')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

options(future.globals.maxSize = 4000 * 1024^5)

##################################  DATA PREPARATION  ############################################
h5_directory_list <- list.files(".//cibersorting/single_cell/GSE188288_RAW/", full.names = T, pattern=NULL, all.files=FALSE)
h5_directory_list

h5_name_list <- gsub(".*//", "", h5_directory_list)
h5_name_list <- gsub(".h5", "", h5_name_list)

h5_name_list

for (i in h5_directory_list){
  print(i)
}


for (k in h5_name_list) {
    assign(paste("counts_",k, sep=""),Read10X_h5(paste0('cibersorting/single_cell/GSE188288_RAW/',k,'.h5'), use.names = TRUE, unique.features = TRUE))
}



sobj1 <- CreateSeuratObject(counts = counts_GSM5676978_Neut_F1_sc_50K, project = 'Neut_F1',min.cells = 3, min.features=200)
sobj2 <- CreateSeuratObject(counts = counts_GSM5676979_Neut_F2_sc_50K, project = 'Neut_F2',min.cells = 3, min.features=200)
sobj3 <- CreateSeuratObject(counts = counts_GSM5676980_Neut_F3_sc_50K, project = 'Neut_F3',min.cells = 3, min.features=200)
sobj4 <- CreateSeuratObject(counts = counts_GSM5676981_Neut_M0_sc_50K, project = 'Neut_M0',min.cells = 3, min.features=200)
sobj5 <- CreateSeuratObject(counts = counts_GSM5676982_Neut_M1_sc_25K, project = 'Neut_M1',min.cells = 3, min.features=200)
sobj6 <- CreateSeuratObject(counts = counts_GSM5676983_Neut_M2_sc_50K, project = 'Neut_M2',min.cells = 3, min.features=200)
sobj7 <- CreateSeuratObject(counts = counts_GSM5676984_Neut_M3_sc_50K, project = 'Neut_M3',min.cells = 3, min.features=200)



# noquote(paste0("counts_",k))
# paste0("counts_",k)
# paste0("project_",k)
# 
# for (k in h5_name_list) {
#   print(k)
#   object_name <- noquote(paste0("counts_",k))
#  assign(paste0("srt_obj_",k),CreateSeuratObject(counts = object_name, project = paste0("project_",k),min.cells = 3, min.features=200))
# }
# 
# # object_name <- noquote(paste0("counts_",k))
# # object_name
# # call(noquote(paste0("counts_",k)))

#Filter out dublets

# Idents(sobj13) = "db.scDblFinder_class"
# sob13_fil <- subset(sobj13, idents = 'singlet')



############################  MERGING ###########################################
projectname <-'Human_neutrophils_'
scrna_merged <- merge(sobj1, y = c(sobj2,
                                   sobj3,
                                   sobj4,
                                   sobj5,
                                   sobj6,
                                   sobj7), add.cell.ids = c('Neut_F1','Neut_F2','Neut_F3','Neut_M0','Neut_M1','Neut_M2','Neut_M3'), project = projectname)




setable(scrna_merged$orig.ident)
table(Idents(scrna_merged))

#raw counts for merged files 

all_counts <- scrna_merged@Dimnames[[1]]
write.xlsx(alldata_counts@Dimnames[[1]], file="Allgenes")

# dataset_Clec7P <- subset(x = scrna_merged, subset = Clec7a > 0.00011)
# write.csv(scrna_merged@assays[["RNA"]]@counts@Dimnames[[1]], file="allgenes")
# genes <- as.data.frame(data1@Dimnames[[1]])


######################### MITOCHONDRIAL REGRESSION ############################################
# Store mitochondrial gene statistics in your Seurat object
scrna_merged[['percent_mt']] <- PercentageFeatureSet(scrna_merged, pattern = '^MT-')

# Basic QC plot to set cutoffs
VlnPlot_nFeature_nCount_mt_unfiltered <- VlnPlot(scrna_merged, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), ncol = 3)
Scatter_nCount_mt_unfiltered <- FeatureScatter(scrna_merged, feature1 = "nCount_RNA", feature2 = "percent_mt")
Scatter_nFeature_mt_unfiltered <- FeatureScatter(scrna_merged, feature1 = "nFeature_RNA", feature2 = "percent_mt")
Scatter_nCount_nFeature_unfiltered <- FeatureScatter(scrna_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

dev.off()
pdf(file = paste0('QC/Scatter_nCount_mt_unfiltered.pdf'), pointsize = 10)
Scatter_nCount_mt_unfiltered
dev.off()

VlnPlot(scrna_merged, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), ncol = 3)

# Filter data to remove unwanted cells
scrna_final <- subset(scrna_merged, subset = nFeature_RNA > 100 & nFeature_RNA < 2500  & percent_mt < 10)

Idents(scrna_final) <- 'orig.ident'
scrna_final$orig.ident<- factor(x = scrna_final$orig.ident, levels = c('Neut_F1','Neut_F2','Neut_F3','Neut_M0','Neut_M1','Neut_M2','Neut_M3'))
sample_counts <- table(scrna_final$orig.ident)
write.xlsx(sample_counts, file = paste0('results/cibersort/Filtered_Sample_counts.xlsx'), overwrite = T)

dev.off()
pdf(file = paste0('plots/cibersort/VlnPlot_nFeature_nCount_mt_filtered'), pointsize = 10)
VlnPlot(scrna_final, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), ncol = 3)
dev.off()

pdf(file = paste0('QC/Scatter_nFeature_nCount_mt_filtered'), pointsize = 10)
scat1 <- FeatureScatter(scrna_final, feature1 = 'nCount_RNA', feature2 = 'percent_mt')
scat2 <- FeatureScatter(scrna_final, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
scat1 + scat2
dev.off()

########################### NORMALIZATION AND SCALING #######################################
# Normalize RNA expression data and scale to variable features
scrna_final <- NormalizeData(scrna_final, normalization.method = "LogNormalize", scale.factor = 10000)
scrna_final <- FindVariableFeatures(scrna_final, selection.method = 'vst', nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scrna_final), 10)
top10
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scrna_final)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

scrna_final <- ScaleData(scrna_final, features = VariableFeatures(scrna_final), vars.to.regress = 'percent_mt')

#################################### DIMENSIONAL REDUCTIONS##############################
# Run principal component analysis
scrna_final <- RunPCA(scrna_final, features = VariableFeatures(object = scrna_final))
scrna_final <- JackStraw(object = scrna_final, num.replicate = 50, prop.freq=0.025, dims = 50)
scrna_final <- ScoreJackStraw(scrna_final, dims = 1:50)
JackStrawPlot <- JackStrawPlot(object = scrna_final, dims = 1:50, xmax = 0.05) + guides(col = guide_legend(ncol = 1)) + theme(legend.text = element_text(size = 6), legend.key.size = unit(0.02, "cm"))

# Visualize PCA to ensure merged samples are comparable
dev.off()
Idents(scrna_final) <- 'orig.ident'
pdf(file = paste0('DimPlot_unprocessed.pdf'), pointsize = 10)
DimPlot(scrna_final, reduction = 'pca')
dev.off()

# Visualize component strengths to decide how many to use
ElbowPlot <- ElbowPlot(object = scrna_final, ndims = 30)
pdf(file = paste0('QC/ElbowPlot.pdf'), pointsize = 10, width = 16)
ElbowPlot
dev.off()


# Determine percent of variation associated with each PC
pct <- scrna_final[["pca"]]@stdev / sum(scrna_final[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pcs

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
elbow_rank_vs_pcs <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

pdf(file = paste0('QC/elbow_rank_vs_pcs' ,projectname,'.pdf'), pointsize = 10, width = 16)
elbow_rank_vs_pcs
dev.off()

print(x = scrna_final[["pca"]], 
      dims = 1:25, 
      nfeatures = 5)


#Heatmap exploration
pdf(file = paste0('QC/heatmap_PCreductions_', projectname, '.pdf'), width = 20, height = 30)
DimHeatmap(scrna_final, dims = 1:20, cells = 1000, balanced = TRUE)
dev.off()

# Cluster cells according to elbow plot dimension choice
scrna_final <- FindNeighbors(scrna_final, dims = 1:10)

#Investigate resolution
scrna_final <- FindClusters(scrna_final, resolution = 1.2)
scrna_final <- FindClusters(scrna_final, resolution = 1.1)
scrna_final <- FindClusters(scrna_final, resolution = 1)
scrna_final <- FindClusters(scrna_final, resolution = 0.9)
scrna_final <- FindClusters(scrna_final, resolution = 0.8)
scrna_final <- FindClusters(scrna_final, resolution = 0.7)
scrna_final <- FindClusters(scrna_final, resolution = 0.6)
scrna_final <- FindClusters(scrna_final, resolution = 0.5)
scrna_final <- FindClusters(scrna_final, resolution = 0.4)
scrna_final <- FindClusters(scrna_final, resolution = 0.3)
scrna_final <- FindClusters(scrna_final, resolution = 0.2)

#Clustree exploration
cluster_tree <- clustree(scrna_final)
pdf(file = paste0('QC/cluster_tree_', projectname, '.pdf'), width = 12, height = 9)
cluster_tree
dev.off()

#Clusterin metrics
projectname <- 'scrna_merged_'
metrics <- 'dim13_res0.9'

# Cluster cells according to elbow plot dimension choice
scrna_final <- FindNeighbors(scrna_final, dims = 1:13)
scrna_final <- FindClusters(scrna_final, resolution = 0.9)

# Run UMAP reduction to visualize clusters
scrna_final <- RunUMAP(scrna_final, dims = 1:10)

# Plot UMAP
Idents(scrna_final) <- 'seurat_clusters'

pdf(file = paste0('plots/seurat_clusters_',projectname, metrics, '.pdf'), width = 12, height = 9)
DimPlot(scrna_final, label = T) 
dev.off()

pdf(file = paste0('plots/seurat_clusters_split_origident_',projectname, metrics, '.pdf'), width = 30, height = 20)
DimPlot(scrna_final, label = T, split.by = 'orig.ident', ncol = 4) 
dev.off()


Idents(scrna_final) <- 'orig.ident'
levels(scrna_final)
new_mappings <- c('1','1','1','1','2','2')
new_mappings <- c('APP/PS1-miR155cKO','APP/PS1-miR155cKO', 'APP/PS1-miR155cKO',  'APP/PS1', 'APP/PS1','APP/PS1')
new_mappings <- c('m','m', 'f',  'f', 'f','m', 'm','f', 'f','f', 'm','m')

names(new_mappings) <- levels(scrna_final)
scrna_final <- RenameIdents(scrna_final, new_mappings)
scrna_final$Batch <- Idents(scrna_final)
scrna_final$Genotype <- Idents(scrna_final)
scrna_final$Sex <- Idents(scrna_final)

Idents(scrna_final) <- 'Batch'
dataset_new_debatch <- RunHarmony(scrna_final, "Batch", plot_convergence = TRUE)
dataset_new_debatch <- RunUMAP(dataset_new_debatch, reduction = "harmony", dims = 1:20)


pdf(file = paste0('plots/DEBATCH_group_batch_',projectname, metrics, '.pdf'), width = 10, height = 10)
DimPlot(dataset_new_debatch, group.by = 'Batch')
dev.off()


#INF
INF_module = list(c('Ifit3', 'Isg15', 'Ifit2', 'Irf7', 'Ifit1', 'Usp18', 'Stat1','Fcgr1'))
INF_module <- AddModuleScore(dataset_new_debatch,
                             features = INF_module,
                             name="INF_module")
INF_module_featureplot <- FeaturePlot(INF_module,
                                      features = "INF_module1", split.by = 'Genotype') 
dev.off()
pdf(file = paste0('plots/umap_Featureplot_Interferon_Module_', projectname,metrics, '.pdf'), width = 12, height = 9)
INF_module_featureplot
dev.off()


#Clusterin metrics
projectname <- 'male_debatched_'
metrics <- 'dim10_res0.9'
# Save final object
saveRDS(dataset_new_debatch, file = paste0('seurat_objects/',projectname, metrics, '.rds'))
