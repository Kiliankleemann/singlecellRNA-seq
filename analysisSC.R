################## -------- INSTALLING AND LOADING PACKAGES ------- ####
list.of.packages <- c("BiocGenerics","tximport","S4Vectors", "DESeq2", "biomaRt",'GseaVis',"OmnipathR",
                      "ggplot2", "ggsignif", "ggpubr", "sva", "devtools", "org.Hs.eg.db", "DoubletFinder",
                      "org.Mm.eg.db", "limma","stringr","KEGGREST","ggrepel", "openxlsx","SingleCellExperiment",
                      "fgsea","clusterProfiler","pheatmap","ggpubr","cowplot",'dplyr',"clustree","decoupleR",
                      "RColorBrewer",'AnnotationDbi', 'tidyverse','pheatmap', 'dendextend','scCustomize',
                      "data.table",  "wesanderson", "AUCell", "GSEABase", "GSVA","HGNChelper","Nebulosa",
                       "SCP", "dittoSeq", "limma", "matrixStats", "sctransform", "scDblFinder", "Seurat")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if(length(new.packages)) BiocManager::install(new.packages)

# Special packages
# devtools::install_github("junjunlab/GseaVis")
# devtools::install_github("zhanghao-njmu/SCP")
# SCP::PrepareEnv()
# devtools::install_github("YosefLab/VISION") 
# devtools::install_github("wu-yc/scMetabolism")
# devtools::install_github('immunogenomics/presto')
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
# devtools::install_github('KChen-lab/METAFlux')

library(SCP)
library(Seurat)
library(DoubletFinder)
library(SingleCellExperiment)
library(scMetabolism)
library(METAFlux)

# Packages loading
invisible(lapply(list.of.packages, library, character.only = TRUE))

#Set computing
options(future.globals.maxSize = 1e9)


######################## LOADING DATA ##########################################
setwd("/Users/")

# Load Seurat data set
projectname = 'projectname'
dataset <- readRDS(paste0('seurat_objects/SeuratObject_', projectname, '.rds'))

######################## QC (Ribosomal, mitochondiral cell cycle) ###############################
data <- CellCycleScoring(dataset, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)

# Basic function to convert human to mouse gene names
human_genes_to_mouse <- function(human_genes, mirror = "www") {
  # Define the mart objects for human and mouse
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
                   host = "https://dec2021.archive.ensembl.org/")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", 
                   host = "https://dec2021.archive.ensembl.org/")
  
  # Retrieve the corresponding mouse genes
  genes.retrieved <- getLDS(attributes = c("hgnc_symbol"), 
                            filters = "hgnc_symbol", 
                            values = human_genes, 
                            mart = human, 
                            attributesL = c("mgi_symbol"), 
                            martL = mouse, 
                            uniqueRows = TRUE)
  
  # Return the list of mouse genes
  mouse_genes <- genes.retrieved$MGI.symbol
  
  return(mouse_genes)
}
m.s.genes <- human_genes_to_mouse(cc.genes.updated.2019$s.genes)
m.g2m.genes <- human_genes_to_mouse(cc.genes.updated.2019$g2m.genes)

dataset <- CellCycleScoring(dataset, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
Cell_cycle_per_cluster <- data@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  dplyr::count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster") +
  scale_fill_manual(values=c('navy', 'cyan3', 'firebrick')) +
  theme_minimal()

dir.create(paste0('plots/CellCyle'))
pdf(file = paste0('plots/CellCyle/cell_cycle_per_cluster', projectname,processing_specs, '.pdf'), width = 18, height = 9)
Cell_cycle_per_cluster
dev.off()

#Umap viz
Cell_cycle_cluster_umap <- CellDimPlot(srt = dataset, group.by = "seurat_clusters", stat.by = "Phase",
  reduction = "UMAP", theme_use = "theme_blank")

pdf(file = paste0('plots/CellCyle/UMAP_cell_cycle_per_cluster_', projectname,processing_specs, '.pdf'), width = 12, height = 9)
Cell_cycle_cluster_umap
dev.off()



######################## METADATA MANIPULATION AND CELLTYPE ASSIGNMENT ##########################
#Seurat Clusters
Idents(dataset) <- 'seurat_clusters'
levels(dataset)
new_mappings <- c("00" , "01" , "02" , "03",  "04",  "05" , "06",  "07",  "08" , "09" , "10",
                  "11", "12", "13" ,"14" ,"15" ,"16" ,"17", "18", "19", "20", "21",
                  "22" ,"23", "24")
names(new_mappings) <- levels(dataset)
dataset <- RenameIdents(dataset, new_mappings)
dataset$seurat_clusters_2<- Idents(dataset)

#Assign Sex metadata
Idents(dataset) <- 'orig.ident'
levels(dataset)
new_mappings <- c("F","M","F","M", "F", "M", "F", "M")
names(new_mappings) <- levels(dataset)
dataset <- RenameIdents(dataset, new_mappings)
dataset$Sex<- Idents(dataset)

#Assign Disease metadata
Idents(dataset) <- 'orig.ident'
levels(dataset)
new_mappings <- c("WT","WT","WT","WT", "AD", "AD", "AD", "AD")
names(new_mappings) <- levels(dataset)
dataset <- RenameIdents(dataset, new_mappings)
dataset$Disease<- Idents(dataset)


#Assign genotype metadata
Idents(dataset) <- 'orig.ident'
levels(dataset)
new_mappings <- c("WT","WT","KO","KO", "WT", "WT", "KO", "KO")
names(new_mappings) <- levels(dataset)
dataset <- RenameIdents(dataset, new_mappings)
dataset$GLS_Genotype<- Idents(dataset)


######## AUTOMATED CELL TYPE ANNOTATION #######
# load auto-detection function
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
#Expadite 
source("https://raw.githubusercontent.com/kris-nader/sc-type/master/R/sctype_wrapper.R"); 
dataset <- run_sctype(dataset,known_tissue_type="Immune system",custom_marker_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",name="sctype_classification",plot=TRUE)
dev.off()

dir.create(paste0('plots/CellAnnotation'))
pdf(file = paste0('plots/CellAnnotation/SCType_Immunesystem_', projectname,processing_specs, '.pdf'), width = 12, height = 9)
DimPlot(dataset, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')        
dev.off()

######## MARKER GENES #######
Idents(dataset) = 'seurat_clusters'
dataset_filtered_seurat_clusters_markers <- FindAllMarkers(dataset, only.pos = T, min.pct = 0.4, logfc.threshold = 0.5)

dir.create(paste0("results/markergenes/"))
write.xlsx(dataset_filtered_seurat_clusters_markers, paste0('results/markergenes/dataset_filtered_seurat_clusters_markers',processing_specs, '.xlsx'))


######## COLORS ##########
#Color Scheme for seurat clusters
colors_seurat_clusters <-c( "#FF0000", "#8B0000", "#FFA500", "#FF8C00", "#FFFF00", "#FFD700",
"#008000", "#32CD32", "#008080", "#00FFFF", "#87CEEB", "#0000FF",
"#00008B",  "#4B0082", "#EE82EE", "#FF00FF", "#800080",
"#FFC0CB", "#FF69B4", "#A52A2A", "#800000", "#F5F5DC", "#808000",
 "#C0C0C0", "#808080")

#Color Scheme for genotypes
colors_genotypes <-c( "#87CEEB", "#4B0082", "#FF8C00", "#8B0000")
colors_genotypes_WT <-c("#87CEEB", "#4B0082")
colors_genotypes_AD <-c("#FF8C00", "#8B0000")

#Color Scheme for Celltypes clusters
colors_CellType <- c("#FF0000", "#8B0000", "#FFA500", "#FF8C00", "#FFD700","#008000","#32CD32", 
                   "#008080", "#00FFFF", "#87CEEB", "#0000FF","#00008B", "#800080","#FF00FF",
                   "#FFC0CB", "#FF69B4","#808000",
                   "#808080")



######## MARKER GENE VIZUALIZATION ########
dir.create(paste0('plots/Marker_genes'))

#### Microglia Genes
Microglia_markers <- c("Tmem119", "P2ry12","Sall1", "Sparc","Crybb1","Gpr34")
dir.create(paste0('plots/Marker_genes/',"Microglia"))

for (i in Microglia_markers) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"Microglia","/FeaturePlot_", i,".pdf"),fplot)
}

pdf(file = paste0('plots/Marker_genes/Microglia/',"VlnPlot_",'microglia_markers' , '.pdf'), width = 10, height = 5)
VlnPlot(dataset, features = Microglia_markers)
dev.off()


#### MGnD genes
Microglia_markers <- c("Cst7", "Lpl","Axl","Clec7a")
dir.create(paste0('plots/Marker_genes/',"MGnD"))
for (i in Microglia_markers) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"MGnD","/FeaturePlot_", i, '.pdf'),fplot )
}

pdf(file = paste0('plots/Marker_genes/MGnD/',"VlnPlot_",Microglia_markers , '.pdf'), width = 10, height = 5)
VlnPlot(dataset, features = Microglia_markers)
dev.off()


##### Interferon responsive Microglia 
Microglia_markers <- c("Ifit3","Clec2d","Usp18","Isg15","Slfn5")
dir.create(paste0('plots/Marker_genes/',"IRM"))
for (i in Microglia_markers) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"IRM","/FeaturePlot_", i, '.pdf'),fplot )
}


##### BAM / CAM Genes
BAM_CAM <- c("Mrc1","Ms4a7","Dab2","Pf4", "Clec12a")
dir.create(paste0('plots/Marker_genes/',"BAM_CAM"))
for (i in BAM_CAM) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"BAM_CAM","/FeaturePlot_", i, '.pdf'),fplot )
}


##### MHC2+ BAM / CAM Genes
BAM_CAM <- c("H2-Aa", "Slamf7", "Stap1", "Ciita", "Cd72","H2-Eb1","H2-Ab1")
dir.create(paste0('plots/Marker_genes/',"BAM_CAM_MHC2"))
for (i in BAM_CAM) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"BAM_CAM_MHC2","/FeaturePlot_", i, '.pdf'),fplot )
}


##### DCs 
cDC_nonClassicalMono <- c("Cd209a","H2-DMb2","P2ry10","Tnip3")
dir.create(paste0('plots/Marker_genes/',"cDC_nonClassicalMono"))

for (i in cDC_nonClassicalMono) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"cDC_nonClassicalMono","/FeaturePlot_", i, '.pdf'),fplot )
}


##### Monocytes
cDC_nonClassicalMono <- c("Plac8", "Napsa" ,"Cytip", "S100a4", "Ifitm6", "Itgal", "Ace")
dir.create(paste0('plots/Marker_genes/',"Monocytes"))

for (i in cDC_nonClassicalMono) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"Monocytes","/FeaturePlot_", i, '.pdf'),fplot )
}



##### Neutrophil genes
genes <- c("Camp","Lcn2","Sell","Hdc", "Hp","Il1r2","Clec4d","Trem1", "Mmp9")
dir.create(paste0('plots/Marker_genes/',"Neutrophils"))
for (i in genes) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"Neutrophils","/FeaturePlot_", i, '.pdf'),fplot )
}


##### Bcell genes
genes <- c("Cd19", "Cd79a", "Cd79b")
dir.create(paste0('plots/Marker_genes/',"Bcells"))
for (i in genes) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"Bcells","/FeaturePlot_", i, '.pdf'),fplot )
}


###### NK-Cells genes
genes <- c("Nkg7", "Gzma", "Gzmb", "Ifng","Il2rb","Prf1")
dir.create(paste0('plots/Marker_genes/',"CD8_TCells"))
for (i in genes) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"CD8_TCells","/FeaturePlot_", i, '.pdf'),fplot )
}

##### Endothelial cell genes
genes <- c("Flt1","Cldn5","Ptprb" )
dir.create(paste0('plots/Marker_genes/',"Endothelial_Cells"))
for (i in genes) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"Endothelial_Cells","/FeaturePlot_", i, '.pdf'),fplot )
}

##### Ependymal cells
genes <- c("Cald1","Myl9","Igfbp7","Tpm2", "Ttr", "Pdgfrb", "Acta2")
dir.create(paste0('plots/Marker_genes/',"Ependymal_Cells"))
for (i in genes) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"Ependymal_Cells","/FeaturePlot_", i, '.pdf'),fplot )
}

##### Smooth muscle cell
genes <- c("Myl9", "Mylk")
dir.create(paste0('plots/Marker_genes/',"Vasc_smooth_musc_cells"))
for (i in genes) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"Vasc_smooth_musc_cells","/FeaturePlot_", i, '.pdf'),fplot )
}


##### Mesenchymal Stem Cells (MSCs)
genes <- c("Islr", "Efemp1")
dir.create(paste0('plots/Marker_genes/',"Mesenchymal_Stem_Cell"))
for (i in genes) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"Mesenchymal_Stem_Cell","/FeaturePlot_", i, '.pdf'),fplot )
}

##### Astrocytes
genes <- c("Slc1a2","Gfap","S100b", "Enpp2")
dir.create(paste0('plots/Marker_genes/',"Astrocytes"))
for (i in genes) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"Astrocytes","/FeaturePlot_", i, '.pdf'),fplot )
}

##### Oligodendrocytes
genes <- c("Cldn11","Mog","Mobp")
dir.create(paste0('plots/Marker_genes/',"Oligodendrocytes"))
for (i in genes) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"Oligodendrocytes","/FeaturePlot_", i, '.pdf'),fplot )
}

##### Cycling cells
genes <- c("Stmn1","Birc5","Top2a","Mki67")
dir.create(paste0('plots/Marker_genes/',"Cycling_cells"))
for (i in genes) {
  fplot <- FeaturePlot(dataset,features = i, order = T)
  ggsave(paste0('plots/Marker_genes/',"Cycling_cells","/FeaturePlot_", i, '.pdf'),fplot )
}




######## ASSIGN FINAL CELL TYPES #######
CellTypes <- read.xlsx(paste0("results/markergenes/dataset_filtered_seurat_clusters_markers",".xlsx")) %>% pull("Celltype") %>% na.omit()
CellTypes

#Major CellTypes
Idents(dataset) <- 'seurat_clusters'
levels(dataset)
new_mappings <- c(CellTypes)
names(new_mappings) <- levels(dataset)
dataset <- RenameIdents(dataset, new_mappings)
dataset$CellType<- Idents(dataset)

levels(dataset)
Idents(dataset) <- 'CellType' 
dir.create(paste0('plots/CellAnnotation/'))
pdf(file = paste0('plots/CellAnnotation/', projectname, '_CellTypes.pdf'), width = 10, height = 9)
DimPlot(dataset, label = T, raster=F, cols = colors_CellType, pt.size = 0.3)
dev.off()

#Sub CellTypes
Idents(dataset) <- 'seurat_clusters'
levels(dataset)
new_mappings <- c(CellTypes)
names(new_mappings) <- levels(dataset)
dataset <- RenameIdents(dataset, new_mappings)
dataset$CellType<- Idents(dataset)

levels(dataset)
Idents(dataset) <- 'CellType' 
pdf(file = paste0('plots/CellAnnotation/', projectname, '_CellTypes.pdf'), width = 10, height = 9)
DimPlot(dataset, label = T, raster=F, cols = cols_CellType, pt.size = 0.3)
dev.off()

saveRDS(dataset, file = paste0('seurat_objects/SeuratObject_', projectname,"_celltype_annotated", '.rds'))


######## SUBSETTING FOR CONDITIONS AND CELL TYPES #######
#Subset WT data
Idents(dataset) <- 'Disease' 

#Subset AD data
Idents(dataset) <- 'Treatment' 
dataset_WTF_5xFAD <- subset(dataset, idents= c('Ctrl', 'AD'))

saveRDS(dataset_WT_5xFAD, file = paste0('seurat_objects/SeuratObject_dataset_WT_5xFAD_', '.rds'))

#Filter out Microglia
Idents(dataset) <- 'CellType' 
levels(dataset)
dataset_microglia <- subset(dataset, idents= c('Microglia','IRM','MGnD'))

#Filter out Macrophages / DCs Bcells
dataset_macrophages<- subset(dataset, idents= c('Macrophage', "MHC2_Macrophages", "non_classical_DCs"))
dataset_neutrophils <- subset(dataset, idents= c('Neutrophils'))

######################## DIFFERENTIAL ANALYSIS #######################################
######## OVERALL DEGS #########
levels(dataset$Treatment)
var1 <- "xx"
var2 <- "yy"
group_comparison <- "Treatment"

dir.create(paste0("results/","DEGs_",group_comparison))

#Set ident
markers <- FindMarkers(dataset, ident.1 = var1, ident.2 = var2, group.by = group_comparison, min.pct = 0.05, logfc.threshold = 0.3, 
                                                     min.diff.pct = 0.0) %>% rownames_to_column('gene')

write.xlsx(markers, paste0("results/","DEGs_",group_comparison,'/', var1 ,'_vs_',var2, '.xlsx'))


######## DEGs PER CELLTYPE #########
Idents(dataset) <- 'CellType'
levels(dataset)

dir.create("results/CellType_DEG/")

var1 <- "xx"
var2 <- "yy"
group_comparison <- "Treatment"
comparison_name <- "xx_vs_yy"

for (i in CellTypes) {
  dir.create(paste0("results/CellType_DEG/", i))
  Celltype_DEGs <- FindMarkers(dataset, ident.1 = var1, ident.2 = var2, group.by = group_comparison, 
    subset.ident = i, only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.4) %>% 
    rownames_to_column('gene')
  write.xlsx(Celltype_DEGs, paste0("results/CellType_DEG/", i,"/",comparison_name, ".xlsx"))
  assign(paste0(i,"_", comparison_name), Celltype_DEGs)
}

#Merge all results
rm(df_total)
df_total <- data.frame()
for (i in CellTypes) {
  DEGs <- get(paste0(i,'_',comparison_name))
  DEGs_padj <- DEGs %>% filter(p_val_adj < 0.05)
  count <- nrow(DEGs_padj)
  df <- data.frame(count)
  df_total <- rbind(df_total, df)
}

df_total <- cbind(df_total,CellTypes)
df_total <- df_total %>% distinct(CellTypes, .keep_all = TRUE)
df_total$CellTypes <- factor(df_total$CellTypes, levels=c('Microglia','Macrophage','MGnD','Monocytes',
                                                          "MHC2_macrophages",'Neutrophils','Endothelial_cells',
                                                          'NK_cells','cDCs','IRM','Cycling_cells',
                                                          "Ependymal_cells", "unkown"))
levels(dataset$CellType)
pointSize = 10 
lineWidth = 0.5

degs_bar <- ggplot(df_total,aes(x = reorder(CellTypes, count), y = count, fill = CellTypes))+
  geom_col(width = 0.7, color = 'black') + scale_fill_manual(values = colors_CellType)+
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black", linewidth = 0.2),
        plot.title  = element_text(color="black", size=10),
        axis.title  = element_text(size = pointSize, colour = "black"),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
        axis.text.y  = element_text(size = pointSize , colour = "black"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.ticks = element_line(linewidth = 0.2, colour = "black"),
        axis.line = element_line(linewidth = 0.2, colour = "black"))+
  labs(y= 'Count of DEGs (Log2FC > 0.2)', x=NULL) +
  ggtitle(paste0('Number of DEGs per CellType',comparison_name,'.xlsx')) +
  scale_y_continuous( expand = c(0,0.1,0,0)) 
degs_bar

ggsave(
  paste0('plots/barplots/','DEGs_per_CellType_',comparison_name,'.pdf'),
  plot = degs_bar,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 150,
  height = 100,
  units = c("mm"),
  dpi = 600)


######## DEGs PER SEURAT CLUS #########
Idents(dataset) <- 'seurat_clusters_2'
levels(dataset)
seurat_clusters_DEGs <- levels(dataset)

dir.create("results/Seurat_clusters_DEG/")

comparison_group <- "Genotype"
ident_1 <- "xx"
ident_2 <- "yy"
comparison_name <- "xx_vs_yy"

for (i in seurat_clusters_DEGs) {
  dir.create(paste0("results/Seurat_clusters_DEG/", i))
  Celltype_DEGs <- FindMarkers(dataset, ident.1 = ident_1, ident.2 = ident_2, group.by = comparison_group, 
                               subset.ident = i, only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.4) %>% 
    rownames_to_column('gene')
  write.xlsx(Celltype_DEGs, paste0("results/Seurat_clusters_DEG/", i,"/",comparison_name, ".xlsx"))
  assign(paste0(i,"_", comparison_name), Celltype_DEGs)
}

#Merge all results
rm(df_total)
df_total <- data.frame()
for (i in seurat_clusters_DEGs) {
  DEGs <- get(paste0(i,'_',comparison_name))
  DEGs_padj <- DEGs %>% filter(p_val_adj < 0.05)
  count <- nrow(DEGs_padj)
  df <- data.frame(count)
  df_total <- rbind(df_total, df)
}

df_total <- cbind(df_total,seurat_clusters_DEGs)


pointSize = 10 
lineWidth = 0.5

degs_bar <- ggplot(df_total,aes(x = reorder(seurat_clusters_DEGs, count), y = count, fill = seurat_clusters_DEGs))+
  geom_col(width = 0.7, color = 'black') + scale_fill_manual(values = colors_seurat_clusters)+
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black", linewidth = 0.2),
        plot.title  = element_text(color="black", size=10),
        axis.title  = element_text(size = pointSize, colour = "black"),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
        axis.text.y  = element_text(size = pointSize , colour = "black"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.ticks = element_line(linewidth = 0.2, colour = "black"),
        axis.line = element_line(linewidth = 0.2, colour = "black"))+
  labs(y= 'Count of DEGs (Log2FC > 0.2)', x=NULL) +
  ggtitle(paste0('Number of DEGs per CellType',comparison_name,'.xlsx')) +
  scale_y_continuous( expand = c(0,0.1,0,0)) 
degs_bar

ggsave(
  paste0('plots/barplots/','Seurat_clusters_DEG_',comparison_name,'.pdf'),
  plot = degs_bar,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 150,
  height = 100,
  units = c("mm"),
  dpi = 600)


######## CHANGES IN CLUSTER PERCENTAGE - BARPLOT  ##########
#Barplot for cluster changes:
barplot_percentage_prefix <- 'CellType'
dataset_name <- 'all_data'
colors_barplot <- colors_CellType
group_barplot <- "Genotype"

Idents(dataset) <- barplot_percentage_prefix
barplot <- dittoBarPlot(dataset, 
                        var = barplot_percentage_prefix,  
                        group.by = group_barplot)
barplot_data <- as.data.frame(barplot$data)

#Select grouping MANUAL!!
barplot_data$label <- factor(barplot_data$label,levels=c('Microglia','Macrophage','MGnD','Monocytes',
                                                         "MHC2_Macrophages",'Neutrophils','Endothelial_cells',
                                                         'NK_cells','non_classical_DCs','IRM','Cycling_cells',
                                                         "Ependymal_cells", "B_cells"))
lineWidth = 1
pointSize = 10
gg_bar <- ggplot(barplot_data, aes(x=grouping, y=percent, fill = label)) +
  geom_bar(position = "fill", stat = "identity",width=0.9) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0.01,0,0)) +
  geom_text(aes(label = paste0(round(percent*100,digits = 2),"%")), 
            position = position_stack(vjust = 0.5), size = 0.5)+
  scale_fill_manual(values = colors_barplot) +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size=20, face="bold.italic"),
        axis.title  = element_text(size = pointSize, colour = "black"),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
        axis.text.y  = element_text(size = pointSize , colour = "black"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.line = element_line(size = lineWidth, colour = "black"))
gg_bar

ggsave(
  paste0('plots/barplots/',barplot_percentage_prefix,'_distribution_group_by_',group_barplot,processing_specs,dataset_name,'.pdf'),
  plot = gg_bar,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 200,
  height = 300,
  units = c("mm"),
  dpi = 600)

#Changes in percentage barplot
barplot_percentage_prefix <- 'CellType'
comparison_name <- 'xx_vs_yy'
colors_barplot <- colors_CellType
group_barplot <- "Genotype"


Idents(dataset) <- barplot_percentage_prefix
barplot <- dittoBarPlot(dataset, 
                        var = barplot_percentage_prefix,  
                        group.by = group_barplot)
barplot_data <- as.data.frame(barplot$data)
CellType <- barplot_data %>% pull(label) %>% unique()

dt_1 <- barplot_data %>% filter(grouping == "xx")
dt_2 <- barplot_data %>% filter(grouping == "yy")


pct_change_df_total <- data.frame(Cluster = character(), Pct_Change = numeric(), stringsAsFactors = FALSE)
for (i in CellType) {
  pct_1 <- dt_1 %>% filter(label == i) %>% pull(percent)
  pct_2 <- dt_2 %>% filter(label == i) %>% pull(percent)
  
  pct_change <- as.numeric(pct_2 - pct_1)*100
  
  pct_change_df <- data.frame(Cluster = i, Pct_Change = pct_change)
  pct_change_df_total <- rbind(pct_change_df_total, pct_change_df)
}

lineWidth = 0.2
pointSize = 10
gg_bar <- ggplot(pct_change_df_total, aes(x=reorder(Cluster,Pct_Change), y=Pct_Change, fill = Pct_Change)) +
  geom_bar(stat = "identity",width=0.9) +
  labs(title = comparison_name, xlab = "", ylab="Percentage (%) change" ) +
  scale_fill_gradient2(limits=c(-3,3), low="navy", mid="whitesmoke", high = "firebrick", na.value = 'navy') +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size=pointSize),
        axis.title  = element_text(size = pointSize, colour = "black"),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
        axis.text.y  = element_text(size = pointSize , colour = "black"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.line = element_line(size = lineWidth, colour = "black"))
gg_bar


ggsave(
  paste0('plots/barplots/',barplot_percentage_prefix,'_percentage_change_for_',comparison_name,processing_specs,'.pdf'),
  plot = gg_bar,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 100,
  height = 100,
  units = c("mm"),
  dpi = 600)









######## TOTAL COUNTS ########
Idents(dataset) <- 'Treatment'
dataset_ctrl <- subset(dataset, ident = 'Ctrl')
dataset_AD <- subset(dataset, ident = 'AD')

# Gene expression
Idents(dataset_filtered_xenon) <- 'SubCellType'
exprTable <- AverageExpression(dataset_filtered_xenon)
exprTable <- exprTable[['RNA']] %>% as.data.frame() 
exprTable <- exprTable %>% select(order(colnames(exprTable))) %>% select(everything())
write.xlsx(exprTable, paste0('results/geneExpression_', 'SubCellType','.xlsx'), rowNames = T, overwrite = T)


######################## VISUALIZATION ##############################################
Idents(dataset) = 'seurat_clusters'

#level active ident 
levels(dataset_filtered_no_macrophage)
dataset_filtered_no_macrophage@active.ident <- factor(dataset_filtered_no_macrophage@active.ident,
                                                      levels=c( "M0",
                                                                "Ribosome-Intermediate",
                                                                "HSP-Intermediate",
                                                                "Cytokine-Intermediate",
                                                                "Interferon",
                                                                "Cycling-S",
                                                                "Cycling-G2M",
                                                                "MGnD",
                                                                "MGnD-Antigen")) # MANUAL EDIT OF CELL TPYPES
levels(dataset_filtered_no_macrophage)

######## GENERAL DIMPLOT FEATUREPLOT ########
dir.create(paste0('plots/DimPlots/'))

levels(dataset)
Idents(dataset) <- 'RNA_snn_res.0.9' 
#Split by DIMPLOT
pdf(file = paste0('plots/DimPlots/',"DimPlot_", "CellType","split_by", "orig_ident", processing_specs, '.pdf'), width = 12, height = 20)
DimPlot(dataset, label = F, label.box = F, cols = colors_CellType, pt.size = 0.3 ,split.by = 'orig.ident', ncol = 2)
dev.off()

#GROUP by DIMPLOT
pdf(file = paste0('plots/DimPlots/',"DimPlot_", "CellType","split_by", "Genotype", processing_specs, '.pdf'), width = 10, height = 10)
DimPlot(dataset, label = F, label.box = F, cols = colors_genotypes, pt.size = 0.3 , group.by = 'Genotype')
dev.off()


##### AD data
levels(dataset_AD)
Idents(dataset_AD) <- 'CellType' 
#Split by DIMPLOT
pdf(file = paste0('plots/DimPlots/',"DimPlot_", "CellType","split_by", "orig_ident", processing_specs, '.pdf'), width = 12, height = 20)
DimPlot(dataset_AD, label = F, label.box = F, cols = colors_genotypes_AD, pt.size = 0.3 , 
        split.by = 'Genotype', ncol = 2)
dev.off()
#GROUP by DIMPLOT
pdf(file = paste0('plots/DimPlots/',"DimPlot_AD_data_", "CellType","split_by", "Genotype", processing_specs, '.pdf'), width = 10, height = 10)
DimPlot(dataset_AD, label = F, label.box = F, cols = colors_genotypes_AD,
        pt.size = 0.1 , group.by = 'Genotype',alpha = 0.2)
dev.off()


##### WT data
levels(dataset_WT)
Idents(dataset_WT) <- 'CellType' 
#Split by DIMPLOT
pdf(file = paste0('plots/DimPlots/',"DimPlot_", "CellType","split_by", "orig_ident", processing_specs, '.pdf'), width = 12, height = 20)
DimPlot(dataset_WT, label = F, label.box = F, cols = colors_genotypes_AD, pt.size = 0.3 , 
        split.by = 'Genotype', ncol = 2)
dev.off()

#GROUP by DIMPLOT
pdf(file = paste0('plots/DimPlots/',"DimPlot_WT_data_","split_by", "Genotype", processing_specs, '.pdf'), width = 10, height = 10)
DimPlot(dataset_WT, label = F, label.box = F, cols = colors_genotypes_WT,
        pt.size = 0.1 , group.by = 'Genotype',alpha = 0.2)
dev.off()



##### FeaturePlots
dir.create(paste0('plots/Featureplots/'))
pdf(file = paste0('plots/Featureplots/','xxx','_','yyy', '.pdf'), width = 10, height = 5)
FeaturePlot(dataset,features = c('xxx','yyy'), pt.size = 0.2, order =T, min.cutoff = 0.5)
dev.off()

#split by 
pdf(file = paste0('plots/Featureplots/','xxx','_', 'yyy',"split_by_","Genotype", '.pdf'), width = 10, height = 5)
FeaturePlot(dataset,features = c('xxx','yyy'),order =T, min.cutoff = 0.5, 
            split.by = "Genotype", ncol = 2)
dev.off()


######## MODULE PLOT ########
#INF
INF_module = list(c("Ifit3","Isg15","Irf7","Ifit2","Ifit1","Stat1","Usp18","Mx1","Ifi204"))
dataset <- AddModuleScore(dataset,features = INF_module,name="INF_module")
INF_module_featureplot <- FeaturePlot(dataset,features = "INF_module1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/Module_plots/Interferon_Module_', projectname, '.pdf'), width = 9, height = 9)
INF_module_featureplot
dev.off()


#M0
M0_module = list(c('P2ry12', 'Fosb', 'Ifngr1', 'Maf', 'Tmem119', 'Jund','Jun', 'Egr1'))
M0_module <- AddModuleScore(dataset,
                            features = M0_module,
                            name="M0_module")
M0_module_featureplot <- FeaturePlot(M0_module,
                                     features = "M0_module1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_M0_Module_', 'dataset_filtered_no_macrophage','.pdf'), width = 9, height = 9)
M0_module_featureplot
dev.off()


#MGnD
MGnD_module = list(c('Clec7a','Axl','Igf1','Lpl','Spp1','Cst7'))
MGnD_module <- AddModuleScore(dataset,
                              features = MGnD_module,
                              name="MGnD_module")
MGnD_module_featureplot <- FeaturePlot(MGnD_module,
                                       features = "MGnD_module1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_MGnD_module_', 'dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
MGnD_module_featureplot
dev.off()

#S-Phase
S_phase_module = list(c('Mcm5',
                     'Mcm6',
                     'Lig1',
                     'Mcm4',
                     'Mcm7',
                     'Stmn1',
                     'Mcm3'))
S_phase_module <- AddModuleScore(dataset,
                              features = S_phase_module,
                              name="S_phase_module")
S_phase_module_featureplot <- FeaturePlot(S_phase_module,
                                       features = "S_phase_module1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_S_phase_module_featureplot_', 'dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
S_phase_module_featureplot
dev.off()

#Ribosome
Ribosome_module = list(Sub_celltype_markers_filtered_top_markers %>% filter(cluster == 'Ribosome-Intermediate') %>% pull(gene))
dataset_filtered_no_macrophage <- AddModuleScore(dataset,
                                 features = Ribosome_module,
                                 name="Ribosome_module")
Ribosome_module_featureplot <- FeaturePlot(dataset_filtered_no_macrophage,
                                          features = "Ribosome_module1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_Ribosome_module_featureplot_', 'dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
Ribosome_module_featureplot
dev.off()

#G2M phase
G2M_module = list(c('Top2a', 'Mki67' ))
G2M_module <- AddModuleScore(dataset,
                                  features = G2M_module,
                                  name="G2M_module")
G2M_module_featureplot <- FeaturePlot(G2M_module,
                                           features = "G2M_module1",label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_G2M_module_featureplot_', 'dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
G2M_module_featureplot
dev.off()

#AntigenMGnD
Antigen_module = list(c('H2-Aa','H2-Eb1', 'H2-Ab1'))
Antigen_module <- AddModuleScore(dataset,
                             features = Antigen_module,
                             name="Antigen_module")
Antigen_module_featureplot <- FeaturePlot(Antigen_module,
                                      features = "Antigen_module1",label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/Featureplot_Antigen_module_','dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
Antigen_module_featureplot
dev.off()

#Macrophages
Macrophage_module = list(c('Mrc1','Clec12a', 'Ms4a7'))
Macrophage_module <- AddModuleScore(dataset,
                                 features = Macrophage_module,
                                 name="Macrophage_module")
Macrophage_module_featureplot <- FeaturePlot(Macrophage_module,
                                          features = "Macrophage_module1", label = F, pt.size = 1, cols = c('grey', 'red')) 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_Macrophage_module_featureplot_', projectname, '.pdf'), width = 12, height = 9)
Macrophage_module_featureplot
dev.off()

#HSP-module
HSP_module = list(Sub_celltype_markers_filtered_top_markers %>% filter(cluster == 'HSP-Intermediate') %>% pull(gene))
dataset_filtered_no_macrophage <- AddModuleScore(dataset,
                                    features = HSP_module,
                                    name="HSP_module")
HSP_module_featureplot <- FeaturePlot(HSP_module,features = "HSP_module1", label = T,  cols = c('grey', 'red')) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_HSP_module_', 'dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
HSP_module_featureplot
dev.off()


#Cytokine-module
Cytokine_module = list(Sub_celltype_markers_filtered_top_markers %>% filter(cluster == 'Cytokine-Intermediate') %>% pull(gene))
dataset_filtered_no_macrophage <- AddModuleScore(dataset,
                             features = Cytokine_module,
                             name="Cytokine_module")
Cytokine_module_featureplot <- FeaturePlot(Cytokine_module,features = "Cytokine_module1", label = T,  cols = c('grey', 'red')) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_Cytokine_module_', 'dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
Cytokine_module_featureplot
dev.off()



######## VIOLIN PLOT ############
dir.create(paste0('plots/Violin_plots'))

#Violin plot for specific genes across conditions
Idents(dataset) = 'Genotype'

pdf(file = paste0('plots/Violin_plots/',"xxx_yyy" ,'_','Genotype',processing_specs, '.pdf'), width = 10, height = 5)
VlnPlot(dataset, features = c("xxx","yyy"), cols = colors_genotypes)
dev.off()



dittoPlot(dataset, var = 'Apoe', group.by = "CellType", split.by = 'Treatment', plots = c('vlnplot', 'jitter'), jitter.size = 0.5, split.ncol = 1)

Module = list(c(Sub_celltype_markers_filtered %>% filter(cluster == 'MGnD-Antigen' & p_val_adj < 0.05) %>% pull(gene)))
dataset <- AddModuleScore(dataset,
                                     features = Module,
                                     name="MGnDAntigen")

Vln_plot <- ExpVlnPlot(
  srt = dataset_filtered_no_macrophage, group.by = "Treatment", 
  features = c("MGnDAntigen1"),
  comparisons = list(
    c("xx", "yy")
  ),
  multiplegroup_comparisons = FALSE
)

Vln_plot

dev.off()
pdf(file = paste0('plots/Violin_plot_', 'Intermediate_combined',':cc_vs_ee', '.pdf'), width = 12, height = 9)
Vln_plot
dev.off()

######## CUSTOM ENRICHMENT VIOLIN PLOT #########
Module_score <- dataset@meta.data$INF_module1
Treatment <- dataset@meta.data$Treatment

vln_data <- map2_dfr(Module_score, Treatment,~ tibble(Module_score = .x, Treatment = .y))

lineWidth = 0.5
pointSize = 15

Vln_plot_custom <- ggplot(vln_data, aes(x = Treatment, y = Module_score, fill = Treatment,group = Treatment)) + 
  geom_violin() +
  stat_summary(fun = "median",geom = "crossbar",width = 0.5,colour = "black")+
  geom_point(aes(x = Treatment), position = position_jitterdodge(jitter.width = 1, jitter.height=0, dodge.width=0.9), size = 0.5, alpha = 0.1) +
  theme(text = element_text(size = pointSize, colour = "black"),
    rect = element_blank(),
    line = element_line(size = lineWidth, colour = "black"),
    plot.title  = element_text(color="black", size=20),
    axis.title  = element_text(size = pointSize, colour = "black"),
    axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
    axis.text.y  = element_text(size = pointSize , colour = "black"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = pointSize, colour = "black"),
    legend.text = element_text(size = pointSize , colour = "black"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    axis.line = element_line(size = lineWidth, colour = "black"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
  labs(title = 'Pre-MGnD Interferon Module') +
  stat_compare_means(comparison =  list(c('atmosphere', 'xenon')), p.adjust.method = "bonferroni",  method = 'wilcox.test', label = "p.signif") +  
  stat_compare_means(label.x = 1.2, label.y = 2)+
  scale_y_continuous(expand = c(0, 0, .05, 0))+
  #scale_fill_manual(values = c(  "#FDBF6F",'#b37627')) #Cytokine Colors
  #scale_fill_manual(values = c("#FF7F00",'#ab5500')) #HSP Colors
  scale_fill_manual(values = c("#33A02C",'#10360e')) #HSP Colors
  
Vln_plot_custom


dev.off()
pdf(file = paste0('plots/Violin_plot_', 'Custom_Module_enrichment','_cc_vs_tt', '.jpeg'), width = 5, height = 9)
Vln_plot_custom
dev.off()


######## HEATMAP #########
DoHeatmap(dataset_filtered,
          features = TOP_dataset_filtered_subcelltype_markers$gene,
          raster = TRUE)

ht <- ExpHeatmap(
  srt = dataset, features = TOP_dataset_filtered_subcelltype_markers$gene, feature_split = TOP_dataset_filtered_subcelltype_markers$cluster, cell_split_by = "SubCellType",
  species = "Mus_musculus", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
  row_title_size = 0, height = 8, width = 10)
ht

pdf(file = paste0('plots/heatmap_Subcelltype_pathways_sc','.pdf'), width = 20, height = 15)
print(ht$plot)
dev.off()

ggsave(paste0('plots/heatmap_pathways_SubCellType' , '.jpeg'),
  plot = ht$plot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 500,
  height = 300,
  units = c("mm"),
  dpi = 300)

######## VOLCANO PLOT ######## 
dir.create('plots/Volcano')

#Set-up
Volcano_prefix <- "xx_vs_tt"

Volcano_data <- read.xlsx(paste0('results/DEGs_Genotype/',Volcano_prefix,'.xlsx'))
Volcano_labels <- read.xlsx(paste0('results/DEGs_Genotype/',Volcano_prefix,'.xlsx')) %>% pull("Labels")

Volcano_data$log_pval_adj <- -log10(Volcano_data$p_val_adj)
Volcano_data$log_pval_adj <- gsub(Inf, 300, Volcano_data$log_pval_adj)
Volcano_data$log_pval_adj <- as.numeric(Volcano_data$log_pval_adj)

sig_data <- Volcano_data %>% filter(p_val_adj < 0.05)
UP_data <- sig_data %>% filter(avg_log2FC > 0) 
DOWN_data <- sig_data %>% filter(avg_log2FC < 0)

##Theme
theme_scatter <- theme(aspect.ratio = 1, 
                       panel.background = element_blank(),
                       panel.border=element_rect(fill=NA),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       strip.background=element_blank(),
                       axis.title  = element_text(size = 8, colour = "black"),
                       axis.text.x=element_text(colour="black",size =8),
                       axis.text.y=element_text(colour="black",size =8),
                       axis.ticks=element_line(colour="black", size = 0.2),
                       plot.margin=unit(c(1,1,1,1),"line"),
                       plot.title = element_text(colour="black",size =10))

#Select data
select_UP_data <- UP_data %>% filter(gene %in% Volcano_labels)
select_DOWN_data <- DOWN_data %>% filter(gene %in% Volcano_labels)

scatter_volcano <- ggplot(Volcano_data, aes(x= avg_log2FC, y= log_pval_adj, label = gene)) +
  geom_point(aes(fill = avg_log2FC),size = 1, shape = 21, color = 'black', stroke=0.01 ) +
  scale_fill_gradient2(limits=c(-3,3), low="navy", mid="whitesmoke", high = "firebrick", na.value = 'firebrick') +
  theme_scatter +
  xlim(-7,7)+
  ylim(0,350)+
  geom_vline(xintercept = 0.1,colour="grey", linetype = "longdash", linewidth = 0.5) +
  geom_vline(xintercept = -0.1,colour="grey", linetype = "longdash", linewidth = 0.5) +
  geom_hline(yintercept = 1.3, colour="grey", linetype = "longdash", linewidth = 0.5) +
  labs(title = str_wrap(paste0(Volcano_prefix),60),
       x = "Log2FC",
       y = "-log(p-adj)") +
  geom_text_repel(data = select_UP_data,
                   color = 'black', size = 3,fontface = 'italic',
                   min.segment.length = 0,
                   # nudge_y = 30,
                   nudge_x = 3,
                   max.overlaps = Inf,
                   point.padding = unit(0.5, 'mm')) +
  geom_text_repel(data = select_DOWN_data,
                   color = 'black',size = 3, fontface = 'italic',
                   min.segment.length =0,
                   # nudge_y = 30,
                   nudge_x = -2,
                   max.overlaps = Inf,
                   point.padding = unit(0.5, 'mm'))
scatter_volcano 
dev.off()
pdf(file = paste0('plots/Volcano/',Volcano_prefix, '.pdf'), width = 5, height = 5)
scatter_volcano
dev.off()



######## DOTPLOT CELLTYPES ##########
Dotplot <- DotPlot(dataset, features = c(
  "Tmem119", "P2ry12","Sall1", "Sparc","Gpr34","Crybb1",                 #Microglia
  "Ifit3","Usp18","Isg15","Irf7",                                               #IRM
  "Cst7", "Lpl","Axl","Clec7a",                                          #MGnD
  "Mrc1","Ms4a7","Dab2","Pf4", "Clec12a",                                #Macrophages
  "H2-Aa", "Slamf7", "Stap1", "Cd72","Ciita", "H2-Eb1","H2-Ab1",         #MHCII_Macrophages
  "Tnip3","Cd209a","H2-DMb2","P2ry10",                                   #Non_classical_DCs
  "Napsa" , "S100a4","Cytip","Plac8",  "Ace",  "Itgal",                  #Monocytes
  "Ifitm6","Hp", "Trem1","Mmp9","Hdc","Camp","Lcn2","Sell" ,             #Neutrophils
  "Cd19", "Cd79a", "Cd79b",                                              #B_cells
  "Nkg7", "Gzma", "Gzmb", "Ifng","Il2rb","Prf1" ,                        #NK_cells
  "Flt1","Cldn5","Ptprb" ,                                               #Endothelial_cells
  "Igfbp7","Cald1","Myl9","Tpm2", "Pdgfrb", "Acta2",                     #Ependymal_cells
  "Stmn1","Birc5","Top2a","Mki67"))                                      #Cycling_cells

Dotplot_df <- Dotplot$data

Dotplot_df$id <- factor(Dotplot_df$id ,levels=c('Microglia','IRM','MGnD','Macrophage',"MHC2_Macrophages",
                                                'non_classical_DCs','Monocytes','Neutrophils', "B_cells",
                                                'NK_cells','Endothelial_cells',"Ependymal_cells",'Cycling_cells'))

pointSize = 12
lineWidth = 0.2

dotplot_markers <- ggplot(Dotplot_df, aes(x = id, y = features.plot)) + 
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size=20),
        axis.title.x =  element_text(size = pointSize , colour = "black"),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 60, hjust=1),
        axis.text.y  = element_text(size = pointSize , colour = "black", face = 'italic'),
        axis.line = element_line(size = lineWidth, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"))+
  scale_colour_gradient2(limits=c(-3,3), low="navyblue", mid="white", high = "firebrick") +
  labs(title = str_wrap("Celltype marker genes", 40))

dotplot_markers

ggsave(
  paste0('plots/Dotplots/', 'Marker_genes_','Celltype' ,processing_specs, '.pdf'),
  plot = dotplot_markers,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 250,
  height = 350,
  units = c( "mm"),
  dpi = 600
)



################################### FUNCTIONAL ASSESSMENT ##############################################
######## TRAJECTORY ASSESSMENT ####################################
Idents(dataset) <- 'SubCellType'
dataset_subset <- subset(dataset, idents = c('M0','MGnD','Cytokine-Intermediate',
                                                               'Intermediate'))
Idents(dataset) <- 'Treatment'
dataset_filtered_subset_xx <- subset(dataset, idents = 'xenon')
dataset_filtered_subset_yy <- subset(dataset, idents = 'atmosphere')

dataset_filtered_subset_xx.slingshot <- RunSlingshot(srt = dataset_filtered_subset_xx, group.by = "SubCellType_PreMGnD", reduction = "umap", start = 'M0')

# dataset_filtered_subset_slingshot@active.ident <- factor(dataset_filtered_subset_slingshot@active.ident,
#                                                         levels=c('M0','Ribosome-Intermediate','HSP-Intermediate',
#                                                                  'Interferon','MGnD-Antigen',
#                                                                  'MGnD','Cytokine-Intermediate'))
# levels(dataset_filtered_subset_slingshot)

pdf(file = paste0('plots/Trajectory_plot_','dataset_filtered_subset_air.slingshot','_SubCellType_PreMGnD' , projectname, '.pdf'), width = 10, height = 10)
Trajectory_plot <- ClassDimPlot(dataset_filtered_subset_xx.slingshot, group.by = "SubCellType_PreMGnD", reduction = "umap", lineages = paste0("Lineage", 1:2), lineages_span = 0.3,lineages_line_bg = "white", palette = "Paired")
dev.off()

Trajectory_plot

ggsave(
  paste0('plots/Trajectory_plot_','dataset_filtered_subset_xx.slingshot','_SubCellType_PreMGnD' , '.pdf'),
  plot = Trajectory_plot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 300,
  height = 300,
  units = c("mm"),
  dpi = 300)


#Dynamic Plot
dynamicplot <- DynamicPlot(
  srt = dataset_filtered_subset_xx.slingshot, lineages = c("Lineage1", "Lineage2"), group.by = "SubCellType_PreMGnD",
  exp_method = 'log1p',
  features = c('Tmem119','P2ry12', 'Stat1','Clec2d'))

dynamicplot

ggsave(
  paste0('plots/dataset_filtered_subset_xx.slingshot','Cst7_Clec7a','PreMGnD.pdf'),
  plot = dynamicplot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 300,
  height = 300,
  units = c( "mm"),
  dpi = 300)

#Expression Plots
ExpDimPlot_trajectory <- ExpDimPlot(dataset_filtered_subset_xx.slingshot, features = paste0("Lineage", 1:2), reduction = "UMAP", theme_use = "theme_blank")
ExpDimPlot_trajectory

ggsave(
  paste0('plots/ExpDimPlot_trajectory_','dataset_filtered_subset_xx.slingshot' ,'.pdf'),
  plot = ExpDimPlot_trajectory,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 300,
  height = 300,
  units = c( "mm"),
  dpi = 300)
