library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(xlsx)

#Liver dataset will be used as an example
Liver_rawdata <- Read10X("~/Documents/Liver/filtered_gene_bc_matrices/mm10_TdTomato")
Liver <- CreateSeuratObject(counts = Liver_rawdata,  project = "Liver", min.cells = 3,  min.features = 200)   

#Read in file with list of dissociation associated genes from O'Flanagan et al. 2019
dag <- read.csv("Dissociation_genes.csv", sep= ",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

#find which DAGs actually overlap with genes present in our datasets
#using "gene" column which contains the 507 upregulated DAGs
overlap_liver <- intersect(rownames(Liver), dag$gene)

#find the percentage of transcripts that map to either mitochondiral genes or to DAGs
Liver[["percent.mito"]] <- PercentageFeatureSet(Liver, pattern = "^mt") 
Liver[["percent.dag"]] <- PercentageFeatureSet(Liver, features = overlap_liver) 

#Show violin plots of various parameters in dataset
VlnPlot(Liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.dag"), ncol = 4)

#Show scatter plots of the same parameters as above
plot1 <- FeatureScatter(Liver, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(Liver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(Liver, feature1 = "nCount_RNA", feature2 = "liver.percent.dag")
plot4 <- FeatureScatter(Liver, feature1 = "percent.dag", feature2 = "percent.mito")
CombinePlots(plots = list(plot1, plot2, plot3, plot4))

#Filter data using cutoffs determined from the above scatterplots using the subset function
liver_filtered <- subset(Liver, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito < 20 & liver.percent.dag < 30)

#Normalize data and find variable features
Liver <- NormalizeData(Liver)
Liver <- FindVariableFeatures(object = Liver, selection.method = 'mean.var.plot', mean.cutoff = c(0.04, 3), dispersion.cutoff = c(0.6, Inf))
length(x = VariableFeatures(object = Liver))
VariableFeaturePlot(Liver)

#Scale data and run PCA
Liver <- ScaleData(object = Liver, features = rownames(x = Liver), vars.to.regress = c("nCount_RNA", "percent.mito"))
Liver <- RunPCA(object = Liver, features = VariableFeatures(object = Liver), verbose = FALSE)
DimPlot(object = Liver, reduction = "pca")

#Choosing the number of PCs for downstream analysis
ElbowPlot(object = Liver)

#Clustering and non-linear dimensionality reduction using tSNE and UMAP
Liver <- FindNeighbors(object = Liver, dims = 1:12)
Liver <- FindClusters(object = Liver, resolution = 1.6)
Liver <- RunUMAP(Liver, dims = 1:12, n.neighbors = 30, min.dist = 0.3)
DimPlot(object = Liver, reduction = 'umap', label = FALSE)

#FindMarkers and Heatmaps
Liver.markers <- FindAllMarkers(Liver, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
Liver_sig <- Liver.markers[which(Liver.markers$p_val_adj < 0.05), ]
Liver_sig <- Liver_sig[!grepl("^mt-", rownames(Liver_sig)), ]
Liver_sig <- Liver_sig[!grepl("^Rp", rownames(Liver_sig)), ]
write.xlsx(Liver.markers, "Liver_DEGs.xlsx")

#Find the top30 based on avg_logFC. In the table, they will be ordered based on pvalue, but they are the top 30 based on avg_logFC
top30_Liver <- Liver_sig %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
write.xlsx(top30_Liver, "Liver_top30.xlsx")

#create Heatmap of the top 30 genes
DoHeatmap(object = subset(Liver, cells = WhichCells(Liver, downsample = 50, seed = 111)), features = top30_Liver$gene) + NoLegend()

#FeaturePlots and VlnPlots
#Create Vln and FeaturePlots for two genes, Fcgr1 and Mafb

VlnPlot(Liver, features = c("Fcgr1", "Mafb"), ncol = 3, pt.size = FALSE)
FeaturePlot(Liver, features = "Fcgr1", order = TRUE, min.cutoff = 1, max.cutoff = 3)

#Subsetting and merging Seurat objects
#in this example, subsetting the liver object to take only 3 clusters, then checked number of cells after using object@assays command
Liver_mac <- subset(x = Liver, idents = c(0, 1, 2))
Liver_mac@assays

#merging 5 Seurat objects
Merged <- merge(x = Liver, y = c(Kidney, Heart, Brain, Lung))




