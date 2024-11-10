#Import libraries
# Seurat is an R package designed for QC, analysis, and exploration of single-cell RNA-seq data. 
# tidyverse is an R programming package assists with data import, tidying, manipulation, and data visualization. 
#ggplot2 is an R package for producing statistical, or data, graphics
library(Seurat)
library(tidyverse)
library(ggplot2)

# Reading Data:
data <- read.csv("Brain_Microglia-counts.csv", row.names = 1)
meta <- read.csv("metadata_FACS.csv")
annotations <- read.csv("annotations_FACS.csv")

# Filteration of the original annotations and meta datasets to include only the entries where the tissue is "Liver",
#The filtered datasets are saved into new variables: annotations_filtered and meta_filtered.
annotations_filtered <- filter(annotations, tissue == "Brain_Microglia")
meta_filtered <- filter(meta, tissue == "Brain_Microglia")



#Creating a Seurat object, which is a specialized data structure that holds and
#organizes the single-cell RNA-seq data and metadata.
s_o_data <- CreateSeuratObject(data)
#calculating the percentage of mitochondrial gene expression for each cell and stores it as a new variable (percent_mito) in the Seurat object (seurate_obj).
s_o_data[["percent_mito"]] = PercentageFeatureSet(s_o_data, pattern = "ˆMT-")

s_o_data@meta.data$IDs <- rownames(s_o_data@meta.data)

#In a Seurat object, the row names of meta.data typically represent the cell identifiers (e.g., cell barcodes or cell names).

VlnPlot(s_o_data, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, )
# Visualization of the distribution of a continuous variable,
#which allows to visually inspect the distribution of these quality control metrics across all the cells in the dataset


##creating a scatter plot that shows the relationship between two features in the dataset.
FeatureScatter(s_o_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

genes <- rownames(s_o_data)
s_o_data1 <- s_o_data %>%  #applying a series of functions sequentially to the Seurat object (s_o_data), modifying it step by step.
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 3000) %>% #common QC step to remove low-quality cells
  NormalizeData() %>% #Ensuring that the data is comparable across cells.
  # Identifying the top 2000 most variable genes in the dataset
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = genes) %>% #making the data mean-centered and variance-stabilized
  RunPCA() #Reducing the dimensionality and extracting the most important features (principal components
rm(genes)

top10 <- head(VariableFeatures(s_o_data1), 10)#This selects the top 10 most variable features.
plot1 <- VariableFeaturePlot(s_o_data1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
plot2

s_o_data1 <- RunPCA(s_o_data1, features = VariableFeatures(object = s_o_data1))
# Examine and visualize PCA results a few different ways
print(s_o_data1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(s_o_data1, dims = 1:2, reduction = "pca")
DimPlot(s_o_data1, reduction = "pca") + NoLegend()
#DimHeatmap isDimHeatmap is used to observe how gene expression varies across the different dimensions of the data.
DimHeatmap(s_o_data1, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(s_o_data1, dims = 1:10, cells = 500, balanced = TRUE)
#Elbow plot visualizes the variance explained by each principal component (PC) in a scree plot format.
ElbowPlot(s_o_data1, ndims = 50)

s_o_data1 <- s_o_data1 %>%
  RunUMAP(dims = 1:24) %>% #visualize high-dimensional data in lower dimensions
  FindNeighbors(dims = 1:24) %>% #calculates the pairwise distances between cells in the high-dimensional space
  FindClusters(resolution = 0.25)
# Look at cluster IDs of the first 5 cells
head(Idents(s_o_data1), 5)
Idents(s_o_data1)

# find all markers of cluster 3
cluster3.markers <- FindMarkers(s_o_data1, ident.1 = 3)
head(cluster2.markers, n = 5)

#4 UMAP-visualization with separation by gender

s_o_data1@meta.data <- s_o_data1@meta.data  %>% mutate(gender = sub(pattern = "\\w\\d+.\\w+\\d+.\\d+_\\d+_", "",x = IDs)) %>% #It replaces part of a string in the IDs column with an empty string ("").
  mutate(gender = sub(pattern = ".1.1", "", x = gender))# processes the gender column by using sub
rownames(s_o_data1@meta.data) <- s_o_data1@meta.data$IDs # to ensure that each row of meta.data is indexed by the unique IDs of the cells
s_o_data1@meta.data
meta_filtered %>% select(c(1))#select the first column of the meta_filtered data frame.
cols = c("#F73C7D", "#1597F7", "red", "green", "orange", "brown", "blue", "darkgray", "black", "pink","yellow")
plot1 <- DimPlot(s_o_data1, reduction = "umap") +
  scale_color_manual(values = cols)
plot2 <- DimPlot(s_o_data1, reduction = "umap", label = F, group.by = "gender")
plot1+plot2
# number of male and female
Malescount <- sum(s_o_data1$gender == "M")
Femalescount <- sum(s_o_data1$gender == "F")
Malescount
Femalescount
# number of all cells
P_o_M <- (Malescount / 4558) * 100 #percentage of male 
P_o_M
P_o_F <- (Femalescount / 4558) * 100 #percentage of female
P_o_F

#5
library(Seurat)
library(dplyr)

# Extracting Gender and Cluster data from metadata
clustergender_data <- s_o_data1@meta.data %>%
  dplyr::select(seurat_clusters, gender)

# calculating the number of cells for each gender in each cluster
gender_counts <- clustergender_data %>%
  group_by(seurat_clusters, gender) %>%
  tally() %>%
  ungroup()

# calculating the total number of cells in each cluster
total_counts <- gender_counts %>%
  group_by(seurat_clusters) %>%
  summarize(total = sum(n))

# Merging Data and calculating the precentages
gender_percentage <- gender_counts %>%
  left_join(total_counts, by = "seurat_clusters") %>%
  mutate(percentage = (n / total) * 100)

head(gender_percentage)

#import ggplot2 library
library(ggplot2)

# stacked bar plot 

ggplot(gender_percentage, aes(x = factor(seurat_clusters), y = percentage, fill = gender)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Clusters", y = "Percentage of Cells (%)", fill = "Gender") +
  scale_fill_manual(values = c("blue", "pink")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

#6

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
s_o_data1.markers <- FindAllMarkers(s_o_data1, only.pos = TRUE)
s_o_data1.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
#print the first 20 rows of the data frame
head(gender_percentage, n = 20)

#Identifying markers with FindMarkers
cluster0.markers <- FindMarkers(s_o_data1, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#Plotting gene expression with VlnPlot
VlnPlot(s_o_data1, features = c("Slc39a1", "0610007L01Rik"))


#Selecting top markers for heatmap visualization
s_o_data1.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Showing with Heatmap that genes-markers separate 2 chosen sub tissues.
  
  ungroup() -> top10
DoHeatmap(s_o_data1, features = top10$gene) + NoLegend()
DoHeatmap(s_o_data1, features = c("Slc39a1", "0610007L01Rik")) + NoLegend()

