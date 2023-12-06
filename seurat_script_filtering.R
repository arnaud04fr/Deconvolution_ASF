# Load required libraries
library(Seurat)
library(dplyr)

# Load single-cell RNA-seq data
counts <- Read10X(data.dir = "./GSE176171") # Read data from the working directory

# Load cell metadata
meta_data <- read.table("./GSE176171/cell_metadata.tsv.gz", header = TRUE, 
                        sep = "\t",
                        quote = "", 
                        row.names = 1, 
                        stringsAsFactors = FALSE)

# Create a Seurat object
sc_data <- CreateSeuratObject(counts = counts,
                              assay="RNA",
                              meta.data = meta_data,
                              min.cells = 3, 
                              min.features = 200,
                              project = "GSE176171_SeuratObj")
#remove counts
rm(counts)

#filter cells that have unique feature counts over 2,500 or less than 200 and >5% mitochondrial counts
sc_data <- subset(sc_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & mitochondrial_percent < 5)

# Filter out ribosomal and mitochondrial genes
ribosomal_genes <- grep("^RPS", rownames(sc_data@assays$RNA@counts), value = TRUE)
mitochondrial_genes <- grep("^MT-", rownames(sc_data@assays$RNA@counts), value = TRUE)
sc_data <- subset(sc_data, features = !rownames(sc_data@assays$RNA) %in% c(ribosomal_genes, mitochondrial_genes))

# Normalize and scale data
sc_data <- NormalizeData(sc_data)
# Select the top variable genes
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 5000)
# Scale data
sc_data <- ScaleData(sc_data)

#linear dimensional reduction
sc_data <- RunPCA(sc_data, features = VariableFeatures(object = sc_data))
DimPlot(sc_data, reduction = "pca")

#determine dimension
ElbowPlot(sc_data) #fast method

#clustering 
sc_data <- FindNeighbors(sc_data, dims = 1:15) #dims -> look at ElbowPlot or paper
sc_data <- FindClusters(sc_data, resolution = 0.5) #resolution 0.4 - 1.2

#UMAP
sc_data <- RunUMAP(sc_data, dims = 1:15)
DimPlot(sc_data, reduction = "umap",
        label = TRUE)

#save file
saveRDS(sc_data, file = "./scdata_GSE176171_filtered-5000feat.rds")

#explore Datasets with gene names
VlnPlot(sc_data, features = c("PRLR","ADIPOQ","PROX1","PDGFRA" ))
FeaturePlot(sc_data, features = c("PRLR", "ADIPOQ","PROX1","IL7R"))


######################################################"
#load from .rds if crash....
#library(Seurat)
#sc_data <- readRDS(file = "scdata_GSE176171.rds")

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
WATmarkers <- FindAllMarkers(sc_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
WATmarkers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
 
write.table(WATmarkers,file="Clusters.tsv",sep = "/", quote = F, col.names = T) 
