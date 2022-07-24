### Spatial Transcriptomics Test
# Re-creating the sample analysis I did using online tutorials

# Load packages
library(tidyverse)
library(readr)
library(Seurat)
library(hdf5r)
library(here)
library(patchwork)

# Common parameters
title <- 'Spatial Transcriptomics Example'
marker_genes <- list(
    'Hypothalamus'=c('Nap1l5', 'Resp18', 'Calb2', 'Zcchc12', 'Sparc'),
    'Cortex.1'=c('Tbr1','Vxn', 'Ttc9b', 'Stx1a'),
    'Cortical Subplate'=c('Nptxr', 'Slc30a3', 'Ccn3', 'Lypd1'),
    'Hippocampal Region'=c('Cabp7', 'Spink8', 'Cnih2', 'Ddn', 'Wipf3'),
    'Thalamus'=c('Tnnt1', 'Prkcd', 'Amotl1', 'Ramp3', 'Rora', 'Adarb1'),
    'Fornix System'=c('Mbp', 'Mobp', 'Cnp', 'Mal', 'Cldn11', 'Apod'),
    'Cortex.2'=c('Lamp5', 'Mef2c', 'Igfbp6', 'Rasgrf2'),
    'sAMY'=c('Penk', 'Tac1', 'Hap1', 'Wfs1'),
    'Pallidum'=c('Ecrg4', 'Clic6', 'Pltp','Dbi'),
    'Stria Medullaris'=c('Igfbp2', 'Igf2', 'Mgp', 'Ptgds', 'Hba-a1', 'Hba-a2'),
    'Ammons Horn'=c('Spink8', 'Hpca', 'Wipf3', 'Neurod6', 'Fibcd1', 'Rprml')
)

cluster_names <- c(
    '0'='Hypothalamus',
    '1'='Cortex.1',
    '2'='Cortical Subplate',
    '3'='Hippocampal Region',
    '4'='Thalamus',
    '5'='Fornix System',
    '6'='Cortex.2',
    '7'='sAMY',
    '8'='Pallidum',
    '9'='Stria Medullaris',
    '10'='Ammons Horn'
)

# Load data
data_dir1 <- here("data")
sample_data <- Load10X_Spatial(data.dir = data_dir1)

###########################################################################################
### Sample QC
###########################################################################################
# Count feature
countfeature_spatial <- VlnPlot(sample_data, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
featcount_spatial_implot <- SpatialFeaturePlot(sample_data, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(countfeature_spatial, count_spatial_implot)



# Count UMI
countumi_spatial <- VlnPlot(sample_data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
umicount_spatial_implot <- SpatialFeaturePlot(sample_data, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(countumi_spatial, umicount_spatial_implot)


# Count percent mitochondrial (Doesn't work)
sample_data <- PercentageFeatureSet(sample_data, '^mt-', col.name='percent.mt')
mtper_spatial <- VlnPlot(sample_data, features = "percent.mt", pt.size = 0.1) + NoLegend()
mtper_spatial_implot <- SpatialFeaturePlot(sample_data, features = "percent.mt") + theme(legend.position = "right")
wrap_plots(mtper_spatial, mtper_spatial_implot)

    # plotting all 3
VlnPlot(sample_data, features = c('nFeature_Spatial', 'nCount_Spatial', 'percent.mt'), ncol = 3)
SpatialFeaturePlot(sample_data, features = c('nFeature_Spatial', 'nCount_Spatial', 'percent.mt'), 
                   ncol = 3) + theme(legend.position = "right")


    # visualizing feature-feature relationships
qc1 <- FeatureScatter(sample_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
qc2 <- FeatureScatter(sample_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
qc1 + qc2



    # normalization
sample_data <- SCTransform(sample_data, assay = "Spatial", verbose = FALSE)

###########################################################################################
### Gene Expression Visualization
###########################################################################################
SpatialFeaturePlot(sample_data, features = marker_genes$Hypothalamus, pt.size.factor = 1.2)
SpatialFeaturePlot(sample_data, features = marker_genes$Cortex.1, pt.size.factor = 1.2)
SpatialFeaturePlot(sample_data, features = marker_genes$`Cortical Subplate`, pt.size.factor = 1.2)
SpatialFeaturePlot(sample_data, features = marker_genes$`Hippocampal Region`, pt.size.factor = 1.2)
SpatialFeaturePlot(sample_data, features = marker_genes$Thalamus, pt.size.factor = 1.2)
SpatialFeaturePlot(sample_data, features = marker_genes$`Fornix System`, pt.size.factor = 1.2)
SpatialFeaturePlot(sample_data, features = marker_genes$Cortex.2, pt.size.factor = 1.2)
SpatialFeaturePlot(sample_data, features = marker_genes$sAMY, pt.size.factor = 1.2)
SpatialFeaturePlot(sample_data, features = marker_genes$Pallidum, pt.size.factor = 1.2)
SpatialFeaturePlot(sample_data, features = marker_genes$`Stria Medullaris`, pt.size.factor = 1.2)
SpatialFeaturePlot(sample_data, features = marker_genes$`Ammons Horn`, pt.size.factor = 1.2)




###########################################################################################
### Clustering and Integration. Dimensionality Reduction
###########################################################################################
sample_data <- RunPCA(sample_data, assay = "SCT", verbose = FALSE)
sample_data <- FindNeighbors(sample_data, reduction = "pca", dims = 1:30)
sample_data <- FindClusters(sample_data, verbose = FALSE)
sample_data <- RunUMAP(sample_data, reduction = "pca", dims = 1:30)

    # plots
p1 <- DimPlot(sample_data, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(sample_data, label = TRUE, label.size = 3)
p1 + p2


    # looking at cells of interest
SpatialDimPlot(
    sample_data, cells.highlight = 
        CellsByIdentities(object = sample_data, idents = c(2, 1, 4, 3, 5, 8)), 
    facet.highlight = TRUE, ncol = 3)

    # Interactive
SpatialDimPlot(sample_data, interactive = TRUE)

LinkedDimPlot(sample_data)




###########################################################################################
### Identification of Spatially Variable Features
###########################################################################################
spatial_markers <- FindMarkers(sample_data, ident.1 = 3, ident.2 = 9)
SpatialFeaturePlot(object = sample_data, features = rownames(spatial_markers)[1:3], 
                   alpha = c(0.1, 1), ncol = 3)


    # using spatially variable
sample_data <- FindSpatiallyVariableFeatures(sample_data, assay = "SCT", 
                                             features = VariableFeatures(sample_data)[1:1000],
                                       selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(sample_data, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(sample_data, features = top.features, ncol = 3, alpha = c(0.1, 1))





###########################################################################################
### Heat Map
###########################################################################################
DimHeatmap(sample_data, dims = 1:15, cells = 500, balanced = TRUE)







###########################################################################################
### GSEA and Volcano Plots
###########################################################################################








