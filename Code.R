### Spatial Transcriptomics Test
# Re-creating the sample analysis I did using online tutorials

# Load packages
library(tidyverse)
library(readr)
library(SeuratObject)
library(here)

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
sample_data <- Load10X_Spatial(data_dir = data_dir1)






