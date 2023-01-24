#### Pre-processing GeoMx DSP data ####
## author: Antonietta Salerno
## date: 06/01/2022

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)

### 1 - Load data ####
setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
datadir <- "TEPA_data"
DCCFiles <- dir(file.path(datadir, "DCC"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(file.path(datadir), pattern = ".pkc$",
                                full.names = TRUE, recursive = TRUE)
SampleAnnotationFile <-
  dir(file.path(datadir), pattern = "annotation.xlsx$",
      full.names = TRUE, recursive = TRUE)

data <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                               pkcFiles = PKCFiles,
                               phenoDataFile = SampleAnnotationFile,
                               phenoDataSheet = "Template",
                               phenoDataDccColName = "Sample_ID",
                               protocolDataColNames = c("aoi", "roi"),
                               experimentDataColNames = c("panel"))
  
### 2 - Study Design ####
  
  
  
  
