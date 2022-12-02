#### Roche Nimblegen data analysis example
# Create working environment
homedir <- getwd()
setwd(homedir)
options(java.parameters = "-Xmx8G")
library(BiocParallel)
library(EnhancedVolcano)
library(enrichR)
library(fgsea)
library(gplots)
library(ggplot2)
library(limma)
library(magrittr)
library(msigdbr)
library(pheatmap)
library(stringr)
require(vulcan)
library(xlsx)
register(SnowParam(12))
dir.create("plots")
dir.create("results")
annotation <- read.table("rawdata/OID39197_Experimental_Report.txt",
                         sep = "\t", as.is = TRUE, header = TRUE)
# Load RMA normalized spot intensities
raw <- read.table("rawdata/processed_data_files/Normalized_Calls_Files/100718_HG18_opt_expr_RMA.calls",
                  sep = "\t", as.is = TRUE, header = TRUE)
