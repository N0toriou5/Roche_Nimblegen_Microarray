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

### Create an annotation file
annotation <- read.table("rawdata/OID39197_Experimental_Report.txt",
                         sep = "\t", as.is = TRUE, header = TRUE)

# Load RMA normalized spot intensities
raw <- read.table("rawdata/processed_data_files/Normalized_Calls_Files/100718_HG18_opt_expr_RMA.calls",
                  sep = "\t", as.is = TRUE, header = TRUE)

# Extract gene name (Entrez ID)
tmp <- strsplit(raw[, 15], ";")
tmp <- sapply(tmp, function(x){
  return(x[2])
})
Gene_id <- gsub(" gene_id ", "", tmp)

# Extract gene symbol
tmp <- strsplit(raw[, 15], ";")
tmp <- sapply(tmp, function(x){
  return(x[4])
})
Gene_symbol <- gsub(" gene_name ", "", tmp)

# Annotate sample names

names <- annotation$CY3_SAMPLE_NAME
names <- gsub(" ", "_", names)
colnames(raw)[3:14] <- names 

# add annotated gene name 

raw <- cbind(Gene_id, Gene_symbol, raw)
subraw <- raw[, -4]
rownames(subraw) <- subraw$SEQ_ID

# load squish function to merge probesets (functions included in the repository)
source("squish.R")
source("geneids.R")

# # create convlist (colnames<-probesets;associate gene_symbols in row)
# same thing but with entrez ids
removed <- subraw[subraw$Gene_symbol!="N/A",]
symbol_list <- as.character(removed[, 1])
names(symbol_list) <- removed[, 3]
input <- as.matrix(removed[, 4:15])
newcounts <- squish(input, symbol_list, method = "highest", verbose = TRUE)
dim(newcounts) # 23720 genes, 12 samples
rownames(newcounts) <- eg2sym(rownames(newcounts))
newcounts <- newcounts[!is.na(rownames(newcounts)),] # 20379    12
save(newcounts, file = "results/rawcounts_highest.rda")
