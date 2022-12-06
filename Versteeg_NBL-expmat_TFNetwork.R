### New Versteeg network
### GSE16476
### 88 samples

setwd("D:/Projects/ABCC3/")
library(affycoretools)
library(hgu133plus2.db)
library(oligo)
library(pd.hg.u133.plus.2)
celFiles <- list.celfiles("GSE16476_RAW/", listGzipped=TRUE)
affyRaw <- read.celfiles(paste0("GSE16476_RAW/",celFiles))
eset <- rma(affyRaw)
# annotate it using the ChipDb package
eset.anno <- annotateEset(eset, hgu133plus2.db)
str(eset.anno)
eset.anno@assayData
expmat <- exprs(eset.anno)

# # create convlist (colnames<-probesets;associate gene_symbols in row)
# same thing but with symbol IDs
source("D:/Archive/squish.R")
symbol_list <- as.character(eset.anno@featureData@data$SYMBOL)
names(symbol_list) <- as.character(eset.anno@featureData@data$PROBEID)
input <- as.matrix(expmat)
newcounts <- squish(input, symbol_list, method = "highest", verbose = TRUE)
dim(newcounts) # 21367 genes, 88 samples
expmat <- newcounts
save(expmat, file = "results/versteeg_NBL-expmat.rda")

### Create TF-centered network
load("D:/Datasets/Expression/Liste/tfs_2022.rda") # A list of TFs
regulon <- corto(expmat, centroids = centroids, nbootstraps = 1000,
                 p = 1e-06,
                 nthreads = 8,
                 verbose = TRUE)
save(regulon, file = "results/versteeg_regulon.rda")
