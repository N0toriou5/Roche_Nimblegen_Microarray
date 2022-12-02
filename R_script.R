#### Roche Nimblegen data analysis example
# Create working environment
homedir <- getwd()
setwd(homedir)
options(java.parameters = "-Xmx8G")
library(BiocParallel)
library(corto)
library(EnhancedVolcano)
library(enrichR)
library(fgsea)
library(gplots)
library(ggplot2)
library(limma)
library(magrittr)
library(matrixStats)
library(msigdbr)
library(pheatmap)
library(stringr)
library(viper)
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

# Boxplot of RMA normalized intensities
png("plots/000_RMA_highest.png", w = 6000, h = 3000, res = 600)
boxplot(newcounts, ylab = "RMA Normalized Intensities", xlab = "Samples",
        cex.axis = 0.7, col = "lightgrey")
dev.off()

# log2 transformation
newcounts <- log2(newcounts)

# Boxplot after log2 transformation
png("plots/000_log2RMA_highest.png", w = 6000, h = 3000, res = 600)
boxplot(newcounts, ylab = "Log2 RMA Normalized Intensities", xlab = "Samples",
        cex.axis = 0.7, col = "lightgrey")
dev.off()
save(newcounts, file = "results/rawcounts_log2.rda")
#Histogram of newcounts
png("plots/000_hist_highest.png", w = 6000, h = 3000, res = 600)
hist(newcounts, xlab = "Log2 RMA Normalized Intensities", col = "lightgrey")
dev.off()

# PCA All samples
png("plots/000_pca_highest.png", w = 2000, h = 2000, p = 45)
pca <- prcomp(t(newcounts))
totvar <- sum(pca$sdev^2)
pcavar <- ((pca$sdev^2) / totvar) * 100
x <- setNames(pca$x[, 1], colnames(newcounts))
y <- setNames(pca$x[, 2], colnames(newcounts))
plot(x, y, pch = 20, main = paste0("ABCC Dataset"),
     xlim = c(min(x)*1.5, max(x)*1.5), cex = 2,
     col = c("red", rep("blue", 2), "green", "red", "grey", 
             "yellow", "green", "yellow", rep("tomato", 2), "grey"),
     xlab = paste0("PC", 1, " (", signif(pcavar[1], 3), "%)"),
     ylab = paste0("PC", 2, " (", signif(pcavar[2], 3), "%)")
)
textplot2(x, y, colnames(newcounts), new = FALSE)
grid()
dev.off()

# Separate the analysis for cell type
# BE analysis
rawbe <- subset(newcounts, select = grepl("Be", colnames(newcounts)))

# PCA only BE cells
png("plots/000_pca_BE_highest.png", w = 2000, h = 2000, p = 45)
pca <- prcomp(t(rawbe))
totvar <- sum(pca$sdev^2)
pcavar <- ((pca$sdev^2) / totvar) * 100
x <- setNames(pca$x[, 1], colnames(rawbe))
y <- setNames(pca$x[, 2], colnames(rawbe))
plot(x, y, pch = 20, main = paste0("SK-N-BE 2(2C)"),
     xlim = c(min(x)*1.5,max(x)*1.5), cex = 2, 
     col = c(rep("red", 2), "blue", rep("green", 2), "blue"),
     xlab = paste0("PC", 1, " (", signif(pcavar[1], 3), "%)"),
     ylab = paste0("PC", 2, " (", signif(pcavar[2], 3), "%)")
)
textplot2(x, y, colnames(rawbe), new = FALSE)
grid()
dev.off()

# Establish array design
trts = factor(c(1, 1, 2, 3, 3, 2))
design <- model.matrix(~ 0 + trts)
rownames(design) <- colnames(rawbe)
colnames(design) <- c("C4mut", "WT", "B12Wt")

# diff analysis block SK-N-BE2
fit <- lmFit(rawbe, design)
contrast.matrix <- makeContrasts(C4mut-WT, B12Wt-WT, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit2$contrasts # 1 C4mut - WT, 2 B12Wt - WT
C4res <- topTable(fit2, coef = 1, adjust = "BH", number = Inf)
BeC4 <- topTable(fit2, sort.by = "none", coef = 1, adjust = "BH", n = Inf)
write.xlsx2(C4res, file = "results/DiffExp_highest.xlsx", sheetName = "Be_C4mut")
B12res <- topTable(fit2, coef = 2, adjust = "BH", number = Inf)
BeB12 <- topTable(fit2, sort.by = "none", coef = 2, adjust = "BH", n = Inf)
write.xlsx2(B12res, file = "results/DiffExp_highest.xlsx", sheetName = "Be_B12Wt", append = TRUE)

summary(decideTests(fit2, lfc = 1)) #list significant genes 2-fold change

#           C4mut - WT B12Wt - WT
# Down          431        607
# NotSig      19612      18797
# Up            336        975

# Volcano Plots
png("plots/000_volcano_BeC4mut_highest.png", w = 2500, h = 2500, res = 300)
res <- BeC4
EnhancedVolcano(res, subtitle = "",
                lab = rownames(res),
                x = 'logFC',
                y = 'adj.P.Val',
                xlim = c(-8, 8),
                ylim = c(0,5),
                title = 'C4mut vs. WT',
                pCutoff = 0.01, #0.01 cutoff
                FCcutoff = 1, # 2-fold change
                labFace = "bold",
                labSize = 3.5,
                caption = paste0('Upregulated = ', nrow(res[res$logFC>1&res$adj.P.Val<=0.01,]), ' genes',"\n",'Downregulated = ',
                                 nrow(res[res$logFC< -1&res$adj.P.Val<=0.01,]), ' genes'))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("plots/000_volcano_BeB12wt_highest.png", w = 2500, h = 2500, res = 300)
res <- BeB12
EnhancedVolcano(res, subtitle = "",
                lab = rownames(res),
                x = 'logFC',
                y = 'adj.P.Val',
                xlim = c(-8,8),
                ylim = c(0,5),
                title = 'B12wt vs. WT',
                pCutoff = 0.01, #0.01 cutoff
                FCcutoff = 1, # 2-fold change
                labFace = "bold",
                labSize = 3.5,
                caption = paste0('Upregulated = ', nrow(res[res$logFC>1&res$adj.P.Val<=0.01,]), ' genes',"\n",'Downregulated = ',
                                 nrow(res[res$logFC< -1&res$adj.P.Val<=0.01,]), ' genes'))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# GSEA per contrast
msigdbr_species()
mdf <- msigdbr(species = "Homo sapiens", category = "C2") # Retrieve all human gene sets
length(mdf$gs_name) # Number of associations 520069
length(length(unique(mdf$gs_name))) # Number of pathways
head(mdf)
mlist <- mdf %>% split(x = .$gene_symbol, f = .$gs_name)
grep("MIGRATION", names(mlist), value = "TRUE")
plength <- sapply(mlist, length)
max(plength) #1932

# C4
signature <- setNames(BeC4[, "t"], rownames(BeC4))
set.seed(1)
gseas <- fgseaMultilevel(pathways = mlist, stats = signature, minSize = 5, maxSize = Inf, eps = 0, nproc = 12)
gseas <- gseas[order(gseas$pval),]
write.xlsx2(gseas, file = "results/000_GSEA.xlsx", sheetName = "C4mut", row.names = FALSE, append = FALSE)
save(gseas, file = "results/000_GSEA_C4mut.rda")

# BeB12
signature <- setNames(BeB12[, "t"], rownames(BeB12))
set.seed(1)
gseas <- fgseaMultilevel(pathways = mlist, stats = signature, minSize = 5, maxSize = Inf, eps = 0, nproc = 12)
gseas <- gseas[order(gseas$pval),]
write.xlsx2(gseas, file = "results/000_GSEA.xlsx", sheetName = "B12Wt",
            row.names = FALSE, append = TRUE)
save(gseas, file = "results/000_GSEA_B12Wt.rda")

# get a quick look at diff expression 
#https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#examining-the-number-of-de-genes
#genes having a log-FC greater than 1 (equivalent to a 2-fold difference on the original scale)
summary(decideTests(fit2, lfc = 1)) 
tfit <- treat(fit2, lfc = 1)
#dt <- decideTests(tfit)
dt <- decideTests(fit2, lfc = 1)
summary(dt)
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
png("plots/000_BE_Venn_highest.png", w = 5000, h = 2500, p = 45)
vennDiagram(dt[,1:2], circle.col = c("blue", "red"), 
            main = "SK-N-BE2 (C2) Differential Expression", cex.main = 2, line = -3, 
            sub = "2-Fold Change, FDR 5%", cex.sub = 1.5)
dev.off()

# 1072 genes are exclusively DE in B12WT
common <- rownames(dt[de.common,])
wt <- which(dt[,2]!=0)
sing <- rownames(dt[wt,])
onlyinbewt <- setdiff(sing, common)
onlysig <- BeB12[onlyinbewt,]
list1 <- rownames(onlysig[onlysig$logFC>0,])
list2 <- rownames(onlysig[onlysig$logFC<0,])
save(list1, list2, file="results/genelists.rda")

### GO Analysis
dbs <- listEnrichrDbs()
head(dbs)
#dbs <- "KEGG_2019_Human"
#dbs <- "ENCODE_TF_ChIP-seq_2015"
dbs <- c("MSigDB_Hallmark_2020","GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018",
         "Reactome_2016","KEGG_2019_Human","BioCarta_2016","WikiPathways_2019_Human","BioPlanet_2019","TF_Perturbations_Followed_by_Expression")
# 
sigB12up <- rownames(onlysig[onlysig$logFC > 0,])
sigB12dn <- rownames(onlysig[onlysig$logFC < 0,])
eup <- enrichr(sigB12up, dbs)
edown <- enrichr(sigB12dn, dbs)

# BioPlanet
up <- eup$BioPlanet_2019
down <- edown$BioPlanet_2019
up$type <- "up"
down$type <- "down"

up <- up[c(1:10),]
up <- up[order(up$Combined.Score),]
down <- down[c(1:10),]
down$Combined.Score <- (-1)*down$Combined.Score
down <- down[order(down$Combined.Score),]
gos <- rbind(down,up)
toplot <- setNames(gos$Combined.Score,gos$Term)
# Format
#names(toplot)<-gsub("GO_","",names(toplot))
#names(toplot)<-gsub("_"," ",names(toplot))
#names(toplot)<-str_to_title(names(toplot))
png("plots/000_be12wt_BioPlanet.png", w = 6000, h = 3000, res = 500)
par(mar = c(4, 1, 3, 1))
bp <- barplot(toplot, horiz = TRUE, xlab = "Combined Score",
            xlim = 1.3*c(-max(abs(toplot)), max(abs(toplot))),
            main = "B12WT top Pathways",
            col = rep(c("cornflowerblue", "salmon"), each = 10),
            yaxt = "n", cex.main = 2 
)
text(0, bp[1:10,1], names(toplot)[1:10], pos = 4)
text(0, bp[11:20,1], names(toplot)[11:20], pos = 2)
dev.off()

#Reactome
up <- eup$Reactome_2016
down <- edown$Reactome_2016
up$type <- "up"
down$type <- "down"
up <- up[c(1:10),]
up <- up[order(up$Combined.Score),]
down <- down[c(1:10),]
down$Combined.Score <- (-1)*down$Combined.Score
down <- down[order(down$Combined.Score),]
gos <- rbind(down, up)
toplot <- setNames(gos$Combined.Score, gos$Term)

# Format
names(toplot) <- gsub("Homo.*", "", names(toplot))
names(toplot)[10] <- "Neurotransmitter Receptor Binding And Downstream Transmission"
#names(toplot)<-gsub("_"," ",names(toplot))
#names(toplot)<-str_to_title(names(toplot))
png("plots/000_be12wt_Reactome.png", w = 6500, h = 3000, res = 500)
par(mar = c(4, 1, 3, 1))
bp <- barplot(toplot, horiz = TRUE, xlab = "Combined Score",
            xlim = 1.3*c(-max(abs(toplot)), max(abs(toplot))),
            main = "B12WT top Pathways",
            col = rep(c("cornflowerblue", "salmon"), each = 10),
            yaxt = "n", cex.main = 2
)
text(0, bp[1:10,1], names(toplot)[1:10], pos = 4)
text(0, bp[11:20,1], names(toplot)[11:20], pos = 2)
dev.off()

# GSEA plots 
gsea <- as.data.frame(gseas)
gseas <- gseas[order(-abs(gseas$NES)),]
let <- grep("MESENCHYMAL|INVASION|MIGRATION", gseas$pathway, value = "TRUE")
pathways <- gseas[grep("MESENCHYMAL|INVASION|MIGRATION", gseas$pathway),]
pathways <- pathways[pathways$padj <= 0.05,]

### Table of top pathways ----
top <- gseas[gseas$NES<0,][1:15,]
top <- rbind(top,gseas[gseas$NES>0,][15:1,])
top <- top[order(top$NES),]
toplot <- setNames(top$NES,top$pathway)

# Format
#names(toplot)<-gsub("GO_","",names(toplot))
names(toplot) <- gsub("_"," ",names(toplot))
names(toplot) <- str_to_title(names(toplot))
png("plots/009_gsea_be12wt.png", w = 6000, h = 3000, res = 500)
par(mar = c(4, 1, 3, 1))
bp <- barplot(toplot, horiz = TRUE, xlab = "Normalized Enrichment Score",
            xlim = 1.3*c(-max(abs(toplot)), max(abs(toplot))),
            main = "BE12WT vs. Ctrl, top Pathways",
            col = rep(c("cornflowerblue", "salmon"), each = 15),
            yaxt = "n", cex.main = 2
)
text(0, bp[1:15,1], names(toplot)[1:15], pos = 4)
text(0, bp[16:30,1], names(toplot)[16:30], pos = 2)
dev.off()

#heatmap of top 100 genes
library(gplots)
Be2topgenes <- topTreat(tfit, coef = 2, n = Inf)
head(Be2topgenes)
Be2topgenes.namegenes <- rownames(Be2topgenes)[1:100]

par(mfrow = c(1, 2))
v <- voom(rawbe, design, plot = TRUE)
v
tmp <- rownames(v$E)
v$genes <- tmp
i <- which(v$genes %in% Be2topgenes.namegenes)

### Plot heatmap
mat_zscore <- t(scale(t(rawbe[i,])))
mat_zscore <- mat_zscore[, c(3, 6, 1, 2, 4, 5)]
#mat_zscore<-mat_zscore[genes,]
png("plots/000b_heatmap_genes_BE_highest.png", w = 2000, h = 5500, res = 350)
pheatmap(mat_zscore, cluster_cols = F, cluster_rows = T)
dev.off()

### Master Regulator Analysis (corto and msViper)
# Create treatment and control expression matrices
load("results/rawcounts_log2.rda")
rawdata <- newcounts[, c(6, 12, 2, 3, 10, 11, 4, 8, 1, 5, 7, 9)]
rawbe <- subset(rawdata, select = grepl("Be", colnames(rawdata)))
rawbe <- rawbe[rowVars(rawbe)>0,]
emt <- c("TWIST1", "SNAI1", "SNAI2", "SNAI3", "ZEB1", "ZEB2", "HIF1A", "PRRX1", "SOX4",
       "SOX9", "FOXC2", "YAP1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "SIRT1", "SREBF2")

# Be B12wt MRA
trt <- rawbe[, 5:6]
ctr <- rawbe[, 1:2]

#### Create a TARGET regulon with corto
load("D:/Datasets/Expression/Neuroblastoma/target_NBL-expmat.rda") #27905   168
load("D:/Datasets/Expression/Liste/tfs_2022.rda")
regulon <- corto(expmat, centroids = centroids, nbootstraps = 1000,
                 p = 1e-06,
                 nthreads = 4,
                 verbose = TRUE)
save(regulon, file = "results/TARGET_regulon.rda")
