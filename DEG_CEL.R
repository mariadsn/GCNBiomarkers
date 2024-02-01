## Preprocesed 2 .CEL

library("BiocManager")
library("affy")
library("gcrma")
library("affyPLM")
library("limma")
library("annaffy")
library("arrayQualityMetrics")

# download files .CEL (raw data)
SDRF <- read.delim("../data_info.csv", sep = ",")
rownames(SDRF) <- SDRF$Ids
SDRF <- AnnotatedDataFrame(SDRF)

microarray.raw.data <- ReadAffy(filenames = SDRF$Array.Data.File, verbose=TRUE, phenoData = SDRF)
microarray.raw.data
cdfName(microarray.raw.data) 
stopifnot(validObject(microarray.raw.data))

# 1. quality analysis
# physical damage microarray 
image(microarray.raw.data[,1], col=rainbow(100))
image(microarray.raw.data[,2], col=rainbow(100))
image(microarray.raw.data[,3], col=rainbow(100))
image(microarray.raw.data[,4], col=rainbow(100))
image(microarray.raw.data[,5], col=rainbow(100))
image(microarray.raw.data[,6], col=rainbow(100))
image(microarray.raw.data[,7], col=rainbow(100))
image(microarray.raw.data[,8], col=rainbow(100))
image(microarray.raw.data[,9], col=rainbow(100))
image(microarray.raw.data[,10], col=rainbow(100))
image(microarray.raw.data[,11], col=rainbow(100))
image(microarray.raw.data[,12], col=rainbow(100))
image(microarray.raw.data[,13], col=rainbow(100))
image(microarray.raw.data[,14], col=rainbow(100))
image(microarray.raw.data[,15], col=rainbow(100))
image(microarray.raw.data[,16], col=rainbow(100))
image(microarray.raw.data[,17], col=rainbow(100))
image(microarray.raw.data[,18], col=rainbow(100))
image(microarray.raw.data[,19], col=rainbow(100))
image(microarray.raw.data[,20], col=rainbow(100))
image(microarray.raw.data[,21], col=rainbow(100))
image(microarray.raw.data[,22], col=rainbow(100))
image(microarray.raw.data[,23], col=rainbow(100))
image(microarray.raw.data[,24], col=rainbow(100))

# quality resumen

arrayQualityMetrics(expressionset = microarray.raw.data,
                    outdir = "../dataArrayQualityMetrics_whAffy",
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = c("organ", "samples"))

apply(probes(microarray.raw.data,"pm"),2,min) 

# 2. Preprocesed data
## Robust Multiarray Average (RMA)
png("../boxplot_Fluorescence_Raw.png", width = 600, height = 600)
boxplot(microarray.raw.data, col = rainbow(24), las=2, ylab="Fluorescence")
dev.off()

png("../hist_Fliorencence_Raw.png", width = 600, height = 600)
hist(microarray.raw.data, col= rainbow(24))
dev.off()

# remove background noise
microarray.processed.data <-rma(microarray.raw.data) #take the log2 as expression data is commonly analyzed on a logaritmic scale

png("../boxplot_Fluorescence_Processed.png", width = 600, height = 600)
boxplot(microarray.processed.data, col = rainbow(24), las=2, ylab="Fluorescence A.U.")
dev.off()

png("../hist_Fliorencence_Processed.png", width = 600, height = 600)
hist(microarray.processed.data, col= rainbow(24))
dev.off()

expression.level <- exprs(microarray.processed.data)
dim(expression.level)

expression.level <- read.csv('../MajorRevision_Enero2024/symbol_expressionLevel_rmaDta.csv', row.names =1)#Enero 2024

sampleID <- c("normal_prostate_1","tumor_prostate_1","normal_prostate_2","tumor_prostate_2",
              "normal_prostate_3","tumor_prostate_3","normal_prostate_4","tumor_prostate_4",
              "normal_prostate_5","tumor_prostate_5","normal_prostate_6","tumor_prostate_6",
              "normal_breast_1","tumor_breast_1","normal_breast_2","tumor_breast_2",
              "normal_breast_3","tumor_breast_3","normal_breast_4","tumor_breast_4",
              "normal_breast_5","tumor_breast_5","normal_breast_6","tumor_breast_6") 
colnames(expression.level) <- sampleID
head(expression.level)

normal.prostate <- (expression.level[,"normal_prostate_1"] + expression.level[,"normal_prostate_2"] + 
                      expression.level[,"normal_prostate_3"] + expression.level[,"normal_prostate_4"] +
                      expression.level[,"normal_prostate_5"] + expression.level[,"normal_prostate_6"])/6
tumor.prostate <- (expression.level[,"tumor_prostate_1"] + expression.level[,"tumor_prostate_2"] + 
                     expression.level[,"tumor_prostate_3"] + expression.level[,"tumor_prostate_4"] +
                     expression.level[,"tumor_prostate_5"] + expression.level[,"tumor_prostate_6"])/6
normal.breast <- (expression.level[,"normal_breast_1"] + expression.level[,"normal_breast_2"] + 
                    expression.level[,"normal_breast_3"] + expression.level[,"normal_breast_4"] +
                    expression.level[,"normal_breast_5"] + expression.level[,"normal_breast_6"])/6
tumor.breast <- (expression.level[,"tumor_breast_1"] + expression.level[,"tumor_breast_2"] + 
                   expression.level[,"tumor_breast_3"] + expression.level[,"tumor_breast_4"] +
                   expression.level[,"tumor_breast_5"] + expression.level[,"tumor_breast_6"])/6


mean.expression <- matrix(c(normal.prostate, tumor.prostate, normal.breast, tumor.breast), ncol=4)
conditions.id <- c("normal_prostate", "tumor_prostate", "normal_breast", "tumor_breast")
rownames(mean.expression) <- names(normal.prostate)
colnames(mean.expression) <- conditions.id
head(mean.expression)

# 3. estimation of expression levels

# comparacion preliminar entre condiciones (tumor vs normal)
png("../PreliminarComparison_tumorProstateVSnormalProstate.png", width = 600, height = 600)
plot(tumor.prostate,normal.prostate, xlab="Normal stromal prostate", ylab= "Tumor stromal prostate", pch=19, cex=0.5)
dev.off()

png("../PreliminarComparison_tumorBreastVSnormalBreast.png", width = 600, height = 600)
plot(tumor.breast,normal.breast, xlab="Normal stromal breast", ylab= "Tumor stromal breast", pch=19, cex=0.5)
dev.off()

# 4. Selection of differentially expressed genes
# limma 
experimental.desing <- model.matrix(~ -1+factor(c(1,2,1,2,1,2,1,2,1,2,1,2,3,4,3,4,3,4,3,4,3,4,3,4)))
colnames(experimental.desing) <- c("normal_prostate", "tumor_prostate", "normal_breast", "tumor_breast")

linear.fit <- lmFit(expression.level, experimental.desing) 

contrast.matrix <-makeContrasts(tumor_prostate-normal_prostate,
                                tumor_breast-normal_breast, 
                                levels = c("normal_prostate", "tumor_prostate", "normal_breast", "tumor_breast"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix) # Fold-Change
contrast.results <- eBayes(contrast.linear.fit) # p-valor/q-valor

nrow(expression.level)
Prostate.results <- topTable(contrast.results, number = 20824, coef = 1)
head(Prostate.results)

Breast.results <- topTable(contrast.results, number = 20824, coef = 2)
head(Breast.results)

# DEGs -- prostate

fold.change.Prostate <- Prostate.results$logFC #FC
adj.Pval.Prostate <- Prostate.results$adj.P.Val # pvalor ajusted
genes.ids.Prostate <- rownames(Prostate.results)

activated.genes.prostate.1 <- genes.ids.Prostate[fold.change.Prostate > 1 & adj.Pval.Prostate < 0.10]
repressed.genes.prostate.1 <- genes.ids.Prostate[fold.change.Prostate < -1 & adj.Pval.Prostate < 0.10]

length(activated.genes.prostate.1)
length(repressed.genes.prostate.1)

# DEGs -- breast

fold.change.Breast <- Breast.results$logFC # FC
adj.Pval.Breast <- Breast.results$adj.P.Val # pvalor ajusted
genes.ids.Breast <- rownames(Breast.results)

# umbral tipico de 2
activated.genes.breast.1 <- genes.ids.Breast[fold.change.Breast > 1 & adj.Pval.Breast < 0.05]
repressed.genes.breast.1 <- genes.ids.Breast[fold.change.Breast < -1 & adj.Pval.Breast < 0.05]

length(activated.genes.breast.1)
length(repressed.genes.breast.1)

# Visualizacion en grafico -- si solo se filtra por fold-change

png("../Comparison_tumorProstateVSnormalProstate_DEG_FC.png", width = 600, height = 600)
plot(tumor.prostate,normal.prostate, xlab="Normal stromal breast", ylab= "Tumor stromal breast", pch=19, cex=0.5,
     col = "grey")
points(tumor.prostate[activated.genes.prostate.1], normal.prostate[activated.genes.prostate.1], pch = 19, cex = 0.5,
      col = "red")
points(tumor.prostate[repressed.genes.prostate.1], normal.prostate[repressed.genes.prostate.1], pch = 19, cex = 0.5,
       col = "blue")
dev.off()

png("../Comparison_tumorBreastVSnormalBreast_DEG_FC.png", width = 600, height = 600)
plot(tumor.breast,normal.breast, xlab="Normal stromal breast", ylab= "Tumor stromal breast", pch=19, cex=0.5,
     col="grey")
points(tumor.breast[activated.genes.breast.1], normal.breast[activated.genes.breast.1], pch = 19, cex = 0.5,
       col = "red")
points(tumor.breast[repressed.genes.breast.1], normal.breast[repressed.genes.breast.1], pch = 19, cex = 0.5,
       col = "blue")
dev.off()

# volvano plot
# Prostate
names(fold.change.Prostate) <- genes.ids.Prostate
log.padj.prostate <- -log10(adj.Pval.Prostate)
names(log.padj.prostate) <- genes.ids.Prostate

png("../Comparison_tumorProstateVSnormalProstate_DEG_volvanoPlot_010.png", width = 600, height = 600)
plot(fold.change.Prostate,log.padj.prostate, ylab="-log10(p value)", xlab= "log2 fold change", pch=19, cex=0.5,
     col = "grey", xlim=c(-6,6))
points(fold.change.Prostate[activated.genes.prostate.1],log.padj.prostate[activated.genes.prostate.1], 
       pch = 19, cex = 0.5, col = "red")
points(fold.change.Prostate[repressed.genes.prostate.1],log.padj.prostate[repressed.genes.prostate.1], 
       pch = 19, cex = 0.5, col = "blue")
dev.off()

# Breast
names(fold.change.Breast) <- genes.ids.Breast
log.padj.breast <- -log10(adj.Pval.Breast)
names(log.padj.breast) <- genes.ids.Breast

png("../Comparison_tumorBreastVSnormalBreast_DEG_volvanoPlot.png", width = 600, height = 600)
plot(fold.change.Breast,log.padj.breast, ylab="-log10(p value)", xlab= "log2 fold change", pch=19, cex=0.5,
     col = "grey", xlim=c(-6,6))
points(fold.change.Breast[activated.genes.breast.1],log.padj.breast[activated.genes.breast.1], 
       pch = 19, cex = 0.5, col = "red")
points(fold.change.Breast[repressed.genes.breast.1],log.padj.breast[repressed.genes.breast.1], 
       pch = 19, cex = 0.5, col = "blue")
dev.off()

# results
BiocManager::install("hgu133plus2.db", force = TRUE)
library("hgu133plus2.db")

totalGenes <- aafTableAnn(rownames(expression.level), "hgu133plus2.db", aaf.handler())
saveText(totalGenes, file="../MajorRevision_Enero2024/aafTableAnn_majorRevision.txt")
write.table(expression.level, file="../MajorRevision_Enero2024/rawAffy_expressionLevel_rmaData.csv", sep = ",", 
            quote = FALSE, row.names = TRUE, col.names = TRUE)

activated.genes.prostate.table <- aafTableAnn(activated.genes.prostate.1, "hgu133plus2.db", aaf.handler())

# guardar datos log2 completos

write.table(exprs(microarray.processed.data), file="../data_log2_prostate_010-breast.txt", sep = "\t", 
            quote = FALSE, row.names = TRUE)

# datasets -- entry co-expression network
prostate.all.DEG <- genes.ids.Prostate[fold.change.Prostate > 1 & adj.Pval.Prostate < 0.10 | fold.change.Prostate < -1 & adj.Pval.Prostate < 0.10]
breast.all.DEG <- genes.ids.Breast[fold.change.Breast > 1 & adj.Pval.Breast < 0.05 | fold.change.Breast < -1 & adj.Pval.Breast < 0.05]

normal.prostate.DEG.table <- expression.level[, c("normal_prostate_1", "normal_prostate_2", 
                                                  "normal_prostate_3", "normal_prostate_4", 
                                                  "normal_prostate_5", "normal_prostate_6")]
normal.prostate.DEG.table <- normal.prostate.DEG.table[rownames(normal.prostate.DEG.table) %in% prostate.all.DEG,]
normal.prostate.DEG.table <- cbind(attr_name = rownames(normal.prostate.DEG.table), normal.prostate.DEG.table)

tumor.prostate.DEG.table <- expression.level[, c("tumor_prostate_1", "tumor_prostate_2", 
                                                 "tumor_prostate_3", "tumor_prostate_4", 
                                                 "tumor_prostate_5", "tumor_prostate_6")]
tumor.prostate.DEG.table <- tumor.prostate.DEG.table[rownames(tumor.prostate.DEG.table) %in% prostate.all.DEG,]
tumor.prostate.DEG.table <- cbind(attr_name = rownames(tumor.prostate.DEG.table), tumor.prostate.DEG.table)

normal.breast.DEG.table <- expression.level[, c("normal_breast_1", "normal_breast_2", 
                                                "normal_breast_3", "normal_breast_4", 
                                                "normal_breast_5", "normal_breast_6")]
normal.breast.DEG.table <- normal.breast.DEG.table[rownames(normal.breast.DEG.table) %in% breast.all.DEG,]
normal.breast.DEG.table <- cbind(attr_name = rownames(normal.breast.DEG.table), normal.breast.DEG.table)

tumor.breast.DEG.table <- expression.level[, c("tumor_breast_1", "tumor_breast_2", 
                                               "tumor_breast_3", "tumor_breast_4", 
                                               "tumor_breast_5", "tumor_breast_6")]
tumor.breast.DEG.table <- tumor.breast.DEG.table[rownames(tumor.breast.DEG.table) %in% breast.all.DEG,]
tumor.breast.DEG.table <- cbind(attr_name = rownames(tumor.breast.DEG.table), tumor.breast.DEG.table)

write.table(normal.prostate.DEG.table, file="../MajorRevision_Enero2024/normalProstate_DEGs_table_MajorRevision.csv", sep = ",", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(tumor.prostate.DEG.table, file="../MajorRevision_Enero2024/tumorProstate_DEGs_table_MajorRevision.csv", sep = ",", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(normal.breast.DEG.table, file="../MajorRevision_Enero2024/normalBreast_DEGs_table_MajorRevision.csv", sep = ",", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(tumor.breast.DEG.table, file="../MajorRevision_Enero2024/tumorBreast_DEGs_table_MajorRevision.csv", sep = ",", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(repressed.genes.breast.1, file="../MajorRevision_Enero2024/repressedBreast_DEGs_table_MajorRevision.txt", sep = ",", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
