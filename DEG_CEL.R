## Preprocesamiento 2 .CEL

library("BiocManager")
library("affy")
library("gcrma")
library("affyPLM")
library("limma")
library("annaffy")
library("arrayQualityMetrics")

# Cargar los ficheros .CEL (datos brutos)
SDRF <- read.delim("../data_info.csv", sep = ",")
rownames(SDRF) <- SDRF$Ids
SDRF <- AnnotatedDataFrame(SDRF)

microarray.raw.data <- ReadAffy(filenames = SDRF$Array.Data.File, verbose=TRUE, phenoData = SDRF)
microarray.raw.data
cdfName(microarray.raw.data) # extrae el diseno de placa
stopifnot(validObject(microarray.raw.data))

# 1. Analisis de calidad
# deteccion de dano fisico en el microarray - imagen que detecto el escaner
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

# resumen de la calidad

arrayQualityMetrics(expressionset = microarray.raw.data,
                    outdir = "../dataArrayQualityMetrics_whAffy",
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = c("organ", "samples"))

apply(probes(microarray.raw.data,"pm"),2,min) #minimos de PM en las distintas muestras
# Y vemos como efectivamente no son nulos. La corrección de fondo
# pretende que aproximadamente sean nulos los valores más pequeños.
#Esto correspondería con aquellos genes que no tienen actividad

# 2. Preprocesamiento de los datos
## Robust Multiarray Average (RMA)
png("../boxplot_Fluorescence_Raw.png", width = 600, height = 600)
boxplot(microarray.raw.data, col = rainbow(24), las=2, ylab="Fluorescence")
dev.off()

png("../hist_Fliorencence_Raw.png", width = 600, height = 600)
hist(microarray.raw.data, col= rainbow(24))
dev.off()

# eliminar ruido de fondo
microarray.processed.data <-rma(microarray.raw.data) #take the log2 as expression data is commonly analyzed on a logaritmic scale
# estos datos estan transformados en log2

png("../boxplot_Fluorescence_Processed.png", width = 600, height = 600)
boxplot(microarray.processed.data, col = rainbow(24), las=2, ylab="Fluorescence A.U.")
dev.off()

png("../hist_Fliorencence_Processed.png", width = 600, height = 600)
hist(microarray.processed.data, col= rainbow(24))
dev.off()

# las muestras ya son comparables entre si

expression.level <- exprs(microarray.processed.data)
dim(expression.level)

sampleID <- c("normal_prostate_1","tumor_prostate_1","normal_prostate_2","tumor_prostate_2",
              "normal_prostate_3","tumor_prostate_3","normal_prostate_4","tumor_prostate_4",
              "normal_prostate_5","tumor_prostate_5","normal_prostate_6","tumor_prostate_6",
              "normal_breast_1","tumor_breast_1","normal_breast_2","tumor_breast_2",
              "normal_breast_3","tumor_breast_3","normal_breast_4","tumor_breast_4",
              "normal_breast_5","tumor_breast_5","normal_breast_6","tumor_breast_6") # renombrar las columnas
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

# 3. Estimacion de los niveles de expresion

# comparacion preliminar entre condiciones (tumor vs normal)
png("../PreliminarComparison_tumorProstateVSnormalProstate.png", width = 600, height = 600)
plot(tumor.prostate,normal.prostate, xlab="Normal stromal prostate", ylab= "Tumor stromal prostate", pch=19, cex=0.5)
dev.off()

png("../PreliminarComparison_tumorBreastVSnormalBreast.png", width = 600, height = 600)
plot(tumor.breast,normal.breast, xlab="Normal stromal breast", ylab= "Tumor stromal breast", pch=19, cex=0.5)
dev.off()

# 4. Seleccion de genes diferencialemnte expresados
# limma - Seleccion de DEGs basada en fold-change(FC) e inferencia estadistica
experimental.desing <- model.matrix(~ -1+factor(c(1,2,1,2,1,2,1,2,1,2,1,2,3,4,3,4,3,4,3,4,3,4,3,4)))
colnames(experimental.desing) <- c("normal_prostate", "tumor_prostate", "normal_breast", "tumor_breast")

linear.fit <- lmFit(expression.level, experimental.desing) # calcula la media de las columnas

contrast.matrix <-makeContrasts(normal_prostate-tumor_prostate,
                                normal_breast-tumor_breast, 
                                levels = c("normal_prostate", "tumor_prostate", "normal_breast", "tumor_breast"))  # especificacion de contrastes

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix) # calculo del Fold-Change
contrast.results <- eBayes(contrast.linear.fit) # calculo p-valor/q-valor

# extraer informacion por comparaciones
nrow(expression.level)
Prostate.results <- topTable(contrast.results, number = 54675, coef = 1)
head(Prostate.results)

Breast.results <- topTable(contrast.results, number = 54675, coef = 2)
head(Breast.results)

# extraer DEGs -- prostate

fold.change.Prostate <- Prostate.results$logFC # para obtener FC
adj.Pval.Prostate <- Prostate.results$adj.P.Val # pvalor ajustado
genes.ids.Prostate <- rownames(Prostate.results)

# umbral tipico de 2
activated.genes.prostate.1 <- genes.ids.Prostate[fold.change.Prostate > 1 & adj.Pval.Prostate < 0.10]
repressed.genes.prostate.1 <- genes.ids.Prostate[fold.change.Prostate < -1 & adj.Pval.Prostate < 0.10]

length(activated.genes.prostate.1)
length(repressed.genes.prostate.1)

# extraer DEGs -- breast

fold.change.Breast <- Breast.results$logFC # para obtener FC
adj.Pval.Breast <- Breast.results$adj.P.Val # pvalor ajustado
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
#text(tumor.prostate["213816_s_at"]+0.3, normal.prostate["213816_s_at"]+0.3, "IRT1", col = "black", cex = 0.7) # para resaltar alguna etiqueta de genes
dev.off()

png("../Comparison_tumorBreastVSnormalBreast_DEG_FC.png", width = 600, height = 600)
plot(tumor.breast,normal.breast, xlab="Normal stromal breast", ylab= "Tumor stromal breast", pch=19, cex=0.5,
     col="grey")
points(tumor.breast[activated.genes.breast.1], normal.breast[activated.genes.breast.1], pch = 19, cex = 0.5,
       col = "red")
points(tumor.breast[repressed.genes.breast.1], normal.breast[repressed.genes.breast.1], pch = 19, cex = 0.5,
       col = "blue")
dev.off()

# Visualizacion en grafico -- si se filtra por fold-change y padj -- volvano plot
# Prostata
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
#text(tumor.prostate["213816_s_at"]+0.3, normal.prostate["213816_s_at"]+0.3, "IRT1", col = "black", cex = 0.7) # para resaltar alguna etiqueta de genes
dev.off()

# Pecho
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
#text(tumor.prostate["213816_s_at"]+0.3, normal.prostate["213816_s_at"]+0.3, "IRT1", col = "black", cex = 0.7) # para resaltar alguna etiqueta de genes
dev.off()

# guardar resultados
BiocManager::install("hgu133plus2.db", force = TRUE)
library("hgu133plus2.db")

activated.genes.prostate.table <- aafTableAnn(activated.genes.prostate.1, "hgu133plus2.db", aaf.handler())

saveHTML(activated.genes.prostate.table, file="../activatedGenes_prostate_table_010.html")
saveText(activated.genes.prostate.table, file="../activatedGenes_prostate_table_010.txt")

repressed.genes.prostate.table <- aafTableAnn(repressed.genes.prostate.1, "hgu133plus2.db", aaf.handler())

saveHTML(repressed.genes.prostate.table, file="../repressedGenes_prostate_table_010.html")
saveText(repressed.genes.prostate.table, file="../repressedGenes_prostate_table_010.txt")

activated.genes.breast.table <- aafTableAnn(activated.genes.breast.1, "hgu133plus2.db", aaf.handler())

saveHTML(activated.genes.breast.table, file="../activatedGenes_breast_table.html")
saveText(activated.genes.breast.table, file="../activatedGenes_breast_table.txt")

repressed.genes.breast.table <- aafTableAnn(repressed.genes.breast.1, "hgu133plus2.db", aaf.handler())

saveHTML(repressed.genes.breast.table, file="../repressedGenes_breast_table.html")
saveText(repressed.genes.breast.table, file="../repressedGenes_breast_table.txt")

# guardar datos log2 completos

write.table(exprs(microarray.processed.data), file="../data_log2_prostate_010-breast.txt", sep = "\t", 
            quote = FALSE, row.names = TRUE)

# praprar datasets por separado -- entrada en pyEnGNet
prostate.all.DEG <- union(activated.genes.prostate.1, repressed.genes.prostate.1)
breast.all.DEG <- union(activated.genes.breast.1, repressed.genes.breast.1)

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

write.table(normal.prostate.DEG.table, file="../normalProstate_DEGs_table_010.csv", sep = ",", 
            quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(tumor.prostate.DEG.table, file="../tumorProstate_DEGs_table_010.csv", sep = ",", 
            quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(normal.breast.DEG.table, file="../normalBreast_DEGs_table.csv", sep = ",", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(tumor.breast.DEG.table, file="../tumorBreast_DEGs_table.csv", sep = ",", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)


