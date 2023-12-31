# Analisis de enriquecimiento funcional
BiocManager::install("clusterProfiler")
BiocManager::install("hgu133plus2.db")
BiocManager::install("DO.db")

library(BiocManager)
library(hgu133plus2.db)
library(clusterProfiler)
library(DO.db)
library(ggplot2)

carpeta <- "hubs_prostate_netwok//"
archivos <- list.files(carpeta)

archivos_txt <- archivos[grep("\\.txt$", archivos)]

# Analisis para GO
for (archivo in archivos_txt) {
  genelist <- readLines(paste0(carpeta,archivo))
  entrez_ids <- mapIds(hgu133plus2.db, keys=genelist, column='ENTREZID', keytype='SYMBOL')
  go <- enrichGO(entrez_ids,  OrgDb = "hgu133plus2.db", ont = "all", pAdjustMethod = "fdr")
  if(length(go@result$ONTOLOGY) != 0){
    write.csv(go, paste0(carpeta,archivo,"all_BH.csv"), row.names = FALSE)
    grafico <- dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    ggsave(filename = paste0(carpeta,archivo,"all_BH.svg"), plot = grafico, device = "svg")
  }
}

# Analisis para KEGG
for (archivo in archivos_txt) {
  genelist <- readLines(paste0(carpeta,archivo))
  entrez_ids <- mapIds(hgu133plus2.db, keys=genelist, column='ENTREZID', keytype='SYMBOL')
  kegg <- enrichKEGG(entrez_ids, pAdjustMethod = "fdr")
  if(!is.null(kegg)){
    if(any(kegg@result$p.adjust <0.05)){
      write.csv(kegg, paste0(carpeta,archivo,"KEGG_BH.csv"), row.names = FALSE)
      grafico <- dotplot(kegg) + facet_grid(scale="free")
      ggsave(filename = paste0(carpeta,archivo,"KEGG_BH.svg"), plot = grafico, device = "svg")
    }
  }
}

 
