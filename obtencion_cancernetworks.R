# Obtener red de cáncer
library(dplyr)
setwd('/Users/juliamf/Desktop/4BTG/tfg/ART/Cancer networks') #guardar output
# CANCER DE PROSTATA

# leer ficheros de las redes - umbral = 0.85
tumor_prostate_0.85 <- read.csv('/Users/juliamf/Desktop/4BTG/tfg/ART/DEGs_ex_data/symbol/prostate/output_majorVoting/DEGs_prostate_tumor_symbol_threshold_0.85.csv', header = TRUE)
normal_prostate_0.85 <- read.csv('/Users/juliamf/Desktop/4BTG/tfg/ART/DEGs_ex_data/symbol/prostate/output_majorVoting/DEGs_prostate_normal_symbol_threshold_0.85.csv',header = TRUE)

# red cancer - red normal = nos quedamos solo con las interacciones pertenecientes a la red de cancer 

prostate_cancer_network <- anti_join(tumor_prostate_0.85, normal_prostate_0.85, by = c('Origen', 'Destino'))
write.csv(prostate_cancer_network, "prostate_cancer_network.csv", row.names = FALSE)

tumor_prostate_0.75 <- read.csv('/Users/juliamf/Desktop/4BTG/tfg/ART/DEGs_ex_data/symbol/prostate/output_majorVoting/DEGs_prostate_tumor_symbol_threshold_0.75.csv', header = TRUE)
normal_prostate_0.75 <- read.csv('/Users/juliamf/Desktop/4BTG/tfg/ART/DEGs_ex_data/symbol/prostate/output_majorVoting/DEGs_prostate_normal_symbol_threshold_0.75.csv',header = TRUE)
prostate_cancer_network_075 <- anti_join(tumor_prostate_0.75, normal_prostate_0.75, by = c('Origen', 'Destino'))
write.csv(prostate_cancer_network_075, "prostate_cancer_network_075.csv", row.names = FALSE)

tumor_prostate_0.7 <- read.csv('/Users/juliamf/Desktop/4BTG/tfg/ART/DEGs_ex_data/symbol/prostate/output_majorVoting/DEGs_prostate_tumor_symbol_threshold_0.7.csv', header = TRUE)
normal_prostate_0.7 <- read.csv('/Users/juliamf/Desktop/4BTG/tfg/ART/DEGs_ex_data/symbol/prostate/output_majorVoting/DEGs_prostate_normal_symbol_threshold_0.7.csv',header = TRUE)
prostate_cancer_network_07 <- anti_join(tumor_prostate_0.7, normal_prostate_0.7, by = c('Origen', 'Destino'))
write.csv(prostate_cancer_network_07, "prostate_cancer_network_07.csv", row.names = FALSE)

# CANCER DE MAMA

# leer ficheros de las redes - umbral = 0.85 
tumor_breast_0.85 <- read.csv('/Users/juliamf/Desktop/4BTG/tfg/ART/DEGs_ex_data/symbol/breast/output_majorVoting/DEGs_breast_tumor_symbol_threshold_0.85.csv', header = TRUE)
normal_breast_0.85 <- read.csv('/Users/juliamf/Desktop/4BTG/tfg/ART/DEGs_ex_data/symbol/breast/output_majorVoting/DEGs_breast_nromal_symbol_threshold_0.85.csv',header = TRUE)

# red cancer - red normal = nos quedamos solo con las interacciones pertenecientes a la red de cancer 

breast_cancer_network <- anti_join(tumor_breast_0.85, normal_breast_0.85, by = c('Origen', 'Destino'))
write.csv(breast_cancer_network, "breast_cancer_network.csv", row.names = FALSE)


