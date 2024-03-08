##### command lines to select the human samples, 2nd version, only for BEAT selection
# done by Jonathan Seguin, group of Prof. Schwaller, DBM, UKBB, Basel, Switzerland
# email: jonathan.seguin@unibas.ch, seguin.jonathan@gmail.com
# Wed Jul 26 09:54:09 2023


# prepare the working environment -----------------------------------------------------

# # update the libraries to prevent bugs ----------------------------------
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data", force = T)
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks", force = T)

## load the librairies -----------------------------------------------------
library(xlsx)
library("TCGAbiolinks")
library(SummarizedExperiment)
library(tidyverse)
library(edgeR)
library(biomaRt)
library(ggpubr)
library(GenomicDataCommons)
library(factoextra)
library(plotly)
library(readxl)
library(pheatmap)
library(PCAtools)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(writexl)
library(rgl)
library(plotly)
library(gridExtra)
library(cluster)
library(survival)
library(ggsurvfit)

## load the functions -----------------------------------------------------
source("Src/DE_functions.R")


##### function to create scatterplot
#  data, data.frame object, contains genes expression
#  x, character value, indicate the name of the column to put values on x-axis
#  y, character value, indicate the name of the column to put values on y-axis
# log2, boolean value, indicate if the log2 values must be used (FALSE by default)
#### 
getScatterplot <- function(data, x, y, log2 = F){
  
  if(log2){
    scatterplot =  ggplot(data = data, aes(x = log2(data[, x]), y = log2(data[, y]), color = Color)) + 
      geom_point() +
      ylab(paste0(y, " (log2(tpm))")) + 
      xlab(paste0(x, " (log2(tpm))")) +
      theme_bw() + 
      ggtitle(paste0("scatterplot, ", x, " vs ", y, ", BEAT database")) 
  } else {
    scatterplot =  ggplot(data = data, aes(x = data[, x], y = data[, y], color = Color)) + 
      geom_point() +
      ylab(paste0(y, " (tpm)")) + 
      xlab(paste0(x, " (tpm)")) + 
      theme_bw() + 
      ggtitle(paste0("scatterplot, ", x, " vs ", y, ", BEAT database")) 
  }
  
  
  # return the scatterplot
  return(scatterplot)
  
}



## read the tables ---------------------------------------------------------

# read the TCGA file
TCGA_table <- read.xlsx2(file = "Dataset/BEAT/combined_study_clinical_data_TCGA-AML.xlsx",
                         sheetIndex = 1)

# read the patient information from cbioportal
BEAT_cbioportal <- read.csv("Dataset/BEAT/combined_study_clinical_data.tsv", header = T, sep = "\t")


# load the table from the TPO-Evi1 RNAseq
fits_tmp <- readRDS(file = "DEA/Rds/RDS_woD/edge_samplesR_fits_results_DEA.rds")

# select the orthologous genes
human_ortho <- getHumanOrthologous()

# select the mouse genes
Mouse_genes <- getMouseGene()


# select only patients from GDC data
GDC_patients_info <- read.table(file = "Dataset/BEAT/clinical.cases_selection.2024-02-09/clinical.tsv", header = T, sep = "\t", quote = "")
GDC_samples_info <- read.table(file = "Dataset/BEAT/biospecimen.cases_selection.2024-02-09/sample.tsv", header = T, sep = "\t", quote = "")
tissue_norm <- GDC_samples_info %>% dplyr::filter(sample_type == "Solid Tissue Normal")


# load the table from the bulk TPO vs PBS
fit_bulkTPOVsnoTPO <- as_tibble(as.data.frame(topTags(fits_tmp$bulk.TPO, 
                                                      n = nrow(fits_tmp$bulk.TPO))))

# add the mouse gene information to the DEA result
fit_bulkTPOVsnoTPO <- merge(x = fit_bulkTPOVsnoTPO, 
                            y = Mouse_genes, 
                            by.x = "SYMBOL", 
                            by.y = "external_gene_name")

# merge with the mouse
fit_bulkTPOVsnoTPO <- merge(x = fit_bulkTPOVsnoTPO, 
                            y = human_ortho, 
                            by.x = "ensembl_gene_id", 
                            by.y = "mmusculus_homolog_ensembl_gene")


# load the table from the bulk TPO vs PBS
fit_GFP_TPOVsnoTPO <- as_tibble(as.data.frame(topTags(fits_tmp$GFP.TPO, 
                                                      n = nrow(fits_tmp$GFP.TPO))))

# add the mouse gene information to the DEA result
fit_GFP_TPOVsnoTPO <- merge(x = fit_GFP_TPOVsnoTPO, 
                            y = Mouse_genes, 
                            by.x = "SYMBOL", 
                            by.y = "external_gene_name")

# merge with the mouse
fit_GFP_TPOVsnoTPO <- merge(x = fit_GFP_TPOVsnoTPO, 
                            y = human_ortho, 
                            by.x = "ensembl_gene_id", 
                            by.y = "mmusculus_homolog_ensembl_gene")



# directories output ------------------------------------------------------
outdir_db <- create_dir(path = "DEA/Human_patients/BEAT/")


# download gene information ------------------------------------------------------

# download ensembl information for human genes
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
                   host = "sep2019.archive.ensembl.org")

# include patients from BEAT database ---------------------
# source (https://www.bioconductor.org/packages/devel/bioc/vignettes/GenomicDataCommons/inst/doc/overview.html#1_What_is_the_GDC)

# source: https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html
projects <- getGDCprojects()$project_id

#Downloading and prepare TCGA CASE
query_BEAT <- GDCquery(project = "BEATAML1.0-COHORT",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")
samplesDown_BEAT <- getResults(query_BEAT, cols = c("cases"))

# download clinical information
#query_clin_BEAT <-  GDCquery_clinic(project = "BEATAML1.0-COHORT", type = "clinical")
BEAT_clinical <- GDCquery_clinic(project = "BEATAML1.0-COHORT", type = "clinical")


# download the patients
GDCdownload(query_BEAT)
GDCdownload(BEAT_clinical)
GDCdownload(query_TCGA_suppClin_BCR)


# get the manifest for the BEAT query
getManifest(query = query_BEAT, save = T)

# load the manifest
#manifest <-  read.table("MANIFEST_old.txt", header = T)
manifest <-  read.table("gdc_manifest.txt", header = T)

# save the information about the patients
write_xlsx(BEAT_clinical, path = "DEA/Human_patients/BEAT/patients_informations.xlsx")

# select the patients who suffer only of myeloid leukemia
BEAT_clinical_ML <- BEAT_clinical[grep(pattern = "Myeloid Leukemias", BEAT_clinical$disease_type),]

# save information about the selected patients
write_xlsx(BEAT_clinical_ML, path = file.path(outdir_db, "selected_patient_BEAT_allMLL.xlsx"))
gc()


# prepare the samples
data_BEAT <- GDCprepare(query_BEAT)

# save in the rds file
saveRDS(data_BEAT, file = "DEA/Human_patients/BEAT/data_BEAT.rds")
saveRDS(BEAT_clinical, file = "DEA/Human_patients/BEAT/BEAT_clinical.rds")

# selected samples
all_samples_ML <- unique(unlist(sapply(BEAT_clinical_ML$submitter_sample_ids, function(x)unlist(strsplit(x, split = ",")))))

# remove patients with modification on t(3) or inv(3)
BEAT_clinical_ML_woEVI1 = BEAT_clinical_ML[intersect(grep(pattern = "t\\(3", BEAT_clinical_ML$primary_diagnosis, invert = T), 
                                              grep(pattern = "inv\\(3", BEAT_clinical_ML$primary_diagnosis, invert = T)), ]
BEAT_clinical_ML_EVI1 = BEAT_clinical_ML[grep(pattern = "t\\(3", BEAT_clinical_ML$primary_diagnosis), ]

# save information about the selected patients
write_xlsx(BEAT_clinical_ML, path = file.path(outdir_db, "selected_patient_BEAT_selected.xlsx"))
gc()

# select only AML patients
BEAT_clinical_allAML = BEAT_clinical_ML[grep(pattern = "Acute myeloid leukemia", BEAT_clinical_ML$primary_diagnosis), ]
BEAT_clinical_allAML_woEVI1 = BEAT_clinical_ML_woEVI1[grep(pattern = "Acute myeloid leukemia", BEAT_clinical_ML_woEVI1$primary_diagnosis), ]
BEAT_clinical_allAML_EVI1 = BEAT_clinical_ML_EVI1[grep(pattern = "Acute myeloid leukemia", BEAT_clinical_ML_EVI1$primary_diagnosis), ]


# save information about the selected patients
write_xlsx(BEAT_clinical_allAML_woEVI1, path = file.path(outdir_db, "selected_patient_BEAT_allAML_woEVI1.xlsx"))
gc()
write_xlsx(BEAT_clinical_allAML_EVI1, path = file.path(outdir_db, "selected_patient_BEAT_allAML_EVI1.xlsx"))
gc()
write_xlsx(BEAT_clinical_allAML, path = file.path(outdir_db, "selected_patient_BEAT_allAML.xlsx"))
gc()

# select only MLL patients
BEAT_clinical_ML_MLL <- BEAT_clinical_ML[grep(pattern = "MLL", BEAT_clinical_ML$primary_diagnosis), ]
BEAT_clinical_allAML_MLL <- BEAT_clinical_allAML[grep(pattern = "MLL", BEAT_clinical_allAML$primary_diagnosis), ]
BEAT_clinical_allAML_woEVI1_MLL <- BEAT_clinical_allAML_woEVI1[grep(pattern = "MLL", BEAT_clinical_allAML_woEVI1$primary_diagnosis), ]
  
# --> all are MLL-AF9 patients (12 patients)

# save information about the selected patients
write_xlsx(BEAT_clinical_ML_MLL, path = file.path(outdir_db, "selected_patient_BEAT_MLL-AF9.xlsx"))
gc()

# get the tpm value for all gene expression
table_BEAT_tpm <- assays(data_BEAT)$tpm_unstrand

# save the rds file
saveRDS(table_BEAT_tpm, file = "DEA/Human_patients/BEAT/table_BEAT_tpm.rds")

# draw the MECOM expression for all patients selected
all_samples_selected <- unique(unlist(sapply(BEAT_clinical_ML$submitter_sample_ids, function(x) strsplit(x, split = ",")))) 
id_patients <-  intersect(all_samples_selected, colnames(table_BEAT_tpm))

# draw the MECOM expression for all AML patients selected without EVI1 translocation/inversion
all_samples_selected_AML <- unique(unlist(sapply(BEAT_clinical_allAML$submitter_sample_ids, function(x) strsplit(x, split = ",")))) 
id_patients_AML <-  intersect(all_samples_selected_AML, colnames(table_BEAT_tpm))

## create a dataframe for the barplot of MECOM expression according to the patients
# select the identifient of the selected patients
id_samples_MLLAF9 <- unique(unlist(sapply(BEAT_clinical_ML_MLL$submitter_sample_ids, function(x) strsplit(x, split = ","))))
id_patients <- sapply(id_samples_MLLAF9, function(x) BEAT_clinical_ML_MLL$submitter_id[grep(pattern = x, BEAT_clinical_ML_MLL$submitter_sample_ids)])
dt_MLLAF9_BEAT <- data.frame("Patient" = id_patients, "Sample" = id_samples_MLLAF9)
# add the Mecom expression
common_samples <- intersect(colnames(table_BEAT_tpm), id_samples_MLLAF9)
dt_MLLAF9_BEAT = dt_MLLAF9_BEAT[dt_MLLAF9_BEAT$Sample %in% common_samples, ]
dt_MLLAF9_BEAT$MECOM = table_BEAT_tpm[mecom_id, dt_MLLAF9_BEAT$Sample]


# perform the differential expression between Evi1 high vs low for BEAT ---------------------

# provide the mecom Ensembl id
mecom_id <- "ENSG00000085276.19"

# extract the mecom expression
mecom_beat_log2tpm <- log2(assays(data_BEAT)$tpm_unstrand[mecom_id,] +1)

# select the rawcount of the dataset
table_beat_unstr <-  assays(data_BEAT)$unstranded

# download supplemental information from human genes
Human_genes <- getHumanGene()



# perform DEA based on k-means clustering result ---------------------

# create the DEA output dir
DEA_dir <- create_dir(path = "DEA/Human_patients/BEAT/allAML")

## find the cluster of samples ---------------------

# create the dataframe for the scatterplot figure
dt_BEAT_ERG_MECOM_log2tpm <- data.frame("MECOM" = log2(table_BEAT_tpm[mecom_id, ]),
                                          "ERG" = log2(table_BEAT_tpm[erg_id, ]),
                                          "Sample" = colnames(table_BEAT_tpm),
                                          "Color" = "none")

# check if there is normal tissues in the selected samples
intersect(dt_BEAT_ERG_MECOM_log2tpm$Sample, tissue_norm$sample_submitter_id) # 0

# check which categories of tissues contains the samples
table(GDC_samples_info$sample_type[GDC_samples_info$sample_submitter_id %in% dt_BEAT_ERG_MECOM_log2tpm$Sample])
# Primary Blood Derived Cancer - Bone Marrow   Primary Blood Derived Cancer - Peripheral Blood      Recurrent Blood Derived Cancer - Bone Marrow 
#  217                                               131                                               149 
# Recurrent Blood Derived Cancer - Peripheral Blood 
# 186 

# select only samples
dt_BEAT_ERG_MECOM_log2tpm = dt_BEAT_ERG_MECOM_log2tpm[dt_BEAT_ERG_MECOM_log2tpm$Sample %in% samples_id, ]


# update the patients information
dt_BEAT_ERG_MECOM_log2tpm$Color[dt_BEAT_ERG_MECOM_log2tpm$Sample %in% common_samples] = "MLL-AF9"
dt_BEAT_ERG_MECOM_log2tpm$Color[dt_BEAT_ERG_MECOM_log2tpm$Sample %in% all_samples_selected_allAML_EVI1] = "EVI1r"

# remove infinite values
dt_BEAT_ERG_MECOM_log2tpm = dt_BEAT_ERG_MECOM_log2tpm[is.finite(dt_BEAT_ERG_MECOM_log2tpm$MECOM), ]
dt_BEAT_ERG_MECOM_log2tpm = dt_BEAT_ERG_MECOM_log2tpm[is.finite(dt_BEAT_ERG_MECOM_log2tpm$ERG), ]


# initialize list to save kmeans results and their corresponding plots
list_kmeans <- NULL
list_kmean_plot <- NULL


# calculate the different k-means clusters
for(k in 2:9){
  
  # calculate the temporary clustering
  tmp_kmeans <-  kmeans(dt_BEAT_ERG_MECOM_log2tpm[, 1:2], centers = k, nstart = 25)
  
  # save the result in the list
  list_kmeans <- c(list_kmeans, list(tmp_kmeans))
  
  # save the plot
  list_kmean_plot <- c(list_kmean_plot, 
                       list(fviz_cluster(tmp_kmeans, data = dt_BEAT_ERG_MECOM_log2tpm[, 1:2], geom = "point", ggtheme = theme_bw()) + ggtitle(paste0("k = ", k))))
  
}

# renames the list
names(list_kmeans) <- paste0("k", 2:9)
names(list_kmean_plot) <- paste0("k", 2:9)


# save k-means object in rds file
saveRDS(list_kmeans, file = "DEA/Human_patients/BEAT/allAML/list_kmeans.rds")
saveRDS(list_kmean_plot, file = "DEA/Human_patients/BEAT/allAML/list_kmean_plot.rds")


# display the result
g <- grid.arrange(grobs = list_kmean_plot, ncol = 4, nrow = 2) 

# save in pdf file
ggsave(filename = "DEA/Human_patients/BEAT/allAML/kmeans_clustering_ERGvsMECOM_BEAT.pdf",
       plot = g, width =297, height = 210, units = "mm")

# determine the optimal number of cluster
set.seed(123)
gap_stat <- clusGap(dt_BEAT_ERG_MECOM_log2tpm[, 1:2], FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)

# save the result in the pdf
pdf(file = "DEA/Human_patients/BEAT/allAML/bestClusterNumber.pdf", paper = "a4r")
set.seed(123)
fviz_nbclust(dt_BEAT_ERG_MECOM_log2tpm[, 1:2], kmeans, method = "wss")
set.seed(123)
fviz_nbclust(dt_BEAT_ERG_MECOM_log2tpm[, 1:2], kmeans, method = "silhouette")
fviz_gap_stat(gap_stat)
dev.off()

# selection of cutoff 4 for the clustering result
# update the scatterplot 
dt_BEAT_ERG_MECOM_log2tpm$cluster = list_kmeans$k4$cluster


## perform the DEA based on clusters  (3 vs 4)---------------------

# get the samples information
cluster1_info <- dt_BEAT_ERG_MECOM_log2tpm[dt_BEAT_ERG_MECOM_log2tpm$cluster == "4", ]

# select the minimal cutoff for ERG
min_cutoff_ERG <- min(cluster1_info$ERG)

# filter the patient based on the ERG expression
patient_kept <- dt_BEAT_ERG_MECOM_log2tpm$Sample[dt_BEAT_ERG_MECOM_log2tpm$ERG >= min_cutoff_ERG]

# select patient high (cluster 4)
patient_high <- intersect(dt_BEAT_ERG_MECOM_log2tpm$Sample[dt_BEAT_ERG_MECOM_log2tpm$cluster == "3"], 
                          patient_kept)

# select patient low (cluster 1)
patient_low <- intersect(dt_BEAT_ERG_MECOM_log2tpm$Sample[dt_BEAT_ERG_MECOM_log2tpm$cluster == "4"], 
                         patient_kept)

# create output directory
tmp_DEA <- create_dir(path = file.path(DEA_dir, "BEAT_ERG_selection_kmeans"))

# select the samples
table_BEAT_unstr_sel <- table_beat_unstr[, colnames(table_beat_unstr) %in% c(patient_high, patient_low)]

# replace NA value in table
table_BEAT_unstr_sel[is.na(table_BEAT_unstr_sel)] = 0

# determine the group
group <- rep("low", ncol(table_BEAT_unstr_sel))

# update the group
group[colnames(table_BEAT_unstr_sel) %in% patient_high] = "high"

# remove the gene lowly expressed
print(paste0("# genes: ", nrow(table_BEAT_unstr_sel)))
# keep gene with more than 10 reads in total
count_10 <- nrow(table_BEAT_unstr_sel[rowSums(table_BEAT_unstr_sel) > 10,])
print(paste0("# kept genes (filter based on 10): ", count_10))
# 48062 kept genes among 60660 genes

# create the dge
dge <- DGEList(counts = table_BEAT_unstr_sel, group = group)

# filter according to the filtering gene expression
keep <- filterByExpr(dge)
count_filter <- length(keep[keep])
print(paste0("# kept gene after filtering by function: ", count_filter))
# 21179 genes kept

# filter the genes
dge <- dge[keep,,keep.lib.sizes=FALSE]

# calculate the normalization factor values
dge <- calcNormFactors(dge)

# create the design object
design <- model.matrix(~ 0 + group)

# change the colnames
colnames(design) <- c("high", "low")

# estimate the dispersion
dge <-  estimateDisp(dge, design = design)

# save the BCV plot
saveFigures("BCVplot_v2", ggplot = plotBCV(dge), dirPlot = tmp_DEA, A4 = T)

# perform the exact test
et <-  exactTest(dge)

# provide all results
allresult <- as.data.frame(topTags(et, n=nrow(et)))

# add ensembl gene id
allresult <-  cbind(allresult, "Ensembl_id" = sapply(rownames(allresult),
                                                     function(x) strsplit(x, split = "\\.")[[1]][1]))

# merge with the human gene information
allresult <- merge(x = allresult, 
                   y = unique(Human_genes[, 1:3]), 
                   by.x = "Ensembl_id", 
                   by.y = "ensembl_gene_id")

# correct the orientation of the FC
if(allresult$logFC[grepl(pattern = "^MECOM$", allresult$external_gene_name)] < 0){
  allresult$logFC = 0-allresult$logFC
}

# add the main result to the list
list_result <- list("Main result" = allresult)

# add to the list of genes for 5%
if(any(allresult$FDR <= 0.05)){
  list_result <- c(list_result, list("DEG, 5%" = allresult[allresult$FDR <= 0.05, ]))
}

# add to the list of genes for 1%
if(any(allresult$FDR <= 0.01)){
  list_result <- c(list_result, list("DEG, 1%" = allresult[allresult$FDR <= 0.01, ]))
}

# save the result in files
write_xlsx(x = list_result, path = file.path(tmp_DEA, "results_DEA_v2.xlsx"), col_names = T)
gc()

# save the selected samples
write.xlsx2(data.frame("Id_sample" = colnames(table_BEAT_unstr_sel), "Group" = group), file = file.path(tmp_DEA, "selected_patients.xlsx"),
            col.names = T, row.names = F)
gc()



### compare with the mouse DEG within a scatterplot ---------------------

# get the output path
outputPath <- tmp_DEA

# load the table with DEA result 
MLLr_BEAT_cutoff_ERG <- read_xlsx(path = file.path(outputPath, "results_DEA.xlsx"), 
                                    sheet = 1)

# merge the tables
merged_table <- merge(x = MLLr_BEAT_cutoff_ERG,
                      y = fit_bulkTPOVsnoTPO,
                      by = "external_gene_name",
                      suffixes = c("_BEAT", "_mouse"))

# create the color for the significativity 
merged_table$color_sign <- "n.s"

# update the color based on the significativity
merged_table$color_sign[merged_table$FDR_mouse <= 0.05] <- "sign. mouse"
merged_table$color_sign[merged_table$FDR_BEAT <= 0.05] <- "sign. BEAT"
merged_table$color_sign[(merged_table$FDR_mouse <= 0.05) & (merged_table$FDR_BEAT <= 0.05)] <- "sign. mouse and BEAT"

# calculate the coefficient correlation
fit <- lm(merged_table$logFC_BEAT ~ merged_table$logFC_mouse, data = merged_table)
coeff_corr <- round(cor(merged_table$logFC_BEAT, merged_table$logFC_mouse, method = "pearson"), 3)

# draw the scatterplot
scatterplot <- ggplot(merged_table, aes(x = logFC_BEAT, y = logFC_mouse, color = color_sign)) +
  geom_point(alpha = 0.75) + theme_bw() + xlab("log2FC BEAT") + ylab("log2FC mouse") +
  xlim(-max(abs(merged_table$logFC_BEAT)), max(abs(merged_table$logFC_BEAT))) +
  ylim(-max(abs(merged_table$logFC_mouse)), max(abs(merged_table$logFC_mouse))) +
  scale_color_manual(values = c("n.s" = "black", "sign. mouse" = "blue", 
                                "sign. BEAT" = "orange", "sign. mouse and BEAT" = "red")) +
  geom_text_repel(data = subset(merged_table[merged_table$external_gene_name %in% merged_table$external_gene_name[merged_table$color_sign == "sign. mouse and BEAT"],]), 
                  color = "black", aes(label = external_gene_name), size = 5) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between BEAT and bulk TPO vs PBS; corr coef (Pearson) = ", coeff_corr))


# save the scatterplot in files
saveFigures("scatterplot_bulk_TPOvsPBS", ggplot = scatterplot, dirPlot = outputPath, A4 = T)


# merge the tables
merged_table <- merge(x = MLLr_BEAT_cutoff_ERG,
                      y = fit_GFP_TPOVsnoTPO,
                      by = "external_gene_name",
                      suffixes = c("_BEAT", "_mouse"))


# create the color for the significativity 
merged_table$color_sign <- "n.s"

# update the color based on the significativity
merged_table$color_sign[merged_table$FDR_mouse <= 0.05] <- "sign. mouse"
merged_table$color_sign[merged_table$FDR_BEAT <= 0.05] <- "sign. BEAT"
merged_table$color_sign[(merged_table$FDR_mouse <= 0.05) & (merged_table$FDR_BEAT <= 0.05)] <- "sign. mouse and BEAT"

# calculate the coefficient correlation
fit <- lm(merged_table$logFC_BEAT ~ merged_table$logFC_mouse, data = merged_table)
coeff_corr <- round(cor(merged_table$logFC_BEAT, merged_table$logFC_mouse, method = "pearson"), 3)

# draw the scatterplot
scatterplot <- ggplot(merged_table, aes(x = logFC_BEAT, y = logFC_mouse, color = color_sign)) +
  geom_point(alpha = 0.75) + theme_bw() + xlab("log2FC BEAT") + ylab("log2FC mouse") +
  xlim(-max(abs(merged_table$logFC_BEAT)), max(abs(merged_table$logFC_BEAT))) +
  ylim(-max(abs(merged_table$logFC_mouse)), max(abs(merged_table$logFC_mouse))) +
  scale_color_manual(values = c("n.s" = "black", "sign. mouse" = "blue", 
                                "sign. BEAT" = "orange", "sign. mouse and BEAT" = "red")) +
  geom_text_repel(data = subset(merged_table[merged_table$external_gene_name %in% merged_table$external_gene_name[merged_table$color_sign == "sign. mouse and BEAT"],]), 
                  color = "black", aes(label = external_gene_name), size = 5) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between BEAT and GFP+ TPO vs PBS; corr coef (Pearson) = ", coeff_corr))


# save the scatterplot in files
saveFigures("scatterplot_GFP+_TPOvsPBS", ggplot = scatterplot, dirPlot = outputPath, A4 = T)



# add list of relevant genes
list_relevant_genes <- c("IL12RB2", "INPP4B", "MECOM", "ADGRG6", "MMP28", "SEMA4F",
                         "SH3BP5", "AIF1", "CDC14A", "FARP1", "GRAP2", "KRI1", 
                         "MLLT3", "MOCOS", "OAS3", "PBX1", "PDLIM2", "PDLIM5", 
                         "SERPINE2", "SLC6A12", "STYK1", "TGM2", "SLC7A2", "ZNF704",
                         "AKAP12", "CACNA2D2", "CALR", "CPSF4", "CST7", "FKBP11",
                         "FNDC3B", "MS4A3", "MST1", "P2RY2", "P4HB", "PHKG1", "PMP22",
                         "PRAG1", "TDRD9", "TRIM9", "UXS1")


# draw the scatterplot
scatterplot <- ggplot(merged_table, aes(x = logFC_BEAT, y = logFC_mouse, color = color_sign)) +
  geom_point(alpha = 0.75) + theme_bw() + xlab("log2FC BEAT") + ylab("log2FC mouse") +
  xlim(-max(abs(merged_table$logFC_BEAT)), max(abs(merged_table$logFC_BEAT))) +
  ylim(-max(abs(merged_table$logFC_mouse)), max(abs(merged_table$logFC_mouse))) +
  scale_color_manual(values = c("n.s" = "black", "sign. mouse" = "blue", 
                                "sign. BEAT" = "orange", "sign. mouse and BEAT" = "red")) +
  geom_text_repel(data = subset(merged_table[merged_table$external_gene_name %in% list_relevant_genes,]), 
                  color = "black", aes(label = external_gene_name), size = 5) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between BEAT and GFP+ TPO vs PBS; corr coef (Pearson) = ", coeff_corr))


# save the scatterplot in files
saveFigures("scatterplot_GFP+_TPOvsPBS_relevantGenes", ggplot = scatterplot, dirPlot = outputPath, A4 = T)

write_xlsx(merged_table, path = "DEA/Human_patients/BEAT/allAML/BEAT_ERG_selection_kmeans/mergedTable_BEAT_GFP+_TPOvsPBS.xlsx")



### perform GSEA and ORA ---------------------


# get the output path
outputPath <- "DEA/Human_patients/BEAT/allAML/BEAT_ERG_selection_kmeans/"

# load the table with DEA result for the cutoff of 0.8
MLLr_BEAT <- read_xlsx(path = file.path(outputPath, "results_DEA.xlsx"), 
                       sheet = 1)


#### perform GSEA -----------------------------------------------------------

# create the output path for the GSEA
outputPath_GSEA <- create_dir(path = file.path(outputPath, "GSEA"))

# select the FC
subtable <- MLLr_BEAT[!is.na(MLLr_BEAT$entrezgene_id),]
id_FC <- sort(subtable$logFC, index.return = T, decreasing = T)
FC <-  subtable$logFC[id_FC$ix]
names(FC) <- subtable$entrezgene_id[id_FC$ix]

# removes duplicate
FC = FC[!duplicated(names(FC))]

# init the list for camera results
list_cameras_human <- NULL

# save the different collection from Msig database for camera
db_cameras_human = c("H", "c2", "c3", "c4", "c5", "c6", "c7")

# start the loop
for(db in db_cameras_human){
  
  # create the filename
  filename = paste0("human_", db, "_v5p2.rdata")
  
  # create the file
  #file = file.path(getwd(), filename)
  
  # download the file
  if(!file.exists(filename)){
    print(paste0("download of ",  filename, " in progress..."))
    download.file(paste0("http://bioinf.wehi.edu.au/software/MSigDB/", filename), filename)
    print(paste0(filename, " downloaded."))
  }
  load(filename)
  print(paste0(filename, " loaded."))
  
}

# save the list of databases
dbs <- paste0("Hs.", c("H", paste0("c", 2:7)))


# analyze in a loop
list_fGSEA_results <- NULL
for(db in dbs){
  
  print(paste0(date(), "; GSEA for ", db, " in progress..."))
  
  # run the fGSEA
  fGSEA_Res <- fgsea(get(db), FC)
  print("GSEA")
  # select the significant GSEA
  fGSEA_Res <- as_tibble(fGSEA_Res)
  fGSEA_Res = dplyr::arrange(fGSEA_Res, pval) %>% filter(padj <= 0.05)
  
  # save in the list
  list_fGSEA_results = c(list_fGSEA_results, list(as.data.frame(fGSEA_Res)))
  
  print(paste0(date(), "; GSEA for ", db, " done."))
  
  
}

# rename the element in the list to provide this name in each sheet within the excel file
names(list_fGSEA_results) = dbs


# save in excel file
write_xlsx(list_fGSEA_results, path = file.path(outputPath_GSEA, "GSEA_Msigdatabase_MLLr_BEAT_kmeans.xlsx"))


# provide the Msigdbr database
msigdbr_species()
hs_tg2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)

# select the hoxa9 pathways

# perform the GSEA for Hoxa9
em2 <- GSEA(FC, TERM2GENE = hs_tg2)
saveFigures("GSEA_C2", ggplot = dotplot(em2, showCategory=30), 
            dirPlot = outputPath_GSEA, A4 = T)



##### perform GSEA with selected marks for custom pathways (source: https://doi.org/10.1182/bloodadvances.2018025866) -----------------------------------------------------------

# download the KEGG pathway
mmu.kegg <- getGeneKEGGLinks(species.KEGG = "mmu")

# download the Jak-Stat pathway
jakStat <- mmu.kegg[grep(pattern = "mmu04630", mmu.kegg$PathwayID), ]

# download the NF_kappaB pathway
NF_kappaB <- mmu.kegg[grep(pattern = "mmu04064", mmu.kegg$PathwayID), ]

# download the MAPK pathway
MAPK <- mmu.kegg[grep(pattern = "mmu04010", mmu.kegg$PathwayID), ]

# download the proteoglycans pathway
proteoglycans <- mmu.kegg[grep(pattern = "mmu05205", mmu.kegg$PathwayID), ]


# update the gene names
jakStat <- merge(x = jakStat, y = Mouse_genes, by.x = "GeneID", by.y  = "entrezgene_id")
NF_kappaB <- merge(x = NF_kappaB, y = Mouse_genes, by.x = "GeneID", by.y  = "entrezgene_id")
MAPK <- merge(x = MAPK, y = Mouse_genes, by.x = "GeneID", by.y  = "entrezgene_id")
proteoglycans <- merge(x = proteoglycans, y = Mouse_genes, by.x = "GeneID", by.y  = "entrezgene_id")

# create the list of genes for the TPO pathway
# source (paper: 10.1007/s12079-018-0480-4)
THPO_pathway <- c("THPO", "MPL", "PM", "BMPR1A", "BMPR2", "PRKCA", "PRKCB", 
                  "PRKCZ", "PRKCD", "PRKCE", "PRKACA", "CTTN", "TNS2", "IRS2", 
                  "GAB2", "GRB2", "SRC", "SYK", "CISH", "GAB1", "PTPN11", "SHC1", 
                  "PIK3R1", "PIK3CA", "AKT1", "MTOR", "RPS6KB1", "EIF4EBP1", 
                  "PLCG1", "SLC2A1", "NFKB1", "BCL2L1", "BAD", "FOXO3", "CDKN1A", 
                  "CDKN1B", "LYN", "HRAS", "BRAF", "RAP1A", "KRAS", "RAF1", "GSK3B", "MAPK8", 
                  "MAPK9", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RPS6KA2", "CDK4",
                  "CDK6", "ELK1", "CREB1", "MEIS1", "HOXA9", "RUNX1", "MAPK14",
                  "SOS1", "PTPN6", "CBL", "PIK3R1", "PIK3R2", "TYK2", "CRKL", 
                  "STAT5A", "STAT5B", "POU2F1", "CRK", "ATXN2L", "VAV1", "TEC",
                  "CASP3", "MYC", "PIM1", "OSM", "HIF1A", "STAT1", "STAT3", 
                  "GATA1", "PDGFRA", "PDGFRB", "NRP1", "BMP4", "IRF2", "PTK2", 
                  "PARP1", "STAT6", "SMAD1", "SMAD5", "SMAD9", "SMAD4", "TGFB1",
                  "SMAD2", "SMAD3", "ID1", "ID2", "ID3", "MX1", "MX2", "P4HA2",
                  "P4HA3", "FN", "COL3A1", "COL4A1", "INPP5D", "JAK2") 



# update the gene names
jakStat <- merge(x = jakStat, y = Mouse_genes, by.x = "GeneID", by.y  = "entrezgene_id")
NF_kappaB <- merge(x = NF_kappaB, y = Mouse_genes, by.x = "GeneID", by.y  = "entrezgene_id")
MAPK <- merge(x = MAPK, y = Mouse_genes, by.x = "GeneID", by.y  = "entrezgene_id")
proteoglycans <- merge(x = proteoglycans, y = Mouse_genes, by.x = "GeneID", by.y  = "entrezgene_id")


# create the custom pathways
custom_pathways <- list("Jak-Stat pathway" = ortho$ensembl_gene_id[ortho$mmusculus_homolog_ensembl_gene %in% jakStat$ensembl_gene_id],
                        "MapK pathway" = ortho$ensembl_gene_id[ortho$mmusculus_homolog_ensembl_gene %in% MAPK$ensembl_gene_id],
                        "NF-Kappa-B pathway" = ortho$ensembl_gene_id[ortho$mmusculus_homolog_ensembl_gene %in% NF_kappaB$ensembl_gene_id],
                        "proteoglycans pathway" = ortho$ensembl_gene_id[ortho$mmusculus_homolog_ensembl_gene %in% proteoglycans$ensembl_gene_id],
                        "THPO_pathway" = ortho$ensembl_gene_id[ortho$external_gene_name %in% THPO_pathway]
)


# perform the GSEA
fgseaRes_custom_pathways <- fgsea(pathways = custom_pathways, stats = geneList_GSEA)

# plot the enrichment
for(enrich in names(custom_pathways)){
  
  # define the filename
  filename <- gsub(pattern = " ", replacement = "_", enrich)
  filename = paste0("GSEA_", filename, "_MLLr_BEAT_kmeans")
  
  # save the enrichment in the pdf
  saveFigures(fileName = filename, ggplot = plotEnrichment(custom_pathways[[enrich]], geneList_GSEA) + 
                labs(title = enrich), dirPlot = outputPath_GSEA, A4 = T)
  
}

# write the excel file
write_xlsx(fgseaRes_custom_pathways, path = file.path(outputPath_GSEA, "GSEA_custom_pathway_Msigdatabase_MLLr_BEAT_kmeans.xlsx"))



#### perform ORA  -----------------------------------------------------------

# create the directory
outputPath_ORA <- create_dir(path = file.path(outputPath, "ORA"))

# select all tested genes
genes <- MLLr_BEAT$Ensembl_id

# select the significant toptags
dge_DEG <- MLLr_BEAT[MLLr_BEAT$FDR <= 0.05, ]

# create the output directory
tmp_dir <- create_dir(path = file.path(outputPath_ORA, "DEG_woUniverse"))

# perform the enrichment analysis
EnrichmentTerm_Analysis(common = dge_DEG$Ensembl_id, outputDir = tmp_dir, 
                        mart = ensembl, isMouse = F)
# create the output directory
tmp_dir <- create_dir(path = file.path(outputPath_ORA, "DEG"))

# perform the enrichment analysis
EnrichmentTerm_Analysis(common = dge_DEG$Ensembl_id, outputDir = tmp_dir, 
                        mart = ensembl, isMouse = F, universe = genes)



# create the survival curve ---------------------

# create the output path
outdir_survival_BEAT <- create_dir(path = "DEA/Human_patients/BEAT/allAML/BEAT_ERG_selection_kmeans/Survival_curve")

# select the samples for the survival curve
samples_id <- read_xlsx(path = "DEA/Human_patients/BEAT/allAML/BEAT_ERG_selection_kmeans/selected_patients.xlsx", sheet = 1)

# select patient high
patient_high <- samples_id$Id_sample[samples_id$Group == "high"]

# select patient low
patient_low <- samples_id$Id_sample[samples_id$Group == "low"]

# add patients id in the sample table
BEAT_cbioportal$Sample.ID_red <- paste0(sapply(BEAT_cbioportal$Sample.ID, function(x)strsplit(x, split = "_")[[1]][5]), "R")
BEAT_cbioportal$Patient.ID_red <- paste0(sapply(BEAT_cbioportal$Sample.ID, function(x)strsplit(x, split = "_")[[1]][4]))

# calculate the days to death for the dead patients
BEAT_cbioportal$days_to_death = round(BEAT_cbioportal$Overall.Survival..Months. * 30.4)


## create the survival curve for selected genes ---------------------
# selection based on the comparison of DEG between the different databases

# create list of Up-Up list
UpUp_DEG <- c("IL12RB2", "INPP4B", "MECOM", "ADGRG6", "MMP28", "SEMA4F", "SH3BP5", "AIF1",
              "CDC14A", "FARP1", "GRAP2", "KRI1", "MLLT3",  "MOCOS", "OAS3","PBX1", "PDLIM2", "PDLIM5", "SERPINE2", "SLC6A12", 
              "STYK1", "TGM2")
DownDown_DEG <- c("AKAP12", "CACNA2D2", "CALR", "CPSF4", "CST7", "FKBP11", "FNDC3B", "MS4A3", 
                  "MST1", "P2RY2", "P4HB", "PHKG1", "PMP22", "PRAG1", "SLC7A2", "TDRD9", "TRIM9", "UXS1", "ZNF704")

# create the dataframe
dt_BEAT_survival_all <- data.frame("Sample_id" = BEAT_cbioportal$Sample.ID_red[BEAT_cbioportal$Sample.ID_red %in% all_samples_ML], 
                               "Patient_id" = BEAT_cbioportal$Patient.ID_red[BEAT_cbioportal$Sample.ID_red %in% all_samples_ML],
                               "days_to_death" = BEAT_cbioportal$days_to_death[BEAT_cbioportal$Sample.ID_red %in% all_samples_ML],
                               "vital_status" = BEAT_cbioportal$Overall.Survival.Status[BEAT_cbioportal$Sample.ID_red %in% all_samples_ML],
                               "days_to_last_follow_up" = BEAT_cbioportal$days_to_death[BEAT_cbioportal$Sample.ID_red %in% all_samples_ML])

# save the Rds
saveRDS(dt_BEAT_survival_all, file = "DEA/Human_patients/BEAT/dt_BEAT_survival_all.rds")

# select only patients from GDC data
GDC_patients_info <- read.table(file = "Dataset/BEAT/clinical.cases_selection.2024-02-07/clinical.tsv", header = T, sep = "\t", quote = "")

# select the samples
samples_id <- dt_BEAT_survival_all$Sample_id[dt_BEAT_survival_all$Patient_id %in% GDC_patients_info$case_submitter_id]

# remove duplicated patients
dt_BEAT_survival_all =  dt_BEAT_survival_all[!duplicated(dt_BEAT_survival_all$Patient_id), ]

# update the dt_BEAT_survival_all table
dt_BEAT_survival_all_sel2 <- dt_BEAT_survival_all[dt_BEAT_survival_all$Sample_id %in% samples_id, ]

# update the vital status
dt_BEAT_survival_all_sel2$vital_status <- gsub(pattern = "0:LIVING", replacement = "Alive", dt_BEAT_survival_all_sel2$vital_status)
dt_BEAT_survival_all_sel2$vital_status <- gsub(pattern = "1:DECEASED", replacement = "Dead", dt_BEAT_survival_all_sel2$vital_status)

# update information about the death of patients
dt_BEAT_survival_all_sel2$days_to_death[grepl(pattern = "Alive", dt_BEAT_survival_all_sel2$vital_status)] = NA
dt_BEAT_survival_all_sel2$days_to_last_follow_up[grepl(pattern = "Dead", dt_BEAT_survival_all_sel2$vital_status)] = NA

# select the data_BEAT
data_BEAT_sel <- data_BEAT[, samples_id]

# draw survival curve for each Up genes
for(gene in c(UpUp_DEG, DownDown_DEG, "HOXA9", "MEIS1")){
  
  # select the ensembl id of the gene
  gene_id_tmp <- rowData(data_BEAT_sel)$gene_id[rowData(data_BEAT_sel)$gene_name == gene] 
  
  # select value of genes expression
  exp <- log2(table_BEAT_tpm[gene_id_tmp, dt_BEAT_survival_all_sel2$Sample_id])
  
  # calculate the median
  median_exp <- median(exp)
  
  # select the samples
  above_sample <- names(exp)[exp >= median_exp]
  below_sample <- names(exp)[exp < median_exp]
  
  # select the patient
  above_patient <- intersect(above_sample, dt_BEAT_survival_all$Sample_id)
  below_patient <- intersect(below_sample, dt_BEAT_survival_all$Sample_id)
  
  # give information about the statistics
  print(paste0(gene, ": above median = ", length(above_patient), "; below median = ", length(below_patient)))
  
  # do the selection for the clinical_TCGA
  dt_BEAT_survival_all_selection = dt_BEAT_survival_all_sel2[dt_BEAT_survival_all_sel2$Sample_id %in% c(above_patient, below_patient), ]
  
  # add information about the expression
  dt_BEAT_survival_all_selection$exp = "below median"
  dt_BEAT_survival_all_selection$exp[dt_BEAT_survival_all_selection$Sample_id %in% above_patient] = "above median"
  
  # draw the survival curve
  TCGAanalyze_survival(
    data = dt_BEAT_survival_all_selection,
    clusterCol = "exp",
    main = paste0(gene, "; BEAT, Kaplan-Meier survival curve"),
    filename = file.path(outdir_survival_BEAT, paste0("survival_curve_", gene, "_allPatients_v3.pdf")),
    color = c("#ED1C24", "#0F75BC"),
    conf.int = F,
    width = 8,
    xlim = c(0, 2000),
    shape = 20
    )
  
}



# create the scatterplot for common DEGs compare to MECOM ---------------------


# get the id for IL12RB2
il12rb2_id = rowData(data_BEAT)$gene_id[rowData(data_BEAT)$gene_name == "IL12RB2"]

# get the id for INPP4B
INPP4B_id = rowData(data_BEAT)$gene_id[rowData(data_BEAT)$gene_name == "INPP4B"]


# create the dataframe for the scatterplot figure
dt_BEAT_ERG_MECOM <- data.frame("MECOM" = table_BEAT_tpm[mecom_id,],
                                "IL12RB2" = table_BEAT_tpm[il12rb2_id, ],
                                "INPP4B" = table_BEAT_tpm[INPP4B_id, ],
                                "Sample" = colnames(table_BEAT_tpm),
                                "Color" = "none")

# select only the samples
dt_BEAT_ERG_MECOM = dt_BEAT_ERG_MECOM[samples_id,]

# update the patients information
dt_BEAT_ERG_MECOM$Color[dt_BEAT_ERG_MECOM$Sample %in% common_samples] = "MLL-AF9"
dt_BEAT_ERG_MECOM$Color[dt_BEAT_ERG_MECOM$Sample %in% all_samples_selected_allAML_EVI1] = "EVI1r"

# calculate the coefficient correlation
dt_BEAT_ERG_MECOM_tmp <- dt_BEAT_ERG_MECOM[dt_BEAT_ERG_MECOM$MECOM != 0, ]
dt_BEAT_ERG_MECOM_tmp <- dt_BEAT_ERG_MECOM_tmp[dt_BEAT_ERG_MECOM_tmp$IL12RB2 != 0, ]
fit <- lm(log2(dt_BEAT_ERG_MECOM_tmp$IL12RB2) ~ log2(dt_BEAT_ERG_MECOM_tmp$MECOM), data = dt_BEAT_ERG_MECOM_tmp)
coeff_corr <- round(cor(log2(dt_BEAT_ERG_MECOM_tmp$MECOM), log2(dt_BEAT_ERG_MECOM_tmp$IL12RB2), method = "pearson"), 3)

# create the scatterplot
scatterplot <- getScatterplot(data = dt_BEAT_ERG_MECOM_tmp, x = "MECOM", y = "IL12RB2", log2 = T)
scatterplot <- scatterplot + geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between MECOM and IL12RB2; corr coef (Pearson) = ", coeff_corr)) +
  scale_color_manual(values = c("none" = "purple", "MLLr" = "cyan3", 
                                "MLL-AF9" = "forestgreen", "EVI1r" = "red")) + geom_point(size=3)


# create the scatterplot for the Erg vs MECOM expression for all patients
saveFigures("scatterplot_IL12RB2vsMECOM_v2", scatterplot, 
            dirPlot = "DEA/Human_patients/BEAT/allAML/", A4 = T)

# calculate the coefficient correlation
dt_BEAT_ERG_MECOM_tmp <- dt_BEAT_ERG_MECOM[dt_BEAT_ERG_MECOM$MECOM != 0, ]
dt_BEAT_ERG_MECOM_tmp <- dt_BEAT_ERG_MECOM_tmp[dt_BEAT_ERG_MECOM_tmp$INPP4B != 0, ]
fit <- lm(log2(dt_BEAT_ERG_MECOM_tmp$INPP4B) ~ log2(dt_BEAT_ERG_MECOM_tmp$MECOM), data = dt_BEAT_ERG_MECOM_tmp)
coeff_corr <- round(cor(log2(dt_BEAT_ERG_MECOM_tmp$MECOM), log2(dt_BEAT_ERG_MECOM_tmp$INPP4B), method = "pearson"), 3)

# create the scatterplot
scatterplot <- getScatterplot(data = dt_BEAT_ERG_MECOM_tmp, x = "MECOM", y = "INPP4B", log2 = T)
scatterplot <- scatterplot + geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between MECOM and INPP4B; corr coef (Pearson) = ", coeff_corr)) +
  scale_color_manual(values = c("none" = "purple", "MLLr" = "cyan3", 
                                "MLL-AF9" = "forestgreen", "EVI1r" = "red")) + geom_point(size=3)


# create the scatterplot for the Erg vs MECOM expression for all patients
saveFigures("scatterplot_INPP4BvsMECOM_v2", scatterplot, 
            dirPlot = "DEA/Human_patients/BEAT/allAML/", A4 = T)
