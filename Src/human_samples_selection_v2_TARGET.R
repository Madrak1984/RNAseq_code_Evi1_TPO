##### command lines to select the human samples, 2nd version, only for TARGET selection
# done by Jonathan Seguin, group of Prof. Schwaller, DBM, UKBB, Basel, Switzerland
# email: jonathan.seguin@unibas.ch, seguin.jonathan@gmail.com
# Wed Jul 26 14:10:27 2023


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
library(RColorBrewer)
library(ggrepel)
library(writexl)
library(plotly)
library(rgl)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(gridExtra)
library(scales)
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
      ggtitle(paste0("scatterplot, ", x, " vs ", y, ", TARGET database")) 
  } else {
    scatterplot =  ggplot(data = data, aes(x = data[, x], y = data[, y], color = Color)) + 
      geom_point() +
      ylab(paste0(y, " (tpm)")) + 
      xlab(paste0(x, " (tpm)")) + 
      theme_bw() + 
      ggtitle(paste0("scatterplot, ", x, " vs ", y, ", TARGET database"))
  }

    
  # return the scatterplot
  return(scatterplot)
  
}



## read the tables ---------------------------------------------------------

# read the TCGA file
TCGA_table <- read.xlsx2(file = "Dataset/TARGET/combined_study_clinical_data_TCGA-AML.xlsx",
                         sheetIndex = 1)


# load the table from the TPO-Evi1 RNAseq
fits_tmp <- readRDS(file = "DEA/Rds/RDS_woD/edge_samplesR_fits_results_DEA.rds")

# select the orthologous genes
ortho <- getHumanOrthologous()

# select the mouse genes
Mouse_genes <- getMouseGene()

# select Human genes
Human_genes <-getHumanGene()



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
                                y = ortho, 
                                by.x = "ensembl_gene_id", 
                                by.y = "mmusculus_homolog_ensembl_gene")

# select the DEG based on FDR <= 5%
DEG_mouse_bulk <- fit_bulkTPOVsnoTPO[fit_bulkTPOVsnoTPO$FDR <= 0.05, ]


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
                            y = ortho, 
                            by.x = "ensembl_gene_id", 
                            by.y = "mmusculus_homolog_ensembl_gene")

# select the DEG based on FDR <= 5%
DEG_mouse_GFP <- fit_GFP_TPOVsnoTPO[fit_GFP_TPOVsnoTPO$FDR <= 0.05, ]



# create the directories output ------------------------------------------------------
outdir_db <- "Human_comparison/Plots/TARGET_v2"
if(!dir.exists(outdir_db)) dir.create(outdir_db, recursive = T)

create_dir(path = "DEA/Human_patients/TARGET/")

# download gene information ------------------------------------------------------

# download ensembl information for human genes
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
                   host = "sep2019.archive.ensembl.org")

# include patients from TARGET database ---------------------
# source (https://www.bioconductor.org/packages/devel/bioc/vignettes/GenomicDataCommons/inst/doc/overview.html#1_What_is_the_GDC)

query_Target <- GDCquery(project = "TARGET-AML",
                         data.category = "Transcriptome Profiling",
                         data.type = "Gene Expression Quantification",
                         workflow.type = "STAR - Counts")

GDCdownload(query_Target)
data_Target <- GDCprepare(query_Target)

# save the data_Target
saveRDS(data_Target, "DEA/Human_patients/TARGET/data_Target.rds")

query_clinical_TARGET <- GDCquery(project = "TARGET-AML",
                                  data.category = "Clinical", 
                                  data.type = "Clinical Supplement")
query_clinical_TARGET_2 <- GDCquery_clinic(project = "TARGET-AML")

# save the query_clinical_TARGET
saveRDS(query_clinical_TARGET_2, file = "DEA/Human_patients/TARGET/query_clinical_TARGET_2.rds")
                                  


# find outlier samples
dataPrep_Target <- TCGAanalyze_Preprocessing(object = data_Target, 
                                             cor.cut = 0.6,
                                             datatype = "unstranded")                      

# perform normalization
dataNorm_Target <- TCGAanalyze_Normalization(tabDF = dataPrep_Target,
                                             geneInfo = geneInfoHT,
                                             method = "gcContent") 

# filter the samples
dataFilt_Target <- TCGAanalyze_Filtering(tabDF = dataNorm_Target,
                                         method = "quantile", 
                                         qnt.cut =  0.25) 

# get the tpm value for all gene expression
table_Target_tpm <- assays(data_Target)$tpm_unstrand

# select the identifier for MECOM
mecom_id = rowData(data_Target)$gene_id[rowData(data_Target)$gene_name == "MECOM"]

# load the patients table
table_AML <- read.xlsx2(file = "Dataset/TARGET/TARGET_AML_ClinicalData_Discovery_20181213.xlsx",
                        header = T, sheetIndex = 1) 

# change table_AML as tibble
table_AML <- as_tibble(table_AML)

# remove the TARGET-21 code
table_AML <- table_AML[grep(pattern = "TARGET-20", table_AML$TARGET.USI), ]


# select the table with the KMT2A patients
table_AML_MF9 <- table_AML[grep(pattern = "KMT2A-MLLT3",
                                  table_AML$Gene.Fusion), ]

# select the information of patients
colData_table <- as.data.frame(colData(data_Target))

# remove control samples
colData_table <- colData_table[colData_table$tumor.code == "20",]


# select the patients
patients <- colData_table$barcode[colData_table$patient %in% table_AML_MF9$TARGET.USI]

# save the patients selection in an excel file
write.xlsx2(as.data.frame(table_AML[table_AML$TARGET.USI %in% intersect(colData_table$patient,
                                                                        table_AML_MF9$TARGET.USI),]), file = file.path(outdir_db, "selected_patient_TARGET-AML.xlsx"), 
            col.names = T, row.names = F)


# select patients on all MLL
patients_allKMT2A <- table_AML[grep(pattern = "KMT2A", table_AML$Gene.Fusion), ]
# select the table with the KMT2A patients
table_AML_MLLr <- table_AML[grep(pattern = "KMT2A",
                                table_AML$Gene.Fusion), ]

# save the patients selection in an excel file
write.xlsx2(as.data.frame(patients_allKMT2A), file = file.path(outdir_db, "selected_patient_allKMT2A_TARGET-AML.xlsx"), 
            col.names = T, row.names = F)


# select patients on all EVI1
table_AML_allEVI1 <- table_AML[grep(pattern = "q26", table_AML$ISCN), ]
table_AML_allEVI1 <- table_AML_allEVI1[grep(pattern = "\\(3", table_AML_allEVI1$ISCN), ]



# perform DEA based on k-means clustering result ---------------------

## find the cluster of samples ---------------------

# create the dataframe for the scatterplot figure
dt_TARGET_ERG_MECOM_log2tpm <- data.frame("MECOM" = log2(table_Target_tpm[mecom_id, colData_table$barcode]),
                                  "ERG" = log2(table_Target_tpm[erg_id, colData_table$barcode]),
                                  "Sample" = colData_table$barcode,
                                  "Color" = "none")




# update the patients information
dt_TARGET_ERG_MECOM_log2tpm$Color[dt_TARGET_ERG_MECOM_log2tpm$Sample %in% colData_table$barcode[colData_table$patient %in% patients_allKMT2A$TARGET.USI]] = "MLLr"
dt_TARGET_ERG_MECOM_log2tpm$Color[dt_TARGET_ERG_MECOM_log2tpm$Sample %in% patients] = "MLL-AF9"
dt_TARGET_ERG_MECOM_log2tpm$Color[dt_TARGET_ERG_MECOM_log2tpm$Sample %in% colData_table$barcode[colData_table$patient %in% table_AML_allEVI1$TARGET.USI]] = "EVI1r"

# remove infinite values
dt_TARGET_ERG_MECOM_log2tpm = dt_TARGET_ERG_MECOM_log2tpm[is.finite(dt_TARGET_ERG_MECOM_log2tpm$MECOM), ]

# initialize list to save kmeans results and their corresponding plots
list_kmeans <- NULL
list_kmean_plot <- NULL


# calculate the different k-means clusters
for(k in 2:9){
  
  # calculate the temporary clustering
  tmp_kmeans <-  kmeans(dt_TARGET_ERG_MECOM_log2tpm[, 1:2], centers = k, nstart = 25)
  
  # save the result in the list
  list_kmeans <- c(list_kmeans, list(tmp_kmeans))
  
  # save the plot
  list_kmean_plot <- c(list_kmean_plot, 
                      list(fviz_cluster(tmp_kmeans, data = dt_TARGET_ERG_MECOM_log2tpm[, 1:2], geom = "point", ggtheme = theme_bw()) + ggtitle(paste0("k = ", k))))
  
}

# update the DEA output directory
DEA_dir <- create_dir("DEA/Human_patients/TARGET/allAML/")

# renames the list
names(list_kmeans) <- paste0("k", 2:9)
names(list_kmean_plot) <- paste0("k", 2:9)

# save in rds file
saveRDS(list_kmeans, file = "DEA/Human_patients/TARGET/allAML/TARGET_ERG_selection_kmeans/list_kmeans_TARGET.rds")
saveRDS(list_kmean_plot, file = "DEA/Human_patients/TARGET/allAML/TARGET_ERG_selection_kmeans/list_kmeans_plot_TARGET.rds")

# display the result
g <- grid.arrange(grobs = list_kmean_plot, ncol = 4, nrow = 2)

# save in pdf file
ggsave(filename = "DEA/Human_patients/TARGET/allAML/kmeans_clustering_ERGvsMECOM_TARGET.pdf",
       plot = g, width =297, height = 210, units = "mm")

# determine the optiimal number of cluster
set.seed(123)
gap_stat <- clusGap(dt_TARGET_ERG_MECOM_log2tpm[, 1:2], FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)

# save the result in the pdf
pdf(file = "DEA/Human_patients/TARGET/allAML/bestClusterNumber.pdf", paper = "a4r")
set.seed(123)
fviz_nbclust(dt_TARGET_ERG_MECOM_log2tpm[, 1:2], kmeans, method = "wss")
set.seed(123)
fviz_nbclust(dt_TARGET_ERG_MECOM_log2tpm[, 1:2], kmeans, method = "silhouette")
fviz_gap_stat(gap_stat)
dev.off()


# selection of cutoff 4 for the clustering result
# update the scatterplot 
dt_TARGET_ERG_MECOM_log2tpm$cluster = list_kmeans$k4$cluster




## perform the DEA based on clusters (4 vs 1) ---------------------

# filter the patient based on the ERG expression
patient_kept <- dt_TARGET_ERG_MECOM_log2tpm$Sample[dt_TARGET_ERG_MECOM_log2tpm$ERG > 0]

# select patient high (cluster 4)
patient_high <- intersect(dt_TARGET_ERG_MECOM_log2tpm$Sample[dt_TARGET_ERG_MECOM_log2tpm$cluster == "4"], 
                          patient_kept)

# select patient low (cluster 1)
patient_low <- intersect(dt_TARGET_ERG_MECOM_log2tpm$Sample[dt_TARGET_ERG_MECOM_log2tpm$cluster == "1"], 
                         patient_kept)


# create output directory
tmp_DEA <- create_dir(path = file.path(DEA_dir, "TARGET_ERG_selection_kmeans"))

# select the samples
table_TARGET_unstr_sel <- table_TARGET_unstr[, colnames(table_TARGET_unstr) %in% c(patient_high, patient_low)]

# replace NA value in table
table_TARGET_unstr_sel[is.na(table_TARGET_unstr_sel)] = 0

# determine the group
group <- rep("low", ncol(table_TARGET_unstr_sel))

# update the group
group[colnames(table_TARGET_unstr_sel) %in% patient_high] = "high"

# remove the gene lowly expressed
print(paste0("# genes: ", nrow(table_TARGET_unstr_sel)))
# keep gene with more than 10 reads in total
count_10 <- nrow(table_TARGET_unstr_sel[rowSums(table_TARGET_unstr_sel) > 10,])
print(paste0("# kept genes (filter based on 10): ", count_10))
# 51224 kept genes among 60660 genes

# create the dge
dge <- DGEList(counts = table_TARGET_unstr_sel, group = group)

# filter according to the filtering gene expression
keep <- filterByExpr(dge)
count_filter <- length(keep[keep])
print(paste0("# kept gene after filtering by function: ", count_filter))
# 24312 genes kept

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
saveFigures("BCVplot", ggplot = plotBCV(dge), dirPlot = tmp_DEA, A4 = T)

# perform the exact test
et <-  exactTest(dge)

# provide all results
allresult <- as.data.frame(topTags(et, n=nrow(et)))

# add ensembl gene id
allresult <-  cbind(allresult, "Ensembl_id" = sapply(rownames(allresult),
                                                     function(x) strsplit(x, split = "\\.")[[1]][1]))

# merge with the human gene information
allresult <- merge(x = allresult, 
                   y = Human_genes, 
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
write_xlsx(x = list_result, path = file.path(tmp_DEA, "results_DEA.xlsx"), col_names = T)
gc()

# save the selected samples
write.xlsx2(data.frame("Id_sample" = colnames(table_TARGET_unstr_sel), "Group" = group), file = file.path(tmp_DEA, "selected_patients.xlsx"),
            col.names = T, row.names = F)
gc()




### perform GSEA and ORA ---------------------


# get the output path
outputPath <- tmp_DEA

# load the table with DEA result for the cutoff of 0.8
MLLr_TARGET_cutoff_ERG <- read_xlsx(path = file.path(outputPath, "results_DEA.xlsx"), 
                                    sheet = 1)


#### perform GSEA -----------------------------------------------------------

# create the output path for the GSEA
outputPath_GSEA <- create_dir(path = file.path(outputPath, "GSEA"))

# select the FC
subtable <- MLLr_TARGET_cutoff_ERG[!is.na(MLLr_TARGET_cutoff_ERG$entrezgene_id),]
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

# rename the element in the list to prived this name in each sheet within the excel file
names(list_fGSEA_results) = dbs


# save in excel file
write_xlsx(list_fGSEA_results, path = file.path(outputPath_GSEA, "GSEA_Msigdatabase_MLLr_TARGET_kmeans.xlsx"))


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
  filename = paste0("GSEA_", filename, "_MLLr_TARGET_kmeans")
  
  # save the enrichment in the pdf
  saveFigures(fileName = filename, ggplot = plotEnrichment(custom_pathways[[enrich]], geneList_GSEA) + 
                labs(title = enrich), dirPlot = outputPath_GSEA, A4 = T)
  
}

# write the excel file
write_xlsx(fgseaRes_custom_pathways, path = file.path(outputPath_GSEA, "GSEA_custom_pathway_Msigdatabase_MLLr_TARGET_kmeans.xlsx"))



#### perform ORA  -----------------------------------------------------------

# create the directory
outputPath_ORA <- create_dir(path = file.path(outputPath, "ORA"))

# select all tested genes
genes <- MLLr_TARGET_cutoff_ERG$Ensembl_id

# select the significant toptags
dge_DEG <- MLLr_TARGET_cutoff_ERG[MLLr_TARGET_cutoff_ERG$FDR <= 0.05, ]

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



### compare with the mouse DEG within a scatterplot ---------------------

# get the output path
outputPath <- tmp_DEA

# load the table with DEA result for the cutoff of 0.8
MLLr_TARGET_cutoff_ERG <- read_xlsx(path = file.path(outputPath, "results_DEA.xlsx"), 
                                    sheet = 1)

# merge the tables
merged_table <- merge(x = MLLr_TARGET_cutoff_ERG,
                      y = fit_bulkTPOVsnoTPO,
                      by = "external_gene_name",
                      suffixes = c("_TARGET", "_mouse"))

# create the color for the significativity 
merged_table$color_sign <- "n.s"

# update the color based on the significativity
merged_table$color_sign[merged_table$FDR_mouse <= 0.05] <- "sign. mouse"
merged_table$color_sign[merged_table$FDR_TARGET <= 0.05] <- "sign. TARGET"
merged_table$color_sign[(merged_table$FDR_mouse <= 0.05) & (merged_table$FDR_TARGET <= 0.05)] <- "sign. mouse and TARGET"

# calculate the coefficient correlation
fit <- lm(merged_table$logFC_TARGET ~ merged_table$logFC_mouse, data = merged_table)
coeff_corr <- round(cor(merged_table$logFC_TARGET, merged_table$logFC_mouse, method = "pearson"), 3)

# draw the scatterplot
scatterplot <- ggplot(merged_table, aes(x = logFC_TARGET, y = logFC_mouse, color = color_sign)) +
  geom_point(alpha = 0.75) + theme_bw() + xlab("log2FC TARGET") + ylab("log2FC mouse") +
  xlim(-max(abs(merged_table$logFC_TARGET)), max(abs(merged_table$logFC_TARGET))) +
  ylim(-max(abs(merged_table$logFC_mouse)), max(abs(merged_table$logFC_mouse))) +
  scale_color_manual(values = c("n.s" = "black", "sign. mouse" = "blue", 
                                "sign. TARGET" = "orange", "sign. mouse and TARGET" = "red")) +
  geom_text_repel(data = subset(merged_table[merged_table$external_gene_name %in% merged_table$external_gene_name[merged_table$color_sign == "sign. mouse and TARGET"],]), 
                  color = "black", aes(label = external_gene_name), size = 5) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between TARGET and bulk TPO vs PBS; corr coef (Pearson) = ", coeff_corr))


# save the scatterplot in files
saveFigures("scatterplot_bulk_TPOvsPBS", ggplot = scatterplot, dirPlot = outputPath, A4 = T)




# merge the tables
merged_table <- merge(x = MLLr_TARGET_cutoff_ERG,
                      y = fit_GFP_TPOVsnoTPO,
                      by = "external_gene_name",
                      suffixes = c("_TARGET", "_mouse"))


# create the color for the significativity 
merged_table$color_sign <- "n.s"

# update the color based on the significativity
merged_table$color_sign[merged_table$FDR_mouse <= 0.05] <- "sign. mouse"
merged_table$color_sign[merged_table$FDR_TARGET <= 0.05] <- "sign. TARGET"
merged_table$color_sign[(merged_table$FDR_mouse <= 0.05) & (merged_table$FDR_TARGET <= 0.05)] <- "sign. mouse and TARGET"

# calculate the coefficient correlation
fit <- lm(merged_table$logFC_mouse ~ merged_table$logFC_TARGET, data = merged_table)
coeff_corr <- round(cor(merged_table$logFC_TARGET, merged_table$logFC_mouse, method = "pearson"), 3)

# draw the scatterplot
scatterplot <- ggplot(merged_table, aes(x = logFC_TARGET, y = logFC_mouse, color = color_sign)) +
  geom_point(alpha = 0.75) + theme_bw() + xlab("log2FC TARGET") + ylab("log2FC mouse") +
  xlim(-max(abs(merged_table$logFC_TARGET)), max(abs(merged_table$logFC_TARGET))) +
  ylim(-max(abs(merged_table$logFC_mouse)), max(abs(merged_table$logFC_mouse))) +
  scale_color_manual(values = c("n.s" = "black", "sign. mouse" = "blue", 
                                "sign. TARGET" = "orange", "sign. mouse and TARGET" = "red")) +
  geom_text_repel(data = subset(merged_table[merged_table$external_gene_name %in% merged_table$external_gene_name[merged_table$color_sign == "sign. mouse and TARGET"],]), 
                  color = "black", aes(label = external_gene_name), size = 5) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between TARGET and GFP+ TPO vs PBS; corr coef (Pearson) = ", coeff_corr))


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
scatterplot <- ggplot(merged_table, aes(x = logFC_TARGET, y = logFC_mouse, color = color_sign)) +
  geom_point(alpha = 0.75) + theme_bw() + xlab("log2FC TARGET") + ylab("log2FC mouse") +
  xlim(-max(abs(merged_table$logFC_TARGET)), max(abs(merged_table$logFC_TARGET))) +
  ylim(-max(abs(merged_table$logFC_mouse)), max(abs(merged_table$logFC_mouse))) +
  scale_color_manual(values = c("n.s" = "black", "sign. mouse" = "blue", 
                                "sign. TARGET" = "orange", "sign. mouse and TARGET" = "red")) +
  geom_text_repel(data = subset(merged_table[merged_table$external_gene_name %in% list_relevant_genes,]), 
                  color = "black", aes(label = external_gene_name), size = 5) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between TARGET and GFP+ TPO vs PBS; corr coef (Pearson) = ", coeff_corr))


# save the scatterplot in files
saveFigures("scatterplot_GFP+_TPOvsPBS_relevantGenes", ggplot = scatterplot, dirPlot = outputPath, A4 = T)


# save the table in an excel file
write_xlsx(merged_table, path = "DEA/Human_patients/TARGET/allAML/TARGET_ERG_selection_kmeans/mergedTable_TARGET_GFP+_TPOvsPBS.xlsx")



# create the survival curve ---------------------

# create the output path
outdir_survival_TARGET <- create_dir(path = "DEA/Human_patients/TARGET/allAML/TARGET_ERG_selection_kmeans/Survival_curve")

# select the samples for the survival curve
samples_id <- read_xlsx(path = "DEA/Human_patients/TARGET/allAML/TARGET_ERG_selection_kmeans/selected_patients.xlsx", sheet = 1)

# select the patients corresponding to the sample_id
patients_id <- unlist(sapply(query_clinical_TARGET_2$submitter_id, function(x) if(any(grepl(pattern = x, samples_id$Id_sample)))return(x)))


## create the survival curve for selected genes ---------------------
# selection based on the comparison of DEG between the different databases

# create list of Up-Up list
UpUp_DEG <- c("IL12RB2", "INPP4B", "MECOM", "ADGRG6", "MMP28", "SEMA4F", "SH3BP5", "AIF1",
              "CDC14A", "FARP1", "GRAP2", "KRI1", "MLLT3",  "MOCOS", "OAS3","PBX1", "PDLIM2", "PDLIM5", "SERPINE2", "SLC6A12", 
              "STYK1", "TGM2")
DownDown_DEG <- c("AKAP12", "CACNA2D2", "CALR", "CPSF4", "CST7", "FKBP11", "FNDC3B", "MS4A3", 
              "MST1", "P2RY2", "P4HB", "PHKG1", "PMP22", "PRAG1", "SLC7A2", "TDRD9", "TRIM9", "UXS1", "ZNF704")




# draw survival curve for selected patients
for(gene in c(UpUp_DEG, DownDown_DEG)){
  
  # select the ensembl id of the gene
  gene_id_tmp <- rowData(data_Target)$gene_id[rowData(data_Target)$gene_name == gene] 
  
  # select the patients
  id_patients_tmp <- sapply(colnames(table_Target_tpm), function(x) paste(strsplit(x, "-")[[1]][1:3], collapse = "-"))
  
  # keep only patients with only one sample
  count_patients <- table(id_patients_tmp)
  kept_patients <- names(count_patients)[count_patients == 1]

  # select the samples 
  samples_id_tmp <- names(id_patients_tmp)[id_patients_tmp %in% intersect(kept_patients, query_clinical_TARGET_2$submitter_id)]
  
  # select value of genes expression
  exp <- log2(table_Target_tpm[gene_id_tmp, samples_id_tmp])
  
  # calculate the median
  median_exp <- median(exp)
  
  # select the samples
  above_sample <- names(exp)[exp >= median_exp]
  below_sample <- names(exp)[exp < median_exp]
  
  # select the patient
  above_patient <- unique(sapply(above_sample, function(x) paste(strsplit(x, "-")[[1]][1:3], collapse = "-")))
  below_patient <- unique(sapply(below_sample, function(x) paste(strsplit(x, "-")[[1]][1:3], collapse = "-")))
  
  # give information about the statistics
  print(paste0(gene, ": above median = ", length(above_patient), "; below median = ", length(below_patient)))
  
  # selection of the patients (# 80 different patients)
  query_clinical_TARGET_2_selection_tmp <- query_clinical_TARGET_2[query_clinical_TARGET_2$submitter_id %in% c(above_patient, below_patient), ]
  
  
  query_clinical_TARGET_2_selection_tmp$exp = "below median"
  query_clinical_TARGET_2_selection_tmp$exp[query_clinical_TARGET_2_selection_tmp$submitter_id %in% above_patient] = "above median"
  
  # change the days to last follow up
  query_clinical_TARGET_2_selection_tmp$days_to_last_follow_up[is.na(query_clinical_TARGET_2_selection_tmp$days_to_death)] = 2500
  
  # draw the survival curve
  TCGAanalyze_survival(
    data = query_clinical_TARGET_2_selection_tmp,
    clusterCol = "exp",
    main = paste0(gene, "; TARGET, Kaplan-Meier survival curve"),
    filename = file.path(outdir_survival_TARGET, paste0("survival_curve_", gene, "_allPatients.pdf")),
    color = c("#ED1C24", "#0F75BC"),
    conf.int = F,
    width = 8,
    xlim = c(0, 2000)
  )
  
}


# create the scatterplot for common DEGs compare to MECOM ---------------------

# get the id for MECOM
mecom_id = rowData(data_Target)$gene_id[rowData(data_Target)$gene_name == "MECOM"]

# get the id for IL12RB2
il12rb2_id = rowData(data_Target)$gene_id[rowData(data_Target)$gene_name == "IL12RB2"]

# get the id for INPP4B
INPP4B_id = rowData(data_Target)$gene_id[rowData(data_Target)$gene_name == "INPP4B"]

# create the dataframe for the scatterplot figure
dt_TARGET_ERG_MECOM <- data.frame("MECOM" = table_Target_tpm[mecom_id, colData_table$barcode],
                                "IL12RB2" = table_Target_tpm[il12rb2_id, colData_table$barcode],
                                "INPP4B" = table_Target_tpm[INPP4B_id, colData_table$barcode],
                                "Sample" = colData_table$barcode,
                                "Color" = "none")

# export the count table
saveRDS(dt_TARGET_ERG_MECOM, file = "DEA/Human_patients/TARGET/dt_TARGET_ERG_MECOM.rds")


# update the patients information
dt_TARGET_ERG_MECOM$Color[dt_TARGET_ERG_MECOM$Sample %in% colData_table$barcode[colData_table$patient %in% patients_allKMT2A$TARGET.USI]] = "MLLr"
dt_TARGET_ERG_MECOM$Color[dt_TARGET_ERG_MECOM$Sample %in% patients] = "MLL-AF9"
dt_TARGET_ERG_MECOM$Color[dt_TARGET_ERG_MECOM$Sample %in% colData_table$barcode[colData_table$patient %in% table_AML_allEVI1$TARGET.USI]] = "EVI1r"



# calculate the coefficient correlation
dt_TARGET_ERG_MECOM_tmp <- dt_TARGET_ERG_MECOM[dt_TARGET_ERG_MECOM$MECOM != 0, ]
dt_TARGET_ERG_MECOM_tmp <- dt_TARGET_ERG_MECOM_tmp[dt_TARGET_ERG_MECOM_tmp$IL12RB2 != 0, ]
fit <- lm(log2(dt_TARGET_ERG_MECOM_tmp$IL12RB2) ~ log2(dt_TARGET_ERG_MECOM_tmp$MECOM), data = dt_TARGET_ERG_MECOM_tmp)
coeff_corr <- round(cor(log2(dt_TARGET_ERG_MECOM_tmp$MECOM), log2(dt_TARGET_ERG_MECOM_tmp$IL12RB2), method = "pearson"), 3)

# create the scatterplot
scatterplot <- getScatterplot(data = dt_TARGET_ERG_MECOM_tmp, x = "MECOM", y = "IL12RB2", log2 = T)
scatterplot <- scatterplot + geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between MECOM and IL12RB2; corr coef (Pearson) = ", coeff_corr)) +
  scale_color_manual(values = c("none" = "purple", "MLLr" = "cyan3", 
                                "MLL-AF9" = "forestgreen", "EVI1r" = "red")) + geom_point(size=3)

# create the scatterplot for the Erg vs MECOM expression for all patients
saveFigures("scatterplot_IL12RB2vsMECOM", scatterplot, 
            dirPlot = "DEA/Human_patients/TARGET/allAML/", A4 = T)

# calculate the coefficient correlation
dt_TARGET_ERG_MECOM_tmp <- dt_TARGET_ERG_MECOM[dt_TARGET_ERG_MECOM$MECOM != 0, ]
dt_TARGET_ERG_MECOM_tmp <- dt_TARGET_ERG_MECOM_tmp[dt_TARGET_ERG_MECOM_tmp$INPP4B != 0, ]
fit <- lm(log2(dt_TARGET_ERG_MECOM_tmp$INPP4B) ~ log2(dt_TARGET_ERG_MECOM_tmp$MECOM), data = dt_TARGET_ERG_MECOM_tmp)
coeff_corr <- round(cor(log2(dt_TARGET_ERG_MECOM_tmp$MECOM), log2(dt_TARGET_ERG_MECOM_tmp$INPP4B), method = "pearson"), 3)

# create the scatterplot
scatterplot <- getScatterplot(data = dt_TARGET_ERG_MECOM_tmp, x = "MECOM", y = "INPP4B", log2 = T)
scatterplot <- scatterplot + geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between MECOM and INPP4B; corr coef (Pearson) = ", coeff_corr)) +
  scale_color_manual(values = c("none" = "purple", "MLLr" = "cyan3", 
                                "MLL-AF9" = "forestgreen", "EVI1r" = "red")) + geom_point(size=3)

# create the scatterplot for the Erg vs MECOM expression for all patients
saveFigures("scatterplot_INPP4BvsMECOM", scatterplot, 
            dirPlot = "DEA/Human_patients/TARGET/allAML/", A4 = T)


