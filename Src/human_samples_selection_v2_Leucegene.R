###### command line to analyse GEO dataset from leucegene database
# done by Jonathan Seguin, group of Prof. Schwaller, DBM, UKBB, Basel, Switzerland
# email: jonathan.seguin@unibas.ch, seguin.jonathan@gmail.com
# created in Thu Jun 22 09:13:41 2023


# set the environment -----------------------------------------------------

## load the libraries -----------------------------------------------------
library(GEOquery)
library(edgeR)
library(writexl)
library(readxl)
library(tidyverse)
library(xlsx)
library(ggrepel)
library(rgl)
library(factoextra)
library(gridExtra)

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
      ylab(paste0(y, " (log2(rpkm))")) + 
      xlab(paste0(x, " (log2(rpkm))")) +
      theme_bw() + 
      ggtitle(paste0("scatterplot, ", x, " vs ", y, ", LEUCEGENE database")) 
  } else {
    scatterplot =  ggplot(data = data, aes(x = data[, x], y = data[, y], color = Color)) + 
      geom_point() +
      ylab(paste0(y, " (rpkm)")) + 
      xlab(paste0(x, " (rpkm)")) + 
      theme_bw() + 
      ggtitle(paste0("scatterplot, ", x, " vs ", y, ", LEUCEGENE database")) 
  }
  
  
  # return the scatterplot
  return(scatterplot)
  
}


## read the tables ---------------------------------------------------------

# load the table from the TPO-Evi1 RNAseq
fits_tmp <- readRDS(file = "DEA/Rds/RDS_woD/edge_samplesR_fits_results_DEA.rds")

# select the orthologous genes
human_ortho <- getHumanOrthologous()

# select the mouse genes
Mouse_genes <- getMouseGene()

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

## create the output directory -----------------------------------------------------
leucegene_dir <- create_dir(path = "Dataset/Leucegene/")

DEA_leucegene <- create_dir(path = "DEA/Human_patients/Leucegene/AllTissues")

# download the data -----------------------------------------------------

# download the GSE file
gse_leuce <- getGEO("GSE67040", GSEMatrix = TRUE, destdir = leucegene_dir)

# download the supplemental file
filePaths = getGEOSuppFiles("GSE67040", baseDir = leucegene_dir)

# list the files
txt_files <- list.files(path = "Dataset/Leucegene/GSE67040/GSE67040_RAW/")

# read the tables
list_tables <- getContentFiles(list_files = txt_files, 
                               directory = "Dataset/Leucegene/GSE67040/GSE67040_RAW/", 
                               header = F)


# prepare the tables -----------------------------------------------------

# check if the MECOM is duplicated in the tables
for(filename in txt_files){
  
  # select the data
  tmp_table <-  table(list_tables[[filename]][, 4])
  
  if(tmp_table["MECOM"] > 1){
    print(filename)
  }
   
}
# no duplicated row for MECOM


# create the table for the RPKM
table_RPKM <- data.frame("Gene" = unique(list_tables$GSM1203305_02H053_RPKM.txt$V4))
for(filename in txt_files){
  
  # select the data
  tmp_table <-  list_tables[[filename]][, 4:5]
  
  # remove duplicated values
  tmp_table = tmp_table[!duplicated(tmp_table$V4),]
  
  # merge the data
  table_RPKM <- merge(x = table_RPKM, y = tmp_table, by.x = "Gene", by.y = "V4", all.x= T)
  
  # update the name
  colnames(table_RPKM)[grep(pattern = "V5", colnames(table_RPKM))] = gsub(pattern = "_RPKM.txt", replacement = "", filename)
  
}
#table_RPKM <- lapply(txt_files, function(x) cbind(table_RPKM, list_tables[[x]]$V5))
# update the table
rownames(table_RPKM) = table_RPKM$Gene
table_RPKM = table_RPKM[, -1]


# create the table for the raw
table_raw <- data.frame("Gene" = unique(list_tables$GSM1203305_02H053_RPKM.txt$V4))
for(filename in txt_files){
  
  # select the data
  tmp_table <-  list_tables[[filename]][, c(4,6)]
  
  # remove duplicated values
  tmp_table = tmp_table[!duplicated(tmp_table$V4),]
  
  # merge the data
  table_raw <- merge(x = table_raw, y = tmp_table, by.x = "Gene", by.y = "V4", all.x= T)
  
  # update the name
  colnames(table_raw)[grep(pattern = "V6", colnames(table_raw))] = gsub(pattern = "_RPKM.txt", replacement = "", filename)
  
}

# update the table
rownames(table_raw) = table_raw$Gene
table_raw = table_raw[, -1]


# load information for the patients
info_patients <- read.table(file = "Dataset/Leucegene/GSE67040_series_matrix.txt", 
                            quote = "", skip = 29, header = T, sep = "\t",
                            fill = T, check.names = F)


# display the MECOM expression -----------------------------------------------------

# select the Mecom id
mecom_count <- table_RPKM["MECOM",]
mecom_log2rpkm <- log2(mecom_count+1)
mecom_log2rpkm <- unlist(mecom_log2rpkm)


# intialize the cut-off values
cutoffs <- 2:4
for(cutoff in cutoffs){
  
  pdf(paste0("Human_comparison/Plots/barplot_Leucegene_rpkm_log2_cutoff_", cutoff, ".pdf"))
  barplot(sort(mecom_log2rpkm, decreasing = T)) 
  abline(h = cutoff)
  dev.off()
  
}


# perform the DEA -----------------------------------------------------

## for all patients -----------------------------------------------------

# read the patients file
leucegene_info <- read_xlsx(path = "Dataset/Leucegene/Leucegene_KMT2A-MLLT3-2023.07.14_final.xlsx", sheet = 1)

# perform statistics based on patient
table(leucegene_info$Tissue)
#       Blood Bone marrow 
#         19          31 
table(leucegene_info$Tissue, leucegene_info$`Cytogenetic information`)
#             EVI1 rearrangements (+EVI1 FISH positive) (Irrespective of additional cytogenetic abnormalities)
# Blood                                                                                                      6
# Bone marrow                                                                                                6
# MLL translocations (+MLL FISH positive) (Irrespective of additional cytogenetic abnormalities)
# Blood                                                                                                   13
# Bone marrow                                                                                             25

# select id_raw
id_raw <- unlist(sapply(colnames(table_raw), function(x) strsplit(x, split = "_")[[1]][2]))

# keep patients with EVI1 rearrangement
table_raw_EVI1 <- table_raw[, id_raw %in% leucegene_info$Sample_id[grepl(pattern = "EVI1", leucegene_info$`Cytogenetic information`)]]

# update id_raw with EVI1 rearrangement
id_raw_EVI1 <- unlist(sapply(colnames(table_raw_EVI1), function(x) strsplit(x, split = "_")[[1]][2]))

# keep patients with MLL rearrangement
table_raw_MLL <- table_raw[, id_raw %in% leucegene_info$Sample_id[grepl(pattern = "MLL", leucegene_info$`Cytogenetic information`)]]

# update id_raw with MLL rearrangement
id_raw_MLL <- unlist(sapply(colnames(table_raw_MLL), function(x) strsplit(x, split = "_")[[1]][2]))


# keep patients with MF9 rearrangement
table_raw_MF9 <- table_raw[, id_raw %in% leucegene_info$Sample_id[grepl(pattern = "1", leucegene_info$`KMT2A::MLLT3 fusion`)]]

# update id_raw with MLL rearrangement
id_raw_MF9 <- unlist(sapply(colnames(table_raw_MF9), function(x) strsplit(x, split = "_")[[1]][2]))



### for all tissues -----------------------------------------------------

# remove patients with EVI1 rearrangement
table_raw_woEVI1 <- table_raw[, !id_raw %in% leucegene_info$Sample_id[grepl(pattern = "EVI1", leucegene_info$`Cytogenetic information`)]]
mecom_log2rpkm_woEVI1 <- mecom_log2rpkm[!id_raw %in% leucegene_info$Sample_id[grepl(pattern = "EVI1", leucegene_info$`Cytogenetic information`)]]

# update id_raw
id_raw_woEVI1 <- unlist(sapply(colnames(table_raw_woEVI1), function(x) strsplit(x, split = "_")[[1]][2]))


### for bone marrows tissues -----------------------------------------------------

# prepare the human genes
Human_genes <- getHumanGene()

# create the output directory
DEA_leucegene_BM <- create_dir(path = "DEA/Human_patients/Leucegene/Bone_Marrow/allAML")

# select only bone marrow samples
id_BM <- names(info_patients[9,grep(pattern = "Bone marrow", info_patients[9,])])
id_BM = gsub(pattern = "\"", replacement = "", id_BM)
all_samples_selected_BM = colnames(table_raw)[id_raw %in% id_BM]

table_LEUCEGENE_rpkm <- table_RPKM[, id_raw %in% id_BM]



# perform DEA based on k-means clustering result ---------------------

## find the cluster of samples ---------------------

# create the dataframe for the scatterplot figure
dt_LEUCEGENE_ERG_MECOM_log2rpkm <- data.frame("MECOM" = unlist(log2(table_LEUCEGENE_rpkm["MECOM", ])),
                                        "ERG" = unlist(log2(table_LEUCEGENE_rpkm["ERG", ])),
                                        "Sample" = colnames(table_LEUCEGENE_rpkm),
                                        "Color" = "none")

# update the patients information
dt_LEUCEGENE_ERG_MECOM_log2rpkm$Color[dt_LEUCEGENE_ERG_MECOM_log2rpkm$Sample %in% names(id_raw_MLL)] = "MLLr"
dt_LEUCEGENE_ERG_MECOM_log2rpkm$Color[dt_LEUCEGENE_ERG_MECOM_log2rpkm$Sample %in% names(id_raw_MF9)] = "MLL-AF9"
dt_LEUCEGENE_ERG_MECOM_log2rpkm$Color[dt_LEUCEGENE_ERG_MECOM_log2rpkm$Sample %in% names(id_raw_EVI1)] = "EVI1r"

# remove infinite values
dt_LEUCEGENE_ERG_MECOM_log2rpkm = dt_LEUCEGENE_ERG_MECOM_log2rpkm[is.finite(dt_LEUCEGENE_ERG_MECOM_log2rpkm$MECOM), ]
dt_LEUCEGENE_ERG_MECOM_log2rpkm = dt_LEUCEGENE_ERG_MECOM_log2rpkm[is.finite(dt_LEUCEGENE_ERG_MECOM_log2rpkm$ERG), ]


# initialize list to save kmeans results and their corresponding plots
list_kmeans <- NULL
list_kmean_plot <- NULL


# calculate the different k-means clusters
for(k in 2:9){
  
  # calculate the temporary clustering
  tmp_kmeans <-  kmeans(dt_LEUCEGENE_ERG_MECOM_log2rpkm[, 1:2], centers = k, nstart = 25)
  
  # save the result in the list
  list_kmeans <- c(list_kmeans, list(tmp_kmeans))
  
  # save the plot
  list_kmean_plot <- c(list_kmean_plot, 
                       list(fviz_cluster(tmp_kmeans, data = dt_LEUCEGENE_ERG_MECOM_log2rpkm[, 1:2], geom = "point", ggtheme = theme_bw()) + ggtitle(paste0("k = ", k))))
  
}

# renames the list
names(list_kmeans) <- paste0("k", 2:9)
names(list_kmean_plot) <- paste0("k", 2:9)

# save k-means object in rds file
saveRDS(list_kmeans, file = "DEA/Human_patients/Leucegene/Bone_Marrow/allAML/list_kmeans.rds")
saveRDS(list_kmean_plot, file = "DEA/Human_patients/Leucegene/Bone_Marrow/allAML/list_kmean_plot.rds")

# display the result
g <- grid.arrange(grobs = list_kmean_plot, ncol = 4, nrow = 2)

# save in pdf file
ggsave(filename = "DEA/Human_patients/Leucegene/Bone_Marrow/allAML/kmeans_clustering_ERGvsMECOM_LEUCEGENE.pdf",
       plot = g, width =297, height = 210, units = "mm")

# determine the optiimal number of cluster
set.seed(123)
gap_stat <- clusGap(dt_LEUCEGENE_ERG_MECOM_log2rpkm[, 1:2], FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)

# save the result in the pdf
pdf(file = "DEA/Human_patients/Leucegene/Bone_Marrow/allAML/bestClusterNumber.pdf", paper = "a4r")
set.seed(123)
fviz_nbclust(dt_LEUCEGENE_ERG_MECOM_log2rpkm[, 1:2], kmeans, method = "wss")
set.seed(123)
fviz_nbclust(dt_LEUCEGENE_ERG_MECOM_log2rpkm[, 1:2], kmeans, method = "silhouette")
fviz_gap_stat(gap_stat)
dev.off()

# selection of cutoff 4 for the clustering result
# update the scatterplot 
dt_LEUCEGENE_ERG_MECOM_log2rpkm$cluster = list_kmeans$k4$cluster


## perform the DEA based on clusters (4 vs 2) ---------------------

# get the samples information
cluster2_info <- dt_LEUCEGENE_ERG_MECOM_log2rpkm[dt_LEUCEGENE_ERG_MECOM_log2rpkm$cluster == "2", ]

# select the minimal cutoff for ERG
min_cutoff_ERG <- min(cluster2_info$ERG)

# filter the patient based on the ERG expression
patient_kept <- dt_LEUCEGENE_ERG_MECOM_log2rpkm$Sample[dt_LEUCEGENE_ERG_MECOM_log2rpkm$ERG >= min_cutoff_ERG]

# select patient high (cluster 4 in figure paper)
patient_high <- intersect(dt_LEUCEGENE_ERG_MECOM_log2rpkm$Sample[dt_LEUCEGENE_ERG_MECOM_log2rpkm$cluster == "4"], 
                          patient_kept)

# select patient low (cluster 1 in figure paper)
patient_low <- intersect(dt_LEUCEGENE_ERG_MECOM_log2rpkm$Sample[dt_LEUCEGENE_ERG_MECOM_log2rpkm$cluster == "2"], 
                         patient_kept)

# create output directory
tmp_DEA <- create_dir(path = file.path(DEA_leucegene_BM, "LEUCEGENE_ERG_selection_kmeans"))

# select the samples
table_raw_BM_sel <- table_raw_BM[, colnames(table_raw_BM) %in% c(patient_high, patient_low)]

# replace NA value in table
table_raw_BM_sel[is.na(table_raw_BM_sel)] = 0

# determine the group
group <- rep("low", ncol(table_raw_BM_sel))

# update the group
group[colnames(table_raw_BM_sel) %in% patient_high] = "high"

# remove the gene lowly expressed
print(paste0("# genes: ", nrow(table_raw_BM_sel)))
# keep gene with more than 10 reads in total
count_10 <- nrow(table_raw_BM_sel[rowSums(table_raw_BM_sel) > 10,])
print(paste0("# kept genes (filter based on 10): ", count_10))
# 21138 kept genes among 21865 genes

# create the dge
dge <- DGEList(counts = table_raw_BM_sel, group = group)

# filter according to the filtering gene expression
keep <- filterByExpr(dge)
count_filter <- length(keep[keep])
print(paste0("# kept gene after filtering by function: ", count_filter))
# 17763 genes kept

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
allresult <-  cbind(allresult, "external_gene_name" = rownames(allresult))

# merge with the human gene information
allresult <- merge(x = allresult, 
                   y = Human_genes, 
                   by = "external_gene_name")

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
write.xlsx2(data.frame("Id_sample" = colnames(table_raw_BM_sel), "Group" = group), file = file.path(tmp_DEA, "selected_patients.xlsx"),
            col.names = T, row.names = F)
gc()



### perform GSEA and ORA ---------------------


# get the output path
outputPath <- tmp_DEA

# load the table with DEA result for the cutoff of 0.8
MLLr_LEUCEGENE <- read_xlsx(path = file.path(outputPath, "results_DEA.xlsx"), 
                       sheet = 1)


#### perform GSEA -----------------------------------------------------------

# create the output path for the GSEA
outputPath_GSEA <- create_dir(path = file.path(outputPath, "GSEA"))

# select the FC
subtable <- MLLr_LEUCEGENE[!is.na(MLLr_LEUCEGENE$entrezgene_id),]
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
write_xlsx(list_fGSEA_results, path = file.path(outputPath_GSEA, "GSEA_Msigdatabase_MLLr_LEUCEGENE_kmeans.xlsx"))


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
  filename = paste0("GSEA_", filename, "_MLLr_LEUCEGENE_kmeans")
  
  # save the enrichment in the pdf
  saveFigures(fileName = filename, ggplot = plotEnrichment(custom_pathways[[enrich]], geneList_GSEA) + 
                labs(title = enrich), dirPlot = outputPath_GSEA, A4 = T)
  
}

# write the excel file
write_xlsx(fgseaRes_custom_pathways, path = file.path(outputPath_GSEA, "GSEA_custom_pathway_Msigdatabase_MLLr_LEUCEGENE_kmeans.xlsx"))



#### perform ORA  -----------------------------------------------------------

# create the directory
outputPath_ORA <- create_dir(path = file.path(outputPath, "ORA"))

# select all tested genes
genes <- MLLr_LEUCEGENE$ensembl_gene_id

# select the significant toptags
dge_DEG <- MLLr_LEUCEGENE[MLLr_LEUCEGENE$FDR <= 0.05, ]

# create the output directory
tmp_dir <- create_dir(path = file.path(outputPath_ORA, "DEG_woUniverse"))

# perform the enrichment analysis
EnrichmentTerm_Analysis(common = dge_DEG$ensembl_gene_id, outputDir = tmp_dir, 
                        mart = ensembl, isMouse = F)
# create the output directory
tmp_dir <- create_dir(path = file.path(outputPath_ORA, "DEG"))

# perform the enrichment analysis
EnrichmentTerm_Analysis(common = dge_DEG$ensembl_gene_id, outputDir = tmp_dir, 
                        mart = ensembl, isMouse = F, universe = genes)




### compare with the mouse DEG within a scatterplot ---------------------

# get the output path
outputPath <- tmp_DEA

# load the table with DEA result for the cutoff of 0.8
MLLr_LEUCEGENE_cutoff_ERG <- read_xlsx(path = file.path(outputPath, "results_DEA.xlsx"), 
                                  sheet = 1)

# merge the tables
merged_table <- merge(x = MLLr_LEUCEGENE_cutoff_ERG,
                      y = fit_bulkTPOVsnoTPO,
                      by = "external_gene_name",
                      suffixes = c("_LEUCEGENE", "_mouse"))

# create the color for the significativity 
merged_table$color_sign <- "n.s"

# update the color based on the significativity
merged_table$color_sign[merged_table$FDR_mouse <= 0.05] <- "sign. mouse"
merged_table$color_sign[merged_table$FDR_LEUCEGENE <= 0.05] <- "sign. LEUCEGENE"
merged_table$color_sign[(merged_table$FDR_mouse <= 0.05) & (merged_table$FDR_LEUCEGENE <= 0.05)] <- "sign. mouse and LEUCEGENE"

# calculate the coefficient correlation
fit <- lm(merged_table$logFC_LEUCEGENE ~ merged_table$logFC_mouse, data = merged_table)
coeff_corr <- round(cor(merged_table$logFC_LEUCEGENE, merged_table$logFC_mouse, method = "pearson"), 3)

# draw the scatterplot
scatterplot <- ggplot(merged_table, aes(x = logFC_LEUCEGENE, y = logFC_mouse, color = color_sign)) +
  geom_point(alpha = 0.75) + theme_bw() + xlab("log2FC LEUCEGENE") + ylab("log2FC mouse") +
  xlim(-max(abs(merged_table$logFC_LEUCEGENE)), max(abs(merged_table$logFC_LEUCEGENE))) +
  ylim(-max(abs(merged_table$logFC_mouse)), max(abs(merged_table$logFC_mouse))) +
  scale_color_manual(values = c("n.s" = "black", "sign. mouse" = "blue", 
                                "sign. LEUCEGENE" = "orange", "sign. mouse and LEUCEGENE" = "red")) +
  geom_text_repel(data = subset(merged_table[merged_table$external_gene_name %in% merged_table$external_gene_name[merged_table$color_sign == "sign. mouse and LEUCEGENE"],]), 
                  color = "black", aes(label = external_gene_name), size = 5) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between LEUCEGENE and bulk TPO vs PBS; corr coef (Pearson) = ", coeff_corr))


# save the scatterplot in files
saveFigures("scatterplot_bulk_TPOvsPBS", ggplot = scatterplot, dirPlot = outputPath, A4 = T)




# merge the tables
merged_table <- merge(x = MLLr_LEUCEGENE_cutoff_ERG,
                      y = fit_GFP_TPOVsnoTPO,
                      by = "external_gene_name",
                      suffixes = c("_LEUCEGENE", "_mouse"))


# create the color for the significativity 
merged_table$color_sign <- "n.s"

# update the color based on the significativity
merged_table$color_sign[merged_table$FDR_mouse <= 0.05] <- "sign. mouse"
merged_table$color_sign[merged_table$FDR_LEUCEGENE <= 0.05] <- "sign. LEUCEGENE"
merged_table$color_sign[(merged_table$FDR_mouse <= 0.05) & (merged_table$FDR_LEUCEGENE <= 0.05)] <- "sign. mouse and LEUCEGENE"

# calculate the coefficient correlation
fit <- lm(merged_table$logFC_mouse ~ merged_table$logFC_LEUCEGENE, data = merged_table)
coeff_corr <- round(cor(merged_table$logFC_LEUCEGENE, merged_table$logFC_mouse, method = "pearson"), 3)

# draw the scatterplot
scatterplot <- ggplot(merged_table, aes(x = logFC_LEUCEGENE, y = logFC_mouse, color = color_sign)) +
  geom_point(alpha = 0.75) + theme_bw() + xlab("log2FC LEUCEGENE") + ylab("log2FC mouse") +
  xlim(-max(abs(merged_table$logFC_LEUCEGENE)), max(abs(merged_table$logFC_LEUCEGENE))) +
  ylim(-max(abs(merged_table$logFC_mouse)), max(abs(merged_table$logFC_mouse))) +
  scale_color_manual(values = c("n.s" = "black", "sign. mouse" = "blue", 
                                "sign. LEUCEGENE" = "orange", "sign. mouse and LEUCEGENE" = "red")) +
  geom_text_repel(data = subset(merged_table[merged_table$external_gene_name %in% merged_table$external_gene_name[merged_table$color_sign == "sign. mouse and LEUCEGENE"],]), 
                  color = "black", aes(label = external_gene_name), size = 5) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between LEUCEGENE and GFP+ TPO vs PBS; corr coef (Pearson) = ", coeff_corr))


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
scatterplot <- ggplot(merged_table, aes(x = logFC_LEUCEGENE, y = logFC_mouse, color = color_sign)) +
  geom_point(alpha = 0.75) + theme_bw() + xlab("log2FC LEUCEGENE") + ylab("log2FC mouse") +
  xlim(-max(abs(merged_table$logFC_LEUCEGENE)), max(abs(merged_table$logFC_LEUCEGENE))) +
  ylim(-max(abs(merged_table$logFC_mouse)), max(abs(merged_table$logFC_mouse))) +
  scale_color_manual(values = c("n.s" = "black", "sign. mouse" = "blue", 
                                "sign. LEUCEGENE" = "orange", "sign. mouse and LEUCEGENE" = "red")) +
  geom_text_repel(data = subset(merged_table[merged_table$external_gene_name %in% list_relevant_genes,]), 
                  color = "black", aes(label = external_gene_name), size = 5) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between LEUCEGENE and GFP+ TPO vs PBS; corr coef (Pearson) = ", coeff_corr))


# save the scatterplot in files
saveFigures("scatterplot_GFP+_TPOvsPBS_relevantGenes", ggplot = scatterplot, dirPlot = outputPath, A4 = T)

# save the table in an excel file
write_xlsx(merged_table, path = "DEA/Human_patients/Leucegene/Bone_Marrow/allAML/LEUCEGENE_ERG_selection_kmeans/mergedTable_LEUCEGENE_GFP+_TPOvsPBS.xlsx")



# create the scatterplot for common DEGs compare to MECOM ---------------------



# create the dataframe for the scatterplot figure
dt_LEUCEGENE_ERG_MECOM <- data.frame("MECOM" = unlist(table_LEUCEGENE_rpkm["MECOM",]),
                                     "IL12RB2" = unlist(table_LEUCEGENE_rpkm["IL12RB2", ]),
                                     "Sample" = colnames(table_LEUCEGENE_rpkm),
                                     "Color" = "none")


# update the patients information
dt_LEUCEGENE_ERG_MECOM$Color[dt_LEUCEGENE_ERG_MECOM$Sample %in% names(id_raw_EVI1)] = "EVI1r"
dt_LEUCEGENE_ERG_MECOM$Color[dt_LEUCEGENE_ERG_MECOM$Sample %in% names(id_raw_MLL)] = "MLLr"
dt_LEUCEGENE_ERG_MECOM$Color[dt_LEUCEGENE_ERG_MECOM$Sample %in% names(id_raw_MF9)] = "MLL-AF9"


# save the rds file
saveRDS(dt_LEUCEGENE_ERG_MECOM, file = "DEA/Human_patients/Leucegene/Bone_Marrow/dt_LEUCEGENE_ERG_MECOM.rds")


# calculate the coefficient correlation
dt_LEUCEGENE_ERG_MECOM_tmp <- dt_LEUCEGENE_ERG_MECOM[dt_LEUCEGENE_ERG_MECOM$MECOM != 0, ]
dt_LEUCEGENE_ERG_MECOM_tmp <- dt_LEUCEGENE_ERG_MECOM_tmp[dt_LEUCEGENE_ERG_MECOM_tmp$IL12RB2 != 0, ]
fit <- lm(log2(dt_LEUCEGENE_ERG_MECOM_tmp$IL12RB2) ~ log2(dt_LEUCEGENE_ERG_MECOM_tmp$MECOM), data = dt_LEUCEGENE_ERG_MECOM_tmp)
coeff_corr <- round(cor(log2(dt_LEUCEGENE_ERG_MECOM_tmp$MECOM), log2(dt_LEUCEGENE_ERG_MECOM_tmp$IL12RB2), method = "pearson"), 3)

# create the scatterplot
scatterplot <- getScatterplot(data = dt_LEUCEGENE_ERG_MECOM_tmp, x = "MECOM", y = "IL12RB2", log2 = T)
scatterplot <- scatterplot + geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between MECOM and IL12RB2; corr coef (Pearson) = ", coeff_corr)) +
  scale_color_manual(values = c("none" = "purple", "MLLr" = "cyan3", 
                                "MLL-AF9" = "forestgreen", "EVI1r" = "red")) + geom_point(size=3)


# create the scatterplot for the Erg vs MECOM expression for all patients
saveFigures("scatterplot_IL12RB2vsMECOM", scatterplot, 
            dirPlot = "DEA/Human_patients/Leucegene/Bone_Marrow/allAML/", A4 = T)

# calculate the coefficient correlation
dt_LEUCEGENE_ERG_MECOM_tmp <- dt_LEUCEGENE_ERG_MECOM[dt_LEUCEGENE_ERG_MECOM$MECOM != 0, ]
dt_LEUCEGENE_ERG_MECOM_tmp <- dt_LEUCEGENE_ERG_MECOM_tmp[dt_LEUCEGENE_ERG_MECOM_tmp$INPP4B != 0, ]
fit <- lm(log2(dt_LEUCEGENE_ERG_MECOM_tmp$INPP4B) ~ log2(dt_LEUCEGENE_ERG_MECOM_tmp$MECOM), data = dt_LEUCEGENE_ERG_MECOM_tmp)
coeff_corr <- round(cor(log2(dt_LEUCEGENE_ERG_MECOM_tmp$MECOM), log2(dt_LEUCEGENE_ERG_MECOM_tmp$INPP4B), method = "pearson"), 3)

# create the scatterplot
scatterplot <- getScatterplot(data = dt_LEUCEGENE_ERG_MECOM_tmp, x = "MECOM", y = "INPP4B", log2 = T)
scatterplot <- scatterplot + geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between MECOM and INPP4B; corr coef (Pearson) = ", coeff_corr)) +
  scale_color_manual(values = c("none" = "purple", "MLLr" = "cyan3", 
                                "MLL-AF9" = "forestgreen", "EVI1r" = "red")) + geom_point(size=3)


# create the scatterplot for the Erg vs MECOM expression for all patients
saveFigures("scatterplot_INPP4BvsMECOM", scatterplot, 
            dirPlot = "DEA/Human_patients/Leucegene/Bone_Marrow/allAML/", A4 = T)
