##### command lines to select the human samples, 1st version, only for St-Jude selection based on Memphis's paper
# done by Jonathan Seguin, group of Prof. Schwaller, DBM, UKBB, Basel, Switzerland
# email: jonathan.seguin@unibas.ch, seguin.jonathan@gmail.com
# Thu Feb 15 13:56:41 2024

# prepare the working environment -----------------------------------------------------

## load the libraries -----------------------------------------------------
library(edgeR)
library(readxl)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(factoextra)
library(gridExtra)
library(cluster)
library(scales)
library(TCGAbiolinks)
edgeRUsersGuide()

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
      ggtitle(paste0("scatterplot, ", x, " vs ", y, ", ST-JUDE database")) 
  } else {
    scatterplot =  ggplot(data = data, aes(x = data[, x], y = data[, y], color = Color)) + 
      geom_point() +
      ylab(paste0(y, " (rpkm)")) + 
      xlab(paste0(x, " (rpkm)")) + 
      theme_bw() + 
      ggtitle(paste0("scatterplot, ", x, " vs ", y, ", ST-JUDE database")) 
  }
  
  
  # return the scatterplot
  return(scatterplot)
  
}

## load the objects -----------------------------------------------------

# read the raw data
data_stJude_Memphis <- readRDS(file = "Dataset/St-Jude/raw_count_Memphis.rds")

# read the samples information
patients_information_Memphis <- read_xlsx(path = "Dataset/St-Jude/patient_info_Memphis_StJude.xlsx")

# select the orthologous genes
ortho <- getHumanOrthologous()

# get all genes information for human
Human_genes <- getHumanGene()

# select all information about genes in human databases
# download ensembl version based on the paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5687824/
# source: https://www.biostars.org/p/83901/
exons.list.per.gene <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,by="gene")
exonic.gene.sizes <- sum(width(GenomicRanges::reduce(exons.list.per.gene)))
Human_genes_red <- unique(Human_genes[Human_genes$entrezgene_id %in% names(exonic.gene.sizes), 1:3])
Human_genes_red$length  = NA
# change as character the gene length
Human_genes_red$entrezgene_id = as.character(Human_genes_red$entrezgene_id)
Human_genes_red$length <- exonic.gene.sizes[Human_genes_red$entrezgene_id]


# calculate the tpm of genes expression ---------------------

# update the rownames for the genes
rownames(data_stJude_Memphis) = data_stJude_Memphis$Gene

# remove the gene colnames
data_stJude_Memphis = data_stJude_Memphis[, -1]

# count the number of samples
ncol(data_stJude_Memphis)
# 349

# count the number of patients
length(unique(sapply(colnames(data_stJude_Memphis), function(x)strsplit(x, split = "_")[[1]][1])))
# 349

# add the rownames information
rownames(patients_information_Memphis) = patients_information_Memphis$Patient

# reorder the table for patients
patients_information_Memphis = patients_information_Memphis[colnames(data_stJude_Memphis), ]

# select only genes where the exons size is known
data_stJude_Memphis_comGenes = data_stJude_Memphis[rownames(data_stJude_Memphis) %in% Human_genes_red$external_gene_name, ]

# create a dataframe for the genes
dt_Genes <- Human_genes_red[Human_genes_red$external_gene_name %in% 
                  rownames(data_stJude_Memphis_comGenes), ] %>% 
  dplyr::arrange(external_gene_name) %>%
  dplyr::select(external_gene_name, length)
dt_Genes = unique(dt_Genes)

# remove the duplicates
dt_Genes = dt_Genes[!duplicated(dt_Genes$external_gene_name), ]

# create the DGEList
dge_Memphis <- DGEList(counts = data_stJude_Memphis_comGenes, samples = patients_information_Memphis, genes = dt_Genes)

# calculate the rpkm counts
rpkm_dge_Memphis <- rpkm(y = dge_Memphis) # need to provide gene.length information

# count the number of patients without MECOM expression
table(data_stJude_Memphis_comGenes["MECOM", ]==0)
# FALSE  TRUE
#   270    79


# perform DEA based on k-means clustering result ---------------------

## find the cluster of samples ---------------------

# create the dataframe for the scatterplot figure
dt_Memphis_ERG_MECOM_log2rpkm <- data.frame("MECOM" = log2(rpkm_dge_Memphis["MECOM", ]),
                                          "ERG" = log2(rpkm_dge_Memphis["ERG",]),
                                          "Sample" = colnames(rpkm_dge_Memphis),
                                          "Color" = "none")

# update the patients information
dt_Memphis_ERG_MECOM_log2rpkm$Color[dt_Memphis_ERG_MECOM_log2rpkm$Sample %in% patients_information_Memphis$Patient[grep(pattern = "KMT2Ar", patients_information_Memphis$`Molecular category`)]] = "MLLr"
dt_Memphis_ERG_MECOM_log2rpkm$Color[dt_Memphis_ERG_MECOM_log2rpkm$Sample %in% patients_information_Memphis$Patient[grep(pattern = "KMT2A::MLLT3", patients_information_Memphis$`Defining alteration`)]] = "MLL-AF9"
dt_Memphis_ERG_MECOM_log2rpkm$Color[dt_Memphis_ERG_MECOM_log2rpkm$Sample %in% patients_information_Memphis$Patient[grep(pattern = "MECOM", patients_information_Memphis$`Molecular category`)]] = "EVI1r"

# remove infinite values
dt_Memphis_ERG_MECOM_log2rpkm = dt_Memphis_ERG_MECOM_log2rpkm[is.finite(dt_Memphis_ERG_MECOM_log2rpkm$MECOM) & is.finite(dt_Memphis_ERG_MECOM_log2rpkm$ERG), ]

# initialize list to save kmeans results and their corresponding plots
list_kmeans <- NULL
list_kmean_plot <- NULL

# calculate the different k-means clusters
for(k in 2:9){
  
  # calculate the temporary clustering
  tmp_kmeans <-  kmeans(dt_Memphis_ERG_MECOM_log2rpkm[, 1:2], centers = k, nstart = 25)
  
  # save the result in the list
  list_kmeans <- c(list_kmeans, list(tmp_kmeans))
  
  # save the plot
  list_kmean_plot <- c(list_kmean_plot, 
                       list(fviz_cluster(tmp_kmeans, data = dt_Memphis_ERG_MECOM_log2rpkm[, 1:2], geom = "point", ggtheme = theme_bw()) + ggtitle(paste0("k = ", k))))
  
}

# renames the list
names(list_kmeans) <- paste0("k", 2:9)
names(list_kmean_plot) <- paste0("k", 2:9)

# create the output directory
output_dir_Memphis <- create_dir(path = "DEA/Human_patients/St-Jude/allAML/")

# save in rds file
saveRDS(list_kmeans, file = file.path(output_dir_Memphis, "list_kmeans.rds"))
saveRDS(list_kmean_plot, file = file.path(output_dir_Memphis, "list_kmeans_plot.rds"))

# display the result
g <- grid.arrange(grobs = list_kmean_plot, ncol = 4, nrow = 2)


# save in pdf file
ggsave(filename = "DEA/Human_patients/St-Jude/allAML/kmeans_clustering_ERGvsMECOM_Memphis.pdf",
       plot = g, width =297, height = 210, units = "mm")

# determine the optimal number of cluster
set.seed(123)
gap_stat <- clusGap(dt_Memphis_ERG_MECOM_log2rpkm[, 1:2], FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)


# save the result in the pdf
pdf(file = "DEA/Human_patients/St-Jude/allAML/bestClusterNumber.pdf", paper = "a4r")
set.seed(123)
fviz_nbclust(dt_Memphis_ERG_MECOM_log2rpkm[, 1:2], kmeans, method = "wss")
set.seed(123)
fviz_nbclust(dt_Memphis_ERG_MECOM_log2rpkm[, 1:2], kmeans, method = "silhouette")
fviz_gap_stat(gap_stat)
dev.off()

# cluster selection based on k = 5
dt_Memphis_ERG_MECOM_log2rpkm$cluster = list_kmeans$k5$cluster



## perform the DEA based on clusters (2 vs 1) for k=5 based on k-means analysis ---------------------


# get the samples information
cluster2_info <- dt_Memphis_ERG_MECOM_log2rpkm[dt_Memphis_ERG_MECOM_log2rpkm$cluster == "2", ]

# select the minimal cutoff for ERG
min_cutoff_ERG <- min(cluster5_info$ERG)

# filter the patient based on the ERG expression
patient_kept <- dt_Memphis_ERG_MECOM_log2rpkm$Sample[dt_Memphis_ERG_MECOM_log2rpkm$ERG >= -1.5] #min_cutoff_ERG

# select patient high (cluster 5 in figure paper)
patient_high <- intersect(dt_Memphis_ERG_MECOM_log2rpkm$Sample[dt_Memphis_ERG_MECOM_log2rpkm$cluster == "2"], 
                          patient_kept)

# select patient low (cluster 1 in figure paper)
patient_low <- intersect(dt_Memphis_ERG_MECOM_log2rpkm$Sample[dt_Memphis_ERG_MECOM_log2rpkm$cluster == "1"], 
                         patient_kept)

# create the output directory
output_dir_Memphis <- create_dir(path = "DEA/Human_patients/St-Jude/allAML/St-Jude_ERG_selection_kmeans_k5/")

# select the samples
data_stJude_Memphis_sel <- data_stJude_Memphis[, colnames(data_stJude_Memphis) %in% c(patient_high, patient_low)]

# replace NA value in table
data_stJude_Memphis_sel[is.na(data_stJude_Memphis_sel)] = 0

# determine the group
group <- rep("low", ncol(data_stJude_Memphis_sel))

# update the group
group[colnames(data_stJude_Memphis_sel) %in% patient_high] = "high"

# remove the gene lowly expressed
print(paste0("# genes: ", nrow(data_stJude_Memphis_sel)))
# keep gene with more than 10 reads in total
count_10 <- nrow(data_stJude_Memphis_sel[rowSums(data_stJude_Memphis_sel) > 10,])
print(paste0("# kept genes (filter based on 10): ", count_10))
# 56060 kept genes among 59050 genes

# create the dge
dge <- DGEList(counts = data_stJude_Memphis_sel, group = group)

# filter according to the filtering gene expression
keep <- filterByExpr(dge)
count_filter <- length(keep[keep])
print(paste0("# kept gene after filtering by function: ", count_filter))
# 22004 genes kept

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
saveFigures("BCVplot", ggplot = plotBCV(dge), dirPlot = output_dir_Memphis, A4 = T)

# perform the exact test
et <-  exactTest(dge)

# save the et result in rds file
saveRDS(et, file = file.path(output_dir_Memphis, "et.rds"))

# provide all results
allresult <- as.data.frame(topTags(et, n=nrow(et)))

# add ensembl gene id
allresult <-  cbind(allresult, "official_gene_symbol" = rownames(allresult))

# merge with the human gene information
allresult <- merge(x = allresult, 
                   y = unique(Human_genes[, c(1,3)]), 
                   by.x = "official_gene_symbol", 
                   by.y = "external_gene_name")

# correct the orientation of the FC
if(allresult$logFC[grepl(pattern = "^MECOM$", allresult$official_gene_symbol)] < 0){
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
write_xlsx(x = allresult, path = file.path(output_dir_Memphis, "results_DEA.xlsx"), col_names = T)
gc()

# save the selected samples
write.xlsx2(data.frame("Id_sample" = colnames(data_stJude_Memphis_sel), "Group" = group), file = file.path(output_dir_Memphis, "selected_patients.xlsx"),
            col.names = T, row.names = F)
gc()

# check the quality of DEA
pdf(file = file.path(output_dir_Memphis, "hist_FDR.pdf"))
hist(allresult$FDR)
dev.off()



# comparison with DEG from GFP+ TPO vs PBS ---------------------

# load the table from the TPO-GFPp vs PBS-GFPn
fits_tmp <- readRDS(file = "DEA/Rds/RDS_woD/edge_samplesR_fits_results_DEA.rds")

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



# load the table with DEA result based on k-means
MLLr_Memphis_cutoff_ERG <- read_xlsx(path = file.path(output_dir_Memphis, "results_DEA.xlsx"), 
                                    sheet = 1)


# merge the tables
merged_table_Memphis <- merge(x = as.data.frame(MLLr_Memphis_cutoff_ERG),
                      y = fit_GFP_TPOVsnoTPO,
                      by.x = "ensembl_gene_id",
                      by.y = "ensembl_gene_id.y", 
                      suffixes = c("_Memphis", "_mouse"))


# create the color for the significativity 
merged_table_Memphis$color_sign <- "n.s"

# update the color based on the significativity
merged_table_Memphis$color_sign[merged_table_Memphis$FDR_mouse <= 0.05] <- "sign. mouse"
merged_table_Memphis$color_sign[merged_table_Memphis$FDR_Memphis <= 0.05] <- "sign. Memphis"
merged_table_Memphis$color_sign[(merged_table_Memphis$FDR_mouse <= 0.05) & (merged_table_Memphis$FDR_Memphis <= 0.05)] <- "sign. mouse and Memphis"

# add list of relevant genes
list_relevant_genes <- c("IL12RB2", "INPP4B", "MECOM", "ADGRG6", "MMP28", "SEMA4F",
                         "SH3BP5", "AIF1", "CDC14A", "FARP1", "GRAP2", "KRI1", 
                         "MLLT3", "MOCOS", "OAS3", "PBX1", "PDLIM2", "PDLIM5", 
                         "SERPINE2", "SLC6A12", "STYK1", "TGM2", "SLC7A2", "ZNF704",
                         "AKAP12", "CACNA2D2", "CALR", "CPSF4", "CST7", "FKBP11",
                         "FNDC3B", "MS4A3", "MST1", "P2RY2", "P4HB", "PHKG1", "PMP22",
                         "PRAG1", "TDRD9", "TRIM9", "UXS1")


# calculate the coefficient correlation
fit <- lm(merged_table_Memphis$logFC_Memphis ~ merged_table_Memphis$logFC_mouse, data = merged_table_Memphis)
coeff_corr <- round(cor(merged_table_Memphis$logFC_Memphis, merged_table_Memphis$logFC_mouse, method = "pearson"), 3)


# draw the scatterplot
scatterplot <- ggplot(merged_table_Memphis, aes(x = logFC_Memphis, y = logFC_mouse, color = color_sign)) +
  geom_point(alpha = 0.75) + theme_bw() + xlab("log2FC Memphis") + ylab("log2FC mouse") +
  xlim(-max(abs(merged_table_Memphis$logFC_Memphis)), max(abs(merged_table_Memphis$logFC_Memphis))) +
  ylim(-max(abs(merged_table_Memphis$logFC_mouse)), max(abs(merged_table_Memphis$logFC_mouse))) +
  scale_color_manual(values = c("n.s" = "black", "sign. mouse" = "blue", 
                                "sign. Memphis" = "orange", "sign. mouse and Memphis" = "red")) +
  geom_text_repel(data = subset(merged_table_Memphis[merged_table_Memphis$external_gene_name %in% list_relevant_genes,]), 
                  color = "black", aes(label = external_gene_name), size = 5) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between Memphis and GFP+ TPO vs PBS; corr coef (Pearson) = ", coeff_corr))


# save the scatterplot in files
saveFigures("scatterplot_GFP+_TPOvsPBS_relevantGenes", ggplot = scatterplot, dirPlot = output_dir_Memphis, A4 = T)


# save the table in an excel file
write_xlsx(merged_table_Memphis, path = "DEA/Human_patients/St-Jude/allAML/St-Jude_ERG_selection_kmeans_k5/mergedTable_Memphis_GFP+_TPOvsPBS.xlsx")


# create the scatterplot for common DEGs compare to MECOM ---------------------


# create the dataframe for the scatterplot figure
dt_Memphis_ERG_MECOM <- data.frame("MECOM" = rpkm_dge_Memphis["MECOM", ],
                                   "IL12RB2" = rpkm_dge_Memphis["IL12RB2", ],
                                   "INPP4B" = rpkm_dge_Memphis["INPP4B", ],
                                   "Sample" = colnames(rpkm_dge_Memphis),
                                   "Color" = "none")

# update the patients information
dt_Memphis_ERG_MECOM$Color[dt_Memphis_ERG_MECOM$Sample %in% patients_information_Memphis$Patient[grep(pattern = "KMT2Ar", patients_information_Memphis$`Molecular category`)]] = "MLLr"
dt_Memphis_ERG_MECOM$Color[dt_Memphis_ERG_MECOM$Sample %in% patients_information_Memphis$Patient[grep(pattern = "KMT2A::MLLT3", patients_information_Memphis$`Defining alteration`)]] = "MLL-AF9"
dt_Memphis_ERG_MECOM$Color[dt_Memphis_ERG_MECOM_log2rpkm$Sample %in% patients_information_Memphis$Patient[grep(pattern = "MECOM", patients_information_Memphis$`Molecular category`)]] = "EVI1r"

# calculate the coefficient correlation
dt_Memphis_ERG_MECOM_tmp <- dt_Memphis_ERG_MECOM[dt_Memphis_ERG_MECOM$MECOM != 0, ]
dt_Memphis_ERG_MECOM_tmp <- dt_Memphis_ERG_MECOM_tmp[dt_Memphis_ERG_MECOM_tmp$IL12RB2 != 0, ]
fit <- lm(log2(dt_Memphis_ERG_MECOM_tmp$IL12RB2) ~ log2(dt_Memphis_ERG_MECOM_tmp$MECOM), data = dt_Memphis_ERG_MECOM_tmp)
coeff_corr <- round(cor(log2(dt_Memphis_ERG_MECOM_tmp$MECOM), log2(dt_Memphis_ERG_MECOM_tmp$IL12RB2), method = "pearson"), 3)

# create the scatterplot
scatterplot <- getScatterplot(data = dt_Memphis_ERG_MECOM_tmp, x = "MECOM", y = "IL12RB2", log2 = T)
scatterplot <- scatterplot + geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between MECOM and IL12RB2; corr coef (Pearson) = ", coeff_corr)) +
  scale_color_manual(values = c("none" = "purple", "MLLr" = "cyan3", 
                                "MLL-AF9" = "forestgreen", "EVI1r" = "red")) + geom_point(size=3)

# create the scatterplot for the Erg vs MECOM expression for all patients
saveFigures("scatterplot_IL12RB2vsMECOM", scatterplot, 
            dirPlot = "DEA/Human_patients/St-Jude/allAML/", A4 = T)

# calculate the coefficient correlation
dt_Memphis_ERG_MECOM_tmp <- dt_Memphis_ERG_MECOM[dt_Memphis_ERG_MECOM$MECOM != 0, ]
dt_Memphis_ERG_MECOM_tmp <- dt_Memphis_ERG_MECOM_tmp[dt_Memphis_ERG_MECOM_tmp$INPP4B != 0, ]
fit <- lm(log2(dt_Memphis_ERG_MECOM_tmp$INPP4B) ~ log2(dt_Memphis_ERG_MECOM_tmp$MECOM), data = dt_Memphis_ERG_MECOM_tmp)
coeff_corr <- round(cor(log2(dt_Memphis_ERG_MECOM_tmp$MECOM), log2(dt_Memphis_ERG_MECOM_tmp$INPP4B), method = "pearson"), 3)

# create the scatterplot
scatterplot <- getScatterplot(data = dt_Memphis_ERG_MECOM_tmp, x = "MECOM", y = "INPP4B", log2 = T)
scatterplot <- scatterplot + geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
  ggtitle(paste0("Comparison between MECOM and INPP4B; corr coef (Pearson) = ", coeff_corr)) +
  scale_color_manual(values = c("none" = "purple", "MLLr" = "cyan3", 
                                "MLL-AF9" = "forestgreen", "EVI1r" = "red")) + geom_point(size=3)

# create the scatterplot for the Erg vs MECOM expression for all patients
saveFigures("scatterplot_INPP4BvsMECOM", scatterplot, 
            dirPlot = "DEA/Human_patients/St-Jude/allAML/", A4 = T)