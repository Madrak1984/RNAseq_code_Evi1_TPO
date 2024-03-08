###### command line to perform the differential expression analysis for the TPO-Evi1 mRNA
# done by Jonathan Seguin, group of Prof. Schwaller, DBM, UKBB, Basel, Switzerland
# email: jonathan.seguin@unibas.ch, seguin.jonathan@gmail.com
# created in Tue Mar 15 10:44:34 2022


# prepare the scripts -----------------------------------------------------


# load the libraries
library(edgeR)
library(DESeq2)
require("ggrepel") # http://www.sthda.com/french/wiki/ggplot2-textes-ajouter-du-texte-a-un-graphique-logiciel-r-et-visualisation-de-donnees
library(xlsx)
library(biomaRt)
library(stringr)
library(writexl)
library(plotly)
library(fgsea)
library(tidyverse)
library(pathview)
library(MetaboSignal)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)
library(readxl)
library(msigdbr)
library(factoextra)
library(rgl)
library(gridExtra)
library(limma)
library(pheatmap)

# Palette with 20 colors (+black and white) that do not conflict so much, adapted from https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
myPalette2 <-  c('#e6194b', '#4363d8', '#3cb44b', '#984EA3', '#f58231', '#ffe119', '#F781BF', '#808080', '#98BFDB', '#bcf60c', '#008080', '#e6beff', '#E5C494', '#000075', '#CD00CD', '#aaffc3', '#808000', '#9a6324', '#fffac8', '#800000', '#000000', '#ffffff')




# load the functions ------------------------------------------------------
source("Src/DE_functions.R")



# create the output directories -------------------------------------------

# create the plot directory
plotDir <- "DEA/Plots"
if(!dir.exists(plotDir)) dir.create(plotDir, recursive = T)

# create the RDS directory
RDSDir <-  "DEA/Rds"
if(!dir.exists(RDSDir)) dir.create(RDSDir)


# create the pathway
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

# convert in list
THPO_pathway <- list("TPO_pathway" = THPO_pathway)

# get the list of genes in mouse
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl")
genes_mouse <-  getBM(attributes = c("ensembl_gene_id","entrezgene_id", "external_gene_name"), mart = mart)



# read the rds files ------------------------------------------------------

# read the rds file from the snakemake pipeline
dge <- readRDS("Snakemake/ensdb_102_dge_list.rds")
se <-  readRDS("Snakemake/ensdb_102_summarized_experiment.rds")

# Data exploration --------------------------------------------------------

# remove the genes without read 
dge = dge[which(rowSums(dge$counts) > 10),]

# replace the "+" and "-" symbol in the SampleGroup to prevent error in the comparison analysis
dge$samples$SampleGroup = gsub(pattern = "GFP-", replacement = "GFPneg", 
                               dge$samples$SampleGroup)
dge$samples$SampleGroup = gsub(pattern = "GFP\\+", replacement = "GFP", 
                               dge$samples$SampleGroup)


# calculate the logCpm
log2cpm <- cpm(dge, log = T, prior.count = 8)
png(file.path(plotDir, "densities.png"))
plotDensities(log2cpm, legend = "topright")
abline(v=1)
dev.off()

# save the PCAs
savePCA(dge_tmp = dge, plotDir = plotDir)




# PCA without the 35088 and 624-2 samples ---------------------------------

# select the identifiers after removing the 35088 and 624-2 samples
id_samples <-  intersect(grep(pattern = "35088", dge$samples$ExternalSampleName,
                            invert = T), grep(pattern = "624-2", dge$samples$ExternalSampleName, 
                                               invert = T))

# remove the samples
dge_samples <- dge[, id_samples]

# update the dge samples
dge_samples$samples$SampleGroup <- gsub(pattern = "noTPO", replacement = "PBS",
                                        dge_samples$samples$SampleGroup)

# save the PCAs
#savePCA(dge_tmp = dge_samples, extensionName = "_wosamples",  plotDir = plotDir)
savePCA(dge_tmp = dge_samples, extensionName = "_wosamples_paper",  plotDir = plotDir, 
        theme_bw = T, text = F)




# Perform the differential expression analysis after removing the 35088 and 624-2 samples----------------------------


# create the RDS dir
#RDSDir_woS <- "DEA/Rds/RDS_wo624-2"
RDSDir_woS <- "DEA/Rds/RDS_woD/"
if(!dir.exists(RDSDir_woS)) dir.create(RDSDir_woS, recursive = T )

# calculate the logCpm
log2cpm <- cpm(dge_samples, log = T, prior.count = 8)

# normalization and estimate dispersion
df.dge_samples <- calcNormFactors(dge_samples)
moma <- model.matrix(~ 0 + SampleGroup, data=dge_samples$samples) # ~ 0 + SampleGroup
colnames(moma) <- sub("SampleGroup", "", colnames(moma))

# prepare the contrasts
contrasts.matrix <- makeContrasts(
  
  # bulk, comparison between TPO vs no_TPO
  bulk.TPO = TPO_bulk - noTPO_bulk,
  
  # GFP+ samples, comparison between TPO vs no_TPO
  GFP.TPO = TPO_GFP - noTPO_GFP,
  
  # GFP- samples, comparison between TPO vs no_TPO
  GFPneg.TPO = TPO_GFPneg - noTPO_GFPneg,
  
  # TPO samples, comparison GFP- vs bulk
  TPO.GFPnegvsbulk = TPO_GFPneg - TPO_bulk,
  
  # TPO samples, comparison GFP+ vs bulk
  TPO.GFPvsbulk = TPO_GFP - TPO_bulk,
  
  # TPO samples, comparison GFP+ vs bulk
  TPO.GFPvsGFPneg = TPO_GFP - TPO_GFPneg,
  
  # TPO samples, comparison GFP- vs bulk
  noTPO.GFPnegvsbulk = noTPO_GFPneg - noTPO_bulk,
  
  # TPO samples, comparison GFP+ vs bulk
  noTPO.GFPvsbulk = noTPO_GFP - noTPO_bulk,
  
  # TPO samples, comparison GFP+ vs bulk
  noTPO.GFPvsGFPneg = noTPO_GFP - noTPO_GFPneg,
  
  # comparison between GFP+ vs GFP-
  GFPvsGFPneg = (TPO_GFP + noTPO_GFP)/2 - (TPO_GFPneg + noTPO_GFPneg)/2,
  
  # comparison between TPO+ VS TPO- with all GFP
  allGFP.TPOvsTPOneg = (TPO_GFP + TPO_GFPneg)/2 - (noTPO_GFP + noTPO_GFPneg)/2,
  
  # special comparison asked by HE, must be similar to late vs early in Vaïa's dataset
  TPO.GFPvsnoTPO.GFPneg = TPO_GFP - noTPO_GFPneg,
  
  # special comparison asked by HE, must be similar to late vs early in Vaïa's dataset
  TPO.GFPvsnoTPO.bulk = TPO_GFP - noTPO_bulk,
  
  # save the levels
  levels = moma
)

## filtering of lowly expressed genes                                                                                       
## based on the figures plotDensity above                                                                                               
dge_samples <- dge_samples[rowSums(log2cpm > 1) >= 2,]                                                                                     
dge_samples$samples$lib.size <- colSums(dge_samples$counts)                                                                                 
dge_samples <- calcNormFactors(dge_samples)                                                                                                
dge_samples <- estimateDisp(dge_samples, design=moma)   
saveRDS(dge_samples, file= file.path(RDSDir_woS,
                                     "EnsemblGenes.dge_samplesList.filtered.rds"))  



############## Perform the DEA

# calculate the fit
fit <- glmQLFit(dge_samples, design = moma)

# create the pdf 
pdf(file.path(plotDir, "DE_analysis_glmQLFit_DEA_wo624-2.pdf"), width = 6, height = 6)

# prepare the variable
decideTestsDGE_glmQLFit = NULL
fits <- list()
for (i in 1:ncol(contrasts.matrix)) {
  
  # select the contrasts
  contr.name <- colnames(contrasts.matrix)[i]
  message("   Working on ... ", contr.name)
  
  ## Subset samples
  dge_samples.temp <- dge_samples[, dge_samples$samples$group %in% names(which(contrasts.matrix[,i] != 0))]
  
  ## Redo TMM normalization
  dge_samples.temp <- calcNormFactors(dge_samples.temp) 
  ## After TMM normalization
  log2cpm.temp <- edgeR::cpm(dge_samples.temp, log=TRUE, normalized.lib.sizes=TRUE)
  
  # perform the differential expression test
  qlf <- glmQLFTest(fit, contrast = contrasts.matrix[, contr.name])
  
  # plot the MD
  plotMD(qlf)
  
  
  # DE genes
  print(table(de <- decideTestsDGE(qlf)))
  decideT =  table(de <- decideTestsDGE(qlf))
  
  if (length(decideT)!= 3) {
    if ('-1'%in%names(decideT)==F) {
      decideT['-1'] = 0
      
    } 
    
    if ('1'%in%names(decideT)==F) {
      decideT['1'] = 0
    } 
    
    decideT = decideT[c('-1','0','1')]
    
  }
  
  decideTestsDGE_glmQLFit = rbind(decideTestsDGE_glmQLFit, decideT)
  
  
  ## MA plot
  detags <- rownames(dge_samples)[as.logical(de)]
  plotSmear(qlf, de.tags=detags, main=contr.name)
  abline(h=c(-1, 1), col=myPalette2[2])
  
  ## add results to list
  fits[[contr.name]] <- qlf
  rm(qlf)
  
  ## add info needed by Shiny app
  fits[[contr.name]]$contrasts <- contrasts.matrix[names(which(contrasts.matrix[,i] != 0)), i, drop=FALSE]
  
}
dev.off()



# save the fits
# Save for shiny-app
saveRDS(fits, file=file.path(RDSDir_woS, "edge_samplesR_fits_results_DEA.rds"))


# Analyze the DEA results -------------------------------------------------

# load the fits 
fits <- readRDS(file = "DEA/Rds/RDS_woD/edge_samplesR_fits_results_DEA.rds")

# list the names of fits
names_fits <- names(fits)

# create the dataframe to save the count of DEGs
data_count <- data.frame("comparison" = names_fits, "down" = 0, "NS" = 0, "up" = 0)

# loop to analyze the fits list
for(name in names_fits){
  print(name)
  
  # select the dge_samples_list
  dge_samples_name <- fits[[name]]
  
  # select the toptags
  dge_samples_DEG <- as.data.frame(topTags(dge_samples_name, n = nrow(dge_samples_name), p.value = 0.05))
  
  # count the DEG down or up
  if(nrow(dge_samples_DEG[dge_samples_DEG$logFC <= 0,]) > 0){
    data_count$down[which(data_count$comparison == name)] = nrow(dge_samples_DEG[dge_samples_DEG$logFC <= 0,])  
  }
  if(nrow(dge_samples_DEG[dge_samples_DEG$logFC >= 0,]) > 0){
    data_count$up[which(data_count$comparison == name)] = nrow(dge_samples_DEG[dge_samples_DEG$logFC >= 0,])  
  }
  print(paste0("count for the dge_samples_DEG ", nrow(dge_samples_name)))
  data_count$NS[which(data_count$comparison == name)] = nrow(dge_samples_name) - 
    (data_count$down[which(data_count$comparison == name)] + 
       data_count$up[which(data_count$comparison == name)])
}


# create the output directory
#outDir <- "DEA/Result/Result_wo624-2"
outDir <- "DEA/Result/Result_woSamples/"

if(!dir.exists(outDir)) dir.create(outDir, recursive = T)

# save the excel files
write.xlsx2(data_count, file = file.path(outDir, "count_DEGs.xlsx"), col.names = T, row.names = F)

# find genes information
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl")

# list the names of fits
names_fits <- names(fits)


# save the DEG in excel files
for(name in names_fits){
  
  # print the message
  print(paste0(date(), ": output for ", name))
  
  # create the tmp_dir
  tmp_dir <- file.path(outDir, name)
  if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = T)
  
  # select the dge_samples_list
  dge_samples_name <- fits[[name]]
  
  # select the toptags
  dge_samples_DEG <- as.data.frame(topTags(dge_samples_name, n = nrow(dge_samples_name)))
  
  # count the DEG down or up
  write.xlsx2(dge_samples_DEG, file = file.path(tmp_dir, paste0("DEGs_", name, ".xlsx")), col.names = T, row.names = F)
  gc()
  
  
  # select the dge_samples_list
  volcano <- getVolcanoplot(dgeobject = fits[[name]], nbGenes = 10)
  saveFigures(fileName = paste0("volcanoplot_", name), dirPlot = tmp_dir, ggplot = volcano)
  
  
  # select the significant toptags
  dge_samples_DEG <- as.data.frame(topTags(dge_samples_name, n = nrow(dge_samples_name), p.value = 0.05))
  
  # select all the analyzed genes
  universe <- rownames(topTags(dge_samples_name, n = nrow(dge_samples_name)))
  
  # perform the enrichment analysis
  print(paste0(date(), ": start enrichment analysis for ", name))
  
  if(!is.null(dge_samples_DEG)){
    if(nrow(dge_samples_DEG) > 0){
      EnrichmentTerm_Analysis(common = rownames(dge_samples_DEG), outputDir = tmp_dir, mart = mart, universe = universe)
    }
  } 
  print(paste0(date(), ": end enrichment analysis for ", name))
  
}

# perform term Analysis for common and specific genes between GFP+ and bulk ---------------------------------

## all DEGs ---------------------------------

# find common genes 
DEG_GFP <- as.data.frame(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO),
                                 p.value = 0.05))
DEG_bulk <- as.data.frame(topTags(fits$bulk.TPO, n = nrow(fits$bulk.TPO),
                                  p.value = 0.05))




# perform the GSEA analysis -----------------------------------------------------------

# use the camera tools as found in the Maria's script
## camera (Gene set enrichment analysis)
## gene annotation is downloaded from http://bioinf.wehi.edu.au/software/MSigDB/, preprocessed as RData lists
## original data from MsigDb http://software.broadinstitute.org/gsea/msigdb/index.jsp


# init the list for camera results
list_cameras <- NULL

# save the different collection from Msig database for camera
db_cameras = c("H", "c2", "c3", "c4", "c5", "c6", "c7")

# start the loop
for(db in db_cameras){
  
  # create the filename
  filename = paste0("mouse_", db, "_v5p2.rdata")
  
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

# update the name
names(list_cameras) = paste0("camera_dge_", db_cameras)

## test with fGSEA -----------------------------------------------------------

### perform GSEA -----------------------------------------------------------

# save the list of databases
dbs <- paste0("Mm.", c("H", paste0("c", 2:7)))

# select the GFP+ TPO vs noTPO results
GSEA_hallmark(fit = fits$GFP.TPO, list_dbs = dbs,
              path = "DEA/Result/Result_woSamples/GFP.TPO/", filename = "GSEA_Msigdatabase_GFP_TPOvsnoTPO.xlsx")



### Redo the analysis for the bulk TPO vs noTPO
# select the GFP+ TPO vs noTPO results
GSEA_hallmark(fit = fits$bulk.TPO, list_dbs = dbs,
              path = "DEA/Result/Result_woSamples/bulk.TPO/", filename = "GSEA_Msigdatabase_bulk_TPOvsnoTPO.xlsx")


### Redo the analysis for the TPO GFP+ vs GFP-
# select the TPO GFP+ vs GFP- results
GSEA_hallmark(fit = fits$TPO.GFPvsGFPneg, list_dbs = dbs,
              path = "DEA/Result/Result_woSamples/TPO.GFPvsGFPneg/", filename = "GSEA_Msigdatabase_TPO_GFPposVsGFPneg.xlsx")


#### perform GSEA with selected marks for custom pathways (source: https://doi.org/10.1182/bloodadvances.2018025866) -----------------------------------------------------------

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


# create the custom pathways
custom_pathways <- list("Jak-Stat pathway" = jakStat$ensembl_gene_id,
                        "MapK pathway" = MAPK$ensembl_gene_id,
                        "NF-Kappa-B pathway" = NF_kappaB$ensembl_gene_id,
                        "proteoglycans pathway" = proteoglycans$ensembl_gene_id
)


# perform the GSEA
fgseaRes_custom_pathways <- fgsea(pathways = custom_pathways, stats = geneList_GSEA)

# plot the enrichment
for(enrich in names(custom_pathways)){
  
  # define the filename
  filename <- gsub(pattern = " ", replacement = "_", enrich)
  filename = paste0("GSEA_", filename, "_GFP_TPOvsnoTPO")
  
  # save the enrichment in the pdf
  saveFigures(fileName = filename, ggplot = plotEnrichment(custom_pathways[[enrich]], geneList_GSEA) + 
                labs(title = enrich), dirPlot = "DEA/Result/Result_woSamples/GFP.TPO/", A4 = T)
  
}

# write the excel file
write_xlsx(fgseaRes_custom_pathways, path = "DEA/Result/Result_woSamples/GFP.TPO/GSEA_custom_pathway_GFP_TPOvsnoTPO.xlsx")




#### perform GSEA with selected marks for custom pathways for Bulk analysis (source: https://doi.org/10.1182/bloodadvances.2018025866) -----------------------------------------------------------

# create the list of genes for GSEA
id_sort = order(fits$bulk.TPO$table$logFC, decreasing = T)
table_GSEA_bulk = fits$bulk.TPO$table[id_sort,]
geneList_GSEA_bulk = table_GSEA_bulk$logFC
names(geneList_GSEA_bulk) = rownames(table_GSEA_bulk)

# perform the GSEA
fgseaRes_custom_pathways <- fgsea(pathways = custom_pathways, stats = geneList_GSEA_bulk)

# plot the enrichment
for(enrich in names(custom_pathways)){
  
  # define the filename
  filename <- gsub(pattern = " ", replacement = "_", enrich)
  filename = paste0("GSEA_", filename, "_bulk_TPOvsnoTPO")
  
  # save the enrichment in the pdf
  saveFigures(fileName = filename, ggplot = plotEnrichment(custom_pathways[[enrich]], geneList_GSEA_bulk) + 
                labs(title = enrich), dirPlot = "DEA/Result/Result_woSamples/bulk.TPO/", A4 = T)
  
}

# write the excel file
write_xlsx(fgseaRes_custom_pathways, path = "DEA/Result/Result_woSamples/bulk.TPO/GSEA_custom_pathway_bulk_TPOvsnoTPO.xlsx")



### deeper analysis of selected term -----------------------------------------------------------

#### GFP+ TPO vs noTPO -----------------------------------------------------------

# select gene list for GFP+ TPO vs noTPO
geneList <- getGeneListforGSEA(fits$GFP.TPO)


# upload the GSEA result with Hugues-Etienne selection
GSEA_res_C2 <- read.xlsx2(file = "DEA/Result/Result_woSamples/GFP.TPO/GSEA_Msigdatabase_GFP_TPOvsnoTPO_HE.xlsx", 
                          header = T, sheetName = "Mm.c2", colClasses = c("character", rep("numeric", 6), "character"))
GSEA_res_H <- read.xlsx2(file = "DEA/Result/Result_woSamples/GFP.TPO/GSEA_Msigdatabase_GFP_TPOvsnoTPO_HE.xlsx",
                         header = T, sheetName = "Mm.H")


# select the HOXA9 pathway
GSEA_res_C2_HOXA9 <- GSEA_res_C2[grep(pattern = "HOXA9", GSEA_res_C2$pathway), ]
GSEA_res_C2_HOXA9 <- rbind(GSEA_res_C2_HOXA9, GSEA_res_C2[grep(pattern = "^KUMAR", GSEA_res_C2$pathway), ])
GSEA_res_C2_HOXA9 <- rbind(GSEA_res_C2_HOXA9, GSEA_res_C2[grep(pattern = "^VERHAAK_AML", GSEA_res_C2$pathway), ])

# count the number of elements
dt_count <- data.frame("Pathway" = names(Mm.c2), "Counts" = unlist(lapply(Mm.c2, length)))

# merge the count
GSEA_res_C2_HOXA9 = merge(x = GSEA_res_C2_HOXA9, y = dt_count, by.x = "pathway",
                          by.y = "Pathway")

# calculate the GeneRatio
GSEA_res_C2_HOXA9$GeneRatio = as.numeric(GSEA_res_C2_HOXA9$size) / GSEA_res_C2_HOXA9$Counts

# reorder according to the GeneRatio
id_GR <-  order(GSEA_res_C2_HOXA9$GeneRatio)

# reorder the table
GSEA_res_C2_HOXA9 = GSEA_res_C2_HOXA9[id_GR,]

# convert vector as factor
GSEA_res_C2_HOXA9$pathway = as.factor(GSEA_res_C2_HOXA9$pathway)

# change the column names
colnames(GSEA_res_C2_HOXA9)[c(3,7)] = c("p.adjust", "Count")

# create the dotplot
ggplot_GSEA <- GSEA_res_C2_HOXA9 %>%
  mutate(pathway = fct_reorder(pathway, GeneRatio)) %>%
  ggplot(aes(x = GeneRatio, y = pathway, colour = p.adjust)) +
  geom_point(aes(size = Count)) +
  scale_colour_gradient(low = "red", high = "blue") +
  theme_bw() + theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black"))

# save the figure in the files
saveFigures(fileName = "GSEA_paper_HOXA9", ggplot = ggplot_GSEA, 
            dirPlot = "DEA/Result/Result_woSamples/GFP.TPO/")


# create figures for paper ------------------------------------------------
 
## perform GSEA based on the TPO pathway for GFP+ ------------------------------------------------
# source (paper: 10.1007/s12079-018-0480-4) 
 
 
# set the table for GFP.TPO
table_GFP_TPO = as.data.frame(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO))) 

# add the EnsemblID
table_GFP_TPO$EnsemblID = rownames(table_GFP_TPO)

# merge with human orthologous gene
table_GFP_TPO = merge(x=table_GFP_TPO, y=ortho, by.x = "EnsemblID", 
                      by.y = "mmusculus_homolog_ensembl_gene")
 
 # reorder the genes 
id_sort = order(table_GFP_TPO$logFC, decreasing = T)
table_GSEA = table_GFP_TPO[id_sort,]
geneList_GSEA = table_GSEA$logFC
names(geneList_GSEA) = table_GSEA$external_gene_name

# perform the GSEA analysis
fgseaRes <- fgsea(pathways = THPO_pathway, stats = geneList_GSEA)

# save the result
write_xlsx(fgseaRes, path = file.path("DEA/Result/Result_woSamples/GFP.TPO/",
                                      "fgseaResult_TPO_pathways_GFP+_TPOvsnoTPO.xlsx"))

# draw the GSEA results
saveFigures(fileName = "fgseaResult_TPO_pathways_GFP+_TPOvsnoTPO", 
            ggplot = plotEnrichment(THPO_pathway$TPO_pathway, geneList_GSEA) + 
              labs(title = "TPO pathway, GFP+ TPO vs PBS"), 
            dirPlot = "DEA/Result/Result_woSamples/GFP.TPO/", A4 = T)

# select the gene names
genes <- unique(table_GFP_TPO$EnsemblID[table_GFP_TPO$external_gene_name %in% names(geneList_GSEA)])
intersect(rownames(mat_norm_count), genes)

# draw the heatmap
drawHeatmap(dge = dge, genes = table_GFP_TPO$EnsemblID[table_GFP_TPO$external_gene_name %in% THPO_pathway$TPO_pathway], 
            mat_norm_count = mat_norm_count, dirPlot = "DEA/Result/Result_woSamples/GFP.TPO/", 
            filename = paste0("heatmap_TPO_pathway"))





## perform GSEA based on the TPO pathway for bulk ------------------------------------------------


# set the table for bulk.TPO
table_bulk_TPO = as.data.frame(topTags(fits$bulk.TPO, n = nrow(fits$bulk.TPO))) 

# add the EnsemblID
table_bulk_TPO$EnsemblID = rownames(table_bulk_TPO)

# merge with human orthologous gene
table_bulk_TPO = merge(x=table_bulk_TPO, y=ortho, by.x = "EnsemblID", 
                      by.y = "mmusculus_homolog_ensembl_gene")

# reorder the genes 
id_sort = order(table_bulk_TPO$logFC, decreasing = T)
table_GSEA = table_bulk_TPO[id_sort,]
geneList_GSEA = table_GSEA$logFC
names(geneList_GSEA) = table_GSEA$external_gene_name

# perform the GSEA analysis
fgseaRes <- fgsea(pathways = THPO_pathway, stats = geneList_GSEA)

# save the result
write_xlsx(fgseaRes, path = file.path("DEA/Result/Result_woSamples/bulk.TPO/",
                                      "fgseaResult_TPO_pathways_bulk_TPOvsnoTPO.xlsx"))

# draw the GSEA results
saveFigures(fileName = "fgseaResult_TPO_pathways_bulk_TPOvsnoTPO", 
            ggplot = plotEnrichment(THPO_pathway$TPO_pathway, geneList_GSEA) + 
              labs(title = "TPO pathway, bulk TPO vs PBS"), 
            dirPlot = "DEA/Result/Result_woSamples/bulk.TPO/", A4 = T)



 
## perform statistics on pathways ------------------------------------------------


# select the HOXA9 pathway
GSEA_res_C2_HOXA9 <- GSEA_res_C2[grep(pattern = "HOXA9", GSEA_res_C2$pathway), ]

# count the number of elements
dt_count <- data.frame("Pathway" = names(Mm.c2), "Counts" = unlist(lapply(Mm.c2, length)))

# merge the count
GSEA_res_C2_HOXA9 = merge(x = GSEA_res_C2_HOXA9, y = dt_count, by.x = "pathway",
                          by.y = "Pathway")

# calculate the GeneRatio
GSEA_res_C2_HOXA9$GeneRatio = as.numeric(GSEA_res_C2_HOXA9$size) / GSEA_res_C2_HOXA9$Counts

# reorder according to the GeneRatio
id_GR <-  order(GSEA_res_C2_HOXA9$GeneRatio)

# reorder the table
GSEA_res_C2_HOXA9 = GSEA_res_C2_HOXA9[id_GR,]

# convert vector as factor
GSEA_res_C2_HOXA9$pathway = as.factor(GSEA_res_C2_HOXA9$pathway)

# change the column names
colnames(GSEA_res_C2_HOXA9)[c(3,7)] = c("p.adjust", "Count")

# create the dotplot
ggplot_GSEA <- GSEA_res_C2_HOXA9 %>%
  mutate(pathway = fct_reorder(pathway, GeneRatio)) %>%
  ggplot(aes(x = GeneRatio, y = pathway, colour = p.adjust)) +
  geom_point(aes(size = Count)) +
  scale_colour_gradient(low = "red", high = "blue") +
  theme_bw() + theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(colour = "black"), 
        axis.text.y = element_text(colour = "black"))

# save the figure in the files
saveFigures(fileName = "GSEA_paper_HOXA9", ggplot = ggplot_GSEA, 
            dirPlot = "DEA/Result/Result_woSamples/GFP.TPO/")



# compare the 3 TPO vs noTPO comparison ------------------------------------------------

# select all DEGs in the GFP+
DEG_GFPp <- as.data.frame(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO), 
                                  p.value = 0.05))

# select all DEGs in the GFP-
DEG_GFPn <- as.data.frame(topTags(fits$GFPneg.TPO, n = nrow(fits$GFPneg.TPO), 
                                  p.value = 0.05))

# select all DEGs in the bulk
DEG_bulk <- as.data.frame(topTags(fits$bulk.TPO, n = nrow(fits$bulk.TPO),
                                  p.value = 0.05))

# find commons DEGs
common_DEGs_TPOcomp <- intersect(intersect(rownames(DEG_bulk), rownames(DEG_GFPp)),
                         rownames(DEG_GFPn))

# find commons up DEGs
common_DEGs_TPOcomp_up <- intersect(intersect(rownames(DEG_bulk[DEG_bulk$logFC > 0,]), 
                                 rownames(DEG_GFPp[DEG_GFPp$logFC > 0,])),
                                 rownames(DEG_GFPn[DEG_GFPn$logFC > 0,]))

# find commons down DEGs
common_DEGs_TPOcomp_down <- intersect(intersect(rownames(DEG_bulk[DEG_bulk$logFC < 0,]), 
                                              rownames(DEG_GFPp[DEG_GFPp$logFC < 0,])),
                                    rownames(DEG_GFPn[DEG_GFPn$logFC < 0,]))

# check if all DEGs are separated according the same orientation
length(common_DEGs_TPOcomp) == length(common_DEGs_TPOcomp_down) + length(common_DEGs_TPOcomp_up)

# create the corresponding Venn Diagram
myCol <- brewer.pal(3, "Pastel2")
draw_VennDiagram(x = list(rownames(DEG_bulk), rownames(DEG_GFPp), rownames(DEG_GFPn)), 
                 category = c("Bulk" , "GFP+" , "GFP-"), 
                 filename = 'DEA/Result/Result_woSamples/allGFP.TPOvsTPOneg/vennDiagram.pdf',
                 color = myCol)
draw_VennDiagram(x = list(rownames(DEG_bulk[DEG_bulk$logFC > 0,]), rownames(DEG_GFPp[DEG_GFPp$logFC > 0,]), 
                          rownames(DEG_GFPn[DEG_GFPn$logFC > 0,])), 
                 category = c("Bulk" , "GFP+" , "GFP-"), 
                 filename = 'DEA/Result/Result_woSamples/allGFP.TPOvsTPOneg/vennDiagram_up.pdf',
                 color = myCol)
draw_VennDiagram(x = list(rownames(DEG_bulk[DEG_bulk$logFC < 0,]), rownames(DEG_GFPp[DEG_GFPp$logFC < 0,]), 
                          rownames(DEG_GFPn[DEG_GFPn$logFC < 0,])), 
                 category = c("Bulk" , "GFP+" , "GFP-"), 
                 filename = 'DEA/Result/Result_woSamples/allGFP.TPOvsTPOneg/vennDiagram_down.pdf',
                 color = myCol)






# creation zscore heatmap for paper ---------------------------------------

# load the fits
fits <- readRDS(file = "DEA/Rds/RDS_woD/edge_samplesR_fits_results_DEA.rds")

# load the dge
dge <- readRDS(file = "DEA/Rds/RDS_woD/EnsemblGenes.dge_samplesList.filtered.rds")

# get the log cpm of dge
logcpm <- cpm(dge, log = T)


## create the heatmap for GFP+ -----------------------------------------------------

# select the EdgeR result
GFP_TPO <- as.data.frame(topTags(object = fits$GFP.TPO, n = nrow(fits$GFP.TPO)))

# select the top100 DEG (based on the FDR)
top100 <- GFP_TPO %>% arrange(FDR) %>% 
  filter(GENEBIOTYPE == "protein_coding" & !grepl("RIKEN", DESCRIPTION)) %>% 
  slice_head(n = 100) %>% select(SYMBOL)

# select the samples
sel_samples <- dge$samples[grep(pattern = "GFP$", dge$samples$SampleGroup),]

# select the genes
logcpm_sel <- logcpm[rownames(top100), rownames(sel_samples)]

# rename the genes
rownames(logcpm_sel) = top100$SYMBOL

# create dataframe for the sample
sample_col <- sel_samples %>% select(SampleGroup)

# draw the zscore heatmap 
heatmap <- pheatmap(logcpm_sel,
                    color = colorRampPalette(c("blue2","white","red"))(100),
                    scale = "row",
                    clustering_method = "ward.D2", 
                    annotation_col = sample_col,
                    fontsize_row = 4)

# save the heatmap in the files
saveFigures(fileName = "heatmap_paper_GFP+_TPOvsPBS", ggplot = heatmap, dirPlot = "DEA/Result/Result_woSamples/GFP.TPO/", A4 = T)

# keep Genes order for the bulk heatmap
top100_order = heatmap$tree_row$labels[heatmap$tree_row$order]




## create the heatmap for bulk+ -----------------------------------------------------

# select the EdgeR result
bulk_TPO <- as.data.frame(topTags(object = fits$bulk.TPO, n = nrow(fits$bulk.TPO)))

# select the samples
sel_samples <- dge$samples[grep(pattern = "bulk$", dge$samples$SampleGroup),]

# select the genes
genes <- dge$genes
sel_genes <- genes[genes$SYMBOL %in% top100_order, ]
sel_genes$ensembl_id <- rownames(sel_genes)
rownames(sel_genes) = sel_genes$SYMBOL
sel_genes = sel_genes[top100_order,]

# select the genes
logcpm_sel <- logcpm[sel_genes$ensembl_id, rownames(sel_samples)]

# rename the genes
rownames(logcpm_sel) = sel_genes$SYMBOL

# create dataframe for the sample
sample_col <- sel_samples %>% select(SampleGroup)

# draw the zscore heatmap 
heatmap <- pheatmap(logcpm_sel[, c("BSSE_QGF_200814", "BSSE_QGF_200813", "BSSE_QGF_200795", "BSSE_QGF_200805", "BSSE_QGF_200812", "BSSE_QGF_200816")],
                    color = colorRampPalette(c("blue2","white","red"))(100),
                    scale = "row",
                    clustering_method = "ward.D2", 
                    annotation_col = sample_col,
                    fontsize_row = 4, 
                    cluster_rows = F, 
                    cluster_cols = F)

# save the heatmap in the files
saveFigures(fileName = "heatmap_paper_bulk_TPOvsPBS", ggplot = heatmap, dirPlot = "DEA/Result/Result_woSamples/bulk.TPO/", A4 = T)



### Which version of R and packages was used
devtools::session_info()
capture.output(devtools::session_info(), file=paste0("session_info_analysis_", 
                                                     format(Sys.Date(), format = "%d_%m_%Y"), ".txt")) 
