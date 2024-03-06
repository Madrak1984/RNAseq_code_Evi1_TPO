###### command line to perform the differential expression analysis for the TPO-Evi1 mRNA
# done by Jonathan Seguin, DBM, UKBB
# created in Tue Mar 15 10:44:34 2022


# prepare the scripts -----------------------------------------------------


# load the libraries
library(edgeR)
library(DESeq2)
library(ggplot2)
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

# Palette with 20 colors (+black and white) that do not conflict so much, adapted from https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
myPalette2 <-  c('#e6194b', '#4363d8', '#3cb44b', '#984EA3', '#f58231', '#ffe119', '#F781BF', '#808080', '#98BFDB', '#bcf60c', '#008080', '#e6beff', '#E5C494', '#000075', '#CD00CD', '#aaffc3', '#808000', '#9a6324', '#fffac8', '#800000', '#000000', '#ffffff')

### Which version of R and packages was used
devtools::session_info()
capture.output(devtools::session_info(), file=paste0("session_info_analysis_", 
                              format(Sys.Date(), format = "%d_%m_%Y"), ".txt")) 


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

# read the rds file
dge <- readRDS("Snakemake/ensdb_102_dge_list.rds")
se <-  readRDS("Snakemake/ensdb_102_summarized_experiment.rds")

# export the tables in the snakemake directory
genes <- cbind("genes" = rownames(as.data.frame(dge$genes)), as.data.frame(dge$genes))
count <- cbind("genes" = rownames(as.data.frame(dge$counts)), as.data.frame(dge$counts))

list_out <- list("count" = count, "genes" = genes, "samples" = as.data.frame(dge$samples))
write_xlsx(list_out, path = "Snakemake/tables_dge.xlsx")

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


# PCA only in bulks,  GFP+ or GFP-  ---------------------------------------

# select the bulk samples
id_bulk <- grep(pattern = "bulk", dge$samples$SampleGroup)

# remove the genes without read 
dge_bulk <- dge[, id_bulk]

# save the PCAs
savePCA(dge_tmp = dge_bulk, extensionName = "_bulk",  plotDir = plotDir)

# select the bulk samples
id_bulk_wosamples <- setdiff(id_bulk, grep(pattern = "624-2", 
                                           dge$samples$ExternalSampleName))
id_bulk_wosamples <- setdiff(id_bulk_wosamples, grep(pattern = "35088",
                                          dge$samples$ExternalSampleName))

# remove the genes without read 
dge_bulk_wosamples <- dge[, id_bulk_wosamples]

# save the PCAs
savePCA(dge_tmp = dge_bulk_wosamples, extensionName = "_bulk_wosamples", 
        plotDir = plotDir)



#### select the GFP samples
id_GFP <- setdiff(grep(pattern = "GFP", dge$samples$SampleGroup), 
                  grep(pattern = "GFPneg", dge$samples$SampleGroup))

# remove the genes without read 
dge_GFP <- dge[, id_GFP]

# save the PCAs
savePCA(dge_tmp = dge_GFP, extensionName = "_GFP",  plotDir = plotDir)




#### select the GFP samples
id_GFPneg <-  grep(pattern = "GFPneg", dge$samples$SampleGroup)

# remove the genes without read 
dge_GFPneg <- dge[, id_GFPneg]

# save the PCAs
savePCA(dge_tmp = dge_GFPneg, extensionName = "_GFPneg",  plotDir = plotDir)




# Perform the differential expression analysis ----------------------------

# create the model matrix
moma <- model.matrix(~ 0 + SampleGroup, data=dge$samples) # ~ 0 + SampleGroup
colnames(moma) <- sub("SampleGroup", "", colnames(moma))

## filtering of lowly expressed genes                                                                                       
## based on the figures plotDensity above                                                                                               
dge <- dge[rowSums(log2cpm > 1) >= 2,]                                                                                     
dge$samples$lib.size <- colSums(dge$counts)                                                                                 
dge <- calcNormFactors(dge)                                                                                                
dge <- estimateDisp(dge, design=moma)   
saveRDS(dge, file= file.path(RDSDir,"EnsemblGenes.DGEList.filtered.rds"))  

# estimating the dispersion
dge$common.dispersion
png(file.path(plotDir, "plotBCV_dispersion.png"))
plotBCV(dge)
dev.off()


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


############## Perform the DEA

# calculate the fit
fit <- glmQLFit(dge, design = moma)

# create the pdf 
pdf(file.path(plotDir, "DE_analysis_glmQLFit_DEA.pdf"), width = 6, height = 6)

# prepare the variable
decideTestsDGE_glmQLFit = NULL
fits <- list()
for (i in 1:ncol(contrasts.matrix)) {
  
  # select the contrasts
  contr.name <- colnames(contrasts.matrix)[i]
  message("   Working on ... ", contr.name)
  ## Subset samples
  dge.temp <- dge[, dge$samples$group %in% names(which(contrasts.matrix[,i] != 0))]
  
  ## Redo TMM normalization
  dge.temp <- calcNormFactors(dge.temp) 
  ## After TMM normalization
  log2cpm.temp <- edgeR::cpm(dge.temp, log=TRUE, normalized.lib.sizes=TRUE)
  
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
  detags <- rownames(dge)[as.logical(de)]
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
saveRDS(fits, file=file.path(RDSDir, "edgeR_fits_results_DEA.rds"))


# Analyze the DEA results -------------------------------------------------

# list the names of fits
names_fits <- names(fits)

# create the dataframe to save the count of DEGs
data_count <- data.frame("comparison" = names_fits, "down" = 0, "NS" = 0, "up" = 0)

# loop to analyze the fits list
for(name in names_fits){
  print(name)
  
  # select the dge_list
  dge_name <- fits[[name]]
  
  # select the toptags
  dge_DEG <- as.data.frame(topTags(dge_name, n = nrow(dge_name), p.value = 0.05))
  
  # count the DEG down or up
  if(nrow(dge_DEG[dge_DEG$logFC <= 0,]) > 0){
    data_count$down[which(data_count$comparison == name)] = nrow(dge_DEG[dge_DEG$logFC <= 0,])  
  }
  if(nrow(dge_DEG[dge_DEG$logFC >= 0,]) > 0){
    data_count$up[which(data_count$comparison == name)] = nrow(dge_DEG[dge_DEG$logFC >= 0,])  
  }
  print(paste0("count for the dge_DEG ", nrow(dge_name)))
  data_count$NS[which(data_count$comparison == name)] = nrow(dge_name) - 
    (data_count$down[which(data_count$comparison == name)] + 
       data_count$up[which(data_count$comparison == name)])
}


# create the output directory
outDir <- "DEA/Result"
if(!dir.exists(outDir)) dir.create(outDir, recursive = T)

# save the excel files
write.xlsx2(data_count, file = file.path(outDir, "count_DEGs.xlsx"),
            col.names = T, row.names = F)

# list the names of fits
names_fits <- names(fits)

# load the DE function
source("Src/DE_functions.R")

# save the DEG in excel files
for(name in names_fits){
  
  # print the message
  print(paste0(date(), ": output for ", name))
  
  # create the tmp_dir
  tmp_dir <- file.path(outDir, name)
  if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = T)
  
  # select the dge_list
  dge_name <- fits[[name]]
  
  # select the toptags
  dge_DEG <- as.data.frame(topTags(dge_name, n = nrow(dge_name)))
  
  # count the DEG down or up
  write.xlsx2(dge_DEG, file = file.path(tmp_dir, paste0("DEGs_", name, ".xlsx")), 
              col.names = T, row.names = F)
  gc()
  
  
  # select the dge_list
  volcano <- getVolcanoplot(dgeobject = fits[[name]], nbGenes = 10)
  saveFigures(fileName = paste0("volcanoplot_", name), dirPlot = tmp_dir,
              ggplot = volcano)
  
  # special command lines only for TPO: GFP+ vs GFP-
  # dge_DEG <- as_tibble(dge_DEG)
  # genes <- dge_DEG %>% arrange(FDR) %>% select(SYMBOL) %>% slice(1:10)
  # genes <- genes$SYMBOL
  # volcano <- getVolcanoplot(dgeobject = fits[[name]], genes = genes)
  # tmp_dir <- "DEA/Result/Result_woSamples/TPO.GFPvsGFPneg/"
  # saveFigures(fileName = "volcanoplot_TPO.GFPvsGFPneg_10genes", dirPlot = tmp_dir,
    #          ggplot = volcano)
  
  # select the significant toptags
  dge_DEG <- as.data.frame(topTags(dge_name, n = nrow(dge_name), p.value = 0.05))
  
  # perform the enrichment analysis
  print(paste0(date(), ": start enrichment analysis for ", name))
  
  if(!is.null(dge_DEG)){
    if(nrow(dge_DEG) > 0){
      EnrichmentTerm_Analysis(common = rownames(dge_DEG), outputDir = tmp_dir, 
                              mart = mart)
    }
  } 
  print(paste0(date(), ": end enrichment analysis for ", name))
  
}




# PCA without the 35088 and 624-2 samples ---------------------------------

# select the identifiers after removing the 35088 and 624-2 samples
id_samples <-  intersect(grep(pattern = "35088", dge$samples$ExternalSampleName,
                            invert = T), grep(pattern = "624-2", dge$samples$ExternalSampleName, 
                                               invert = T))
#id_samples <- grep(pattern = "624-2", dge$samples$ExternalSampleName, 
                                               # invert = T)

# remove the genes without read 
dge_samples <- dge[, id_samples]

# update the dge samples
dge_samples$samples$SampleGroup <- gsub(pattern = "noTPO", replacement = "PBS",
                                        dge_samples$samples$SampleGroup)
#dge_samples$samples$SampleGroup <- gsub(pattern = "GFP$", replacement = "GFP+",
                                        #dge_samples$samples$SampleGroup)
#dge_samples$samples$SampleGroup <- gsub(pattern = "GFPneg", replacement = "GFP-",
 #                                       dge_samples$samples$SampleGroup)

# save the PCAs
#savePCA(dge_tmp = dge_samples, extensionName = "_wosamples",  plotDir = plotDir)
savePCA(dge_tmp = dge_samples, extensionName = "_wosamples_paper",  plotDir = plotDir, 
        theme_bw = T, text = F)

# see to create the 3d PCA for the PCA
dge_tmp = dge_samples[which(rowSums(dge_samples$counts) > 10), ]

# normalization and estimate dispersion
df.dge <- calcNormFactors(dge_tmp)
moma <- model.matrix(~ 0 + SampleGroup, data=dge_tmp$samples) # ~ 0 + SampleGroup
colnames(moma) <- sub("SampleGroup", "", colnames(moma))
dge_tmp <- calcNormFactors(dge_tmp)
dge_tmp <- estimateDisp(dge_tmp)

# plot a PCA with Deseq2
dds = DESeqDataSetFromMatrix(countData = round(df.dge$counts),
                             colData = df.dge$samples, design = ~SampleGroup)
transform.dds <- rlog(dds, blind = T)

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(transform.dds)))

# save the barplot with the observed batch effect 
saveFigures(fileName = "PC_dimension_HE_paper", ggplot = fviz_eig((pca)), dirPlot = plotDir,
            A4 = T)

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

# assembly the data for the plot
d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], group=dge_samples$samples$SampleGroup,  
                name=colnames(transform.dds))

# open 3d window
open3d()

# resize window
par3d(windowRect = c(100, 100, 612, 612))

# create the dataframe for the colours
dt_colours <-  data.frame("group" = unique(d$group),
                          "colour" = brewer.pal(length(unique(d$group)), 
                                                name = "Set1"))

# merge with the colour
d = merge(x=d, y = dt_colours, by.x = "group", by.y = "group")

# plot in 3D
plot3d(x=d$PC1, y=d$PC2, z=d$PC3, size = 1, type="s", col = d$colour,
       xlab = "PC1", ylab = "PC2", zlab = "PC3")

# add the legend
legend3d("topright", legend = dt_colours$group, pch = 16, col= dt_colours$colour,
         cex= 1, inset=c(0.02))

# add treatment information
d$treatment = "PBS"
d$treatment[grep(pattern = "TPO", d$group)] = "TPO"  

# draw the pcaplot
pcaplot <- ggplot(data=d, aes_string(x="PC2", y="PC3", color="group", shape = "treatment")) + 
  geom_point(size=3) + 
  xlab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) +
  ylab(paste0("PC3: ",round(percentVar[3]*100),"% variance")) +
  coord_fixed() + scale_colour_brewer(palette = "Set1") +
  theme_bw()

# save the PCAplot with the observed batch effect 
saveFigures(fileName = "PCA_PC2-PC3_group_treatment_HE_paper", 
            ggplot = pcaplot, dirPlot = "DEA/Plots/")


# draw the pcaplot
pcaplot <- ggplot(data=d, aes_string(x="PC1", y="PC2", color="group", shape = "treatment")) + 
  geom_point(size=3) + 
  xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) +
  coord_fixed() + scale_colour_brewer(palette = "Set1") +
  theme_bw()

# save the PCAplot with the observed batch effect 
saveFigures(fileName = "PCA_PC1-PC2_group_treatment_HE_paper", 
            ggplot = pcaplot, dirPlot = "DEA/Plots/")


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

# load the DE function
source("Src/DE_functions.R")

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
common_DEG_bulk_GFP <- intersect(rownames(DEG_GFP), rownames(DEG_bulk))
spe_DEG_GFP <- setdiff(rownames(DEG_GFP), rownames(DEG_bulk))
spe_DEG_bulk <- setdiff(rownames(DEG_bulk), rownames(DEG_GFP))

# create the tmp_dir
tmp_dir <- file.path(outDir, "common_DEG_GFP+_bulk")
if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = T)

# count the DEG down or up
write.xlsx2(genes_mouse[genes_mouse$ensembl_gene_id %in% common_DEG_bulk_GFP, ], file = file.path(tmp_dir, paste0("common_DEG_GFP+_bulk.xlsx")), col.names = T, row.names = F)
gc()


# select all the analyzed genes
universe <- rownames(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO)))

# perform the term analysis
EnrichmentTerm_Analysis(common = common_DEG_bulk_GFP, outputDir = tmp_dir, mart = mart, universe = universe)


# create the tmp_dir
tmp_dir <- file.path(outDir, "spe_DEG_GFP+_bulk")
if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = T)

# count the DEG down or up
write.xlsx2(genes_mouse[genes_mouse$ensembl_gene_id %in% spe_DEG_GFP, ], file = file.path(tmp_dir, paste0("spe_DEG_GFP+_bulk.xlsx")), col.names = T, row.names = F)
gc()


# select all the analyzed genes
universe <- rownames(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO)))

# perform the term analysis
EnrichmentTerm_Analysis(common = spe_DEG_GFP, outputDir = tmp_dir, mart = mart, universe = universe)

# create the tmp_dir
tmp_dir <- file.path(outDir, "spe_DEG_bulk")
if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = T)

# count the DEG down or up
write.xlsx2(genes_mouse[genes_mouse$ensembl_gene_id %in% spe_DEG_bulk, ], file = file.path(tmp_dir, paste0("spe_DEG_bulk.xlsx")), col.names = T, row.names = F)
gc()


# select all the analyzed genes
universe <- rownames(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO)))

# perform the term analysis
EnrichmentTerm_Analysis(common = spe_DEG_bulk, outputDir = tmp_dir, mart = mart, universe = universe)


## only common up DEGs ---------------------------------

# find common genes 
common_DEG_bulk_GFP_up <- intersect(rownames(DEG_GFP[DEG_GFP$logFC > 0, ]), 
                                 rownames(DEG_bulk[DEG_bulk$logFC > 0, ]))
spe_DEG_GFP_up <- setdiff(rownames(DEG_GFP[DEG_GFP$logFC > 0, ]), 
                            rownames(DEG_bulk[DEG_bulk$logFC > 0, ]))
spe_DEG_bulk_up <- setdiff(rownames(DEG_bulk[DEG_bulk$logFC > 0, ]),
                           rownames(DEG_GFP[DEG_GFP$logFC > 0, ]))

# create the tmp_dir
tmp_dir <- file.path(outDir, "common_DEG_GFP+_bulk_up")
if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = T)

# count the DEG down or up
write.xlsx2(genes_mouse[genes_mouse$ensembl_gene_id %in% common_DEG_bulk_GFP_up, ], file = file.path(tmp_dir, paste0("common_DEG_GFP+_bulk_up.xlsx")), col.names = T, row.names = F)
gc()


# select all the analyzed genes
universe <- rownames(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO)))

# perform the term analysis
EnrichmentTerm_Analysis(common = common_DEG_bulk_GFP_up, outputDir = tmp_dir, mart = mart, universe = universe)
#redo without universe
#EnrichmentTerm_Analysis(common = common_DEG_bulk_GFP_up, outputDir = file.path(tmp_dir, "WoUniverse"), mart = mart)
#redo without universe
#EnrichmentTerm_Analysis(common = common_DEG_bulk_GFP_up, outputDir = file.path(tmp_dir, "WoSimplify"), 
#                        mart = mart, universe = universe, simplify = F)

# create the tmp_dir
tmp_dir <- file.path(outDir, "spe_DEG_GFP+_bulk_up")
if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = T)

# count the DEG down or up
write.xlsx2(genes_mouse[genes_mouse$ensembl_gene_id %in% spe_DEG_GFP_up, ], file = file.path(tmp_dir, paste0("spe_DEG_GFP+_bulk_up.xlsx")), col.names = T, row.names = F)
gc()


# select all the analyzed genes
universe <- rownames(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO)))

# perform the term analysis
EnrichmentTerm_Analysis(common = spe_DEG_GFP_up, outputDir = tmp_dir, mart = mart, universe = universe)


# create the tmp_dir
tmp_dir <- file.path(outDir, "spe_DEG_bulk_up")
if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = T)

# count the DEG down or up
write.xlsx2(genes_mouse[genes_mouse$ensembl_gene_id %in% spe_DEG_bulk_up, ], file = file.path(tmp_dir, paste0("spe_DEG_bulk_up.xlsx")), col.names = T, row.names = F)
gc()


# select all the analyzed genes
universe <- rownames(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO)))

# perform the term analysis
EnrichmentTerm_Analysis(common = spe_DEG_bulk_up, outputDir = tmp_dir, mart = mart, universe = universe)
  




## only common down DEGs ---------------------------------

# find common genes 
common_DEG_bulk_GFP_down <- intersect(rownames(DEG_GFP[DEG_GFP$logFC < 0, ]), 
                                    rownames(DEG_bulk[DEG_bulk$logFC < 0, ]))
spe_DEG_GFP_down <- setdiff(rownames(DEG_GFP[DEG_GFP$logFC < 0, ]), 
                               rownames(DEG_bulk[DEG_bulk$logFC < 0, ]))
spe_DEG_bulk_down <- setdiff(rownames(DEG_bulk[DEG_bulk$logFC < 0, ]),
                             rownames(DEG_GFP[DEG_GFP$logFC < 0, ]))

# create the tmp_dir
tmp_dir <- file.path(outDir, "common_DEG_GFP+_bulk_down")
if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = T)

# count the DEG down or up
write.xlsx2(genes_mouse[genes_mouse$ensembl_gene_id %in% common_DEG_bulk_GFP_down, ], file = file.path(tmp_dir, paste0("common_DEG_GFP+_bulk_down.xlsx")), col.names = T, row.names = F)
gc()


# select all the analyzed genes
universe <- rownames(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO)))

# perform the term analysis
EnrichmentTerm_Analysis(common = common_DEG_bulk_GFP_down, outputDir = tmp_dir, mart = mart, universe = universe)

# create the tmp_dir
tmp_dir <- file.path(outDir, "spe_DEG_GFP+_bulk_down")
if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = T)

# count the DEG down or up
write.xlsx2(genes_mouse[genes_mouse$ensembl_gene_id %in% spe_DEG_GFP_down, ], file = file.path(tmp_dir, paste0("spe_DEG_GFP+_bulk_down.xlsx")), col.names = T, row.names = F)
gc()


# select all the analyzed genes
universe <- rownames(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO)))

# perform the term analysis
EnrichmentTerm_Analysis(common = spe_DEG_GFP_down, outputDir = tmp_dir, mart = mart, universe = universe)



# create the tmp_dir
tmp_dir <- file.path(outDir, "spe_DEG_bulk_down")
if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = T)

# count the DEG down or up
write.xlsx2(genes_mouse[genes_mouse$ensembl_gene_id %in% spe_DEG_bulk_down, ], file = file.path(tmp_dir, paste0("spe_DEG_bulk_down.xlsx")), col.names = T, row.names = F)
gc()


# select all the analyzed genes
universe <- rownames(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO)))

# perform the term analysis
EnrichmentTerm_Analysis(common = spe_DEG_bulk_down, outputDir = tmp_dir, mart = mart, universe = universe)



# PCA without the 624-2 samples ---------------------------------

# select the identifiers after removing the 35088 and 624-2 samples
id_wo624 <-  grep(pattern = "624-2", dge$samples$ExternalSampleName, 
                                                invert = T)

# remove the genes without read 
dge_wo624 <- dge[, id_wo624]

# save the PCAs
savePCA(dge_tmp = dge_wo624, extensionName = "_wo624",  plotDir = plotDir)



# PCA with only TPO ---------------------------------

# select the identifiers after removing the 35088 and 624-2 samples
id_TPO <-  grep(pattern = "noTPO", dge$samples$SampleGroup, 
                  invert = T)

# remove the genes without read 
dge_TPO <- dge[, id_TPO]

# save the PCAs
savePCA(dge_tmp = dge_TPO, extensionName = "_TPO",  plotDir = plotDir)

# PCA with TPO and without 624
id_TPO_wo624 <- setdiff(id_TPO, grep(pattern = "624-2", 
                                     dge$samples$ExternalSampleName))

# remove the genes without read 
dge_TPO_wo624 <- dge[, id_TPO_wo624]

# save the PCAs
savePCA(dge_tmp = dge_TPO_wo624, extensionName = "_TPO_wo624",  plotDir = plotDir)




# PCA without TPO ---------------------------------

# select the identifiers after removing the 35088 and 624-2 samples
id_woTPO <-  grep(pattern = "noTPO", dge$samples$SampleGroup)

# remove the genes without read 
dge_woTPO <- dge[, id_woTPO]

# save the PCAs
savePCA(dge_tmp = dge_woTPO, extensionName = "_woTPO",  plotDir = plotDir)

# select the identifiers after removing the 35088 and 624-2 samples
id_woTPO_wo35088 <-  setdiff(id_woTPO, grep(pattern = "35088", 
                                            dge$samples$ExternalSampleName))

# remove the genes without read 
dge_woTPO_wo35088 <- dge[, id_woTPO_wo35088]

# save the PCAs
savePCA(dge_tmp = dge_woTPO_wo35088, extensionName = "_woTPO_wo35088",
        plotDir = plotDir)




# PCA with  GFP+ and GFP- ---------------------------------

# select the identifiers after removing the 35088 and 624-2 samples
id_allGFP <-  grep(pattern = "GFP", dge$samples$ExternalSampleName)

# remove the genes without read 
dge_allGFP <- dge[, id_allGFP]

# save the PCAs
savePCA(dge_tmp = dge_allGFP, extensionName = "_allGFP",  plotDir = plotDir)

# PCA with TPO and without 624
id_allGFP_wo624 <- setdiff(id_allGFP, grep(pattern = "624-2", 
                                          dge$samples$ExternalSampleName))

# remove the genes without read 
dge_allGFP_wo624 <- dge[, id_allGFP_wo624]

# save the PCAs
savePCA(dge_tmp = dge_allGFP_wo624, extensionName = "_allGFP_wo624",  
        plotDir = plotDir)

# PCA with TPO and without 624
id_allGFP_wosamples <- setdiff(id_allGFP_wo624, grep(pattern = "35088",
                                        dge$samples$ExternalSampleName))

# remove the genes without read 
dge_allGFP_wosamples <- dge[, id_allGFP_wosamples]

# save the PCAs
savePCA(dge_tmp = dge_allGFP_wosamples, extensionName = "_allGFP_wosamples", 
        plotDir = plotDir)




# draw the plot for the MECOM expression ----------------------------------

# create the barplot
createPlots(pattern = "mecom", dge = dge, plotDir = plotDir)




# draw the plot for the CD52 expression ----------------------------------

# create the barplot
createPlots(pattern = "cd52", dge = dge, plotDir = plotDir)



# metatable creation ------------------------------------------------------

fits <-  readRDS(file = "DEA/Rds/RDS_woD/edge_samplesR_fits_results_DEA.rds")
# check the count of DEG
for(name in names(fits)){
  
  print(name)
  
  # select the count of fits
  count <- as.data.frame(topTags(fits[[name]], n = nrow(fits[[name]]), p.value = 0.05))
  print(nrow(count))
  
  
}



# write the metatable for fdr <= 5% in the excel file
write_xlsx(unique(createMetatable(fits = fits)), path = file.path(outDir, "metatable_DEGs_FDR5%.xlsx"),
           col_name = T)
#write.xlsx2(metatable, file = file.path(outDir, "metatable_DEGs_FDR5%.xlsx"),
            #col.names = T, row.names = F) --> problem to save the folder

# write the metatable for fdr <= 10% in the excel file
write_xlsx(unique(createMetatable(fits = fits, FDR = 0.1)), path = file.path(outDir,
                                  "metatable_DEGs_FDR10%.xlsx"), col_name = T)

# write the metatable for pval <= 5% in the excel file
write_xlsx(unique(createMetatable_pval(fits = fits)), path = file.path(outDir,
                               "metatable_DEGs_pval5%.xlsx"), col_name = T)


# Term enrichment analyis for noTPO: GFP+ vs bulk ---------------------------------------------

# create the tmp_dir for pval <= 5%
outDir = "DEA/Result/Result_woSamples"
pval_noTPO <- file.path(outDir, "noTPO.GFPvsbulk/pval")
if(!dir.exists(pval_noTPO)) dir.create(pval_noTPO, recursive = T)

# select the dge_samples_list
dge_samples_name <- fits$noTPO.GFPvsbulk
  
# select the toptags
dge_samples_DEG <- as.data.frame(topTags(dge_samples_name, n = nrow(dge_samples_name)))
  
# select the DEG based on pval <= 5%
dge_samples_DEG = dge_samples_DEG[dge_samples_DEG$PValue <= 0.05,]

# perform the enrichmentTerm analysis  
EnrichmentTerm_Analysis(common = rownames(dge_samples_DEG), outputDir = pval_noTPO, mart = mart)



# scatterplot FC GFP+ vs FC GFP- for TPO+ vs noTPO comparison -------------------

# create a temporary table for the DEGs to use for the scatterplot
tmp_data = metatable[,c(2, 11:13, 18:20)]

# change color according to the FDR
tmp_data$color = "black"
tmp_data$color[tmp_data$`FDR_GFP+: TPO vs noTPO` <= 0.05] = "blue"
tmp_data$color[tmp_data$`FDR _GFP-: TPO vs noTPO` <= 0.05] = "red"
tmp_data$color[(tmp_data$`FDR _GFP-: TPO vs noTPO` <= 0.05) &
                 (tmp_data$`FDR_GFP+: TPO vs noTPO` <= 0.05)] = "purple"


# create the scatterplot
scatterplot <- ggplot(data = tmp_data, aes(x=`logFC _GFP-: TPO vs noTPO`, 
               y=`logFC_GFP+: TPO vs noTPO`, text = SYMBOL)) + 
              geom_point(aes(color = color)) + theme_bw() +
              scale_color_manual(values = c("black", "blue", "purple", "red"),
                                name = "TPO vs noTPO comparison\n(FDR <= 5%)", 
                                breaks = c("black", "blue", "purple", "red"), 
                                labels = c("N.S.", "significant in GFP+", 
                                    "significant in GFP+ & GFP-" ,
                                    "significant in GFP-"))+ xlab("logFC in GFP-") +
            ylab("logFC in GFP+") + theme(legend.title = element_text(face = "bold"))
         # geom_text_repel(data=subset(tmp_data[tmp_data$SYMBOL %in% c("Mecom", "Dsp", "Mpl", "Stat5a", "Osm", "Cish"),]), 
          #                color="darkgreen", aes(label=SYMBOL), size =5)

# display interactive scatterplot
ggplotly(scatterplot)

# save the plot in files
saveFigures(fileName = "scatterplot_TPO_comp_GFP", ggplot = scatterplot, 
            dirPlot = "DEA/Result/Result_woSamples/", A4 = T)


# select genes to remove immunoglobulin and predicted genes
genes <- c(metatable$SYMBOL[grep(pattern = "^immunoglobulin", metatable$DESCRIPTION)], 
           metatable$SYMBOL[grep(pattern = "^predicted", metatable$DESCRIPTION)])
genes <- c(genes, c("Cxcl12", "Gimap3", "Gimap4", "Gimap6", "Thy1", "Lck", 
                    "Ablim1", "Il2rb"))


# remove the selected genes
tmp_data <- tmp_data[!tmp_data$SYMBOL %in% genes, ]


# create the scatterplot
scatterplot <- ggplot(data = tmp_data, aes(x=`logFC _GFP-: TPO vs noTPO`, 
                                           y=`logFC_GFP+: TPO vs noTPO`, text = SYMBOL)) + 
  geom_point(aes(color = color)) + theme_bw() +
  scale_color_manual(values = c("black", "blue", "purple", "red"),
                     name = "TPO vs noTPO comparison\n(FDR <= 5%)", 
                     breaks = c("black", "blue", "purple", "red"), 
                     labels = c("N.S.", "significant in GFP+", 
                                "significant in GFP+ & GFP-" ,
                                "significant in GFP-"))+ xlab("logFC in GFP-") +
  ylab("logFC in GFP+") + theme(legend.title = element_text(face = "bold"))
# geom_text_repel(data=subset(tmp_data[tmp_data$SYMBOL %in% c("Mecom", "Dsp", "Mpl", "Stat5a", "Osm", "Cish"),]), 
#                color="darkgreen", aes(label=SYMBOL), size =5)

# display interactive scatterplot
ggplotly(scatterplot)

# save the plot in files
saveFigures(fileName = "scatterplot_TPO_comp_GFP_woIG_GM", ggplot = scatterplot, 
            dirPlot = "DEA/Result/Result_woSamples/", A4 = T)


# create the scatterplot for Evi1 according to the Hoxa9
scatterplot <- ggplot(data = tmp_data, aes(x=`logFC _GFP-: TPO vs noTPO`, 
                                           y=`logFC_GFP+: TPO vs noTPO`, text = SYMBOL)) + 
  geom_point(aes(color = color)) + theme_bw() +
  scale_color_manual(values = c("black", "blue", "purple", "red"),
                     name = "TPO vs noTPO comparison\n(FDR <= 5%)", 
                     breaks = c("black", "blue", "purple", "red"), 
                     labels = c("N.S.", "significant in GFP+", 
                                "significant in GFP+ & GFP-" ,
                                "significant in GFP-"))+ xlab("logFC in GFP-") +
  ylab("logFC in GFP+") + theme(legend.title = element_text(face = "bold"))
# geom_text_repel(data=subset(tmp_data[tmp_data$SYMBOL %in% c("Mecom", "Dsp", "Mpl", "Stat5a", "Osm", "Cish"),]), 
#                color="darkgreen", aes(label=SYMBOL), size =5)

# display interactive scatterplot
ggplotly(scatterplot)

# save the plot in files
saveFigures(fileName = "scatterplot_TPO_comp_GFP_woIG_GM", ggplot = scatterplot, 
            dirPlot = "DEA/Result/Result_woSamples/", A4 = T)



# scatterplot FC TPO vs FC noTPO for GFP+ vs GFP- comparison -------------------

# read the TPO GFP+ vs GFP- comparison (file created in differential_expression_analysis_TPO)
TPO_comp <- read.xlsx2(file = "DEA/Result/TPO_wo624/TPO.GFPvsGFPneg/DEGs_TPO.GFPvsGFPneg.xlsx", 
                       header = T, as.data.frame = T, sheetIndex = 1)

# read the TPO GFP+ vs GFP- comparison (file created in differential_expression_analysis_noTPO)
noTPO_comp <- read.xlsx2(file = "DEA/Result/noTPO/noTPO.GFPvsGFPneg/DEGs_noTPO.GFPvsGFPneg.xlsx", 
                       header = T, as.data.frame = T, sheetIndex = 1)

# create a temporary table for the DEGs to use for the scatterplot
tmp_data = merge(x = TPO_comp, y = noTPO_comp, by = "SYMBOL", suffixes = c("_TPO", "_noTPO"))

# change color according to the FDR
tmp_data$color = "black"
tmp_data$color[tmp_data$FDR_TPO <= 0.05] = "blue"
tmp_data$color[tmp_data$FDR_noTPO <= 0.05] = "red"
tmp_data$color[(tmp_data$FDR_TPO <= 0.05) &
                 (tmp_data$FDR_noTPO <= 0.05)] = "purple"

# reduce the tmp_data table
tmp_data = tmp_data[, c(1, 8, 19, 24)]
tmp_data[,2:3] = apply(tmp_data[,2:3], 2, as.numeric)

# create the scatterplot
scatterplot <- ggplot(data = tmp_data, aes(x=logFC_noTPO, y=logFC_TPO, 
  text = SYMBOL)) + geom_point(aes(color = color)) + theme_bw() +
  scale_color_manual(values = c("black", "blue", "purple", "red"),
                     name = "GFP+ vs GFP- comparison\n(FDR <= 5%)", 
                     breaks = c("black", "blue", "purple", "red"), 
                     labels = c("N.S.", "significant in TPO", 
                                "significant in TPO & no TPO" ,
                                "significant in no TPO"))+ xlab("logFC in no TPO") +
  ylab("logFC in TPO") + theme(legend.title = element_text(face = "bold"))
# geom_text_repel(data=subset(tmp_data[tmp_data$SYMBOL %in% c("Mecom", "Dsp", "Mpl", "Stat5a", "Osm", "Cish"),]), 
#                color="darkgreen", aes(label=SYMBOL), size =5)

ggplotly(scatterplot)

# save the plot in files
saveFigures(fileName = "scatterplot_GFP_comp_TPO", ggplot = scatterplot, 
            dirPlot = "DEA/Result/", A4 = T)




# scatterplot FC TPO vs FC noTPO for GFP+ vs GFP- comparison  without immunoglobulin genes -------------------

# read the TPO GFP+ vs GFP- comparison (file created in differential_expression_analysis_TPO)
TPO_comp_woIg <- read.xlsx2(file = "DEA/Result/TPO_wo624_woIg/TPO.GFPvsGFPneg/DEGs_TPO.GFPvsGFPneg.xlsx", 
                       header = T, as.data.frame = T, sheetIndex = 1)

# read the TPO GFP+ vs GFP- comparison (file created in differential_expression_analysis_noTPO)
noTPO_comp_woIg <- read.xlsx2(file = "DEA/Result/noTPO_woIg/noTPO.GFPvsGFPneg/DEGs_noTPO.GFPvsGFPneg.xlsx", 
                         header = T, as.data.frame = T, sheetIndex = 1)

# create a temporary table for the DEGs to use for the scatterplot
tmp_data = merge(x = TPO_comp_woIg, y = noTPO_comp_woIg, by = "SYMBOL", suffixes = c("_TPO", "_noTPO"))

# change color according to the FDR
tmp_data$color = "black"
tmp_data$color[tmp_data$FDR_TPO <= 0.05] = "blue"
tmp_data$color[tmp_data$FDR_noTPO <= 0.05] = "red"
tmp_data$color[(tmp_data$FDR_TPO <= 0.05) &
                 (tmp_data$FDR_noTPO <= 0.05)] = "purple"

# reduce the tmp_data table
tmp_data = tmp_data[, c(1, 8, 19, 24)]
tmp_data[,2:3] = apply(tmp_data[,2:3], 2, as.numeric)


# remove the selected genes
tmp_data <- tmp_data[!tmp_data$SYMBOL %in% genes, ]

# create the scatterplot
scatterplot <- ggplot(data = tmp_data, aes(x=logFC_noTPO, y=logFC_TPO, 
                                           text = SYMBOL)) + geom_point(aes(color = color)) + theme_bw() +
  scale_color_manual(values = c("black", "blue", "purple", "red"),
                     name = "GFP+ vs GFP- comparison\n(FDR <= 5%)", 
                     breaks = c("black", "blue", "purple", "red"), 
                     labels = c("N.S.", "significant in TPO", 
                                "significant in TPO & no TPO" ,
                                "significant in no TPO"))+ xlab("logFC in no TPO") +
  ylab("logFC in TPO") + theme(legend.title = element_text(face = "bold"))
# geom_text_repel(data=subset(tmp_data[tmp_data$SYMBOL %in% c("Mecom", "Dsp", "Mpl", "Stat5a", "Osm", "Cish"),]), 
#                color="darkgreen", aes(label=SYMBOL), size =5)

ggplotly(scatterplot)

# save the plot in files
saveFigures(fileName = "scatterplot_GFP_comp_TPO_woIg", ggplot = scatterplot, 
            dirPlot = "DEA/Result/", A4 = T)




# scatterplot TPO comparisons vs GFP comparisons --------------------------

# select DEG
DEG_GFP_comp <- as.data.frame(topTags(fits$GFPvsGFPneg, n = nrow(fits$GFPvsGFPneg)))
DEG_TPO_comp <- as.data.frame(topTags(fits$allGFP.TPOvsTPOneg,
                                      n = nrow(fits$allGFP.TPOvsTPOneg)))

# create a temporary table for the DEGs to use for the scatterplot
tmp_data = merge(x = DEG_GFP_comp, y = DEG_TPO_comp, by = "SYMBOL", suffixes = c("_GFP", "_TPO"))

# change color according to the FDR
tmp_data$color = "black"
tmp_data$color[tmp_data$FDR_TPO <= 0.05] = "blue"
tmp_data$color[tmp_data$FDR_GFP <= 0.05] = "red"
tmp_data$color[(tmp_data$FDR_TPO <= 0.05) &
                 (tmp_data$FDR_GFP <= 0.05)] = "purple"

# reduce the tmp_data table
tmp_data = tmp_data[, c(1, 8, 19, 24)]
tmp_data[,2:3] = apply(tmp_data[,2:3], 2, as.numeric)

# create the scatterplot
scatterplot <- ggplot(data = tmp_data, aes(x=logFC_GFP, y=logFC_TPO, 
                                           text = SYMBOL)) + geom_point(aes(color = color)) + theme_bw() +
  scale_color_manual(values = c("black", "blue", "purple", "red"),
                     name = "comparison\n(FDR <= 5%)", 
                     breaks = c("black", "blue", "purple", "red"), 
                     labels = c("N.S.", "significant in TPO+ vs TPO-", 
                                "significant in all comparisons" ,
                                "significant in GFP+ vs GFP-"))+ xlab("logFC in GFP+ vs GFP-") +
  ylab("logFC in TPO+ vs TPO-") + theme(legend.title = element_text(face = "bold"))
# geom_text_repel(data=subset(tmp_data[tmp_data$SYMBOL %in% c("Mecom", "Dsp", "Mpl", "Stat5a", "Osm", "Cish"),]), 
#                color="darkgreen", aes(label=SYMBOL), size =5) 

ggplotly(scatterplot)

# save the plot in files
saveFigures(fileName = "scatterplot_TPOvsGFP_comp", ggplot = scatterplot, 
            dirPlot = "DEA/Result/", A4 = T)


# remove the selected genes
tmp_data <- tmp_data[!tmp_data$SYMBOL %in% genes, ]

# create the scatterplot
scatterplot <- ggplot(data = tmp_data, aes(x=logFC_GFP, y=logFC_TPO, 
                                           text = SYMBOL)) + geom_point(aes(color = color)) + theme_bw() +
  scale_color_manual(values = c("black", "blue", "purple", "red"),
                     name = "comparison\n(FDR <= 5%)", 
                     breaks = c("black", "blue", "purple", "red"), 
                     labels = c("N.S.", "significant in TPO+ vs TPO-", 
                                "significant in all comparisons" ,
                                "significant in GFP+ vs GFP-"))+ xlab("logFC in GFP+ vs GFP-") +
  ylab("logFC in TPO+ vs TPO-") + theme(legend.title = element_text(face = "bold"))
# geom_text_repel(data=subset(tmp_data[tmp_data$SYMBOL %in% c("Mecom", "Dsp", "Mpl", "Stat5a", "Osm", "Cish"),]), 
#                color="darkgreen", aes(label=SYMBOL), size =5) 

ggplotly(scatterplot)

# save the plot in files
saveFigures(fileName = "scatterplot_TPOvsGFP_comp_woIg_Gm", ggplot = scatterplot, 
            dirPlot = "DEA/Result/", A4 = T)




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

  
  # find the index
  #idx <- ids2indices(get(paste0("Mm.", db)),id=metatable$ENTREZID) # see to run voom if needed
  
  # run the camera function
  #list_cameras <- c(list_cameras, list(camera(dge_samples, idx, moma)))
}

# update the name
names(list_cameras) = paste0("camera_dge_", db_cameras)

# extract information for all camera lists
for(id in names(list_cameras)){
  
  # set camera
  camera = list_cameras[[id]]
  
  # extract the significant term
  select_camera = camera[camera$FDR <= 0.05, ]
  
  # write the table
  write.table(select_camera, file = file.path("DEA/Result/Result_woSamples/", paste0(id, ".txt")), col.names = T, row.names = T, sep = "\t", quote = F)
}


# provide the test 
test$genes$ENTREZID
list_db$C2@geneSets$ABBUD_LIF_SIGNALING_1_DN

id <- which(test$genes$ENTREZID %in% list_db$C2@geneSets$ABBUD_LIF_SIGNALING_1_DN)


test_list <- list_db$C2@geneSets
for(id_list in 1:length(test_list)){
 test_list[[id_list]] = which(test$genes$ENTREZID %in% test_list[[id_list]])  
}


dge_samples <- readRDS("DEA/Rds/RDS_woD/EnsemblGenes.dge_samplesList.filtered.rds")
cpm_dge <- cpm(y = dge_samples, normalized.lib.sizes = T, log = T)
#test <- fits$GFP.TPO
# cam.BasalvsLP <- camera(v,idx,design,contrast=contr.matrix[,1])
test_list@Data
count_log <- cpm(test_list)
result_camera <- camera(cpm_dge, test_list, design = dge_samples$design, 
                        contrast = contrasts.matrix[,1])


result_camera <- camera(cpm_dge, test_list, design = dge_samples$design, 
                        contrast = contrasts.matrix[,2])


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


### perform the GSEA analysis for the C2 database
# create the gene list
FC <- getGeneListforGSEA(fit = fits$GFP.TPO)

# provide the Msigdbr database
msigdbr_species()
m_tg2 <- msigdbr(species = "Mus musculus", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)

# select the hoxa9 pathways

# perform the GSEA for Hoxa9
em2 <- GSEA(FC, TERM2GENE = m_tg2)
dotplot(em2, showCategory=30)



#### perform GSEA with selected marks from from Slany's paper (source: https://doi.org/10.1182/bloodadvances.2018025866) -----------------------------------------------------------

# load the mouse ChIPseq
mm_Chip <- read_excel(path = "External_Data/Slany_paper/ba025866-suppl2.xlsx")
gc()

## correct delay in the column by taking all ENSMUG id
# select columns with ENSEMBL Id
id_ens <- apply(mm_Chip, 2, function(x)length(grep(pattern = "ENSMUSG", x = x)) > 0)

# concatenate the gene
genes_mm_Chip <- as.vector(unlist(mm_Chip[, id_ens]))
# keep only the ensembl ID
genes_mm_Chip = unique(genes_mm_Chip[grep(pattern = "ENSMUSG", genes_mm_Chip)])

# load the human ChIPseq
#mm_Chip <- read.xlsx2(file = "External_Data/Slany_paper/ba025866-suppl2.xlsx", 
#                     sheetIndex = 1, header = T)
hs_Chip <- read_excel(path = "External_Data/Slany_paper/ba025866-suppl3.xlsx")
gc()

## correct delay in the column by taking all ENSMUG id
# select columns with ENSEMBL Id
id_ens <- apply(hs_Chip, 2, function(x)length(grep(pattern = "ENSG", x = x)) > 0)

# concatenate the gene
genes_hs_Chip <- as.vector(unlist(hs_Chip[, id_ens]))
# keep only the ensembl ID
genes_hs_Chip = unique(genes_hs_Chip[grep(pattern = "ENSG", genes_hs_Chip)])
# find orthologous genes
ortho <-  getHumanOrthologous()

# select known orthologous genes
ortho_genes <-  ortho[ortho$ensembl_gene_id %in% genes_hs_Chip, ]

# select human genes close to the promoter
genes_hs_Chip_prom <- unique(hs_Chip$`Nearest Ensembl`[grep(pattern = "promoter", hs_Chip$Annotation)])
ortho_genes_prom <-  ortho[ortho$ensembl_gene_id %in% genes_hs_Chip_prom, ]

# load the RNAseq results
hoxa9_activated <- read_excel(path = "External_Data/Slany_paper/ba025866-suppl4.xlsx", 
                              sheet = "Hoxa9_activated")
hoxa9_activated_ensID <- unique(ortho$mmusculus_homolog_ensembl_gene[ortho$mmusculus_homolog_associated_gene_name %in% hoxa9_activated$gene_short_name])
hoxa9_repressed <- read_excel(path = "External_Data/Slany_paper/ba025866-suppl4.xlsx", 
                              sheet = "Hoxa9_repressed")
hoxa9_repressed_ensID <- unique(ortho$mmusculus_homolog_ensembl_gene[ortho$mmusculus_homolog_associated_gene_name %in% hoxa9_repressed$gene_short_name])


# load the RNAseq result for interesection of genes  between Hoxa9 and MLL-ENL
inter_hoxa9 <- read_excel(path = "External_Data/Slany_paper/ba025866-suppl4.xlsx",
                          sheet = "Hoxa9vsMLLENL")
inter_hoxa9_activated = inter_hoxa9$...1[inter_hoxa9$Hoxa9 < 0]
inter_hoxa9_repressed = inter_hoxa9$...1[inter_hoxa9$Hoxa9 > 0]

# get all genes in mouse
mouse_genes <- getMouseGene()

# convert to ensemblID
inter_hoxa9_activated_ensID = unique(mouse_genes$ensembl_gene_id[mouse_genes$external_gene_name %in% inter_hoxa9_activated])
inter_hoxa9_repressed_ensID = unique(mouse_genes$ensembl_gene_id[mouse_genes$external_gene_name %in% inter_hoxa9_repressed])

# create the custom pathways
Slany_pathways <- list("Slany Hoxa9 activated" = hoxa9_activated_ensID,
                       "Slany Hoxa9 repressed" = hoxa9_repressed_ensID,
                       "Slany Hoxa9-MLLENL activated" = inter_hoxa9_activated_ensID,
                       "Slany Hoxa9-MLLENL repressed" = inter_hoxa9_repressed_ensID,
                       "Slany mRNAseq" = unique(c(hoxa9_activated_ensID, hoxa9_repressed_ensID)),
                       "Slany mRNAseq Hoxa9-MLLENL" = unique(c(inter_hoxa9_activated_ensID, inter_hoxa9_repressed_ensID)),
                       "Slany ChIPseq mouse" = genes_mm_Chip,
                       "Slany ChIPseq mouse promoter" = unique(mm_Chip$`Nearest Ensembl`[grep(pattern = "promoter", mm_Chip$Annotation)]),
                       "Slany ChIPseq human" = unique(ortho_genes$mmusculus_homolog_ensembl_gene),
                       "Slany ChIPseq human promoter" = unique(ortho_genes_prom$mmusculus_homolog_ensembl_gene)
                       )

# create the list of genes for GSEA
id_sort = order(fits$GFP.TPO$table$logFC, decreasing = T)
table_GSEA = fits$GFP.TPO$table[id_sort,]
geneList_GSEA = table_GSEA$logFC
names(geneList_GSEA) = rownames(table_GSEA)


# perform the GSEA analysis
fgseaRes <- fgsea(pathways = Slany_pathways, stats = geneList_GSEA)

# plot the enrichment
for(enrich in names(Slany_pathways)){
  
  # define the filename
  filename <- gsub(pattern = " ", replacement = "_", enrich)
  filename = paste0("GSEA_", filename, "_GFP_TPOvsnoTPO")
  
  # save the enrichment in the pdf
  saveFigures(fileName = filename, ggplot = plotEnrichment(Slany_pathways[[enrich]], geneList_GSEA) + 
                labs(title = enrich), dirPlot = plotDir, A4 = T)
  
}

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


#### perform ORA for Up_Genes for GFP+ TPO vs PBS (source: https://doi.org/10.1182/bloodadvances.2018025866) -----------------------------------------------------------

# select all tested genes
genes <- rownames(as.data.frame(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO))))

# select the significant toptags
dge_DEG <- as.data.frame(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO), p.value = 0.05))
dge_DEG_up <- rownames(dge_DEG)[dge_DEG$logFC > 0]
dge_DEG_down <- rownames(dge_DEG)[dge_DEG$logFC < 0]


# create the output directory
tmp_dir <- create_dir(path = "DEA/Result/Result_woSamples/GFP.TPO/Up_DEG_woUniverse")

# perform the enrichment analysis
EnrichmentTerm_Analysis(common = dge_DEG_up, outputDir = tmp_dir, 
                            mart = mart)

# create the output directory
tmp_dir <- create_dir(path = "DEA/Result/Result_woSamples/GFP.TPO/Up_DEG")

# perform the enrichment analysis
EnrichmentTerm_Analysis(common = dge_DEG_up, outputDir = tmp_dir, 
                        mart = mart, universe = genes)


# create the output directory
tmp_dir <- create_dir(path = "DEA/Result/Result_woSamples/GFP.TPO/Down_DEG_woUniverse")

# perform the enrichment analysis
EnrichmentTerm_Analysis(common = dge_DEG_down, outputDir = tmp_dir, 
                        mart = mart)

# create the output directory
tmp_dir <- create_dir(path = "DEA/Result/Result_woSamples/GFP.TPO/Down_DEG")

# perform the enrichment analysis
EnrichmentTerm_Analysis(common = dge_DEG_down, outputDir = tmp_dir, 
                        mart = mart, universe = genes)



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


#### perform ORA for Up_Genes for bulk TPO vs PBS (source: https://doi.org/10.1182/bloodadvances.2018025866) -----------------------------------------------------------

# select all tested genes
genes <- rownames(as.data.frame(topTags(fits$bulk.TPO, n = nrow(fits$bulk.TPO))))

# select the significant toptags
dge_DEG <- as.data.frame(topTags(fits$bulk.TPO, n = nrow(fits$bulk.TPO), p.value = 0.05))
dge_DEG_up <- rownames(dge_DEG)[dge_DEG$logFC > 0]
dge_DEG_down <- rownames(dge_DEG)[dge_DEG$logFC < 0]


# create the output directory
tmp_dir <- create_dir(path = "DEA/Result/Result_woSamples/bulk.TPO/Up_DEG_woUniverse")

# perform the enrichment analysis
EnrichmentTerm_Analysis(common = dge_DEG_up, outputDir = tmp_dir, 
                        mart = mart)

# create the output directory
tmp_dir <- create_dir(path = "DEA/Result/Result_woSamples/bulk.TPO/Up_DEG")

# perform the enrichment analysis
EnrichmentTerm_Analysis(common = dge_DEG_up, outputDir = tmp_dir, 
                        mart = mart, universe = genes)

# create the output directory
tmp_dir <- create_dir(path = "DEA/Result/Result_woSamples/bulk.TPO/Down_DEG_woUniverse")

# perform the enrichment analysis
EnrichmentTerm_Analysis(common = dge_DEG_down, outputDir = tmp_dir, 
                        mart = mart)


# create the output directory
tmp_dir <- create_dir(path = "DEA/Result/Result_woSamples/bulk.TPO/Down_DEG")

# perform the enrichment analysis
EnrichmentTerm_Analysis(common = dge_DEG_down, outputDir = tmp_dir, 
                        mart = mart, universe = genes)


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


# select results selected
GSEA_res_C2_selected <- subset(GSEA_res_C2, category != "")
GSEA_res_H_selected <- subset(GSEA_res_H, category != "")

# select the KEGG term
KEGG_term <- GSEA_res_C2_selected$pathway[grep(pattern = "KEGG", GSEA_res_C2_selected$pathway)]
# find in the KEGG database the id for the HEMATOPOIETIC_CELL_LINEAGE term

# display the pathway 
#pathview(gene.data = geneList, kegg.dir = "DEA/Result/Result_woSamples/GFP.TPO/", 
#         pathway.id = "04640", species = "ko", limit = list(gene=max(abs(geneList)), cpd=1))

test <- MS_convertGene(names(geneList), organism_code = "mmu", organism_name = "mouse", orthology = T)
geneList2 = geneList
names(geneList2) = test
setwd("DEA/Result/Result_woSamples/GFP.TPO/")
pathview(gene.data = geneList2, pathway.id = "04640", species = "ko", limit = list(gene=round(max(abs(geneList2))), cpd=1),
         gene.idtype = "KEGG")

pathview(gene.data = geneList2, pathway.id = "05200", species = "ko", limit = list(gene=round(max(abs(geneList2)))/2, cpd=1),
         gene.idtype = "KEGG")
pathview(gene.data = geneList2, pathway.id = "04144", species = "ko", limit = list(gene=round(max(abs(geneList2)))/2, cpd=1),
         gene.idtype = "KEGG")
pathview(gene.data = geneList2, pathway.id = "04072", species = "ko", limit = list(gene=round(max(abs(geneList2)))/2, cpd=1),
         gene.idtype = "KEGG")
pathview(gene.data = geneList2, pathway.id = "04218", species = "ko", limit = list(gene=round(max(abs(geneList2)))/2, cpd=1),
         gene.idtype = "KEGG")

setwd("../../../../")

## draw the heatmap

# read the dge file
dge <- readRDS("DEA/Rds/RDS_woD/EnsemblGenes.dge_samplesList.filtered.rds")

# create the matrix for the heatmap
mat_norm_count <- cpm(dge, normalized.lib.sizes = T)

# set the information about the matrix
annot_col <- data.frame("group" = dge$samples[colnames(mat_norm_count), "SampleGroup"])
rownames(annot_col) <- rownames(dge$samples[colnames(mat_norm_count),])


# function to draw the heatmaps
drawHeatmap <- function(dge, genes, mat_norm_count, dirPlot, filename){
  
  
  # reduce the matrix to the gene of interest
  tmp_mat <- mat_norm_count[genes,]
  rownames(tmp_mat) = dge$genes[rownames(tmp_mat),]$SYMBOL
  #colnames(tmp_mat) <- dge$samples[colnames(tmp_mat),]$ExternalSampleName
  
  
  # create the matrix with pheatmap
  heatmap <- pheatmap(tmp_mat, scale = "row", annotation_col = annot_col, fontsize_row = 4) #, annotation_row = )
  saveFigures(ggplot = heatmap, fileName = filename,
              dirPlot = dirPlot, A4 = T)
  
  
}


 
 # draw all the heatmap in the loop for C2 database
 for(pathway in GSEA_res_C2_selected$pathway){
   
   print(paste0(date(),", heatmap for ", pathway, " in progress..."))
   
   # select genes of interest
   genes <- rownames(dge$genes)[dge$genes$ENTREZID %in% Mm.c2[[pathway]]]
   
   # draw the heatmap
   drawHeatmap(dge = dge, genes = genes, mat_norm_count = mat_norm_count,
               dirPlot = "DEA/Result/Result_woSamples/GFP.TPO/", filename = paste0("heatmap_", pathway))
   print(paste0(date(),", heatmap for ", pathway, " done."))
   
 }
 
 
 # draw all the heatmap in the loop for H database
 for(pathway in GSEA_res_H_selected$pathway){
   
   print(paste0(date(),", heatmap for ", pathway, " in progress..."))
   
   # select genes of interest
   genes <- rownames(dge$genes)[dge$genes$ENTREZID %in% Mm.H[[pathway]]]
   
   # draw the heatmap
   drawHeatmap(dge = dge, genes = genes, mat_norm_count = mat_norm_count,
               dirPlot = "DEA/Result/Result_woSamples/GFP.TPO/", filename = paste0("heatmap_", pathway))
   print(paste0(date(),", heatmap for ", pathway, " done."))
   
 }
 
 
 
 # try to find common gene in Jürg's selection
 GSEA_res_C2_Jurg <- read.xlsx2(file = "DEA/Result/Result_woSamples/GFP.TPO/HE-190822-GSEA_Msigdatabase_GFP_TPOvsnoTPO.xlsx", 
                           header = T, sheetName = "Mm.c2")
 GSEA_res_H_Jurg <- read.xlsx2(file = "DEA/Result/Result_woSamples/GFP.TPO/HE-190822-GSEA_Msigdatabase_GFP_TPOvsnoTPO.xlsx", 
                                header = T, sheetName = "Mm.H")
 
 # check if the pathway contains Hoxa9 (entrez Id 15405)
 GSEA_res_C2_Jurg$'contains Hoxa9' = sapply(GSEA_res_C2_Jurg$pathway, function(x)any(Mm.c2[[x]] == "15405"))
 GSEA_res_C2_Jurg$'contains Hoxa11' = sapply(GSEA_res_C2_Jurg$pathway, function(x)any(Mm.c2[[x]] == "15396"))
 GSEA_res_H_Jurg$'contains Hoxa9' = sapply(GSEA_res_H$pathway, function(x)any(Mm.H[[x]] == "15405"))
 GSEA_res_H_Jurg$'contains Hoxa11' = sapply(GSEA_res_H$pathway, function(x)any(Mm.H[[x]] == "15396"))
 
 
 # select HOXA9 terms
 HOXA9_terms <- GSEA_res_C2_Jurg$pathway[GSEA_res_C2_Jurg$category == "Hoxa9_J"]
 
 # write the table in the excel file
 write.xlsx2(GSEA_res_C2_Jurg, file = "DEA/Result/Result_woSamples/GFP.TPO/HE-190822-GSEA_Msigdatabase_GFP_TPOvsnoTPO_Hoxa9.xlsx", 
             col.names = T, row.names = F, sheetName = "Mm.c2")
 gc()
 
 # write the subtable in the excel file
 write.xlsx2(GSEA_res_C2_Jurg[GSEA_res_C2_Jurg$`contains Hoxa9`,], file = "DEA/Result/Result_woSamples/GFP.TPO/HE-190822-GSEA_Msigdatabase_GFP_TPOvsnoTPO_Hoxa9.xlsx", 
             col.names = T, row.names = F, append = T, sheetName = "Mm.c2; Hoxa9")
 gc()
 
 # select list of genes
 HOXA9_terms_genes <- Mm.c2[HOXA9_terms]
 
 # select common genes
 HESS_terms_commonGenes <- unique(Reduce(c, HOXA9_terms_genes[grep(pattern = "HESS", names(HOXA9_terms_genes))]))
 TAKEDA_terms_commonGenes <- unique(Reduce(c, HOXA9_terms_genes[grep(pattern = "TAKEDA", names(HOXA9_terms_genes))]))
 KUMAR_terms_commonGenes <- unique(Reduce(c, HOXA9_terms_genes[grep(pattern = "KUMAR", names(HOXA9_terms_genes))]))
 VERHAAK_terms_commonGenes <- unique(Reduce(c, HOXA9_terms_genes[grep(pattern = "VERHAAK", names(HOXA9_terms_genes))]))
 common_genes <- intersect(intersect(intersect(HESS_terms_commonGenes, TAKEDA_terms_commonGenes), KUMAR_terms_commonGenes), VERHAAK_terms_commonGenes)
 
 # function to draw the Venn diagram
 draw_VennDiagram <- function(x, category, filename, color){
   venn.diagram(
     x = x,
     category.names = category,
     filename = filename,
     output=TRUE,
     
     # Output features
     imagetype="png" ,
     height = 480 , 
     width = 480 , 
     resolution = 300,
     compression = "lzw",
     
     # Circles
     lwd = 2,
     lty = 'blank',
     fill = color,
     
     # Numbers
     cex = .6,
     fontface = "bold",
     fontfamily = "sans",
     
     # Set names
     cat.cex = 0.3,
     cat.fontface = "bold",
     cat.default.pos = "outer")
 }
 
 
 # create the corresponding Venn Diagram
 myCol <- brewer.pal(4, "Pastel2")
 draw_VennDiagram(x = list(HESS_terms_commonGenes, TAKEDA_terms_commonGenes, KUMAR_terms_commonGenes, VERHAAK_terms_commonGenes), 
                  category = c("HESS" , "TAKEDA" , "KUMAR", "VERHAAK"), 
                  filename = 'DEA/Result/Result_woSamples/GFP.TPO/HOXA9_venn_diagramm.png',
                  color = myCol)
 
 
 # provide list of gene for TAKEDA, HESS and VERHAAK comparison
 genes_24 <- intersect(intersect(TAKEDA_terms_commonGenes, HESS_terms_commonGenes), VERHAAK_terms_commonGenes)
 
 # convert entrezId in ensembl_genes
 genes_24 <- getBM(filters = "entrezgene_id",
                attributes = c("ensembl_gene_id","entrezgene_id", "external_gene_name"),
                values = genes_24,
                mart = mart)
 
 # save the table in a excel file
 write.xlsx2(genes_24, file = "DEA/Result/Result_woSamples/GFP.TPO/HOXA9_24genes_Takeda_Hess_Verhaak.xlsx", 
             col.names = T, row.names = F)
 gc()
 
 
 # select interferon terms
 interferon_terms <- GSEA_res_C2_Jurg$pathway[GSEA_res_C2_Jurg$category == "interferon"]
 
 # select list of genes
 interferon_terms_genes <- Mm.c2[interferon_terms]
 
 # select common genes
 REACTOME_terms_commonGenes <- unique(Reduce(c, interferon_terms_genes[grep(pattern = "REACTOME", names(interferon_terms_genes))]))
 BOSCO_terms_commonGenes <- unique(Reduce(c, interferon_terms_genes[grep(pattern = "BOSCO", names(interferon_terms_genes))]))
 BROWNE_terms_commonGenes <- unique(Reduce(c, interferon_terms_genes[grep(pattern = "BROWNE", names(interferon_terms_genes))]))
 ZHANG_terms_commonGenes <- unique(Reduce(c, interferon_terms_genes[grep(pattern = "ZHANG", names(interferon_terms_genes))]))
 common_genes_interferon <- intersect(intersect(intersect(REACTOME_terms_commonGenes, BOSCO_terms_commonGenes), BROWNE_terms_commonGenes), ZHANG_terms_commonGenes)
 
 # draw the Venn diagram
 draw_VennDiagram(x = list(REACTOME_terms_commonGenes, BOSCO_terms_commonGenes, BROWNE_terms_commonGenes, ZHANG_terms_commonGenes), 
                  category = c("REACTOME" , "BOSCO" , "BROWNE", "ZHANG"), 
                  filename = 'DEA/Result/Result_woSamples/GFP.TPO/interferon_venn_diagramm.png',
                  color = myCol)
 
 
 ##### comparison with Slany -----------------------------------------------------------
 
 # get all genes in mouse
 mouse_genes <- getMouseGene()
 
 # hoxa list genes
 hoxa9_activated <- merge(x = hoxa9_activated, y = mouse_genes, all.x = T, by.x = "gene_short_name", by.y = "external_gene_name")
 hoxa9_repressed <- merge(x = hoxa9_repressed, y = mouse_genes, all.x = T, by.x = "gene_short_name", by.y = "external_gene_name")
 
 # draw the Venn diagram for Takeda
 draw_VennDiagram(x = list(TAKEDA_terms_commonGenes, na.omit(hoxa9_activated$entrezgene_id), na.omit(hoxa9_repressed$entrezgene_id)), 
                  category = c("TAKEDA" , "Slany - Hoxa9 activated" , "Slany - Hoxa9 repressed"), 
                  filename = 'DEA/Result/Result_woSamples/GFP.TPO/Takeda_Slany_venn_diagramm.png',
                  color = myCol[1:3])
 
 # selected terms
 HESS_terms_up <- HOXA9_terms_genes[["HESS_TARGETS_OF_HOXA9_AND_MEIS1_UP"]]
 
 # draw the Venn diagram for Takeda
 draw_VennDiagram(x = list(HESS_terms_up, na.omit(hoxa9_activated$entrezgene_id), na.omit(hoxa9_repressed$entrezgene_id)), 
                  category = c("HESS" , "Slany - Hoxa9 activated" , "Slany - Hoxa9 repressed"), 
                  filename = 'DEA/Result/Result_woSamples/GFP.TPO/Hess_Slany_venn_diagramm.png',
                  color = myCol[1:3])
 
 
 # selected terms
 KUMAR_MLLAF9_terms <- HOXA9_terms_genes[["KUMAR_TARGETS_OF_MLL_AF9_FUSION"]]
 
 # draw the Venn diagram for Takeda
 draw_VennDiagram(x = list(KUMAR_MLLAF9_terms, na.omit(hoxa9_activated$entrezgene_id), na.omit(hoxa9_repressed$entrezgene_id)), 
                  category = c("KUMAR" , "Slany - Hoxa9 activated" , "Slany - Hoxa9 repressed"), 
                  filename = 'DEA/Result/Result_woSamples/GFP.TPO/Kumar_Slany_venn_diagramm.png',
                  color = myCol[1:3])
 
 # do the count based on intersection
 inter_hess_activated <- intersect(HESS_terms_up, na.omit(hoxa9_activated$entrezgene_id))
 inter_hess_repressed <- intersect(HESS_terms_up, na.omit(hoxa9_repressed$entrezgene_id))
 
 
 
 #### bulk TPO vs noTPO -----------------------------------------------------------
 
 
 # select gene list for bulk TPO vs noTPO
 geneList_bulk <- getGeneListforGSEA(fits$bulk.TPO)

 # upload the GSEA result with Hugues-Etienne selection
 GSEA_res_C2_bulk <- read.xlsx2(file = "DEA/Result/Result_woSamples/bulk.TPO/GSEA_Msigdatabase_bulk_TPOvsnoTPO_HE.xlsx", 
                           header = T, sheetName = "Mm.c2")
 GSEA_res_H_bulk <- read.xlsx2(file = "DEA/Result/Result_woSamples/bulk.TPO/GSEA_Msigdatabase_bulk_TPOvsnoTPO_HE.xlsx",
                          header = T, sheetName = "Mm.H")
 
 # select results selected
 GSEA_res_C2_bulk_selected <- subset(GSEA_res_C2_bulk, category != "")
 GSEA_res_H_bulk_selected <- subset(GSEA_res_H_bulk, category != "")
 
 
 # draw all the heatmap in the loop for C2 database
 for(pathway in GSEA_res_C2_bulk_selected$pathway){
   
   print(paste0(date(),", heatmap for ", pathway, " in progress..."))
   
   # select genes of interest
   genes <- rownames(dge$genes)[dge$genes$ENTREZID %in% Mm.c2[[pathway]]]
   
   # draw the heatmap
   drawHeatmap(dge = dge, genes = genes, mat_norm_count = mat_norm_count,
               dirPlot = "DEA/Result/Result_woSamples/bulk.TPO/", filename = paste0("heatmap_", pathway))
   print(paste0(date(),", heatmap for ", pathway, " done."))
   
 }
 
 
 # draw all the heatmap in the loop for H database
 for(pathway in GSEA_res_H_bulk_selected$pathway){
   
   print(paste0(date(),", heatmap for ", pathway, " in progress..."))
   
   # select genes of interest
   genes <- rownames(dge$genes)[dge$genes$ENTREZID %in% Mm.H[[pathway]]]
   
   # draw the heatmap
   drawHeatmap(dge = dge, genes = genes, mat_norm_count = mat_norm_count,
               dirPlot = "DEA/Result/Result_woSamples/bulk.TPO/", filename = paste0("heatmap_", pathway))
   print(paste0(date(),", heatmap for ", pathway, " done."))
   
 }
 
 
 # try to find common gene in Jürg's selection
 GSEA_res_C2_bulk_Jurg <- read.xlsx2(file = "DEA/Result/Result_woSamples/bulk.TPO/HE-190822-GSEA_Msigdatabase_bulk_TPOvsnoTPO.xlsx", 
                                header = T, sheetName = "Mm.c2")
 
 # select HOXA9 terms
 HOXA9_bulk_terms <- GSEA_res_C2_bulk_Jurg$pathway[GSEA_res_C2_bulk_Jurg$category == "Hoxa9_J"]
 
 # select list of genes
 HOXA9_bulk_terms_genes <- Mm.c2[HOXA9_bulk_terms]
 
 # select common genes
 HESS_bulk_terms_commonGenes <- unique(Reduce(c, HOXA9_bulk_terms_genes[grep(pattern = "HESS", names(HOXA9_bulk_terms_genes))]))
 TAKEDA_bulk_terms_commonGenes <- unique(Reduce(c, HOXA9_bulk_terms_genes[grep(pattern = "TAKEDA", names(HOXA9_bulk_terms_genes))]))
 KUMAR_bulk_terms_commonGenes <- unique(Reduce(c, HOXA9_bulk_terms_genes[grep(pattern = "KUMAR", names(HOXA9_bulk_terms_genes))]))
 VERHAAK_bulk_terms_commonGenes <- unique(Reduce(c, HOXA9_bulk_terms_genes[grep(pattern = "VERHAAK", names(HOXA9_bulk_terms_genes))]))
 common_bulk_genes <- intersect(intersect(intersect(HESS_bulk_terms_commonGenes, TAKEDA_bulk_terms_commonGenes), KUMAR_bulk_terms_commonGenes), VERHAAK_bulk_terms_commonGenes)
 
 # create the corresponding Venn Diagram
 draw_VennDiagram(x = list(HESS_bulk_terms_commonGenes, TAKEDA_bulk_terms_commonGenes, KUMAR_bulk_terms_commonGenes, VERHAAK_bulk_terms_commonGenes), 
                  category = c("HESS" , "TAKEDA" , "KUMAR", "VERHAAK"), 
                  filename = 'DEA/Result/Result_woSamples/bulk.TPO/HOXA9_venn_diagramm.png',
                  color = myCol)
 
 
 # select interferon terms
 interferon_bulk_terms <- GSEA_res_C2_bulk_Jurg$pathway[GSEA_res_C2_bulk_Jurg$category == "interferon"]
 
 # select list of genes
 interferon_bulk_terms_genes <- Mm.c2[interferon_bulk_terms]
 
 # select common genes
 REACTOME_bulk_terms_commonGenes <- unique(Reduce(c, interferon_bulk_terms_genes[grep(pattern = "REACTOME", names(interferon_bulk_terms_genes))]))
 BROWNE_bulk_terms_commonGenes <- unique(Reduce(c, interferon_bulk_terms_genes[grep(pattern = "BROWNE", names(interferon_bulk_terms_genes))]))
 DER_bulk_terms_commonGenes <- unique(Reduce(c, interferon_bulk_terms_genes[grep(pattern = "DER", names(interferon_bulk_terms_genes))]))
 HECKER_terms_commonGenes <- unique(Reduce(c, interferon_bulk_terms_genes[grep(pattern = "HECKER", names(interferon_bulk_terms_genes))]))
 common_genes_bulk_interferon <- intersect(intersect(intersect(REACTOME_bulk_terms_commonGenes, BROWNE_bulk_terms_commonGenes), DER_bulk_terms_commonGenes), HECKER_terms_commonGenes)
 
 # draw the Venn diagram
 draw_VennDiagram(x = list(REACTOME_bulk_terms_commonGenes, BROWNE_bulk_terms_commonGenes, DER_bulk_terms_commonGenes, HECKER_terms_commonGenes), 
                  category = c("REACTOME" , "BROWNE" , "DER", "HECKER"), 
                  filename = 'DEA/Result/Result_woSamples/bulk.TPO/interferon_venn_diagramm.png',
                  color = myCol)
 
 
 
 # heatmap for top 100 gene expression from Slany's paper -----------------------------------------------------------
 
 # load the RNAseq result with the expression value from Slany's paper
Slany_RNAseq <- read_excel(path = "External_Data/Slany_paper/ba025866-suppl4.xlsx", sheet = "Expressed genes")
 
# add EnsemblId for genes
Slany_RNAseq = merge(x= mouse_genes, y = Slany_RNAseq, 
                                by.x = "external_gene_name",
                                by.y = "gene_short_name")
 
# select genes already present in the matrix
Slany_RNAseq = Slany_RNAseq[Slany_RNAseq$ensembl_gene_id %in% rownames(mat_norm_count), ]
 
 # reorder the Slany_RNAseq object
 Slany_RNAseq = Slany_RNAseq %>% arrange(`sum logfold 24-72`) 

 # select top 100 genes hoxa9_activated
 hoxa9_activated_top100 <- head(Slany_RNAseq, n=100)
 
 # select top 100 genes hoxa9_repressed
 hoxa9_repressed_top100 <- tail(Slany_RNAseq, n=100)
 
 

 # perform the heatmaps
 drawHeatmap(dge = dge, genes = hoxa9_activated_top100$ensembl_gene_id, mat_norm_count = mat_norm_count,
             dirPlot = "DEA/Result/Result_woSamples/GFP.TPO/", filename = paste0("heatmap_top100_Hoxa9_Activated"))
 drawHeatmap(dge = dge, genes = hoxa9_repressed_top100$ensembl_gene_id, mat_norm_count = mat_norm_count,
             dirPlot = "DEA/Result/Result_woSamples/GFP.TPO/", filename = paste0("heatmap_top100_Hoxa9_Repressed"))
 

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




## create matplot based on the TPO pathway for DEG in GFP+ and bulk ------------------------------------------------

# find orthologous genes of THPO pathway
THPO_mouse_genes = ortho[ortho$external_gene_name %in% THPO_pathway$TPO_pathway, ]

# create the dataframe to store the values
THPO_mouse_genes = THPO_mouse_genes[, c("mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name")]
# add the GFP+ results
tmp_table <- as.data.frame(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO)))
THPO_mouse_genes = merge(x = THPO_mouse_genes, y = tmp_table, 
                         by.x = "mmusculus_homolog_associated_gene_name",
                         by.y = "SYMBOL", suffixes = c("", "_GFP+"))
# add the bulk results
tmp_table <- as.data.frame(topTags(fits$bulk.TPO, n = nrow(fits$bulk.TPO)))
THPO_mouse_genes = merge(x = THPO_mouse_genes, y = tmp_table, 
                         by.x = "mmusculus_homolog_associated_gene_name",
                         by.y = "SYMBOL", suffixes = c( "_GFP+", "_bulk"))

# remove the tmp_table
rm(tmp_table)
gc()

# add the rownames for the THPO_mouse_gene
rownames(THPO_mouse_genes) = THPO_mouse_genes$mmusculus_homolog_associated_gene_name

# create the dataframe for the heatmap
dt_THPO_plot <- data.frame("Cells" = rep(c("GFP+", "bulk"), nrow(THPO_mouse_genes)),
                           "Genes" = rep(THPO_mouse_genes$mmusculus_homolog_associated_gene_name, each = 2),
                           "FC" = NA)

# add information about the FC
id_FC <- THPO_mouse_genes$mmusculus_homolog_associated_gene_name[THPO_mouse_genes$`FDR_GFP+` <= 0.05]
dt_THPO_plot$FC[(dt_THPO_plot$Cells == "GFP+") & (dt_THPO_plot$Genes %in% id_FC)] = THPO_mouse_genes[id_FC, "logFC_GFP+"]
id_FC <- THPO_mouse_genes$mmusculus_homolog_associated_gene_name[THPO_mouse_genes$FDR_bulk <= 0.05]
dt_THPO_plot$FC[(dt_THPO_plot$Cells == "bulk") & (dt_THPO_plot$Genes %in% id_FC)] = THPO_mouse_genes[id_FC, "logFC_bulk"]


# create the parameters for the matplot
param_plot <-  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 60, vjust = 0.5, face = "bold", color = "black"),
        axis.title = element_lank(), plot.background = element_rect(fill = "white"),
        axis.text.y = element_text(size = 12, face="bold", color = "black"),
        legend.text = element_text(size = 10, face = "bold", color = "black"),
        legend.title = element_text(size = 12, face = "bold", color = "black")) 

# create the first heatmap
matplot1 <- ggplot(dt_THPO_plot[1:(43*2), ], aes(x = Genes , y = Cells, fill = FC)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "white") + param_plot
  

# create the second heatmap
matplot2 <- ggplot(dt_THPO_plot[((43*2)+1):nrow(dt_THPO_plot), ], aes(x = Genes , y = Cells, fill = FC)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "white") + param_plot

# combine the 2 matplot
ggsave(filename="DEA/Result/Result_woSamples/common_DEG_GFP+_bulk/THPO_matplot.pdf", 
       plot = grid.arrange(matplot1, matplot2, ncol=1, nrow= 2), 
       device = cairo_pdf, 
       width = 297, 
       height = 210, 
       units = "mm")
saveFigures(fileName = "THPO_matplot2", ggplot = grid.arrange(matplot1, matplot2, ncol=1, nrow= 2),
            dirPlot = "DEA/Result/Result_woSamples/common_DEG_GFP+_bulk/", A4 = T)

 
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


# select commons genes between HESS and TAKEDA
intersect_commonGenes <- intersect(HESS_bulk_terms_commonGenes, 
                             TAKEDA_terms_commonGenes)



# see to do a heatmap or dotplot with specific terms 
DEG_GFP <- as.data.frame(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO),
                                 p.value = 0.05))

# change the rownames of the table
DEG_GFP$ENTREZID = as.character(DEG_GFP$ENTREZID)

# select DEG in common with HESS and TAKEDA
final_genes_heatplot <- intersect(DEG_GFP$ENTREZID, intersect_commonGenes)

# select the pathways
pathways_heatmap <-  c(names(Mm.c2)[grep(pattern = "TAKEDA", names(Mm.c2))],
                       names(Mm.c2)[grep(pattern = "^HESS_", names(Mm.c2))])

# select the count of genes
mat_heatmap <- matrix(data = NA, nrow = length(pathways_heatmap),
                      ncol = length(final_genes_heatplot))
colnames(mat_heatmap) = final_genes_heatplot
rownames(mat_heatmap) = pathways_heatmap  


# update the table in loop
for(pathway in pathways_heatmap){
  
  # select the genes in the pathways
  tmp_genes <- Mm.c2[[pathway]]
  
  # select the subtable for the DEG_GFP
  tmp_DEG_GFP = DEG_GFP[DEG_GFP$ENTREZID %in% tmp_genes, ]
  
  # change the rownames
  rownames(tmp_DEG_GFP) = tmp_DEG_GFP$ENTREZID
  
  # select the genes
  inter_genes <- intersect(dt_heatmap[dt_heatmap$Pathway %in% pathway, "Genes"], tmp_DEG_GFP$ENTREZID)
  
  # update the table
  print(paste0("pathway: ", pathway, "; inter_gene: ", length(inter_genes)))
  if(length(inter_genes) > 0){
    
    mat_heatmap[pathway, inter_genes] = tmp_DEG_GFP[inter_genes,"logFC"]

  }
  
  
}

# change the colnames of the matrix
tmp_DEG_GFP = DEG_GFP[DEG_GFP$ENTREZID %in% colnames(mat_heatmap), ]
rownames(tmp_DEG_GFP) = tmp_DEG_GFP$ENTREZID
colnames(mat_heatmap) = tmp_DEG_GFP[colnames(mat_heatmap),"SYMBOL"]

# remove pathways without values
mat_heatmap = mat_heatmap[!apply(mat_heatmap, 1, function(x) all(is.na(x))), ]

# create the data.frame based on matrix
dt_heatmap <-  as.data.frame(mat_heatmap)
dt_heatmap <- stack(dt_heatmap)
dt_heatmap <- cbind(dt_heatmap, rownames(mat_heatmap))

# change the colnames of the dataframe
colnames(dt_heatmap) = c("values", "gene", "pathway")

# create heatmap with ggplot2
matplot <- ggplot(dt_heatmap, aes(x = gene , y = pathway, fill = values)) +
  geom_tile(colour = "white") + theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 60, vjust = 0.5),
        axis.title = element_blank(), plot.background = element_rect(fill = "white")) +
  scale_fill_gradient(low = "blue", high = "red", na.value = "white")


# save the figure in the files
saveFigures(fileName = "matplot_GSEA_paper_HOXA9", ggplot = matplot, 
            dirPlot = "DEA/Result/Result_woSamples/GFP.TPO/", A4 = T)


# create the heatmap for the common genes between Valk and Target ------------------------------------------------

# read the table of common genes
table_human_gene_common_Valk_Target <- read.xlsx2(file = "DEA/Human_patients/TARGET_allMLL/common_patients_Valk_TARGET_ortho.xlsx", 
                                                  sheetIndex = 1)

# draw the heatmap
drawHeatmap(dge = dge, genes = intersect(rownames(mat_norm_count), table_human_gene_common_Valk_Target$mmusculus_homolog_ensembl_gene), 
            mat_norm_count = mat_norm_count,
            dirPlot = "DEA/Human_patients/TARGET_allMLL/", filename = paste0("heatmap_common_patients_Valk_TARGET_ortho"))


# create the heatmap for the common genes between Valk, Target and BEAT ------------------------------------------------

# read the table of common genes
table_human_gene_common_Valk_Target_BEAT <- read.xlsx2(file = "DEA/Human_patients/BEAT_MLLAF9/common_patients_Valk_TARGET_BEAT_ortho.xlsx", 
                                                  sheetIndex = 1)

# draw the heatmap
drawHeatmap(dge = dge, genes = intersect(rownames(mat_norm_count), table_human_gene_common_Valk_Target_BEAT$mmusculus_homolog_ensembl_gene), 
            mat_norm_count = mat_norm_count,
            dirPlot = "DEA/Human_patients/BEAT_MLLAF9/", filename = "heatmap_common_patients_Valk_TARGET_BEAT_ortho")




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



# compare the DEG for Bulk and GFP+ TPO vs noTPO comparison with Slany ------------------------------------------------


# create the corresponding Venn Diagram with Slany
myCol <- brewer.pal(3, "Pastel2")
draw_VennDiagram(x = list(common_DEG_bulk_GFP, hoxa9_activated_ensID, hoxa9_repressed_ensID), 
                 category = c("DEG bulk\n & GFP+" , "Hoxa9\n activated" , "Hoxa9\n repressed"), 
                 filename = 'DEA/Result/Result_woSamples/common_DEG_GFP+_bulk/vennDiagram.png',
                 color = myCol)

# create the corresponding Venn Diagram with Slany (down)
myCol <- brewer.pal(3, "Pastel2")
draw_VennDiagram(x = list(common_DEG_bulk_GFP_down, hoxa9_activated_ensID, hoxa9_repressed_ensID), 
                 category = c("DEG bulk\n & GFP+" , "Hoxa9\n activated" , "Hoxa9\n repressed"), 
                 filename = 'DEA/Result/Result_woSamples/common_DEG_GFP+_bulk_down/vennDiagram_down.png',
                 color = myCol)

# create the corresponding Venn Diagram with Slany (up)
myCol <- brewer.pal(3, "Pastel2")
draw_VennDiagram(x = list(common_DEG_bulk_GFP_up, hoxa9_activated_ensID, hoxa9_repressed_ensID), 
                 category = c("DEG bulk\n & GFP+" , "Hoxa9\n activated" , "Hoxa9\n repressed"), 
                 filename = 'DEA/Result/Result_woSamples/common_DEG_$GFP+_bulk_up/vennDiagram_up.png',
                 color = myCol)


## create the GSEA analysis for the comparison with Slany ------------------------------------------------

# select all DEGs in the GFP+
table_GFPp <- as.data.frame(topTags(fits$GFP.TPO, n = nrow(fits$GFP.TPO)))

# select all DEGs in the bulk
table_bulk <- as.data.frame(topTags(fits$bulk.TPO, n = nrow(fits$bulk.TPO)))

# update with the genes
table_GFPp$ensembl_ID = rownames(table_GFPp)
table_bulk$ensembl_ID = rownames(table_bulk)

# prepare the gene list 
id_sort_GFPp <- order(table_GFPp$logFC, decreasing = T)
GFPp <- table_GFPp$logFC[id_sort_GFPp]
names(GFPp) <- table_GFPp$ensembl_ID[id_sort_GFPp]

id_sort_bulk <- order(table_bulk$logFC, decreasing = T)
bulk <- table_bulk$logFC[id_sort_bulk]
names(bulk) <- table_bulk$ensembl_ID[id_sort_bulk]

# perform the GSEA
fgseaRes <- fgsea(pathways = Slany_pathways, stats = GFPp)

# create the output directory for the GSEA analysis
plotDir_GSEA <- "DEA/Plots/GFP+_HE"
if(!dir.exists(plotDir_GSEA)) dir.create(plotDir_GSEA, recursive = T)

# plot the enrichment
for(enrich in names(Slany_pathways)){
  
  # define the filename
  filename <- gsub(pattern = " ", replacement = "_", enrich)
  filename = paste0("GSEA_", filename, "_GFP_TPOvsnoTPO")
  
  # save the enrichment in the pdf
  saveFigures(fileName = filename, ggplot = plotEnrichment(Slany_pathways[[enrich]], GFPp) + 
                labs(title = enrich), dirPlot = plotDir_GSEA, A4 = T)

}

# write the excel file
write_xlsx(fgseaRes, path = file.path(plotDir_GSEA, "GSEA_result_GFP+_HE.xlsx"))



# create the output directory for the GSEA analysis
plotDir_GSEA <- "DEA/Plots/bulk_HE"
if(!dir.exists(plotDir_GSEA)) dir.create(plotDir_GSEA, recursive = T)

# perform the GSEA
fgseaRes <- fgsea(pathways = Slany_pathways, stats = bulk)


# plot the enrichment
for(enrich in names(Slany_pathways)){
  
  # define the filename
  filename <- gsub(pattern = " ", replacement = "_", enrich)
  filename = paste0("GSEA_", filename, "_bulk_TPOvsnoTPO")
  
  # save the enrichment in the pdf
  saveFigures(fileName = filename, ggplot = plotEnrichment(Slany_pathways[[enrich]], GFPp) + 
                labs(title = enrich), dirPlot = plotDir_GSEA, A4 = T)
  
}

# write the excel file
write_xlsx(fgseaRes, path = file.path(plotDir_GSEA, "GSEA_result_bulk_HE.xlsx"))

