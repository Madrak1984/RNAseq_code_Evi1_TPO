##### command lines to create z-score heatmap for paper
# done by Jonathan Seguin, group of Prof. Schwaller, DBM, UKBB, Basel, Switzerland
# email: jonathan.seguin@unibas.ch, seguin.jonathan@gmail.com
# Tue Mar 12 11:24:35 2024

# prepare the working environment -----------------------------------------------------

## load the libraries -----------------------------------------------------
library(pheatmap)
library(tidyverse)
library(edgeR)

## create the virtual object -----------------------------------------------------

# load the fits
fits <- readRDS(file = "DEA/Rds/RDS_woD/edge_samplesR_fits_results_DEA.rds")


# load the dge
dge <- readRDS(file = "DEA/Rds/RDS_woD/EnsemblGenes.dge_samplesList.filtered.rds")

# get the log cpm of dge
logcpm <- cpm(dge, log = T)


# create the heatmap for GFP+ -----------------------------------------------------

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

# draw the zscore heatmap (source: https://bioinformatics.stackexchange.com/questions/18359/in-a-rna-seq-heatmap-should-you-do-z-score-standardisation-before-clustering-the)
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




# create the heatmap for bulk+ -----------------------------------------------------

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

# draw the zscore heatmap (source: https://bioinformatics.stackexchange.com/questions/18359/in-a-rna-seq-heatmap-should-you-do-z-score-standardisation-before-clustering-the)
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
