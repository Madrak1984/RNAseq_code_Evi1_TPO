###### command lines to analyze the RNA-seq from Vaïa's paper (PMID:27344946) &
# GSE65383 from the Gene expression dataset
# done by Jonathan Seguin, group of Prof. Schwaller, DBM, UKBB, Basel, Switzerland
# email: jonathan.seguin@unibas.ch, seguin.jonathan@gmail.com
# Tue Apr  5 10:58:17 2022



# load the libraries ------------------------------------------------------
library(edgeR)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(factoextra)
library(plotly)

# create the output directories -------------------------------------------

# create the plot directory
plotDir <- "DEA/Plots"
if(!dir.exists(plotDir)) dir.create(plotDir, recursive = T)

# create the RDS directory
RDSDir <-  "DEA/Rds"
if(!dir.exists(RDSDir)) dir.create(RDSDir)

# prepare the DGE object -------------------------------------------------

# load the count table
count_table <- read.table("Dataset/GSE65383/GSE65383_raw_read_counts_mm9_Refseq.txt", 
                          header = T, sep = "\t", row.names = "refseq_name")

# create the dataframe for the sample information
sample_name <- gsub(pattern = "AP_SK_MllAF9_BMTs_", replacement = "", 
                    colnames(count_table))
colnames(count_table) = sample_name
dt_group <- data.frame("sampleName" = sample_name, "group" = c(rep("WT-GMP", 4), 
            rep("LTearly", 4), rep("LTlate", 4), rep("AML-GMP", 4))) 

# add gene information
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
annot <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
              "refseq_mrna"), filters = "refseq_mrna", 
              values = rownames(count_table), mart = ensembl)
annot <- merge(x=rownames(count_table), y=annot)


# create the DGEList object
dgevaia <- DGEList(counts = count_table, group = as.factor(dt_group$group))



# create the PCA to compare with the PCA in the paper ---------------------


# remove genes without more than 10 reads in all samples
dgevaia_tmp = dgevaia[which(rowSums(dgevaia$counts) > 10),]

# normalization and estimate dispersion
df.dgevaia_tmp <- calcNormFactors(dgevaia_tmp)
moma_vaia <- model.matrix(~ 0 + group, data=dgevaia_tmp$samples) # ~ 0 + SampleGroup
colnames(moma_vaia) <- sub("SampleGroup", "", colnames(moma_vaia))
dgevaia_tmp <- calcNormFactors(dgevaia_tmp)
dgevaia_tmp <- estimateDisp(dgevaia_tmp, design=moma_vaia)


# plot a PCA with Deseq2
dds_vaia = DESeqDataSetFromMatrix(countData = round(dgevaia_tmp$counts),
                             colData = dgevaia_tmp$samples, design = ~group)
transform.dds_vaia <- rlog(dds_vaia, blind = T)

# create the PCA data
pcaData <- DESeq2::plotPCA(transform.dds_vaia, intgroup = c("group"), returnData =T, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar")) # must be before to keep the percentVar information
pcaData <- cbind(pcaData, "SampleName" = SampleName)

# do the PCA plot
pcaplot <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  #coord_fixed() +
  #scale_color_manual(values = rainbow(6))
  scale_fill_brewer(palette = "Set1") +
  geom_text_repel(data = pcaData, color="black", aes(label=name))



# problem with the sample 15, redo without this sample --------------------

# remove the sample 15
dgevaia = dgevaia[,-15]

# remove genes without more than 10 reads in all samples
dgevaia_tmp = dgevaia[which(rowSums(dgevaia$counts) > 10),]

# normalization and estimate dispersion
df.dgevaia_tmp <- calcNormFactors(dgevaia_tmp)
moma_vaia <- model.matrix(~ 0 + group, data=dgevaia_tmp$samples) # ~ 0 + SampleGroup
colnames(moma_vaia) <- sub("SampleGroup", "", colnames(moma_vaia))
dgevaia_tmp <- calcNormFactors(dgevaia_tmp)
dgevaia_tmp <- estimateDisp(dgevaia_tmp, design=moma_vaia)


# plot a PCA with Deseq2
dds_vaia = DESeqDataSetFromMatrix(countData = round(dgevaia_tmp$counts),
                                  colData = dgevaia_tmp$samples, design = ~group)
transform.dds_vaia <- rlog(dds_vaia, blind = T)

# create the PCA data
pcaData <- DESeq2::plotPCA(transform.dds_vaia, intgroup = c("group"), returnData =T, ntop = nrow(dds_vaia))
percentVar <- round(100 * attr(pcaData, "percentVar")) # must be before to keep the percentVar information
pcaData <- cbind(pcaData, "SampleName" = name)

# do the PCA plot
pcaplot <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  #coord_fixed() +
  #scale_color_manual(values = rainbow(6))
  scale_fill_brewer(palette = "Set1") +
  geom_text_repel(data = pcaData, color="black", aes(label=name))


# save the dge object in one rds file
saveRDS(dgevaia, "DEA/Rds/dge_vaia_wo15.rds")


# combine the PCAs --------------------------------------------------------

# read the rds file from the TPO-Evi1 dataset
dge <- readRDS("Snakemake/ensdb_102_dge_list.rds")
se <-  readRDS("Snakemake/ensdb_102_summarized_experiment.rds")

# remove the genes without read 
dge = dge[which(rowSums(dge$counts) > 10),]

# reduce the count table
count_table_red <- count_table[annot$refseq_mrna,]

# select the count gene
count_table_dge <-  as.data.frame(cpm(dge))
count_table_dge$gene = rownames(dge)

# merge the count_table
count_table_dge <- merge(count_table_dge, annot, by.x = "gene",
                         by.y = "ensembl_gene_id")

# get the count table 
count_table_vaia <- as.data.frame(cpm(dgevaia))
count_table_vaia$gene = rownames(count_table_vaia)

# merge the tables
merge_tables <- merge(x = count_table_dge, y = count_table_vaia, 
                      by.x = "refseq_mrna", by.y = "gene")

# prepare the count table for the PCA
dt_merge_PCA <- merge_tables[,-c(1,2,grep(pattern = "external",
                                          colnames(merge_tables)))]
rownames(dt_merge_PCA) = merge_tables$refseq_mrna

# remove genes without more than 10 reads in all samples
dt_merge_PCA = dt_merge_PCA[which(rowSums(dt_merge_PCA) > 10),]

# prepare the samples
tmp_dge <- dge$samples[,1:3]
tmp_dge$group = dge$samples$SampleGroup
merge_dge <- rbind(tmp_dge, dgevaia_tmp$samples)

# create the DGEList object
dgeMerge <- DGEList(counts = dt_merge_PCA, group = as.factor(merge_dge$group))



# normalization and estimate dispersion
df.dgeMerge <- calcNormFactors(dgeMerge)
moma_vaia <- model.matrix(~ 0 + group, data=merge_dge) # ~ 0 + SampleGroup
#colnames(moma_vaia) <- sub("SampleGroup", "", colnames(moma_vaia))
dgeMerge_tmp <- calcNormFactors(dgeMerge)
dgeMerge_tmp <- estimateDisp(dgeMerge_tmp, design=moma_vaia)


# create the deseq2 object
dds_merge = DESeqDataSetFromMatrix(countData = round(dgeMerge_tmp$counts),
                                  colData = dgeMerge_tmp$samples, design = ~group)

transform.dds_merge <- rlog(dds_merge, blind = T)

# create the PCA data
pcaData <- DESeq2::plotPCA(transform.dds_merge, intgroup = c("group"), 
                           returnData =T, ntop = nrow(dds_merge))
percentVar <- round(100 * attr(pcaData, "percentVar")) # must be before to keep the percentVar information
pcaData <- cbind(pcaData, "SampleName" = name)


# do the PCA plot (https://github.com/mikelove/DESeq2/blob/master/R/plots.R)
pcaplot <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  #coord_fixed() +
  #scale_color_manual(values = rainbow(6))
  scale_colour_brewer(palette = "Paired") 


# save the PCAplot with the observed batch effect 
pdf(file = "DEA/Plots/PCA_vaia_HE_batch_effect.pdf", paper = "a4r")
pcaplot
dev.off()

# change the LTlate to LTearly for the sample 9 from Vaïa dataset
pcaData["09_LTlate_1_M2grLT3", "group"] = "LTearly"
pcaData["09_LTlate_1_M2grLT3", "name"] = "09_LTearly_1_M2grLT3"

# remove the duplicated columns
pcaData = pcaData[, -4]

# redo the pca without the GFPneg
pcaData_woGFPneg = pcaData[grep(pattern = "GFPneg", pcaData$group, invert = T),]

# plot the PCA
pcaplot <- ggplot(pcaData_woGFPneg, aes(PC1, PC2, color = group)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  #coord_fixed() +
  #scale_color_manual(values = rainbow(6))
  scale_colour_brewer(palette = "Set1") 


# save the PCAplot with the observed batch effect 
pdf(file = "DEA/Plots/PCA_vaia_HE_batch_effect_woGFPneg.pdf", paper = "a4r")
pcaplot
dev.off()


# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(transform.dds_merge)))

# save the barplot with the observed batch effect 
pdf(file = "DEA/Plots/barplot_PC_dimension.pdf", paper = "a4r")
fviz_eig(pca)
dev.off()

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

# create the group
intgroup.df <- as.data.frame(colData(transform.dds_merge)[, intgroup, drop=FALSE])

# add the intgroup factors together to create a new grouping factor
group <- if (length(intgroup) > 1) {
  factor(apply( intgroup.df, 1, paste, collapse=":"))
} else {
  colData(transform.dds_merge)[[intgroup]]
}

# assembly the data for the plot
d <- data.frame(PC1=pca$x[,2], PC2=pca$x[,3], group=group, intgroup.df, 
                name=colnames(transform.dds_merge))


attr(d, "percentVar") <- percentVar[2:3]

pcaplot <- ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + 
  geom_point(size=3) + 
  xlab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) +
  ylab(paste0("PC3: ",round(percentVar[3]*100),"% variance")) +
  coord_fixed() + scale_colour_brewer(palette = "Paired") 

# save the PCAplot with the observed batch effect 
pdf(file = "DEA/Plots/PCA_vaia_HE_batch_effect_PC2_PC3.pdf", paper = "a4r")
pcaplot
dev.off()



# correct batch effect ----------------------------------------------------

# source of code https://www.biostars.org/p/403053/

# add batch information in the design
dgeMerge_tmp$samples$batch = "Vaia"
dgeMerge_tmp$samples$batch[grep(pattern = "BSSE_QGF_", 
                                rownames(dgeMerge_tmp$samples))] = "HE"

# redo the conversion to DESeq object
dds_merge = DESeqDataSetFromMatrix(countData = round(dgeMerge_tmp$counts),
                                   colData = dgeMerge_tmp$samples, design = ~ batch + group)
# error message, it is not possible to compare the samples

transform.dds_merge <- rlog(dds_merge, blind = T)




# scatterplot Erg vs Evi1 -------------------------------------------------

# find the Erg and Evi1 genes
genes_selection <- annot[annot$external_gene_name %in% c("Erg", "Mecom", "Gapdh", "Dsp", "Hoxa9"),]

# normalized the read count
dgeMerge <- calcNormFactors(dgeMerge)

# cpm count
cpm_count <- cpm(dgeMerge, normalized.lib.sizes = T)

# select the genes of interest
cpm_count <- cpm_count[genes_selection$refseq_mrna[c(1,3,5,6)],]# 2 not 3
# Mecom: NM_007963
# genes_selection[c(1,3,5,6),]

# create the dataframe for the scatterplot
dt_plot <- as.data.frame(t(cpm_count))
dt_plot$sample <- rownames(dt_plot)
colnames(dt_plot) <- c("Mecom", "Hoxa9", "Dsp", "Erg", "sample")
dt_plot$group <- unlist(lapply(dt_plot$sample, function(x)strsplit(x, split= "_")[[1]][2]))
dt_plot$group[grep(pattern = "BMT", dt_plot$group)] = "AML-GMP"
dt_plot$group[grep(pattern = "QGF", rownames(dt_plot))] = as.vector(dgeMerge$samples$group)[grep(pattern = "QGF"
                                                , rownames(dgeMerge$samples))]

# change the sample 9 from Vaïa to LTearly
dt_plot["09_LTlate_1_M2grLT3", "group"] = "LTearly"
dt_plot["09_LTlate_1_M2grLT3", "sample"] = "09_LTearly_1_M2grLT3"


# remove the GFP neg
dt_plot <- dt_plot[grep(pattern = "GFP-", dt_plot$group, invert = T),]

# add the shape information
dt_plot$data = "TPO"
dt_plot$data[grep(pattern = "GMP", dt_plot$group)] = "Vaïa"
dt_plot$data[grep(pattern = "LT", dt_plot$group)] = "Vaïa"
dt_plot$data[grep(pattern = "noTPO", dt_plot$group)] = "noTPO"


# remove the wtGMP
dt_plot <- dt_plot[grep(pattern = "wt", dt_plot$group, invert = T),]

# create the scatterplot
scatterplot <- ggplot(data = dt_plot, aes(x = Mecom, y = Erg, color = group, shape = data, ext=sample)) + 
  geom_point(size=3) + theme_bw() + 
  #geom_text_repel(data = dt_plot, color="black", aes(label=sample)) +
  xlab("Mecom norm cpm") + ylab(paste0("Erg norm cpm")) +
  scale_colour_brewer(palette = "Set1") 

saveFigures(fileName = "scatterplot_HE_Vaia_Evi1_Erg_norm_woGFPneg_paper", 
            ggplot = scatterplot, dirPlot = "DEA/Plots/")
