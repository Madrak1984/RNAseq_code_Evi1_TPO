# laod the libraries
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(xlsx)
library(biomaRt)
library("msigdbr")
library(edgeR)
require("ggrepel") # http://www.sthda.com/french/wiki/ggplot2-textes-ajouter-du-texte-a-un-graphique-logiciel-r-et-visualisation-de-donnees
library(fgsea)
library(tidyverse)
library(writexl)
library(VennDiagram)
library(svglite)

## getContentFiles ---------------------------------------------------------

######## function to create a list object which content all tables from the files found in the list_files argument
#   list_files, vector object of character values, contains names of all files to be read
#   directory, character value, contains the path oft the directory which contains all files, NULL by default
#   header, boolean value, indicates if there is a header in the file (TRUE by default)
#########
getContentFiles <- function(list_files, directory = NULL, header = T){
  
  
  # check the arguments
  if(!is.vector(list_files)) stop("list_files is not a vector object")
  
  # initialize list of tables
  list_tables <- NULL
  
  # list the complete path of files
  if(is.null(directory)){
    files <- list_files
  } else  files <- file.path(directory, list_files)
  
  
  
  
  # read all the files
  for(file in files){
    
    # display the message
    print(paste0(date(), ": reading of ", file, " in progress..."))
    
    if(length(grep(pattern = "csv$", file)) != 0){
      list_tables <- c(list_tables, list(read.csv(file = file, header = T, fill = T)))
    } else if(length(grep(pattern = "bed$", file)) != 0){
      list_tables <- c(list_tables, list(read.csv(file = file, header = F, fill = T, sep = "\t")))
    } else if((length(grep(pattern = "narrowPeak$", file)) != 0) | (length(grep(pattern = "broadPeak$", file)) != 0)){
      list_tables <- c(list_tables, list(read.table(file = file, header = F, sep = "\t", quote = "")))
    } else list_tables <- c(list_tables, list(read.table(file = file, header = header, sep = "\t", quote = "", fill = T)))
    
    # display the message
    print(paste0(date(), ": reading of ", file, " done."))
  } 
  
  # add the name of files within the list object
  names(list_tables) <- list_files
  
  return(list_tables)
}



# DE_analysis -------------------------------------------------------------
# create the function for the Differential analysis

DE_analysis <- function(Final_sce, day, cluster, output){
  
  
  # create the output directory
  path_stats_dir = file.path(output, paste0("Cluster",cluster,"_day_", day))
  if(!dir.exists(path_stats_dir)){
    dir.create(path_stats_dir, recursive = T)
  }
  
  # reduce the table
  Final_sce = Final_sce[, grep(pattern = day, colData(Final_sce)$days)]
  
  print(paste0(date(), " - Sum the umi per samples"))
  print(Final_sce)
  
  ## Sum of all cell from each sample --> Pseudobulk analysis based on Athimed's command lines
  # define names of group to compare
  splitCells <- factor(paste0(Final_sce$sampleName, "_", Final_sce$MLLAF9_expr))
  # count the cells per group
  num_cells <- table(splitCells)
  # exclude the group with less than 25 cells
  excludedCluster = names(num_cells[num_cells < 25]) 
  # remove the cells excluded
  if(!isEmpty(excludedCluster)){
    Final_sce.cell = Final_sce[, -which(splitCells %in% excludedCluster)]
    # redo the groups to compare
    splitCells <- factor(paste0(Final_sce.cell$sampleName, "_", Final_sce.cell$MLLAF9_expr))
    num_cells = table(splitCells)
  } else  Final_sce.cell = Final_sce
  
  if(ncol(Final_sce.cell) > 2){
    # sum the gene expression per cell in each group 
    sum_by_cluster <- vapply(split(colnames(Final_sce.cell), splitCells), function(x){ Matrix::rowSums(counts(Final_sce.cell)[, x]) }, numeric(length = nrow(counts(Final_sce.cell)))) 
    
    # do the count of cells per samples
    pheno <- data.frame(summed_sample=colnames(sum_by_cluster))
    row.names(pheno) <- pheno$summed_sample
    pheno$num_cells <- num_cells
    pheno$condition <- gsub(pattern = "^[0-9]*_MF9_d[0-9]*_", replacement = "", pheno$summed_sample)
    dim(pheno) # 118 samples
    sum(pheno$num_cells) # total number of cell: 43'670
    
    # draw the barplot
    png(filename = file.path(path_stats_dir, paste0("barplot_MLL-AF9_cluster", cluster, "_", day,"_countCells.png")))
    par(mar=c(5,15,1,1))
    barplot(pheno$num_cells,las=2,cex.names=0.8,names=pheno$summed_sample,horiz = T)
    dev.off()
    
    print(paste0(date(), " - Create the DGEList object"))
    
    # create the DGEList object
    dge <- DGEList(counts = sum_by_cluster, samples = pheno, group = pheno$condition, genes = as.data.frame(rowData(Final_sce)[,c('ID', 'Symbol', 'Entrez_ID', 'Symbol')]), remove.zeros = F)
    
    print(dim(dge))
    
    # normalisation
    dge <- calcNormFactors(dge)
    
    # Keep only the genes that are expressed in at least 2 samples:
    keep <- rowSums(edgeR::cpm(dge) > 1) >= 2 ## Sometimes we have only 2 replicates
    table(keep)
    
    # subset, while updating total counts
    dge <- dge[keep, , keep.lib.sizes=FALSE] 
    # Redo TMM normalization
    dge <- calcNormFactors(dge) 
    
    #Visualize CPMs after filtering
    log2cpm.raw <- edgeR::cpm(dge, log=TRUE, prior.count=8, normalized.lib.sizes=FALSE)
    # After TMM normalization
    log2cpm <- edgeR::cpm(dge, log=TRUE, prior.count=8, normalized.lib.sizes=TRUE) 
    # log	==> logical, if TRUE then log2 values are returned.
    # prior.count	==> average count to be added to each observation to avoid taking log of zero. Used only if log=TRUE.
    
    # plot the densities according to the groups or samples
    png(filename = file.path(path_stats_dir, paste0("plotDensities_MLL-AF9_cluster", cluster, "_", day,"_countCells_group.png")))
    plotDensities(log2cpm, group=dge$samples$group, col=rep(myPalette2[1:12]), legend="topright")
    dev.off()
    
    png(filename = file.path(path_stats_dir, paste0("plotDensities_MLL-AF9_cluster", cluster, "_", day,"_countCells_samples.png")))
    plotDensities(log2cpm, group=dge$samples$summed_sample, col=rep(myPalette2[1:12]), legend="topright")
    dev.off()
    
    ## PCA
    summary(apply(log2cpm, 1, var) == 0) ## Any gene with no variability?
    pca1 <- prcomp(t(log2cpm), scale = T)
    summary(pca1)
    
    # visualization for the PCA
    png(filename = file.path(path_stats_dir, paste0("plotPCA_MLL-AF9_cluster", cluster, "_", day,"_axis.png")))
    print(fviz_eig(pca1))
    dev.off()
    
    png(filename = file.path(path_stats_dir, paste0("plotPCA_MLL-AF9_cluster", cluster, "_", day,"_1vs2.png")))
    print(fviz_pca_ind(pca1,
                 col.ind = "cos2", # Colorer par le cos2
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     
    ))
    dev.off()
    png(filename = file.path(path_stats_dir, paste0("plotPCA_MLL-AF9_cluster", cluster, "_", day,"_1vs3.png")))
    print(fviz_pca_ind(pca1,
                 col.ind = "cos2", # Colorer par le cos2
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE,
                 axes = c(1,3)
    ))
    dev.off()
    png(filename = file.path(path_stats_dir, paste0("plotPCA_MLL-AF9_cluster", cluster, "_", day,"_2vs3.png")))
    print(fviz_pca_ind(pca1,
                 col.ind = "cos2", # Colorer par le cos2
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE,
                 axes = c(2,3)
    ))
    dev.off()
    #fviz_pca_ind(pca1, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE )
    
    # design matrix 
    moma <- model.matrix(~ 0 + group, data=dge$samples) 
    colnames(moma) <- gsub("group", "", colnames(moma)) 
    colnames(moma) <- make.names(colnames(moma))
    colnames(moma)
    
    print(paste0(date(), " - Define the contrasts"))
    
    test = makeContrasts(paste(colnames(moma), collapse = "-"), levels = moma)
    
    
    
    # Contrast matrix
    if(ncol(moma) < 3){
      contrasts.matrix <- makeContrasts(paste(colnames(moma), collapse = "-"), levels = moma)
      colnames(contrasts.matrix) = gsub(pattern = "-", replacement = "vs", colnames(contrasts.matrix))
    } else {
      contrasts.matrix <- makeContrasts(
        HighvsLow = high - low,
        HighvsNotDetected = high - Not_detected,
        LowvsNotDetected = low - Not_detected,
        levels = moma
      )
    }
    
    
    
    ## Tell which samples are used in each contrast
    fitPhenoData <- list()
    for (i in 1:ncol(contrasts.matrix)) {
      fitPhenoData[[colnames(contrasts.matrix)[i]]] <- list(
        group=levels(dge$samples$group)
      )
    }
    
    
    # perform the DE with glmQLFit based on the Athimed's code
    print(paste0(date(), " - Perform the differential analysis"))
    
    decideTestsDGE_glmQLFit = NULL
    fits <- list()
    for (i in 1:ncol(contrasts.matrix)) {
      contr.name <- colnames(contrasts.matrix)[i]
      message("   Working on ... ", contr.name)
      
      # initialize the comparison
      comp = paste(names(which(contrasts.matrix[,i] != 0)), collapse = "VS")
      
      ## Subset samples
      dge.temp <- dge[, dge$samples$group %in% names(which(contrasts.matrix[,i] != 0))]
      
      ## Subset and renormalize
      keep <- rowSums(edgeR::cpm(dge.temp) > 1) >= 2  & !dge$genes$ID %in% rowData(Final.sce.n14)$ID[rowData(Final.sce.n14)$high_in_ambient]
      dge.temp <- dge.temp[keep, , keep.lib.sizes=FALSE] 
      
      ## Add a filtering step on the min number of cells with expression!
      ## keep only the cell belonging to the contrast we are making
      sce.temp <- Final_sce[
      rowData(Final_sce)$ID %in% dge.temp$genes$ID, 
      Final_sce$MLLAF9_expr %in% names(which(contrasts.matrix[,i] != 0))]
      
      ## Should be detected in at least 5% of the cells
      keep <- rowSums(counts(sce.temp) > 0 ) > dim(sce.temp)[2]/100*5
      dge.temp <- dge.temp[keep, , keep.lib.sizes=FALSE] 
      
      print(dim(dge.temp)) # total 7857 genes for cluster3 ", day,"
      
      ## Redo TMM normalization
      dge.temp <- calcNormFactors(dge.temp) 
      ## After TMM normalization
      log2cpm.temp <- edgeR::cpm(dge.temp, log=TRUE, prior.count=8, normalized.lib.sizes=TRUE)
      
      ## Plot density 
      par(mfrow = c(1, 1))
      png(filename = file.path(path_stats_dir, paste0("plotDensities_MLL-AF9_cluster", cluster, "_", day,"_group_", comp,".png")))
      plotDensities(log2cpm.temp, group=dge.temp$samples$group, col=myPalette2[1:2], legend="topright")
      dev.off()
      png(filename = file.path(path_stats_dir, paste0("plotDensities_MLL-AF9_cluster", cluster, "_", day,"_name_", comp,".png")))
      plotDensities(log2cpm.temp, group=dge.temp$samples$summed_sample, legend="topright")
      dev.off()
      
      ## design matrix 
      moma.temp <- model.matrix(~ 0 + group, data=dge.temp$samples) 
      colnames(moma.temp) <- gsub("group", "", colnames(moma.temp)) 
      colnames(moma.temp) <- make.names(colnames(moma.temp))
      colnames(moma.temp)
      
      ## estimate dispersion (TMM norm object
      dge.temp <- estimateDisp(dge.temp, moma.temp)
      plotBCV(dge.temp, pch=1) ## Not nice looking!!
      #plotMDS(dge.temp)
      
      fit1 <- glmQLFit(dge.temp, moma.temp, prior.count=8) 
      lrt <- glmQLFTest(fit1, contrast=contrasts.matrix[names(which(contrasts.matrix[,i] != 0)), i])
      png(filename = file.path(path_stats_dir, paste0("plotMD_MLL-AF9_cluster", cluster, "_", day,"_name_", comp,".png")))
      plotMD(lrt)
      dev.off()
      
      # DE genes
      print(table(de <- decideTestsDGE(lrt)))
      decideT =  table(de <- decideTestsDGE(lrt))
      
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
      detags <- rownames(dge.temp)[as.logical(de)]
      plotSmear(lrt, de.tags=detags, main=contr.name)
      abline(h=c(-1, 1), col=myPalette2[2])
      
      ## boxplot of 12 top genes (if at least 12 below FDR 5%)
      par(mfrow = c(3, 4))
      for (gene in row.names(topTags(lrt, p.value = 0.05, n=12))){
        boxplot(log2cpm.temp[gene, ] ~ dge.temp$samples$group, 
                col=rep(myPalette2[i], each=2), main=dge.temp$genes[gene,]$SYMBOL, xlab="", ylab="") 
      }
      
      ## add results to list
      fits[[contr.name]] <- lrt
      rm(lrt)
      ## Save filtered DGEList and other info useful for later 
      fits[[contr.name]]$dge <- dge.temp
      fits[[contr.name]]$moma <- moma.temp
      
      ## add info needed by Shiny app
      fits[[contr.name]]$contrasts <- contrasts.matrix[names(which(contrasts.matrix[,i] != 0)), i, drop=FALSE]
      fits[[contr.name]]$fitPhenoData <- fitPhenoData[[contr.name]]
      
      
      
    }
  } else fits = "not enough samples for the pseudo-bulk comparison"
  
  
  
  return(fits)
  
}



# function to perform a differential enrichment analysis
###### function to perform enrichment analysis based on one gene list, saved the result in excel files
# common, vector object, contains EnsemblID of genes to analyze
# outputDir, character value, path to save output files
# mart, mart object, provide the mart object to obtain all names of genes
# universe, character vector, Ensembl ID of analyzed genes 
# simplify, boolean value, indicates if you want to simplify the result (TRUE by default)
# isMouse, boolean value, indicate if the function must use the mouse genome or the human genome (TRUE by default)
#####
EnrichmentTerm_Analysis <- function(common, outputDir, mart, universe = NULL, simplify = T, isMouse = T){
  
  
  # create the output directory
  dir_plot <- outputDir
  if(!dir.exists(dir_plot)){
    dir.create(dir_plot, recursive = T)
  }
  
  
  
  # create the table of genes
  genes <- getBM(filters = "ensembl_gene_id",
                 attributes = c("ensembl_gene_id","entrezgene_id", "external_gene_name"),
                 values = common,
                 mart = mart)
  if(!is.null(universe)){
    genes_uni <- getBM(filters = "ensembl_gene_id",
                   attributes = c("ensembl_gene_id","entrezgene_id", "external_gene_name"),
                   values = universe,
                   mart = mart)
  }
  
  
  # set the different object for the genome
  if(isMouse){
    OrgDb = org.Mm.eg.db
    species = "Mus musculus"
    organism  = "mmu"
  } else {
    OrgDb = org.Hs.eg.db
    species = "Homo sapiens"
    organism  = "hsa"
  }
  
  
  
  # https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html
  # documentation for clusterprofiler: https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
  # initialize the different enrichment analysis
  enrichments <- c("CC", "BP", "MF", "KEGG")
  
  # initialize the list to save the enrichment plots
  list_Enrich_term <- NULL
  enrich_name <- NULL
  
  # perform the enrichments
  for(enrich in enrichments){
    
    # print initial message
    print(paste0(date(), ": ", enrich, " enrichment analysis in progress..."))
    
    # check if it is for the KEGG analysis
    if(enrich != "KEGG"){
      
      # find the enrichment
      if(is.null(universe)){
        tmp_enrich <- enrichGO(gene = na.omit(genes$entrezgene_id[genes$ensembl_gene_id %in% common]), 
                               OrgDb = OrgDb, ont = enrich, readable = T)
      } else {
        tmp_enrich <- enrichGO(gene = na.omit(genes$entrezgene_id[genes$ensembl_gene_id %in% common]), 
                               OrgDb = OrgDb, ont = enrich, readable = T,
                               universe = as.character(na.omit(genes_uni$entrezgene_id[genes_uni$ensembl_gene_id %in% universe])))
      }
      
      
      # simplify the result
      if(simplify) tmp_enrich = clusterProfiler::simplify(tmp_enrich)

      # write the excel files
      if(!any(dim(as.data.frame(tmp_enrich)) == 0)){
        write.xlsx(as.data.frame(tmp_enrich), file = file.path(dir_plot, "GO_term_analysis.xlsx"), col.names = T, row.names = F, sheetName = enrich, append = T)
        gc()
      }
      
      # save the enrichment in the list
      if(!is.null(tmp_enrich)){
        list_Enrich_term <- c(list_Enrich_term, tmp_enrich)
        enrich_name <- c(enrich_name, enrich)
      }
    } else{
      
      # perform the KEGG enrichment analysis
      if(is.null(universe)){
        eKegg <- enrichKEGG(gene = as.character(as.vector(na.omit(genes$entrezgene_id[genes$ensembl_gene_id %in% common]))),
                            organism  = organism)
        } else {
        eKegg <- enrichKEGG(gene = as.character(as.vector(na.omit(genes$entrezgene_id[genes$ensembl_gene_id %in% common]))),
                            organism  = organism, universe = as.character(na.omit(genes_uni$entrezgene_id[genes_uni$ensembl_gene_id %in% universe])))
      }
      
      if(!is.null(eKegg)){
        eKegg <- setReadable(eKegg, OrgDb = OrgDb, keyType = "ENTREZID")
        # write the excel files
        if(!any(dim(as.data.frame(eKegg)) == 0)){
          write.xlsx(as.data.frame(eKegg), file = file.path(dir_plot, "KEGG_term_analysis.xlsx"), col.names = T, row.names = F, sheetName = "KEGG")
          gc()
        }
      
      }
      # save the enrichment in the list
      if(!is.null(eKegg)){
        list_Enrich_term <- c(list_Enrich_term, eKegg)
        enrich_name <- c(enrich_name, enrich)
        
      }
      
    } # end if KEGG
    
    # print final message
    print(paste0(date(), ": ", enrich, " enrichment analysis in done."))
    
  } # end loop
  
  # save the names of list
  names(list_Enrich_term) = enrich_name
  
  
  
  # perform the msig database enrichment analysis
  dbs <- c("H", paste0("C", 1:8))
  list_db <- NULL
  tmp_db <- NULL
  for(db in dbs){
    
    # print initial message
    print(paste0(date(), ": Msig ", db, " enrichment analysis in progress..."))
    
    # select the list gene for the database
    m_t2g <- msigdbr(species = species, category = db) %>% dplyr::select(gs_name, entrez_gene)
    
    # perform the enrichment analysis
    if(is.null(universe)){
      enr <- enricher(na.omit(genes$entrezgene_id[genes$ensembl_gene_id %in% common]), TERM2GENE = m_t2g)
    } else {
      enr <- enricher(na.omit(genes$entrezgene_id[genes$ensembl_gene_id %in% common]), 
                      TERM2GENE = m_t2g, universe = as.character(na.omit(genes_uni$entrezgene_id[genes_uni$ensembl_gene_id %in% universe])))
    }
    if(!is.null(enr)){
      enr <- setReadable(enr, OrgDb = OrgDb, keyType = "ENTREZID") # convert the genesymbol name
      list_db <- c(list_db, enr)
      tmp_db <- c(tmp_db, db)
    }
    
    # print final message
    print(paste0(date(), ": Msig ", db, " enrichment analysis in done."))
    
  }
  
  # save the list name
  names(list_db) = tmp_db
  
  
  
  #### draw the plot for the enrichment analysis
  
  #save the all dotplot for GO and KEGG
  for(plot in names(list_Enrich_term)){
    
    # check if there is a plot to draw
    if(!any(dim(list_Enrich_term[[plot]]) == 0)){
      
      # name for the title
      if(plot == "KEGG"){
        title = "dotplot for KEGG term"
      } else {
        title = paste0("dotplot for GO-", plot ," term")
      }
      
      # print the plot in file
      png(filename = file.path(dir_plot, paste0(plot, "_dotplot.png")), width = 800, height = 800)
      print(dotplot(list_Enrich_term[[plot]], showCategory=20)  + ggtitle(title) + theme(title = element_text(size = 20, face = "bold")))
      dev.off()
      ggsave(filename = file.path(dir_plot, paste0(plot, "_dotplot.svg")), plot = dotplot(list_Enrich_term[[plot]], showCategory=20)  + 
               ggtitle(title) + theme(title = element_text(size = 20, face = "bold")), 
             device = "svg",)
      #svg(filename = file.path(dir_plot, paste0(plot, "_dotplot.svg")), width = 800, height = 800)
      #dotplot(list_Enrich_term[[plot]], showCategory=20)  + ggtitle(title) + theme(title = element_text(size = 20, face = "bold"))
      #dev.off()
      #pdf(file = file.path(dir_plot, paste0(plot, "_dotplot.pdf")), width = 800, height = 800)
      #dotplot(list_Enrich_term[[plot]], showCategory=20)  + ggtitle(title) + theme(title = element_text(size = 20, face = "bold"))
      #dev.off()
      ggsave(filename = file.path(dir_plot, paste0(plot, "_dotplot.pdf")), plot = dotplot(list_Enrich_term[[plot]], showCategory=20)  + 
            ggtitle(title) + theme(title = element_text(size = 20, face = "bold")), 
            device = "pdf")
    }
    
    
  }
  
  
  # save all the plot for all the Msig databases
  for(db in names(list_db)){
    
    # test if there are result within the object
    if(nrow(as.data.frame(list_db[[db]])) != 0){
      
      # name for the title
      title = paste0("dotplot for the database ", db," in Msig")
      
      # print the plot in file
      png(file.path(dir_plot, paste0("Msig_", db, "_dotplot.png")), width = 800, height = 800)
      print(dotplot(list_db[[db]], showCategory=20)  + ggtitle(title) + theme(title = element_text(size = 20, face = "bold")))
      dev.off()
      
      # save the excel file
      write.xlsx(as.data.frame(list_db[[db]]), file = file.path(dir_plot, "MsigDatabase_term_analysis.xlsx"), col.names = T, row.names = F, append = T, sheetName = paste0("Msigdatabase", db))
      gc()
      
    } # end if
  }
  
  
  # extract list of terms only based on MLL target in C2 MsigDatabase
  subtable <- as.data.frame(list_db[["C2"]])
  subtable <- subtable[grep(pattern = "MLL", subtable$ID),]
  if(nrow(subtable) > 0){
    write.xlsx(subtable, file = file.path(dir_plot, "MsigDatabase_term_analysis.xlsx"),
               col.names = T, row.names = F, append = T, sheetName = "Msigdatabase C2, MLL")
    
    # create a subenrich table
    subEnrichtable = list_db[["C2"]]
    subEnrichtable@result = subEnrichtable@result[grep(pattern = "MLL", subEnrichtable$ID),]
    png(file.path(dir_plot, paste0("heatplot_MLL.png")), width = 1600, height = 400)
    print(heatplot(subEnrichtable, showCategory = 8) + theme(axis.text.x = element_text(face="bold", size = 11)))
    dev.off()
    
  }
  
  
}





# getVolcanoplot ----------------------------------------------------------


######   function to draw the volcanoplot
#   dgeobject, DGELRT or dataframe object, contains the Differential expression analysis done by EdgeR
#   FC, numeric value, indicate a threshold for the FC (NULL by default, and not implemented yet)
#   FDR, numeric value, threshold to change the color according to the FDR (0.05 by default)
#   genes, vector of character, genes names (SYMBOL) to display on the volcanoplot (NULL by default, see addGenes parameters)
#   addGenes, boolean value, indicate if genes must be added (TRUE by default, add significant genes by default)
#   setAxis, boolean value, indicate if we want to set the FC axis to be symetrical (FALSE by defaults)
#   nbGenes, numeric value, indicate the number of significant genes to display (NULL by default)
#   yaxis, numeric values, vector object, must indicate the limit for the y-axis with 2 values (minimal and maximal), NULL by default
#####
getVolcanoplot <-  function(dgeobject, FC = NULL, FDR= 0.05, genes = NULL,  addGenes = T, setAxis = F, nbGenes = NULL, yaxis = NULL){
  
  #### check the parameters
  
  # check dgeobject
  if(class(dgeobject) == "DGELRT"){
    dge_DEG <- as.data.frame(topTags(dgeobject, n = nrow(dgeobject)))
  } else if(is.data.frame(dgeobject)){
    dge_DEG <- dgeobject
  } else stop("dgeobject is not a dataframe or a DGELRT object (see EdgeR package documentation for more details")
  
  # check FDR
  if(!is.numeric(FDR)) stop("FDR must be numeric")
  
  # check genes
  if(!is.logical(addGenes)) stop("addGenes must be boolean")
  if(addGenes & (!is.null(genes))){
    if(is.character(genes)){
      if(!is.vector(genes)) stop("genes must be a vector")
    } else stop("genes must contain character")
  }
  
  # check setAxis
  if(!is.logical(setAxis)) stop("setAxis must be boolean")
  
  # check the nbGenes
  if(addGenes & (!is.null(nbGenes))){
    if(!is.numeric(nbGenes)) stop("nbGenes must be numeric")
  }
  
  # check the y-axis
  if((!is.numeric(yaxis)) & (!is.null(yaxis))){
    stop("y-axis is not a nemurical values!")
  } else if((is.numeric(yaxis)) & (!is.null(yaxis))){
    if(length(yaxis) != 2){
      stop("yaxis must be a vector of length 2.")
    }
  }
  
  # add the colors
  dge_DEG$color <- "darkblue"
  dge_DEG$color[dge_DEG$FDR <= FDR] <- "red"
  
  # select genes to display
  if(addGenes){
    
    # select the genes
    if(is.null(genes)){
      all_genes <- dge_DEG$SYMBOL[dge_DEG$FDR <= FDR]  
    }
    
    # select genes according to the count
    if(!is.null(nbGenes)){
      
      # select the dataframe according to the FDR
      tmp_data <- dge_DEG[dge_DEG$FDR <= FDR,]  
      
      # select the ordered genes
      tmp_genes <- tmp_data$SYMBOL[order(tmp_data$FDR)]
      
      # remove genes already selected
      if(!is.null(genes)) tmp_genes <- tmp_genes[!(tmp_genes %in% genes)]
      
      # select the genes according to the nbGenes
      if(nbGenes < length(tmp_genes)){
        tmp_genes <- head(tmp_genes, n= nbGenes)
      }
      
      # concatenate the genes
      if(!is.null(genes)){
        genes <- c(genes, tmp_genes)
      } else genes <- tmp_genes
    } else if(is.null(genes)) genes <- all_genes
  }
  
  volcano <- ggplot(dge_DEG, aes(x=logFC, y=-log10(FDR), color=color)) + 
    theme_bw() + geom_point(alpha=0.5, color=dge_DEG$color)
  
  # add genes  
  if(addGenes){
    volcano <- volcano +
      geom_text_repel(data=subset(dge_DEG[dge_DEG$SYMBOL %in% genes,]), color="black", aes(label=SYMBOL), size =5)
  } 
  
  
  # set the x-axis
  if(setAxis){
    
    # get the max values
    max_x <- ceiling(max(na.omit(abs(dge_DEG$logFC))))
    volcano <- volcano + xlim(-max_x, max_x)
  }
  #volcano <- volcano + ggtitle(title)
  
  # set the y-axis
  if(!is.null(yaxis)){
    volcano <- volcano + ylim(yaxis)
  }
  
  # return the plot
  return(volcano)
}





# save the figures --------------------------------------------------------

##### function to save figures for paper in png, svg and pdf formats
# filename, character value, name of the file to save without extension
# ggplot, ggplot object, contains the figure to save in the file
# dirPlot, character value (current working directory by default), indicate the path were the file must be saved
# A4, boolean value (FALSE by default), indicate if the figure must be saved in A4 landscape format
#####
saveFigures <- function(fileName, ggplot, dirPlot = getwd(), A4 = F){
  
  if(!dir.exists(dirPlot)){
    stop(paste0(dirPlot, " does not exist!"))
  }
  
  # save the figures in png, svg and pdf format
  # save in the svg format
  if(A4){
    svg(filename = file.path(dirPlot, paste0(fileName, ".svg")), width = 11.7, height = 8.3)
  } else {
    svg(filename = file.path(dirPlot, paste0(fileName, ".svg")))
  }
  print(ggplot)
  dev.off()
  
  # save in the png format
  if(A4){
    png(filename = file.path(dirPlot, paste0(fileName, ".png")), width = 297, height = 210, units = "mm", res = 300)
  } else {
    png(filename = file.path(dirPlot, paste0(fileName, ".png")))
  }
  print(ggplot)
  dev.off()
  
  # save in the pdf format
  if(A4){
    pdf(file = file.path(dirPlot, paste0(fileName, ".pdf")), paper = "a4r", width = 11.7, height = 8.3)
  } else {
    pdf(file = file.path(dirPlot, paste0(fileName, ".pdf")))
  }
  print(ggplot)
  dev.off()
  
}



# scatterplots ------------------------------------------------------------



# create the function to save the scatterplot
draw_scatterplot_cluster <- function(genes, table_DEG_cl3, table_DEG_cl6, filename, xlab ="logFC cluster 6", ylab="logFC cluster 3"){
  
  # create the dataframe for the scatterplot
  dt_d10 <- data.frame("genes" = genes, 
                       "symbol" = table_DEG_cl3[genes,]$Symbol, 
                       "logFC_cl3" = table_DEG_cl3[genes,]$log2FC,
                       "logFC_cl6" = table_DEG_cl6[genes,]$log2FC,
                       "Adj-pvalue <= 5%" = "N.S.", check.names = F)
  
  # update the color according to the significant
  clust3_sign = table_DEG_cl3[table_DEG_cl3$adj.P.Val <= 0.05,]$GeneId
  clust6_sign = table_DEG_cl6[table_DEG_cl6$adj.P.Val <= 0.05,]$GeneId  
  common_sign = intersect(clust3_sign, clust6_sign)
  dt_d10[dt_d10$genes %in% clust3_sign,]$`Adj-pvalue <= 5%` = "cluster 3"
  dt_d10[dt_d10$genes %in% clust6_sign,]$`Adj-pvalue <= 5%` = "cluster 6"
  dt_d10[dt_d10$genes %in% common_sign,]$`Adj-pvalue <= 5%` = "common clusters"
  
  # get the max lim
  max_xlim <- ceiling(max(na.omit(abs(dt_d10$logFC_cl6))))
  max_ylim <- ceiling(max(na.omit(abs(dt_d10$logFC_cl3))))
  
  # create the ggplot object
  gg <- ggplot(data = dt_d10, aes(x=logFC_cl6, y=logFC_cl3, color = `Adj-pvalue <= 5%`)) +
    geom_point(alpha=0.5) + 
    scale_colour_manual(values=c("N.S." = "gray50", "cluster 3" = "red", "cluster 6" = "blue", "common clusters" = "purple")) +
    theme_bw() +
    xlab(xlab) +
    ylab(ylab) +
    geom_text_repel(data=subset(dt_d10[dt_d10$symbol %in% c( "MLL-AF9","Meis1","Eya1","Hoxa9","Hoxa10", "Dsp", "Ikzf2", "Chd9", "Epha7", "Igf2bp3", "Slpi"),]), color="black", aes(label=symbol), size =5) +
    xlim(-max_xlim, max_xlim) +
    ylim(-max_ylim, max_ylim)
  
  
  # draw the scatterplot
  png(paste0(filename, ".png"))
  print(gg)
  dev.off()
  pdf(paste0(filename, ".pdf"))
  print(gg)
  dev.off()
  
}





# create the function to save the scatterplot
draw_scatterplot_day <- function(genes, table_DEG_d10, table_DEG_d60, filename){
  
  # create the dataframe for the scatterplot
  dt_d10 <- data.frame("genes" = genes, 
                       "symbol" = table_DEG_d10[genes,]$Symbol, 
                       "logFC_d10" = table_DEG_d10[genes,]$log2FC,
                       "logFC_d60" = table_DEG_d60[genes,]$log2FC,
                       "Adj-pvalue <= 5%" = "N.S.", check.names = F)
  
  # update the color according to the significant
  d10_sign = table_DEG_d10[table_DEG_d10$adj.P.Val <= 0.05,]$GeneId
  d60_sign = table_DEG_d60[table_DEG_d60$adj.P.Val <= 0.05,]$GeneId  
  common_sign = intersect(d10_sign, d60_sign)
  dt_d10[dt_d10$genes %in% d10_sign,]$`Adj-pvalue <= 5%` = "day 10"
  dt_d10[dt_d10$genes %in% d60_sign,]$`Adj-pvalue <= 5%` = "day 60"
  dt_d10[dt_d10$genes %in% common_sign,]$`Adj-pvalue <= 5%` = "common DEG"
  
  # get the max lim
  max_xlim <- ceiling(max(na.omit(abs(dt_d10$logFC_d10))))
  max_ylim <- ceiling(max(na.omit(abs(dt_d10$logFC_d60))))
  
  # create the ggplot object
  gg <- ggplot(data = dt_d10, aes(x=logFC_d10, y=logFC_d60, color = `Adj-pvalue <= 5%`)) +
    geom_point(alpha=0.5) + 
    scale_colour_manual(values=c("N.S." = "grey50", "day 10" = "red", "day 60" = "blue", "common DEG" = "purple")) +
    theme_bw() +
    xlab("logFC day 10") +
    ylab("logFC day 60") +
    geom_text_repel(data=subset(dt_d10[dt_d10$symbol %in% c( "MLL-AF9","Meis1","Eya1","Hoxa9","Hoxa10", "Dsp", "Ikzf2", "Chd9", "Epha7", "Igf2bp3", "Slpi"),]), color="black", aes(label=symbol), size =5) +
    xlim(-max_xlim, max_xlim) +
    ylim(-max_ylim, max_ylim)
  
  
  # draw the scatterplot
  png(paste0(filename, ".png"))
  print(gg)
  dev.off()
  pdf(paste0(filename, ".pdf"))
  print(gg)
  dev.off()
  
}


# get human gene annotation for orthology with mouse ----------------------

##### get human gene annotation for orthology with mouse, code provided by Julien Roux, must be updated for the host
#
#####
getHumanOrthologous <- function(){
  
  # load the ensembl gene information
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "sep2019.archive.ensembl.org")
  orthology <- getBM(attributes=c("ensembl_gene_id",
                                  "external_gene_name",
                                  "mmusculus_homolog_associated_gene_name",
                                  "mmusculus_homolog_ensembl_gene",
                                  "mmusculus_homolog_orthology_type",
                                  "mmusculus_homolog_subtype",
                                  "mmusculus_homolog_orthology_confidence"),
                     mart = ensembl)
  
  
  # Some filtering to keep only 1-to-1 orthologs
  orthology1to1 <- orthology[orthology$mmusculus_homolog_ensembl_gene != "",]
  orthology1to1 <- orthology1to1[orthology1to1$mmusculus_homolog_orthology_type %in% c("ortholog_one2one"),]
  
  
  ## Filter on taxonomic level. We could be more stringent here but would lose some orthologs identified as "dubious duplications"
  orthology1to1 <- orthology1to1[orthology1to1$mmusculus_homolog_subtype %in% c("Boreoeutheria", "Eutheria", "Euarchontoglires", "Theria", "Amniota"), ]
  row.names(orthology1to1) <- orthology1to1$ensembl_gene_id
  
  return(orthology1to1)
  
}


# get mouse gene annotations ----------------------

##### get mouse gene annotation
#
#####
getMouseGene <- function(withProtein = F){
  
  # load the ensembl gene information
  mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  
  # list all genes
  if(withProtein){
    genes_mm <- getBM(attributes = c("ensembl_gene_id", 
                                     "entrezgene_id", 
                                     "external_gene_name", 
                                     "description", 
                                     "uniprotswissprot"), 
                      mart = mart)
  } else {
    genes_mm <- getBM(attributes = c("ensembl_gene_id", 
                                     "entrezgene_id", 
                                     "external_gene_name"), 
                      mart = mart)
    
  }
  
  return(genes_mm)
}



# get human gene annotations ----------------------

##### get human gene annotation
#
#####
getHumanGene <- function(withProtein = F){
  
  # load the ensembl gene information
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  # list all genes
  if(withProtein){
    genes_mm <- getBM(attributes = c("ensembl_gene_id", 
                                     "entrezgene_id", 
                                     "external_gene_name", 
                                     "description", 
                                     "uniprotswissprot"), 
                      mart = mart)
  } else {
    genes_mm <- getBM(attributes = c("ensembl_gene_id", 
                                     "entrezgene_id", 
                                     "external_gene_name"), 
                      mart = mart)
    
  }
  
  return(genes_mm)
}



# internal functions ------------------------------------------------------

##### function to save the PCA in files
# dge, DGEList object, result of differential expression analysis done with edgeR
# extensionName, character value, extension of the file name, "" by default
# plotDir, character value, path where saved the files
# theme_bw, boolean value, TRUE to draw the black and white background (FALSE by default)
# text, boolean value, TRUE to add name of point, (TRUE by default)
#####
savePCA <- function(dge_tmp, extensionName = "", plotDir, theme_bw = F, text = T){
  
  # remove genes without more than 10 reads in all samples
  dge = dge_tmp[which(rowSums(dge_tmp$counts) > 10),]
  
  # normalization and estimate dispersion
  df.dge <- calcNormFactors(dge)
  moma <- model.matrix(~ 0 + SampleGroup, data=dge$samples) # ~ 0 + SampleGroup
  colnames(moma) <- sub("SampleGroup", "", colnames(moma))
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design=moma)
  
  # plot a PCA with Deseq2
  dds = DESeqDataSetFromMatrix(countData = round(df.dge$counts),
                               colData = df.dge$samples, design = ~SampleGroup)
  transform.dds <- rlog(dds, blind = T)
  
  # plotPCA for the top 1000
  plotPCA(filename = paste0("PCA_top1000", extensionName), 
          transform.dds_pca = transform.dds, ntop = 1000, plotDir = plotDir, 
          SampleName = colData(dds)$ExternalSampleName, theme_bw = theme_bw, 
          text = text)
  
  # plotPCA for the top 1000
  plotPCA(filename = paste0("PCA_top500", extensionName), 
          transform.dds_pca = transform.dds, ntop = 500, plotDir = plotDir,
          SampleName = colData(dds)$ExternalSampleName, theme_bw = theme_bw, 
          text = text)
  
  # plotPCA for all genes 
  plotPCA(filename = paste0("PCA_all", extensionName), 
          transform.dds_pca = transform.dds, ntop = nrow(transform.dds), 
          plotDir = plotDir, SampleName = colData(dds)$ExternalSampleName,  
          theme_bw = theme_bw, text = text)
  
}


##### functions which create the PCAs
# filename, character value, file name
# transform.dds_pca, DESeqTransform object
# ntop, numeric value, indicate the number of top genes to select for the PCA plot
# plotDir, character value, path where saved the files
# theme_bw, boolean value, TRUE to draw the black and white background (FALSE by default)
# text, boolean value, TRUE to add name of point, (TRUE by default)
##### save the PCA in png and pdf file
plotPCA <- function(filename, transform.dds_pca, ntop, plotDir, SampleName, 
                    theme_bw = F, text = T){
  
  # create the PCA data
  pcaData <- DESeq2::plotPCA(transform.dds_pca, intgroup = c("SampleGroup"), returnData =T, ntop = ntop)
  percentVar <- round(100 * attr(pcaData, "percentVar")) # must be before to keep the percentVar information
  pcaData <- cbind(pcaData, "SampleName" = SampleName)
  
  # do the PCA plot
  pcaplot <- ggplot(pcaData, aes(PC1, PC2, color = group)) +
    geom_point(size = 5) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
    #coord_fixed() +
    #scale_color_manual(values = rainbow(6))
    scale_fill_brewer(palette = "Set1")
  
  # add the text
  if(text){
    pcaplot <- pcaplot + geom_text_repel(data = pcaData, color="black", 
                                         aes(label=SampleName))
  }
  
  # add the theme black and white
  if(theme_bw){
    pcaplot <- pcaplot + theme_bw()
  }
  
  # save the PCA in the png and pdf files
  png(file.path(plotDir, paste0(filename,".png")))
  print(pcaplot)
  dev.off()
  pdf(file.path(plotDir, paste0(filename,".pdf")))
  print(pcaplot)
  dev.off()
  svg(file.path(plotDir, paste0(filename,".svg")))
  print(pcaplot)
  dev.off()
}


##### functions to save the barplot in files
# data, dataframe object, dataframe used by the ggplot2
# plotDir, character value, path of the directory to save the figure
# fileName, character value, prefix of name file
##### save the barplot in the file
savebarplot <- function(data, plotDir, fileName){
  
  # create the barplot
  bar_cpm <- ggplot(data = data, aes(x = samples, y = norm_cpm, fill = group)) + 
    geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ylab("normalized cpm")
  
  # save the barplot in files
  saveFigures(fileName = fileName, ggplot = bar_cpm, dirPlot = plotDir)
  
}


##### functions to save the lineplot in files
# data, dataframe object, dataframe used by the ggplot2
# plotDir, character value, path of the directory to save the figure
# fileName, character value, prefix of name file
##### save the lineplot in the file
savelineplot <- function(data, plotDir, fileName){
  
  # create the barplot
  line_cpm <- ggplot(data = data, aes(x = mouse, y = norm_cpm, group = group,
                                      colour = group, shape = treatment)) + geom_line() + geom_point() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ylab("normalized cpm")
  
  # save the barplot in files
  saveFigures(fileName = fileName, ggplot = line_cpm, dirPlot = plotDir)
  
}


##### functions to create the barplot and save it in files
# pattern, character value, gene to select
# dge, dge object
# plotDir, character value, path of the directory to save the figure
##### save the barplot in the file
createPlots <- function(pattern, dge, plotDir){
  
  # calculate the normalized read count
  norm_count_dge <- cpm(dge)
  
  # select the mecom count
  count <- norm_count_dge[grep(pattern = paste0(pattern, "$"), dge$genes$SYMBOL, 
                               ignore.case = T),]
  
  # change the sample name
  names(count) <- paste(dge$samples$ExternalSampleName, 
                        unlist(lapply(dge$samples$SampleGroup, 
                                      function(x) strsplit(x, "_")[[1]][1])), sep = "_")
  
  # reorder the sample name
  count <-  count[order(count, decreasing = T)]
  
  # create the dataframe for the mecom gene
  dt_count <- data.frame("samples" = paste(unlist(lapply(names(count), 
                                                         function(x) strsplit(x, "_")[[1]][1])), unlist(lapply(names(count), 
                                                                                                               function(x) strsplit(x, "_")[[1]][2])), sep="_"), 
                         "norm_cpm" = count, 
                         "group" = unlist(lapply(names(count), 
                                                 function(x) strsplit(x, "_")[[1]][2])),
                         "mouse" = unlist(lapply(names(count), 
                                                 function(x) strsplit(x, "_")[[1]][1])),
                         "treatment" = unlist(lapply(names(count), 
                                                     function(x) strsplit(x, "_")[[1]][3]))) 
  
  # keep the order of samples
  dt_count$samples <- factor(dt_count$samples, levels = dt_count$samples)
  dt_count$mouse <- factor(dt_count$mouse, levels = unique(dt_count$mouse))
  dt_count$treatment <- factor(dt_count$treatment, levels = unique(dt_count$treatment))
  
  # create the file names
  filenames <- paste(c("barplot", "lineplot"), str_replace(pattern, "^\\w{1}",
                                                           toupper), "normCPM", sep = "_")
  
  # create the barplot
  savebarplot(data = dt_count, plotDir = plotDir, fileName = filenames[1])
  
  # create the lineplot
  savelineplot(data = dt_count, plotDir = plotDir, fileName = filenames[2])
  
  # remove the indexes for the 35088 & 624-2 samples
  id_sample <- intersect(grep(pattern = "35088", dt_count$samples, invert = T), 
                         grep(pattern = "624-2", dt_count$samples, invert = T))
  
  # remove the 35088 & 624-2 samples
  dt_count_wosamples <- dt_count[id_sample,]
  
  # create the barplot
  savebarplot(data = dt_count_wosamples, plotDir = plotDir, 
              fileName = paste0(filenames[1], "_wo35088--624-2"))
  
  # create the lineplot
  savelineplot(data = dt_count_wosamples, plotDir = plotDir,
               fileName = paste0(filenames[2], "_wo35088--624-2"))
  
}

##### function to create a metatable
# fits, list of edgeR results
# FDR, numerical value, threshold for the FDR (0.05 by default)
##### return a metatable which summarize the DEA results
createMetatable <- function(fits, FDR = 0.05){
  
  # initialize the metatable
  DEG_table <- fits[[1]]
  dge_samples_DEG <- as.data.frame(topTags(DEG_table, n = nrow(DEG_table)))
  metatable <- dge_samples_DEG[,c(1:7)]
  
  # colnames to select
  col_name <- c("EnsemblID", "logFC", "PValue", "FDR")
  
  # add ensembl_ID
  metatable$EnsemblID = rownames(metatable)
  
  # add DEG for GFP+.TPOvsnoTPO
  DGE_tmp = as.data.frame(topTags(fits[["GFP.TPO"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'GFP+: TPO vs noTPO' = F
  metatable$`GFP+: TPO vs noTPO`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`GFP+: TPO vs noTPO`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC < 0)]] = "down"
  
  # add DEG for Bulk.TPOvsnoTPO
  DGE_tmp = as.data.frame(topTags(fits[["bulk.TPO"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'Bulk: TPO vs noTPO' = F
  metatable$`Bulk: TPO vs noTPO`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`Bulk: TPO vs noTPO`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC < 0)]] = "down"
  
  # merge the statistics values
  X = as.data.frame(topTags(fits[["GFP.TPO"]], n = nrow(DEG_table)))
  X$EnsemblID = rownames(X)
  Y = as.data.frame(topTags(fits[["bulk.TPO"]], n = nrow(DEG_table)))
  Y$EnsemblID = rownames(Y)
  tmp_stat <- merge(x=X[,col_name],
                    y=Y[,col_name], 
                    by="EnsemblID", suffixes = c('_GFP+: TPO vs noTPO', '_Bulk: TPO vs noTPO'))
  metatable <- merge(x = metatable, y = tmp_stat, by = "EnsemblID")
  
  
  # add DEG for GFP-.TPOvsnoTPO
  DGE_tmp = as.data.frame(topTags(fits[["GFPneg.TPO"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'GFP-: TPO vs noTPO' = F
  metatable$`GFP-: TPO vs noTPO`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`GFP-: TPO vs noTPO`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC < 0)]] = "down"
  
  # merge the statistics values
  Y = as.data.frame(topTags(fits[["GFPneg.TPO"]], n = nrow(DEG_table)))
  Y$EnsemblID = rownames(Y)
  metatable <- merge(x = metatable, y = Y[,col_name], by = "EnsemblID")
  # change the last column name
  colnames(metatable)[c((ncol(metatable)-2) : ncol(metatable))] = paste(colnames(metatable)[c((ncol(metatable)-2) : ncol(metatable))], "_GFP-: TPO vs noTPO")
  
  # add DEG for TPO.GFP+vsbulk
  DGE_tmp = as.data.frame(topTags(fits[["TPO.GFPvsbulk"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'TPO: GFP+ vs bulk' = F
  metatable$`TPO: GFP+ vs bulk`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`TPO: GFP+ vs bulk`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC < 0)]] = "down"
  
  # add DEG for noTPO.GFP+vsbulk
  DGE_tmp = as.data.frame(topTags(fits[["noTPO.GFPvsbulk"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'noTPO: GFP+ vs bulk' = F
  metatable$`noTPO: GFP+ vs bulk`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`noTPO: GFP+ vs bulk`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC < 0)]] = "down"
  
  # merge the statistics values
  X=as.data.frame(topTags(fits[["TPO.GFPvsbulk"]], n = nrow(DEG_table)))
  X$EnsemblID = rownames(X)
  Y=as.data.frame(topTags(fits[["noTPO.GFPvsbulk"]], n = nrow(DEG_table)))
  Y$EnsemblID = rownames(Y)
  tmp_stat <- merge(x=X[,col_name],
                    y=Y[,col_name], 
                    by="EnsemblID", suffixes = c('_TPO: GFP+ vs bulk', '_noTPO: GFP+ vs bulk'))
  metatable <- merge(x = metatable, y = tmp_stat, by = "EnsemblID")
  
  
  
  # add DEG for TPO.GFP+vsGFP-
  DGE_tmp = as.data.frame(topTags(fits[["TPO.GFPvsGFPneg"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'TPO: GFP+ vs GFP-' = F
  metatable$`TPO: GFP+ vs GFP-`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`TPO: GFP+ vs GFP-`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC < 0)]] = "down"
  
  
  # add DEG for noTPO.GFP+vsGFP-
  DGE_tmp = as.data.frame(topTags(fits[["noTPO.GFPvsGFPneg"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'noTPO: GFP+ vs GFP-' = F
  metatable$`noTPO: GFP+ vs GFP-`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`noTPO: GFP+ vs GFP-`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$FDR <= FDR) & (DGE_tmp$logFC < 0)]] = "down"
  
  
  # merge the statistics values
  X=as.data.frame(topTags(fits[["TPO.GFPvsGFPneg"]], n = nrow(DEG_table)))
  X$EnsemblID = rownames(X)
  Y=as.data.frame(topTags(fits[["noTPO.GFPvsGFPneg"]], n = nrow(DEG_table)))
  Y$EnsemblID = rownames(Y)
  tmp_stat <- merge(x=X[,col_name],
                    y=Y[,col_name], 
                    by="EnsemblID", suffixes = c('_TPO: GFP+ vs GFP-', '_noTPO: GFP+ vs GFP-'))
  metatable <- merge(x = metatable, y = tmp_stat, by = "EnsemblID")
  
  # return the metatable
  return(metatable)
}


##### function to create a metatable (must be factorize with previous function)
# fits, list of edgeR results
# Pvalue, numerical value, threshold for the Pvalue (0.05 by default)
##### return a metatable which summarize the DEA results
createMetatable_pval <- function(fits, Pvalue = 0.05){
  # initialize the metatable
  DEG_table <- fits[[1]]
  dge_samples_DEG <- as.data.frame(topTags(DEG_table, n = nrow(DEG_table)))
  metatable <- dge_samples_DEG[,c(1:7)]
  
  # colnames to select
  col_name <- c("EnsemblID", "logFC", "PValue", "FDR")
  
  # add ensembl_ID
  metatable$EnsemblID = rownames(metatable)
  
  # add DEG for GFP+.TPOvsnoTPO
  DGE_tmp = as.data.frame(topTags(fits[["GFP.TPO"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'GFP+: TPO vs noTPO' = F
  metatable$`GFP+: TPO vs noTPO`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`GFP+: TPO vs noTPO`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC < 0)]] = "down"
  
  # add DEG for Bulk.TPOvsnoTPO
  DGE_tmp = as.data.frame(topTags(fits[["bulk.TPO"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'Bulk: TPO vs noTPO' = F
  metatable$`Bulk: TPO vs noTPO`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`Bulk: TPO vs noTPO`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC < 0)]] = "down"
  
  # merge the statistics values
  X = as.data.frame(topTags(fits[["GFP.TPO"]], n = nrow(DEG_table)))
  X$EnsemblID = rownames(X)
  Y = as.data.frame(topTags(fits[["bulk.TPO"]], n = nrow(DEG_table)))
  Y$EnsemblID = rownames(Y)
  tmp_stat <- merge(x=X[,col_name],
                    y=Y[,col_name], 
                    by="EnsemblID", suffixes = c('_GFP+: TPO vs noTPO', '_Bulk: TPO vs noTPO'))
  metatable <- merge(x = metatable, y = tmp_stat, by = "EnsemblID")
  
  
  # add DEG for GFP-.TPOvsnoTPO
  DGE_tmp = as.data.frame(topTags(fits[["GFPneg.TPO"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'GFP-: TPO vs noTPO' = F
  metatable$`GFP-: TPO vs noTPO`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`GFP-: TPO vs noTPO`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC < 0)]] = "down"
  
  # merge the statistics values
  Y = as.data.frame(topTags(fits[["GFPneg.TPO"]], n = nrow(DEG_table)))
  Y$EnsemblID = rownames(Y)
  metatable <- merge(x = metatable, y = Y[,col_name], by = "EnsemblID")
  # change the last column name
  colnames(metatable)[c((ncol(metatable)-2) : ncol(metatable))] = paste(colnames(metatable)[c((ncol(metatable)-2) : ncol(metatable))], "_GFP-: TPO vs noTPO")
  
  # add DEG for TPO.GFP+vsbulk
  DGE_tmp = as.data.frame(topTags(fits[["TPO.GFPvsbulk"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'TPO: GFP+ vs bulk' = F
  metatable$`TPO: GFP+ vs bulk`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`TPO: GFP+ vs bulk`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC < 0)]] = "down"
  
  # add DEG for noTPO.GFP+vsbulk
  DGE_tmp = as.data.frame(topTags(fits[["noTPO.GFPvsbulk"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'noTPO: GFP+ vs bulk' = F
  metatable$`noTPO: GFP+ vs bulk`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`noTPO: GFP+ vs bulk`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC < 0)]] = "down"
  
  # merge the statistics values
  X=as.data.frame(topTags(fits[["TPO.GFPvsbulk"]], n = nrow(DEG_table)))
  X$EnsemblID = rownames(X)
  Y=as.data.frame(topTags(fits[["noTPO.GFPvsbulk"]], n = nrow(DEG_table)))
  Y$EnsemblID = rownames(Y)
  tmp_stat <- merge(x=X[,col_name],
                    y=Y[,col_name], 
                    by="EnsemblID", suffixes = c('_TPO: GFP+ vs bulk', '_noTPO: GFP+ vs bulk'))
  metatable <- merge(x = metatable, y = tmp_stat, by = "EnsemblID")
  
  
  
  # add DEG for TPO.GFP+vsGFP-
  DGE_tmp = as.data.frame(topTags(fits[["TPO.GFPvsGFPneg"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'TPO: GFP+ vs GFP-' = F
  metatable$`TPO: GFP+ vs GFP-`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`TPO: GFP+ vs GFP-`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC < 0)]] = "down"
  
  
  # add DEG for noTPO.GFP+vsGFP-
  DGE_tmp = as.data.frame(topTags(fits[["noTPO.GFPvsGFPneg"]], n = nrow(DEG_table)))
  DGE_tmp$EnsemblID = rownames(DGE_tmp)
  metatable$'noTPO: GFP+ vs GFP-' = F
  metatable$`noTPO: GFP+ vs GFP-`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC > 0)]] = "up"
  metatable$`noTPO: GFP+ vs GFP-`[metatable$EnsemblID %in% DGE_tmp$EnsemblID[(DGE_tmp$PValue <= Pvalue) & (DGE_tmp$logFC < 0)]] = "down"
  
  
  # merge the statistics values
  X=as.data.frame(topTags(fits[["TPO.GFPvsGFPneg"]], n = nrow(DEG_table)))
  X$EnsemblID = rownames(X)
  Y=as.data.frame(topTags(fits[["noTPO.GFPvsGFPneg"]], n = nrow(DEG_table)))
  Y$EnsemblID = rownames(Y)
  tmp_stat <- merge(x=X[,col_name],
                    y=Y[,col_name], 
                    by="EnsemblID", suffixes = c('_TPO: GFP+ vs GFP-', '_noTPO: GFP+ vs GFP-'))
  metatable <- merge(x = metatable, y = tmp_stat, by = "EnsemblID")
  
  # return the metatable
  return(metatable)
  
}



# getGeneListforGSEA(fit) -----------------------------------------------------------

###### function to a vector which contains a list genes for a GSEA analysis based on fit obejct (provided by EdgeR)
###### return a vector object, with decreased rank FC values named with the EntrezId of corresponding genes
# fit, DGELRT object, contains DEA results provided by EdgeR object (must contains entrezID names)
#####
getGeneListforGSEA <- function(fit){
  
  # check if the fit contains a entrezID column
  if(!any(colnames(fit$genes) == "ENTREZID")) stop("No ENTREZID column in the Gene table of the fit object!")
  
  # select genes with entrezID
  genes <- fit$genes[!is.na(fit$genes$ENTREZID),]
  
  # select the FC
  subtable <- fit[!is.na(fit$genes$ENTREZID),]
  id_FC <- sort(subtable$table$logFC, index.return = T, decreasing = T)
  FC <-  subtable$table$logFC[id_FC$ix]
  names(FC) <- subtable$genes$ENTREZID[id_FC$ix]
  
  # removes duplicate
  FC = FC[!duplicated(names(FC))]
  
  # return the vector
  return(FC)
}


# GSEA_hallmark(fit, list_dbs, filename, path = getwd()) -----------------------------------------------------------



# function to perform a GSEA based on hallmark database
###### function to perform GSEA based on fit object (provided by EdgeR), saved the result in excel files
# fit, DGELRT object, contains DEA results provided by EdgeR object (must contains entrezID names)
# list_dbs, vector of character values, must contain of database available in the environment (i.e c("Mm.H", "Mm.c2", etc...))
# filename, character value, names of the excel file
# path, character value, name of the directory which mut contain the excel file 
#####
GSEA_hallmark <- function(fit, list_dbs, filename, path = getwd()){
  
  # check if the output directory exist
  if(!dir.create(path)) dir.create(path, recursive = T)
  
  # create the gene list
  FC <- getGeneListforGSEA(fit = fit)
  
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
  write_xlsx(list_fGSEA_results, path = file.path(path, filename))
}





draw_VennDiagram <- function(x, category, filename, color = NULL, height = 480, width = 480, lty = 'blank', col = "black", imagetype="png"){
  venn.diagram(
    x = x,
    category.names = category,
    filename = filename,
    output=TRUE,
    
    # Output features
    imagetype=imagetype ,
    height = height , 
    width = width , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 1,
    lty = lty,
    fill = color,
    col = col,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.3,
    cat.fontface = "bold",
    cat.default.pos = "outer")
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# create_dir(path) -----------------------------------------------------------

# function to create a directory
###### function to create a directory
# path, character value, indicate the path to create if the folders did not exist yet
##### return path
create_dir <- function(path){
  
  if(is.character(path)){
    
    if(is.vector(path)){
      
      if(!dir.exists(paths = path)) dir.create(path = path, recursive = T)
      
    } else stop("path must be a vector!")
    
  } else stop("path must be a character!")
  
  return(path)
}


