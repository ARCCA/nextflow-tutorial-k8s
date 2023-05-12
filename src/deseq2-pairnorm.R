#!/usr/bin/env Rscript

## read command line arguments

args = commandArgs(trailingOnly=TRUE)
inputFilePath <- args[1]
targetsFile <- args[2]

## read in targets file and create all possible comparisons

targets <- read.csv(targetsFile, header=T)
sampleGroups <- as.character(unique(targets$sampleGroup))
comparisonTable <- data.frame(t(combn(sampleGroups, 2)))
allComparisons <- paste(comparisonTable$X1, "-vs-", comparisonTable$X2, sep="")

## load libraries

library(DESeq2)  
library(ape) 
library(ggplot2) 
library(dplyr)  
library(ggrepel)

## loop over the genecount and transcriptcount files

files = c("markdup.genecount", "markdup.transcriptcount")

for (file in files) {

  ## loop over each comparison

  for (i in allComparisons) {

    c1 <- unlist(strsplit(i, "-vs-"))
    myComparisons <- i

    ## read-in full targets and raw counts

    targets <- read.csv(targetsFile, header=T)
    rawCountsFile <- paste(inputFilePath, "/all.", file, ".txt", sep="")
    rawCounts <- read.delim(rawCountsFile, header=TRUE, check.names=F)
    allGeneLengths <- rawCounts[,c(1,2)]
    colnames(allGeneLengths) = c("gene_id","max.tx_len")

    ## reduce the targets and rawCounts down to the samples in the comparison
  
    targets2 = targets %>% filter(sampleGroup != c1[1]) %>% filter(sampleGroup != c1[2])
    outliers <- as.vector(targets2$analysisID)

    if (length(outliers) > 0) {
      targets <- dplyr::filter(targets, !analysisID %in% outliers)
      rawCounts <- dplyr::select(rawCounts, -outliers)
    }

    ## remove the annotation from the rawCounts to give only counts

    if (file == "markdup.genecount") {
      annotationNames <- rawCounts[,c(1:4)]
      rownames(rawCounts) = rawCounts[,1]
      rawCounts[,c(1:4)] <- NULL
    } else if (file == "markdup.transcriptcount") {
      annotationNames <- rawCounts[,c(1:5)]
      colnames(annotationNames)[1] = "ensemblGeneID"
      rownames(rawCounts) = rawCounts[,1]
      rawCounts[,c(1:5)] <- NULL
    }

    if (!identical(as.character(targets$analysisID), colnames(rawCounts))) {
      stop("colnames in sheets are not in the correct order")
    }  
 
    ## remove all row where there are zero counts in all samples

    rawCounts = rawCounts[rowSums(rawCounts > 0) >= 1, ]

    ## create the analysis objects

    sampleGroups <- targets$sampleGroup

    exptDesign = data.frame(
      row.names = colnames(rawCounts),
      condition = sampleGroups)

    exptObject <- DESeqDataSetFromMatrix(
      countData = rawCounts,
      colData = exptDesign,
      design = ~ condition)

    ## perform DESeq    

    analysisObject = DESeq(exptObject)

    ## add gene lengths to the output object

    geneLengths <- allGeneLengths[allGeneLengths$gene_id %in% rownames(analysisObject),]
    mcols(analysisObject)$basepairs <- geneLengths[match(rownames(analysisObject), geneLengths$gene_id),]$max.tx_len

    ## create the raw, normalised and FPKM objects
    
    rawCounts <- counts(analysisObject, normalized = FALSE)
    normalisedCounts <- counts(analysisObject, normalized = TRUE)
    fpkmNormalisedCounts <- fpkm(analysisObject, robust = TRUE)
  
    colnames(rawCounts) = gsub("^", "raw.counts.", colnames(rawCounts))
    colnames(normalisedCounts) = gsub("^", "norm.counts.", colnames(normalisedCounts))
    colnames(fpkmNormalisedCounts) = gsub("^", "fpkm.norm.counts.", colnames(fpkmNormalisedCounts))
  
    ## create the first part of the output file

    tempData = merge(rawCounts, normalisedCounts, by="row.names", all=T)
    rownames(tempData) <- tempData[,1]
    tempData[,1] <- NULL
    finalData = merge(tempData, fpkmNormalisedCounts, by="row.names", all=T)
    rownames(finalData) <- finalData[,1]
    finalData[,1] <- NULL
    printFPKMCounts <- data.frame(ensembleID=rownames(fpkmNormalisedCounts), fpkmNormalisedCounts)
     
    ## loop over all posible comparisons
 
    for (one in unique(targets$sampleGroup)) {
  	     
      for (two in unique(targets$sampleGroup)) {
  	       
        if (any(myComparisons == paste(one, "-vs-", two, sep=""))) {
    	                        
    	  result = as.data.frame(results(analysisObject, contrast=c("condition", one, two), independentFiltering=TRUE, pAdjustMethod="BH", alpha = 0.1))
    	            
    	  slimData = as.data.frame(cbind(rownames(finalData), finalData))
    	  colnames(slimData)[1] <-  "tracking_id"
    	         
    	  slimData = dplyr::select(slimData, tracking_id, colnames(dplyr::select(slimData, paste("raw.counts", unlist(subset(targets, sampleGroup == one, analysisID)), sep="."), paste("raw.counts", unlist(subset(targets, sampleGroup == two, analysisID)), sep="."), paste("norm.counts", unlist(subset(targets, sampleGroup == one, analysisID)), sep="."), paste("norm.counts", unlist(subset(targets, sampleGroup == two, analysisID)), sep="."), paste("fpkm.norm.counts", unlist(subset(targets, sampleGroup == one, analysisID)), sep="."), paste("fpkm.norm.counts", unlist(subset(targets, sampleGroup == two, analysisID)), sep="."))))
    	                
    	  printData = merge(slimData, result, by.x = "tracking_id", by.y = "row.names", all = T)
    	        
          if (file == "markdup.genecount") {
            tt <- merge(annotationNames, printData, by.x="ensemblGeneID", by.y="tracking_id")
            tt = tt[order(tt$pvalue), ]
            write.table(tt, file=paste("deseq2-pairnorm/output/", one, "_vs_", two, ".", file, ".pairnorm.txt", sep=""), row.names=F, sep="\t", quote=F)
          } else if (file == "markdup.transcriptcount") {  
            colnames(annotationNames)[1] = "ensemblTranscriptID"
            tt <- merge(annotationNames, printData, by.x="ensemblTranscriptID", by.y="tracking_id")
            tt = tt[order(tt$pvalue), ]
            write.table(tt, file=paste("deseq2-pairnorm/output/", one, "_vs_", two, ".", file, ".pairnorm.txt", sep=""), row.names=F, sep="\t", quote=F)
          }
        }
      }
    }
  }
}
