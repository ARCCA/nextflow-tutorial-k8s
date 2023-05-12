#!/usr/bin/env Rscript

## read command line arguments

args = commandArgs(trailingOnly=TRUE)
inputFilePath <- args[1]
targetsFile <- args[2]

## read in targets file and create all possible comparisons

targets <- read.csv(targetsFile, header=T)
sampleGroups <- as.character(unique(targets$sampleGroup))
comparisonTable <- data.frame(t(combn(sampleGroups, 2)))
myComparisons <- paste(comparisonTable$X1, "-vs-", comparisonTable$X2, sep="")

## load libraries

library(DESeq2)
library(ape)   
library(ggplot2)
library(dplyr)
library(ggrepel)

## loop over the genecount and transcriptcount files

files = c("markdup.genecount", "markdup.transcriptcount")
  
for (file in files) {

  ## read counts file
  
  rawCountsFile <- paste(inputFilePath, "/all.", file, ".txt", sep="")
  rawCounts <- read.delim(rawCountsFile, header=TRUE, check.names=F)

  ## get gene lengths from rawCounts 

  allGeneLengths <- rawCounts[,c(1,2)]
  colnames(allGeneLengths) = c("gene_id", "max.tx_len")
    
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
    
  significantGenesOne <- NULL
  significantGenesFive <- NULL
  
  ## loop over all posible comparisons
 
  for (one in unique(targets$sampleGroup)) {

      for (two in unique(targets$sampleGroup)) {

        if (any(myComparisons == paste(one, "-vs-", two, sep=""))) {

            result = as.data.frame(results(analysisObject, contrast=c("condition", one, two), independentFiltering=TRUE, pAdjustMethod="BH"))

            slimData = as.data.frame(cbind(rownames(finalData), finalData))

            colnames(slimData)[1] <-  "tracking_id"

            slimData = dplyr::select(slimData, tracking_id, colnames(dplyr::select(slimData, paste("raw.counts", unlist(subset(targets, sampleGroup == one, analysisID)), sep="."), paste("raw.counts", unlist(subset(targets, sampleGroup == two, analysisID)), sep="."), paste("norm.counts", unlist(subset(targets, sampleGroup == one, analysisID)), sep="."), paste("norm.counts", unlist(subset(targets, sampleGroup == two, analysisID)), sep="."), paste("fpkm.norm.counts", unlist(subset(targets, sampleGroup == one, analysisID)), sep="."), paste("fpkm.norm.counts", unlist(subset(targets, sampleGroup == two, analysisID)), sep="."))))

            printData <- merge(slimData, result, by.x = "tracking_id", by.y = "row.names", all = T)

            if (file == "markdup.genecount") {
              tt <- merge(annotationNames, printData, by.x="ensemblGeneID", by.y="tracking_id")
              tt <- tt[order(tt$pvalue), ]

              ## write data to output

              write.table(tt, file=paste("deseq2-allnorm/output/", one, "_vs_", two, ".", file, ".allnorm.txt", sep=""), row.names=F, sep="\t", quote=F)
              geneList <- as.vector(unlist(subset(tt, padj <= 0.05, select=c("ensemblGeneID"))))
              significantGenesFive <- c(significantGenesFive, geneList)
              geneList <- as.vector(unlist(subset(tt, padj <= 0.01, select=c("ensemblGeneID"))))
              significantGenesOne <- c(significantGenesOne, geneList)
   	    } else if (file == "markdup.transcriptcount") {
              colnames(annotationNames)[1] = "ensemblTranscriptID"
    	      tt <- merge(annotationNames, printData, by.x="ensemblTranscriptID", by.y="tracking_id")
              tt <- tt[order(tt$pvalue), ]
              write.table(tt, file=paste("deseq2-allnorm/output/", one, "_vs_", two, ".", file, ".allnorm.txt", sep=""), row.names=F, sep="\t", quote=F)
              geneList <- as.vector(unlist(subset(tt, padj <= 0.05, select=c("ensemblTranscriptID"))))
              significantGenesFive <- c(significantGenesFive, geneList)
              geneList <- as.vector(unlist(subset(tt, padj <= 0.01, select=c("ensemblTranscriptID"))))
              significantGenesOne <- c(significantGenesOne, geneList)
 	    } 
          }
       }
    }

    pointSize <- 0.4

    ## make qc

    pdf(paste("deseq2-allnorm/qc/PCA.", file, ".0-05.pdf", sep=""), onefile=T)

    uniqueSignificantGenes <- unique(significantGenesFive)
    pcaDD <- subset(fpkmNormalisedCounts, rownames(rawCounts) %in% uniqueSignificantGenes)
    pFPKMS <- cbind(rownames(pcaDD), pcaDD)

    if (file == "markdup.genecount") {
      colnames(pFPKMS)[1] <- "ensemblGeneID"
    } else if (file == "markdup.transcriptcount") {
      colnames(pFPKMS)[1] <- "ensemblTranscriptID"
    }

    ## write data for clustering

    write.table(pFPKMS, file=paste("deseq2-allnorm/output/clusterFPKMS.",file,".0-05.txt",sep=""), row.names=F, sep="\t", quote=F)

    pca <- prcomp(t(pcaDD), center=TRUE, scale=TRUE)
    scores <- data.frame(targets$analysisID, pca$x[,1:2])
    pointSize <- 0.4
    print(qplot(x=PC1, y=PC2, main = file, data=scores, label=paste("         ", factor(targets$analysisID), sep=""), colour=factor(targets$sampleGroup))  + geom_text_repel(size=2.5) + scale_colour_discrete(name="sampleGroup") )
    dev.off()


    pdf(paste("deseq2-allnorm/qc/PCA.", file, ".0-01.pdf", sep=""), onefile=T)
    uniqueSignificantGenes <- unique(significantGenesOne)
    pcaDD <- subset(fpkmNormalisedCounts, rownames(rawCounts) %in% uniqueSignificantGenes)
    pFPKMS <- cbind(rownames(pcaDD), pcaDD)
    if(file == "markdup.genecount") {
      colnames(pFPKMS)[1] <- "ensemblGeneID"
    } else if(file == "markdup.transcriptcount") {
      colnames(pFPKMS)[1] <- "ensemblTranscriptID"
    }

    write.table(pFPKMS, file=paste("deseq2-allnorm/output/clusterFPKMS.",file,".0-01.txt",sep=""), row.names=F, sep="\t", quote=F)
    pca <- prcomp(t(pcaDD), center=TRUE, scale=TRUE)
    scores <- data.frame(targets$analysisID, pca$x[,1:2])
    pointSize <- 0.4
    print(qplot(x=PC1, y=PC2, main = file, data=scores, label=paste("         ", factor(targets$analysisID), sep=""), colour=factor(targets$sampleGroup))  + geom_text_repel(size=2.5) + scale_colour_discrete(name="sampleGroup"))
    dev.off()
}


