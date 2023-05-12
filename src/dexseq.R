#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(DEXSeq)

library(parallel)

library(BiocParallel)

library(stringr)

library(dplyr)

RunAnalysis <- function() {

	# gff file name
	gfffile = args[1]

	# gff file path
	gfffilePath = args[2]

	# targets filename with path 
	targetsFile = args[3]

	# Directory path for count files
	workDir = args[4]

	# Number of workers for BPPARAM
	wk = args[5]

	ffile <- paste(gfffilePath,"/",gfffile,sep="")

	# Read in targets file
	targets <- read.csv(targetsFile,header=T,check.names=F)

	# Substitute '-' with '_' in targets file columns
	targets$suppliedID = gsub("-", "_", targets$suppliedID)
	targets$analysisID = gsub("-", "_", targets$analysisID)

	# Read in analysisID and sampleGroup columns as character vectors
	sampleList <- as.character(targets$analysisID)
	sampleList1 <- as.character(targets$sampleGroup)

	# Split the targets file dataframe into separate lists based on sampleGroup values
	targetsGroup <- targets %>% group_split(sampleGroup)

	# Get unique sampleGroup values from targets file
	samplesU <- as.character(unique(targets$sampleGroup))

	# Generate of different combinations of two values from samplesU vector. Store these combinations in comparisonTable dataframe 
	comparisonTable <- data.frame(t(combn(samplesU,2)))

	# Join the combinations by '-vs-' and store it in myComparisons
	myComparisons1 <- paste(comparisonTable$X1,"-vs-",comparisonTable$X2,sep="")
	myComparisons <- myComparisons1[1]

	# Run the loop over all combinations of sampleGroup values
	for(i in myComparisons) {
    		# Split the string and store the values as list in c1
    		c1 <- unlist(strsplit(i, "-vs-"))
    		print(paste("Analysing ",c1[1]," vs ",c1[2],sep=""))
        	# Run a loop over all the elements of targetsGroup list 
		for(j in 1:length(targetsGroup)) {
			# if first element of c1 is present in any of the sampleGroup column inside targetsGroup lists, store the values of analysisID column from jth targetsGroup list into c11. If the value is not present, jump to next iteration  
			if(c1[1] %in% targetsGroup[[j]]$sampleGroup) {
				c11 <- as.character(targetsGroup[[j]]$analysisID)
				print(as.character(targetsGroup[[j]]$analysisID))
			} else {
				next
			}
		}

		for(j in 1:length(targetsGroup)) {
			# if second element of c1 is present in any of the sampleGroup column inside targetsGroup lists, store the values of analysisID column from jth targetsGroup list into c22. If the value is not present, jump to next iteration
			if(c1[2] %in% targetsGroup[[j]]$sampleGroup) {
				c22 <- as.character(targetsGroup[[j]]$analysisID)
				print(as.character(targetsGroup[[j]]$analysisID))
			} else {
				next
			}
		}
	
		# create a list of tr and wt values for c11 and c22 vectors. The c11class contains only 'tr' values over the length of c11. The c22class contains only 'wt' values over the length of c22.  
		c11class <- rep("tr",length(c11))
		c22class <- rep("wt",length(c22))
	
		# Combine both the c11class and c22class vectors and store the combined list as class vector
		class <- c(c11class,c22class)
	
		# Combine both the c11 and c22 vectors containing analysisID values and store it as s3 vector
		s3 = c(c11,c22)
	
		# Construct sampleTable dataframe with values of s3 (analysisID) as rownames and values of class as condition column
		sampleTable = data.frame(row.names = s3,condition = class)
        
        	print(paste("c11: ",c11," c22: ",c22,sep=""))

		# Replace '_' with '-' in c11 and c22 and store it in c111 and c222. This is only for finding the dexseqcount files inside the output directory.
        	c111 <- gsub("_", "-", c11)
        	c222 <- gsub("_", "-", c22)

		# Construct a cfiles list having full path to dexseqcount files inside sample names in output folder
        	cfiles = c(paste(workDir,"/",c111,"/",c111,".markdup.dexseqcount.txt",sep=""),paste(workDir,"/",c222,"/",c222,".markdup.dexseqcount.txt",sep=""))
	
		# Construct DEXSeqDataSet object from this data. The dxda1 object holds all the input data and it will be passed along the stages of the subsequent analysis. The third argument is a formula of the form ~ sample + exon + condition:exon that specifies the contrast with of a variable from the sample table columns and the ‘exon’ variable.
		dxda1 = DEXSeqDataSetFromHTSeq(
	     		cfiles,
	     		sampleData=sampleTable,
	     		design= ~ sample + exon + condition:exon,
	     		flattenedfile=ffile )

		# For parallelization of the script, MulticoreParam function is used from the BiocParallel package
		BPPARAM = MulticoreParam(workers = wk)

		# This function estimates the size factors using the "median ratio method" introduced in DESeq.
		dxda2 = estimateSizeFactors( dxda1 )

		print("Running estimate dispersion")	
		# This function obtains dispersion estimates for Negative Binomial distributed data
		dxda3 = estimateDispersions( dxda2, BPPARAM=BPPARAM )

		print("Running test for DEU")
		# The testForDEU tests for differential exon usage. It perform these tests for each exon in each gene.
		dxda4 = testForDEU( dxda3, BPPARAM=BPPARAM )

		print("Running estimate exon fold changes")
		# This function estimates relative exon usage fold changes. It is calculated based on the coefficients of a GLM that uses the formula: count ~ condition + exon + condition:exon. The resulting coefficients allow the estimation of changes on the usage of exons across different conditions.
		dxda5 = estimateExonFoldChanges( dxda4, fitExpToVar="condition", BPPARAM=BPPARAM )
	
		# In order to summarize the results without showing the values from intermediate steps, DEXSeqResults() function is called
		dxr1 = DEXSeqResults( dxda5 )
	
		print(dxr1)

		# Convert the dxr1 into dataframe for writing the data to file
		dxr2 <- data.frame(dxr1)

		# Sort the results by the lowest p-value
		dxr3 <- dxr2[order(dxr2$pvalue),]

		# Give the output filenames
		comparisonName = paste("dexseq/output/",c1[1],"_vs_",c1[2],".txt",sep="")
		# Write the data to file
		write.table(dxr3,quote=F,row.names=F,sep="\t",file=comparisonName)

	}
}

errorMessage <- try(RunAnalysis())

