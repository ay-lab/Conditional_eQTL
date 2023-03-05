#!/usr/bin/env Rscript

#===============================
# Conditional Matrix eQTL implementation
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# and the paper 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6075455/

# Sourya Bhattacharyya
# La Jolla Institute for Immunology
#===============================

library(MatrixEQTL)
library(data.table)
library(optparse)
library(plyr)
library(dplyr)
library(ggplot2)
require(gridExtra)
library(viridis)

options(scipen=10)
options(datatable.fread.datatable=FALSE)

# threhold of eQTL significance
FDR_THR_eQTL_SIG <- 0.05

# maximum no of iterations
MAX_ITERATION_COUNT <- 20

#===========================================================
option_list = list(

	make_option(c("--GenotypeFile"), type="character", default=NULL, help="Genotype file name (containing SNP ID and genotype information for different samples). Default = NULL. Mandatory parameter."),
	make_option(c("--SNPLocFile"), type="character", default=NULL, help="Genotype location file name (containing SNP ID and SNP coordinates). Default = NULL. Mandatory parameter."),
	make_option(c("--ExprFile"), type="character", default=NULL, help="Gene expression file (containing gene ID and gene expression for different samples). Default = NULL. Mandatory parameter."),
	make_option(c("--GeneLocFile"), type="character", default=NULL, help="Gene location file name (containing gene ID and gene coordinates). Default = NULL. Mandatory parameter."),
	
	make_option(c("--OutDir"), type="character", default=NULL, help="Output directory name to store all the results. Default = current working directory."),
	
	make_option(c("--CovariateFile"), type="character", default=NULL, help="File containing covariate information. If provided, covariates are taken into account for measuring QTL associations. Default = NULL, means no covariate file is provided."),

	# using absolute / log transformed gene expression
	make_option(c("--logGeneExpr"), type="integer", default=0, help="Integer value, if 1, signifies that input gene expression values are log2 transformed. Otherwise, absolute gene expression values are provided. Default = 0."),

	# previously this parameter's default value was 1
	make_option(c("--PValTrans"), type="numeric", default=1e-8, help="p-value threshold for trans-eQTLs. Default = 1e-8 (1e-8 is recommended in main package). This value was selected specifically for this application."),

	# previously this parameter's default value was 0
	make_option(c("--PValCis"), type="numeric", default=1e-2, help="p-value threshold for cis-eQTLs. Default = 1e-2 (1e-2 is recommended in main package). This value was selected specifically for this application."),

	make_option(c("--CisDist"), type="numeric", default=1e6, help="maximum distance at which gene-SNP pair is considered local. Default = 1e6 (i.e. 1 Mb). For iQTL, gene is replaced by interacting bin (interval) and so this distance may be varied to say, even 5000.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#=========
# parse input parameters
#=========
if (is.null(opt$CovariateFile)) {
	cat(sprintf("\n covariate file input is not provided - check the option --CovariateFile"))
	quit(save = "no")
}

ModelName <- 'modelLINEAR'
useModel <- modelLINEAR

SNP_file_name <- opt$GenotypeFile
snps_location_file_name <- opt$SNPLocFile
expression_file_name <- opt$ExprFile
gene_location_file_name <- opt$GeneLocFile

if (is.null(opt$CovariateFile)) {
	covariates_file_name <- NULL
} else {
	covariates_file_name <- opt$CovariateFile
}

# Only associations significant at this level will be saved
pvOutputThreshold_cis <- as.numeric(opt$PValCis)
pvOutputThreshold_tra <- as.numeric(opt$PValTrans)

# Distance for local gene-SNP pairs
cisDist = as.numeric(opt$CisDist)

OUTDIR <- opt$OutDir
system(paste("mkdir -p", OUTDIR))

#==================
# a few other internal parameters
#==================
# Building null distribution
nullDist <- FALSE

# Error covariance matrix
errorCovariance = numeric();

# Minimum MAF
minMAF = 0.05;

# Maximum missingness
maxMiss = 0.05;

# now decide the output file names
output_file_name_all <- paste0(OUTDIR, '/out_all_eqtl.bed')
output_file_name_cis <- paste0(OUTDIR, '/out_cis_filt.bed')
output_file_name_tra <- paste0(OUTDIR, '/out_trans_filt.bed')
output_file_name_qqplot <- paste0(OUTDIR, '/out_qqplot.png')
output_file_name_all_cis <- paste0(OUTDIR, '/out_cis_all.bed')

#================
# now write this parameters in a output file
#================
Parameter_File <- paste0(OUTDIR, '/Parameters.txt')
fp_Param <- file(Parameter_File, "w")
outtext <- paste0("\n *** Parameters involved in Matrix EQTL current execution *** \n GenotypeFile : ", SNP_file_name, "\n Gene expression (or loop) file name: ", expression_file_name, "\n Base output directory specified : ", OUTDIR, "\n Covariares file (NULL means not used): ", covariates_file_name, "\n ModelName : ", ModelName, "\n p-value threshold for cis-eQTLs: ", pvOutputThreshold_cis, "\n p-value threshold for trans-eQTLs: ", pvOutputThreshold_tra, "\n Using identity matrix as covariance \n Using 0.05 as MAF threshold \n Using 0.05 as maximum missingness threshold ")
writeLines(outtext, con=fp_Param, sep="\n")
close(fp_Param)

# also define an output log file
logfile <- paste0(OUTDIR, '/Complete_Output_', gsub("-","_",gsub(" ","_",format(Sys.time(), "%F %H-%M"))), '.log')
outtext <- paste0("\n ********* Complete output log ****** \n")
cat(outtext, file=logfile, append=FALSE, sep="\n")

#================
## Load gene expression data
cat(sprintf("\n before loading gene expression data \n"))
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load genotype data
cat(sprintf("\n before loading genotype data \n"))
if(nrow(fread(SNP_file_name, header = TRUE, nrows = 10)) == 0) {
	# Touch the output file and exit if there are no suitable snps
	try(system(paste("touch", output_file_name_cis), intern = TRUE, ignore.stderr = TRUE))
	try(system(paste("touch", output_file_name_tra), intern = TRUE, ignore.stderr = TRUE))
	try(system(paste("touch", output_file_name_qqplot), intern = TRUE, ignore.stderr = TRUE))
	try(system(paste("touch", output_file_name_all_cis), intern = TRUE, ignore.stderr = TRUE))
	quit(save = "no")
}

snps = SlicedData$new();
snps$fileDelimiter = " ";	#"\t"      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000    	# read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load covariates
if (!is.null(opt$CovariateFile)) {
	cov <- data.table::fread(covariates_file_name)
	# covMat <- as.matrix(t(cov[, .(PC1,PC2)]))
	tempDF <- cbind.data.frame(cov$PC1, cov$PC2)
	covMat <- as.matrix(t(tempDF))
	# colnames(covMat) <- cov[, IID]
	colnames(covMat) <- as.vector(cov$IID)
}

if (!is.null(opt$CovariateFile)) {
	target_colnames <- intersect(intersect(gene$columnNames, snps$columnNames), as.vector(cov$IID))
	cat(sprintf("\n\n ===>> number of donors in the gene expression file : %s \n ===>> number of donors in the genotype file : %s  \n ===>> number of donors in the covariate file : %s \n ==>> final number of common donors : %s ", length(gene$columnNames), length(snps$columnNames), length(as.vector(cov$IID)), length(target_colnames)))
	cat(sprintf("\n\n target_colnames : %s ", paste(target_colnames, collapse=" ")))

	# now subset snps, gene and covariate information using these columns
	cat(sprintf("\n before re-ordering SNP columns \n"))
	cols <- sapply(target_colnames, function(x) which(x == snps$columnNames))
	cat(sprintf("\n\n cols : %s ", paste(cols, collapse=" ")))
	cat(sprintf("\n number of entries in cols : %s ", length(cols)))
	snps$ColumnSubsample(cols)

	cat(sprintf("\n before re-ordering gene columns \n"))		
	cols <- c()
	for (x in target_colnames) {
		cols <- c(cols, which(x == gene$columnNames)[1])
	}
	cat(sprintf("\n\n cols : %s ", paste(cols, collapse=" ")))
	cat(sprintf("\n number of entries in cols : %s ", length(cols)))
	cat(sprintf("\n\n starting ColumnSubsample for gene"))
	gene$ColumnSubsample(cols)

	cat(sprintf("\n selecting covariate information \n"))
	covMat <- covMat[, as.vector(as.character(target_colnames))]
	cvrt = SlicedData$new()
	cvrt$CreateFromMatrix(covMat)

} else {
	cat(sprintf("\n before re-ordering gene columns \n"))
	cols <- sapply(gene$columnNames, function(x) which(x == snps$columnNames))
	snps$ColumnSubsample(cols)
}

## Permutate samples if we are calculating the null distribution
if (nullDist) {
	random <- sample(1:snps$nCols(), size = snps$nCols());
	snps$ColumnSubsample(random);
}

## Filter SNPs with high missingness
miss.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
	slice = snps[[sl]];
	miss.list[[sl]] = unlist(apply(slice, 1,
		function(x) sum(is.na(x))))/ncol(slice);
}
miss = unlist(miss.list)

outtext <- paste0("\n SNPs before missingness filtering: ", nrow(snps))
cat(outtext, file=logfile, append=TRUE, sep="\n")
snps$RowReorder(miss < maxMiss);
outtext <- paste0("\n SNPs after missingness filtering: ", nrow(snps))
cat(outtext, file=logfile, append=TRUE, sep="\n")

## Filter SNPs with low MAF
maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
	slice = snps[[sl]];
	maf.list[[sl]] = rowMeans(slice, na.rm=TRUE)/2;
	maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

## Look at the distribution of MAF
# hist(maf[maf<0.1],seq(0,0.1,length.out=100))
outtext <- paste0("\n SNPs before MAF filtering: ", nrow(snps))
cat(outtext, file=logfile, append=TRUE, sep="\n")
if(any(maf > minMAF)) {
	snps$RowReorder(maf>minMAF);
	outtext <- paste0("\n SNPs after MAF filtering: ", nrow(snps))
	cat(outtext, file=logfile, append=TRUE, sep="\n")	
} else {
	# Touch the output file and exit if there are no suitable snps
	try(system(paste("touch", output_file_name_cis), intern = TRUE, ignore.stderr = TRUE))
	try(system(paste("touch", output_file_name_tra), intern = TRUE, ignore.stderr = TRUE))
	try(system(paste("touch", output_file_name_qqplot), intern = TRUE, ignore.stderr = TRUE))
	try(system(paste("touch", output_file_name_all_cis), intern = TRUE, ignore.stderr = TRUE))
	quit(save = "no")
}

# reading gene and SNP coordinates
cat(sprintf("\n before reading SNP location file \n"))
snpspos <- data.table::fread(snps_location_file_name, header = T)
# snpspos <- read.table(snps_location_file_name, header = T, stringsAsFactors = FALSE);
colnames(snpspos) <- c('SNP', 'snpchr', 'snppos')
cat(sprintf("\n before reading gene location file \n"))
genepos <- data.table::fread(gene_location_file_name, header = T);
# genepos <- read.table(gene_location_file_name, header = T, stringsAsFactors = FALSE);
colnames(genepos) <- c('gene', 'genechr', 'genestart', 'geneend')

# store the raw gene expression data into a matrix
# (before quantile normalization)
# all the residual gene expression computation would be done in this matrix
Orig_GeneExpr_Mat <- as.matrix(gene)

#======================================
# iteration to find eQTL
#======================================
itn_count <- 0

# also convert the SNPs (genotype data) into a matrix
z <- as.matrix(snps)

while(1) {

	itn_count <- itn_count + 1
	outtext <- paste0("\n\n\n  ==> now calling main matrix QTL function - iteration count : ", itn_count, "\n\n\n")
	cat(outtext, file=logfile, append=TRUE, sep="\n")	

	curroutdir <- paste0(OUTDIR, '/iteration_', itn_count)
	system(paste("mkdir -p", curroutdir))	

	output_file_name_all <- paste0(curroutdir, '/out_all_eqtl.bed')
	output_file_name_cis <- paste0(curroutdir, '/out_cis_eqtl.bed')
	output_file_name_cis_detailed <- paste0(curroutdir, '/out_cis_eqtl_detailed.bed')

	# Quantile normalization of the gene expression values
	for( sl in 1:length(gene) ) {
	  mat = gene[[sl]];
	  mat = t(apply(mat, 1, rank, ties.method = "average"));
	  mat = qnorm(mat / (ncol(gene)+1));
	  gene[[sl]] = mat;
	}
	rm(sl, mat);

	if ((file.exists(output_file_name_cis_detailed) == FALSE) | (file.exists(output_file_name_all) == FALSE) | (file.exists(output_file_name_cis) == FALSE)) {

		# covariates is defined
		me <- Matrix_eQTL_main(snps = snps, gene = gene, cvrt = cvrt, output_file_name = output_file_name_all, pvOutputThreshold = pvOutputThreshold_tra, useModel = useModel, errorCovariance = errorCovariance, verbose = TRUE, output_file_name.cis = output_file_name_cis, pvOutputThreshold.cis = pvOutputThreshold_cis, snpspos = snpspos, genepos = genepos, cisDist = cisDist, pvalue.hist = "qqplot", min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE);

		## Results:
		outtext <- paste0("\n Analysis done in: ", me$time.in.sec, " seconds")
		cat(outtext, file=logfile, append=TRUE, sep="\n")	

		# now read this all eQTL file
		# 1st column stores the SNP: column names = snps
		# 2nd column stores the gene: column names = gene
		#alleqtldata <- data.table::fread(output_file_name_all, header=T)
		alleqtldata <- data.table::fread(output_file_name_cis, header=T)
		mergeDF <- dplyr::inner_join(alleqtldata, snpspos) %>% dplyr::inner_join(genepos)
		write.table(mergeDF, output_file_name_cis_detailed, row.names=F, col.names=T, quote=F, sep="\t", append=F)

	} else {		
		# eQTL for this iteration is already computed - read the computed file
		outtext <- paste0("\n\n ===>> iteration : ", itn_count, " eQTL is already computed !!! ")
		cat(outtext, file=logfile, append=TRUE, sep="\n")	
		alleqtldata <- data.table::fread(output_file_name_cis, header=T)
		mergeDF <- data.table::fread(output_file_name_cis_detailed, header=T)
	}

	# for the first iteration, get the significant gene list
	# which have at least one significantly associated SNP (with respect to the given FDR threshold)
	# only this set of genes will be used in the subsequent iterations
	if (itn_count == 1) {
		Significant_GeneList <- as.vector(unique(mergeDF[which(mergeDF[, 6] < FDR_THR_eQTL_SIG), 2]))
		outtext <- paste0("\n\n ===>> iteration : ", itn_count, " number of genes with significant QTL association: ", length(Significant_GeneList))
		cat(outtext, file=logfile, append=TRUE, sep="\n")		
		
		# if no significant gene-SNP association at the first iteration, we break from the loop
		if (length(Significant_GeneList) == 0) {			
			break
		}

		# also check if any gene has maximum expression < 1 TPM for any genotype
		# discard those genes
		# assuming gene expression is provided in terms of log2 expression
		discard_gene_list <- c()
		for (i in 1:length(Significant_GeneList)) {
			currgene <- Significant_GeneList[i]
			r <- which(rownames(Orig_GeneExpr_Mat) == currgene)
			yvec <- as.numeric(Orig_GeneExpr_Mat[r, 1:ncol(Orig_GeneExpr_Mat)])
			yvec <- yvec[!is.na(yvec)]
			if (max(yvec) < 1) {
				discard_gene_list <- c(discard_gene_list, currgene)
			}
		}
		if (length(discard_gene_list) == length(Significant_GeneList)) {
			# all genes are discarded
			break
		}
		Significant_GeneList <- setdiff(Significant_GeneList, discard_gene_list)
		outtext <- paste0("\n\n ===>> iteration : ", itn_count, " number of genes with significant QTL association (after filtering out very low expression genes - outliers: ", length(Significant_GeneList))
		cat(outtext, file=logfile, append=TRUE, sep="\n")		

		# also find out the reference p-value threshold corresponding to the FDR < FDR_THR_eQTL_SIG
		# as suggested in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5993513/
		REF_PVAL_THR <- max(as.numeric(mergeDF[which(mergeDF[, 6] < FDR_THR_eQTL_SIG), 5]))
		outtext <- paste0("\n\n ===>> ALSO FINALIZING THE P-value threshold : ", REF_PVAL_THR, "\n\n")
		cat(outtext, file=logfile, append=TRUE, sep="\n")

		# initialize a separate vector for the first iteration
		# this vector will be used
		Target_Significant_GeneList <- Significant_GeneList
	}

	# for the 2nd iteration onwards, check the genes with p-value threshold
	# and also check the low gene expression
	if (itn_count > 1) {
		Target_Significant_GeneList <- intersect(as.vector(unique(mergeDF[which(mergeDF[, 5] < REF_PVAL_THR), 2])), Significant_GeneList)
		outtext <- paste0("\n\n ===>> iteration : ", itn_count, " number of genes with significant QTL association: ", length(Target_Significant_GeneList))
		cat(outtext, file=logfile, append=TRUE, sep="\n")
		
		# if no significant gene-SNP association at the first iteration, we break from the loop
		if (length(Target_Significant_GeneList) == 0) {			
			break
		}

		# also check if any gene has maximum expression < 1 TPM for any genotype
		# discard those genes
		# assuming gene expression is provided in terms of log2 expression
		discard_gene_list <- c()
		for (i in 1:length(Target_Significant_GeneList)) {
			currgene <- Target_Significant_GeneList[i]
			r <- which(rownames(Orig_GeneExpr_Mat) == currgene)
			yvec <- as.numeric(Orig_GeneExpr_Mat[r, 1:ncol(Orig_GeneExpr_Mat)])
			yvec <- yvec[!is.na(yvec)]
			if (max(yvec) < 1) {
				discard_gene_list <- c(discard_gene_list, currgene)
			}
		}
		if (length(discard_gene_list) == length(Target_Significant_GeneList)) {
			# all genes are discarded
			break
		}
		Target_Significant_GeneList <- setdiff(Target_Significant_GeneList, discard_gene_list)
		outtext <- paste0("\n\n ===>> iteration : ", itn_count, " number of genes with significant QTL association (after filtering out very low expression genes - outliers: ", length(Target_Significant_GeneList))
		cat(outtext, file=logfile, append=TRUE, sep="\n")
	}

	# sanity check - may not be required
	# for the second iteration onwards, discard the entries from "mergeDF"
	# (the gene-SNP pairs) which are already present in the "Complete_SignGeneSNP_DF"
	if (itn_count > 1) {
		curroutdata <- mergeDF[, 1:2]
		refdata <- Complete_SignGeneSNP_DF[,1:2]
		colnames(curroutdata) <- c('snp', 'gene')
		colnames(refdata) <- c('snp', 'gene')
		m <- plyr::match_df(curroutdata, refdata)
		Match_RowIdx <- as.integer(rownames(m))
		if (length(Match_RowIdx) > 0) {
			outtext <- paste0("\n\n **** DISCARD ENTRIES -- number of entries in the significant gene-SNP pair of mergeDF -- ", nrow(curroutdata), "\n number of entries / gene-SNP pairs already in the final output (duplicates) : ", length(Match_RowIdx))
			cat(outtext, file=logfile, append=TRUE, sep="\n")
			# filtered set of entries in the merged data frame
			idx <- setdiff(seq(1, nrow(curroutdata)), Match_RowIdx)
			mergeDF <- mergeDF[idx, ]
		}
	}	

	#==================
	# boolean variable indicating construction of data frame
	bool_DF <- FALSE
	
	# then scan through individual genes	
	for (i in 1:length(Target_Significant_GeneList)) {
		currgene <- Target_Significant_GeneList[i]
		if (itn_count == 1) {			
			# for the first iteration, get the gene-SNP with lowest FDR
			# for multiple mminimum FDR, take the first entry - it'll have the lowest p-value
			minFDR <- min(mergeDF[which(mergeDF[, 2] == currgene), 6])	
			outtext <- paste0("\n itn_count : ", itn_count, " processing gene : ", currgene, "  minimum FDR : ", minFDR)
			cat(outtext, file=logfile, append=TRUE, sep="\n")			
			rowidx <- which((mergeDF[, 2] == currgene) & (mergeDF[, 6] == minFDR))[1]
			beta_val <- mergeDF[rowidx, 3]
			currSNP <- mergeDF[rowidx, 1]
			outtext <- paste0("  --- matching row : ", rowidx, " significantly associated SNP : ", currSNP, "  beta val : ", beta_val)
			cat(outtext, file=logfile, append=TRUE, sep="\n")
		} else {
			# follow https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5993513/
			# for subsequent iterations, use the REF_PVAL_THR for finding significant association
			rowidxset <- which((mergeDF[, 2] == currgene) & (mergeDF[, 5] < REF_PVAL_THR))
			if (length(rowidxset) > 0) {
				# get the most significant SNP for this gene - lowest p-value
				min_p_val <- min(mergeDF[rowidxset, 5])
				rowidx <- which((mergeDF[, 2] == currgene) & (mergeDF[, 5] == min_p_val))[1]
				beta_val <- mergeDF[rowidx, 3]
				currSNP <- mergeDF[rowidx, 1]
				outtext <- paste0("  --- matching row : ", rowidx, " p-value of the current associated SNP : ", mergeDF[rowidx, 5], "  reference p-value threshold (global) : ", REF_PVAL_THR, "  SNP ID :  ", currSNP, "  beta val : ", beta_val)
				cat(outtext, file=logfile, append=TRUE, sep="\n")
			} else {
				rowidx <- 0
				outtext <- paste0(" --- no significantly associated SNP")
				cat(outtext, file=logfile, append=TRUE, sep="\n")
			}
		}		
			
		# if the current gene has at least one significant SNP, rowidx > 0: otherwise 0
		if (rowidx > 0) {
			# store the significant gene-SNP pair information in a data frame
			# called "currItn_SignGeneSNP_DF"			
			if (bool_DF == FALSE) {
				currItn_SignGeneSNP_DF <- mergeDF[rowidx, ]
				bool_DF <- TRUE
			} else {
				currItn_SignGeneSNP_DF <- rbind.data.frame(currItn_SignGeneSNP_DF, mergeDF[rowidx, ])
			}

			#===================
			# most important step
			# here we perform residual gene expression, by minimizing the effect for the current significant SNP
			# moticated by the paper
			# https://academic.oup.com/hmg/article/26/8/1444/2970473
			#===================
			r <- which(rownames(Orig_GeneExpr_Mat) == currgene)
			outtext <- paste0("\n\n --->> updating the gene expression matrix ---- \n the gene : ", currgene, "  is placed in row : ", r, " of the gene expression matrix ")
			cat(outtext, file=logfile, append=TRUE, sep="\n")
			r1 <- which(rownames(z) == currSNP)
			outtext <- paste0("\n ---- the associated SNP : ", currSNP, "  is placed in row : ", r1, " of the reference genotype data")
			cat(outtext, file=logfile, append=TRUE, sep="\n")
			# first get the gene expression for all donors for this gene
			# Note: there can be donors for which gene expression may be NA
			# so filter out those entries
			# and then compute the standard deviation of the remaining gene expression values
			yvec <- as.numeric(Orig_GeneExpr_Mat[r, 1:ncol(Orig_GeneExpr_Mat)])
			outtext <- paste0("\n number of donors (gene expression data) : ", length(yvec))
			cat(outtext, file=logfile, append=TRUE, sep="\n")
			yvec <- yvec[!is.na(yvec)]
			outtext <- paste0(" ---- number of donors with non-NA gene expression : ", length(yvec))
			cat(outtext, file=logfile, append=TRUE, sep="\n")
			mean_yvec <- mean(yvec)
			sd_yvec <- sd(yvec)
			outtext <- paste0(" ---- mean of gene expression : ", mean_yvec, "  standard deviation of gene expression : ",  sd_yvec)
			cat(outtext, file=logfile, append=TRUE, sep="\n")
			
			# the residual gene expression is computed as
			# original_gene_expr - beta * stdev * genotype_val
			# where genotype_val is 0, 1, or 2 for different donors
			# check https://en.wikipedia.org/wiki/Effect_size
			genotype_currSNP <- as.integer(as.numeric(z[r1, ]))
			orig_gene_expr_values_vec <- as.numeric(Orig_GeneExpr_Mat[r, ])
			for (j in 1:ncol(Orig_GeneExpr_Mat)) {
				if ((!is.na(Orig_GeneExpr_Mat[r, j])) & (!is.na(z[r1, j]))) {
					Orig_GeneExpr_Mat[r, j] <- (Orig_GeneExpr_Mat[r, j] - beta_val * sd_yvec * z[r1, j])
				}
			}
			modified_gene_expr_values_vec <- as.numeric(Orig_GeneExpr_Mat[r, ])

			# plot the genotype specific variation of gene expression for this SNP 
			# before and after this adjustment
			bool_plot_geno <- FALSE
			for (genotype_val in c(0,1,2)) {			
				idx <- which(genotype_currSNP == genotype_val)
				outtext <- paste0("\n genotype_val : ", genotype_val, "  number of donors with this genotype : ", length(idx))
				cat(outtext, file=logfile, append=TRUE, sep="\n")
				if (length(idx) > 0) {
					prevexprvec <- orig_gene_expr_values_vec[idx]
					currexprvec <- modified_gene_expr_values_vec[idx]
					if (opt$logGeneExpr == 1) {
						# applicable only if input gene expression values are log2 transformed
						prevexprvec_abs <- (2^orig_gene_expr_values_vec[idx])
						currexprvec_abs <- (2^modified_gene_expr_values_vec[idx])
					}
					if (bool_plot_geno == FALSE) {
						prevexprvec_DF <- data.frame(Genotype=paste0(genotype_val, ' (', length(idx), ')'), expr=prevexprvec)
						currexprvec_DF <- data.frame(Genotype=paste0(genotype_val, ' (', length(idx), ')'), expr=currexprvec)
						if (opt$logGeneExpr == 1) {
							# applicable only if input gene expression values are log2 transformed
							prevexprvec_abs_DF <- data.frame(Genotype=paste0(genotype_val, ' (', length(idx), ')'), expr=prevexprvec_abs)
							currexprvec_abs_DF <- data.frame(Genotype=paste0(genotype_val, ' (', length(idx), ')'), expr=currexprvec_abs)
						}
						bool_plot_geno <- TRUE
					} else {
						prevexprvec_DF <- rbind.data.frame(prevexprvec_DF, data.frame(Genotype=paste0(genotype_val, ' (', length(idx), ')'), expr=prevexprvec))
						currexprvec_DF <- rbind.data.frame(currexprvec_DF, data.frame(Genotype=paste0(genotype_val, ' (', length(idx), ')'), expr=currexprvec))
						if (opt$logGeneExpr == 1) {
							# applicable only if input gene expression values are log2 transformed		
							prevexprvec_abs_DF <- rbind.data.frame(prevexprvec_abs_DF, data.frame(Genotype=paste0(genotype_val, ' (', length(idx), ')'), expr=prevexprvec_abs))
							currexprvec_abs_DF <- rbind.data.frame(currexprvec_abs_DF, data.frame(Genotype=paste0(genotype_val, ' (', length(idx), ')'), expr=currexprvec_abs))
						}
					}
				}
			}	# end genotype loop

			# plot these gene expression
			plotfile <- paste0(curroutdir, '/plot_variation_geneExpr_', currgene, '_temp.pdf')
			currplot1 <- ggplot2::ggplot(prevexprvec_DF, aes(x=Genotype, y=expr, fill=Genotype)) + geom_boxplot() + scale_fill_viridis(discrete = TRUE, alpha=0.6) + geom_jitter(position=position_jitter(width=.1, height=0), color="black", size=0.4, alpha=0.9) + labs(title=paste0(currgene, "\n", currSNP, "\n beta : ", beta_val, " BEFORE"), x="Genotype", y="GeneExpr (log2)") + theme(panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=16), legend.position="none")
			currplot2 <- ggplot2::ggplot(currexprvec_DF, aes(x=Genotype, y=expr, fill=Genotype)) + geom_boxplot() + scale_fill_viridis(discrete = TRUE, alpha=0.6) + geom_jitter(position=position_jitter(width=.1, height=0), color="black", size=0.4, alpha=0.9) + labs(title=paste0(currgene, "\n", currSNP, "\n AFTER MOD"), x="Genotype", y="GeneExpr (log2)") + theme(panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=16), legend.position="none")	
			if (opt$logGeneExpr == 1) {
				# applicable only if input gene expression values are log2 transformed	
				currplot3 <- ggplot2::ggplot(prevexprvec_abs_DF, aes(x=Genotype, y=expr, fill=Genotype)) + geom_boxplot() + scale_fill_viridis(discrete = TRUE, alpha=0.6) + geom_jitter(position=position_jitter(width=.1, height=0), color="black", size=0.4, alpha=0.9) + labs(title=paste0(currgene, "\n", currSNP, "\n beta : ", beta_val, "\n BEFORE"), x="Genotype", y="GeneExpr") + theme(panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=16), legend.position="none")
				currplot4 <- ggplot2::ggplot(currexprvec_abs_DF, aes(x=Genotype, y=expr, fill=Genotype)) + geom_boxplot() + scale_fill_viridis(discrete = TRUE, alpha=0.6) + geom_jitter(position=position_jitter(width=.1, height=0), color="black", size=0.4, alpha=0.9) + labs(title=paste0(currgene, "\n", currSNP, "\n AFTER MOD"), x="Genotype", y="GeneExpr") + theme(panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=16), legend.position="none")
			}

			if (opt$logGeneExpr == 1) {
				# applicable only if input gene expression values are log2 transformed	
				ggsave(plotfile, arrangeGrob(currplot1, currplot2, currplot3, currplot4), width=10, height=8)
			} else {
				ggsave(plotfile, arrangeGrob(currplot1, currplot2), width=10, height=8)
			}

			#===================
			# end update and validate the change in gene expression
			#===================

		}	# end presence of significant SNP condition

	}	# end significant gene loop

	if (bool_DF == FALSE) {
		# no significant SNP at the current iteration
		break
	}

	# otherwise, one or more genes had changed expression
	# so create a combined plot file
	outplotfile <- paste0(curroutdir, '/plot_variation_geneExpr.pdf')
	system(paste("pdfunite", paste0(curroutdir, '/plot_variation_geneExpr_*.pdf'), outplotfile))
	system(paste("rm", paste0(curroutdir, '/plot_variation_geneExpr_*.pdf')))

	# we also dump the significant gene-SNP pairs
	# along with the iteration number - append one column named "itn"
	currItn_SignGeneSNP_DF$itn <- rep(itn_count, nrow(currItn_SignGeneSNP_DF))

	# write the significant gene-SNP pair for the current iteration
	write.table(currItn_SignGeneSNP_DF, paste0(curroutdir, '/significant_gene_SNP_pairs.bed'), row.names=F, col.names=T, quote=F, sep="\t", append=F)

	if (itn_count == 1) {		
		Complete_SignGeneSNP_DF <- currItn_SignGeneSNP_DF
	} else {	
		Complete_SignGeneSNP_DF <- rbind.data.frame(Complete_SignGeneSNP_DF, currItn_SignGeneSNP_DF)
	}

	#===========================
	# write back the updated gene expression data for the next iteration
	# keep row.names=T option
	#===========================
	modified_gene_expr_file <- paste0(curroutdir, '/modified_gene_expr.bed')
	write.table(Orig_GeneExpr_Mat, modified_gene_expr_file, row.names=T, col.names=T, quote=F, sep="\t", append=F)

	# again assign the converted gene expression data back to the genes structure, for the next iteration
	gene$Clear()
	# then load from the matrix
	# row and column names will be copied
	gene$CreateFromMatrix( Orig_GeneExpr_Mat );

	# condition - sourya
	# continue upto MAX_ITERATION_COUNT iterations
	if (itn_count == MAX_ITERATION_COUNT) {
		break
	}

}	# end while loop - iteration QTL

# write the complete set of significant gene-SNP pair for the current iteration
if (exists("Complete_SignGeneSNP_DF")) {
	write.table(Complete_SignGeneSNP_DF, paste0(OUTDIR, '/complete_set_significant_gene_SNP_pairs.bed'), row.names=F, col.names=T, quote=F, sep="\t", append=F)
}


