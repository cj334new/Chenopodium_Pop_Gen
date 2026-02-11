library("MatrixEQTL")
base.dir = ('eQTL_work_PATH'); # change path

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(base.dir, "/hapmap_for_eQTL.txt", sep="");
snps_location_file_name = paste(base.dir, "/snpsloc.txt", sep="");

# Gene expression file name
expression_file_name = paste(base.dir, "/quinoa_rnaseq_230.fpkm", sep="");
gene_location_file_name = paste(base.dir, "/geneloc.txt", sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste(base.dir, "/Covariates.txt", sep="");
#covariates_file_name = character();

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-5;
pvOutputThreshold_tra = 1e-5;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

# Distance for local gene-SNP pairs
cisDist = 3e4;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
#cat('Detected local eQTLs:', '\n');
#show(me$cis$eqtls)
write.table(me$cis$eqtls, file="local_eQTLs.txt", quote=FALSE, sep = "\t")
#cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)
write.table(me$trans$eqtls, file="distant_eQTLs.txt", quote=FALSE, sep = "\t")

## Plot the Q-Q plot of local and distant p-values

plot(me)
