# MatrixEQTL cis/trans eQTL Analysis Pipeline

This document describes a standardized, publication-ready pipeline for
performing **cis- and trans-eQTL analysis** using the `MatrixEQTL` R
package.

The pipeline integrates genotype, gene expression, and covariate data
and applies Bonferroni correction to determine statistically significant
associations.

------------------------------------------------------------------------

## Overview

-   **MatrixEQTL framework**\
    Efficient linear modeling approach for large-scale eQTL mapping.

-   **Population structure correction**\
    The first three principal components (PC1--PC3), derived from
    genotype data using PLINK, are included as covariates.

-   **Cis vs Trans definition**

    -   **cis-eQTL**: SNP located within 30 kb of gene TSS or TTS\
    -   **trans-eQTL**: SNP located outside 30 kb window

------------------------------------------------------------------------

## Requirements

### Input Files

-   `hapmap_for_eQTL.txt`\
    Genotype matrix (SNP × individuals)

-   `snpsloc.txt`\
    SNP genomic positions (chr, position)

-   `quinoa_rnaseq_230.fpkm`\
    Gene expression matrix (genes × individuals)

-   `geneloc.txt`\
    Gene genomic coordinates (chr, start, end)

-   `Covariates.txt`\
    Covariate file (e.g., PC1, PC2, PC3)

------------------------------------------------------------------------

## Preprocessing Steps

### 1. Sample Selection

-   230 tetraploid accessions selected from 558 accessions with SNP
    data.

### 2. Gene Filtering

-   Genes with **median FPKM = 0** were excluded.

### 3. Expression Normalization

Gene expression values were normalized using:

``` r
qqnorm(expression_values)
```

------------------------------------------------------------------------

## Statistical Framework

-   Linear model: `modelLINEAR`
-   Error covariance: Identity matrix
-   cis window size: **30 kb (3e4 bp)**
-   Significance threshold: **Bonferroni correction (α = 0.05)**

Bonferroni-corrected P-value threshold:

    1.90e-13
    -log10(P) = 12.72

------------------------------------------------------------------------

## Pipeline

### 1. Load Required Library

``` r
library("MatrixEQTL")
```

------------------------------------------------------------------------

### 2. Set Parameters

``` r
useModel = modelLINEAR
cisDist = 3e4
pvOutputThreshold_cis = 1e-5
pvOutputThreshold_tra = 1e-5
```

------------------------------------------------------------------------

### 3. Load Genotype Data

``` r
snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$LoadFile("hapmap_for_eQTL.txt")
```

------------------------------------------------------------------------

### 4. Load Expression Data

``` r
gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$LoadFile("quinoa_rnaseq_230.fpkm")
```

------------------------------------------------------------------------

### 5. Load Covariates

``` r
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$LoadFile("Covariates.txt")
```

------------------------------------------------------------------------

### 6. Run MatrixEQTL

``` r
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = "trans_output.txt",
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  output_file_name.cis = "cis_output.txt",
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot"
)
```

------------------------------------------------------------------------

## Output Files

-   `local_eQTLs.txt`\
    Significant cis-eQTL associations

-   `distant_eQTLs.txt`\
    Significant trans-eQTL associations

------------------------------------------------------------------------

## Notes

-   Bonferroni correction is conservative but reduces false positives.
-   Population stratification correction is essential to avoid spurious
    associations.
-   The 30 kb cis window is commonly used in plant eQTL studies but may
    be adjusted depending on LD decay.
-   Q-Q plots can be used to evaluate inflation of test statistics.

------------------------------------------------------------------------