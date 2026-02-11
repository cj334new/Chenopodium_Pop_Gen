# Chenopodium_Pop_Gen

## Overview

This project investigates genomic diversity, introgression, selective
sweeps, demographic history, and regulatory variation (eQTLs) across
*Chenopodium* accessions spanning 20 species, with a focus on quinoa
environmental adaptation.

------------------------------------------------------------------------

## Directory Structure

    001.selection/           # Selective sweep analysis (π ratio, FST)
    002.IBD/                 # Identity-by-descent analysis
    003.introgression/       # ABBA-BABA (fd statistics) introgression analysis
    004.population_history/  # MSMC demographic inference
    005.eQTL/                # eQTL mapping and LD block detection
    README.md

------------------------------------------------------------------------

## 1. Nucleotide Diversity and FST
-   Window size: 50 kb
-   Step size: 5 kb

Example:

``` bash
vcftools --vcf input.vcf          --window-pi 50000          --window-pi-step 5000          --out pi_output

vcftools --vcf input.vcf          --weir-fst-pop group1.txt          --weir-fst-pop group2.txt          --fst-window-size 50000          --fst-window-step 5000          --out fst_output
```

------------------------------------------------------------------------

## 2. Selective Sweep Identification

Selective sweeps were identified using:

-   π ratio (π_group1 / π_group2)
-   FST (top 5% threshold)
-   Window size: 50 kb
-   Step size: 5 kb

Regions in the top 5% of both statistics were considered candidate
selective sweeps.

------------------------------------------------------------------------


## 3. Demographic History (MSMC)

-   Tool: **MSMC**
-   Mutation rate: 7e-9 per bp per generation
-   Generation time: 1 year

Effective population size (Ne) was scaled using:

Ne = 1 / (2 × mutation rate × coalescent rate)

------------------------------------------------------------------------

## 4. Introgression Analysis (ABBA-BABA)

-   Tool: genomics_general (fd statistics)
-   Window size: 1 Mb
-   Step size: 100 kb
-   Minimum SNPs per window: 100

``` bash
python ABBABABA.py        --window 1000000        --step 100000        --minSites 100
```

------------------------------------------------------------------------

## 5. eQTL Mapping

-   Tool: **MatrixEQTL**
-   Covariates: First 3 genotype PCs (Plink)
-   Expression normalization: `qqnorm()` in R
-   Bonferroni threshold: 1.90e-13

### cis/trans Definition

-   cis-eQTL: within 30 kb of gene TSS or TTS
-   trans-eQTL: outside 30 kb

------------------------------------------------------------------------

------------------------------------------------------------------------

## Keywords

Chenopodium · Quinoa · Population genomics · eQTL · Selective sweep ·
Introgression · Polyploid · Environmental adaptation · Flowering time

------------------------------------------------------------------------