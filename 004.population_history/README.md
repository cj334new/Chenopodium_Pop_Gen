# MSMC2 Population History Analysis Pipeline

This document describes a standardized, publication-ready pipeline for
performing **MSMC2 (Multiple Sequentially Markovian Coalescent)**
analysis to infer historical effective population size (Ne) trajectories
from phased genome data.

The pipeline includes: - Individual-specific phased VCF extraction -
Mask generation based on sequencing depth - MSMC input file
construction - MSMC2 execution - Downstream visualization in R

------------------------------------------------------------------------

## Overview

**MSMC2** estimates changes in effective population size over time using
coalescent theory and phased whole-genome variation data.

Typical applications:

-   Inferring demographic history
-   Detecting bottlenecks or expansions
-   Comparing population size dynamics
-   Complementing selection sweep and introgression analyses

------------------------------------------------------------------------

## Requirements

### Software

-   `vcftools`
-   `bcftools`
-   `msmc2`
-   `msmc-tools`
-   `bgzip`
-   `python3`
-   `R`

### Input Data

-   Phased chromosome VCF files\
-   Individual BAM files\
-   Reference genome FASTA\
-   Selected individuals for MSMC analysis

------------------------------------------------------------------------

## Pipeline

------------------------------------------------------------------------

### 1. Extract Individual Phased VCF

For each chromosome (Chr01--Chr18) and each individual:

``` bash
vcftools   --gzvcf Chr${i}.phased.vcf.gz   --indv ${listname}   --recode --recode-INFO-all --stdout |   bgzip -c > ${listname}/Chr${i}.${listname}.phased.vcf.gz
```

Purpose: - Extract phased SNPs for each individual - Generate
chromosome-specific VCF files

------------------------------------------------------------------------

### 2. Calculate Sequencing Depth and Generate Mask

``` bash
bcftools mpileup -r Chr${i} -f reference.fa ${bam} | bcftools call -c -V indels | python3 bamCaller.py 2 ${listname}/Chr${i}.${listname}.mask.gz | bgzip -c > ${listname}/Chr${i}.${listname}.vcf.gz
```

Purpose: - Calculate depth information - Generate mask files - Exclude
low-quality and low-coverage regions

------------------------------------------------------------------------

### 3. Generate MSMC Input Files

``` bash
python3 generate_multihetsep.py   --mask ${listname}/Chr${i}.${listname}.mask.gz   ${listname}/Chr${i}.${listname}.phased.vcf.gz   > ${listname}/Chr${i}.${listname}.msmcinput
```

Purpose: - Convert phased VCF to MSMC input format - Apply mask
filtering

------------------------------------------------------------------------

### 4. Run MSMC2

``` bash
msmc2   -o ${listname}/${listname}   Chr01.msmcinput   Chr02.msmcinput   ...
  Chr18.msmcinput
```

Purpose: - Estimate coalescent rate (lambda) - Infer historical
effective population size

Output:

-   `${listname}.final.txt`
-   Contains time boundaries and lambda values

------------------------------------------------------------------------

## Important Parameters

-   Mutation rate (mu)
-   Generation time (gen)
-   Mask coverage threshold

------------------------------------------------------------------------

## Notes

-   High-quality phasing is critical for accurate MSMC inference.
-   Mask files help avoid bias from low-depth regions.
-   Mutation rate and generation time strongly influence scaling.
-   Results represent relative historical trends; scaling depends on
    parameter choice.

------------------------------------------------------------------------