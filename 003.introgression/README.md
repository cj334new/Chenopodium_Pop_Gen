# ABBA-BABA (D-statistic) Sliding Window Analysis Pipeline

This document describes a standardized, publication-ready pipeline for
performing **ABBA-BABA (D-statistic)** analysis in sliding windows
across the genome using phased genotype data.

The pipeline is suitable for detecting **introgression and gene flow**
among populations and complements population genetic analyses such as
nucleotide diversity (π), Fst, XP-CLR, and IBD analyses.

------------------------------------------------------------------------

## Overview

-   **ABBA-BABA (D-statistic)**\
    A four-taxon test used to detect excess allele sharing between
    non-sister populations, indicating gene flow.

-   **Sliding window framework**\
    Calculates D-statistics across genomic windows to identify localized
    introgression signals.

-   **Typical applications**

    -   Detection of introgressed genomic regions\
    -   Identifying asymmetric gene flow\
    -   Complementing selection sweep analyses\
    -   Studying domestication or adaptive introgression

------------------------------------------------------------------------

## Requirements

### Input Files

-   `ALL_samples.vcf.gz`\
    Genotype file (gzipped), formatted for ABBABABAwindows.py.

-   `sample_pop.txt`\
    Population assignment file mapping individuals to populations.

------------------------------------------------------------------------

## Pipeline

### 1. ABBA-BABA Sliding Window Analysis

Genome-wide D-statistics are calculated in sliding windows using phased
genotype data.


`https://github.com/simonhmartin/genomics_general/blob/master/ABBABABAwindows.py`

``` bash
python ABBABABAwindows.py     -g ALL_samples.vcf.gz     -f phased     -o Q3toQ2.csv     -w 1000000     -m 100     -s 100000     -P1 Q1     -P2 Q2     -P3 Q3     -O Outgroup     -T 10     --minData 0.5     --popsFile sample_pop.txt     --writeFailedWindows
```

------------------------------------------------------------------------

## Key Parameters

-   `-g` : Input genotype file (gzipped).
-   `-f phased` : Specifies that input data are phased.
-   `-o` : Output CSV file containing window-based D-statistic results.
-   `-w 1000000` : Window size (1 Mb).
-   `-s 100000` : Step size (100 kb).
-   `-m 100` : Minimum number of informative sites per window.
-   `-P1` : Population 1 (reference population).
-   `-P2` : Population 2.
-   `-P3` : Test population (candidate introgressed population).
-   `-O` : Outgroup population.
-   `-T 10` : Number of threads.
-   `--minData 0.5` : Minimum proportion of non-missing data required
    per window.
-   `--popsFile` : Population assignment file.
-   `--writeFailedWindows` : Outputs windows that fail filtering
    criteria.

------------------------------------------------------------------------

## Output Files

-   `Q3toQ2.csv`\
    Contains window-based D-statistic values and associated statistics.

-   Optional failed windows file\
    Lists windows excluded due to insufficient data or filtering
    thresholds.

------------------------------------------------------------------------

## Notes

-   Input genotype data should be high quality and properly phased when
    using `-f phased`.
-   Window size and step size can be adjusted depending on marker
    density and genome size.
-   Interpretation of D-statistic:
    -   **D \> 0**: Excess allele sharing between P3 and P2\
    -   **D \< 0**: Excess allele sharing between P3 and P1\
-   Significance assessment may require block-jackknife or Z-score
    estimation depending on downstream implementation.

------------------------------------------------------------------------

