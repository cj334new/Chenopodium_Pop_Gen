
# Nucleotide Diversity (π) and Genetic Differentiation (Fst) Analysis Pipeline

This document describes a standardized, publication-ready pipeline for calculating **nucleotide diversity (π)** within populations and **genetic differentiation (Fst)** between populations using sliding-window approaches.

---

## Overview

- **π (Nucleotide Diversity)**: Measures genetic variation within a population.  
- **Fst (Fixation Index)**: Quantifies genetic differentiation between populations.  
- **Sliding-window strategy**:
  - Window size: **50 kb**
  - Step size: **5 kb**

This pipeline is suitable for population genomics studies and downstream selection-scan analyses.

---

## Requirements

### Input Files

- `ALL_samples.vcf.gz`  
  BGZF-compressed and indexed VCF file containing all samples.
- `${pop}.list`  
  Plain-text files listing sample IDs (one per line) for each population.

---

## Pipeline

### 1. Define Population Identifiers

```shell
pop1="pop1"   # Population 1 identifier
pop2="pop2"   # Population 2 identifier
```

---

### 2. Calculate Nucleotide Diversity (π)

Nucleotide diversity is calculated independently for each population using a sliding-window approach.

```shell
for pop in $pop1 $pop2; do
    vcftools --gzvcf ALL_samples.vcf.gz \           # Input VCF file (gzipped)
        --keep ${pop}.list \                         # Sample list for current population
        --window-pi 50000 \                          # Window size: 50,000 base pairs
        --window-pi-step 5000 \                      # Step size: 5,000 base pairs (sliding window)
        --out ${pop}.win50Kstep5K.pi                 # Output file prefix
done
```

**Output files**
- `${pop}.win50Kstep5K.pi.windowed.pi`

**Key columns**
- `CHROM`: Chromosome
- `BIN_START`: Window start position
- `BIN_END`: Window end position
- `PI`: Nucleotide diversity estimate

---

### 3. Calculate Genetic Differentiation (Fst)

Genetic differentiation between populations is estimated using Weir and Cockerham’s Fst.

```shell
vcftools --gzvcf ALL_samples.vcf.gz \                # Input VCF file (gzipped)
    --weir-fst-pop ${pop1}.list \                    # Sample list for population 1
    --weir-fst-pop ${pop2}.list \                    # Sample list for population 2
    --fst-window-size 50000 \                        # Window size: 50,000 base pairs
    --fst-window-step 5000 \                         # Step size: 5,000 base pairs (sliding window)
    --out ${pop1}_${pop2}.win50Kstep5k.fst           # Output file prefix 
```

**Output file**
- `${pop1}_${pop2}.win50Kstep5k.fst.windowed.weir.fst`

**Key columns**
- `CHROM`
- `BIN_START`
- `BIN_END`
- `MEAN_FST`

---

### 4. Integrate π and Fst Results

π values from both populations and Fst values are merged by genomic window for downstream analyses.

```shell
python pi_fst_integrator.py \
    ${pop1} \                                          # Population 1 identifier
    ${pop1}.win50Kstep5K.pi.windowed.pi \              # π file for population 1 (50kb window, 5kb step)
    ${pop2} \                                          # Population 2 identifier 
    ${pop2}.win50Kstep5K.pi.windowed.pi \              # π file for population 2 (50kb window, 5kb step)
    ${pop1}_${pop2}.win50Kstep5k.fst.windowed.weir.fst \  # Fst file between populations (50kb window, 5kb step)
    ${pop1}_${pop2}.win50Kstep5k.pi.fst.tab            # Output merged file (tab-delimited)
```

**Final output**
- `${pop1}_${pop2}.win50Kstep5k.pi.fst.tab`

The integrated table can be directly used for visualization (e.g., Manhattan plots) and selection-scan analyses.

---

## Notes

- Consistent window and step sizes are required for correct data integration.
- This pipeline can be extended with additional statistics (e.g., XP-CLR, Tajima’s D).

---
