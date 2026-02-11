# Haplotype Phasing and Identity-by-Descent (IBD) Analysis Pipeline

This document describes a standardized, publication-ready pipeline for **haplotype phasing** and **identity-by-descent (IBD)** detection based on a genome-wide VCF file.  
The pipeline is designed to complement population genetic analyses such as nucleotide diversity (π) and genetic differentiation (Fst), and is suitable for downstream selection and demographic inference.

---

## Overview

- **Haplotype phasing**  
  Infers phased haplotypes from unphased genotype data using a statistical model.

- **Identity-by-Descent (IBD)**  
  Detects shared haplotype segments inherited from a recent common ancestor.

- **Typical applications**
  - Demographic history inference
  - Detection of recent shared ancestry
  - Complementary evidence for selection and population structure analyses

---

## Requirements

### Input Files

- `ALL_samples.vcf.gz`  
  BGZF-compressed and indexed VCF file containing all samples.

---

## Pipeline

### 1. Haplotype Phasing with Beagle

Genotype phasing is performed using **Beagle**, which statistically infers haplotypes for each individual based on linkage disequilibrium patterns.

```shell
java -Xmx50g -jar beagle.25Nov19.28d.jar \
    gt=ALL_samples.vcf.gz \
    out=ALL_samples.phased \
    nthreads=20
```

**Key parameters**
- `-Xmx50g`  
  Maximum Java heap memory allocation.
- `gt`  
  Input VCF file.
- `out`  
  Output prefix for phased VCF.
- `nthreads`  
  Number of computational threads.

**Output files**
- `ALL_samples.phased.vcf.gz`  
  Phased genotype data.
- `ALL_samples.phased.log`  
  Beagle log file.

---

### 2. Identity-by-Descent (IBD) Detection with Refined IBD

IBD segments are identified from the phased haplotypes using **Refined IBD**, which improves detection accuracy by refining candidate shared segments.

```shell
java -jar refined-ibd.17Jan20.102.jar \
    gt=ALL_samples.phased.vcf.gz \
    length=0.002 \
    out=ALL_samples.phased_0.002Mb
```
---

## Notes

- Accurate IBD detection requires high-quality phased data; therefore, phasing should be performed prior to any IBD analysis.
---
