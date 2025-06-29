# radiation-mouse

**RNA-Seq Analysis of Radiation Effects in A20 Lymphoma Mouse Model**

This project was completed as part of the *Data Mining Modeling Biostatistics* course at Humber Polytechnic using R. It explores the transcriptional response to radiation exposure in a murine lymphoma model using RNA-seq count data.

## Summary

- Preprocessed raw count data using **DESeq2** with variance-stabilizing transformation.
- Performed **PCA** for sample clustering and quality assessment.
- Identified differentially expressed genes (DEGs) using DESeq2.
- Conducted pathway enrichment analysis using **fgsea**.
- Observed delayed radiation-induced upregulation of metabolic stress, immune modulation, and neuronal pathways at 7 days post-treatment.

## Files

- `annotated.R`: Main R script used for analysis.
- `R Final Report.pdf`: Summary of methods and key results.
- `GSE281695_raw_counts_all_samples.xlsx`: Raw count matrix from GEO dataset GSE281695.

## Tools & Libraries

- R (tidyverse, DESeq2, fgsea, pheatmap)

---

*For academic use and learning purposes only.*

