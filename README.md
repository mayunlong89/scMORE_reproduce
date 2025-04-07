# scMORE_reproduce

scMORE (single cell MultiOmics Regulon Enrichment) is designed to identify disease-relevant regulons by leveraging GWAS summary statistics and multimodal single-cell measurements, where gene expression and chromatin accessibility are either measured for each individual cell or integrated into metacell or clusters to capture both modalities. Please find the scMORE R package in [Github] (https://github.com/mayunlong89/scMORE).

Here, we provide all the necessary scripts for reproducing scMORE analyses, including those used for analyzing and visualizing ground-truth datasets, real single-cell multiomics data, and various GWAS summary statistics for complex diseases.

### Citations
Ma et al., Integrating polygenic signals and single-cell multiomics identifies cell type-specific regulomes critical for immune- and aging-related diseases and traits, `Nature Aging` (Under review), 2025


## scHOB Database
Human organoids are advanced three-dimensional structures that accurately recapitulate key characteristics of human organ development and functions. Unlike two-dimensional cultures lacking critical cell-cell communications, organoids provides a powerful model for recovering complex cellular dynamics involved in developmental and homeostatic processes. Organoids also allow genetic and pharmacological manipulation in a more physiologically relevant context compared to animal models. Although single-cell sequencing advancements have accelerated their biological and therapeutic use, there has been no systematic platform for unified processing and analysis of organoid-based single-cell multiomics data. 

We thus established scHOB (single-cell Human Organoid Bank), a multi-omic single-cell database, consisting of both scRNA-seq and scATAC-seq data on 10 types of widely-adopted human organoids (i.e., brain, lung, heart, eye, liver & bile duct, pancreas, intestine, kidney, and skin) spanning more than 1.5 million cells with 67 main cell types in 385 samples across 83 distinct protocols. see [Github code](https://github.com/mayunlong89/scHOB/tree/main); see [scHOB Website](https://schob.su-lab.org/).


### Application examples:
1. Ma et al., Integration of human organoids single-cell transcriptomic profiles and human genetics repurposes critical cell type-specific drug targets for severe COVID-19. [Cell Proliferation](https://onlinelibrary.wiley.com/doi/full/10.1111/cpr.13558),2024, and see related [Github codes](https://github.com/mayunlong89/scHuman_organoids_COVID19).
2. Ma et al., Sytematic dissection of pleiotropic loci and critical regulons in excitatory neurons and microglia relevant to neuropsychiatric and ocular diseases, [Translational Psychiatry](https://rdcu.be/d7qof), 2025. preprint version see: [Research Square](https://www.researchsquare.com/article/rs-4514542/v1), 2024.


### Other references:
1. [scPagwas](https://www.cell.com/cell-genomics/pdf/S2666-979X(23)00180-5.pdf)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8137370.svg)](https://doi.org/10.5281/zenodo.8137370)

2. Development of novel polygenic regression method scPagwas for integrating scRNA-seq data with GWAS on complex diseases. see [Ma et al. Cell Genomics, 2023](https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00180-5), and see related [Github codes](https://github.com/mayunlong89/scPagwas_main)




