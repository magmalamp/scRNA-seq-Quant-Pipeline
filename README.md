# scRNA-seq Quantitative Pipeline

This is a Snakemake port of an HPC-based alevin-fry single-cell RNA-seq lab workflow, reimplemented to run locally and reproducibly on 10x Chromium data. The pipeline runs end-to-end, from mapping with salmon, to cell-barcode/UMI resolution and quantification with alevin-fry. As implemented, the workflow is demonstrated on the 10x Genomics "1k PBMCs from a Healthy Donor" dataset.

![Pipeline Directed Acyclic Graph](plots/quant_dag_full.svg)
