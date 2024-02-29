
# README for Network Enrichment Analysis Pipeline

## Overview
This repository contains a Nextflow pipeline for Network Enrichment Analysis, primarily designed to process proteomics and genomics data. The pipeline focuses on integrating differential expression data with network data to identify significant pathways and network enrichments under various conditions.

## Usage
To run the pipeline, you need to specify several parameters:

1. `--meta_file`: Path to the metadata file, which includes sample information like experimental conditions.
2. `--count_file`: Path to the count file, which contains the expression data for analysis.
3. `--network_file`: Path to the network file, providing the interaction network information.
4. `--logFC`: Boolean flag to apply a log fold change (logFC) threshold (default: `true`).
5. `--logFC_up`: Upper log2 fold change threshold for upregulated elements (default: `1`).
6. `--logFC_down`: Lower log2 fold change threshold for downregulated elements (default: `-1`).
7. `--p_adj`: Boolean flag to use adjusted p-values (default: `true`).
8. `--alpha`: Significance threshold for p-values (default: `0.05`).
9. `--output`: Directory for storing output files (default: `"./output/"`).

Example command:
```
./nextflow run main.nf --meta_file path/to/meta.txt --count_file path/to/count.txt --network_file path/to/network.txt --output output_directory
```

## Output Description
The Network Enrichment Analysis pipeline generates several types of outputs, including:

- **Enrichment Analysis Results**: Comprehensive results of the network enrichment analysis, highlighting key pathways and networks.
- **Statistical Summaries**: Statistical analysis summaries, including p-values, fold changes, and other relevant metrics.
- **Visualization Files**: Graphs and charts that visually represent the analysis results, facilitating easier interpretation and discussion.

## Getting Started
To get started, clone this repository and ensure that Nextflow is installed on your system. Prepare your metadata, count, and network files in the required format. Then execute the pipeline using the command with the correct paths to your input files. Results will be stored in the specified output directory.

---
**Note**: It's important to adjust paths, filenames, and parameters according to your specific project requirements and data structure.