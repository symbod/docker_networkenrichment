#!/usr/bin/env nextflow

params.output = "./output/"
params.meta_file = '' // Path to meta.txt
params.count_file = '' // Path to count.txt
params.gene_column = '' // Name of gene column with human genes

params.logFC = true
params.logFC_up = 1
params.logFC_down = -1
params.p_adj = true
params.alpha = 0.05

// data
meta_file = file(params.meta_file)
count_file = file(params.count_file)

// scripts
netenrich_analysis_script = file("${projectDir}/NetworkEnrichment.R")

process network_enrichment {
    container 'kadam0/netenrichment:0.0.2'
    publishDir params.output, mode: "copy"

    input:
    path script_file
    path meta_file
    path count_file

    output:                                
    path "*"

    script:
    """
    Rscript $script_file --meta_file ${meta_file} --count_file ${count_file} --gene_column ${params.gene_column} --out_dir ${params.output} --logFC ${params.logFC} --logFC_up ${params.logFC_up} --logFC_down ${params.logFC_down} --p_adj ${params.p_adj} --alpha ${params.alpha}
    """
}

workflow {
  network_enrichment(netenrich_analysis_script, meta_file, count_file)
}
