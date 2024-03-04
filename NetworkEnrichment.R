#!/usr/bin/env Rscript

## Script name: NetworkEnrichment.R
##
## Purpose of script: Network Enrichment and Drug Repurposing
##
## Author: Klaudia Adamowicz
## Email: klaudia.adamowicz@uni-hamburg.de
##
## Date Created: 2024-03-04
##
## Copyright (c) Dr. Tanja Laske, 2024

# Load required libraries --------------------------
suppressPackageStartupMessages({
  required_packages <- c("optparse","argparse","data.table", "tibble", "igraph", "tidyverse", "rjson", "reticulate")
  for(package in required_packages){
    if(!require(package,character.only = TRUE)) install.packages(package)
    library(package, character.only = TRUE)
  }
  if(!require("limma",character.only = TRUE, quietly = TRUE)) BiocManager::install("limma")
  library("limma", character.only = TRUE)
  
  reticulate::use_virtualenv("/home/python_env", required = TRUE)
  ds <- reticulate::import("drugstone")
  ds$print_license()
  ds$accept_license()
})


# Define Methods --------------------------

#' Method performing limma
#'
#' @param set_condition Vector of experimental design specifying the condition(s)
#' @param set_counts Counts (rows = proteins, columns = samples)
#' @param set_comparisons Vector of comparisons
#'
#' @return Results of the limma analysis
#'
#' @export
perform_de <- function(set_condition, set_counts, set_comparisons){
  # create design matrix
  groupsM <- as.factor(set_condition)
  designM <- model.matrix(~0+groupsM) 
  colnames(designM) <-levels(groupsM) 
  fit <- lmFit(set_counts, designM)
  
  # create contrasts
  contr <- makeContrasts(contrasts = set_comparisons, levels = colnames(coef(fit)))
  
  fit2 <- contrasts.fit(fit, contr)
  ebfit <- eBayes(fit2, trend = TRUE)
  return(ebfit)
}

#' Method to extract results of fit object
#'
#' @param set_fit Fit object of the perform_de method
#' @param set_comparisons Vector of comparisons
#' @param out_dir Output directory
#' @param lfc_up Log fold change threshold (upregulated)
#' @param lfc_down Log fold change threshold (downregulated)
#' @param alpha p-value or p.adjust threshold (depending on padj)
#' @param strict TRUE if < and > should be used for threshold comparison, FALSE if <= and >=
#' @param padj TRUE if alpha should be applied to padj, FALSE if alpha applied to p-value
#' @param logFC_thr TRUE if a threshold for logFC should be set, FALSE if not
#'
#' @return Extracted results from the fit object
#'
#' @export
compare_de_expr <- function(set_fit, set_comparisons, out_dir, lfc_up = args$logFC_up, lfc_down = args$logFC_down, alpha = args$alpha, 
                            strict = FALSE, padj = args$p_adj, logFC_thr=args$logFC, write_output = TRUE){
  de_results <- list()
  for (i in 1:length(set_comparisons)){
    # get results for each comparison
    top.table <- topTable(set_fit, sort.by = "P", number=Inf, coef=c(i))
    gene_reg <- setDT(top.table, keep.rownames = "gene") # save row names as column
    
    # different threshold settings
    if (logFC_thr){
      if (strict){
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$logFC > lfc_up & gene_reg$adj.P.Val < alpha, 1, ifelse(gene_reg$logFC < lfc_down & gene_reg$adj.P.Val < alpha, 1, 0))
        } else {
          gene_reg$Change <- ifelse(gene_reg$logFC > lfc_up & gene_reg$P.Value < alpha, 1, ifelse(gene_reg$logFC < lfc_down & gene_reg$P.Value < alpha, 1, 0))
        }
      } else {
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$logFC >= lfc_up & gene_reg$adj.P.Val <= alpha, 1, ifelse(gene_reg$logFC <= lfc_down & gene_reg$adj.P.Val <= alpha, 1, 0))
        } else {
          gene_reg$Change <- ifelse(gene_reg$logFC >= lfc_up & gene_reg$P.Value <= alpha, 1, ifelse(gene_reg$logFC <= lfc_down & gene_reg$P.Value <= alpha, 1, 0))
        }
      }
    } else {
      if (strict){
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$adj.P.Val < alpha, 1, 0)
        } else {
          gene_reg$Change <- ifelse(gene_reg$P.Value < alpha, 1, 0)
        }
      } else {
        if (padj){
          gene_reg$Change <- ifelse(gene_reg$adj.P.Val <= alpha, 1, 0)
        } else {
          gene_reg$Change <- ifelse(gene_reg$P.Value <= alpha, 1, 0)
        }
      }
    }
    if (write_output){
      write.table(gene_reg, file = paste0(out_dir,"de_data/","DEdata.", set_comparisons[[i]],".txt"),sep=" ",row.names = FALSE) 
    }
    de_results[[set_comparisons[[i]]]] <- gene_reg
  }
  return(de_results)
}

#' Prepare matrix for differential expression analysis
#'
#' @param counts Counts (rows = proteins, columns = samples)
#' @param md Metadata
#' @param condition_name Name of the condition in metadata
#' @param comparisons Vector of comparisons
#' @param out_dir Output directory
#' @param network_proteins Proteins in the network
#' @param lfc_up Log fold change threshold for upregulated proteins
#' @param lfc_down Log fold change threshold for downregulated proteins
#' @param alpha p-value or p.adjust threshold (depending on padj)
#' @param strict TRUE if < and > should be used for threshold comparison, FALSE if <= and >=
#' @param padj TRUE if alpha should be applied to padj, FALSE if alpha applied to p-value
#' @param logFC_thr TRUE if a threshold for logFC should be set, FALSE if not
#' @param plot TRUE to generate plots, FALSE otherwise
#' @param write_output TRUE to write output files, FALSE otherwise
#'
#' @return Prepared matrix for differential expression analysis
#'
#' @export
prepare_matrix <- function(counts, md, condition_name, comparisons, out_dir,
                           lfc_up = args$logFC_up, lfc_down = args$logFC_down, alpha = args$alpha, 
                           strict = FALSE, padj = args$p_adj, logFC_thr=args$logFC, plot = TRUE, write_output = TRUE){
  if (write_output){
    # create output directories
    dir.create(out_dir, showWarnings = FALSE) #stops warnings if folder already exists
    dir.create(file.path(out_dir,"de_data"), showWarnings = FALSE) #stops warnings if folder already exists
    
  }
  
  # get condition vector
  condition <- md[, condition_name]
  
  # perform DE analysis pairwise
  pairwise_de <- perform_de(set_condition = condition, set_counts = counts, set_comparisons = comparisons)
  # extract results
  de_results <- compare_de_expr(set_fit = pairwise_de, set_comparisons = comparisons, out_dir = out_dir, lfc_up = lfc_up, 
                                lfc_down = lfc_down, alpha = alpha, strict = strict, padj = padj, logFC_thr = logFC_thr, write_output = write_output)
  mat = data.frame(matrix(ncol=1,nrow=0, dimnames=list(NULL, c("gene"))))
  for (comp in comparisons) {
    cur_df <- de_results[[comp]][,c("gene","Change"), drop = FALSE]
    cur_df <- cur_df[!is.na(cur_df$Change), ]
    colnames(cur_df) <- c("gene", comp)
    mat = merge(mat, cur_df, all=TRUE)
  }
  
  # Fill NA with 0 
  mat[is.na(mat)] <- 0
  
  # filter for only expressed rows 
  mat <- mat[rowSums(mat==0) != (ncol(mat)-1), ]
  # explode id column
  mat <- mat %>% mutate(gene=strsplit(gene, ";")) %>% unnest(gene)
  return(mat)
}


protein_to_gene_mapping <- function(count_data, ind_mat, gene_col) {
  # Split and explode the Protein.IDs column
  long_format <- count_data[, .(Protein.ID = unlist(strsplit(Protein.IDs, ";"))), by = .(Gene = get(gene_col))]
  
  # Remove duplicate rows
  long_format <- unique(long_format)
  
  # Rename column in ind_mat
  setnames(ind_mat, old = "gene", new = "Protein.ID")
  
  # Merge with ind_mat
  merged_data <- merge(ind_mat, long_format, by = "Protein.ID", all.x = TRUE)
  
  return(merged_data)
}


save_meta_data <- function(comp, cond, meta_data, out_dir, filename_prefix) {
  meta_all <- list()
 
  sample_one <- strsplit(comp, split="-")[[1]][1]
  sample_two <- strsplit(comp, split="-")[[1]][2]
  
  # save meta data into json 
  sample_group <- if (cond == "TimeCond") list("Timepoint", "Condition") else if (cond == "TimeCondLoc") list("Timepoint", "Condition", "Location") else list("Condition")
  meta_all[[comp]] = list(
    samples_group_A = meta_data[meta_data[[cond]] %in% c(sample_one)]$Column_name,
    samples_group_B = meta_data[meta_data[[cond]] %in% c(sample_two)]$Column_name,
    group_A = str_split(sample_one,"_")[[1]],
    group_B = str_split(sample_two,"_")[[1]],
    sample_group = sample_group
  )
 
  # Add links to the entire meta_all list
  meta_all_links <- list(
    druglist = file.path(out_dir, paste0(filename_prefix, "_drugs.csv")),
    genelist = file.path(out_dir, paste0(filename_prefix, "_genelist.csv")),
    graph_df = file.path(out_dir, paste0(filename_prefix, "_graph.csv")),
    graph_file = file.path(out_dir, paste0(filename_prefix, "_network.graphml")))
  
  # Combine metadata with links
  meta_all_combined <- list(meta_data = meta_all, links = meta_all_links)
  # Write to JSON
  write(toJSON(meta_all_combined), file.path(out_dir, paste0(filename_prefix,".txt.json")))
}

## ------------- Prepare Network Enrichment ------------------

### Parse arguments -----

# set up arguments
parser <- OptionParser()
parser <- add_option(parser, c("-m","--meta_file"), help="Meta data description file")
parser <- add_option(parser, c("-c","--count_file"), help="Preprocessed count file")
parser <- add_option(parser, c("-o","--out_dir"), help="Directory for output files", default="")

# Adding option for column selection

parser <- add_option(parser, c("-g","--gene_column"), help="Gene column with human gene names.")

# Adding new options for thresholds with defaults
parser <- add_option(parser, c("--logFC"), help = "Boolean specifying whether to apply a logFC threshold (TRUE) or not (FALSE)", type = "logical", default = TRUE)
parser <- add_option(parser, c("--logFC_up"), help = "Upper log2 fold change threshold (dividing into upregulated)", type = "numeric", default = 1)
parser <- add_option(parser, c("--logFC_down"), help = "Lower log2 fold change threshold (dividing into downregulated)", type = "numeric", default = -1)
parser <- add_option(parser, c("--p_adj"), help = "Boolean specifying whether to apply a threshold on adjusted p-values (TRUE) or on raw p-values (FALSE)", type = "logical", default = TRUE)
parser <- add_option(parser, c("--alpha"), help = "Threshold for adjusted p-values or p-values", type = "numeric", default = 0.05)

# get command line options, if help option encountered print help and exit
args <- parse_args(parser)
# check if mandatory options are set
check_options <- function(tags){
  for (tag in tags){
    if (is.null(args[tag])){
      print_help(parser)
      stop("Missing mandatory option.", call.=FALSE)
    }
  }
}
check_options(c('meta_file','count_file','gene_column'))

# save arguments
meta_file_path <- "proteomics-genevention/example_data/plasma/meta_data.csv" #args$meta_file
count_file_path <- "proteomics-genevention/example_data/plasma/IRS_on_Median_normalized_data.csv"#args$count_file
out_dir <- "proteomics-genevention/example_data/test/netenrich"#args$out_dir
gene_column <- "Orthologs"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE) #stops warnings if folder already exists

#### Load data --------
meta_data <- fread(meta_file_path)
count_data <- fread(count_file_path)

## Prepare data ----

#### Correct data ----
ID_column <- "Animal"

# remove ref
meta_data <- meta_data[meta_data[[ID_column]] != "ref",]

# convert timepoint column
meta_data[, Timepoint := sapply(Timepoint, function(tp) ifelse(grepl("pre", tp), -as.numeric(gsub("pre", "", tp)), ifelse(grepl("post", tp), as.numeric(gsub("post", "", tp)), as.numeric(tp))))]
meta_data[, Sample_name := if ("Sample_name" %in% names(meta_data)) Sample_name else if ("Label" %in% names(meta_data)) Label else NULL]
meta_data[, Column_name := if ("Column_name" %in% names(meta_data)) Column_name else if ("Column" %in% names(meta_data)) Column else NULL]
#meta_data[, Timepoint := as.numeric(Timepoint)]

### Rename Columns ---

# rename from file name to sample name
names(count_data) <- plyr::mapvalues(names(count_data), from = meta_data$Column_name, to = meta_data$Sample_name, warn_missing=FALSE) 
meta_data <- subset(meta_data, meta_data$Sample_name %in% names(count_data))

# remove columns that are not in meta_data
columns_to_keep <- c("Protein.IDs", gene_column, meta_data$Sample_name)
existing_columns <- columns_to_keep[columns_to_keep %in% names(count_data)]
count_data <- count_data[, existing_columns, with = FALSE]

#### Encode Columns ----

mappings <- list(
  Condition = list(setNames(paste0("C", seq_along(unique(meta_data$Condition))), unique(meta_data$Condition))),
  Location = if("Location" %in% names(meta_data) && length(unique(meta_data$Location)) > 1) 
    list(setNames(paste0("L", seq_along(unique(meta_data$Location))), unique(meta_data$Location))),
  Timepoint = list(setNames(paste0("T", seq_along(sort(unique(meta_data$Timepoint)))), sort(unique(meta_data$Timepoint))))
)

write(toJSON(mappings), file.path(out_dir, "Metadata_encoding.txt.json")) 

reverse_mappings <- list(
  Condition = if (!is.null(mappings$Condition)) setNames(names(mappings$Condition[[1]]), unlist(mappings$Condition[[1]])) else NULL,
  Location = if (!is.null(mappings$Location)) setNames(names(mappings$Location[[1]]), unlist(mappings$Location[[1]])) else NULL,
  Timepoint = if (!is.null(mappings$Timepoint)) setNames(names(mappings$Timepoint[[1]]), unlist(mappings$Timepoint[[1]])) else NULL
)

write(toJSON(reverse_mappings), file.path(out_dir, "Metadata_reverse_encoding.txt.json")) 

# rename condition
meta_data[, Condition := mappings$Condition[[1]][Condition]]
# rename timepoint
meta_data[, Timepoint := as.character(Timepoint)]
meta_data[, Timepoint := mappings$Timepoint[[1]][Timepoint]]
# rename location
if(!is.null(mappings$Location)) {
  meta_data[, Location := mappings$Location[[1]][Location]]
}

### Create indicator matrix ----
#### Prepare input data ----

# create count matrix 
counts <- as.data.frame(count_data[,c(meta_data$Sample_name), with=F], )
rownames(counts) <- count_data$Protein.IDs
counts <- counts[rowSums(is.na(counts)) != ncol(counts), ] # remove where full row is NA


#### Assign active and inactive cases ----
##### Conditions ----
###### Save comparisons ----
comparisons <- c()
elements <- sort(unique(meta_data$Condition))
for (index_a in 1:(length(elements)-1)){
  for (index_b in (index_a+1):length(elements)){
    comparisons <- c(comparisons, paste0(elements[index_a], "-", elements[index_b]))
  }
}
ind_mat_cond <- prepare_matrix(counts=counts, md=as.data.frame(meta_data), condition_name="Condition", 
                               comparisons=comparisons, out_dir = out_dir, write_output = FALSE)
ind_mat_cond <- protein_to_gene_mapping(count_data = count_data, ind_mat = ind_mat_cond, gene_col = gene_column)
write.table(ind_mat_cond, file.path(out_dir, "indicator_matrix_cond.tsv"), sep="\t", row.names=FALSE)


##### Conditions and Timepoints ----
# save sample groups
target_column <- ifelse("Location" %in% names(meta_data) && length(unique(meta_data$Location)) > 1, 
                        "TimeCondLoc", 
                        "TimeCond")
# Create the target column based on the presence of 'Location'
meta_data[, (target_column) := if(target_column == "TimeCondLoc") {
  paste(meta_data$Timepoint, meta_data$Condition, meta_data$Location, sep = "_")} else {
    paste(meta_data$Timepoint, meta_data$Condition, sep = "_")}]

###### Save comparisons ----
comparisons <- c()
elements <- sort(unique(meta_data[[target_column]]))
for (index_a in 1:(length(elements)-1)){
  for (index_b in (index_a+1):length(elements)){
    # only if time point or condition are the same!
    split_a = str_split(elements[index_a],"_")[[1]]
    split_b = str_split(elements[index_b],"_")[[1]]
    # Count the number of differences
    differences <- sum(split_a != split_b)
    # Only add to comparisons if exactly one of timepoint, condition, or location is different
    if (differences == 1){
      comparisons <- c(comparisons, paste0(elements[index_a], "-", elements[index_b]))
    }
  }
}

ind_mat_timecond <- prepare_matrix(counts=counts, md=as.data.frame(meta_data), condition_name=target_column, 
                                   comparisons=comparisons, out_dir = out_dir, write_output = FALSE)
ind_mat_timecond <- protein_to_gene_mapping(count_data = count_data, ind_mat = ind_mat_timecond, gene_col = gene_column)
write.table(ind_mat_timecond, file.path(out_dir,"indicator_matrix_timecond.tsv"), sep="\t", row.names=FALSE)


# Run Network Enrichment -----------------------

## Functions ----

run_drugstone_ami <- function(genes){
  parameters = list("algorithm" = "multisteiner", "target" = "drug-target", 
                    "tolerance" = 5, "hubPenalty" = 0.5, "num_trees" = 5, 
                    "ppiDataset"= 'nedrex', "identifier" = "symbol")
  ds_out <- ds$new_task(genes, parameters)
  ds_results <- ds_out$get_result()
  ds_genes <- ds_results$get_genes()
  ds_edges <- do.call(rbind, lapply(ds_results$get_edges(), function(x) { data.frame(from = x$from, to = x$to) }))
  ds_graph <- graph_from_edgelist(as.matrix(ds_edges), directed = FALSE)
  return(list("ds_genes" = names(ds_genes), "ds_graph" = ds_graph, "ds_edges" = ds_edges))
}

run_drugstone_drugs <- function(genes, edges, num_drugs = 20){
  ds_out_2 <- ds$new_task(genes, list("algorithm" = "trustrank", "target" = "drug","pdiDataset"= 'nedrex', "includeNonApprovedDrugs" = TRUE, "resultSize" = num_drugs))
  ds_results_2 <- ds_out_2$get_result()
  ds_drugs_2 <- ds_results_2$get_drugs()
  ds_genes <- ds_results_2$get_genes()
  ds_edges_2 <- do.call(rbind, lapply(ds_results_2$get_edges(),  function(x) { data.frame(from = x$from, to = x$to) }))
  ds_edges_2 <- rbind(ds_edges_2, edges)
  ds_graph_2 <- graph_from_edgelist(as.matrix(ds_edges_2), directed = FALSE)
  return(list("ds_drugs" = ds_drugs_2, "ds_gene" = names(ds_genes), "ds_graph" = ds_graph_2))
}

save_results <- function(ds_result, filename){
  if ( "ds_drugs" %in% names(ds_result) ) {
    ds_drugs <- do.call(rbind, lapply(ds_result$ds_drugs, as.data.frame))
    ds_drugs <- ds_drugs %>% group_by(across(-hasEdgesTo)) %>%
      summarize(hasEdgesTo = paste(hasEdgesTo, collapse = ";"))
    write.table(ds_drugs, file.path(out_dir, paste0(filename, "_drugs.csv")), quote=FALSE, sep=",", row.names = FALSE)
  }
  if ( "ds_gene" %in% names(ds_result) ) {
    ds_gene <- ds_result$ds_gene
    write.table(ds_gene, file.path(out_dir, paste0(filename, "_genelist.csv")), quote=FALSE, sep=",", row.names = FALSE)
  }
  if ( "ds_graph" %in% names(ds_result) ) {
    ds_graph <- ds_result$ds_graph
    # save as graphml
    write.graph(ds_graph, file.path(out_dir, paste0(filename, "_graph.graphml")), format = "graphml")
    write.table(get.data.frame(ds_graph), file.path(out_dir, paste0(filename, "_graph.csv")), quote=FALSE, sep=",", row.names = FALSE)
  }
}


## condition vs condition ----

net_results <- list()

if (nrow(ind_mat_cond) > 0) {
  unique_genes <- unique(unlist(strsplit(na.omit(gsub("^$|^NA$", "", ind_mat_cond$Gene)), ";")))
  ds_results <- run_drugstone_ami(genes = unique_genes)
  ds_results_2 <- run_drugstone_drugs(ds_results$ds_genes, ds_results$ds_edges,
                                      num_drugs = 100)
  
  net_results[["cond"]][[names(ind_mat_cond)[2]]] <- list("gene_net" = ds_results,
                                                          "drug_net" = ds_results_2)
  save_results(ds_result = ds_results_2, filename = "cond")
  save_meta_data(comp = names(ind_mat_cond)[2], cond = "Condition", meta_data = meta_data,
                 out_dir = out_dir, filename_prefix = "cond")
  
} else {
  warning("No data found. Skipping condition vs condition analysis.")
}
## time point condition vs time point condition ----

net_results[["tp"]] <- list()

for ( case in setdiff(names(ind_mat_timecond), c("Protein.ID", "Gene"))){
  if (sum(ind_mat_timecond[[case]]) > 0) {
    df <- ind_mat_timecond[c("Protein.ID", "Gene", case)]  %>% filter(.data[[case]] != 0)
    unique_genes <- unique(unlist(strsplit(na.omit(gsub("^$|^NA$", "", df$Gene)), ";")))
    if ( length(unique_genes) > 0) {
      ds_results <- run_drugstone_ami(genes = unique_genes)
      ds_results_2 <- run_drugstone_drugs(ds_results$ds_genes, ds_results$ds_edges,
                                          num_drugs = 100)
      net_results[["tp"]][[case]] <- list("gene_net" = ds_results,
                                          "drug_net" = ds_results_2)
      
      ### save data
      if (!is.null(net_results[["tp"]][[case]])){
        save_results(ds_result = ds_results_2, filename = paste0("timepoint_",case))
        save_meta_data(comp = case, cond = "Condition", meta_data = meta_data,
                       out_dir = out_dir, filename_prefix = paste0("timepoint_",case))
      }
    }
  } else {
    warning(paste0("No data found. Skipping ", case, " analysis."))
  }
}

# Overview statistics ----

# Initialize maps for counting unique values and drug IDs
top_comparisons <- list()
drugID_info <- list()

# Loop through each main part and sub-part of kpm_results
for (main_part in names(net_results)) {
  for (sub_part in names(net_results[[main_part]])) {
    # Convert sub-part to dataframe
    drug_df <- do.call(rbind, lapply(net_results[[main_part]][[sub_part]]$drug_net$ds_drugs, 
                                     function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
    
    # Count unique values in hasEdgesTo column
    unique_edges <- unique(unlist(drug_df$hasEdgesTo))
    # Update top_comparisons with the number of unique edges for the sub_part
    top_comparisons[[sub_part]] <- length(unique_edges)
    
    # Count and store drugID and label
    for (row in 1:nrow(drug_df)) {
      drug_id <- drug_df$drugId[[row]]
      label <- drug_df$label[[row]]
      status <- ifelse(drug_df$status[[row]] == "approved", "approved", "unapproved")
      edges_to <- unique(drug_df$hasEdgesTo[[row]])
      
      if (is.null(drugID_info[[drug_id]])) {
        drugID_info[[drug_id]] <- list(count = 1, name = label, status = status, edges_to = edges_to)
        } else {
        drugID_info[[drug_id]]$count <- drugID_info[[drug_id]]$count + 1
        drugID_info[[drug_id]]$edges_to <- unique(c(drugID_info[[drug_id]]$edges_to, edges_to))
      }
    }
  }
}

## Top comparisons ----

# Convert top_comparisons list into a dataframe
top_comparisons_df <- data.frame(comparison = names(top_comparisons),
                                 drugtargets = unlist(top_comparisons))

# Sort the dataframe by drugtargets in descending order
top_comparisons_df <- top_comparisons_df[order(-top_comparisons_df$drugtargets),]

write.table(top_comparisons_df, file.path(out_dir,"top-comparisons.tsv"), sep="\t", row.names=FALSE)

## Top IDs ----

# Convert drugID_info into a dataframe
drug_info_df <- do.call(rbind, lapply(names(drugID_info), function(drug_id) {
  info <- drugID_info[[drug_id]]
  num_drugtargets <- length(info$edges_to)             # Count of drug targets
  drugtargets <- paste(info$edges_to, collapse = ";")  # Concatenate the edges_to vector into a single string
  data.frame(
    drugId = drug_id, 
    label = info$name, 
    status = info$status, 
    count = info$count, 
    drugtargets = drugtargets, 
    num_drugtargets = num_drugtargets, 
    stringsAsFactors = FALSE
  )
}))

# Convert the result into a dataframe if it's not already
drug_info_df <- as.data.frame(drug_info_df)


# Sort the dataframe (example: sorting by count in descending order)
drug_info_df <- drug_info_df[order(-drug_info_df$num_drugtargets),]

write.table(drug_info_df, file.path(out_dir,"top-drugs.tsv"), sep="\t", row.names=FALSE)

# Extracting all drug targets and counting occurrences
all_drugtargets <- unlist(strsplit(drug_info_df$drugtargets, ";"))
target_counts <- table(all_drugtargets)

# Check presence in 'Gene' column of ind_mat_cond or ind_mat_timecond
is_seed <- (names(target_counts) %in% ind_mat_cond$Gene) | (names(target_counts) %in% ind_mat_timecond$Gene)

# Create the final dataframe
drugtarget_info_df <- data.frame(ID = names(target_counts), 
                                 seed = is_seed, 
                                 drugs = as.integer(target_counts), 
                                 stringsAsFactors = FALSE)
drugtarget_info_df <- drugtarget_info_df[order(-drugtarget_info_df$drugs),]

write.table(drugtarget_info_df, file.path(out_dir,"top-ids.tsv"), sep="\t", row.names=FALSE)
