#!/usr/bin/env Rscript

# parameters --------------------------------------------------------------------
#' @issue how to transfer parameters like linux command
args <- commandArgs(trailingOnly= TRUE)
workdir <- args[1]
group_file <- args[2]
SPECIES <- args[3]
FOLD_CHANGE <- as.numeric(args[4])
FDRVALUE <- as.numeric(args[5])
PVALUE <- as.numeric(args[6])
#OUTPUT_DIR <- args[6]

CODE_DIR <- "/data3/circRNA/bin"

# configure ---------------------------------------------------------------------
options(stringsAsFactors= FALSE)
source(paste(CODE_DIR, "data.R", sep= "/"))
source(paste(CODE_DIR, "functions.R", sep= "/"))
library(pheatmap)

# perform differential expression analysis --------------------------------------

#' any two conditions among samples
#' e.g there are 3 conditions 1, 2 and 3, then
#'     we can get three condition combinations, 
#'     that is (2, 1), (3, 1) and (3, 2).
uniq_condition <- unique(sample_info$condition)
condition_combn <- combn(uniq_condition, 2)

for(combn_index in seq_len(ncol(condition_combn))){
  
  # extract one from all combination
  condition_ctrl <- min(condition_combn[, combn_index])
  condition_treat <- max(condition_combn[, combn_index])
  # subset of sample information
  col_data <- sample_info[sample_info$condition %in% c(condition_ctrl, condition_treat), ]
  rownames(col_data) <- col_data$sample
  
  # output directory
  result_dir <- paste0(workdir,"circRNA/Diff_Result")
  if(!dir.exists(result_dir))dir.create(result_dir)
  
  reslt_dir <- paste0(result_dir, "/condition_", condition_treat, "vs", condition_ctrl)
  if(!dir.exists(reslt_dir))dir.create(reslt_dir)
  
  go_dir <- paste0(reslt_dir, "/Go_analysis")
  if(!dir.exists(go_dir))dir.create(go_dir)
  
  heatmap_dir <- paste0(reslt_dir, "/heatmap")
  if(!dir.exists(heatmap_dir))dir.create(heatmap_dir)
  
  
  # DEA
  # 1. annotated
  circ_list <- getCircList(col_data$sample, type= "annotated")
  count_data <- getCountMatrix(circ_list, col_data)
  reslt <- doDeseq(count_data, col_data)
  anno_info <- getAnno(circ_list, col_data)
  
  diff_circ <- degAssess(reslt$results, fc= FOLD_CHANGE, fdr= FDRVALUE, pvalue= PVALUE)
  
  out_tbl <- outPut(diff_circ,
                    reslt$normalized_counts,
                    col_data,
                    output_dir  = reslt_dir, 
                    anno_matrix = anno_info,
                    type= "annotated")
  
  
  #' @issue this should be a new function
  ### go pathway -----------------
  GO_dir <- go_dir
  source(paste(CODE_DIR, "Go-pathway_mainscript.R", sep= "/"))
  
  # heatmap ----------------------
  heatmap_input <- out_tbl[, col_data$sample]
  heatmap_labels <- pheatmap(heatmap_input,
                             scale = "row", 
                             clustering_distance_rows = "euclidean",
                             clustering_distance_cols = "euclidean",
                             color = colorRampPalette(c("blue", "white", "red"))(50),
                             cluster_rows = TRUE,
                             cluster_cols = TRUE, 
                             border_color = FALSE,
                             #main = "ssssssssss",
                             annotation_col = col_data["condition"],
                             show_rownames = TRUE, 
                             show_colnames = TRUE, 
                             number_color = "blue",
                             height=8,width=9,filename=paste0(heatmap_dir,"/annotated_heatmap.pdf"),fontsize = 10)
  
  
  
  # 2. novel
  
  circ_list <- getCircList(col_data$sample, type= "novel")
  count_data <- getCountMatrix(circ_list, col_data)
  reslt <- doDeseq(count_data, col_data)
  
  diff_circ <- degAssess(reslt$results, fc= FOLD_CHANGE, fdr= FDRVALUE, pvalue= PVALUE)
  out_tbl <- outPut(diff_circ, 
                    reslt$normalized_counts, 
                    col_data,
                    output_dir  = result_dir, 
                    anno_matrix = NULL,
                    type= "novel")  
  
  # heatmap ----------------------
 
  heatmap_input <- out_tbl[, col_data$sample]
  heatmap_labels <- pheatmap(heatmap_input,
                             scale = "row", 
                             clustering_distance_rows = "euclidean",
                             clustering_distance_cols = "euclidean",
                             color = colorRampPalette(c("blue", "white", "red"))(50),
                             cluster_rows = TRUE,
                             cluster_cols = TRUE, 
                             border_color = FALSE,
                             #main = "ssssssssss",
                             annotation_col = col_data["condition"],
                             show_rownames = TRUE, 
                             show_colnames = TRUE, 
                             number_color = "blue",
                             height=8,width=9,filename=paste0(heatmap_dir,"/annotated_heatmap.pdf"),fontsize = 10)
}

