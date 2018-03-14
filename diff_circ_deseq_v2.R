#!/usr/bin/env Rscript

# parameters --------------------------------------------------------------------
args <- commandArgs(trailingOnly= TRUE)
workdir <- args[1]
group_file <- args[2]
FOLD_CHANGE <- args[3]
FDRVALUE <- args[4]
PVALUE <- args[5]
OUTPUT_DIR <- args[6]

# configure ---------------------------------------------------------------------
options(stringsAsFactors= FALSE)


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
  
  # DEA
  # 1. annotated
  circ_list <- getCircList(col_data$sample, type= "annotated")
  count_data <- getCountMatrix(circ_list, col_data)
  reslt <- doDeseq(count_data, col_data)
  anno_info <- getAnno(circ_list, col_data)
  
  diff_circ <- degAssess(reslt$results)
  
  outPut(diff_circ, reslt$normalized_counts, anno_info, type= "annotated")
  
  # 2. novel
  
  circ_list <- getCircList(col_data$sample, type= "novel")
  count_data <- getCountMatrix(circ_list, col_data)
  reslt <- doDeseq(count_data, col_data)
  
  diff_circ <- degAssess(reslt$results)
  outPut(diff_circ, reslt$normalized_counts, type= "novel")
  
}

