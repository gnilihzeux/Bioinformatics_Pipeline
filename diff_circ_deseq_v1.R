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
UNIVERSE_GENES <- NULL
#OUTPUT_DIR <- args[6]

SPECIES <- switch(SPECIES,
                  human = "hsa",
                  hsa   = "hsa",
                  mouse = "mmu",
                  mmu   = "mmu",
                  rat   = "rno", 
                  rno   = "rno"
                  )

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
uniq_condition <- unique(sample_info$sampleCondition)
condition_combn <- combn(uniq_condition, 2)

for(combn_index in seq_len(ncol(condition_combn))){
  
  # extract one from all combination
  condition_ctrl <- min(condition_combn[, combn_index])
  condition_treat <- max(condition_combn[, combn_index])
  # subset of sample information
  col_data <- sample_info[sample_info$sampleCondition %in% c(condition_ctrl, condition_treat), ]
  rownames(col_data) <- col_data$sampleName
  
  # output directory
  result_dir <- paste0(workdir,"/circRNA/2-Diff_Result")
  if(!dir.exists(result_dir))dir.create(result_dir)
  
  reslt_dir <- paste0(result_dir, "/condition_", condition_treat, "vs", condition_ctrl)
  if(!dir.exists(reslt_dir))dir.create(reslt_dir)
  
  go_dir <- paste0(reslt_dir, "/Go_analysis")
  if(!dir.exists(go_dir))dir.create(go_dir)
  
  heatmap_dir <- paste0(reslt_dir, "/heatmap")
  if(!dir.exists(heatmap_dir))dir.create(heatmap_dir)
  
  
  # DEA
  # 1. annotated
  circ_list <- getCircList(col_data$sampleName, type= "annotated")
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
  # go pathway -----------------
  source(paste(CODE_DIR, "outputGoPathway.R", sep= "/"))

  # heatmap ----------------------
  e_tbl <- out_tbl[out_tbl$diffState != "maintain", col_data$sampleName]
  rownames(e_tbl) <- out_tbl$circName
  
  source(paste(CODE_DIR, "HeatMap.R", sep= "/"))
  HeatMap(e_tbl, annotation_col = col_data["sampleCondition"], path = heatmap_dir, Figurename = "annotated_heatmap.pdf")
  
  # 2. novel
  
  circ_list <- getCircList(col_data$sampleName, type= "novel")
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
  e_tbl <- out_tbl[out_tbl$diffState != "maintain", col_data$sampleName]
  rownames(e_tbl) <- out_tbl$circName
  
  source(paste(CODE_DIR, "HeatMap.R", sep= "/"))
  HeatMap(e_tbl, annotation_col = col_data["sampleCondition"], path = heatmap_dir, Figurename = "novel_heatmap.pdf")
  
}

