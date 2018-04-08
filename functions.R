#' function collections

# DEA(differential expressionn analysis) input --------------------------------

#' get circRNA results from CIRCexplorer
getCircList <- function(sample_names, type= "annotated"){
  
  circ_list <- NULL
  for(sample_name in sample_names){
    
    # annotated or denovo circRNA
    if(type == "annotated"){
      #' @issue
      #' where are the files
      sample_dir <- paste0(workdir,"/circRNA/", sample_name, "/circ_out/annotate/circ_fusion.txt")
    }else{
      sample_dir <- paste0(workdir,"/circRNA/", sample_name, "/circ_out/denovo/novel_circ.txt")
    }
    
    sample_file <-  read.table(sample_dir, header= FALSE, sep= "\t", stringsAsFactors= FALSE)
    # name each circRNA with 'chr:start:end:strand' format
    circ_name <- paste(sample_file[, 1], sample_file[, 2], sample_file[, 3], sample_file[, 6], sep= ":")
    # get circRNA count
    circ_count <- sub("circular_RNA/", "", sample_file[, 4])
    
    # tidy circRNA table
    if(type == "annotated"){
      circ_tbl <- data.frame(circ_name, circ_count, sample_file[, c(1:3, 6, 10, 14:16, 18)])
      colnames(circ_tbl) <- c("circName", "count", "chrom","start","end", "strand", "exonCount", "circType","geneName","isoformName", "flankIntron")
    }else{
      circ_tbl <- data.frame(circ_name, circ_count, sample_file[, c(1:3, 6)])
      colnames(circ_tbl) <- c("circName", "count", "chrom","start","end", "strand")
    }
    
    circ_list <- c(circ_list, list(circ_tbl))
    
  }
  
  names(circ_list) <- sample_names
  circ_list
}

#' @section get circRNA count matrix
#' common circRNA among samples
#'     we should consider intra-condition and inter-condition
#'         intra-condition: biological repeats or technique replicates
getCountMatrix <- function(circ_list, col_data){
  
  samples <- col_data$sampleName
  count_list <- lapply(circ_list[samples], "[", c("circName", "count"))
  for(i in seq_len(length(samples))){
    colnames(count_list[[i]]) <- c("circName", samples[i])
  }
  count_matrix <- Reduce(function(x, y)merge(x, y, by= "circName", all= TRUE, stringsAsFactors= FALSE),
                         count_list[samples]
  )
  
  # Convert data.frame columns from factors to characters
  circ_name <- count_matrix$circName
  col_name <- colnames(count_matrix)
  count_matrix <- count_matrix[, - 1]
  count_matrix <- data.frame(lapply(count_matrix, function(x)as.integer(as.character(x))), stringsAsFactors=FALSE)
  count_matrix[is.na(count_matrix)] <- 0L
  rownames(count_matrix) <- circ_name
  colnames(count_matrix) <- col_name[-1]
  
  #' @issue
  #' expressed in at least one sample of each condition
  indx <- apply(count_matrix, 1, function(x){
                                   ctrl_col <- col_data$sampleCondition == min(col_data$sampleCondition)
                                   (sum(x[ctrl_col] > 0) >=1) & (sum(x[!ctrl_col] > 0) >=1) 
                                 }
  )
  count_matrix[indx, ]
}

#' get annotated circRNA's annotation
getAnno <- function(circ_list, col_data){
  
  samples <- col_data$sampleName
  circ_list <- unname(
                 lapply(circ_list[samples], function(x)subset(x, select= - count))
  )
  circ_tbl <- unique(do.call(rbind, circ_list))
  
  circ_matrix <- circ_tbl[, -1]
  rownames(circ_matrix) <- circ_tbl$circName
  circ_matrix
}

# DEA main --------------------------------------------------------------------

#' differential expression analysis using DESeq2
doDeseq <- function(count_data, col_data){
  
  # DEG analysis from count matrix using DESeq2
  if(!require(DESeq2)){
    source("https://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
  }else{
    
    #' @issue how to suppress message
    require(DESeq2)
  }
  
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = col_data,
                                design = ~ sampleCondition
                                )
  dds$sampleCondition <- factor(dds$sampleCondition, levels = sort(unique(col_data$sampleCondition)))
  
  # parameters depend on sample counts
  # @section minReplicatesForReplace - outlier removal
  # 
  # \code{DESeq} function will automatically replace counts with large Cook’s distance with
  # the trimmed mean over all samples, scaled up by the size factor or normalization
  # factor for that sample. This outlier replacement only occurs when there are 7 or
  # more replicates, and can be turned off with \code{DESeq(dds, minReplicatesForReplace=Inf)}.
  # 
  # @issue
  dds <- DESeq(dds, fitType= "parametric", minReplicatesForReplace = 7, quiet= TRUE)
  
  # tidy results
  #' @section cooksCutoff - filter results
  #' 
  #' The \code{results} function automatically flags genes which contain a Cook’s distance 
  #' above a cutoff for samples which have 3 or more replicates. The p values and 
  #' adjusted p values for these genes are set to NA. At least 3 replicates are 
  #' required for flagging, as it is difficult to judge which sample might be an 
  #' outlier with only 2 replicates. This filtering can be turned off with 
  #' \code{results(dds, cooksCutoff=FALSE)}.
  #' 
  #' @section independentFiltering - filter results
  reslt <- results(dds, cooksCutoff= FALSE, independentFiltering= FALSE)
  
  # normalized counts
  normalized_counts <- counts(dds, normalized= TRUE)
  
  list(results=as.data.frame(reslt), normalized_counts= normalized_counts)
}

degAssess <- function(reslt, fc= 2, fdr= 0.05, pvalue= 0.001){
  
  state <- rep("maintain", nrow(reslt))
  state[(reslt$log2FoldChange > log2(fc)) & (reslt$pvalue < pvalue) & (reslt$padj < fdr)] <- "up"
  state[(reslt$log2FoldChange < log2(1 / fc)) & (reslt$pvalue < pvalue) & (reslt$padj < fdr)] <- "down"
  
  tbl <- data.frame(reslt, diffState= state, stringsAsFactors= FALSE)
  tbl[order(tbl$diffState), ]
}

# output ----------------------------------------------------------------------
outPut <- function(diff_matrix, 
                   count_matrix,
                   col_data,
                   output_dir,
                   anno_matrix = NULL,
                   type        = "annotated"
                   ){
  
  ### #####
  circ_names <- rownames(diff_matrix)
  
  if(is.null(anno_matrix)){
    
    tbl <- data.frame(circName= circ_names, diff_matrix, count_matrix[circ_names, ], check.names= FALSE)
  }else{
    tbl <- data.frame(circName= circ_names, diff_matrix, count_matrix[circ_names, ], anno_matrix[circ_names, ], check.names= FALSE)
    
    #' @issue colnames
    #' @issue one host but many circRNA
    # only contain two columns: circRNA host, diffState
    gn_tbl <- tbl[tbl$diffState != "maintain", c("geneName", "diffState")]
    colnames(gn_tbl) <- c("geneSymbol", "diffState")
  }
  rownames(tbl) <- tbl$circName

  ### write #####
  if(type == "annotated"){
    
    #' @issue output file name
    write.table(tbl, 
                paste0(output_dir,
                       "/annotated-circRNA_result", 
                       #"_(FC",  FOLD_CHANGE, "-Padj", FOLD_CHANGE, "-Pvalue", PVALUE, ")", 
                       #"_condition-", condition_treat, "-vs-", condition_ctrl,
                       ".txt"),
                sep= "\t", quote= F, row.names= FALSE, col.names= TRUE)
    
    write.table(gn_tbl, 
                paste0(output_dir, "/Go_analysis",
                       "/annotated_diff-circ", 
                       #"_(FC",  FOLD_CHANGE, "-Padj", FOLD_CHANGE, "-Pvalue", PVALUE, ")", 
                       #"_condition-", condition_treat, "-vs-", condition_ctrl,
                       ".txt"),
                sep= "\t", quote= F, row.names= FALSE, col.names= TRUE)
    
    
  }else{
    
    write.table(tbl, 
                paste0(reslt_dir,
                       "/novel-circRNA_result", 
                       #"_(FC",  FOLD_CHANGE, "-Padj", FOLD_CHANGE, "-Pvalue", PVALUE, ")", 
                       #"_condition-", condition_treat, "-vs-", condition_ctrl,
                       ".txt"),
                sep= "\t", quote= F, row.names= FALSE, col.names= TRUE)
    
  }

  # return -----------------------
  tbl
  
}

