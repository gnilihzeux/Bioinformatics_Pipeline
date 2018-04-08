# DESeq2 --------------------------------------------------
doDeseq <- function(countMatrix, sampleInfo){
  # differential expression analysis using DESeq2
  #
  # Args:
  #   countMatrix: read count matrix with colnames and rownames
  #                    row - gene, column - sample, value - count
  #                class - matrix/data.frame
  #   sampleInfo : sample information matrix contains at least tow columns
  #                1.sampleID - usually is the sample name, class - character
  #                2.sampleConditon - groups presents by 0,1,2,3..., class - integer
  #                class - data.frame
  # Returns:
  #   
  
  # header ------------------------------------
  if (!require(DESeq2)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
  } else {
    require(DESeq2)
  }
  
  options(stringsAsFactors = FALSE)
  options(digits = 7)
  
  # sampleInfo --------------------------------
  sample_info <- sampleInfo
  rownames(sample_info) <- sample_info$sampleID
  
  # DESeq -------------------------------------
  
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData   = sample_info,
                                design    = ~ sampleCondition
  )
  dds$sampleCondition <- factor(dds$sampleCondition, 
                                levels = sort(unique(sample_info$sampleCondition))
                                )
  
  # parameters depend on sample counts
  #' @section minReplicatesForReplace - outlier removal
  # 
  #' \code{DESeq} function will automatically replace counts with large Cook’s distance with
  #' the trimmed mean over all samples, scaled up by the size factor or normalization
  #' factor for that sample. This outlier replacement only occurs when there are 7 or
  #' more replicates, and can be turned off with \code{DESeq(dds, minReplicatesForReplace=Inf)}.
  #' 
  #' @issue
  dds <- DESeq(dds,
               test                    = "Wald",
               fitType                 = "parametric",
               quiet                   = TRUE,
               minReplicatesForReplace = 7
               )
  
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
  reslt <- results(dds,
                   cooksCutoff          = FALSE,
                   independentFiltering = FALSE,
                   alpha                = 0.1,
                   pAdjustMethod        = "BH",
                   tidy                 = TRUE
                   )
  rslt <- reslt[, c("row", "log2FoldChange", "pvalue", "padj")]
  colnames(rslt) <- c("geneID", "log2FoldChange", "pvalue", "padj")
  
  # normalized counts
  normalized_counts <- DESeq2::counts(dds, normalized= TRUE)
  normalized_counts <- data.frame(geneID = rownames(normalized_counts), normalized_counts)
  
  merge(rslt, normalized_counts)
}
