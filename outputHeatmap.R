# main ------------------------------------------------------------------
CODE_dir <- sub("/$", "", CODE_dir)					
reslt_dir <- sub("/$", "", reslt_dir)

# 
expression_tbl <- read.table(paste0(reslt_dir), "/annotated-circRNA_result.txt",
                             header = TRUE,
                             sep = "\t",
                             check.names = FALSE,
                             stringsAsFactors = FALSE
)
e_tbl <- expression_tbl[expression_tbl$diffState != "maintain", col_data$sampleName]
rownames(e_tbl) <- e_tbl$circName

HeatMap(e_tbl, annotation_col = col_data$sampleCondition, path = heatmap_dir)
