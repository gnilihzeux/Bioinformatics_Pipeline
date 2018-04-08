# main ------------------------------------------------------------------
CODE_DIR <- sub("/$", "", CODE_DIR)					
source(paste0(CODE_DIR, "/goPathwayEnrichment.R"))					
source(paste0(CODE_DIR, "/appendXlsxResult.R"))
source(paste0(CODE_DIR, "/pictureGoPath.R"))

# differential expression gene list -------
temp_files <- list.files(go_dir, pattern = ".txt", full.names = TRUE)			
files_index <- sapply(temp_files, function(x){
  all(c("geneSymbol", "diffState") %in% 
        strsplit(readLines(x, n = 1)[[1]], split = "\t")[[1]]
  )
}
)
diff_genes <- read.table(temp_files[files_index], head = TRUE, sep = "\t", stringsAsFactors = FALSE)

if(!nrow(diff_genes))stop("No differential expression genes, please loose your threholds!")

up_genes <- diff_genes$geneSymbol[diff_genes$diffState == "up"]
down_genes <- diff_genes$geneSymbol[diff_genes$diffState == "down"]

# go/pathway enrichment ---------------------
up_enrich <- goPathwayEnrichment(degList      = up_genes,
                                 universeList = UNIVERSE_GENES,
                                 species      = SPECIES,
                                 minGeneNum   = 5,
                                 maxGeneNum   = 500
)

down_enrich <- goPathwayEnrichment(degList      = down_genes,
                                 universeList = UNIVERSE_GENES,
                                 species      = SPECIES,
                                 minGeneNum   = 5,
                                 maxGeneNum   = 500
)

# output ------------------------------------
# table
appendXlsxResult(up_enrich, 
                 workdir = go_dir,
                 Analysis_Method = "up_gene-analysis"
)
appendXlsxResult(down_enrich, 
                 workdir = go_dir,
                 Analysis_Method = "down_gene-analysis"
)
# go/pathway plots
pictureGoPath(go_dir)

