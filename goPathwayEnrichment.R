# go and kegg enrichment analysis -----------------------------
goPathwayEnrichment <- function(degList,
                                universeList = NULL,
                                species      = "hsa",
                                minGeneNum   = 5,
                                maxGeneNum   = 500
){
  # GO and KEGG pathway enrichment analysis with fisher's test based on local database
  #
  # Args:
  #   degList     : differential expression gene list
  #                 class-character
  #   universeList: background gene list
  #                 class-character
  #   species     : human, mouse or rat, correspongding to 'hsa', 'mmu' or 'rno'
  #                 class-character
  #   minGeneNum  : go/pathway item contains at least 5 genes
  #                 class-numeric
  #   maxGeneNum  : go/pathway item contains at most 500 genes
  #                 class-numeric
  #
  # Returns:
  #   4 tables with columns: go/pahID, go/pathDescription, goType, geneRatio, bgRatio, pvalue, padj, overlapGeneList, overlapGeneCount
  #                 class-list
  #           1. go:bp class-data.frame
  #           2. go:cc class-data.frame
  #           3. go:mf class-data.frame
  #           4. kegg  class-data.frame
  # @issue package version
  #   there are some differences between 'v3.4' and 'v3.6' of 'clusterProfiler': 
  #           1. Bioconductor version
  #                Bioc v3.4 --> v3.4
  #                Bioc v3.6 --> v3.6
  #           2. R version
  #                R v3.3.x --> v3.4
  #                R >= v3.4.2 --> v3.6
  #           3. enrichGO
  #                keytype --> v3.4
  #                keyType --> v3.6
  #  @strategy
  #    omit the argue name 'keytype/keyType' but input parameters in order
  
  # header --------------------------------------
  if (!require(clusterProfiler)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("clusterProfiler")
  } else {
    
    #' @issue how to suppress message
    require(clusterProfiler)
  }
  
  options(stringsAsFactors = FALSE)
  options(digits = 7)
  
  # input ---------------------------------------
  # species
  species <- switch(species,
                    human = "hsa",
                    hsa   = "hsa",
                    mouse = "mmu",
                    mmu   = "mmu",
                    rat   = "rno",
                    rno   = "rno"
  )
  species_db <- switch(species,
                       hsa = "org.Hs.eg.db",
                       mmu = "org.Mm.eg.db",
                       rno = "org.Rn.eg.db"
                       )
  if (!require(species_db, character.only = TRUE)) {
    # @issue
    source("http://bioconductor.org/biocLite.R")
    biocLite(species_db)
  } else {
    
    #' @issue how to suppress message
    require(species_db, character.only = TRUE)
  }
  
  
  # main  -----------------------------------
  enrichScore <- function(enrichReslt){
    # enrichment score = overlapGeneCount*bgGeneNum / (diffGeneNum*termGeneNum)
    #
    # Args:
    #   enrichReslt: enrichGO's or enrichKEGG's result, class-enrichReslt
    #
    # Returns:
    #   enrichment score, class-numeric
    overlapGeneCount <- as.numeric(sapply(strsplit(enrichReslt$GeneRatio, "/"), "[", 1))
    diffGeneNum      <- as.numeric(sapply(strsplit(enrichReslt$GeneRatio, "/"), "[", 2))
    bgGeneNum        <- as.numeric(sapply(strsplit(enrichReslt$BgRatio, "/"), "[", 2))
    termGeneNum      <- as.numeric(sapply(strsplit(enrichReslt$BgRatio, "/"), "[", 1))
    
    overlapGeneCount*bgGeneNum / (diffGeneNum*termGeneNum)
  }
   
  
  if (!is.null(universeList)) {
    enrich_bp <- enrichGO(gene          = degList,
                          OrgDb         = species_db,
                          "SYMBOL", # keytype or keyType
                          ont           = "BP",
                          pvalueCutoff  = 1,
                          pAdjustMethod = "BH",
                          universe      = universeList,
                          qvalueCutoff  = 1,
                          minGSSize     = minGeneNum,
                          maxGSSize     = maxGeneNum,
                          readable      = FALSE
    )@result
    format_bp <- data.frame(goID          = enrich_bp$ID,
                            goDescription = enrich_bp$Description,
                            goType        = "biological pathway",
                            geneRatio     = enrich_bp$GeneRatio,
                            bgRatio       = enrich_bp$BgRatio,
                            pvalue        = enrich_bp$pvalue,
                            padj          = enrich_bp$p.adjust,
                            qvalue        = enrich_bp$qvalue,
                            enrichScore   = enrichScore(enrich_bp),
                            overlapGeneList = enrich_bp$geneID,
                            overlapGeneCount = enrich_bp$Count
                            )
    
    enrich_cc <- enrichGO(gene          = degList,
                          OrgDb         = species_db,
                          "SYMBOL", # keytype or keyType
                          ont           = "CC",
                          pvalueCutoff  = 1,
                          pAdjustMethod = "BH",
                          universe      = universeList,
                          qvalueCutoff  = 1,
                          minGSSize     = minGeneNum,
                          maxGSSize     = maxGeneNum,
                          readable      = FALSE
    )@result
    format_cc <- data.frame(goID          = enrich_cc$ID,
                            goDescription = enrich_cc$Description,
                            goType        = "cellular component",
                            geneRatio     = enrich_cc$GeneRatio,
                            bgRatio       = enrich_cc$BgRatio,
                            pvalue        = enrich_cc$pvalue,
                            padj          = enrich_cc$p.adjust,
                            qvalue        = enrich_cc$qvalue,
                            enrichScore   = enrichScore(enrich_cc),
                            overlapGeneList = enrich_cc$geneID,
                            overlapGeneCount = enrich_cc$Count
    )
    
    enrich_mf <- enrichGO(gene          = degList,
                          OrgDb         = species_db,
                          "SYMBOL", # keytype or keyType
                          ont           = "MF",
                          pvalueCutoff  = 1,
                          pAdjustMethod = "BH",
                          universe      = universeList,
                          qvalueCutoff  = 1,
                          minGSSize     = minGeneNum,
                          maxGSSize     = maxGeneNum,
                          readable      = FALSE
    )@result
    format_mf <- data.frame(goID          = enrich_mf$ID,
                            goDescription = enrich_mf$Description,
                            goType        = "molecular function",
                            geneRatio     = enrich_mf$GeneRatio,
                            bgRatio       = enrich_mf$BgRatio,
                            pvalue        = enrich_mf$pvalue,
                            padj          = enrich_mf$p.adjust,
                            qvalue        = enrich_mf$qvalue,
                            enrichScore   = enrichScore(enrich_mf),
                            overlapGeneList = enrich_mf$geneID,
                            overlapGeneCount = enrich_mf$Count
                            )
    
    entrez_list <- bitr(degList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = species_db)$ENTREZID
    entrez_universe <- bitr(universeList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = species_db)$ENTREZID
    enrich_kegg <- enrichKEGG(gene          = entrez_list,
                              organism      = species,
                              "kegg", # keytype or keyType
                              pvalueCutoff  = 1,
                              pAdjustMethod = "BH",
                              universe      = entrez_universe,
                              minGSSize     = minGeneNum,
                              maxGSSize     = maxGeneNum,
                              qvalueCutoff  = 1,
                              use_internal_data = FALSE
                            )@result
    format_kegg <- data.frame(goID          = enrich_kegg$ID,
                              goDescription = enrich_kegg$Description,
                              geneRatio     = enrich_kegg$GeneRatio,
                              bgRatio       = enrich_kegg$BgRatio,
                              pvalue        = enrich_kegg$pvalue,
                              padj          = enrich_kegg$p.adjust,
                              qvalue        = enrich_kegg$qvalue,
                              enrichScore   = enrichScore(enrich_kegg),
                              overlapGeneList = enrich_kegg$geneID,
                              overlapGeneCount = enrich_kegg$Count
                              )
  } else {
    enrich_bp <- enrichGO(gene          = degList,
                          OrgDb         = species_db,
                          "SYMBOL", # keytype or keyType
                          ont           = "BP",
                          pvalueCutoff  = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff  = 1,
                          minGSSize     = minGeneNum,
                          maxGSSize     = maxGeneNum,
                          readable      = FALSE
    )@result
    format_bp <- data.frame(goID          = enrich_bp$ID,
                            goDescription = enrich_bp$Description,
                            goType        = "biological pathway",
                            geneRatio     = enrich_bp$GeneRatio,
                            bgRatio       = enrich_bp$BgRatio,
                            pvalue        = enrich_bp$pvalue,
                            padj          = enrich_bp$p.adjust,
                            qvalue        = enrich_bp$qvalue,
                            enrichScore   = enrichScore(enrich_bp),
                            overlapGeneList = enrich_bp$geneID,
                            overlapGeneCount = enrich_bp$Count
    )
    
    enrich_cc <- enrichGO(gene          = degList,
                          OrgDb         = species_db,
                          "SYMBOL", # keytype or keyType
                          ont           = "CC",
                          pvalueCutoff  = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff  = 1,
                          minGSSize     = minGeneNum,
                          maxGSSize     = maxGeneNum,
                          readable      = FALSE
    )@result
    format_cc <- data.frame(goID          = enrich_cc$ID,
                            goDescription = enrich_cc$Description,
                            goType        = "cellular component",
                            geneRatio     = enrich_cc$GeneRatio,
                            bgRatio       = enrich_cc$BgRatio,
                            pvalue        = enrich_cc$pvalue,
                            padj          = enrich_cc$p.adjust,
                            qvalue        = enrich_cc$qvalue,
                            enrichScore   = enrichScore(enrich_cc),
                            overlapGeneList = enrich_cc$geneID,
                            overlapGeneCount = enrich_cc$Count
    )
    
    enrich_mf <- enrichGO(gene          = degList,
                          OrgDb         = species_db,
                          "SYMBOL", # keytype or keyType
                          ont           = "MF",
                          pvalueCutoff  = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff  = 1,
                          minGSSize     = minGeneNum,
                          maxGSSize     = maxGeneNum,
                          readable      = FALSE
    )@result
    format_mf <- data.frame(goID          = enrich_mf$ID,
                            goDescription = enrich_mf$Description,
                            goType        = "molecular function",
                            geneRatio     = enrich_mf$GeneRatio,
                            bgRatio       = enrich_mf$BgRatio,
                            pvalue        = enrich_mf$pvalue,
                            padj          = enrich_mf$p.adjust,
                            qvalue        = enrich_mf$qvalue,
                            enrichScore   = enrichScore(enrich_mf),
                            overlapGeneList = enrich_mf$geneID,
                            overlapGeneCount = enrich_mf$Count
    )
    
    entrez_list <- bitr(degList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = species_db)$ENTREZID
    enrich_kegg <- enrichKEGG(gene          = entrez_list,
                              organism      = species,
                              "kegg", # keytype or keyType
                              pvalueCutoff  = 1,
                              pAdjustMethod = "BH",
                              minGSSize     = minGeneNum,
                              maxGSSize     = maxGeneNum,
                              qvalueCutoff  = 1,
                              use_internal_data = FALSE
    )@result
    format_kegg <- data.frame(pathID          = enrich_kegg$ID,
                              pathDescription = enrich_kegg$Description,
                              geneRatio     = enrich_kegg$GeneRatio,
                              bgRatio       = enrich_kegg$BgRatio,
                              pvalue        = enrich_kegg$pvalue,
                              padj          = enrich_kegg$p.adjust,
                              qvalue        = enrich_kegg$qvalue,
                              enrichScore   = enrichScore(enrich_kegg),
                              overlapGeneList = enrich_kegg$geneID,
                              overlapGeneCount = enrich_kegg$Count
    )
  }
  
  list(format_bp, format_cc, format_mf, format_kegg)
}

