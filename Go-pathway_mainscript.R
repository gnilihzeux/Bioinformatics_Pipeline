#!/usr/bin/Rscript
##path

path1 <- GO_dir

GO_DB_PATH <- "GO-PATH-DB"

#setwd(path1)
#the parameters for screening differentially expressed genes

#the parameters for Go-Pathway 
#species : "hsa"(default),"mmu","rno"

species = SPECIES

## source_method : "Up_Down Gene","Target Up_Down Gene","Target Differential Gene","Differential Gene"
source_method = "Up_Down Gene"
######################################################################
#ctrl + A
#
options(java.parameters = "-Xmx15000m")
source(paste(CODE_DIR, 'Go-pathway-analysis.R', sep= "/"), encoding = 'UTF-8')
library(XLConnect)
library(xlsx)

Analysis_Method <<- "Path-analysis"
Analysis_Species <<- species
Res_Source <<- source_method
#Ana_Wid =  paste(path1,"/gopath",sep="")
Ana_Wid =  path1
outputdir <<- paste(path1,"/output",sep="")
if (!file.exists(outputdir)) { dir.create(outputdir) }

Res_Source = source_method

Pathway_Analysis_Run()

 
#Ana_Wid =  paste(path1,"/gopath",sep="")

Analysis_Method <<- "GO-analysis"


Pathway_Analysis_Run()

library(ggplot2)
Go_picture(path=path1,method = source_method)
Path_picture(path=path1,method = source_method)


