#' 1. read in sample annotaion information with conditions
#' 2. read in sample counts of all
#' 

# sample groups
group_file <- sub("/$", "", group_file)
sample_info <- read.table(paste(workdir, group_file, sep = "/"), header= TRUE, sep= "\t", stringsAsFactors= FALSE)
rownames(sample_info) <- sample_info$sampleName
