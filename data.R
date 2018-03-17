#' 1. read in sample annotaion information with conditions
#' 2. read in sample counts of all
#' 

# sample groups
sample_info <- read.table(paste0(workdir, group_file), header= F, sep= "\t", stringsAsFactors= FALSE)
colnames(sample_info) <- c("sample", "label", "condition")

