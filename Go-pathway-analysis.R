

# Setup color gradation function
# green to black to red (default)
# Also can be #FFFFEE, #0037BA ...
filled_cols <- function(lowi = "green", highi = "red", ncolors = 10) {
  low <- col2rgb(lowi)/255
  high <- col2rgb("black")/255
  col1 <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
  low <- col2rgb("black")/255
  high <- col2rgb(highi)/255
  col2 <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
  cols_distribution <- c(col1[1:(ncolors-1)],col2)
  return(cols_distribution)
}
# Comprehensive normalization method for data visualization.
comprehensive_normalization <- function(x,control_method="median",standardize=FALSE,column_norm=FALSE,
                                        calibration=FALSE,calibrationFile=NULL,...) {
  if (mode(x) != "numeric") x <- data.matrix(x)
  if (calibration && file.exists(calibrationFile)) {
    calibrationData <- read.table(calibrationFile,header=TRUE,sep="\t")
    print(paste0("calibration to: ", calibrationFile))
    x <- overlapSets(x,calibrationData)
  } else {
    if (control_method=="median" && standardize && column_norm) {
      x <- t(scale(t(x),center=apply(x,1,median)))
      x <- scale(x,center=apply(x,2,median))
    } else if (control_method=="median" && standardize) {
      x <- t(scale(t(x),center=apply(x,1,median)))
    } else if (control_method=="median") {
      x <- x - apply(x,1,median,na.rm=T)
    } else if (control_method=="mean" && standardize & column_norm) {
      x <- t(scale(t(x),center=rowMeans(x)))
      x <- scale(x,center=colMeans(x))
    } else if (control_method=="mean" && standardize) {
      x <- t(scale(t(x),center=rowMeans(x)))
    } else if (control_method=="mean") {
      x <- x - rowMeans(x)
    } else {
      x
    }
  }
  return(x)
}

Pathway_Analysis_Run <- function(...) {
if (Res_Source == "Up_Down Gene") {
  if (Analysis_Method == "Path-analysis") { Table_1_Name <<- "上调基因显著性Pathway"; Table_2_Name <<- "显著性Pathway所包含上调基因"; Table_3_Name <<- "下调基因显著性Pathway"; Table_4_Name <<- "显著性Pathway所包含下调基因" } else { Table_1_Name <<- "上调基因显著性功能"; Table_2_Name <<- "显著性功能所包含上调基因"; Table_3_Name <<- "下调基因显著性功能"; Table_4_Name <<- "显著性功能所包含下调基因"}
} else if (Res_Source == "Differential Gene") {
  if (Analysis_Method == "Path-analysis") { Table_1_Name <<- "差异基因显著性Pathway"; Table_2_Name <<- "显著性Pathway所包含差异基因" } else { Table_1_Name <<- "差异基因显著性功能"; Table_2_Name <<- "显著性功能所包含差异基因" }
} else if (Res_Source == "Target Up_Down Gene") {
  if (Analysis_Method == "Path-analysis") { Table_1_Name <<- "上调MicroRNA预测靶基因显著性Pathway"; Table_2_Name <<- "显著性Pathway所包含上调MicroRNA预测的靶基因"; Table_3_Name <<- "下调MicroRNA预测靶基因显著性Pathway"; Table_4_Name <<- "显著性Pathway所包含下调MicroRNA预测的靶基因" } else {Table_1_Name <<- "上调MicroRNA预测靶基因显著性功能"; Table_2_Name <<- "显著性功能所包含上调MicroRNA预测的靶基因"; Table_3_Name <<- "下调MicroRNA预测靶基因显著性功能"; Table_4_Name <<- "显著性功能所包含下调MicroRNA预测的靶基因" }
} else {
  if (Analysis_Method == "Path-analysis") { Table_1_Name <<- "差异MicroRNA预测靶基因显著性Pathway"; Table_2_Name <<- "显著性Pathway所包含差异MicroRNA预测的靶基因" } else { Table_1_Name <<- "差异MicroRNA预测靶基因显著性功能"; Table_2_Name <<- "显著性功能所包含差异MicroRNA预测的靶基因" }
}




  
  Significant_Analysis_Method <- function(X) {                
    p_value <- 2*fisher.test(matrix(c(diff_gene_in_sets[X],diff_gene_amount-diff_gene_in_sets[X],gene_in_sets[X]-diff_gene_in_sets[X],residual_gene[X]),nr=2))$p.value
    return(p_value)       
  }
  
  # library(RMySQL)
  #drv <- dbDriver("MySQL")
  # con <- dbConnect(drv,dbname="go_path_analysis",user="devuser",password="111111",host="192.168.2.10") 
  gene_list <<- read.csv(paste(CODE_DIR, "/", GO_DB_PATH, "/pathway/",Analysis_Species,"_gene_list.csv",sep = ""))
  #gene_list <- dbGetQuery(con,paste("select ",Analysis_Species,"_gene.gene_id,",Analysis_Species,"_gene.gene_symbol ","from ",Analysis_Species,"_gene",sep=""))
  gene_amount <- nrow(gene_list)
  
  if (Analysis_Method == "Path-analysis") {
    column_description_1_id <<- c("path_id","path_name","category","path_diffgene_count","path_gene_count","enrichment","pvalue","FDR")
    column_description_2_id <<- c("path_id","path_name","category","enrichment","pvalue","FDR","gene_id","gene_name")      
    pathway_list <<- read.csv(paste(CODE_DIR, "/", GO_DB_PATH, "/pathway/",Analysis_Species,"_pathway_list.csv",sep = ""))
    relation_list<<- read.csv(paste(CODE_DIR, "/", GO_DB_PATH, "/pathway/",Analysis_Species,"_relation_list.csv",sep = ""))
    
    
  } else {
    
    
    
    column_description_1_id <<- c("go_id","go_name","type","go_diffgene_count","go_gene_count","enrichment","pvalue","FDR")
    column_description_2_id <<- c("go_id","go_name","type","enrichment","pvalue","FDR","gene_id","gene_name")
    relation_list <<- read.csv(paste(CODE_DIR, "/", GO_DB_PATH, "/GO/",Analysis_Species,"_relation_list.csv",sep = ""),stringsAsFactor=F)
    pathway_list <- read.csv(paste(CODE_DIR, "/", GO_DB_PATH, "/GO/",Analysis_Species,"_go_list.csv",sep = ""),stringsAsFactor=F)
    
  }
  
  
  
  
  Analysis_File_Lists <- list.files(Ana_Wid,pattern=".csv$|.txt$|.xls$|.xlsx",full.names=TRUE)
  for (single_file in Analysis_File_Lists) {
    diff_gene_list <- suffix_file_read(single_file)
    up_diff_gene_list <<- diff_gene_list[which(diff_gene_list$style == "up"),]
    down_diff_gene_list <<- diff_gene_list[which(diff_gene_list$style == "down"),]
    New_File_Directory <<- paste(outputdir,"/",substr(single_file,max(unlist(gregexpr('[/]',single_file)))+1,nchar(single_file)-4),sep="")
    # if (!file.exists(New_File_Directory)) {  dir.create(New_File_Directory)  }
    Create_File_Description <- Create_File_Doc(File_Details_Directory = paste(New_File_Directory,"_",Analysis_Method,".xlsx",sep=""))
    if (nrow(up_diff_gene_list) != 0) {
      comput_diff_gene <- match(up_diff_gene_list[,1],gene_list$gene_symbol)[!is.na(match(up_diff_gene_list[,1],gene_list$gene_symbol))]
      if (identical(comput_diff_gene[1],integer(0))) { svalue(status_bar) <- sprintf("Zero Matched Genes !"); stop() };
      comput_diff_gene_id <- gene_list$gene_id[comput_diff_gene]
      diff_gene_amount <- length(comput_diff_gene_id)  
      required_relation_list <- relation_list[relation_list$gene_id %in% comput_diff_gene_id,]
      if (Analysis_Method == "Path-analysis") {
        new_relation_list <- relation_list[relation_list$path_id %in% required_relation_list$path_id,]
        diff_gene_in_sets <- as.vector(table(required_relation_list$path_id))
        gene_in_sets <- as.vector(table(new_relation_list$path_id))  
        residual_gene <- gene_amount-diff_gene_amount-gene_in_sets+diff_gene_in_sets
        #print("Diff Gene Total:", diff_gene_in_sets,"All Gene Total:", gene_in_sets, residual_gene);
        p_value <- unlist(lapply(X=1:length(residual_gene),FUN=Significant_Analysis_Method)) # p_value <- unlist(multicore:::mclapply(X=1:length(residual_gene),FUN=Significant_Analysis_Method))
        required_pathway_list <- pathway_list[match(sort(unique(required_relation_list$path_id)),pathway_list$path_id),]
      } else if (Analysis_Method == "GO-analysis") {
        new_relation_list <- relation_list[relation_list$go_id %in% required_relation_list$go_id & relation_list$go_id %in% pathway_list$go_id,]
        required_relation_list <- required_relation_list[required_relation_list$go_id %in% new_relation_list$go_id,]
        diff_gene_in_sets <- as.vector(table(required_relation_list$go_id))
        gene_in_sets <- as.vector(table(new_relation_list$go_id))  
        residual_gene <- gene_amount-diff_gene_amount-gene_in_sets+diff_gene_in_sets
        p_value <- unlist(lapply(X=1:length(residual_gene),FUN=Significant_Analysis_Method))
        # p_value <- unlist(multicore:::mclapply(X=1:length(residual_gene),FUN=Significant_Analysis_Method))
        required_pathway_list <- pathway_list[match(sort(unique(new_relation_list$go_id)),pathway_list$go_id),]
      }
      Enrichment_Score <- diff_gene_in_sets*gene_amount/(diff_gene_amount*gene_in_sets)
      Ori_Ana_Result <- cbind(required_pathway_list,diff_gene_in_sets,gene_in_sets,Enrichment_Score,p_value)
      Ana_Result <- Ori_Ana_Result[order(Ori_Ana_Result$p_value),]
      Ana_Result$p_value <- ifelse(2*Ana_Result$p_value <= 1,2*Ana_Result$p_value,1)
      FDR <- p.adjust(Ana_Result$p_value,method="fdr")
      Ana_Result <- cbind(Ana_Result,FDR)
      colnames(Ana_Result) <- column_description_1_id           
      # write.table(Ana_Result,paste(New_File_Directory,"_",Analysis_Method,".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
      # write.xlsx2(Ana_Result,paste(New_File_Directory,"_",Analysis_Method,".xlsx",sep=""),sheetName=Table_1_Name,row.names=FALSE,append=TRUE)
      Append_XLSX_Result(File_Name=Table_1_Name,Res_Data=Ana_Result,Fill_Row_Number=sum(Ana_Result$pvalue < 0.05))
      result_gene_list <- merge(Ana_Result,required_relation_list,sort=FALSE)
      gene_symbol <- gene_list$gene_symbol[match(result_gene_list$gene_id,gene_list$gene_id)]
      pathway_gene_list <- cbind(result_gene_list,gene_symbol)[,-c(4,5)]
      colnames(pathway_gene_list) <- column_description_2_id
      # write.xlsx2(pathway_gene_list,paste(New_File_Directory,"_",Analysis_Method,".xlsx",sep=""),sheetName=Table_2_Name,row.names=FALSE,append=TRUE)
      Append_XLSX_Result(File_Name=Table_2_Name,Res_Data=pathway_gene_list,Fill_Row_Number=sum(pathway_gene_list$pvalue < 0.05))
      
    }
    
    if (nrow(down_diff_gene_list) != 0) {
      comput_diff_gene <- match(down_diff_gene_list[,1],gene_list$gene_symbol)[!is.na(match(down_diff_gene_list[,1],gene_list$gene_symbol))]
      if (identical(comput_diff_gene[1],integer(0))) { print("Zero Matched Genes !"); stop() }
      comput_diff_gene_id <- gene_list$gene_id[comput_diff_gene]
      diff_gene_amount <- length(comput_diff_gene_id)
      required_relation_list <- relation_list[relation_list$gene_id %in% comput_diff_gene_id,]
      if (Analysis_Method == "Path-analysis") {
        new_relation_list <- relation_list[relation_list$path_id %in% required_relation_list$path_id,]
        diff_gene_in_sets <- as.vector(table(required_relation_list$path_id))
        gene_in_sets <- as.vector(table(new_relation_list$path_id))
        residual_gene <- gene_amount-diff_gene_amount-gene_in_sets+diff_gene_in_sets
        p_value <- unlist(lapply(X=1:length(residual_gene),FUN=Significant_Analysis_Method)) 
        # p_value <- unlist(multicore:::mclapply(X=1:length(residual_gene),FUN=Significant_Analysis_Method))
        required_pathway_list <- pathway_list[match(sort(unique(required_relation_list$path_id)),pathway_list$path_id),]
      } else if (Analysis_Method == "GO-analysis") {
        new_relation_list <- relation_list[relation_list$go_id %in% required_relation_list$go_id & relation_list$go_id %in% pathway_list$go_id,]
        required_relation_list <- required_relation_list[required_relation_list$go_id %in% new_relation_list$go_id,]
        diff_gene_in_sets <- as.vector(table(required_relation_list$go_id))
        gene_in_sets <- as.vector(table(new_relation_list$go_id))
        residual_gene <- gene_amount-diff_gene_amount-gene_in_sets+diff_gene_in_sets
        p_value <- unlist(lapply(X=1:length(residual_gene),FUN=Significant_Analysis_Method))
        # p_value <- unlist(multicore:::mclapply(X=1:length(residual_gene),FUN=Significant_Analysis_Method))
        required_pathway_list <- pathway_list[match(sort(unique(new_relation_list$go_id)),pathway_list$go_id),]
      }
      Enrichment_Score <- diff_gene_in_sets*gene_amount/(diff_gene_amount*gene_in_sets)
      Ori_Ana_Result <- cbind(required_pathway_list,diff_gene_in_sets,gene_in_sets,Enrichment_Score,p_value)
      Ana_Result <- Ori_Ana_Result[order(Ori_Ana_Result$p_value),]
      Ana_Result$p_value <- ifelse(2*Ana_Result$p_value <= 1,2*Ana_Result$p_value,1)
      FDR <- p.adjust(Ana_Result$p_value,method="fdr")
      Ana_Result <- cbind(Ana_Result,FDR)
      colnames(Ana_Result) <- column_description_1_id
      # write.xlsx2(Ana_Result,paste(New_File_Directory,"_",Analysis_Method,".xlsx",sep=""),sheetName=Table_3_Name,row.names=FALSE,append=TRUE)
      Append_XLSX_Result(File_Name=Table_3_Name,Res_Data=Ana_Result,Fill_Row_Number=sum(Ana_Result$pvalue < 0.05))
      result_gene_list <- merge(Ana_Result,required_relation_list,sort=FALSE)
      gene_symbol <- gene_list$gene_symbol[match(result_gene_list$gene_id,gene_list$gene_id)]
      pathway_gene_list <- cbind(result_gene_list,gene_symbol)[,-c(4,5)]
      colnames(pathway_gene_list) <- column_description_2_id
      # write.xlsx2(pathway_gene_list,paste(New_File_Directory,"_",Analysis_Method,".xlsx",sep=""),sheetName=Table_4_Name,row.names=FALSE,append=TRUE)
      Append_XLSX_Result(File_Name=Table_4_Name,Res_Data=pathway_gene_list,Fill_Row_Number=sum(pathway_gene_list$pvalue < 0.05))
      
    }
  }
}





suffix_file_read <- function(File_Directory) {
  parts <- strsplit(File_Directory,".",fixed=TRUE)
  # substr(svalue(tbl[1,2]),nchar(svalue(tbl[1,2]))-2,nchar(svalue(tbl[1,2])))
  if (parts[[1]][length(parts[[1]])] == "txt") {
    Table.Content <- read.table(File_Directory,sep="\t",header=TRUE,stringsAsFactors=FALSE)
  } else if (parts[[1]][length(parts[[1]])] == "csv") {
    Table.Content <- read.delim(File_Directory,sep=",",header=TRUE,stringsAsFactors=FALSE)
  } else {
    library(gdata)
    Table.Content <- read.xls(File_Directory,sheet=1,header=TRUE,stringsAsFactors=FALSE)
  }
}

Create_File_Doc <- function(File_Details_Directory) {
  
  # Create a new workbook
  wb <- createWorkbook();
  # wb = loadWorkbook(paste(New_File_Directory,"_",Analysis_Method,".xlsx",sep=""));
  # Create a new sheet with a name
  sheet1 = createSheet(wb, "文件说明");    
  cells = createCell(createRow(sheet1, 1:33), 1:3);
  Header_head_doc <- c("表名称","列头","列头的含义");
  Column_Header_Style <- CellStyle(wb,alignment = Alignment(h="ALIGN_CENTER"),font = Font(wb,isBold = TRUE),
                                   border = Border(position=c("BOTTOM","LEFT", "TOP", "RIGHT"),pen="BORDER_THIN"),fill=Fill(backgroundColor="gray"))
  tmp_0 = mapply(function(x) setCellValue(cells[[1,x]], Header_head_doc[x]), 1:3);    
  tmp_0 <- mapply(function(x) setCellStyle(cells[[1,x]],Column_Header_Style), 1:3)
  
  # Merge the first column into one cell
  addMergedRegion(sheet1, 2, 9, 1, 1);
  addMergedRegion(sheet1, 10, 17, 1, 1);
  # Create the style for title cell
  title_cell_style = CellStyle(wb,
                               alignment = Alignment(h="ALIGN_CENTER",v="VERTICAL_CENTER"),
                               font = Font(wb, color="blue",name="Times New Roman",isBold = TRUE));
  title_cell_1 = getCells(getRows(sheet1, 2), 1)[[1]];
  every_cell_frame <- CellStyle(wb,font=Font(wb,name="Times New Roman"),border = Border(position=c("BOTTOM","LEFT", "TOP", "RIGHT"),pen="BORDER_THIN"))
    if (Analysis_Method == "Path-analysis") {
    column_description_1 <- c("pathway索引号，与KEGG数据库的pathway索引号一致","pathway名称，与KEGG数据库中pathway的命名方式一致","pathway 所属的类别","差异基因计数，表示所有差异基因中属于某一pathway的差异基因数量","基因计数，表示数据库中属于某一pathway的基因数量","富集度，若p值相同，富集度越大的pathway，表示该pathway受到实验的影响越大","p值，评估pathway的显著性水平，(p<0.05表示pathway具有显著性差异,用黄色标注)","误判率，对p值准确率的判断，对pathway显著性水平的再判断")
    column_description_2 <- c(column_description_1[-c(4,5)],"基因索引号，与NCBI GenBank的基因索引号一致","基因名称，与NCBI GenBank的基因标识一致")

    if (nrow(up_diff_gene_list) != 0 & nrow(down_diff_gene_list) != 0) {
      tmp0 <- mapply(function(x) setCellStyle(cells[[x,1]],every_cell_frame), 2:33)
      setCellValue(title_cell_1, Table_1_Name);
      setCellStyle(title_cell_1, title_cell_style);          
      title_cell_2 = getCells(getRows(sheet1, 10), 1)[[1]];
      setCellValue(title_cell_2, Table_2_Name);
      setCellStyle(title_cell_2, title_cell_style);
      addMergedRegion(sheet1, 18, 25, 1, 1);title_cell_3 = getCells(getRows(sheet1, 18), 1)[[1]];
      setCellValue(title_cell_3, Table_3_Name);
      setCellStyle(title_cell_3, title_cell_style);
      addMergedRegion(sheet1, 26, 33, 1, 1);title_cell_4 = getCells(getRows(sheet1, 26), 1)[[1]];
      setCellValue(title_cell_4, Table_4_Name);
      setCellStyle(title_cell_4, title_cell_style);
      setColumnWidth(sheet1,1,30)
      tmp_1 <- mapply(function(x,y) setCellValue(cells[[x,2]], column_description_1_id[y]), 2:9, 1:8)
      tmp_2 <- mapply(function(x,y) setCellValue(cells[[x,2]], column_description_2_id[y]), 10:17, 1:8)
      tmp_3 <- mapply(function(x,y) setCellValue(cells[[x,2]], column_description_1_id[y]), 18:25, 1:8)
      tmp_4 <- mapply(function(x,y) setCellValue(cells[[x,2]], column_description_2_id[y]), 26:33, 1:8)
      setColumnWidth(sheet1,2,25)                 
      tmp_5 = mapply(function(x,y) setCellValue(cells[[x,3]], column_description_1[y]), 2:9, 1:8);
      tmp_6 = mapply(function(x,y) setCellValue(cells[[x,3]], column_description_2[y]), 10:17, 1:8);
      tmp_7 = mapply(function(x,y) setCellValue(cells[[x,3]], column_description_1[y]), 18:25, 1:8);
      tmp_8 = mapply(function(x,y) setCellValue(cells[[x,3]], column_description_2[y]), 26:33, 1:8);
      setColumnWidth(sheet1,3,90)              
      tmp9 <- mapply(function(x) setCellStyle(cells[[x,2]],every_cell_frame), 2:33)
      tmp10 <- mapply(function(x) setCellStyle(cells[[x,3]],every_cell_frame), 2:33)
    } else {
      tmp0 <- mapply(function(x) setCellStyle(cells[[x,1]],every_cell_frame), 2:17)
      setCellValue(title_cell_1, Table_1_Name);
      setCellStyle(title_cell_1, title_cell_style);    
      title_cell_2 = getCells(getRows(sheet1, 10), 1)[[1]];
      setCellValue(title_cell_2, Table_2_Name);
      setCellStyle(title_cell_2, title_cell_style);          
      setColumnWidth(sheet1,1,30)
      tmp_1 <- mapply(function(x,y) setCellValue(cells[[x,2]], column_description_1_id[y]), 2:9, 1:8)
      tmp_2 <- mapply(function(x,y) setCellValue(cells[[x,2]], column_description_2_id[y]), 10:17, 1:8)
      setColumnWidth(sheet1,2,25)
      tmp_5 = mapply(function(x,y) setCellValue(cells[[x,3]], column_description_1[y]), 2:9, 1:8);
      tmp_6 = mapply(function(x,y) setCellValue(cells[[x,3]], column_description_2[y]), 10:17, 1:8);          
      setColumnWidth(sheet1,3,90)         
      tmp9 <- mapply(function(x) setCellStyle(cells[[x,2]],every_cell_frame), 2:17)
      tmp10 <- mapply(function(x) setCellStyle(cells[[x,3]],every_cell_frame), 2:17)
    } 
    
  } else {
    column_description_1 <- c("GO索引号，与Gene Ontology数据库的GO索引号一致","GO名称，与Gene Ontology数据库中GO的命名方式一致","GO所属的type","差异基因计数，表示所有差异基因中属于某一GO的差异基因数量","基因计数，表示数据库中属于某一GO的基因数量","富集度，若p值相同，富集度越大的GO，表示该GO受到实验的影响越大","p值，评估GO的显著性水平，(p<0.05表示GO具有显著性差异,用黄色标注)","误判率，对p值准确率的判断，对GO显著性水平的再判断")
    column_description_2 <- c(column_description_1[-c(4,5)],"基因索引号，与NCBI GenBank的基因索引号一致","基因名称，与NCBI GenBank的基因标识一致")

    if (nrow(up_diff_gene_list) != 0 & nrow(down_diff_gene_list) != 0) {
      tmp0 <- mapply(function(x) setCellStyle(cells[[x,1]],every_cell_frame), 2:33)
      setCellValue(title_cell_1, Table_1_Name);
      setCellStyle(title_cell_1, title_cell_style);    
      title_cell_2 = getCells(getRows(sheet1, 10), 1)[[1]];
      setCellValue(title_cell_2, Table_2_Name);
      setCellStyle(title_cell_2, title_cell_style);
      addMergedRegion(sheet1, 18, 25, 1, 1);title_cell_3 = getCells(getRows(sheet1, 18), 1)[[1]];
      setCellValue(title_cell_3, Table_3_Name);
      setCellStyle(title_cell_3, title_cell_style);
      addMergedRegion(sheet1, 26, 33, 1, 1);title_cell_4 = getCells(getRows(sheet1, 26), 1)[[1]];
      setCellValue(title_cell_4, Table_4_Name);
      setCellStyle(title_cell_4, title_cell_style);
      setColumnWidth(sheet1,1,30)
      tmp_1 <- mapply(function(x,y) setCellValue(cells[[x,2]], column_description_1_id[y]), 2:9, 1:8)
      tmp_2 <- mapply(function(x,y) setCellValue(cells[[x,2]], column_description_2_id[y]), 10:17, 1:8)
      tmp_3 <- mapply(function(x,y) setCellValue(cells[[x,2]], column_description_1_id[y]), 18:25, 1:8)
      tmp_4 <- mapply(function(x,y) setCellValue(cells[[x,2]], column_description_2_id[y]), 26:33, 1:8)
      setColumnWidth(sheet1,2,25)    
      tmp_5 = mapply(function(x,y) setCellValue(cells[[x,3]], column_description_1[y]), 2:9, 1:8);
      tmp_6 = mapply(function(x,y) setCellValue(cells[[x,3]], column_description_2[y]), 10:17, 1:8);
      tmp_7 = mapply(function(x,y) setCellValue(cells[[x,3]], column_description_1[y]), 18:25, 1:8);
      tmp_8 = mapply(function(x,y) setCellValue(cells[[x,3]], column_description_2[y]), 26:33, 1:8);
      setColumnWidth(sheet1,3,90)
      tmp9 <- mapply(function(x) setCellStyle(cells[[x,2]],every_cell_frame), 2:33)
      tmp10 <- mapply(function(x) setCellStyle(cells[[x,3]],every_cell_frame), 2:33)       
    } else {
      tmp0 <- mapply(function(x) setCellStyle(cells[[x,1]],every_cell_frame), 2:17)
      setCellValue(title_cell_1, Table_1_Name);
      setCellStyle(title_cell_1, title_cell_style);    
      title_cell_2 = getCells(getRows(sheet1, 10), 1)[[1]];
      setCellValue(title_cell_2, Table_2_Name);
      setCellStyle(title_cell_2, title_cell_style);          
      setColumnWidth(sheet1,1,30)
      tmp_1 <- mapply(function(x,y) setCellValue(cells[[x,2]], column_description_1_id[y]), 2:9, 1:8)
      tmp_2 <- mapply(function(x,y) setCellValue(cells[[x,2]], column_description_2_id[y]), 10:17, 1:8)
      setColumnWidth(sheet1,2,25)
      tmp_5 = mapply(function(x,y) setCellValue(cells[[x,3]], column_description_1[y]), 2:9, 1:8);
      tmp_6 = mapply(function(x,y) setCellValue(cells[[x,3]], column_description_2[y]), 10:17, 1:8);          
      setColumnWidth(sheet1,3,90)         
      tmp9 <- mapply(function(x) setCellStyle(cells[[x,2]],every_cell_frame), 2:17)
      tmp10 <- mapply(function(x) setCellStyle(cells[[x,3]],every_cell_frame), 2:17)
    }
  }    
  # Save the workbook into a file
  saveWorkbook(wb,File_Details_Directory);    
}


Append_XLSX_Result <- function(File_Name,Res_Data,Fill_Row_Number) {
  
  wb = loadWorkbook(paste(New_File_Directory,"_",Analysis_Method,".xlsx",sep=""))
  sheets = createSheet(wb, File_Name)
  Column_Header_Style <- CellStyle(wb,font = Font(wb,name="Times New Roman",isBold = TRUE))
  S_Font_Style <- CellStyle(wb,font = Font(wb,name="Times New Roman"),fill=Fill(foregroundColor = "YELLOW"))
  US_Font_Style <- CellStyle(wb,font = Font(wb,name="Times New Roman"))
  setColumnWidth(sheets,1,15); setColumnWidth(sheets,2,25); setColumnWidth(sheets,3:7,15)
  addDataFrame(Res_Data[1:Fill_Row_Number,],sheets,row.names=FALSE,colnamesStyle=Column_Header_Style,colStyle=list('1'=S_Font_Style,'2'=S_Font_Style,'3'=S_Font_Style,'4'=S_Font_Style,'5'=S_Font_Style,'6'=S_Font_Style,'7'=S_Font_Style,'8'=S_Font_Style))
  addDataFrame(Res_Data[(Fill_Row_Number+1):nrow(Res_Data),],sheets,col.names=FALSE,row.names=FALSE,startRow=Fill_Row_Number+2,colnamesStyle=Column_Header_Style,colStyle=list('1'=US_Font_Style,'2'=US_Font_Style,'3'=US_Font_Style,'4'=US_Font_Style,'5'=US_Font_Style,'6'=US_Font_Style,'7'=US_Font_Style,'8'=US_Font_Style))
  saveWorkbook(wb,paste(New_File_Directory,"_",Analysis_Method,".xlsx",sep=""))
}

# Append_XLSX_Result(File_Name=Table_1_Name,Res_Data=Ana_Result,Fill_Row_Number=sum(Ana_Result$pvalue < 0.05))
library(ggplot2)
ggplot2_image <- function(path_plot,path_index) {
  library(ggplot2)
  
  path_plot <- data.frame(cbind(rep(path_plot[,2],2),c(rep("pvalue",nrow(path_plot)),rep("FDR",nrow(path_plot))),c(-log10(path_plot[,6]),-log10(path_plot[,7]))),stringsAsFactors=FALSE)
  colnames(path_plot) <- c("path_name","legend","values")
  
  p <- ggplot(path_plot,aes(x=reorder(path_name,as.numeric(values)),y=values,fill=legend))+xlab(NULL)+ylab(NULL)+scale_y_discrete(breaks=range(path_plot$values),labels=signif(as.numeric(range(path_plot$values)),3))
  p+geom_bar(stat="identity",width=0.5,position="dodge")+labs(title="Significant Pathway")+coord_flip()
  ggsave(paste(New_File_Directory,"_",path_index,"_",Analysis_Method,".png",sep=""),width=10,height=20)
  #dev.off()
  
  wb = loadWorkbook(paste(New_File_Directory,"_",Analysis_Method,".xlsx",sep=""))
  sheets = createSheet(wb,paste(path_index,Analysis_Method,"Image",sep=""))
  addPicture(paste(New_File_Directory,"_",Analysis_Method,".png",sep=""),sheets,startRow=3, startColumn=2)
  saveWorkbook(wb,paste(New_File_Directory,"_",Analysis_Method,".xlsx",sep=""))
}




#set your paramters
Go_picture <- function(path,method = "Up_Down Gene"){

    path_GO = paste0(path,"/output")
    #path_Pathway = paste0(path,"/output")
    filename = list.files(paste0(path,"/output"))
    
    GO_file_name = filename [grepl("GO-analysis.xlsx", filename )]
    #Path_file_name = filename [grepl("Path-analysis.xlsx", filename )]
    num_go = length(GO_file_name)
    samlpenames= gsub(pattern = "_GO-analysis.xlsx", replacement = "", x = GO_file_name)
        for ( i in 1:num_go){

  if(method == "Up_Down Gene"|method == "Target Up_Down Gene"){
    #Ctrl+A run
    ############################################################
    #GO-Analysis
    setwd(path_GO)
    library("xlsx")
    go_upGene_fun =  read.xlsx(GO_file_name[i],sheetIndex = 2)
        
    go_downGene_fun =  read.xlsx(GO_file_name[i],sheetIndex = 4)
        
    sheet1 = go_upGene_fun[which(go_upGene_fun$pvalue<0.05),c("go_name","pvalue")]
    sheet3 = go_downGene_fun[which(go_downGene_fun$pvalue<0.05),c("go_name","pvalue")] 
    sheet2_up = go_upGene_fun[which(go_upGene_fun$pvalue<0.05),c("go_name","enrichment")]
    sheet2_down = go_downGene_fun[which(go_downGene_fun$pvalue<0.05),c("go_name","enrichment")]
    
    
    
    if (nrow(sheet1)>=30 & nrow(sheet3)>=30){
      sheet1$pvalue = -log10(sheet1$pvalue) 
      colnames(sheet1) = c("go_name","(-LgP)")
      sheet1_plot = sheet1[1:30,]
      colnames(sheet1_plot) = c("go_name","LgP")
      
      
      
      sheet3$pvalue = -log10(sheet3$pvalue) 
      colnames(sheet3) = c("go_name","(-LgP)")
      sheet3_plot = sheet3[1:30,]
      colnames(sheet3_plot) = c("go_name","LgP")
      
      
      
      sheet2_down$enrichment = -sheet2_down$enrichment
      sheet2_up_plot = sheet2_up[1:30,]
      sheet2_down_plot = sheet2_down[1:30,]
      
      
      
    }else if (nrow(sheet1)>=30 & nrow(sheet3) < 30){
      sheet1$pvalue = -log10(sheet1$pvalue) 
      colnames(sheet1) = c("go_name","(-LgP)")
      sheet1_plot = sheet1[1:30,]
      colnames(sheet1_plot) = c("go_name","LgP")
      
      
      
      sheet3$pvalue = -log10(sheet3$pvalue) 
      colnames(sheet3) = c("go_name","(-LgP)")
      sheet3_plot = sheet3
      colnames(sheet3_plot) = c("go_name","LgP")
      
      
      
      sheet2_down$enrichment = -sheet2_down$enrichment
      sheet2_up_plot = sheet2_up[1:30,]
      sheet2_down_plot = sheet2_down
      
      
    }else if (nrow(sheet1)<30 & nrow(sheet3) >= 30){
      sheet1$pvalue = -log10(sheet1$pvalue) 
      colnames(sheet1) = c("go_name","(-LgP)")
      sheet1_plot = sheet1
      colnames(sheet1_plot) = c("go_name","LgP")
      
      
      
      sheet3$pvalue = -log10(sheet3$pvalue) 
      colnames(sheet3) = c("go_name","(-LgP)")
      sheet3_plot = sheet3[1:30,]
      colnames(sheet3_plot) = c("go_name","LgP")
      
      
      
      sheet2_down$enrichment = -sheet2_down$enrichment
      sheet2_up_plot = sheet2_up
      sheet2_down_plot = sheet2_down[1:30,]
    }else{
      sheet1$pvalue = -log10(sheet1$pvalue) 
      colnames(sheet1) = c("go_name","(-LgP)")
      sheet1_plot = sheet1
      colnames(sheet1_plot) = c("go_name","LgP")
      
      
      
      sheet3$pvalue = -log10(sheet3$pvalue) 
      colnames(sheet3) = c("go_name","(-LgP)")
      sheet3_plot = sheet3
      colnames(sheet3_plot) = c("go_name","LgP")
      
      
      
      sheet2_down$enrichment = -sheet2_down$enrichment
      sheet2_up_plot = sheet2_up
      sheet2_down_plot = sheet2_down
    }
    
    sheet2_down_plot$enrichment = -sheet2_down_plot$enrichment
    
    #sheet2 = rbind(sheet2_up,sheet2_down)
    #sheet2 = sheet2[order(-sheet2$enrichment),]
    
    #file <- paste(path, "/GO-analysis-picture2.xlsx", sep="")
    #write.xlsx(sheet1,file ,sheetName = "Go-(-LgP) of significant up-regulated genes",row.names=FALSE)
    #write.xlsx(sheet2,file ,sheetName = "Go Enrichment",row.names=FALSE,append = TRUE)
    #write.xlsx(sheet3,file ,sheetName = "Go-(-LgP) of significant down-regulated genes",row.names=FALSE,append = TRUE)
    
    
    #plot 4 pictures
    
    #
    sheet1_plot = sheet1_plot[order(sheet1_plot$LgP),]
    sheet3_plot = sheet3_plot[order(sheet3_plot$LgP),]
    sheet2_up_plot = sheet2_up_plot[order(sheet2_up_plot$enrichment),]
    sheet2_down_plot = sheet2_down_plot[order(sheet2_down_plot$enrichment),]
    
    #pdf("mutation_type.pdf")
     
    
    ggplot(sheet1_plot, aes(x = go_name,y=LgP)) + 
      geom_bar(stat = "identity",fill ="#007bba")+ scale_x_discrete(limits=sheet1_plot$go_name)+
      labs(title = "Difgene Sig Go", y ="(-LgP)",x="Gene ontology category")+
      coord_flip()
    
    
    ggsave(paste0(samlpenames[i],"UPgene_Sig_Go.pdf"),width = 9.18, height = 6.53)
    #ggsave("UPgene_Sig_Go.png",width = 9.18, height = 6.53)
    
    
    
    ggplot(sheet3_plot, aes(x = go_name,y=LgP)) + 
      geom_bar(stat = "identity",fill ="#007bba")+ scale_x_discrete(limits=sheet3_plot$go_name)+
      labs(title = "Difgene Sig Go", y ="(-LgP)",x="Gene ontology category")+
      coord_flip()
    
    
    ggsave(paste0(samlpenames[i],"Downgene_Sig_Go.pdf"),width = 9.18, height = 6.53)
    #ggsave("Downgene_Sig_Go.png",width = 9.18, height = 6.53)
    
    
    ggplot(sheet2_up_plot, aes(x = go_name,y = enrichment)) + 
      geom_bar(stat = "identity",fill ="#007bba")+ scale_x_discrete(limits=sheet2_up_plot$go_name)+
      labs(title = "Difgene Sig Go", y ="enrichment",x="Gene ontology category")+
      coord_flip()
    ggsave(paste0(samlpenames[i],"Upgene_Sig_GoEnrichment.pdf"),width = 9.18, height = 6.53)
    #ggsave("Upgene_Sig_GoEnrichment.png",width = 9.18, height = 6.53)
    
    ggplot(sheet2_down_plot, aes(x = go_name,y = enrichment)) + 
      geom_bar(stat = "identity",fill ="#007bba")+ scale_x_discrete(limits=sheet2_down_plot$go_name)+
      labs(title = "Difgene Sig Go", y ="enrichment",x="Gene ontology category")+
      coord_flip()
    ggsave(paste0(samlpenames[i],"Downgene_Sig_GoEnrichment.pdf"),width = 9.18, height = 6.53)
   # ggsave("Downgene_Sig_GoEnrichment.png",width = 9.18, height = 6.53)
    
    
    
  }else if(method == "Target Differential Gene"|method =="Differential Gene"){
        
    
    #Ctrl+A run
    ############################################################
    #GO-Analysis
    setwd(path_GO)
    library("xlsx")
    go_upGene_fun =  read.xlsx(GO_file_name[i],sheetIndex = 2)
    #go_Fun_upGene =  read.xlsx("GO-analysis.xlsx",sheetIndex = 3)
    
    #go_downGene_fun =  read.xlsx(GO_file_name,sheetIndex = 3)
     
    sheet1 = go_upGene_fun[which(go_upGene_fun$pvalue<0.05),c("go_name","pvalue")]
    #sheet3 = go_downGene_fun[which(go_downGene_fun$pvalue<0.05),c(2,6)] 
    sheet2_up = go_upGene_fun[which(go_upGene_fun$pvalue<0.05),c("go_name","enrichment")]
    #sheet2_down = go_downGene_fun[which(go_downGene_fun$pvalue<0.05),c(2,5)]
    
    
    
    if (nrow(sheet1)>=30 ){
      sheet1$pvalue = -log10(sheet1$pvalue) 
      colnames(sheet1) = c("go_name","(-LgP)")
      sheet1_plot = sheet1[1:30,]
      colnames(sheet1_plot) = c("go_name","LgP")

      #sheet2_down$enrichment = -sheet2_down$enrichment
      sheet2_up_plot = sheet2_up[1:30,]
      #sheet2_down_plot = sheet2_down[1:30,]
      
      
      
    
    }else{
      sheet1$pvalue = -log10(sheet1$pvalue) 
      colnames(sheet1) = c("go_name","(-LgP)")
      sheet1_plot = sheet1
      colnames(sheet1_plot) = c("go_name","LgP")
      
      
      
    
      
      #sheet2_down$enrichment = -sheet2_down$enrichment
      sheet2_up_plot = sheet2_up
      #sheet2_down_plot = sheet2_down
    }
    
    #sheet2_down_plot$enrichment = -sheet2_down_plot$enrichment
    
    #sheet2 = rbind(sheet2_up,sheet2_down)
    #sheet2 = sheet2[order(-sheet2$enrichment),]
    
    #file <- paste(path, "/GO-analysis-picture2.xlsx", sep="")
    #write.xlsx(sheet1,file ,sheetName = "Go-(-LgP) of significant up-regulated genes",row.names=FALSE)
    #write.xlsx(sheet2,file ,sheetName = "Go Enrichment",row.names=FALSE,append = TRUE)
    #write.xlsx(sheet3,file ,sheetName = "Go-(-LgP) of significant down-regulated genes",row.names=FALSE,append = TRUE)
    
    
    #plot 4 pictures
    
    #
    sheet1_plot = sheet1_plot[order(sheet1_plot$LgP),]
    #sheet3_plot = sheet3_plot[order(sheet3_plot$LgP),]
    sheet2_up_plot = sheet2_up_plot[order(sheet2_up_plot$enrichment),]
    #sheet2_down_plot = sheet2_down_plot[order(sheet2_down_plot$enrichment),]
    
    #pdf("mutation_type.pdf")
    
    
    ggplot(sheet1_plot, aes(x = go_name,y=LgP)) + 
      geom_bar(stat = "identity",fill ="#007bba")+ scale_x_discrete(limits=sheet1_plot$go_name)+
      labs(title = "Difgene Sig Go", y ="(-LgP)",x="Gene ontology category")+
      coord_flip()
    
    ggsave(paste0(samlpenames[i],"Diff_gene_Sig_Go.pdf"),width = 9.18, height = 6.53)
    #ggsave("Diff_gene_Sig_Go.png",width = 9.18, height = 6.53)
    
    

    ggplot(sheet2_up_plot, aes(x = go_name,y = enrichment)) + 
      geom_bar(stat = "identity",fill ="#007bba")+ scale_x_discrete(limits=sheet2_up_plot$go_name)+
      labs(title = "Difgene Sig Go", y ="enrichment",x="Gene ontology category")+
      coord_flip()
    ggsave(paste0(samlpenames[i],"Diff_gene_Sig_GoEnrichment.pdf"),width = 9.18, height = 6.53)
      
    
  }else{
    print("please input the right source_method")
  }



        }




  
  
  
  
  
  ##########################################
  
  
  
}



#Pathway-analysis

Path_picture <- function(path,method = "Up_Down Gene"){


    
    path_Pathway = paste0(path,"/output")
    filename = list.files(paste0(path,"/output"))
    Path_file_name = filename [grepl("Path-analysis.xlsx", filename )]

    num_path =  length(Path_file_name)
    sample_path=gsub(pattern = "_Path-analysis.xlsx", replacement = "", x = Path_file_name)

   for ( j in 1:num_path){

if(method == "Up_Down Gene"|method == "Target Up_Down Gene"){
 
  setwd(path_Pathway)
    library("xlsx")
    path_upGene_fun =  read.xlsx(Path_file_name[j],sheetIndex = 2)
       
    path_downGene_fun =  read.xlsx(Path_file_name[j],sheetIndex = 4)
        
    sheet1 = path_upGene_fun[which(path_upGene_fun$pvalue<0.05),c("path_name","pvalue")]
    sheet3 = path_downGene_fun[which(path_downGene_fun$pvalue<0.05),c("path_name","pvalue")] 
    
    
    
    if (nrow(sheet1)>=30 & nrow(sheet3)>=30){
      sheet1$pvalue = -log10(sheet1$pvalue) 
      colnames(sheet1) = c("path_name","(-LgP)")
      sheet1_plot = sheet1[1:30,]
      colnames(sheet1_plot) = c("path_name","LgP")
      
      
      
      sheet3$pvalue = -log10(sheet3$pvalue) 
      colnames(sheet3) = c("path_name","(-LgP)")
      sheet3_plot = sheet3[1:30,]
      colnames(sheet3_plot) = c("path_name","LgP")
      
      
      
      
    }else if (nrow(sheet1)>=30 & nrow(sheet3) < 30){
      sheet1$pvalue = -log10(sheet1$pvalue) 
      colnames(sheet1) = c("path_name","(-LgP)")
      sheet1_plot = sheet1[1:30,]
      colnames(sheet1_plot) = c("path_name","LgP")
      
      
      
      sheet3$pvalue = -log10(sheet3$pvalue) 
      colnames(sheet3) = c("path_name","(-LgP)")
      sheet3_plot = sheet3
      colnames(sheet3_plot) = c("path_name","LgP")
      
      
      
    }else if (nrow(sheet1)<30 & nrow(sheet3) >= 30){
      sheet1$pvalue = -log10(sheet1$pvalue) 
      colnames(sheet1) = c("path_name","(-LgP)")
      sheet1_plot = sheet1
      colnames(sheet1_plot) = c("path_name","LgP")
      
      
      
      sheet3$pvalue = -log10(sheet3$pvalue) 
      colnames(sheet3) = c("path_name","(-LgP)")
      sheet3_plot = sheet3[1:30,]
      colnames(sheet3_plot) = c("path_name","LgP")
      
      
    }else{
      sheet1$pvalue = -log10(sheet1$pvalue) 
      colnames(sheet1) = c("path_name","(-LgP)")
      sheet1_plot = sheet1
      colnames(sheet1_plot) = c("path_name","LgP")
      
      
      
      sheet3$pvalue = -log10(sheet3$pvalue) 
      colnames(sheet3) = c("path_name","(-LgP)")
      sheet3_plot = sheet3
      colnames(sheet3_plot) = c("path_name","LgP")
      
      
      
    }
    
    
    
    
    #
    sheet1_plot = sheet1_plot[order(sheet1_plot$LgP),]
    sheet3_plot = sheet3_plot[order(sheet3_plot$LgP),]
    
    ggplot(sheet1_plot, aes(x = path_name,y=LgP)) + 
      geom_bar(stat = "identity",fill ="#007bba")+ scale_x_discrete(limits=sheet1_plot$path_name)+
      labs(title = "Difgene Sig Pathway", y ="(-LgP)",x="Pathway")+
      coord_flip()
    

    ggsave(paste0(sample_path[j],"UPgene_Sig_Pathway.pdf"),width = 9.18, height = 6.53)
   # ggsave("UPgene_Sig_Pathway.png",width = 9.18, height = 6.53)
    
    
    
    ggplot(sheet3_plot, aes(x = path_name,y=LgP)) + 
      geom_bar(stat = "identity",fill ="#007bba")+ scale_x_discrete(limits=sheet3_plot$path_name)+
      labs(title = "Difgene Sig Pathway", y ="(-LgP)",x="Pathway")+
      coord_flip()
    
    
    ggsave(paste0(sample_path[j],"Downgene_Sig_Pathway.pdf"),width = 9.18, height = 6.53)
   # ggsave("Downgene_Sig_Pathway.png",width = 9.18, height = 6.53)
    
    
    
  }else if (method == "Target Differential Gene"|method == "Differential Gene") {
    
       
    setwd(path_Pathway)
    library("xlsx")
    path_upGene_fun =  read.xlsx(Path_file_name[j],sheetIndex = 2)
       
    #path_downGene_fun =  read.xlsx(Path_file_name,sheetIndex = 4)
    #go_Fun_downGene =  read.xlsx("GO-analysis.xlsx",sheetIndex = 5)
    
    sheet1 = path_upGene_fun[which(path_upGene_fun$pvalue<0.05),c("path_name","pvalue")]
    #sheet3 = path_downGene_fun[which(path_downGene_fun$pvalue<0.05),c(2,6)] 
    
    
    
    if (nrow(sheet1)>=30){
      sheet1$pvalue = -log10(sheet1$pvalue) 
      colnames(sheet1) = c("path_name","(-LgP)")
      sheet1_plot = sheet1[1:30,]
      colnames(sheet1_plot) = c("path_name","LgP")
      
      
      
      #sheet3$pvalue = -log10(sheet3$pvalue) 
      #colnames(sheet3) = c("path_name","(-LgP)")
      #sheet3_plot = sheet3[1:30,]
     # colnames(sheet3_plot) = c("path_name","LgP")
      
      
      
      
    }else{
      sheet1$pvalue = -log10(sheet1$pvalue) 
      colnames(sheet1) = c("path_name","(-LgP)")
      sheet1_plot = sheet1
      colnames(sheet1_plot) = c("path_name","LgP")
      

    }
    
    
    
    
    #
    sheet1_plot = sheet1_plot[order(sheet1_plot$LgP),]
    #sheet3_plot = sheet3_plot[order(sheet3_plot$LgP),]
    
    ggplot(sheet1_plot, aes(x = path_name,y=LgP)) + 
      geom_bar(stat = "identity",fill ="#007bba")+ scale_x_discrete(limits=sheet1_plot$path_name)+
      labs(title = "Difgene Sig Pathway", y ="(-LgP)",x="Pathway")+
      coord_flip()
    
    
    ggsave(paste0(sample_path[j],"Diff_gene_Sig_Pathway.pdf"), width = 9.18, height = 6.53)
   # ggsave("Diff_gene_Sig_Pathway.png",width = 9.18, height = 6.53)
    
    
    

    
  }else{
    print("please input right source_method")
  }
  
  






     }




  
  
  
  
  
  
}




