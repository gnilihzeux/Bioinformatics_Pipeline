appendXlsxResult <- function(Go_pathway_result, workdir, Analysis_Method="up_gene-analysis") {
  
  if (!require(xlsx)) {
    install.packages("xlsx")
  } else {
    
    #' @issue how to suppress message
    require(xlsx)
  }
  
  
  
  #sheet1 文件说明
  wb <- createWorkbook()
  #creat a rew sheet with name "文件说明"
  sheet1 = createSheet(wb, "文件说明") 
  #这个张表22行３列
  cells = createCell(createRow(sheet1, 1:22), 1:3)
  Header_head_doc <- c("表名称","列头","列头的含义")
  Column_Header_Style <- CellStyle(wb,
                                   alignment = Alignment(h="ALIGN_CENTER"),
                                   font = Font(wb,isBold = TRUE),
                                   border = Border(position=c("BOTTOM","LEFT", "TOP", "RIGHT"),pen="BORDER_THIN"),
                                   fill = Fill(backgroundColor="gray"))
  #设置列头分别的文本和style
  tmp_0 = mapply(function(x) setCellValue(cells[[1,x]], Header_head_doc[x]), 1:3)    
  tmp_0 = mapply(function(x) setCellStyle(cells[[1,x]],Column_Header_Style), 1:3)
  
  
  # Merge the first column into one cell
  addMergedRegion(sheet1, 2, 12, 1, 1);
  addMergedRegion(sheet1, 13, 22, 1, 1);
  # Create the style for title cell
  title_cell_style = CellStyle(wb,
                               alignment = Alignment(h="ALIGN_CENTER",v="VERTICAL_CENTER"),
                               font = Font(wb, color="blue",name="Times New Roman",isBold = TRUE));
  
  
  title_cell_1 = getCells(getRows(sheet1, 2), 1)[[1]];
  every_cell_frame <- CellStyle(wb,font=Font(wb,name="Times New Roman"),border = Border(position=c("BOTTOM","LEFT", "TOP", "RIGHT"),pen="BORDER_THIN"))
  
  if (Analysis_Method == "up_gene-analysis") {
    
    File_Name=c("上调基因GO-BP富集","上调基因GO-CC富集","上调基因GO-MF富集","上调基因KEGG富集","上调基因KEGG富集")
    column1_description_1 = c("上调基因GO富集")
    column1_description_2 = c("上调基因KEGG富集")
    
    column2_description_1 = c("goID",
                              "goDescription",
                              "goType",
                              "geneRatio",
                              "bgRatio",
                              "pvalue",
                              "padj",
                              "qvalue",
                              "enrichScore",
                              "overlapGeneList",
                              "overlapCount")
    column2_description_2 = c("pathID",
                              "pathDescription",
                              "geneRatio",
                              "bgRatio",
                              "pvalue",
                              "padj",
                              "qvalue",
                              "enrichScore",
                              "overlapGeneList",
                              "overlapCount")
    
    
    
    column3_description_1 = c("GO term 索引号，与Gene Ontology数据库的GO索引号一致",
                              "GO term 的描述信息",
                              "GO term 的分类，包括 Biological Pathway/Cell Component/Molecular Function",
                              "富集到 GO term 中的基因占所有表达上调基因的比例，展示形式为：“富集基因数目/上调基因数目”",
                              "background ratio，GO term 基因数目与背景基因数目的比率，表示为“GO term 基因数目/背景基因数目”",
                              "p 值，评估上调基因富集到 GO term 的统计显著性水平",
                              "p 值校正值 FDR，即对相同条件下的试验的差异显著性进行多重检验校正，常用的方法包括 Bonfferoni，BH等",
                              "q 值是 FDR 值的一种扩展，也可以简单地理解为另外一种 FDR 计算方法",
                              "上调基因在GO term中的富集分数",
                              "上调基因集合与 GO term 基因集合的交集",
                              "上调基因集合与 GO term 基因集合交集的基因数目")
    column3_description_2 = c("KEGG term 索引号，与KEGG数据库的索引号一致",
                              "KEGG term 的描述信息",
                              "富集到 KEGG term 中的基因占所有表达上调基因的比例，展示形式为：“富集基因数目/上调基因数目”",
                              "background ratio，KEGG term 基因数目与背景基因数目的比率，表示为“KEGG term 基因数目/背景基因数目”",
                              "p 值，评估上调基因富集到 KEGG term 的统计显著性水平",
                              "p 值校正值 FDR，即对相同条件下的试验的差异显著性进行多重检验校正，常用的方法包括 Bonfferoni，BH等",
                              "q 值是 FDR 值的一种扩展，也可以简单地理解为另外一种 FDR 计算方法",
                              "上调基因集合与 KEGG term 基因集合的交集",
                              "上调基因在KEGG中的富集分析",
                              "上调基因集合与 KEGG term 基因集合交集的基因数目")
    
    
    
  }else{
    File_Name=c("下调基因GO-BP富集","下调基因GO-CC富集","下调基因GO-MF富集","下调基因KEGG富集","下调基因KEGG富集")
    
    column1_description_1 = c("下调基因GO富集")
    column1_description_2 = c("下调基因KEGG富集")
    
    column2_description_1 = c("goID",
                              "goDescription",
                              "goType",
                              "geneRatio",
                              "bgRatio",
                              "pvalue",
                              "padj",
                              "qvalue",
                              "enrichScore",
                              "overlapGeneList",
                              "overlapCount")
    column2_description_2 = c("pathID",
                              "pathDescription",
                              "geneRatio",
                              "bgRatio",
                              "pvalue",
                              "padj",
                              "qvalue",
                              "enrichScore",
                              "overlapGeneList",
                              "overlapCount")
    
    
    
    column3_description_1 = c("GO term 索引号，与Gene Ontology数据库的GO索引号一致",
                              "GO term 的描述信息",
                              "GO term 的分类，包括 Biological Pathway/Cell Component/Molecular Function",
                              "富集到 GO term 中的基因占所有表达上调基因的比例，展示形式为：“富集基因数目/下调基因数目”",
                              "background ratio，GO term 基因数目与背景基因数目的比率，表示为“GO term 基因数目/背景基因数目”",
                              "p 值，评估上调基因富集到 GO term 的统计显著性水平",
                              "p 值校正值 FDR，即对相同条件下的试验的差异显著性进行多重检验校正，常用的方法包括 Bonfferoni，BH等",
                              "q 值是 FDR 值的一种扩展，也可以简单地理解为另外一种 FDR 计算方法",
                              "下调基因在GO term中的富集分数",
                              "下调基因集合与 GO term 基因集合的交集",
                              "下调基因集合与 GO term 基因集合交集的基因数目")
    column3_description_2 = c("KEGG term 索引号，与KEGG数据库的索引号一致",
                              "KEGG term 的描述信息",
                              "富集到 KEGG term 中的基因占所有表达上调基因的比例，展示形式为：“富集基因数目/下调基因数目”",
                              "background ratio，KEGG term 基因数目与背景基因数目的比率，表示为“KEGG term 基因数目/背景基因数目”",
                              "p 值，评估上调基因富集到 KEGG term 的统计显著性水平",
                              "p 值校正值 FDR，即对相同条件下的试验的差异显著性进行多重检验校正，常用的方法包括 Bonfferoni，BH等",
                              "q 值是 FDR 值的一种扩展，也可以简单地理解为另外一种 FDR 计算方法",
                              "下调基因集合与 KEGG term 基因集合的交集",
                              "下调基因在KEGG中的富集分析",
                              "下调基因集合与 KEGG term 基因集合交集的基因数目")
    
    
  }
  
  
  
  tmp0 <- mapply(function(x) setCellStyle(cells[[x,1]],every_cell_frame), 2:22)
  
  setCellValue(title_cell_1, column1_description_1);
  setCellStyle(title_cell_1, title_cell_style);          
  title_cell_2 = getCells(getRows(sheet1, 13), 1)[[1]];
  setCellValue(title_cell_2, column1_description_2);
  setCellStyle(title_cell_2, title_cell_style);
  
  tmp1 <- mapply(setCellValue, cells[2:22,2], c(column2_description_1,column2_description_2))
  tmp2 <- mapply(setCellValue, cells[2:22,3], c(column3_description_1,column3_description_2))
  tmp3 <- mapply(function(x) setCellStyle(cells[[x,2]],every_cell_frame), 2:22)
  tmp4 <- mapply(function(x) setCellStyle(cells[[x,3]],every_cell_frame), 2:22)
  
  # 设置列宽
  setColumnWidth(sheet1,1,30)
  setColumnWidth(sheet1,2,25)                 
  setColumnWidth(sheet1,3,90)  
  
  for(i in 1:4){
    sheets = createSheet(wb, File_Name[i])
    Column_Header_Style = CellStyle(wb,font = Font(wb,name="Times New Roman",isBold = TRUE))
    S_Font_Style = CellStyle(wb,font = Font(wb,name="Times New Roman"),fill=Fill(foregroundColor = "YELLOW"))
    US_Font_Style = CellStyle(wb,font = Font(wb,name="Times New Roman"))
    setColumnWidth(sheets,1,15); 
    setColumnWidth(sheets,2,25); 
    setColumnWidth(sheets,3:11,15);
    Fill_Row_Number=sum(Go_pathway_result[[i]]$pvalue < 0.05)
	Res_Data <- Go_pathway_result[[i]]
    addDataFrame(Res_Data[1:Fill_Row_Number,],
                 sheets,row.names=FALSE,
                 colnamesStyle=Column_Header_Style,
                 colStyle=list('1'=S_Font_Style,
                               '2'=S_Font_Style,
                               '3'=S_Font_Style,
                               '4'=S_Font_Style,
                               '5'=S_Font_Style,
                               '6'=S_Font_Style,
                               '7'=S_Font_Style,
                               '8'=S_Font_Style,
                               '9'=S_Font_Style,
                               '10'=S_Font_Style,
                               '11'=S_Font_Style))
    addDataFrame(Res_Data[(Fill_Row_Number+1):nrow(Go_pathway_result[[i]]),],
                 sheets,
                 col.names=FALSE,
                 row.names=FALSE,
                 startRow=Fill_Row_Number+2,
                 colnamesStyle=Column_Header_Style,
                 colStyle=list('1'=US_Font_Style,
                               '2'=US_Font_Style,
                               '3'=US_Font_Style,
                               '4'=US_Font_Style,
                               '5'=US_Font_Style,
                               '6'=US_Font_Style,
                               '7'=US_Font_Style,
                               '8'=US_Font_Style,
                               '9'=US_Font_Style,
                               '10'=US_Font_Style,
                               '11'=US_Font_Style))
  }
  
  workdir <- sub("/$", "", workdir)
  if(Analysis_Method == "up_gene-analysis"){
    
    saveWorkbook(wb,paste0(workdir,"/up_gene_goPathway.xlsx")); 
  }else{
    saveWorkbook(wb,paste0(workdir,"/down_gene_goPathway.xlsx")); 
  }
  
}
