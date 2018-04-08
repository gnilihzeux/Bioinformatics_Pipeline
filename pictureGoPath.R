# GO/pathway plots ---------------------
pictureGoPath <- function(inputDir){
  #
  #
  library("xlsx")
  library(ggplot2)
  inputDir <- sub("/$", "", inputDir)
  
  go_upGene_bp =  read.xlsx(paste0(inputDir, "/up_gene_goPathway.xlsx"),sheetIndex = 2)
  go_upGene_cc =  read.xlsx(paste0(inputDir, "/up_gene_goPathway.xlsx"),sheetIndex = 3)
  go_upGene_mf =  read.xlsx(paste0(inputDir, "/up_gene_goPathway.xlsx"),sheetIndex = 4)
  go_downGene_bp =  read.xlsx(paste0(inputDir, "/down_gene_goPathway.xlsx"),sheetIndex = 2)
  go_downGene_cc =  read.xlsx(paste0(inputDir, "/down_gene_goPathway.xlsx"),sheetIndex = 3)
  go_downGene_mf =  read.xlsx(paste0(inputDir, "/down_gene_goPathway.xlsx"),sheetIndex = 4)
  
  
  sheet1 = rbind(go_upGene_bp[which(go_upGene_bp$pvalue<0.05),c("goDescription","goType","pvalue")],go_upGene_cc[which(go_upGene_cc$pvalue<0.05),c("goDescription","goType","pvalue")],go_upGene_mf[which(go_upGene_mf$pvalue<0.05),c("goDescription","goType","pvalue")])
  
  sheet3 = rbind(go_downGene_bp[which(go_downGene_bp$pvalue<0.05),c("goDescription","goType","pvalue")],go_downGene_cc[which(go_downGene_cc$pvalue<0.05),c("goDescription","goType","pvalue")],go_downGene_mf[which(go_downGene_mf$pvalue<0.05),c("goDescription","goType","pvalue")])
  
  sheet2_up = rbind(go_upGene_bp[which(go_upGene_bp$pvalue<0.05),c("goDescription","goType","enrichScore")],go_upGene_cc[which(go_upGene_cc$pvalue<0.05),c("goDescription","goType","enrichScore")],go_upGene_mf[which(go_upGene_mf$pvalue<0.05),c("goDescription","goType","enrichScore")])
  
  sheet2_down = rbind(go_downGene_bp[which(go_downGene_bp$pvalue<0.05),c("goDescription","goType","enrichScore")],go_downGene_cc[which(go_downGene_cc$pvalue<0.05),c("goDescription","goType","enrichScore")],go_downGene_mf[which(go_downGene_mf$pvalue<0.05),c("goDescription","goType","enrichScore")])
  
  
  #sheet3 = go_downGene_fun[which(go_downGene_fun$pvalue<0.05),c("go_name","goType","pvalue")] 
  #sheet2_up = go_upGene_fun[which(go_upGene_fun$pvalue<0.05),c("go_name","type","enrichScore")]
  #sheet2_down = go_downGene_fun[which(go_downGene_fun$pvalue<0.05),c("go_name","type","enrichment")]
  
  sheet1$pvalue = -log10(sheet1$pvalue)
  sheet3$pvalue = -log10(sheet3$pvalue)
  colnames(sheet1) = c("go_name","type","(-LgP)")
  colnames(sheet3) = c("go_name","type","(-LgP)")
  colnames(sheet2_up) = c("go_name","type","enrichment")
  colnames(sheet2_down) = c("go_name","type","enrichment")
  
  plot_sheet<-function(Go_data,enrichment_data,...){
    
    n_cc=length(which(Go_data$type=="cellular component"))
    n_mf=length(which(Go_data$type=="molecular function"))
    n_bp=length(which(Go_data$type=="biological pathway"))
    
    
    if (n_cc>=10 & n_mf>=10 & n_bp>=10){
      gene_Go = rbind(Go_data[which(Go_data$type=="cellular component"),][1:10,],Go_data[which(Go_data$type=="molecular function"),][1:10,],Go_data[which(Go_data$type=="biological pathway"),][1:10,])
    }else if (n_cc>=10 & n_mf>=10 & n_bp<10){
      gene_Go = rbind(Go_data[which(Go_data$type=="cellular component"),][1:10,],Go_data[which(Go_data$type=="molecular function"),][1:10,],Go_data[which(Go_data$type=="biological pathway"),])
      
    }else if (n_cc>=10 & n_mf<10 & n_bp>=10) {
      gene_Go = rbind(Go_data[which(Go_data$type=="cellular component"),][1:10,],Go_data[which(Go_data$type=="molecular function"),],Go_data[which(Go_data$type=="biological pathway"),][1:10,])
      
      
    }else if (n_cc<10 & n_mf>=10 & n_bp>=10) {
      
      gene_Go = rbind(Go_data[which(Go_data$type=="cellular component"),],Go_data[which(Go_data$type=="molecular function"),][1:10,],Go_data[which(Go_data$type=="biological pathway"),][1:10,])
      
    }else if (n_cc<10 & n_mf<10 & n_bp>=10){
      gene_Go = rbind(Go_data[which(Go_data$type=="cellular component"),],Go_data[which(Go_data$type=="molecular function"),],Go_data[which(Go_data$type=="biological pathway"),][1:10,])
      
      
    }else if (n_cc<10 & n_mf>=10 & n_bp<10){
      gene_Go = rbind(Go_data[which(Go_data$type=="cellular component"),],Go_data[which(Go_data$type=="molecular function"),][1:10,],Go_data[which(Go_data$type=="biological pathway"),])
      
      
      
    }else if (n_cc>=10 & n_mf<10 & n_bp<10){
      gene_Go = rbind(Go_data[which(Go_data$type=="cellular component"),][1:10,],Go_data[which(Go_data$type=="molecular function"),],Go_data[which(Go_data$type=="biological pathway"),])
      
      
    }else{
      gene_Go = rbind(Go_data[which(Go_data$type=="cellular component"),],Go_data[which(Go_data$type=="molecular function"),],Go_data[which(Go_data$type=="biological pathway"),])
      
      
    }
    
    colnames(gene_Go) = c("go_name","type","LgP")
    gene_enrichment=enrichment_data[which(enrichment_data$go_name %in% gene_Go$go_name),]
    return(list(gene_Go,gene_enrichment))
  }
  
  
  up_gene=plot_sheet(Go_data=sheet1,enrichment_data=sheet2_up)
  down_gene=plot_sheet(Go_data=sheet3,enrichment_data=sheet2_down)
  sheet1_plot=up_gene[[1]]
  sheet2_up_plot=up_gene[[2]]
  sheet3_plot=down_gene[[1]]
  sheet2_down_plot=down_gene[[2]]
  
  
  
  sheet2_up_bp=sheet2_up_plot[sheet2_up_plot$type=="biological pathway",]
  sheet2_up_cc=sheet2_up_plot[sheet2_up_plot$type=="cellular component",]
  sheet2_up_mf=sheet2_up_plot[sheet2_up_plot$type=="molecular function",]
  sheet2_up_bp= sheet2_up_bp[order(sheet2_up_bp$enrichment),]
  sheet2_up_cc= sheet2_up_cc[order(sheet2_up_cc$enrichment),]
  sheet2_up_mf= sheet2_up_mf[order(sheet2_up_mf$enrichment),]
  sheet2_up_plot=rbind(sheet2_up_bp,sheet2_up_cc,sheet2_up_mf)
  
  sheet2_down_bp=sheet2_down_plot[sheet2_down_plot$type=="biological pathway",]
  sheet2_down_cc=sheet2_down_plot[sheet2_down_plot$type=="cellular component",]
  sheet2_down_mf=sheet2_down_plot[sheet2_down_plot$type=="molecular function",]
  sheet2_down_bp= sheet2_down_bp[order(sheet2_down_bp$enrichment),]
  sheet2_down_cc= sheet2_down_cc[order(sheet2_down_cc$enrichment),]
  sheet2_down_mf= sheet2_down_mf[order(sheet2_down_mf$enrichment),]
  sheet2_down_plot=rbind(sheet2_down_bp,sheet2_down_cc,sheet2_down_mf)
  
  
  
  #pdf("mutation_type.pdf")
  sheet1_names = as.vector(sheet1_plot[,1])
  sheet1_nchr = nchar(sheet1_names)
  #the gene names 
  if (is.na(which(sheet1_nchr > 70)[1])){
    
  }else{
    sheet1_names[which(sheet1_nchr > 70)]  = paste0(substring( sheet1_names[which(sheet1_nchr > 70)],0,70),"...")
    sheet1_plot[,1] = sheet1_names
  }
  
  sheet3_names = as.vector(sheet3_plot[,1])
  sheet3_nchr = nchar(sheet3_names)
  #the gene names 
  if (is.na(which(sheet3_nchr > 70)[1])){
    
  }else{
    sheet3_names[which(sheet3_nchr > 70)]  = paste0(substring( sheet3_names[which(sheet3_nchr > 70)],0,70),"...")
    sheet3_plot[,1] = sheet3_names
  }
  
  
  sheet2up_names = as.vector(sheet2_up_plot[,1])
  sheet2up_nchr = nchar(sheet2up_names)
  #the gene names 
  if (is.na(which(sheet2up_nchr > 70)[1])){
    
  }else{
    sheet2up_names[which(sheet2up_nchr> 70)]  = paste0(substring( sheet2up_names[which(sheet2up_nchr > 70)],0,70),"...")
    sheet2_up_plot[,1] = sheet2up_names
    
  }
  
  sheet2down_names = as.vector(sheet2_down_plot[,1])
  sheet2down_nchr = nchar(sheet2down_names)
  #the gene names 
  if (is.na(which(sheet2down_nchr > 70)[1])){
    
  }else{
    sheet2down_names[which(sheet2down_nchr> 70)]  = paste0(substring( sheet2down_names[which(sheet2down_nchr > 70)],0,70),"...")
    sheet2_down_plot[,1] = sheet2down_names
  }
  
  
  p1 <- ggplot(sheet1_plot,aes(x=go_name,y=LgP,fill=type,group=type))+geom_bar(position="identity",stat="identity") +coord_flip()+scale_x_discrete(limits=sheet1_plot$go_name)
  p1+xlab("Gene ontology category") + ylab("(-LgP)") + labs(fill="biotype")
  ggsave(paste0(inputDir, "/UPgene_Sig_Go.pdf"),width = 9.18, height = 6.53)
  p1 <- ggplot(sheet1_plot,aes(x=go_name,y=LgP,fill=type,group=type))+geom_bar(position="identity",stat="identity") +coord_flip()+scale_x_discrete(limits=sheet1_plot$go_name)
  p1+xlab("Gene ontology category") + ylab("(-LgP)") + labs(fill="biotype")
  
  ggsave(paste0(inputDir, "/UPgene_Sig_Go.png"),width = 9.18, height = 6.53)
  
  p3 <- ggplot(sheet3_plot,aes(x=go_name,y=LgP,fill=type)) +geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=sheet3_plot$go_name)
  p3+labs(title = "Difgene Sig Go", y ="(-LgP)",x="Gene ontology category",labs(fill="biotype")) 
  ggsave(paste0(inputDir, "/Downgene_Sig_Go.pdf"),width = 9.18, height = 6.53)
  p3 <- ggplot(sheet3_plot,aes(x=go_name,y=LgP,fill=type)) +geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=sheet3_plot$go_name)
  p3+labs(title = "Difgene Sig Go", y ="(-LgP)",x="Gene ontology category",labs(fill="biotype")) 
  ggsave(paste0(inputDir, "/Downgene_Sig_Go.png"),width = 9.18, height = 6.53)
  
  
  p2_up <- ggplot(sheet2_up_plot,aes(x=go_name,y=enrichment,fill=type)) +geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=sheet2_up_plot$go_name)
  p2_up+labs(title = "Difgene Sig Go", y ="enrichment",x="Gene ontology category",labs(fill="biotype")) 
  ggsave(paste0(inputDir, "/Upgene_Sig_GoEnrichment.pdf"),width = 9.18, height = 6.53)
  p2_up <- ggplot(sheet2_up_plot,aes(x=go_name,y=enrichment,fill=type)) +geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=sheet2_up_plot$go_name)
  p2_up+labs(title = "Difgene Sig Go", y ="enrichment",x="Gene ontology category",labs(fill="biotype")) 
  ggsave(paste0(inputDir, "/Upgene_Sig_GoEnrichment.png"),width = 9.18, height = 6.53)
  
  
  p2_down <- ggplot(sheet2_down_plot,aes(x=go_name,y=enrichment,fill=type)) +geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=sheet2_down_plot$go_name)
  p2_down+labs(title = "Difgene Sig Go", y ="enrichment",x="Gene ontology category",labs(fill="biotype")) 
  ggsave(paste0(inputDir, "/Downgene_Sig_GoEnrichment.pdf"),width = 9.18, height = 6.53)
  p2_down <- ggplot(sheet2_down_plot,aes(x=go_name,y=enrichment,fill=type)) +geom_bar(stat="identity")+coord_flip()+scale_x_discrete(limits=sheet2_down_plot$go_name)
  p2_down+labs(title = "Difgene Sig Go", y ="enrichment",x="Gene ontology category",labs(fill="biotype")) 
  ggsave(paste0(inputDir, "/Downgene_Sig_GoEnrichment.png"),width = 9.18, height = 6.53)
  
  #-----------------------------------------
  path_upGene =  read.xlsx(paste0(inputDir, "/up_gene_goPathway.xlsx"),sheetIndex = 5)
  
  
  path_downGene =  read.xlsx(paste0(inputDir, "/down_gene_goPathway.xlsx"),sheetIndex = 5)
  
  
  sheet1 = path_upGene[which(path_upGene$pvalue<0.05),c("pathID","overlapGeneCount","pvalue","bgRatio")]
  sheet3 = path_downGene[which(path_downGene$pvalue<0.05),c("pathID","overlapGeneCount","pvalue","bgRatio")] 
  
  sheet1_genenum=as.data.frame(apply(as.data.frame(sheet1$bgRatio),2,strsplit,"/"))[1,]
  sheet1_genenum=t(sheet1_genenum)
  sheet1$bgRatio=sheet1_genenum
  
  sheet3_genenum=as.data.frame(apply(as.data.frame(sheet3$bgRatio),2,strsplit,"/"))[1,]
  sheet3_genenum=t(sheet3_genenum)
  sheet3$bgRatio=sheet3_genenum
  
  sheet1$pvalue = -log10(sheet1$pvalue) 
  colnames(sheet1) = c("path_name","path_diffgene_count","LgP","bg_count")
  sheet3$pvalue = -log10(sheet3$pvalue) 
  colnames(sheet3) = c("path_name","path_diffgene_count","LgP","bg_count")
  
  
  if (nrow(sheet1)>=30 & nrow(sheet3)>=30){
    
    sheet1_plot = sheet1[1:30,]
    sheet3_plot = sheet3[1:30,]
    
  }else if (nrow(sheet1)>=30 & nrow(sheet3) < 30){
    sheet1_plot = sheet1[1:30,]
    sheet3_plot = sheet3
  }else if (nrow(sheet1)<30 & nrow(sheet3) >= 30){
    
    sheet1_plot = sheet1
    sheet3_plot = sheet3[1:30,]
  }else{
    sheet1_plot = sheet1
    sheet3_plot = sheet3
    
  }
  
  
  
  
  #
  sheet1_plot = sheet1_plot[order(sheet1_plot$LgP),]
  sheet3_plot = sheet3_plot[order(sheet3_plot$LgP),]
  
  
  ggplot(sheet1_plot,aes(LgP,path_name))+ geom_point()  + geom_point(aes(size=path_diffgene_count,colour=LgP))+scale_color_gradient(low="green",high = "red")+labs(color=expression(LgP),size="Diffgene_count",x="(-LgP)",y="Pathway",title="Pathway enrichment")
  
  ggsave(paste0(inputDir, "/UPgene_Sig_Pathway.pdf"),width = 8, height = 8)
  
  ggplot(sheet1_plot,aes(LgP,path_name))+ geom_point()  + geom_point(aes(size=path_diffgene_count,color=LgP))+scale_color_gradient(low="green",high = "red")+labs(color=expression(LgP),size="Diffgene_count",x="(-LgP)",y="Pathway",title="Pathway enrichment")
  
  ggsave(paste0(inputDir, "/UPgene_Sig_Pathway.png"),width = 8, height = 8)
  
  ggplot(sheet3_plot,aes(LgP,path_name))+ geom_point()  + geom_point(aes(size=path_diffgene_count,color=LgP))+scale_color_gradient(low="green",high = "red")+labs(color=expression(LgP),size="Diffgene_count",x="(-LgP)",y="Pathway",title="Pathway enrichment")
  
  ggsave(paste0(inputDir, "/Downgene_Sig_Pathway.pdf"),width = 8, height = 8)
  
  ggplot(sheet3_plot,aes(LgP,path_name))+ geom_point()  + geom_point(aes(size=path_diffgene_count,color=LgP))+scale_color_gradient(low="green",high = "red")+labs(color=expression(LgP),size="Diffgene_count",x="(-LgP)",y="Pathway",title="Pathway enrichment")
  
  ggsave(paste0(inputDir, "/Downgene_Sig_Pathway.png"),width = 8, height = 8)
  
}

    