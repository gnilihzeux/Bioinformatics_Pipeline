HeatMap<-function(data,
                  annotation_col=NA,
                  path=getwd(),
                  adjust_data=FALSE,
                  data_log=FALSE,
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  color = colorRampPalette(c("blue", "white", "red"))(50),
                  height=8,
                  width=6,
                  Figurename="Diff_heatmap.pdf",
                  fontsize = 10,...
  
){
  
  
  library(pheatmap)
  data[is.na(data)] <- 0
  
  if(data_log){
    data = log2(data + 1 )
  }
  
 
  
  
  if (adjust_data){
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
    
    data <- comprehensive_normalization(x=data,control_method="median",standardize=FALSE);
    contrast_threshold <- signif(quantile(data,probs=0.9,names=FALSE),2)
    imageVals <- data;
    imageVals[data > contrast_threshold] <- contrast_threshold;
    imageVals[data < -1*contrast_threshold] <- -1*contrast_threshold;
    data <- imageVals
  }
  
  
  
  
  if ( nrow (data) > 35) {
    show_genenames = FALSE
  }else{
    show_genenames = TRUE
  }
 
  if ( ncol (data) > 35) {
    show_colnames = FALSE
  }else{
    show_colnames = TRUE
  }
  
  
  
  pheatmap(data,
           scale                    = "row", 
           clustering_distance_rows = clustering_distance_rows,
           clustering_distance_cols = clustering_distance_cols,
           color                    = color,
           cluster_rows             = TRUE,
           cluster_cols             = TRUE, 
           border_color             = FALSE,
           annotation_col           = annotation_col,
           show_rownames            = show_genenames, 
           show_colnames            = show_colnames, 
           number_color             = "blue",
           height                   = height,
           width                    = width,
           filename                 = paste0(path,"/",Figurename),
           fontsize                 = fontsize
           )
  
}