#' @title 读取731免疫细胞获取SNP的参考面板
#' @param inputfile 输入参考文件的位置
#' @export
read_ref_data<-function(inputfile){
library(tidyr)
path<-paste0(inputfile,"ref.txt.gz")
data<-data.table::fread(path)%>%data.frame()
return(data)
}
