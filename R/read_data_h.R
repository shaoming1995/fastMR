#' @title 读取工具变量用于多变量孟德尔随机化分析
#' @param outfile 输入单变量文件夹存放工具变量的文件名
#' @export
read_data_h<-function(outfile){
  path0<-paste0(getwd(),"/",outfile)
  path00<-paste0(path0,"/IV.csv")
  Atemp<- read.csv(path00,header = T)
  Atemp<-Atemp[,c("SNP","effect_allele.exposure","other_allele.exposure", "eaf.exposure",
                  "beta.exposure","se.exposure", "pval.exposure","id.exposure","exposure","mr_keep")]

  Atemp<-subset(Atemp,mr_keep=="TRUE")
  return(Atemp)
}


