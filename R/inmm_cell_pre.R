#' @title 免疫细胞数据预处理
#' @param inputfile 免疫细胞数据文件位置
#' @param savefile 预处理后数据保存位置
#' @param exp_or_out 处理成暴露还是结局，默认是T表示处理成暴露
#' @param P_exp 设定处理暴露用于潜在工具变量的P值
#' @param P_out 设定处理结局用于潜在工具变量的P值
#' @export
inmm_cell_pre<-function(inputfile,savefile,exp_or_out=T,P_exp=1e-05,P_out=5e-08){
library(tidyr)
file<-dir(inputfile)%>%data.frame()

for(i in 1:nrow(file)){
  dir.create(savefile)
  path<-paste0(inputfile,"/",file[i,1])
  data<-data.table::fread(path)%>%data.frame()
  data$mer<-paste0(data$chromosome,":",data$base_pair_location)
  colnames(data)[c(3,4,5,7,8,9,10)]<-c("effect_allele.exposure","other_allele.exposure", "samplesize.exposure",
                                       "eaf.exposure","beta.exposure","se.exposure","pval.exposure")
  data$id.exposure<-file[i,]
  data$exposure<-data$id.exposure
  data<-data[,c("effect_allele.exposure","other_allele.exposure","samplesize.exposure",
                "beta.exposure","se.exposure" ,"eaf.exposure","pval.exposure","mer","id.exposure","exposure")]
  if(exp_or_out==T){
    data$pval.exposure<-as.numeric(data$pval.exposure)
    data<-subset(data,pval.exposure<P_exp)
    data1<-merge(ref_data,data,by="mer",all = F)
    pathe<-paste0(savefile,"/",file[i,1],".csv")
    write.csv(data1,pathe,row.names = F,quote = F)
    pv<-round(i/nrow(file),4)*100
    cat("已完成暴露数据转化",pv,"%")}else{
      data2<-data[,c("effect_allele.exposure","other_allele.exposure","samplesize.exposure",
                     "beta.exposure","se.exposure" ,"eaf.exposure","pval.exposure","mer","id.exposure","exposure")]

      colnames(data2)<-c("effect_allele.outcome","other_allele.outcome","samplesize.outcome",
                         "beta.outcome","se.outcome" ,"eaf.outcome","pval.outcome","mer","id.outcome","outcome")
      data2<-merge(ref_data,data2,by="mer",all = F)
      pathe<-paste0(savefile,"/",file[i,1])
      data2<-subset(data2,pval.outcome>P_out)
     write.table(data2, file = gzfile(pathe), row.names = F)
      pv<-round(i/nrow(file),4)*100
      cat("已完成结果数据转化",pv,"%")
    }
}
}
