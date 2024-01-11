#' @title 炎症因子数据预处理
#' @param inputfile 炎症因子数据文件位置
#' @param savefile 预处理后数据保存位置
#' @param exp_or_out 处理成暴露还是结局，默认是T表示处理成暴露
#' @param P_exp 设定处理暴露用于潜在工具变量的P值
#' @param P_out 设定处理结局用于潜在工具变量的P值
#' @export
infla_factor_pre<-function(inputfile,savefile,exp_or_out=T,P_exp=1e-05,P_out=5e-08){
  library(tidyr)
  file<-dir(inputfile)%>%data.frame()
  for(i in 1:nrow(file)){
    dir.create(savefile)
    path<-paste0(inputfile,"/",file[i,1])
    data<-data.table::fread(path)%>%data.frame()
    colnames(data)[c(3,4,11,5,6,7,8,10)]<-c("effect_allele.exposure","other_allele.exposure","samplesize.exposure",
                                            "beta.exposure","se.exposure" ,"eaf.exposure","pval.exposure","SNP")
    data<-data[,c("effect_allele.exposure","other_allele.exposure","samplesize.exposure",
                  "beta.exposure","se.exposure" ,"eaf.exposure","pval.exposure","SNP")]
    data$id.exposure<-file[i,]
    data$exposure<-data$id.exposure
    if(exp_or_out==T){
      data<-subset(data,pval.exposure<P_exp)
      pathe<-paste0(savefile,"/",file[i,1],".csv")
      write.csv(data,pathe,row.names = F,quote = F)
      pv<-round(i/nrow(file),4)*100
      cat("已完成暴露数据转化",pv,"%")}else{
        data2<-data[,c("effect_allele.exposure","other_allele.exposure","samplesize.exposure",
                       "beta.exposure","se.exposure" ,"eaf.exposure","pval.exposure","SNP","id.exposure","exposure")]
        colnames(data2)<-c("effect_allele.outcome","other_allele.outcome","samplesize.outcome",
                           "beta.outcome","se.outcome" ,"eaf.outcome","pval.outcome","SNP","id.outcome","outcome")
        pathe<-paste0(savefile,"/",file[i,1])
        data2<-subset(data2,pval.outcome>P_out)
        write.table(data2,file = pathe,row.names = F)
        pv<-round(i/nrow(file),4)*100
        cat("已完成结果数据转化",pv,"%")
      }


  }
}



