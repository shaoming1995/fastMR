#' @title 代谢数据预处理
#' @param inputfile 代谢数据文件位置
#' @param savefile 预处理后数据保存位置
#' @param exp_or_out 处理成暴露还是结局，默认是T表示处理成暴露
#' @param P_exp 设定处理暴露用于潜在工具变量的P值
#' @param P_out 设定处理结局用于潜在工具变量的P值
#' @export
metb_pre<-function(inputfile,savefile,exp_or_out=T,P_exp=1e-05,P_out=5e-08){
  library(tidyr)
  file<-dir(inputfile)%>%data.frame()
  for(i in 1:nrow(file)){
    #dir.create(savefile)
    path<-paste0(inputfile,"/",file[i,1])
    data<-data.table::fread(path)%>%data.frame()

    head(data)

    colnames(data)[c(3:9)]<-c("effect_allele.exposure","other_allele.exposure","eaf.exposure",
                              "beta.exposure","se.exposure" ,"pval.exposure","SNP")
    data<-data[,c("effect_allele.exposure","other_allele.exposure",
                  "beta.exposure","se.exposure" ,"eaf.exposure","pval.exposure","SNP")]

    data$samplesize.exposure<-8264
    data$id.exposure<-file[i,]
    data$exposure<-data$id.exposure

    if(exp_or_out==T){
      data$pval.exposure<-as.numeric(data$pval.exposure)
      data<-subset(data,pval.exposure<P_exp)
      data[data==""]<-NA
      data<-na.omit(data)
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
        data2[data2==""]<-NA
        data2<-na.omit(data2)
         write.table(data2, file = gzfile(pathe), row.names = F)
        pv<-round(i/nrow(file),4)*100
        cat("已完成结果数据转化",pv,"%")
      }


  }
}


