#' @title 来源芬兰数据库： https://www.finngen.fi/fi/hyodynna_tuloksia
#' @param file 下载的文件名称
#' @param name 输入表型名称
#' @param exp_or_out 是否作为暴露数据集，默认是T
#' @export
finn_read_data<-function(file,name,exp_or_out=T,N){

  A<-fread(file)%>%data.frame()
  if(exp_or_out==T){
    colnames(A)[c(4,3,5,7,9,10,11)]<-c("effect_allele.exposure","other_allele.exposure", "SNP","pval.exposure","beta.exposure","se.exposure",
                                       "eaf.exposure")
    A$id.exposure<-name
    A$exposure<-name
    A$samplesize.exposure<-N
    return(A)}else{
      colnames(A)[c(4,3,5,7,9,10,11)]<-c("effect_allele.outcome","other_allele.outcome", "SNP","pval.outcome","beta.outcome","se.outcome",
                                         "eaf.outcome")
      A$id.outcome<-name
      A$outcome<-name
      A$samplesize.outcome<-N
      return(A)
    }
  }
