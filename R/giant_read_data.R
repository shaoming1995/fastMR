#' @title 来源GIANT数据库： https://portals.broadinstitute.org/collaboration/giant/index.php
#' @param file 下载的文件名称
#' @param name 输入表型名称
#' @param exp_or_out 是否作为暴露数据集，默认是T
#' @export
giant_read_data<-function(file,name,exp_or_out=T){

  A<-fread(file)%>%data.frame()
  if(exp_or_out==T){
    colnames(A)[c(3,5,4,6,7,8,9,10)]<-c("SNP",
                                         "effect_allele.exposure",
                                         "other_allele.exposure",
                                         "eaf.exposure",
                                         "beta.exposure",
                                         "se.exposure",
                                         "pval.exposure",
                                         "samplesize.exposure")
    A$id.exposure<-name
    A$exposure<-name
    return(A)}else{
      colnames(A)[c(3,5,4,6,7,8,9,10)]<-c("SNP",
                                          "effect_allele.outcome",
                                          "other_allele.outcome",
                                          "eaf.outcome",
                                          "beta.outcome",
                                          "se.outcome",
                                          "pval.outcome",
                                          "samplesize.outcome")
      A$id.outcome<-name
      A$outcome<-name
      return(A)
    }

}
