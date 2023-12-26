#' @title 来源烟草数据库： https://conservancy.umn.edu/handle/11299/201564
#' @param file 下载的文件名称
#' @param name 输入表型名称
#' @param exp_or_out 是否作为暴露数据集，默认是T
#' @export
yancao_read_data<-function(file,name,exp_or_out=T){

  A<-fread(file)%>%data.frame()
  if(exp_or_out==T){
  colnames(A)[c(3,5,4,6,8,9,10,12)]<-c("SNP",
                                         "effect_allele.exposure",
                                         "other_allele.exposure",
                                         "eaf.exposure",
                                         "pval.exposure",
                                         "beta.exposure",
                                         "se.exposure","samplesize.exposure")
  A$id.exposure<-name
  A$exposure<-name
  return(A)}else{
    colnames(A)[c(3,5,4,6,8,9,10,12)]<-c("SNP",
                                         "effect_allele.outcome",
                                         "other_allele.outcome",
                                         "eaf.outcome",
                                         "pval.outcome",
                                         "beta.outcome",
                                         "se.outcome","samplesize.outcome")
    A$id.outcome<-name
    A$outcome<-name

  }

}
