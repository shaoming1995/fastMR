#' @title 适用根据Z与eaf转换beta和se
#' @param data 输入需要转换的数据
#' @param n_col 输入样本量列名
#' @param eaf_col 输入eaf列名
#' @param Z_col 输入Z列名
#' @export
trans_from_Z<-function(data,n_col,eaf_col,Z_col,exp_or_out=T){
  if(exp_or_out==T){
  data1<-data[,c(n_col,eaf_col,Z_col)]
  data1$beta.exposure<-2*(data1[,eaf_col]*(1-data1[,eaf_col]))
  data1$beta.exposure<-data1$beta.exposure*(data1[,Z_col]*data1[,Z_col]+data1[,n_col])
  data1$beta.exposure<-sqrt(data1$beta.exposure)
  data1$beta.exposure<- data1[,Z_col]/data1$beta.exposure
  data1$se.exposure<-data1$beta.exposure/data1[,Z_col]
  data2<-cbind(data1[,c("beta.exposure","se.exposure")],data)
  return(data2)}else{
    data1<-data[,c(n_col,eaf_col,Z_col)]
    data1$beta.outcome<-2*(data1[,eaf_col]*(1-data1[,eaf_col]))
    data1$beta.outcome<-data1$beta.outcome*(data1[,Z_col]*data1[,Z_col]+data1[,n_col])
    data1$beta.outcome<-sqrt(data1$beta.outcome)
    data1$beta.outcome<- data1[,Z_col]/data1$beta.outcome
    data1$se.outcome<-data1$beta.outcome/data1[,Z_col]
    data2<-cbind(data1[,c("beta.outcome","se.outcome")],data)
    return(data2)
  }
}
