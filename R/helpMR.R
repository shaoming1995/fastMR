#' @title 获取安装标准MR分析的密钥
#' @param Yourname 输入您的用户名
#' @export
helpMR<-function(Yourname=NULL){
  if(class(Yourname)!="NULL"){
  if(!require("keyring",quietly=T))
    result1 <-try(install.packages("keyring"))
  if(!require("blastula",quietly=T))
    result2 <-try(install.packages("blastula"))
  if(!require("rmarkdown",quietly=T))
    result3 <- try(install.packages("rmarkdown"))
  result4 <-try(library(keyring))
  result5 <-try(library(blastula))
  result6 <-try(library(rmarkdown))
  #if(!require("sofmClump",quietly=T))
    #try(devtools::install_github("shaoming1995/sofmClump"))
  #if(!require("sofmstandMR",quietly=T))
  #devtools::install_github("shaoming1995/sofmstandMR")
  id_temp<-paste0(round(runif(1,min = 1,max = 100000),0))
  result7 <-try(create_smtp_creds_key(
    id = id_temp,
    user = "sofm_teacher@163.com",
    host = "smtp.163.com",
    port = 25,
    use_ssl = T))
  subject = "please usering sofmClump"
  from <- "sofm_teacher@163.com"
  to <- "shaoming0928@stu.njmu.edu.cn"
  email <-compose_email(body = md(paste0(Yourname,"的密钥是",Sys.info()["nodename"])))
  smtp_send(
    from = from,
    to = to,
    subject =subject,
    email = email,
    credentials = creds_key(id = id_temp)
  )}else{message("必须输入DK开头的用户名才可以激活安装
抖音ID：793742981（医小研）
微信:SFM19950928
感谢您的使用,我们会不断为您提供优质的R包服务
")}
}








