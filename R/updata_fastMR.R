#' @title 更新fastMR包
#' @export
updata_fastMR<-function(){
  x<-try(detach("package:fastMR", unload = TRUE))
  if(class(x)!="try-error"){
    remove.packages("fastMR")
    devtools::install_github("shaoming1995/fastMR")
    message("fastMR已经更新到最新版本")}else{
      remove.packages("fastMR")
      devtools::install_github("shaoming1995/fastMR")
      message("fastMR已经更新到最新版本")}
}
