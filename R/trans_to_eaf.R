#' @title 适用根据样本量与beta和se转换eaf
#' @param data 输入需要转换的数据
#' @param n_col 输入样本量列名
#' @param beta_col 输入效应系数列名
#' @param se_col 输入标准误列名
#' @export
trans_to_eaf<-function(data, n_col, beta_col, se_col){
   data1 <- data[, c(beta_col, se_col, n_col)]
  data1$d2 <- ((1/data1[, 2]) * (1/data1[, 2]))/2
  data1$z <- data1[, 1]/data1[, 2]
  data1$z2n <- data1$z * data1$z + data1[, 3]
  data1$zf1 <- data1$d2/data1$z2n
  data1$eaf <- (1 + sqrt(1 -4 * data1$zf1))/2
  return(data1$eaf)
}


