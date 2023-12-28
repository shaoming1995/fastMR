#' @title 基于暴露本地,结局是在线数据的标准多变量孟德尔随机化
#' @param keyssh 密钥
#' @param exposure_dat_iv 输入合并的暴露工具变量文件
#' @param outgwas 数据结局的GWAS summary数据
#' @param ivw_MVMR 是否进行IVW方法估计
#' @param lasso_MVMR 是否进行lasso方法估计
#' @param egger_MVMR 是否进行egger方法估计
#' @param median_MVMR 是否进行median方法估计
#' @export
stand_MVMR_local_local<-function(keyssh,exposure_dat_iv,outgwas,
                               ivw_MVMR=T,lasso_MVMR=F,egger_MVMR=F,median_MVMR=F){
  if (Sys.info()["nodename"] == keyssh) {
    exposure_dat <- exposure_dat_iv
    #获取结局数据
    outcome_dat <- merge(exposure_dat,outgwas,by="SNP",all=F)
    #多变量MR分析
    outcome_dat<-outcome_dat[!duplicated(outcome_dat$SNP),]
    mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)
    #转换数据格式讲S3变成S4数据格式
    library(MendelianRandomization)
    MRMVInput <- mr_mvinput(bx = mvdat$exposure_beta,
                            bxse = mvdat$exposure_se,
                            by = mvdat$outcome_beta,
                            byse = mvdat$outcome_se,
                            correlation =matrix())
    if(ivw_MVMR==T){
      res <- mv_multiple(mvdat)
      #查看结果
      res<-generate_odds_ratios(res$result)
      return(res)
    }
    if(lasso_MVMR==T){
      #LASSO回归
      mv_lasso<-mr_mvlasso(MRMVInput)
      return(mv_lasso)}
    if(egger_MVMR==T){
      mvergger<-mr_mvegger(MRMVInput)#egger法
      return(mvergger)}
    if(median_MVMR==T){
      #加权法
      mv_median<-mr_mvmedian(MRMVInput)
      #查看加权法结果
      return(mv_median)}
  }else{
    cat("请联系管理员获取账户学号和密码或微信联系SFM19950928或DKYXS666")}
}
