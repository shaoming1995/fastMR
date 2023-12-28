#' @title 多变量暴露工具变量合并
#' @description 适用于0个来自IEU的暴露，2个来自非IEU的暴露工具变量合并
#' @param data1IV 来自非IEU暴露1工具变量文件
#' @param data1GWAS 来自非IEU暴露1的GWAS summary数据文件
#' @param data2IV 来自非IEU暴露2工具变量文件
#' @param keyssh 密钥
#' @param data2GWAS 来自非IEU暴露2的GWAS summary数据文件
#' @export
IEU_0_to_local_2<-function(keyssh,data1IV,data1GWAS,data2IV,data2GWAS){
  if (!require(tidyfst)) install.packages("tidyfst")
  library(tidyfst)
  if (Sys.info()["nodename"] == keyssh) {
  #假设data1IV的暴露来自IEU data2IV和data2IV的暴露来自非IEU
  exp_name<-c("SNP","effect_allele.exposure","other_allele.exposure", "eaf.exposure", "beta.exposure","se.exposure", "pval.exposure","id.exposure","exposure")
  out_name<-c("SNP","effect_allele.outcome","other_allele.outcome", "eaf.outcome", "beta.outcome","se.outcome","pval.outcome","id.outcome","outcome")
  #去EXP2中淘EXP1 IV
  data_h_SNP_steiger_mv1_2<-merge(data1IV[,c("SNP","outcome")],data2GWAS,by="SNP",all=F)
  #去EXP1中淘EXP2 IV
  data_h_SNP_steiger_mv2_1<-merge(data2IV[,c("SNP","outcome")],data1GWAS,by="SNP",all=F)
  #开始合并暴露工具变量
  MVIV1_1<-data1IV[,exp_name]#自己本身的IVs
  MVIV1_2<-data_h_SNP_steiger_mv2_1[,exp_name]#EXP2在自己GWAS中的IVs


  MVIV2_2<-data2IV[,exp_name]#自己本身的IVs
  MVIV2_1<-data_h_SNP_steiger_mv1_2[,exp_name]#EXP1在自己GWAS中的IVs
  #合并工具变量
  exposure_dat_temp<-rbind(MVIV1_1,MVIV1_2,
                           MVIV2_2,MVIV2_1)
  #写出工具变量剔除唯一的SNP
  # 对添加列的命名
  library(tidyr)
  exposure_dat_temp<-exposure_dat_temp %>% add_count_dt(exposure_dat_temp,.name = "N")
  exposure_dat_temp<-subset(exposure_dat_temp,N==2)
  return(exposure_dat_temp)
  #write.csv(exposure_dat_temp,"exposure_dat.csv",quote = F,row.names = F)
  #cat("已完成多变量暴露工具变量合并，请前往文件夹下exposure_dat.csv文件中进行人工筛选")
  message("已完成多变量暴露工具变量合并")
  }else{
    cat("请联系管理员服获取账户学号和密码或微信联系SFM19950928或DKYXS666")}
}




