#' @title 多变量暴露工具变量合并
#' @description 适用于1个来自IEU的暴露，4个来自非IEU的暴露工具变量合并
#' @param keyssh 密钥
#' @param data1IV 来自非IEU暴露1工具变量文件
#' @param data1GWAS 来自非IEU暴露1的GWAS summary数据文件
#' @param data2IV 来自非IEU暴露2工具变量文件
#' @param data2GWAS 来自非IEU暴露2的GWAS summary数据文件
#' @param data3IV 来自非IEU暴露3工具变量文件
#' @param data3GWAS 来自非IEU暴露3的GWAS summary数据文件
#' @param data4IV 来自非IEU暴露4工具变量文件
#' @param data4GWAS 来自非IEU暴露4的GWAS summary数据文件
#' @param data5IV 来自非IEU暴露4工具变量文件
#' @param data5GWAS 来自非IEU暴露4的GWAS summary数据文件
#' @export
IEU_1_to_local_4<-function(keyssh,data1IV,GWASID,data2IV,data2GWAS,data3IV,data3GWAS,data4IV,data4GWAS,data5IV,data5GWAS){
  if (!require(tidyfst)) install.packages("tidyfst")
  library(tidyfst)
  if (Sys.info()["nodename"] == keyssh) {
    #假设data1IV的暴露来自IEU data2IV和data2IV的暴露来自非IEU
    exp_name<-c("SNP","effect_allele.exposure","other_allele.exposure", "eaf.exposure", "beta.exposure","se.exposure", "pval.exposure","id.exposure","exposure")
    out_name<-c("SNP","effect_allele.outcome","other_allele.outcome", "eaf.outcome", "beta.outcome","se.outcome","pval.outcome","id.outcome","outcome")
    #去EXP2和EXP3和 EXP4 GWAS summary中去淘暴露1的工具变量
    data_h_SNP_steiger_mv1_2<-merge(data1IV[,c("SNP","outcome")],data2GWAS,by="SNP",all=F)
    data_h_SNP_steiger_mv1_3<-merge(data1IV[,c("SNP","outcome")],data3GWAS,by="SNP",all=F)
    data_h_SNP_steiger_mv1_4<-merge(data1IV[,c("SNP","outcome")],data4GWAS,by="SNP",all=F)
    data_h_SNP_steiger_mv1_5<-merge(data1IV[,c("SNP","outcome")],data5GWAS,by="SNP",all=F)
    #去EXP1和EXP3和EXP4GWAS summary中去淘暴露2的工具变量
    data_h_SNP_steiger_mv2_1<-extract_outcome_data(snps=data2IV$SNP,outcomes=GWASID,proxies=T,maf_threshold = 0.01,access_token = NULL)
    data_h_SNP_steiger_mv2_1<-data_h_SNP_steiger_mv2_1[,out_name]
    colnames(data_h_SNP_steiger_mv2_1)<-exp_name
    data_h_SNP_steiger_mv2_3<-merge(data2IV[,c("SNP","outcome")],data3GWAS,by="SNP",all=F)
    data_h_SNP_steiger_mv2_4<-merge(data2IV[,c("SNP","outcome")],data4GWAS,by="SNP",all=F)
    data_h_SNP_steiger_mv2_5<-merge(data2IV[,c("SNP","outcome")],data5GWAS,by="SNP",all=F)
    #去EXP1和EXP2 和EXP4 GWAS summary中去淘暴露3的工具变量
    data_h_SNP_steiger_mv3_1<-extract_outcome_data(snps=data3IV$SNP,outcomes=GWASID,proxies=T,maf_threshold = 0.01,access_token = NULL)
    data_h_SNP_steiger_mv3_1<-data_h_SNP_steiger_mv3_1[,out_name]
    colnames(data_h_SNP_steiger_mv3_1)<-exp_name
    data_h_SNP_steiger_mv3_2<-merge(data3IV[,c("SNP","outcome")],data2GWAS,by="SNP",all=F)
    data_h_SNP_steiger_mv3_4<-merge(data3IV[,c("SNP","outcome")],data4GWAS,by="SNP",all=F)
    data_h_SNP_steiger_mv3_5<-merge(data3IV[,c("SNP","outcome")],data5GWAS,by="SNP",all=F)
    #去EXP1和EXP2 和EXP3 GWAS summary中去淘暴露3的工具变量
    data_h_SNP_steiger_mv4_1<-extract_outcome_data(snps=data4IV$SNP,outcomes=GWASID,proxies=T,maf_threshold = 0.01,access_token = NULL)
    data_h_SNP_steiger_mv4_1<-data_h_SNP_steiger_mv4_1[,out_name]
    colnames(data_h_SNP_steiger_mv4_1)<-exp_name
    data_h_SNP_steiger_mv4_2<-merge(data4IV[,c("SNP","outcome")],data2GWAS,by="SNP",all=F)
    data_h_SNP_steiger_mv4_3<-merge(data4IV[,c("SNP","outcome")],data3GWAS,by="SNP",all=F)
    data_h_SNP_steiger_mv4_5<-merge(data4IV[,c("SNP","outcome")],data5GWAS,by="SNP",all=F)

    data_h_SNP_steiger_mv5_1<-extract_outcome_data(snps=data5IV$SNP,outcomes=GWASID,proxies=T,maf_threshold = 0.01,access_token = NULL)
    data_h_SNP_steiger_mv5_1<-data_h_SNP_steiger_mv5_1[,out_name]
    colnames(data_h_SNP_steiger_mv5_1)<-exp_name
    data_h_SNP_steiger_mv5_2<-merge(data5IV[,c("SNP","outcome")],data2GWAS,by="SNP",all=F)
    data_h_SNP_steiger_mv5_3<-merge(data5IV[,c("SNP","outcome")],data3GWAS,by="SNP",all=F)
    data_h_SNP_steiger_mv5_4<-merge(data5IV[,c("SNP","outcome")],data4GWAS,by="SNP",all=F)
    #开始合并暴露工具变量
    MVIV1_1<-data1IV[,exp_name]#自己本身的IVs
    MVIV1_2<-data_h_SNP_steiger_mv2_1[,exp_name]#EXP2在自己GWAS中的IVs
    MVIV1_3<-data_h_SNP_steiger_mv3_1[,exp_name]#EXP3在自己GWAS中的IVs
    MVIV1_4<-data_h_SNP_steiger_mv4_1[,exp_name]#EXP4在自己GWAS中的IVs
    MVIV1_5<-data_h_SNP_steiger_mv5_1[,exp_name]#EXP4在自己GWAS中的IVs

    MVIV2_2<-data2IV[,exp_name]#自己本身的IVs
    MVIV2_1<-data_h_SNP_steiger_mv1_2[,exp_name]#EXP1在自己GWAS中的IVs
    MVIV2_3<-data_h_SNP_steiger_mv3_2[,exp_name]#EXP3在自己GWAS中的IVs
    MVIV2_4<-data_h_SNP_steiger_mv4_2[,exp_name]#EXP4在自己GWAS中的IVs
    MVIV2_5<-data_h_SNP_steiger_mv5_2[,exp_name]#EXP4在自己GWAS中的IVs

    MVIV3_3<-data3IV[,exp_name]#自己本身的IVs
    MVIV3_1<-data_h_SNP_steiger_mv1_3[,exp_name]#EXP1在自己GWAS中的IVs
    MVIV3_2<-data_h_SNP_steiger_mv2_3[,exp_name]#EXP2在自己GWAS中的IVs
    MVIV3_4<-data_h_SNP_steiger_mv4_3[,exp_name]#EXP2在自己GWAS中的IVs
    MVIV3_5<-data_h_SNP_steiger_mv5_3[,exp_name]#EXP2在自己GWAS中的IVs

    MVIV4_4<-data4IV[,exp_name]#自己本身的IVs
    MVIV4_1<-data_h_SNP_steiger_mv1_4[,exp_name]#EXP1在自己GWAS中的IVs
    MVIV4_2<-data_h_SNP_steiger_mv2_4[,exp_name]#EXP2在自己GWAS中的IVs
    MVIV4_3<-data_h_SNP_steiger_mv3_4[,exp_name]#EXP2在自己GWAS中的IVs
    MVIV4_5<-data_h_SNP_steiger_mv5_4[,exp_name]#EXP2在自己GWAS中的IVs

    MVIV5_5<-data5IV[,exp_name]#自己本身的IVs
    MVIV5_1<-data_h_SNP_steiger_mv1_5[,exp_name]#EXP1在自己GWAS中的IVs
    MVIV5_2<-data_h_SNP_steiger_mv2_5[,exp_name]#EXP2在自己GWAS中的IVs
    MVIV5_3<-data_h_SNP_steiger_mv3_5[,exp_name]#EXP3在自己GWAS中的IVs
    MVIV5_4<-data_h_SNP_steiger_mv4_5[,exp_name]#EXP4在自己GWAS中的IVs
    #合并工具变量
    exposure_dat_temp<-rbind(MVIV1_1,MVIV1_2,MVIV1_3,MVIV1_4,MVIV1_5,
                             MVIV2_2,MVIV2_1,MVIV2_3,MVIV2_4,MVIV2_5,
                             MVIV3_3,MVIV3_1,MVIV3_2,MVIV3_4,MVIV3_5,
                             MVIV4_4,MVIV4_1,MVIV4_2,MVIV4_3,MVIV4_5,
                             MVIV5_5,MVIV5_1,MVIV5_2,MVIV5_3,MVIV5_4)

    #写出工具变量剔除唯一的SNP
    library(tidyr)
    exposure_dat_temp<-exposure_dat_temp %>% add_count_dt(exposure_dat_temp,.name = "N")
    exposure_dat_temp<-subset(exposure_dat_temp,N==5)
    return(exposure_dat_temp)
    #write.csv(exposure_dat_temp,"exposure_dat.csv",quote = F,row.names = F)
    #cat("已完成多变量暴露工具变量合并，请前往文件夹下exposure_dat.csv文件中进行人工筛选")
    message("已完成多变量暴露工具变量合并")
  }else{
      cat("请联系管理员获取账户学号和密码或微信联系SFM19950928或DKYXS666")}
}
