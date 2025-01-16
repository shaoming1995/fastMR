#' @title 用于两表型的GWAS荟萃分析
#' @param data1 输入标准MR格式的暴露列名文件，额外包含必要的两列，chr.exposure,pos.exposure
#' @param data2 输入第二个表型的GWAS文件，需要包含必要的三列，SNP,BETA,SE
#' @param snp_col 输入SNP列名
#' @param beta_col 表型2的效应值列名
#' @param se_col 表型2的效应值标准误列名
#' @param a1_col 表型2的效应等位基因列名
#' @param a2_col 表型2的参考等位基因列名
#' @param filename1 输入临时储存文件夹1
#' @param filename2 输入临时储存文件夹2
#' @param Nsample 输入合计的样本量
#' @param triat 输入表型的名称
#' @export

GWAS_meta<-function(data1,data2,snp_col="SNP",beta_col="BETA",se_col="SE",a1_col="effect_allele",a2_col="other_allele",filename1="芬兰肺癌GWAS",filename2="catlog肺癌GWAS",Nsample=1000000,triat="Lung cancer"){
  RegistID_dat <- RegistID_dat
  RegistID_u <- subset(RegistID_dat, IK == keyssh)
  tempid <- paste0(keyssh, "_", Sys.info()["nodename"], "_",
                   RegistID_u$RegistID)
  if (RegistID_u$FINN %in% tempid) {
  if (!requireNamespace("progress", quietly = TRUE))
    install.packages("progress")
  library(tidyr)
  library(data.table)
  library(progress)
  #合并交集
  message("正在获取有效SNP...")
  d1_d2<-merge(data1,
               data2[,c(snp_col,beta_col,se_col,a1_col,a2_col)],
               by.x = "SNP",
               by.y = snp_col,
               all = F)
  #查看有效SNP数目
  message(paste0("获取完成，一共获取到",dim(d1_d2)[1]),"个SNP进行meta分析")
  message("开始进行切割数据...")
  #写出备用文件
  d11<-d1_d2[,c("SNP","beta.exposure","se.exposure","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","pval.exposure")]
  colnames(d11)<-c("SNP","OR","SE","CHR","BP","A1","A2","P")
  d11$OR<-exp(d11$OR)
  N<-22
  pb<-progress_bar$new(total=N)
  for(i in 1:N){
    d_chr<-subset(d11,CHR==i)
    dir.create(filename1)
    write.table(d_chr,paste0(getwd(),"/",filename1,"/芬兰肺癌GWAS",i,".txt"),quote=F,row.names=F)
    pb$tick()
  }

  message("数据集1已经切割完成，准备切割数据集2...")
  d22<-d1_d2[,c("SNP",beta_col,se_col,"chr.exposure","pos.exposure",a1_col,a2_col,"pval.exposure")]
  colnames(d22)<-c("SNP","OR","SE","CHR","BP","A1","A2","P")
  d22$OR<-exp(d22$OR)
  N<-22
  pb<-progress_bar$new(total=N)
  for(i in 1:N){
    d_chr<-subset(d22,CHR==i)
    dir.create(filename2)
    write.table(d_chr,paste0(getwd(),"/",filename2,"/catlog肺癌GWAS",i,".txt"),quote=F,row.names=F)
    pb$tick()

  }
  message("数据集2已经切割完成，准备进入GWAS的meta分析...")

  #GWAS的meta分析
  N<-22
  pb<-progress_bar$new(total=N)
  for (i in 1:N) {
    dir.create("GWAS_meta")
    pl<-paste0("plink/plink --meta-analysis ",
               paste0(getwd(),"/",filename1,"/芬兰肺癌GWAS",i,".txt "),
               paste0(getwd(),"/",filename2,"/catlog肺癌GWAS",i,".txt"),paste0(" --out GWAS_meta/gwas_meta",i))
    system(pl)
    pb$tick()
  }

  message("GWAS的meta分析已经完成，细节文件保存至GWAS_meta文件夹下,开始汇总分文件...")
  file<-dir(paste0(getwd(),"/GWAS_meta"))%>%data.frame()
  #输出GWAS meta
  N<-22
  pb<-progress_bar$new(total=N)
  ET<-c()
  for (i in 1:N) {
    ET0<-fread(paste0(getwd(),"/GWAS_meta/gwas_meta",i,".meta"))%>%data.frame()
    ET<-rbind(ET,ET0)
    pb$tick()
  }

  message("汇总分文件已经完成，正在输出结果....")
  ET$P.fin<-ifelse(ET$I>0.5,ET$P.R.,ET$P)
  ET$OR.fin<-ifelse(ET$I>0.5,ET$OR.R.,ET$OR)
  ET$EFFECT.fin<-log(ET$OR.fin)
  ET$SE.fin<-sqrt(((ET$EFFECT.fin)^2)/qchisq(ET$P.fin,1,lower.tail=F))
  ET_finn<-ET[,c("CHR","BP","SNP","A1","A2","P.fin","EFFECT.fin","SE.fin")]
  head(ET_finn)

  ET_finn2<-merge(ET_finn,data1[,c("SNP","eaf.exposure")])
  colnames(ET_finn2)[9]<-"EAF.fin"
  ET_finn2$N.fin<-Nsample
  ET_finn2$id.fin<-triat
  dir.create("最终GWAS荟萃文件")
  write.table(ET_finn2,"最终GWAS荟萃文件/肺癌GWAS_META.txt",quote=F,row.names=F)
  message("结果已经输出在最终GWAS荟萃文件中")
    }else {
    warning("keyssh不正确,请联系管理员微信SFM19950928或DKYXS666获取密钥")
  }
}
