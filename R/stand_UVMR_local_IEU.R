#' @title 适用暴露来自本地数据结局来自IEU的标准单变量孟德尔随机化分析
#' @param keyssh 密钥
#' @param expgwas 输入暴露的GWAS摘要数据
#' @param GWASID 输入结局的GWAS ID号
#' @param samplesize_outcome 输入结局数据的样本量，默认100000
#' @param name_outcome 输入结局的名称，默认outcome
#' @param clump_p1 输入工具变量的选择P值,默认5e-08
#' @param clump_r2 输入工具变量的选择的r2,默认0.001
#' @param clump_kb 输入工具变量的选择的距离,默认10000
#' @param pop 输入工具变量的选择的人群,默认EUR
#' @param outfile 输入分析结果的文件夹
#' @param steiger 是否进行反向过滤，默认是TURE
#' @param Fvalue 是否计算F值，默认是TURE
#' @param confounding_SNP 是否查询混杂因素，默认是TURE
#' @param local_clump 是否启动本地聚类，默认是FALSE
#' @param pt 是否进行绘图，默认是TURE
#' @param presso 是否进行MRPRESSO，默认是FALSE
#' @export

stand_UVMR_local_IEU<-function(keyssh,expgwas,GWASID,samplesize_outcome=100000,name_outcome="outcome",
                                 local_clump=F,confounding_SNP=NULL,clump_p1=5e-08,clump_r2=0.001,clump_kb=10000,pop="EUR",outfile="MR结果",presso=F,
                                 steiger=T,Fvalue=T,pt=T){

   library(tidyr)
  RegistID_dat <- RegistID_dat
  RegistID_u <- subset(RegistID_dat, IK == keyssh)
  tempid <- paste0(keyssh, "_", Sys.info()["nodename"], "_",RegistID_u$RegistID)
  if (RegistID_u$FINN %in% tempid) {
  dir.create(outfile)
  EXP<-expgwas[,c("SNP",
                  "effect_allele.exposure",
                  "other_allele.exposure",
                  "eaf.exposure",
                  "beta.exposure",
                  "se.exposure",
                  "pval.exposure",
                  "id.exposure",
                  "exposure",
                  "samplesize.exposure"
  )]
  expiv<-subset(EXP,pval.exposure<clump_p1)
  if(local_clump==F){
    expiv<- clump_data(expiv,clump_kb = clump_kb,clump_r2 = clump_r2,clump_p1 = clump_p1,clump_p2 = 1,pop = pop)}else{
      local_clump_data1<-function(temp_dat,pop,clump_kb,clump_r2){
        temp_dat$rsid <- temp_dat$SNP
        temp_dat$id <- temp_dat$id.exposure
        temp_dat$pval <- temp_dat$pval.exposure
        ld_sofm1<- function (dat, clump_kb, clump_r2, clump_p, bfile, plink_bin)
        {
          shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd",
                          "sh")
          fn <- tempfile()
          write.table(data.frame(SNP = dat[["rsid"]], P = dat[["pval"]]),
                      file = fn, row.names = F, col.names = T, quote = F)
          fun2 <- paste0(shQuote(plink_bin, type = shell), " --bfile ",
                         shQuote(bfile, type = shell), " --clump ", shQuote(fn,
                                                                            type = shell), " --clump-p1 ", clump_p, " --clump-r2 ",
                         clump_r2, " --clump-kb ", clump_kb, " --threads 20 --out ", shQuote(fn,
                                                                                             type = shell))
          system(fun2)
          res <- read.table(paste(fn, ".clumps", sep = ""), header = F)
          unlink(paste(fn, "*", sep = ""))
          y <- subset(dat, !dat[["rsid"]] %in% res[["V3"]])
          if (nrow(y) > 0) {
            message("Removing ", length(y[["rsid"]]), " of ", nrow(dat),
                    " variants due to LD with other variants or absence from LD reference panel")
          }
          return(subset(dat, dat[["rsid"]] %in% res[["V3"]]))
        }

        ld_sofm2<-function (dat = NULL, clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99,
                            pop = "EUR", access_token = NULL, bfile = NULL, plink_bin = NULL)
        {
          stopifnot("rsid" %in% names(dat))
          stopifnot(is.data.frame(dat))
          if (is.null(bfile)) {
            message("Please look at vignettes for options on running this locally if you need to run many instances of this command.")
          }
          if (!"pval" %in% names(dat)) {
            if ("p" %in% names(dat)) {
              warning("No 'pval' column found in dat object. Using 'p' column.")
              dat[["pval"]] <- dat[["p"]]
            }
            else {
              warning("No 'pval' column found in dat object. Setting p-values for all SNPs to clump_p parameter.")
              dat[["pval"]] <- clump_p
            }
          }
          if (!"id" %in% names(dat)) {
            dat$id <- random_string(1)
          }
          if (is.null(bfile)) {
            access_token = check_access_token()
          }
          ids <- unique(dat[["id"]])
          res <- list()
          for (i in 1:length(ids)) {
            x <- subset(dat, dat[["id"]] == ids[i])
            if (nrow(x) == 1) {
              message("Only one SNP for ", ids[i])
              res[[i]] <- x
            }
            else {
              message("Clumping ", ids[i], ", ", nrow(x), " variants, using ",
                      pop, " population reference")
              if (is.null(bfile)) {
                res[[i]] <- ld_clump_api(x, clump_kb = clump_kb,
                                         clump_r2 = clump_r2, clump_p = clump_p, pop = pop,
                                         access_token = access_token)
              }
              else {
                res[[i]] <- ld_sofm1(x, clump_kb = clump_kb,
                                     clump_r2 = clump_r2, clump_p = clump_p, bfile = bfile,
                                     plink_bin = plink_bin)
              }
            }
          }
          res <- dplyr::bind_rows(res)
          return(res)
        }
        filepath1<-paste0(getwd(),"/1kg.v3/",pop)# /home/Refrence/EUR/ 改成你参考文件的路径
        filepath2<-paste0(getwd(),"/1kg.v3/plink2_win64_20231212/plink2")
        filepath3<-paste0(getwd(),"/1kg.v3/plink2_mac_20231212/plink2")
        if(Sys.info()["sysname"] == "Windows"){
        temp_dat <- ld_sofm2(temp_dat,
                             plink_bin = filepath2, # /home/plink2.0/plink2改成你plink文件的路径
                             bfile = filepath1,
                             clump_kb = clump_kb, clump_r2 = clump_r2,pop=pop)}else{
                              temp_dat <- ld_sofm2(temp_dat,
                                                    plink_bin = filepath3, # /home/plink2.0/plink2改成你plink文件的路径
                                                    bfile = filepath1,
                                                    clump_kb = clump_kb, clump_r2 = clump_r2,pop=pop)

                             }

        temp_dat$rsid <- NULL
        temp_dat$pval <- NULL
        return(temp_dat)
      }
      expiv<- local_clump_data1(expiv,clump_kb = clump_kb,clump_r2 = clump_r2,pop = pop)
    }
  if(Fvalue==T){
   if(class(expiv$eaf.exposure[1])!="logical"){
        expiv$R2<-expiv$beta.exposure*expiv$beta.exposure*2*(expiv$eaf.exposure)*(1-expiv$eaf.exposure)
        expiv$Fvalue<-(expiv$samplesize.exposure-2)*expiv$R2/(1-expiv$R2)
        expiv<-subset(expiv,Fvalue>10)}else{
          expiv$R2<-NA
          expiv$Fvalue<-(expiv$beta.exposure/expiv$se.exposure)*(expiv$beta.exposure/expiv$se.exposure)
          expiv<-subset(expiv,Fvalue>10)
        }}
  if(dim(expiv)[[1]]!=0){
    OUT<-extract_outcome_data(snps=expiv$SNP,outcomes=GWASID,proxies=T,maf_threshold = 0.01)
    OUT$id.outcome=name_outcome
  OUT$outcome=name_outcome
  OUT$samplesize.outcome=samplesize_outcome
  OUT<-subset(OUT,pval.outcome>5e-08)
  OUT<-OUT[!duplicated(OUT$SNP),]
  total1<-merge(OUT,expiv,by.x="SNP",by.y="SNP",all = F)
      if(dim(total1)[[1]]!=0){
        #confonding_name<-c("Whole body fat mass","Arm fat mass left")
        # Atemp<- read.csv(path00,header = T,row.names = 1)  #PhenoScanSNP1(dim(total1)[[1]])
        # Atemp0<-Atemp%>%filter(trait %in% confonding_name)
        # Atemp0<-Atemp0[!duplicated(Atemp0$snp),]
        total1<-total1%>% filter(!SNP %in%confounding_SNP)
        #分别取出暴露与结局的数据
        EXP1<-total1[,c("SNP","effect_allele.exposure","other_allele.exposure", "eaf.exposure",
                        "beta.exposure","se.exposure", "pval.exposure","id.exposure","exposure",
                        "samplesize.exposure")]
        OUT1<-total1[,c("SNP","effect_allele.outcome","other_allele.outcome", "eaf.outcome",
                        "beta.outcome","se.outcome","pval.outcome","id.outcome","outcome",
                        "samplesize.outcome")]
        #去除回文
        dat1<-harmonise_data(exposure_dat=EXP1,outcome_dat=OUT1,action=2)
        #Steiger过滤
        if(steiger==T){
          dat1<-steiger_filtering(dat1)
          dat1<-subset(dat1,steiger_dir==TRUE)}
        res <- mr(dat1)
        mr_OR<-generate_odds_ratios(res)
        mr_OR$or<-round(mr_OR$or,3)
        mr_OR$or_lci95<-round(mr_OR$or_lci95,3)
        mr_OR$or_uci95 <- round(mr_OR$or_uci95,3)
        mr_OR$OR_CI <- paste0(mr_OR$or,"(",mr_OR$or_lci95,"-",mr_OR$or_uci95,")")
        het <- mr_heterogeneity(dat1)
        ple <- mr_pleiotropy_test(dat1)
        data_h_TableS1 <- dat1
        data_h_TableS1$R2<-data_h_TableS1$beta.exposure*data_h_TableS1$beta.exposure*2*(data_h_TableS1$eaf.exposure)*(1-data_h_TableS1$eaf.exposure)
        data_h_TableS1$Fvalue<-(data_h_TableS1$samplesize.exposure-2)*data_h_TableS1$R2/(1-data_h_TableS1$R2)
        path1<-paste0(getwd(),"/",outfile,"/MR.csv")
        path2<-paste0(getwd(),"/",outfile,"/het.csv")
        path3<-paste0(getwd(),"/",outfile,"/ple.csv")
        path4<-paste0(getwd(),"/",outfile,"/IV.csv")
        write.csv(mr_OR, path1, row.names = F)
        write.csv(het, path2, row.names = F)
        write.csv(ple, path3, row.names = F)
        write.csv(data_h_TableS1, path4, row.names = F)
        if(pt==T){
          #散点图
          path5<-paste0(getwd(),"/",outfile,"/散点图.pdf")
          pdf(path5, width = 10, height = 10)
          p1 <- mr_scatter_plot(res[1:5,], dat1)
          print(p1[[1]])
          dev.off()
          #敏感性分析Leave-one-out plot
          path6<-paste0(getwd(),"/",outfile,"/留一法.pdf")
          pdf(path6, width = 10, height = 10)
          mr_outcome_loo <- mr_leaveoneout(dat1)
          write.csv(mr_outcome_loo, paste0(getwd(),"/",outfile,"/留一法.csv"))
          p3 <- mr_leaveoneout_plot(mr_outcome_loo)
          print(p3[[1]])
          dev.off()
          #森林图
          path7<-paste0(getwd(),"/",outfile,"/森林图.pdf")
          pdf(path7, width = 10, height = 10)
          mr_outcome_single <- mr_singlesnp(dat1)
          write.csv(mr_outcome_single, paste0(getwd(),"/",outfile,"/森林图.csv"))
          p2 <- mr_forest_plot(mr_outcome_single)
          print(p2[[1]])
          dev.off()
          #漏斗图
          path8<-paste0(getwd(),"/",outfile,"/漏斗图.pdf")
          pdf(path8,width = 10,height = 10)
          mr_outcome_single <- mr_singlesnp(dat1)
          write.csv(mr_outcome_single, paste0(getwd(),"/",outfile,"/漏斗图.csv"))
          p4 <- mr_funnel_plot(mr_outcome_single)
          print(p4[[1]])
          dev.off()}
        if(presso==T){
          mrpresso_data<-MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",BetaExposure="beta.exposure",
                                             SdOutcome = "se.outcome", OUTLIERtest = T, DISTORTIONtest = T,SdExposure = "se.exposure",data=dat1)
          res_presso_main=mrpresso_data[["Main MR results"]]
          pathpre1<-paste0(getwd(),"/",outfile,"/res_presso.csv")
          write.csv(res_presso_main, pathpre1, row.names = F)
          #Residual Sum of Squares of the observed data观察数据的残差平方和
          res_mrpresso=data.frame(RSSobs=mrpresso_data[["MR-PRESSO results"]][["Global Test"]][["RSSobs"]],
                                  Pvalue=mrpresso_data[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]])
          pathpre2<-paste0(getwd(),"/",outfile,"/RSSobs.csv")
          write.csv(res_mrpresso, pathpre2, row.names = F)}
      }else{cat("当前阈值可能严格，未找到工具变量")}

  }
  else{cat("当前阈值可能严格，未找到工具变量")}
  warning("此R包由作者邵明编制，请关注抖音号793742981或者顶刊研习社公众号")
  }
  else {
    warning("keyssh不正确,请联系管理员微信SFM19950928或DKYXS666获取密钥")
  }
}
