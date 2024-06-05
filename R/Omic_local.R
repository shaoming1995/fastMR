#' @title 组学孟德尔随机化
#' @param keyssh 密钥
#' @param fac_cell_met 选择运行那种组学,1,2,3代表炎症因子,免疫细胞,代谢
#' @param savefile 结果输出保存的文件夹
#' @param omicfile 输入组学数据存放的文件夹
#' @param finish_out 输入结局的GWAS summary
#' @param local_clump 是否启动本地聚类,默认不启动
#' @param clump_p1 输入工具变量的选择P值,默认1e-05
#' @param clump_r2 输入工具变量的选择的r2,默认0.001
#' @param clump_kb 输入工具变量的选择的距离,默认10000
#' @param presso 是否启动MRPRESSO,默认F
#' @param pop 输入工具变量的选择的人群,默认EUR
#' @export
Omic_local<-function(keyssh,fac_cell_met=1,savefile="MR结果",omicfile,finish_out,local_clump=F,clump_p1=1e-05,clump_r2=0.001,clump_kb=10000,pop="EUR",presso=F){
  A_temp <- c()#
  B_temp <- c()#
  C_temp <- c()#
  D_temp <- c()
  E_temp <- c()#
  F_temp <- c()#
  G_temp <- c()#
  G_temp1 <- c()#
  H_temp1 <- c()#
  H_temp2 <- c()#
  H_temp3 <- c()#
  H_temp4 <- c()#
  if(fac_cell_met==1){
    N=91
    nm<-"炎症因子"
  }
  if(fac_cell_met==2){
    N=731
    nm<-"免疫细胞"
  }
  if(fac_cell_met==3){
    N=1400
    nm<-"血浆代谢"
  }
  if(fac_cell_met==4){
    N=211
    nm<-"肠道菌群"
  }
  if(presso==F){
  if (Sys.info()["nodename"] == keyssh) {
    dir.create(savefile)
    file<-dir(omicfile)
    file<-data.frame(file)
    #file<-file[-1,]
    file<-data.frame(file)
    for (id in file[,1]){
      exppath<-paste0(omicfile,"/",id)
      exp<-fread(exppath,header=T)%>%data.frame()
      expiv<-subset(exp,pval.exposure<clump_p1)
      if(local_clump==F){
        expiv<- clump_data(expiv,clump_kb = clump_kb,clump_r2 = clump_r2,clump_p1 = 1,clump_p2 = 1,pop = pop)}else{
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
          expiv<- local_clump_data(expiv,clump_kb = clump_kb,clump_r2 = clump_r2,pop = pop)
        }
        expiv$R2<-expiv$beta.exposure*expiv$beta.exposure*2*(expiv$eaf.exposure)*(1-expiv$eaf.exposure)
        expiv$Fvalue<-(expiv$samplesize.exposure-2)*expiv$R2/(1-expiv$R2)
        expiv<-subset(expiv,Fvalue>10)
      if(dim(expiv)[[1]]!=0){
        #在结局GWAS summary中寻找与暴露对应的SNPs
        total1<-merge(finish_out,expiv,by.x="SNP",by.y="SNP",all = F)
        #去除与结局有gwas显著性的SNPs以及可能重复的SNP
        total1<-subset(total1,pval.outcome>5e-08)
        total1<-total1[!duplicated(total1$SNP),]
        if(dim(total1)[[1]]!=0){
          #分别取出暴露与结局的数据
          EXP1<-total1[,c("SNP","effect_allele.exposure","other_allele.exposure", "eaf.exposure",
                          "beta.exposure","se.exposure", "pval.exposure","id.exposure","exposure",
                          "samplesize.exposure")]
          OUT1<-total1[,c("SNP","effect_allele.outcome","other_allele.outcome", "eaf.outcome",
                          "beta.outcome","se.outcome","pval.outcome","id.outcome","outcome",
                          "samplesize.outcome")]
          #去除回文
          dat1<-harmonise_data(exposure_dat=EXP1,outcome_dat=OUT1,action=2)
          test1<-try(data_h_F10<-dat1%>%subset(dat1$mr_keep==TRUE))#F10+steiger
          if(class(test1)!="try-error"&dim(data_h_F10)[[1]]!=0){
            #Steiger过滤
            data_h_F10_steiger<-steiger_filtering(data_h_F10)
            data_h_F10_steiger<-subset(data_h_F10_steiger,steiger_dir==TRUE)
            data_h_F10_steiger<-steiger_filtering(data_h_F10)
            data_h_F10_steiger<-subset(data_h_F10_steiger,steiger_dir==TRUE)
            if(dim(data_h_F10_steiger)[[1]]!=0){
              res <- mr(data_h_F10_steiger)
              res$steiger<- "yes"
              res$F10<-"yes"
              res$id<-id
              res$IV<-paste0("threshold value <",clump_p1)
              mr_OR<-generate_odds_ratios(res)
              mr_OR$or<-round(mr_OR$or,3)
              mr_OR$or_lci95<-round(mr_OR$or_lci95,3)
              mr_OR$or_uci95 <- round(mr_OR$or_uci95,3)
              mr_OR$OR_CI <- paste0(mr_OR$or,"(",mr_OR$or_lci95,"-",mr_OR$or_uci95,")")
              het <- mr_heterogeneity(data_h_F10_steiger)
              #het$file<-id
              ple <- mr_pleiotropy_test(data_h_F10_steiger)
              #ple$file<-id
              data_h_TableS1 <- data_h_F10_steiger
              data_h_TableS1$R2<-data_h_TableS1$beta.exposure*data_h_TableS1$beta.exposure*2*(data_h_TableS1$eaf.exposure)*(1-data_h_TableS1$eaf.exposure)
              data_h_TableS1$Fvalue<-(data_h_TableS1$samplesize.exposure-2)*data_h_TableS1$R2/(1-data_h_TableS1$R2)
              data_h_TableS1$steiger<-"yes"
              data_h_TableS1$F10<-"yes"
              data_h_TableS1$file<-id
              A_temp <- rbind(mr_OR, A_temp)
              B_temp <- rbind(het, B_temp)
              C_temp <- rbind(ple, C_temp)
              D_temp <- rbind(data_h_TableS1, D_temp)

              pt1<-paste0(getwd(),"/",savefile,"/MR.csv")
              pt2<-paste0(getwd(),"/",savefile,"/het.csv")
              pt3<-paste0(getwd(),"/",savefile,"/ple.csv")
              pt4<-paste0(getwd(),"/",savefile,"/IV.csv")
              write.csv(A_temp,pt1, row.names = F)
              write.csv(B_temp,pt2, row.names = F)
              write.csv(C_temp,pt3, row.names = F)
              write.csv(D_temp,pt4, row.names = F) }else{H_temp1 <- rbind(H_temp1, id)%>%data.frame()
              H_temp1$reason<-paste0("因不满足反向过滤被排除的",nm)
              pt5<-paste0(getwd(),"/",savefile,"/NOsteigerid.csv")
              write.csv(H_temp1,pt5, row.names = F)}
          }else{H_temp2 <- rbind(H_temp2, id)%>%data.frame()
          H_temp2$reason<-paste0("因不满足回文等被排除的",nm)
          pt6<-paste0(getwd(),"/",savefile,"/NOclumpid.csv")
          write.csv(H_temp2, pt6, row.names = F)}
          }else{H_temp3 <- rbind(H_temp3, id)%>%data.frame()
          H_temp3$reason<-paste0("因无法匹配到工具变量被排除的",nm)
          pt7<-paste0(getwd(),"/",savefile,"/NOmergeid.csv")
          write.csv(H_temp3, pt7, row.names = F)}
        }else{H_temp4 <- rbind(H_temp4, id)%>%data.frame()
          H_temp4$reason<-paste0("因不满足F值大于10被排除的",nm)
          pt8<-paste0(getwd(),"/",savefile,"/NOF10id.csv")
          write.csv(H_temp4, pt8, row.names = F)}
          print(id)
          row_numbers <- which(file$file == id)
          pv<-round((row_numbers/N)*100,4)
          cat("已完成",pv,"%")
        }}else{
    warning("keyssh不正确,请联系管理员微信SFM19950928或DKYXS666获取密钥")
        }
    }

else{
      library(MRPRESSO)
      #多效性偏差
      pt4<-paste0(getwd(),"/",savefile,"/IV.csv")
      exp<-read.csv(pt4)
      data_class<-unique(exp$id.exposure)
      PPT<-paste0("切分好的",nm,"暴露文件")
      dir.create(PPT)
      for (i in data_class){
        A<-subset(exp,id.exposure==i)
        cfilename = paste0(getwd(),"/切分好的",nm,"暴露文件/",i,".csv")
        write.csv(A,cfilename,row.names = F)
        }
      message("暴露已经完成切分,presso准备中...")
      pt4<-paste0(getwd(),"/切分好的",nm,"暴露文件/")

      file<-dir(pt4)
      file<-data.frame(file)
      #file<-file[-1,]
      file<-data.frame(file)
      for(ic in file[,1]){
        presso_path<-paste0(getwd(),"/切分好的",nm,"暴露文件/",ic)
        data_h_F10_steiger<-read.csv(presso_path,header = T)
       if (dim(data_h_F10_steiger)[[1]]>3){
      mrpresso_data<-mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                               SdOutcome = "se.outcome", SdExposure = "se.exposure",
                               data =data_h_F10_steiger,
                               OUTLIERtest = TRUE, DISTORTIONtest = F,SignifThreshold = 1)
      res_presso_main=mrpresso_data[["Main MR results"]]
      res_presso_main$id<-ic
      E_temp <- rbind(res_presso_main, E_temp)
      pt9<-paste0(getwd(),"/",savefile,"/result_presso.csv")
      write.csv(E_temp, pt9, row.names = F)
      #Residual Sum of Squares of the observed data观察数据的残差平方和
      res_mrpresso=data.frame(RSSobs=mrpresso_data[["MR-PRESSO results"]][["Global Test"]][["RSSobs"]],
                              Pvalue=mrpresso_data[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]])
      res_mrpresso$id<-ic
      F_temp <- rbind(res_mrpresso, F_temp)

      pt10<-paste0(getwd(),"/",savefile,"/Global_Test_RSSobs.csv")
      write.csv(F_temp,pt10, row.names = F)
      Outliersnp<-mrpresso_data$`MR-PRESSO results`$`Outlier Test`
      data_h_F10_steiger_snp<-data_h_F10_steiger[,c("SNP","id.exposure")]
      Outliersnp<-cbind(Outliersnp,data_h_F10_steiger_snp)
      G_temp1 <- rbind(Outliersnp,G_temp1)
      pt11<-paste0(getwd(),"/",savefile,"/single_Test_RSSobs.csv")
      write.csv(G_temp1, pt11, row.names = F)
      row_numbers <- which(file$file == ic)
      pv<-round((row_numbers/nrow(file))*100,4)
      cat("已完成",pv,"%")
      }else{
        G_temp <- rbind(G_temp, id)%>%data.frame()
        G_temp$reason<-"因工具变量少于4"
        write.csv(G_temp, "NOpressoid.csv", row.names = F)}
  }
  }
}

