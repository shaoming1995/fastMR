#' @title 肠道菌群MR工具
#' @param name 学号
#' @param key 密码
#' @param savefile 设置一个保存结果的文件夹名字
#' @param PATH 切分的肠道菌群暴露文件位置
#' @param GWASsummay 本地结局的预处理好的GWAS summay
#' @param pop 输入工具变量的选择的人群,默认EUR
#' @param local_clump 是否启动本地聚类,默认不启动
#' @param kb 聚类距离
#' @param r2 相关系数
#' @param outname 本地结局的名称
#' @export
Gut_local<-function(name,key,savefile,PATH,GWASsummay,outname,local_clump=F,kb,r2,pop="EUR"){
  A <- name
  A <- as.numeric(gsub("DK", "00", A))
  C <- A + key
  if (C == 2310000) {
    library(TwoSampleMR)
    dir.create(savefile)
    filename <- data.frame(dir(PATH))
    A_temp <- c()
    B_temp <- c()
    C_temp <- c()
    D_temp <- c()
    E_temp <-c()
    F_temp <-c()
    for (i in filename[, 1]) {
      ipath <- paste0(PATH, "/", i)
      exp_temp <- read.csv(ipath, header = T)
      if(local_clump==F){
      test2 <- (try(exp_temp <- clump_data(exp_temp, clump_kb = kb,
                                           clump_r2 = r2)))}else{
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

                                             test2 <- (try(exp_temp <- local_clump_data1(exp_temp, clump_kb = kb,
                                                                                  clump_r2 = r2,pop=pop)))
                                           }

      if (class(test2) == "try-error") {
        cat(i,"由于网络502问题未完成clump")
        total <- merge(GWASsummay, exp_temp, by = "SNP")
        total$eaf.exposure <- NA   #修改这里，让从本地和IEU数据库获得的菌群数据一致
        exp3 <- total[, c("SNP", "effect_allele.exposure",
                          "other_allele.exposure", "beta.exposure", "se.exposure",
                          "pval.exposure", "id.exposure", "exposure",
                          "eaf.exposure")]
        out3 <- total[, c("SNP", "effect_allele.outcome",
                          "other_allele.outcome", "eaf.outcome", "beta.outcome",
                          "se.outcome", "pval.outcome", "id.outcome",
                          "outcome", "eaf.outcome")]
        dat <- harmonise_data(exposure_dat = exp3, outcome_dat = out3,
                              action = 2)
        data_h<-dat%>%subset(dat$mr_keep==TRUE)
        data_h$Fvalue <- (data_h$beta.exposure/data_h$se.exposure)*(data_h$beta.exposure/data_h$se.exposure)
        data_h_TableS1 <- data_h[, c("exposure","SNP","effect_allele.exposure", "other_allele.exposure",
                                     "beta.exposure", "se.exposure","Fvalue","pval.exposure",
                                     "beta.outcome","se.outcome", "pval.outcome")] #导出SNP表格用
        data_h_TableS1$cluster <- 0  #0表达没有富集成功
        res <- mr(data_h)
        res$cluster <- 0    #0表达没有富集成功
        mr_OR<-generate_odds_ratios(res)
        mr_OR$or<-round(mr_OR$or,3)
        mr_OR$or_lci95<-round(mr_OR$or_lci95,3)
        mr_OR$or_uci95 <- round(mr_OR$or_uci95,3)
        mr_OR$OR_CI <- paste0(mr_OR$or,"(",mr_OR$or_lci95,"-",mr_OR$or_uci95,")") #生成成一个类似1.475(0.806-2.701)的格式，方便复制用
      }
      else {
        total <- merge(GWASsummay, exp_temp, by = "SNP")
        if(class(total)=="NULL"){
          Ename <- paste0(savefile, "/", "肠道菌群与",
                        outname, "未匹配到工具变量的菌群ID.csv")
          E_temp <- rbind(i, E_temp)
          write.csv(E_temp, Ename, row.names = F)}else{
            total$eaf.exposure <- NA
        exp3 <- total[, c("SNP", "effect_allele.exposure",
                          "other_allele.exposure", "beta.exposure", "se.exposure",
                          "pval.exposure", "id.exposure", "exposure",
                          "eaf.exposure")]
        out3 <- total[, c("SNP", "effect_allele.outcome",
                          "other_allele.outcome", "eaf.outcome", "beta.outcome",
                          "se.outcome", "pval.outcome", "id.outcome",
                          "outcome", "eaf.outcome")]
        dat <- harmonise_data(exposure_dat = exp3, outcome_dat = out3,
                              action = 2)
        data_h<-dat%>%subset(dat$mr_keep==TRUE)
          print(data_h)
        data_h$Fvalue <- (data_h$beta.exposure/data_h$se.exposure)*(data_h$beta.exposure/data_h$se.exposure)
        data_h_TableS1 <- data_h[, c("exposure","SNP","effect_allele.exposure", "other_allele.exposure",
                                     "beta.exposure", "se.exposure","Fvalue","pval.exposure",
                                     "beta.outcome","se.outcome", "pval.outcome")]
        data_h_TableS1$cluster <- 1
        res<-mr(data_h)
        res$cluster <- 1
        mr_OR<-generate_odds_ratios(res)
        mr_OR$or<-round(mr_OR$or,3)
        mr_OR$or_lci95<-round(mr_OR$or_lci95,3)
        mr_OR$or_uci95 <- round(mr_OR$or_uci95,3)
        mr_OR$OR_CI <- paste0(mr_OR$or,"(",mr_OR$or_lci95,"-",mr_OR$or_uci95,")")
      }}
      if (dim(res)[[1]] != 0) {
        het <- mr_heterogeneity(dat)
        ple <- mr_pleiotropy_test(dat)
        A_temp <- rbind(mr_OR, A_temp)
        B_temp <- rbind(het, B_temp)
        C_temp <- rbind(ple, C_temp)
        D_temp <- rbind(data_h_TableS1, D_temp)
        print(paste0("当前运行到", i, "文件"))
        Aname <- paste0(savefile, "/", "肠道菌群与",
                        outname, "的MR结果.csv")
        Bname <- paste0(savefile, "/", "肠道菌群与",
                        outname, "的异质性结果.csv")
        Cname <- paste0(savefile, "/", "肠道菌群与",
                        outname, "的多效性结果.csv")
        Dname <- paste0(savefile, "/", "肠道菌群与",
                        outname, "的SNPs情况.csv")
        write.csv(A_temp, Aname, row.names = F)
        write.csv(B_temp, Bname, row.names = F)
        write.csv(C_temp, Cname, row.names = F)
        write.csv(D_temp, Dname, row.names = F)
      }
      else {
        cat("请前往切分好的肠道菌群暴露文件下删除",
            i, "文件再次重新运行")
      }
    }
    cat("当前分析已全部完成！请前往", savefile,
        "文件夹下查看结果")
  }
  else {
    cat("请联系客服获取账户学号和密码或微信联系SFM19950928")
  }
}
