#' @title 肠道菌群数据预处理
#' @param keyssh 密钥
#' @param inputfile 肠道菌群数据文件位置
#' @param savefile 预处理后数据保存位置
#' @param exp_or_out 处理成暴露还是结局，默认是T表示处理成暴露
#' @param P_exp 设定处理暴露用于潜在工具变量的P值
#' @param P_out 设定处理结局用于潜在工具变量的P值
#' @export
gut_pre<-function (keyssh,inputfile, savefile, exp_or_out = T, P_exp = 1e-05,
                   P_out = 5e-08)
{library(tidyr)
  RegistID_dat <- RegistID_dat
  RegistID_u <- subset(RegistID_dat, IK == keyssh)
  tempid <- paste0(keyssh, "_", Sys.info()["nodename"], "_",
                   RegistID_u$RegistID)
  if (RegistID_u$FINN %in% tempid) {
    library(tidyr)
    file <- dir(inputfile) %>% data.frame()
    for (i in 1:nrow(file)) {
      dir.create(savefile)
      path <- paste0(inputfile, "/", file[i, 1])
      data <- data.table::fread(path) %>% data.frame()
      colnames(data)[c( 4, 6,5,11, 7, 8,10)] <- c("SNP","effect_allele.exposure",
                                                  "other_allele.exposure", "samplesize.exposure",
                                                  "beta.exposure", "se.exposure", "pval.exposure")
      data$eaf.exposure<-NA
      data$id.exposure <- data$bac
      data$exposure <- data$id.exposure
      data <- data[, c("effect_allele.exposure", "other_allele.exposure",
                       "samplesize.exposure", "beta.exposure", "se.exposure",
                       "eaf.exposure", "pval.exposure", "SNP", "id.exposure",
                       "exposure")]
      if (exp_or_out == T) {
        data <- subset(data, pval.exposure < P_exp)
        pathe <- paste0(savefile, "/", file[i, 1], ".csv")
        write.csv(data, pathe, row.names = F, quote = F)
        pv <- round(i/nrow(file), 4) * 100
        message("已完成暴露数据转化", pv, "%")
      }
      else {
        data2 <- data[, c("effect_allele.exposure", "other_allele.exposure",
                          "samplesize.exposure", "beta.exposure", "se.exposure",
                          "eaf.exposure", "pval.exposure", "SNP", "id.exposure",
                          "exposure")]
        colnames(data2) <- c("effect_allele.outcome", "other_allele.outcome",
                             "samplesize.outcome", "beta.outcome", "se.outcome",
                             "eaf.outcome", "pval.outcome", "SNP", "id.outcome",
                             "outcome")

        pathe <- paste0(savefile, "/", file[i, 1])
        data2 <- subset(data2, pval.outcome > P_out)
        write.table(data2, file = pathe, row.names = F)
        pv <- round(i/nrow(file), 4) * 100
        cat("已完成结果数据转化", pv, "%")
      }
    }
  }else {
    warning("keyssh不正确,请联系管理员微信SFM19950928或DKYXS666获取密钥")
  }
}
