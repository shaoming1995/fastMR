#' @title 适用qtl与表型的SMR分析
#' @param keyssh 学号
#' @param pop 参考人种，默认EUR
#' @param maGWASfile 输入表型gwas数据
#' @param eqtlfile 输入qtl数据
#' @param single 是否top_SNP计算SMR，默认T
#' @param kb top_SNP范围
#' @param r2 top_SNP与区域LD的r2
#' @param outfile 保存文件名称
#' @param num 计算机启动线程最大数，默认10
#' @export
SMR_qtl_GWAS<-function(keyssh,pop="EUR",maGWASfile,eqtlfile,single=T,kb=500,r2=0.1,outfile,num=10){
  RegistID_dat <- RegistID_dat
  RegistID_u <- subset(RegistID_dat, IK == keyssh)
  tempid <- paste0(keyssh, "_", Sys.info()["nodename"], "_",
                   RegistID_u$RegistID)
  if (RegistID_u$FINN %in% tempid) {
    if(Sys.info()["sysname"] == "Windows"){
      if(single==T){
        shell<-paste0(getwd(),"/SMR/smr_Win/","smr-1.3.1-win.exe --bfile ",getwd(),"/1kg.v3/",pop," --gwas-summary ",maGWASfile," --beqtl-summary ",eqtlfile," --out ",outfile," --thread-num ",num)
        system(shell)}else{
          shell<-paste0(getwd(),"/SMR/smr_Win/","smr-1.3.1-win.exe --bfile ",getwd(),"/1kg.v3/",pop," --gwas-summary ",maGWASfile," --beqtl-summary ",eqtlfile," --out ",outfile," --smr-multi"," --set-wind ",kb," --ld-multi-snp ",r2," --thread-num ",num)
          system(shell)}}else{
            if(single==T){
              shell<-paste0(getwd(),"/SMR/smr_Mac/","smr-1.3.1-macos-arm64 --bfile ",getwd(),"/1kg.v3/",pop," --gwas-summary ",maGWASfile," --beqtl-summary ",eqtlfile," --out ",outfile," --thread-num ",num)
              system(shell)}else{
                shell<-paste0(getwd(),"/SMR/smr_Mac/","smr-1.3.1-macos-arm64 --bfile ",getwd(),"/1kg.v3/",pop," --gwas-summary ",maGWASfile," --beqtl-summary ",eqtlfile," --out ",outfile," --smr-multi"," --set-wind ",kb," --ld-multi-snp ",r2," --thread-num ",num)
                system(shell)}}}else {
                  warning("keyssh不正确,请联系管理员微信SFM19950928或DKYXS666获取密钥")}
}
