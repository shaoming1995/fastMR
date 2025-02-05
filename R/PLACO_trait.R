#' @title 用于两表型的多效性位点分析
#' @param keyssh 学号
#' @param exp_gwas 暴露的GWAS summary数据,包含必要的四列:SNP，beta.exposure,se.exposure,pval.exposure
#' @param out_gwas 结局的GWAS summary数据，包含必要的四列:SNP，beta.outcome,se.outcome,pval.outcome
#' @param p.threshold 过滤P值，默认5e-08
#' @param save_file 输出保存的文件
#' @export
PLACO_triat<-function(keyssh,exp_gwas,out_gwas,p.threshold=5e-08,save_file="PLACO"){
  RegistID_dat <- RegistID_dat
  RegistID_u <- subset(RegistID_dat, IK == keyssh)
  tempid <- paste0(keyssh, "_", Sys.info()["nodename"], "_",
                   RegistID_u$RegistID)
  if (RegistID_u$FINN %in% tempid) {
  ############################################
  #---------------- Function for normal product based tail probability calculation
  # (Using modified Bessel function of the 2nd kind with order 0)
  .pdfx<-function(x) besselK(x=abs(x),nu=0)/pi
  .p.bessel<-function(z, varz, AbsTol=1e-13){
    p1<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[1])),Inf, abs.tol=AbsTol)$value)
    p2<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[2])),Inf, abs.tol=AbsTol)$value)
    p0<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]),Inf, abs.tol=AbsTol)$value)
    pval.compnull<-p1+p2-p0
    return(pval.compnull)
  }

  #---------------- Function for estimating the variances for PLACO
  var.placo<-function(Z.matrix, P.matrix, p.threshold=1e-4){
    # Here Z.matrix is the pxk matrix of Z-scores where p is total no. of variants in the dataset, and k is the no. of traits
    # Similarly, P.matrix is the corresponding pxk matrix of p-values where p is total no. of variants in the dataset, and k is the no. of traits
    # p.threshold determines which variants are deemed to have no association marginally (default: 1e-4)
    # checks
    k<-ncol(Z.matrix)
    if(k!=2) stop("This method is meant for 2 traits only. Columns correspond to traits.")
    ZP<-cbind(Z.matrix,P.matrix)
    ZP<-na.omit(ZP)

    rows.alt<-which(ZP[,3]<p.threshold & ZP[,4]<p.threshold)
    if(length(rows.alt)>0){
      ZP<-ZP[-rows.alt,]
      if(nrow(ZP)==0) stop(paste("No 'null' variant left at p-value threshold",p.threshold))
      if(nrow(ZP)<30) warning(paste("Too few 'null' variants at p-value threshold",p.threshold))
    }
    varz<-diag(var(ZP[,c(1,2)]))
    return(varz)
  }

  #---------------- Function for estimating correlation matrix of the Z's
  cor.pearson<-function(Z.matrix, P.matrix, p.threshold=1e-4){
    # Here Z.matrix is the pxk matrix of Z-scores where p is total no. of variants in the dataset, and k is the no. of traits
    # Similarly, P.matrix is the corresponding pxk matrix of p-values where p is total no. of variants in the dataset, and k is the no. of traits
    # p.threshold determines which variants are deemed to have no association marginally (default: 1e-4)
    # checks
    k<-ncol(Z.matrix)
    if(k!=2) stop("This method is meant for 2 traits only.")
    # estimating correlation
    row.exclude<-which( apply(P.matrix, MARGIN = 1, function(x) any(x < p.threshold)) == TRUE )
    if(length(row.exclude)>0) Z.matrix<-Z.matrix[-row.exclude,]
    R<-cor(Z.matrix)
    return(R)
  }

  ############################################
  placo<-function(Z, VarZ, AbsTol=.Machine$double.eps^0.8){
    # Z: vector of Z-scores of size k=2 (i.e., collection of Z-scores of a particular SNP for k=2 traits)
    # VarZ: vector of variances of Z-scores (covariance assumed 0; so need to be independent traits)
    # AbsTol: absolute tolerance (accuracy paramater) for numerical integration.
    # checks
    k<-length(Z)
    if(k!=2) stop("This method is meant for 2 traits only.")
    if(length(VarZ)!=k) stop("Provide variance estimates for 2 traits as obtained using var.placo() function.")

    # test of pleiotropy: PLACO
    pvalue.b=.p.bessel(z=Z, varz=VarZ, AbsTol=AbsTol)
    return(list(T.placo=prod(Z), p.placo=pvalue.b))
  }
message("完成配置，正在进行PLACO分析，耗时较长...."）
  combin_dat<-merge(exp_gwas,out_gwas,by="SNP",all=F)
  combin_dat$Zexp<-combin_dat$beta.exposure/combin_dat$se.exposure
  combin_dat$Zout<-combin_dat$beta.outcome/combin_dat$se.outcome
  Z.matrix<-combin_dat[,c("SNP","Zexp","Zout")]
  rownames(Z.matrix)<-Z.matrix$SNP
  Z.matrix<-Z.matrix[,-1]
  Z.matrix<-as.matrix(Z.matrix)
  P.matrix<-combin_dat[,c("SNP","pval.exposure","pval.outcome")]
  rownames(P.matrix)<-P.matrix$SNP
  P.matrix<-P.matrix[,-1]
  P.matrix<-as.matrix(P.matrix)

  # Step 1: Obtain the variance parameter estimates (only once)
  VarZ <- var.placo(Z.matrix, P.matrix, p.threshold=p.threshold)
  # Step 2: Apply test of pleiotropy for each variant
  out <- sapply(1:nrow(Z.matrix), function(i) placo(Z=Z.matrix[i,], VarZ=VarZ))
  # Check the output for say variant 100
  out1<-c()
  for (i in nrow(Z.matrix):1){
    out0<-cbind(out[,i]$T.placo,out[,i]$p.placo)
    out1<-rbind(out0,out1)
  }
  out1<-data.frame(out1)
  colnames(out1)<-c("T.placo","p.placo")
  out1$SNP<-rownames(Z.matrix)
  out1$FDR_placo<-p.adjust(out1$p.placo)
  write.table(out1,paste0(getwd(),"/",save_file,"PLACO.txt"),row.names = F)
  }
else {
  warning("keyssh不正确,请联系管理员微信SFM19950928或DKYXS666获取密钥")}
}

