#' @title 用于查找SNP附件的基因
#' @param data 输入查找的SNP，包含必要的三列SNP,chr.exposure,pos.exposure
#' @param flanking 输入SNP距离基因的大小
#' @param build 输入SNP参考的基因组位置，默认hg19，其他hg18或者hg38
#' @param filename 输出结果的文件夹名称
#' @export
#' @examples
Find_nearest_gene<-function(data,flanking = 0, build = "hg19",filename="匹配基因文件"){

  data1<-data[,c("SNP","chr.exposure","pos.exposure")]#染色体 位置 SNP
  result<-find_nearest_gene(data1, flanking = flanking, build = build,

                            collapse = TRUE, snp = "SNP", chr = "chr.exposure",

                            bp = "pos.exposure")

  D<-merge(result,data,by.x = "SNP",by.y="SNP",all = F)
  dir.create(filename)
  path<-paste0(getwd(),"/",filename,"/SNP_gene.csv")
  write.csv(D,path,row.names = F)
}




