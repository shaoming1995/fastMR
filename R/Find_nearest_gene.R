#' @title 用于查找SNP附件的基因
#' @param data 输入查找的SNP，包含必要的三列SNP,chr.exposure,pos.exposure
#' @param flanking 输入SNP距离基因的大小
#' @param build 输入SNP参考的基因组位置，默认hg19，其他hg18或者hg38
#' @param snp_col 输入SNP列名
#' @param chr_col 输入染色体列名
#' @param bp_col 输入基因组位置列名
#' @param filename 输出结果的文件夹名称
#' @export
Find_nearest_gene<-function(data,flanking = 0, snp_col='rsid',chr_col='chromosome', bp_col='position',build = "hg19",filename="匹配基因文件"){
  library(dplyr)
  data1<-data[,c(snp_col,chr_col,bp_col)]#染色体 位置 SNP
  find_nearest_gene1 <-function(data, flanking=100, build='hg19', collapse=TRUE, snp='rsid', chr='chromosome', bp='position'){

    data <- data # not sure this is needed

    if(build == 'hg18'){
      genelist <- hg18genelist
    }

    if(build == 'hg19'){
      genelist <- hg19genelist
    }

    if(build == 'hg38'){
      genelist <- hg38genelist
    }

    if(flanking != 100){
      flanking = as.numeric(flanking)
    }

    if(snp != 'rsid'){
      data <- data %>% rename(rsid = snp)
    }

    if(chr != 'chromosome'){
      data <- data %>% rename(chromosome = chr)
    }

    if(bp != 'position'){
      data <- data %>% rename(position = bp)
    }

    data<-sqldf::sqldf(sprintf("select A.*,B.* from
              data A left join genelist B
              ON (A.chromosome == B.CHR and
              A.position >= B.START - %1$s and
              A.position <= B.STOP + %1$s)", flanking*1000)
    )
    if (collapse==TRUE){
      data %>%
        dplyr::group_by(rsid, chromosome, position) %>%
        dplyr::summarise(GENES = paste(GENE, collapse=',')) %>%
        data.frame
    } else {
      data <- data %>%
        dplyr::rename(geneSTART = START, geneSTOP = STOP) %>%
        dplyr::select(rsid, chromosome, position, geneSTART,geneSTOP,GENE)

      data$distance <- apply(data,1, FUN=function(x){
        ifelse(
          !is.na(x['GENE']) & x['position'] < x['geneSTART'], -(as.numeric(x['geneSTART']) - as.numeric(x['position'])),
          ifelse(!is.na(x['GENE']) & x['position'] > x['geneSTOP'], as.numeric(x['position']) - as.numeric(x['geneSTOP']),
                 ifelse(!is.na(x['GENE']) & (x['position'] > x['geneSTART']) & (x['position'] < x['geneSTOP']),'intergenic',NA))
        )
      })
      data
    }
}
 result<-find_nearest_gene1(data1, flanking = flanking, build = build,

                            collapse = TRUE, snp = snp_col, chr = chr_col,
                            
                            bp = bp_col)
  D<-merge(data,result,by.x="SNP",by.y = "rsid",all = F)
  dir.create(filename)
  path<-paste0(getwd(),"/",filename,"/SNP_gene.csv")
  write.csv(D,path,row.names = F)
}

