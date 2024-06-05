#' @title 组学留一法
#' @param filepath IV.csv的绝对文件路径例如G:/MR结果
#' @param savefilepath 保存留一法的绝对路径例如G:/MR结果
#' @export
omic_LOO<-function(filepath,savefilepath=filepath){
F_temp <- c()#
G_temp <- c()#
library(TwoSampleMR)
pt4<-paste0(filepath,"/IV.csv")
exp<-read.csv(pt4)
data_class<-unique(exp$id.exposure)
PPT<-paste0(filepath,"/切分好的临时文件")
dir.create(PPT)
for (i in data_class){
  A<-subset(exp,id.exposure==i)
  cfilename = paste0(filepath,"/切分好的临时文件/",i,".csv")
  write.csv(A,cfilename,row.names = F)
}
message("暴露已经完成切分,留一法准备中...")
pt4<-paste0(getwd(),"/切分好的临时文件/")
file<-dir(pt4)
file<-data.frame(file)
#file<-file[-1,]
file<-data.frame(file)
for(ic in file[,1]){
  presso_path<-paste0(getwd(),"/切分好的临时文件/",ic)
  data_h_F10_steiger<-read.csv(presso_path,header = T)
  if (dim(data_h_F10_steiger)[[1]]>2){
    Fout<-mr_leaveoneout(data_h_F10_steiger)
    F_temp <- rbind(Fout, F_temp)
    pt9<-paste0(savefilepath,"/LOO.csv")
    write.csv(F_temp, pt9, row.names = F)
    row_numbers <- which(file$file == ic)
    pv<-round((row_numbers/nrow(file))*100,4)
    cat("已完成",pv,"%")}else{
      G_temp <- rbind(G_temp, id)%>%data.frame()
      G_temp$reason<-"Num_iv<3"
      write.csv(G_temp, "NOLOO.csv", row.names = F)}}
  }
