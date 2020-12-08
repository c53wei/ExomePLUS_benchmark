args=commandArgs(F)

## Rscript happy_metrics.R <INPUT_DIR> <OUTPUT_DIR>

exome<-T
rep_1<-F

if(exome && rep_1 ) 
{
  args[1]="/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_sample/illumina/happy_WES"
  args[2]="/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_sample/illumina/happy_WES"
  
  IDs<-c("HG002_WGS_30X", "HG002_XTV2_4plex", "HG002_XTV2_1plex", "HG002_XTV2_4plex", "HG002_XT", 
         "HG003_WGS_30X", "HG003_XTV2_1plex", "HG003_XTV2_4plex", "HG003_XT",
         "HG004_WGS_30X", "HG004_XTV2_1plex", "HG004_XTV2_4plex", "HG004_XT")
  
  protocol<-c("XT", "XTV2_1plex", "XTV2_4plex", 
              "WGS_30X")
  
  
  
  
} else if(!exome && rep_1) {
  args[1]="/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_sample/illumina/happy_WGS"
  args[2]="/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_sample/illumina/happy_WGS"
  
  IDs<-c("HG002_2X_DS", "HG002_30X", "HG002_30X_DS", "HG002_5X_DS", "HG002_2X", "HG002_5X",
         "HG003_2X_DS", "HG003_30X", "HG003_30X_DS", "HG003_5X_DS", "HG003_2X", "HG003_5X",
         "HG004_2X_DS", "HG004_30X", "HG004_30X_DS", "HG004_5X_DS", "HG004_2X", "HG004_5X")
  
  protocol<-c("2X_DS", "2X", "5X_DS", "5X", "30X_DS", "30X" )
  
} else if(exome && !rep_1) {
  args[1]<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_sample/illumina/technical_replicate/happy_WES"
  args[2]<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_sample/illumina/technical_replicate/happy_WES"
  
  IDs<-c("HG002_WGS_30X", "HG002_IDT", "HG002_IDT", "HG002_XTV2_4plex", "HG002_XT", "HG002_XTV2_1plex", "HG002_XTV2_4plex",   
         "HG003_WGS_30X", "HG003_IDT", "HG003_XT", "HG003_XTV2_1plex", "HG003_XTV2_4plex", 
         "HG004_WGS_30X", "HG004_IDT", "HG004_XT", "HG004_XTV2_1plex", "HG004_XTV2_4plex" )
  
  protocol<-c("IDT", "XT", "XTV2_1plex", "XTV2_4plex", 
              "WGS_30X")
  
} else {
  args[1]<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_sample/technical_replicate/happy_WGS"
  args[2]<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_sample/technical_replicate/happy_WGS"
  
  IDs<-c("HG002_30X", "HG002_2X", "HG002_5X",
         "HG003_30X", "HG003_2X", "HG003_5X",
         "HG004_30X", "HG004_2X", "HG004_5X")
  
  protocol<-c("2X", "5X","30X" )
}


input_dir<-args[1]
output_dir<-args[2]


library(ggplot2)
library(stringr)

filenames<-list.files(input_dir, full.names = T, pattern = '.csv')
files<-vector(mode="list", length=length(filenames))

snp_all<-vector(mode="list", length=length(filenames))
indel_all<-vector(mode="list", length=length(filenames))
snp_pass<-vector(mode="list", length=length(filenames))
indel_pass<-vector(mode="list", length=length(filenames))


metrics<-c("Recall", "Precision", "Fraction NA", "F1 Score")

for( i in 1:length(files))
{
  files[[i]]<-read.csv(filenames[[i]])
  temp<-split(files[[i]], paste(files[[i]]$Type, files[[i]]$Filter))
  snp_all[[i]]<-temp[["SNP ALL"]][c(11:14)]
  indel_all[[i]]<-temp[["INDEL ALL"]][c(11:14)]
  snp_pass[[i]]<-temp[["SNP PASS"]][c(11:14)]
  indel_pass[[i]]<-temp[["INDEL PASS"]][c(11:14)]
}
snp_all<-do.call(rbind, snp_all)
snp_all$IDs<-IDs
snp_all$Sample<-substr(snp_all$IDs, 1, 5)
snp_all$X<-str_split_fixed(snp_all$IDs, '_', 2)[,2]
snp_all$X<-factor(snp_all$X, levels=protocol)
rownames(snp_all)<-basename(filenames)

indel_all<-do.call(rbind, indel_all)
indel_all$IDs<-IDs
indel_all$Sample<-substr(indel_all$IDs, 1, 5)
indel_all$X<-str_split_fixed(indel_all$IDs, '_', 2)[,2]
indel_all$X<-factor(snp_all$X, levels=protocol)
rownames(indel_all)<-basename(filenames)

snp_pass<-do.call(rbind, snp_pass)
snp_pass$IDs<-IDs
snp_pass$Sample<-substr(snp_pass$IDs, 1, 5)
snp_pass$X<-str_split_fixed(snp_pass$IDs, '_', 2)[,2]
snp_pass$X<-factor(snp_pass$X, levels=protocol)
rownames(snp_pass)<-basename(filenames)

indel_pass<-do.call(rbind, indel_pass)
indel_pass$IDs<-IDs
indel_pass$Sample<-substr(indel_pass$IDs, 1, 5)
indel_pass$X<-str_split_fixed(indel_pass$IDs, '_', 2)[,2]
indel_pass$X<-factor(snp_pass$X, levels=protocol)
rownames(indel_pass)<-basename(filenames)


for(i in seq(1, 2, 1))
{
  graphics.off()
  print(ggplot(snp_all) +
          geom_point(aes(x=X, y=snp_all[,i], col=Sample), size=5) +
          labs(x="", y=metrics[i]) +
          ylim(0,1) +
          theme_bw() +
          theme(axis.text.x =element_text(size=18,face="bold"),
                axis.text.y =element_text(size=20),
                axis.title.y = element_text(size=22,face="bold"),
                legend.title=element_text(size=22,face="bold"),
                legend.text=element_text(size=20))
  )
  ggsave(paste(output_dir, '/SNP_ALL_', metrics[i], '.pdf', sep=''), device = 'pdf', width=11, height=8.5)
  
  graphics.off()
  print(ggplot(indel_all) +
          geom_point(aes(x=X, y=indel_all[,i], col=Sample), size=5) +
          labs(x="", y=metrics[i]) +
          ylim(0,1) +
          theme_bw() +
          theme(axis.text.x =element_text(size=18,face="bold"),
                axis.text.y =element_text(size=20),
                axis.title.y = element_text(size=22,face="bold"),
                legend.title=element_text(size=22,face="bold"),
                legend.text=element_text(size=20))
  )
  ggsave(paste(output_dir, '/INDEL_ALL_', metrics[i], '.pdf', sep=''), device = 'pdf', width=11, height=8.5)
  
  graphics.off()
  print(ggplot(snp_pass) +
          geom_point(aes(x=X, y=snp_pass[,i], col=Sample), size=5) +
          labs(x="", y=metrics[i]) +
          ylim(0,1) +
          theme_bw() +
          theme(axis.text.x =element_text(size=18,face="bold"),
                axis.text.y =element_text(size=20),
                axis.title.y = element_text(size=22,face="bold"),
                legend.title=element_text(size=22,face="bold"),
                legend.text=element_text(size=20))
  )
  ggsave(paste(output_dir, '/SNP_PASS_', metrics[i], '.pdf', sep=''), device = 'pdf', width=11, height=8.5)
  
  graphics.off()
  print(ggplot(indel_pass) +
          geom_point(aes(x=X, y=indel_pass[,i], col=Sample), size=5) +
          labs(x="", y=metrics[i]) +
          ylim(0,1) +
          theme_bw() +
          theme(axis.text.x =element_text(size=18,face="bold"),
                axis.text.y =element_text(size=20),
                axis.title.y = element_text(size=22,face="bold"),
                legend.title=element_text(size=22,face="bold"),
                legend.text=element_text(size=20))
  )
  ggsave(paste(output_dir, '/INDEL_PASS_', metrics[i], '.pdf', sep=''), device = 'pdf', width=11, height=8.5)
}

