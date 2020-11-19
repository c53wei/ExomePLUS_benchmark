input_dir<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/WGS_all"
output_dir<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/ExomePLUS_benchmark/metrics_by_sample/metric_data/illumina_WGS"
  
filenames<-list.files(input_dir, full.names = T, pattern='.summary.csv')
files<-vector(mode="list", length=length(filenames))

metrics<-vector(mode='list', length=4)
names(metrics)<-c('snp_all', 'indel_all', 'snp_pass', 'indel_pass')

for(i in 1:length(metrics))
{
  metrics[[i]]<-vector(mode="list", length=length(filenames))
}

for( i in 1:length(files))
{
  files[[i]]<-read.csv(filenames[[i]])
  temp<-split(files[[i]], paste(files[[i]]$Type, files[[i]]$Filter))
  metrics[['snp_all']][[i]]<-temp[["SNP ALL"]][c(11:14)]
  metrics[['indel_all']][[i]]<-temp[["INDEL ALL"]][c(11:14)]
  metrics[['snp_pass']][[i]]<-temp[["SNP PASS"]][c(11:14)]
  metrics[['indel_pass']][[i]]<-temp[["INDEL PASS"]][c(11:14)]
}
mean_coverage<-c(41.7, 2.6, 6.6, 3.1, 8.8, 
                 30, 2.0, 6.4, 3.4, 12.4,
                 30, 3.2, 5.5, 3.9, 9.6)
median_coverage<-c(42.0, 2.0, 6.6, 3.0, 8.0,
                   30, 2.0, 6.0, 3.0, 12.0,
                   30, 3.0, 5.0, 4.0, 9.0)

for(i in 1:length(metrics))
{
  metrics[[i]]<-do.call(rbind, metrics[[i]])
  metrics[[i]]$mean_coverage<-mean_coverage
  metrics[[i]]$median_coverage<-median_coverage
  metrics[[i]]$Sample<-substr(basename(filenames), 1, 5)
  rownames(metrics[[i]])<-basename(filenames)
  write.csv(metrics[[i]], paste(output_dir, '/', names(metrics)[i], '.csv', sep = '') )
}

