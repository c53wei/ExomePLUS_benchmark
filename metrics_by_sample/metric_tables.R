make_table<-function(input_dir, output_dir, X_axis)
{
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
  
  for(i in 1:length(metrics))
  {
    metrics[[i]]<-do.call(rbind, metrics[[i]])
    metrics[[i]]$Sample<-substr(basename(filenames), 1, 5)
    metrics[[i]]$X_axis<-X_axis
    rownames(metrics[[i]])<-basename(filenames)
    write.csv(metrics[[i]], paste(output_dir, '/', names(metrics)[i], '.csv', sep = '') )
  }
}

input_dir<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/WGS_all"
output_dir<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/ExomePLUS_benchmark/metrics_by_sample/metric_data/illumina_WGS"

X_axis<-c(41.7, 2.6, 6.6, 1.6, 1.6, 4.4, 4.4,
          30, 2.0, 6.4, 1.7, 1.7, 6.2, 6.2,
          30, 3.2, 5.5, 1.9, 2.0, 4.8, 4.8)

# input_dir<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_sample/illumina/illumina_round2"
# output_dir<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/ExomePLUS_benchmark/metrics_by_sample/metric_data/illumina_WES"
#   
# X_axis<-c("4plex_1", "4plex_A", "4plex_B", "8plex_A", "8plex_A", "8plex_B", "8plex_B",
#            "4plex_1", "4plex_A", "4plex_B", "8plex_A", "8plex_A", "8plex_B", "8plex_B",
#            "4plex_1", "4plex_A", "4plex_B", "8plex_A", "8plex_A", "8plex_B", "8plex_B", 
#            "4plex_1", "4plex_A", "4plex_B", "8plex_A", "8plex_A", "8plex_B", "8plex_B")

make_table(input_dir, output_dir, X_axis)


