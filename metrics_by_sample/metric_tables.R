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
    metrics[['snp_all']][[i]]<-cbind(temp[["SNP ALL"]][c(11:14, 4:5)], temp[["SNP ALL"]][6]-temp[["SNP ALL"]][7]-temp[["SNP ALL"]][8], temp[["SNP ALL"]]$QUERY.FP)
    metrics[['indel_all']][[i]]<-cbind(temp[["INDEL ALL"]][c(11:14, 4:5)], temp[["INDEL ALL"]][6]-temp[["INDEL ALL"]][7]-temp[["INDEL ALL"]][8], temp[["INDEL ALL"]]$QUERY.FP)
    metrics[['snp_pass']][[i]]<-cbind(temp[["SNP PASS"]][c(11:14, 4:5)], temp[["SNP PASS"]][6]-temp[["SNP PASS"]][7]-temp[["SNP PASS"]][8], temp[["SNP PASS"]]$QUERY.FP)
    metrics[['indel_pass']][[i]]<-cbind(temp[["SNP ALL"]][c(11:14, 4:5)], temp[["INDEL PASS"]][6]-temp[["INDEL PASS"]][7]-temp[["INDEL PASS"]][8], temp[["INDEL PASS"]]$QUERY.FP)
  }
  
  for(i in 1:length(metrics))
  {
    metrics[[i]]<-do.call(rbind, metrics[[i]])
    metrics[[i]]$Sample<-substr(basename(filenames), 1, 5)
    metrics[[i]]$X_axis<-X_axis
    rownames(metrics[[i]])<-basename(filenames)
    colnames(metrics[[i]])[c(7,8)]<-c("QUERY.TP", "QUERY.FP")
    write.csv(metrics[[i]], paste(output_dir, '/', names(metrics)[i], '.csv', sep = '') )
  }
}


## -----------------------------------------------------------------------------------------------------------------


# input_dir<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/WGS_all"
# output_dir<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/ExomePLUS_benchmark/metrics_by_sample/metric_data/illumina_WGS"
# 
# X_axis<-c(41.7, 2.6, 6.6, 1.6, 1.6, 4.4, 4.4,
#           30, 2.0, 6.4, 1.7, 1.7, 6.2, 6.2,
#           30, 3.2, 5.5, 1.9, 2.0, 4.8, 4.8)

input_dir<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_sample/mgi/happy_WES"
output_dir<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/ExomePLUS_benchmark/metrics_by_sample/metric_data/mgi_WES"

X_axis<-c( "4plex_A", "4plex_B", "8plex_A", "8plex_B", "8plex_B",
           "4plex_A", "4plex_B", "8plex_A", "8plex_A", "8plex_B", "8plex_B", 
           "4plex_A", "4plex_B", "8plex_A", "8plex_A", "8plex_B", "8plex_B", 
           "4plex_A", "4plex_B", "8plex_A", "8plex_A", "8plex_B", "8plex_B")

# X_axis<-c("4plex_1", "4plex_A", "4plex_B", "8plex_A", "8plex_A", "8plex_B", "8plex_B",
#            "4plex_1", "4plex_A", "4plex_B", "8plex_A", "8plex_A", "8plex_B", "8plex_B", "XT",
#            "4plex_1", "4plex_A", "4plex_B", "8plex_A", "8plex_A", "8plex_B", "8plex_B", "XT",
#            "4plex_1", "4plex_A", "4plex_B", "8plex_A", "8plex_A", "8plex_B", "8plex_B", "XT")

make_table(input_dir, output_dir, X_axis)


