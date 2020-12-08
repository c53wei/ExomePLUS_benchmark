args<-commandArgs(T)

## Rscript depth_bins_count_tables.R <INPUT_VCF> <OUTPUT_DIR>

# args[1]<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_depth_bins/test/HG002_NA24385_Son_100_exome_T.vcf.gz"
# args[2]<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_depth_bins/test"

output_folder<-args[2]


library(stringr)

vcf<-read.table(args[1], header=F, sep='\t', fill=T)
vcf<-vcf[c(1,2,4,5,8, 10,11)]
colnames(vcf)<-c("CHROM", "POS", "REF", "ALT", "INFO", "TRUTH", "QUERY")

vcf$TYPE_TRUTH<-str_split_fixed(vcf$TRUTH, ":", 7)[,5]
vcf$TRUTH<-str_split_fixed(vcf$TRUTH, ":", 7)[,2]
vcf$TYPE_QUERY<-str_split_fixed(vcf$QUERY, ":", 7)[,5]
vcf$QUERY<-str_split_fixed(vcf$QUERY, ":", 7)[,2]
vcf$DP_BIN<-as.character(str_split_fixed(vcf$INFO, "DP_BIN=", 2)[,2])

vcf<-vcf[nchar(vcf$DP_BIN)>0,]

# Total count by coverage bins
coverage_tables<-as.data.frame(table(vcf$DP_BIN))


tables_query<-table(vcf$TYPE_QUERY, vcf$QUERY, vcf$DP_BIN, exclude = c(".", " ", "", "NOCALL"))
tables_query<-as.data.frame(tables_query)

snp_tables_query<-split(tables_query, tables_query$Var1 == "INDEL")[['FALSE']]
indel_tables_query<-split(tables_query, tables_query$Var1 == "INDEL")[['TRUE']]
snp_metric_tables_query<-split(snp_tables_query, f=snp_tables_query$Var2, drop=T)
indel_metric_tables_query<-split(indel_tables_query, f=indel_tables_query$Var2, drop=T)

tables_truth<-table(vcf$TYPE_TRUTH, vcf$TRUTH, vcf$DP_BIN, exclude = c(".", " ", "", "NOCALL"))
tables_truth<-as.data.frame(tables_truth)


snp_tables_truth<-split(tables_truth, tables_truth$Var1 == "INDEL")[['FALSE']]
indel_tables_truth<-split(tables_truth, tables_truth$Var1 == "INDEL")[['TRUE']]
snp_metric_tables_truth<-split(snp_tables_truth, f=snp_tables_truth$Var2, drop=T)
indel_metric_tables_truth<-split(indel_tables_truth, f=indel_tables_truth$Var2, drop=T)

# SNP_tables
write.csv2(snp_metric_tables_query[["TP"]], file = paste(output_folder, '/', str_split_fixed(basename(args[1]), '.vcf', 3)[1], "_snp_query_TP.csv", sep = "" ))
write.csv2(snp_metric_tables_query[["FP"]], file = paste(output_folder, '/', str_split_fixed(basename(args[1]), '.vcf', 3)[1], "_snp_query_FP.csv", sep = "" ))
write.csv2(snp_metric_tables_query[["UNK"]], file = paste(output_folder, '/', str_split_fixed(basename(args[1]), '.vcf', 3)[1], "_snp_query_UNK.csv", sep = "" ))

write.csv2(snp_metric_tables_truth[["TP"]], file = paste(output_folder, '/', str_split_fixed(basename(args[1]), '.vcf', 3)[1], "_snp_truth_TP.csv", sep = "" ))
write.csv2(snp_metric_tables_truth[["FN"]], file = paste(output_folder, '/', str_split_fixed(basename(args[1]), '.vcf', 3)[1], "_snp_truth_FN.csv", sep = "" ))

# INDEL_tables

write.csv2(indel_metric_tables_query[["TP"]], file = paste(output_folder, '/', str_split_fixed(basename(args[1]), '.vcf', 3)[1], "_indel_query_TP.csv", sep = "" ))
write.csv2(indel_metric_tables_query[["FP"]], file = paste(output_folder, '/', str_split_fixed(basename(args[1]), '.vcf', 3)[1], "_indel_query_FP.csv", sep = "" ))
write.csv2(indel_metric_tables_query[["UNK"]], file = paste(output_folder, '/', str_split_fixed(basename(args[1]), '.vcf', 3)[1], "_indel_query_UNK.csv", sep = "" ))

write.csv2(indel_metric_tables_truth[["TP"]], file = paste(output_folder, '/', str_split_fixed(basename(args[1]), '.vcf', 3)[1], "_indel_truth_TP.csv", sep = "" ))
write.csv2(indel_metric_tables_truth[["FN"]], file = paste(output_folder, '/', str_split_fixed(basename(args[1]), '.vcf', 3)[1], "_indel_truth_FN.csv", sep = "" ))