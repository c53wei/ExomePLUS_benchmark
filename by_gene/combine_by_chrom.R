prep_data<-function(to_combine)
{
  combined<-do.call(rbind, to_combine)
  colnames(combined)[1]<-"gene"
  
  by_gene<-split.data.frame(combined, combined$gene)
  result<-vector(mode="list", length=length(by_gene))
  
  for(i in 1:length(by_gene))
  {
    
    new_variant_count<-sum(by_gene[[i]]$variant_count)
    new_average_depth<-0
    
    if(new_variant_count > 0)
    {
      new_average_depth<-sum(by_gene[[i]]$average_depth*by_gene[[i]]$variant_count)/new_variant_count
    }
    new_depth_as_char<-paste(by_gene[[i]]$depth_as_char[by_gene[[i]]$depth_as_char!="0"], collapse=" ")
    result[[i]]<-as.character(c(by_gene[[i]][1,][,c(1:6)], new_depth_as_char, new_average_depth, new_variant_count))
  }
  result<-as.data.frame(do.call(rbind, result))
  names(result)<-c("gene", "chrom", "start", "end", "region_length", "gc", "depth_as_char", "average_depth", "variant_count")
  result<-result[order(as.numeric(result[,2]), as.numeric(result[,3])), ]
  
  if(length(result[nchar(result$depth_as_char)<1,]$depth_as_char) > 0)
    result[nchar(result$depth_as_char)<1,]$depth_as_char<-'0'
  
  result$gc_content<-as.character(as.numeric(result$gc)/as.numeric(result$region_length))
  return(result)
}

args=commandArgs(T)

## Rscript combine_by_chrom.R <INPUT_DIR> <OUTPUT_DIR>
# args[1]<-"~/Documents/BME/Co-Op/McGill/R/ExomePLUS/CDS/bam_vcfs/9"
# args[2]<-"~/Documents/BME/Co-Op/McGill/R/ExomePLUS/CDS/bam_vcfs/output"

input_dir<-args[1]
output_dir<-args[2]

files<-list.files(args[1], ".csv", full.names = T)

to_combine<-vector(mode="list", length=length(files))

for(i in 1:length(files))
  to_combine[[i]]<-read.csv(files[i])

df<-prep_data(to_combine)

write.csv(df, paste(output_dir, '/', strsplit(basename(input_dir), '\\.')[[1]][1], ".csv", sep=''), row.names = F)
