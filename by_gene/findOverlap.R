numbers_only <- function(x) !grepl("\\D", x)

bed_to_granges <- function(vcf){
  df <- read.table(vcf,
                   header=F,
                   stringsAsFactors=F)
  
  
  df$GC<-df[,ncol(df)-2] + df[,ncol(df)-3]
  
  df<-df[,c(1:3, ncol(df)-6, ncol(df))]
  
  names(df)<-c('chr','start','end','score','gc')
  
  df<-df[numbers_only(df$chr),]
  
  gr<-with(df, GRanges(chr, IRanges(start, end), score=score, gc=gc))

  return(gr)
}

vcf_to_granges<-function(vcf)
{
  df <- read.table(vcf,
                   header=F,
                   stringsAsFactors=F)
  names(df)<-c("chr", "start", "DP")
  df<-df[numbers_only(df$chr),]
  
  gr <- with(df, GRanges(chr, IRanges(start, start), score=DP))
  
  return(gr)
}

args=commandArgs(T)

## Rscript findOverlap.R <VCF> <BED> <OUTPUTDIR>

# args[1]="~/Documents/BME/Co-Op/McGill/R/ExomePLUS/CDS/bam_vcfs/04.txt"
# args[2]="~/Documents/BME/Co-Op/McGill/R/ExomePLUS/CDS/bam_vcfs/1.bed"
# args[3]="~/Documents/BME/Co-Op/McGill/R/ExomePLUS/CDS/bam_vcfs/output"


library(BiocManager)
library(VariantAnnotation)
library(GenomicFeatures)

vcf<-vcf_to_granges(args[1])
bed<-bed_to_granges(args[2])
outputDir<-args[3]

bed_by_genes<-split(bed, score(bed))
region_length<-vector(mode="integer", length=length(bed_by_genes))
chrom<-vector(mode="integer", length=length(bed_by_genes))
start<-vector(mode="integer", length=length(bed_by_genes))
end<-vector(mode="integer", length=length(bed_by_genes))
gc<-vector(mode="integer", length=length(bed_by_genes))
gene<-vector(mode="integer", length=length(bed_by_genes))

# computing length of each gene
for(i in 1:length(bed_by_genes))
{
  gene[i]<-score(bed_by_genes[[i]][1])
  region_length[i]<-sum(width(bed_by_genes[[i]]))
  gc[i]<-sum(mcols(bed_by_genes[[i]])$gc)
  chrom[i]<-seqlevelsInUse(bed_by_genes[[i]])
  start[i]<-start(bed_by_genes[[i]])[1]
  end[i]<-end(bed_by_genes[[i]])[length(bed_by_genes[[i]])]
}

gene_info<-as.data.frame(cbind(gene,chrom, start, end, region_length, gc))
desired_order<-unique(score(bed))
gene_info$gene <- factor( as.character(gene_info$gene), levels=desired_order )
gene_info <- gene_info[order(gene_info$gene),]

overlaps <- findOverlaps(vcf, bed)

# Features frooverlaps vcf with overlaps in bed
# Note: The same feature frooverlaps vcf can overlap with mulitple features frooverlaps bed
vcf_matched <- vcf[queryHits(overlaps)]

# Add the metadata frooverlaps bed
mcols(vcf_matched) <- cbind.data.frame(
  mcols(vcf_matched),
  mcols(bed[subjectHits(overlaps)]));


depth<-vector(mode="list", length=length(unique(score(bed))))
names(depth)<-unique(score(bed))

unique_genes<-sort(unique(mcols(vcf_matched)[2])@listData[["score"]])

vcf_matched<-split(vcf_matched, mcols(vcf_matched)[2])


for(i in 1:length(vcf_matched))
{
  index<-unique_genes[i]
  
  if(length(strsplit(index, " ")[[1]])>1)
  {
    for(j in index)
    {
      depth[[j]]<-append(depth[[j]], score(vcf_matched[[i]]))
    }
  }
  else{
    depth[[index]]<-score(vcf_matched[[i]])
  }
}

average_depth<-unlist(lapply(depth, mean))
average_depth[is.na(average_depth)] <- 0

variant_count<-unlist(lapply(depth, length))

depth_as_char<-sapply(depth, paste, collapse=' ')

             
df<-as.data.frame(cbind(gene_info, depth_as_char, average_depth, variant_count))
df$depth_as_char[unlist(lapply(depth, is.null))]<-0
rownames(df)<-unique(score(bed))
                     
write.csv(df, paste(outputDir, '/', strsplit(basename(args[1]), '\\.')[[1]][1], '_', strsplit(basename(args[2]), '\\.')[[1]][1], ".csv", sep=''), row.names = F)

  