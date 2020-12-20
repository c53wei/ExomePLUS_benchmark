median_depth<-function(to_plot)
{
  for( index in 1:length(to_plot))
  {
    to_plot[[index]]$depth_as_char<-lapply(sapply(to_plot[[index]]$depth_as_char, strsplit, " "), as.numeric)
    to_plot[[index]]$median_depth<-as.numeric(lapply(to_plot[[index]]$depth_as_char, median))
  }
  return(to_plot)
}

my_boxplot<-function(to_plot, median)
{
  stat_type<-"Average"
  
  if(median)
  {
    to_plot<-median_depth(to_plot)
    stat_type<-"Median"
  }
  
  boxplot<-do.call(rbind, to_plot)
  
  
  boxplot<-boxplot[boxplot$variant_count>0,]
  boxplot$ID <- factor(boxplot$ID, levels=IDs_order)
  boxplot$Sample<-substr(boxplot$ID, 1, 5)
  stat_to_plot<-boxplot$average_depth
  
  if(median)
    stat_to_plot<-boxplot$median_depth
    
  
  graphics.off()
  result<-print(ggplot(boxplot, aes(x=ID, y=stat_to_plot, col=Sample)) +
          geom_boxplot() +
          geom_hline(yintercept=100) +
          labs(title = paste(boxplot$type[1],  "Coverage by Gene", sep=' '), x = " ", y = paste(stat_type, " Depth"), sep='') +
          theme_bw() +
          theme(plot.title = element_text(size=16,face="bold", hjust=0.5),
                axis.text.x =element_text(angle = 20, hjust = 1, face='bold'),
                axis.text.y =element_text(size=14),
                axis.title.y = element_text(size=14,face="bold"),
                axis.title.x  = element_text(size=14,face="bold"),
                legend.position='none')
  ) 
  ggsave(paste(outputDir, '/boxplot_', stat_type, ".pdf", sep=''), device = "pdf", width=11, height=8.5)
  
  return(to_plot)
}

coverage_by_gc<-function(to_plot, i)
{
  to_plot<-to_plot[to_plot$average_depth > 0 & to_plot$gc_content > 0,]
  
    
  graphics.off()
  print(ggplot(to_plot) +
          geom_point(aes(x=gc_content, y=average_depth), alpha=0.07) +
          labs(title = i, x = "GC Content", y = "Average Depth") +
          theme_bw() +
          theme(plot.title = element_text(size=14,face="bold", hjust=0.5),
                axis.title.y = element_text(size=12,face="bold"),
                axis.title.x  = element_text(size=12,face="bold")) +
          xlim(0,0.8) +
          ylim(0,200)
          
  )
  ggsave(paste(outputDir, '/', i, '_gc_content.pdf', sep='' ), device = "pdf", width=11, height=8.5 )
}

stats<-function(to_plot)
{
  median<-stats_summ(to_plot, T)
  average<-stats_summ(to_plot, F)
  summary<-cbind(average, median[,-1])
  summary$ID<-factor(summary$ID, levels=IDs_order)
  summary<-summary[order(summary$ID),]
  summary[,2:5]<-round(summary[,2:5], 2)
  
  write.csv(summary, paste(outputDir, "/stats_summary.csv", sep=''), row.names = F)
  
}

stats_summ<-function(to_plot, median)
{
  combined<-do.call(rbind, to_plot)
  
  if(median)
  {
    summ <-  combined %>%
      group_by(ID) %>%
      summarize(mean_median = mean(median_depth), median_median = median(median_depth))
  }
  else
  {
    summ <-  combined %>%
      group_by(ID) %>%
      summarize(mean_average = mean(average_depth), median_average = median(average_depth))
  }
    
  
  
  return(summ)
}
## --------------------------------------------------------------------------------------------------------

args=commandArgs(T)

## Rscript plot_by_gene <INPUT_DIR> <OUTPUT_DIR>

args[1]<-"/Users/clarewei/Documents/BME/Co-Op/McGill/R/ExomePLUS/CDS/bam_vcfs/UTR_combined_new"
args[2]<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/by_gene/UTR_figures"

inputDir<-args[1]
outputDir<-args[2]

library(ggpubr)
library(ggplot2)
library(dplyr)

IDs<-c("HG002_XTV2_4plex_New_Index", "HG002_XTV2_1plex", "HG002_XTV2_4plex", 
       "HG003_XTV2_1plex", "HG003_XTV2_4plex",
       "HG004_XTV2_1plex", "HG004_XTV2_4plex",
       "HG002_XT", "HG003_XT", "HG004_XT")

IDs_order<-c("HG002_XT", "HG003_XT", "HG004_XT",
             "HG002_XTV2_1plex", "HG003_XTV2_1plex", "HG004_XTV2_1plex",
             "HG002_XTV2_4plex", "HG002_XTV2_4plex_New_Index", "HG003_XTV2_4plex", "HG004_XTV2_4plex")

dirnames<-list.dirs(inputDir, full.names=T, recursive = F)
to_plot<-vector(mode="list", length=length(dirnames))
names(to_plot)<-IDs

for(d in 1:length(dirnames))
{
  index<-IDs[d]
  
  filenames<-list.files(dirnames[d], full.names = T)
  files<-vector(mode="list", length=length(filenames))
  
  for(i in 1:length(filenames))
  {
    files[[i]]<-read.csv(filenames[i])
    colnames(files[[i]])[1]<-"gene"
  }
  
  to_plot[[index]]<-do.call(rbind, files)
  to_plot[[index]]$ID<-IDs[d]
  to_plot[[index]]$type<-substr(basename(inputDir), 1, 3)
  
  to_plot[[index]][,c(2:6,8:10)]<-lapply(to_plot[[index]][,c(2:6,8:10)], as.numeric)
}

to_plot<-my_boxplot(to_plot, F);
to_plot<-my_boxplot(to_plot, T)

stats(to_plot)

# gc_plot<-vector(mode="list", length=10)
# names(gc_plot)<-IDs_order

for(i in 1:length(to_plot))
{
  index<-IDs[i]
  coverage_by_gc(to_plot[[index]], index)
}

# figure<-ggarrange(gc_plot[[1]], gc_plot[[2]], gc_plot[[3]], nrow=2, ncol=2)
# ggexport( figure, filename=paste(outputDir, "/V6UTR_1plex_gc_content.pdf", sep=''))
# 
# figure<-ggarrange(gc_plot[[4]], gc_plot[[5]], gc_plot[[6]], nrow=2, ncol=2)
# ggexport( figure, filename=paste(outputDir, "/XTV2_1plex_gc_content.pdf", sep=''))
# 
# figure<-ggarrange(gc_plot[[7]], gc_plot[[8]], gc_plot[[9]], gc_plot[[10]], nrow=2, ncol=2)
# ggexport( figure, filename=paste(outputDir, "/XTV2_4plex_gc_content.pdf", sep=''))
