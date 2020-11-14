library(ggplot2)
input_dir<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/ExomePLUS_benchmark/metrics_by_sample/metric_data/illumina_WGS"
output_dir<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_sample/illumina/WGS_figures"

filelist<-list.files(input_dir, full.names=T, pattern='.csv')
metrics<-c('Recall', 'Precision', 'Frac_UNK', 'F1_Score')

for(i in 1:length(filelist))
{
  to_plot<-as.data.frame(read.csv(filelist[[i]], row.names = 1))
  for(j in 1:4 )
  {
    graphics.off()
    print(ggplot(to_plot) +
            geom_point(aes(x=mean_coverage, y=to_plot[,j], col=Sample), size=5) +
            labs(x="", y=metrics[j]) +
            ylim(0,1) +
            theme_bw() +
            theme(axis.text.x =element_text(size=18,face="bold"),
                  axis.text.y =element_text(size=20),
                  axis.title.y = element_text(size=22,face="bold"),
                  legend.title=element_text(size=22,face="bold"),
                  legend.text=element_text(size=20))
    )
    ggsave(paste(output_dir, '/', strsplit(basename(filelist), '\\.')[[i]][1], '_', metrics[j], '.pdf', sep=''), device = 'pdf', width=11, height=8.5)
  }
  
}
