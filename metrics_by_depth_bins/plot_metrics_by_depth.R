read.files<-function(input_folder, intervals)
{
  filenames <- list.files(input_folder, pattern = ".csv", full.names = T)
  
  indel_snp<-grepl('indel', filenames)
  
  filenames<-filenames[indel_snp]
  
  shortnames <- list.files(input_folder, pattern = ".csv", full.names = F)[indel_snp]
  
  filelist  <- vector(mode="list", length=length(filenames))
  
  metric_tables<-vector(mode="list", length=length(filenames))
  names(metric_tables)<-shortnames
  
  
  for(i in 1:length(metric_tables))
  {
    metric_tables[[shortnames[i]]]<-vector(mode='list', length=length(intervals))
    names(metric_tables[[shortnames[i]]])<-intervals
    
    for(j in 1:length(metric_tables[[shortnames[[i]]]]))
      metric_tables[[shortnames[[i]]]][[j]]<-0
  }
  
  for(i in 1:length(filenames))
  {
    
    file<-tryCatch({
      if (file.size(filenames[i]) > 0){
        read.csv2(filenames[i])
      }
    }, error = function(err) {
      # error handler picks up where error was generated
      print(paste("Read.table didn't work!:  ",err))
    })
    
    tryCatch({
      file$Var3<-gsub(',', '', file$Var3)
      file$Var3<-gsub(' ', ', ', file$Var3)
      for(j in 1:nrow(file))
      {
        interval<-file$Var3[j]
        freq<-file$Freq[j]
        
        metric_tables[[shortnames[[i]]]][[interval]]<-metric_tables[[shortnames[[i]]]][[interval]]+freq
      }
    }, error = function(err) {
      # error handler picks up where error was generated
      print(paste("indexing frequencies error—no table",err))
    })
    
    metric_tables[[shortnames[i]]]<-as.data.frame((do.call(rbind, metric_tables[[shortnames[i]]])))
    rownames(metric_tables[[shortnames[i]]])<-intervals
    colnames(metric_tables[[shortnames[i]]])<-"Freq"
  }
  
  return(metric_tables)
}

library(ggplot2)
library(stringr)


exome<-T

if(exome)
{
  intervals<-c("(0.0, 10.0]", "(10.0, 20.0]", "(20.0, 40.0]", "(40.0, 60.0]", "(60.0, 80.0]", "(80.0, 100.0]", 
               "(100.0, 200.0]", ">200.0")
  
  IDs<-c("HG002_WGS_30X", "HG002_XTV2_4plex", "HG002_XTV2_1plex", "HG002_XTV2_4plex", "HG002_XT", 
         "HG003_WGS_30X", "HG003_XTV2_1plex", "HG003_XTV2_4plex", "HG003_XT",
         "HG004_WGS_30X", "HG004_XTV2_1plex", "HG004_XTV2_4plex", "HG004_XT")
  
  IDs_order<-c("HG002_XT", "HG003_XT", "HG004_XT",
              "HG002_XTV2_1plex", "HG003_XTV2_1plex", "HG004_XTV2_1plex",
              "HG002_XTV2_4plex", "HG003_XTV2_4plex", "HG004_XTV2_4plex",
              "HG002_WGS_30X", "HG003_WGS_30X", "HG004_WGS_30X")
  
  custom<-c(brewer.pal(3, 'Reds'), brewer.pal(3, 'Greens'), brewer.pal(3, 'Blues'), brewer.pal(3, 'Greys') )
  
  type<-"WES"
  
  input_folder<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_depth_bins/WES_count_tables_PASS/"
  output_folder<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_depth_bins/WES_figures_PASS/"
} else {
  
  intervals<-c("(0.0, 1.0]", "(1.0, 5.0]", "(5.0, 10.0]", "(10.0, 20.0]", "(20.0, 30.0]", ">30.0")
  
  IDs<-c("HG002_30X", "HG002_2X", "HG002_5X",
              "HG003_30X", "HG003_2X", "HG003_5X",
              "HG004_30X", "HG004_2X", "HG004_5X")
  
  IDs_order<-c("HG002_2X", "HG003_2X", "HG004_2X",
               "HG002_5X", "HG003_5X", "HG004_5X",
               "HG002_30X", "HG003_30X", "HG004_30X")
  
  custom<-c(brewer.pal(3, 'Reds'), brewer.pal(3, 'Greens'), brewer.pal(3, 'Blues') )
  
  type<-"WGS"
  
  input_folder<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_depth_bins/WGS_count_tables_PASS/"
  output_folder<-"/Users/clarewei/Documents/BME/Co-Op/McGill/ExomePLUS/metrics_by_depth_bins/WGS_figures_PASS/"
}

filelist<-read.files(input_folder, intervals)

count_tables<-vector(mode="list", length=length(filelist)/5)
metric_tables<-vector(mode="list", length=length(filelist)/5)
names(metric_tables)<-IDs
samples<-names(metric_tables)
metrics<-c("Recall", "Precision", "Frac_NA", "F1_Score")

for( i in 1:length(count_tables) )
{
  count_tables[[i]]<-as.data.frame(cbind(filelist[[5*i-4]]$Freq, filelist[[5*i-3]]$Freq, filelist[[5*i-2]]$Freq, filelist[[5*i-1]]$Freq, filelist[[5*i-0]]$Freq ))
  colnames(count_tables[[i]])<-c("FP_Query", "TP_Query", "UNK_Query", "FN_Truth", "TP_Truth")
  
  rownames(count_tables[[i]])<-intervals
  
  metric_tables[[samples[i]]]$Recall <- count_tables[[i]]$TP_Truth/(count_tables[[i]]$TP_Truth+count_tables[[i]]$FN_Truth)
  metric_tables[[samples[i]]]$Precision <- count_tables[[i]]$TP_Query/(count_tables[[i]]$TP_Query+count_tables[[i]]$FP_Query)
  metric_tables[[samples[i]]]$Frac_NA <- count_tables[[i]]$UNK_Query/(count_tables[[i]]$TP_Query+count_tables[[i]]$FP_Query+count_tables[[i]]$UNK_Query)
  metric_tables[[samples[i]]]$F1_Score <- 2*metric_tables[[samples[i]]]$Recall*metric_tables[[samples[i]]]$Precision/(metric_tables[[samples[i]]]$Recall+metric_tables[[samples[i]]]$Precision)
  
  metric_tables[[samples[i]]]<-as.data.frame(cbind(metric_tables[[samples[i]]]$Recall, metric_tables[[samples[i]]]$Precision, metric_tables[[samples[i]]]$Frac_NA, metric_tables[[samples[i]]]$F1_Score))
  colnames(metric_tables[[samples[i]]])<-metrics
  rownames(metric_tables[[samples[i]]])<-intervals
  
  metric_tables[[samples[i]]]$Index <- seq(from=1, to=nrow(metric_tables[[samples[i]]]), by = 1)
  metric_tables[[samples[i]]]$Legend <- samples[i]
}

to_plot<-do.call("rbind", metric_tables)
to_plot$Legend<-factor(to_plot$Legend, levels=IDs_order)

for( i in 1:length(metrics) )
{
  graphics.off()
  print(ggplot(to_plot, aes(x=Index, y=to_plot[[i]], col=Legend)) +
          geom_point(size=7) + 
          scale_colour_manual(values=custom) +
          labs(x="Depth Bins", y = metrics[i]) +
          theme(plot.title = element_text(hjust=0.5)) +
          theme_bw() +
          theme(axis.text.x =element_text(size=18, angle=20, hjust=1, face='bold'),
                axis.text.y =element_text(size=18),
                axis.title.x = element_text(size=22,face="bold"),
                axis.title.y = element_text(size=22,face="bold"),
                legend.position = 'none') +
          ylim(0,1) +
    
          scale_x_discrete(limits=intervals)
  )
  ggsave(paste(output_folder, metrics[i],"_indel.pdf",sep=""), device = "pdf", width=11, height=8.5)
}

  



  
  

