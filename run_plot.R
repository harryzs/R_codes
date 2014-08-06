##6.8.2014

run_plot_for_FASTQC_summary<-function(path_to_file, file){
  
  cat ("\n run_plot_for_FASTQC_summary start\n")
  
  setwd(path_to_file)
  
  data <- read.table(file = file, header = TRUE, sep = "\t",stringsAsFactors  = FALSE )
  
  ## --------------- define a function to get Lanes from sample_names
  get_Lanes <- function(string, pattern = "_") {
    
    index = length(unlist(strsplit(string, pattern, )))
    
    return(unlist(strsplit(string, pattern, ))[index-1])
  }

  ## --------------- get Lanes
  Lanes <- unlist(lapply(data$sample_name, get_Lanes))
  
  data$lanes <- Lanes
  
  ##--------load library ggplots2
  library(ggplot2)
  library("grid")  
  
  ##-------set parameters for ggplots
  ##plot width
  width = 10
  
  ##plot height
  height = 10
  
  ##themes
  my_themes = theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "none") +
    theme(axis.text = element_text( size = rel(2), face = "bold", colour = "blue")) +
    theme(axis.title = element_text( size = rel(2), face = "bold" )) +
    theme(plot.title = element_text( size = rel(3), face = "bold", colour = "red")) 
    
  ##--------using ggplots2
  p <- ggplot(data = data , aes(x = lanes, y = as.numeric(reads_number)/1e6, fill = lanes)) + 
    geom_bar(stat = "identity", colour = "black", size = 1) 

    # scale_y_continuous(limits=c(0, 400))
                       
  p <- p + labs(title = "Total number of reads in each lane \n", x = "\n Lanes", y = "Number of reads (Million) \n") +
    my_themes

  ##Save plot to a PDF file
  ggsave(filename = "Plot_FastQC.pdf", plot = p, width = width, height = height)

  
  cat ("\n run_plot_for_FASTQC_summary finished\n")
}



run_plot_for_EstimateLibraryComplexity <- function(path_to_file, file){
  
  cat ("\n run_plot_for_EstimateLibraryComplexity start\n")
  
  setwd(path_to_file)
  
  data <- read.table(file = file, header = TRUE, sep = "\t")
  ##--------load library ggplots2
  library(ggplot2)
  library("grid")  
  
  ##-------set parameters for ggplots
  ##plot width
  width = 10
  
  ##plot height
  height = 10
  
  ##themes
  my_themes = theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none") +
    theme(axis.text = element_text( size = rel(1.5), face = "bold", colour = "blue")) +
    theme(axis.title = element_text( size = rel(1.5), face = "bold" )) +
    theme(plot.title = element_text( size = rel(2), face = "bold", colour = "red")) 
  
  
  ##--------PERCENT_DUPLICATION
  data$sample_name <- factor(data$sample_name)
  p <- ggplot(data = data , aes(x = sample_name, y = PERCENT_DUPLICATION, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "The percentage of mapped reads\n that is marked as duplicate \n",
                x = "\nSamples", y = "Perecntage\n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_PERCENT_DUPLICATION.pdf", plot = p, width = width, height = height)
  
  
  ##--------ESTIMATED_LIBRARY_SIZE
  p <- ggplot(data = data , aes(x = sample_name, y = ESTIMATED_LIBRARY_SIZE/1e6, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "The estimated number of \n unique molecules in the library  \n", 
                x = "\nSamples", y = "Unique molecules (1e6)\n") + my_themes
  p
  ##Save plot to a PDF file
  ggsave(filename = "Plot_ESTIMATED_LIBRARY_SIZE.pdf", plot = p, width = width, height = height)
  
  cat ("\n run_plot_for_EstimateLibraryComplexity finished\n")
}


run_plot_for_STAR_summary<-function(path_to_file, file){
  cat ("\n run_plot_for_STAR_summary start\n")
  setwd(path_to_file)
  
  data <- read.table(file = file, header = TRUE, sep = "\t")
#   colnames(data)
#   [1] "sample_name"                               "Number.of.input.reads"                    
#   [3] "Average.input.read.length"                 "UNIQUE.READS."                            
#   [5] "Uniquely.mapped.reads.number"              "Uniquely.mapped.reads.."                  
#   [7] "Average.mapped.length"                     "Number.of.splices..Total"                 
#   [9] "Number.of.splices..Annotated..sjdb."       "Number.of.splices..GT.AG"                 
#   [11] "Number.of.splices..GC.AG"                  "Number.of.splices..AT.AC"                 
#   [13] "Number.of.splices..Non.canonical"          "Mismatch.rate.per.base..."                
#   [15] "Deletion.rate.per.base"                    "Deletion.average.length"                  
#   [17] "Insertion.rate.per.base"                   "Insertion.average.length"                 
#   [19] "MULTI.MAPPING.READS."                      "Number.of.reads.mapped.to.multiple.loci"  
#   [21] "X..of.reads.mapped.to.multiple.loci"       "Number.of.reads.mapped.to.too.many.loci"  
#   [23] "X..of.reads.mapped.to.too.many.loci"       "UNMAPPED.READS."                          
#   [25] "X..of.reads.unmapped..too.many.mismatches" "X..of.reads.unmapped..too.short"          
#   [27] "X..of.reads.unmapped..other"   
  
  head(data)
  data$total_mapping_rate <- 100-data[,23]-data[,25]-data[,26]-data[,27]
  data$PEC_of_annotated_splices <- (data[,9]/data[,8])*100
  ##--------load library ggplots2
  library(ggplot2)
  library("grid")  
  ##-------set parameters for ggplots
  ##plot width
  width = 12
  
  ##plot height
  height = 10
  
  ##themes
  my_themes = theme(axis.text.x = element_text(angle = 0, hjust = 1)) + theme(legend.position = "none") +
              theme(axis.text = element_text( size = rel(1.5), face = "bold", colour = "blue")) +
              theme(axis.title = element_text( size = rel(1.5), face = "bold" )) +
              theme(plot.title = element_text( size = rel(2), face = "bold", colour = "red")) 
      
  
  ##--------Number.of.input.reads
  p <- ggplot(data = data , aes(x = sample_name, y = (Number.of.input.reads)/1e6 , fill = sample_name)) + 
  geom_bar(stat = "identity",  colour = "black", size = 1) 
  
  p <- p + labs(title = "Number of raw reads for each sample \n ", x = "Samples", y = "Number of reads (Million)\n") + my_themes
  p <- p + coord_flip()

  ##Save plot to a PDF file
  ggsave(filename = "Plot_STAR_input_reads.pdf", plot = p, width = width, height = height)
  

  #----------PEC of annotated splices
  p <- ggplot(data = data , aes(x = sample_name, y = PEC_of_annotated_splices, fill = sample_name)) + 
  geom_bar(stat = "identity",  colour = "black", size = 1) 

  p <- p + labs(title = "Percent of annotated splieces \n ", x = "Samples \n", y = "Percent \n") + my_themes
  p <- p + coord_flip()
  ##Save plot to a PDF file
  ggsave(filename = "Plot_STAR_PEC_of_annotated_splices.pdf", plot = p, width = width, height = height)


  ##--------total mapping rates
  p <- ggplot(data = data , aes(x = sample_name, y = total_mapping_rate, fill = sample_name)) + 
  geom_bar(stat = "identity",  colour = "black", size = 1) 

  p <- p + labs(title = "Total mapping rates \n ", x = "Samples", y = "Mapping rates (%)\n") + my_themes
  p <- p + coord_flip()
  ##Save plot to a PDF file
  ggsave(filename = "Plot_STAR_total_mapping_rates.pdf", plot = p, width = width, height = height)

  ##--------total unmapping rates
  p <- ggplot(data = data , aes(x = sample_name, y = 100-total_mapping_rate, fill = sample_name)) + 
   geom_bar(stat = "identity",  colour = "black", size = 1) 

  p <- p + labs(title = "Total unmapping rates \n ", x = "Samples", y = "Unmapping rates (%)\n") + my_themes
  p <- p + coord_flip()
  ##Save plot to a PDF file
  ggsave(filename = "Plot_STAR_total_unmapping_rates.pdf", plot = p, width = width, height = height)

  
  ##--------Uniquely.mapped.reads.number
  p <- ggplot(data = data , aes(x = sample_name, y = Uniquely.mapped.reads.number/1e6 , fill = sample_name)) + 
  geom_bar(stat = "identity",  colour = "black", size = 1) 
  
  p <- p + labs(title = "Uniquely.mapped.reads.number \n", x = "Samples", y = "Number of reads (Million)\n") +  my_themes
  p <- p + coord_flip()
  ##Save plot to a PDF file
  ggsave( filename = "Plot_STAR_Uniquely_mapped_reads.pdf", plot = p,  width = width, height = height)

 
  ##--------percent of Uniquely.mapped.reads.number
  p <- ggplot(data = data , aes(x = sample_name, y = Uniquely.mapped.reads.., fill = sample_name)) + 
  geom_bar(stat = "identity",  colour = "black", size = 1) 
  
  p <- p + labs(title = "Percentage_of_uniquely.mapped.reads \n ", x = "Samples", y = "Percentage (%)\n") + 
    my_themes
  p <- p + coord_flip()
  ##Save plot to a PDF file
  ggsave( filename = "Plot_STAR_percent_of_Uniquely_mapped_reads.pdf", plot = p, width = width, height = height)

  
  ##--------Average mapped length
  p <-ggplot(data = data , aes(x = sample_name, y = Average.mapped.length , fill = sample_name)) + 
  geom_bar(stat = "identity",  colour = "black", size = 1) 
  
  p <- p + labs(title = "Average.mapped.length \n", x = "Samples", y = "Length (bp)\n") + 
    my_themes
  p <- p + coord_flip()
  ##Save plot to a PDF file
  ggsave(filename = "Plot_STAR_Average_mapped_length.pdf", plot = p, width = width, height = height)
  
  ##--------Number.of.splices..Total
  p <- ggplot(data = data , aes(x = sample_name, y = Number.of.splices..Total/1e6, fill = sample_name)) + 
  geom_bar(stat="identity",  colour = "black", size = 1)

  p <- p + labs(title = "Total splices number \n", x = "Samples", y = "Number (Million)\n") + 
  my_themes
  p <- p + coord_flip()
  ##Save plot to a PDF file
  ggsave(filename = "Plot_STAR_Number_of_Total_splices.pdf", plot = p, width = width, height = height)

  
  ##--------Mismatch.rate.per.base...
  p <- ggplot(data = data , aes(x = sample_name, y = Mismatch.rate.per.base..., fill = sample_name)) + 
  geom_bar(stat="identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "Mismatch.rate.per.base \n", x = "Samples", y = "Mismatch Rate (%)\n") + 
    my_themes
  p <- p + coord_flip()
  ##Save plot to a PDF file
  ggsave(filename = "Plot_STAR_Mismatch_rate_per_base.pdf", plot = p, width = width, height = height)
  

  ##--------X..of.reads.mapped.to.multiple.loci
  p <- ggplot(data = data , aes(x = sample_name, y = X..of.reads.mapped.to.multiple.loci, fill = sample_name)) + 
  geom_bar(stat="identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "Perecentage of Reads.mapped.to.multiple.loci \n", x = "Samples", y = "Perecentage of Reads (%)\n") + 
    my_themes 
  p <- p + coord_flip()
  ##Save plot to a PDF file
  ggsave(filename ="Plot_STAR_reads_mapped_to_multiple_loci.pdf", plot = p, width = width, height = height )

  ##--------X..of.reads.mapped.to.too.many.loci
  p <- ggplot(data = data , aes(x = sample_name, y = X..of.reads.mapped.to.too.many.loci, fill = sample_name)) + 
  geom_bar(stat="identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "Perecentage of Reads.mapped.to.too.many.loci \n", x = "Samples", y = "Perecentage of Reads (%)\n") + 
    my_themes  
  p <- p + coord_flip()
  ##Save plot to a PDF file
  ggsave( filename = "Plot_STAR_reads_mapped_to_too_many_loci.pdf", plot = p, width = width, height = height)
    


  ##-------Plot_STAR_Uniquely.mapped.reads_V.S._Number.of.Total.Splices
  
  p <- ggplot(data = data, aes(x = Uniquely.mapped.reads.number/1e6, y= Number.of.splices..Total/1e6)) + 
    geom_smooth(method = "lm")  + geom_point()
  p <- p + 
    labs(title = "Uniquely.mapped.reads V.S. Number.of.Total.Splices \n", 
         x = "Uniquely.mapped.reads (Million)\n", 
         y = "Number.of.Total.Splices (Million)\n") +
    theme(axis.text = element_text( size = rel(1.5), face = "bold", colour = "blue")) +
    theme(axis.title = element_text( size = rel(1.5), face = "bold" )) +
    theme(plot.title = element_text( size = rel(2), face = "bold", colour = "red")) +
    theme(legend.text = element_text(size = rel(1.5), face = "bold")) +
    theme(legend.title = element_text(size = rel(1.5), face = "bold"))
    
  x = max((data$Uniquely.mapped.reads.number)/1e6)-2
  y = max((data$Number.of.splices..Total)/1e6)-2
  p <-p + geom_text(data = data, 
  aes(label = paste("r = ", round(cor(Uniquely.mapped.reads.number, Number.of.splices..Total), 2), sep=" ")),
  x = max((data$Uniquely.mapped.reads.number)/1e6)-2, y = max((data$Number.of.splices..Total)/1e6)-2, size = 10)

  
  ##Save plot to a PDF file
  ggsave( filename = "Plot_STAR_Uniquely_mapped_reads_VS_Number_of_Total_Splices.pdf", 
          plot = p, width = width, height = height)
    
  cat ("\n run_plot_for_STAR_summary finished\n")
 
  
}


run_plot_for_CollectRnaSeqMetrics_summary <-function(path_to_file, file){
  
  cat ("\n run_plot_for_CollectRnaSeqMetrics_summary start\n")
  
  setwd(path_to_file)
  
  data <- read.table(file = file, header = TRUE, sep = "\t")
  dim(data)
  colnames(data)
  ##--------load library ggplots2
  library(ggplot2)
  library("grid")  
  
  ##-------set parameters for ggplots
  ##plot width
  width = 10
  
  ##plot height
  height = 10
  
  ##themes
  my_themes = theme(axis.text.x= element_text(angle = 45, hjust = 1)) + theme(legend.position = "none") +
    theme(axis.text = element_text( size = rel(1.5), face = "bold", colour = "blue")) +
    theme(axis.title = element_text( size = rel(1.5), face = "bold" )) +
    theme(plot.title = element_text( size = rel(2), face = "bold", colour = "red")) 

  ##--------PCT_RIBOSOMAL_BASES
  p <- ggplot(data = data , aes(x = sample_name, y = PCT_RIBOSOMAL_BASES*1e5, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "Perecntage_of_Ribosomal_Reads \n", x = "\nSamples", y = "PCT_RIBOSOMAL_READS (1e-5)\n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_PCT_RIBOSOMAL_BASES.pdf", plot = p, width = width, height = height)
  
  ##----------PCT_USABLE_BASES
  p <- ggplot(data = data , aes(x = sample_name, y = PCT_USABLE_BASES*100, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "Perecntage_of_usabel_reads \n", x = "\nSamples", y = "PCT_USABL_Reads\n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_PCT_USABLE_BASES.pdf", plot = p, width = width, height = height)
    
  
  ##--------PCT_MRNA_BASES
  p <- ggplot(data = data , aes(x = sample_name, y = PCT_MRNA_BASES*100, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "Perecntage_of_mRNA_Reads \n", x = "\nSamples", y = "PCT_MRNA_Reads\n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_PCT_MRNA_BASES.pdf", plot = p, width = width, height = height)
  
  
  ##--------MEDIAN_5PRIME_TO_3PRIME_BIAS
  p <- ggplot(data = data , aes(x = sample_name, y = MEDIAN_5PRIME_TO_3PRIME_BIAS, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(
    title = "The ratio of coverage at the 5' end of to the 3' end \n based on the 1000 most highly expressed transcripts\n", 
    x = "\nSamples", y = "MEDIAN_5PRIME_TO_3PRIME_BIAS\n") + my_themes
  
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_MEDIAN_5PRIME_TO_3PRIME_BIAS.pdf", plot = p, width = width, height = height)
  
  
  ##--------PCT_CODING_BASES
  p <- ggplot(data = data , aes(x = sample_name, y = PCT_CODING_BASES*100, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "PCT_CODING_Reads \n", x = "\nSamples", y = "Percentage\n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_PCT_CODING_BASES.pdf", plot = p, width = width, height = height)
  
  
  ##--------PCT_UTR_BASES
  p <- ggplot(data = data , aes(x = sample_name, y = PCT_UTR_BASES*100, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "PCT_UTR_Reads \n", x = "\nSamples", y = "Percentage\n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_PCT_UTR_BASES.pdf", plot = p, width = width, height = height)
  
  
  
  
  ##--------PCT_INTRONIC_BASES
  p <- ggplot(data = data , aes(x = sample_name, y = PCT_INTRONIC_BASES*100, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "PCT_INTRONIC_Reads\n", x = "\nSamples", y = "Percentage\n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_PCT_INTRONIC_BASES.pdf", plot = p, width = width, height = height)
  
  
  ##--------PCT_INTERGENIC_BASES
  p <- ggplot(data = data , aes(x = sample_name, y = PCT_INTERGENIC_BASES*100, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "PCT_INTERGENIC_Reads\n", x = "\nSamples", y = "Percentage\n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_PCT_INTERGENIC_BASES.pdf", plot = p, width = width, height = height)
  
  
  ##--------MEDIAN_CV_COVERAGE
  p <- ggplot(data = data , aes(x = sample_name, y = MEDIAN_CV_COVERAGE, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "Median CV of coverage of the 1000 \n most highly expressed transcripts\n", 
                x = "\nSamples", y = "median CV\n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_MEDIAN_CV_COVERAGE.pdf", plot = p, width = width, height = height)
  
  
  ##--------MEDIAN_5PRIME_BIAS
  p <- ggplot(data = data , aes(x = sample_name, y = MEDIAN_5PRIME_BIAS, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "The median 5 prime bias of the 1000 \n most highly expressed transcripts\n", 
                x = "\nSamples", y = "5 prime bias\n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_MEDIAN_5PRIME_BIAS.pdf", plot = p, width = width, height = height)
  
  
  ##--------MEDIAN_3PRIME_BIAS
  p <- ggplot(data = data , aes(x = sample_name, y = MEDIAN_3PRIME_BIAS, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "The median 3 prime bias of the \n1000 most highly expressed transcripts\n", 
                x = "\nSamples", y = "3 prime bias\n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_MEDIAN_3PRIME_BIAS.pdf", plot = p, width = width, height = height)
  
    
  cat ("\n run_plot_for_CollectRnaSeqMetrics_summary finished\n")
    
}


run_plot_for_Infer_experiment_summary <- function(path_to_file, file){
  cat ("\n run_plot_for_Infer_experiment_summary start\n")
  
  setwd(path_to_file)
  
  data <- read.table(file = file, header = TRUE, sep = "\t")
  dim(data)
  colnames(data)
  ##--------load library ggplots2
  library(ggplot2)
  library("grid")  
  
  ##-------set parameters for ggplots
  ##plot width
  width = 10
  
  ##plot height
  height = 10
  
  ##themes
  my_themes = theme(axis.text.x = element_text(angle = 0, hjust = 1)) + theme(legend.position = "none") +
    theme(axis.text = element_text( size = rel(1.5), face = "bold", colour = "blue")) +
    theme(axis.title = element_text( size = rel(1.5), face = "bold" )) +
    theme(plot.title = element_text( size = rel(2), face = "bold", colour = "red")) +
    theme(legend.text = element_text( size = rel(1.5), face = "bold")) +
    theme(legend.title = element_text( size = rel(1.5), face = "bold"))
  
  ##--------Infer_experiment
  p <- ggplot(data = data, aes(x = sample_names, y = type_one_reads/type_two_reads, fill = sample_names)) +
  geom_bar(stat = "identity",  colour = "black", size = 1)

  p <- p + labs(title = "Infer_experiment \n", x = "Samples", y = "Ratio of different types of reads \n") + my_themes
  p <- p + coord_flip()
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_Infer_experiment.pdf", plot = p, width = width, height = height)
  
  cat ("\n run_plot_for_Infer_experiment_summary finished\n")
}


run_plot_for_preseq <- function(path_to_file, sample_name, c_curve = TRUE, lc_extrap = TRUE){
  
  cat ("\n run_plot_for_preseq  start\n")
  
  setwd(path_to_file)
  
  file_to_c_curve <- paste(sample_name,"_c_curve.txt", sep = "")
  file_to_c_curve 
  
  file_to_lc_extrap <- paste(sample_name,"_lc_extrap.txt", sep ="")
  file_to_lc_extrap
  
  ##--------load library ggplots2
  library(ggplot2)
  library("grid")  
  
  
  ##-------set parameters for ggplots
  ##plot width
  width = 10
  
  ##plot height
  height = 10
  
  ##themes
  my_themes = theme(axis.text = element_text( size = rel(1.5), face = "bold", colour = "blue")) +
    theme(axis.title = element_text( size = rel(1.5), face = "bold" )) +
    theme(plot.title = element_text( size = rel(2), face = "bold", colour = "red")) +
    theme(legend.text = element_text( size = rel(1.5), face = "bold")) +
    theme(legend.title = element_text( size = rel(1.5), face = "bold"))
  
  
  
  ##---------------------------------------------------------------------------------------
  ##----------------------------for lc_extrap
  if(lc_extrap == TRUE){
  data <- read.table(file = file_to_lc_extrap, header = TRUE, sep = "\t")

  ##---------------------------------------------------------------------------------------
  ##---------------for all data in lc_extrap 
  p <- ggplot(data = data, aes(x = TOTAL_READS/1e6, y = EXPECTED_DISTINCT/1e6)) +
    geom_line()
  
  p <- p + labs(title = paste("lc_extrap for ",sample_name, "\n", sep = ""), 
                x = "\n TOTAL_READS (Mio)", y = "EXPECTED_DISTINCT (Mio) \n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = paste(sample_name,"_Plot_lc_extrap.pdf", sep = ""), plot = p, width = width, height = height)
  
  
  ##---------------------------------------------------------------------------------------
  ##---------------for first xxx lines in lc_extrap
  data_sub <- data[1:500,]
  p <- ggplot(data = data_sub, aes(x = TOTAL_READS/1e6, y = EXPECTED_DISTINCT/1e6)) +
    geom_line()
  
  p <- p + labs(title = paste("lc_extrap for ",sample_name, "\n", sep = ""), 
                x = "\n TOTAL_READS (Mio) ", y = "EXPECTED_DISTINCT (Mio) \n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = paste(sample_name,"_Plot_lc_extrap_first_","500", ".pdf", sep = ""), 
         plot = p, width = width, height = height)
  
  }
  
 
  ##---------------------------------------------------------------------------------------
  ##---------------for c_curve
  if(c_curve == TRUE){  
    
  data <- read.table(file = file_to_c_curve, header = TRUE, sep = "\t")
  #dim(data)
  #head(data)
  #colnames(data)

  ##--------Infer_experiment
  p <- ggplot(data = data, aes(x = total_reads/1e6, y = distinct_reads/1e6)) +
    geom_line()
  
  p <- p + labs(title = paste( "c_curve for ",sample_name, "\n", sep = ""), 
                x = "\n total_reads (Mio) ", y = "distinct_reads (Mio) \n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = paste( sample_name,"_Plot_c_curv.pdf",sep = ""), plot = p, width = width, height = height)
  
  }
  
  cat ("\n run_plot_for_preseq finished \n")
    
}


run_plot_for_Btrim_summary <- function(path_to_file, file){
  cat ("\n run_plot_for_Btrim_summary start\n")
  
  setwd(path_to_file)
  
  data <- read.table(file = file, header = TRUE, sep = "\t")
  dim(data)
  colnames(data)
  ##--------load library ggplots2
  library(ggplot2)
  library("grid")  
  
  ##-------set parameters for ggplots
  ##plot width
  width = 10
  
  ##plot height
  height = 10
  
  ##themes
  my_themes = theme(axis.text.x = element_text(angle = 0, hjust = 1)) + theme(legend.position = "none") +
    theme(axis.text = element_text( size = rel(0.8), face = "bold", colour = "blue")) +
    theme(axis.title = element_text( size = rel(1.5), face = "bold" )) +
    theme(plot.title = element_text( size = rel(2), face = "bold", colour = "red")) +
    theme(legend.text = element_text( size = rel(1.5), face = "bold")) +
    theme(legend.title = element_text( size = rel(1.5), face = "bold"))
  
  ##--------PEC_of_Pattern_trimmed
  p <- ggplot(data = data, aes(x = sample_names, y = PEC_of_Pattern_trimmed *100, fill = sample_names)) +
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "Percentage of reads with adapters \n", x = "Samples", y = "Percentage \n") + my_themes
  p <- p + coord_flip()
  ##Save plot to a PDF file
  ggsave(filename = "Plot_Btrim_PEC_of_Pattern_trimmed.pdf", plot = p, width = width, height = height)
  
  
  ##--------PEC_of_Qual_trimmed
  p <- ggplot(data = data, aes(x = sample_names, y = PEC_of_Qual_trimmed, fill = sample_names)) +
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "Percentage of Qual_trimmed_reads  \n", x = "Samples", y = "Percentage \n") + my_themes
  p <- p + coord_flip()
  ##Save plot to a PDF file
  ggsave(filename = "Plot_Btrim_PEC_of_Qual_trimmed.pdf", plot = p, width = width, height = height)
  
  cat ("\n run_plot_for_Btrim_summary\n")
  
}


run_plot_for_HTSeq <- function(path_to_file, file){
  cat ("\n run_plot_for_HTSeq start\n")
  
  setwd(path_to_file)
  
  data <- read.table(file = file, header = TRUE, sep = "\t")
  dim(data)
  colnames(data)
  ##--------load library ggplots2
  library(ggplot2)
  library("grid")  
  
  ##-------set parameters for ggplots
  ##plot width
  width = 10
  
  ##plot height
  height = 10
  
  ##themes
  my_themes = theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.position = "none") +
    theme(axis.text = element_text( size = rel(1.5), face = "bold", colour = "blue")) +
    theme(axis.title = element_text( size = rel(1.5), face = "bold" )) +
    theme(plot.title = element_text( size = rel(2), face = "bold", colour = "red")) +
    theme(legend.text = element_text( size = rel(1.5), face = "bold")) +
    theme(legend.title = element_text( size = rel(1.5), face = "bold"))
  
  ##--------percentage_of_reads_mapped_to_features
  p <- ggplot(data = data, aes(x = sample_name, y = percentage_of_reads_mapped_to_features*100, fill = sample_name)) +
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "Reads_mapped_to_features \n", x = "Samples", y = "Percentage \n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_HTSeq_percentage_of_reads_mapped_to_features.pdf", plot = p, width = width, height = height)
  
  
  
  ##--------reads_mapped_to_features
  p <- ggplot(data = data, aes(x = sample_name, y = reads_mapped_to_features/1e6, fill = sample_name)) +
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "Reads_mapped_to_features \n", x = "Samples", y = "Raw reads number (Million) \n") + my_themes
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_HTSeq_reads_mapped_to_features.pdf", plot = p, width = width, height = height)
  
  cat ("\n run_plot_for_HTSeq finished\n")
}


run_plot_for_alignment_summary_metrics <-function(path_to_file, file){
  
  cat ("\n run_plot_for_alignment_summary_metrics start\n")
  
  setwd(path_to_file)
  
  data <- read.table(file = file, header = TRUE, sep = "\t")
  dim(data)
  colnames(data)
  #   [1] "sample_name"                "CATEGORY"                   "TOTAL_READS"                "PF_READS"                  
  #   [5] "PCT_PF_READS"               "PF_NOISE_READS"             "PF_READS_ALIGNED"           "PCT_PF_READS_ALIGNED"      
  #   [9] "PF_ALIGNED_BASES"           "PF_HQ_ALIGNED_READS"        "PF_HQ_ALIGNED_BASES"        "PF_HQ_ALIGNED_Q20_BASES"   
  #   [13] "PF_HQ_MEDIAN_MISMATCHES"    "PF_MISMATCH_RATE"           "PF_HQ_ERROR_RATE"           "PF_INDEL_RATE"             
  #   [17] "MEAN_READ_LENGTH"           "READS_ALIGNED_IN_PAIRS"     "PCT_READS_ALIGNED_IN_PAIRS" "BAD_CYCLES"                
  #   [21] "STRAND_BALANCE"             "PCT_CHIMERAS"               "PCT_ADAPTER"                "SAMPLE"                    
  #   [25] "LIBRARY"                    "READ_GROUP"     
  
  
  data$sample_name <- paste(data$sample_name,data$CATEGORY, sep ="_")
  head(data)
  ##--------load library ggplots2
  library(ggplot2)
  library("grid")  
  
  ##-------set parameters for ggplots
  ##plot width
  width = 10
  
  ##plot height
  height = 10
  
  ##themes
  my_themes = theme(axis.text.x = element_text(angle = 0, hjust = 1)) + theme(legend.position = "none") +
    theme(axis.text = element_text( size = rel(1), face = "bold", colour = "blue")) +
    theme(axis.title = element_text( size = rel(1.5), face = "bold" )) +
    theme(plot.title = element_text( size = rel(2), face = "bold", colour = "red")) 
  
  ##--------TOTAL_READS
  p <- ggplot(data = data , aes(x = sample_name, y = TOTAL_READS/1e6, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "TOTAL_READS \n", x = "Samples\n", y = "\n raw number (Million)") + my_themes + coord_flip()
  
  ##Save plot to a PDF file
  ggsave(filename = "TOTAL_READS.pdf", plot = p, width = width, height = height)
  
  
  ##--------PF_READS_ALIGNED
  p <- ggplot(data = data , aes(x = sample_name, y = PF_READS_ALIGNED/1e6, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "READS_ALIGNED \n", x = "Samples\n", y = "\n raw number (Million)") + my_themes + coord_flip()
  
  ##Save plot to a PDF file
  ggsave(filename = "PF_READS_ALIGNED.pdf", plot = p, width = width, height = height)
  
  
  ##--------PCT_PF_READS_ALIGNED
  p <- ggplot(data = data , aes(x = sample_name, y = PCT_PF_READS_ALIGNED*100, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "PCT_PF_READS_ALIGNED \n", x = "Samples\n", y = "\nPercentage") + my_themes + coord_flip()
  
  ##Save plot to a PDF file
  ggsave(filename = "PCT_PF_READS_ALIGNED.pdf", plot = p, width = width, height = height)
  
  ##--------PCT_READS_ALIGNED_IN_PAIRS
  p <- ggplot(data = data , aes(x = sample_name, y = PCT_READS_ALIGNED_IN_PAIRS*100, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "PCT_READS_ALIGNED_IN_PAIRS \n", x = "Samples\n", y = "\nPercentage") + my_themes + coord_flip()
  
  ##Save plot to a PDF file
  ggsave(filename = "PCT_READS_ALIGNED_IN_PAIRS.pdf", plot = p, width = width, height = height)
  
  ##--------STRAND_BALANCE
  p <- ggplot(data = data , aes(x = sample_name, y = STRAND_BALANCE, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "Reads mapped to positive strand/total mapped reads\n", x = "Samples\n", y = "\nSTRAND_BALANCE") + my_themes
  p <- p + coord_flip()
  ##Save plot to a PDF file
  ggsave(filename = "STRAND_BALANCE.pdf", plot = p, width = width, height = height)
  
  ##--------PCT of High quality aligned PF reads (high quality == mapping quality >= 20)
  p <- ggplot(data = data , aes(x = sample_name, y = (PF_HQ_ALIGNED_READS/TOTAL_READS)*100, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1)
  
  p <- p + labs(title = "High quality aligned PF reads \n (high_quality:mapping quality >= 20)\n", 
                x = "Samples\n", y = "\nPercent") + my_themes
  p <- p + coord_flip()
  
  ##Save plot to a PDF file
  ggsave(filename = "PCT_HQ_ALIGNED_READS.pdf", plot = p, width = width, height = height)
  
  cat ("\n run_plot_for_alignment_summary_metrics finished\n")
  
}


run_plot_for_DuplicationMetrics <-function(path_to_file, file){
  
  cat ("\n run_plot_for_DuplicationMetrics start\n")
  
  setwd(path_to_file)
  
  data <- read.table(file = file, header = TRUE, sep = "\t")
  dim(data)
  colnames(data)
  head(data)
  
  #   [1] "sample_name"                  "LIBRARY"                      "UNPAIRED_READS_EXAMINED"     
  #   [4] "READ_PAIRS_EXAMINED"          "UNMAPPED_READS"               "UNPAIRED_READ_DUPLICATES"    
  #   [7] "READ_PAIR_DUPLICATES"         "READ_PAIR_OPTICAL_DUPLICATES" "PERCENT_DUPLICATION"         
  #   [10] "ESTIMATED_LIBRARY_SIZE"  
  
  ##--------load library ggplots2
  library(ggplot2)
  library("grid")  
  
  ##-------set parameters for ggplots
  ##plot width
  width = 10
  
  ##plot height
  height = 10
  
  ##themes
  my_themes = theme(axis.text.x = element_text(angle = 0, hjust = 1)) + theme(legend.position = "none") +
    theme(axis.text = element_text( size = rel(1.5), face = "bold", colour = "blue")) +
    theme(axis.title = element_text( size = rel(1.5), face = "bold" )) +
    theme(plot.title = element_text( size = rel(2), face = "bold", colour = "red")) 
  
  ##--------PERCENT_DUPLICATION
  p <- ggplot(data = data , aes(x = sample_name, y = PERCENT_DUPLICATION*100, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1) 
  
  p <- p + labs(title = "Duplication Rate\n", x = "Samples \n", y = "\nPercentage") + my_themes + coord_flip()
  
  ##Save plot to a PDF file
  ggsave(filename = "PERCENT_DUPLICATION.pdf", plot = p, width = width, height = height)
  
  
  ##--------ESTIMATED_LIBRARY_SIZE
  p <- ggplot(data = data , aes(x = sample_name, y = ESTIMATED_LIBRARY_SIZE/1e6, fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1) 
  
  p <- p + labs(title = "ESTIMATED_LIBRARY_SIZE\n", x = "Samples \n", y = "\n Library Size (Million)") + my_themes + coord_flip()
  
  ##Save plot to a PDF file
  ggsave(filename = "ESTIMATED_LIBRARY_SIZE.pdf", plot = p, width = width, height = height)
  
  
  
  cat ("\n run_plot_for_DuplicationMetrics finished\n")
  
}


run_plot_for_abnormal_reads<-function(path_to_file, file){
  
  cat ("\n run_plot_for_abnormal_reads start\n")
  
  setwd(path_to_file)
  
  data <- read.table(file = file, header = TRUE, sep = "\t")
  colnames(data)
  #   [1] "sample_name"                  "Total_STAR_mapped_reads"      "Reads_mapeed_to_normal"      
  #   [4] "Reads_mapeed_to_abnormal"     "Percentage_of_normal_reads"   "Percentage_of_abnormal_reads"
  
  head(data)
  
  ##--------load library ggplots2
  library(ggplot2)
  library("grid")  
  ##-------set parameters for ggplots
  ##plot width
  width = 12
  
  ##plot height
  height = 10
  
  ##themes
  my_themes = theme(axis.text.x = element_text(angle = 0, hjust = 1)) + theme(legend.position = "none") +
    theme(axis.text = element_text( size = rel(1.5), face = "bold", colour = "blue")) +
    theme(axis.title = element_text( size = rel(1.5), face = "bold" )) +
    theme(plot.title = element_text( size = rel(2), face = "bold", colour = "red")) 
  
  
  ##--------Number.of.input.reads
  p <- ggplot(data = data , aes(x = sample_name, y = Percentage_of_abnormal_reads*100 , fill = sample_name)) + 
    geom_bar(stat = "identity",  colour = "black", size = 1) 
  
  p <- p + labs(title = "Percentage of reads mapped to scaffolds \n ", x = "Samples \n", y = " \n Percentage") + my_themes
  p <- p + coord_flip()
  
  ##Save plot to a PDF file
  ggsave(filename = "Plot_STAR_Percentage_of_reads_mapped_to_scaffolds.pdf", plot = p, width = width, height = height)
  
  cat ("\n run_plot_for_abnormal_reads finished\n")
  
}


