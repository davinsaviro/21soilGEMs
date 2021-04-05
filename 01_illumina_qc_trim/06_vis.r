library(tidyverse)
library(ggplot2)
library(patchwork)

rawdata_readLength <- read_csv("C:/Users/Davin/Desktop/steps/003_raw_filter_clean/vis scripts/002004/20181119_1648_002-004_raw_readLength.txt",  col_names = "read_length") %>% arrange(., read_length)

number_reads <- nrow(rawdata_readLength)
mean_read <- round(mean(rawdata_readLength$read_length),digits = 0)
longest_read <- max(rawdata_readLength$read_length)

print(paste0("The sample has ",number_reads, " reads. The mean read length is ",mean_read, " and the longest read is ",longest_read))

rawdata_readLength$Cumsum_forward <- cumsum(rawdata_readLength$read_length)
rawdata_readLength$Cumsum_reverse <- max(rawdata_readLength$Cumsum_forward) - cumsum(rawdata_readLength$read_length)
Half_sequencing_output <- max(rawdata_readLength$Cumsum_forward)/2

print(paste0("The total number of sequenced bases is ",max(rawdata_readLength$Cumsum_forward)," and half of that is ",Half_sequencing_output))

histoPlot_01 <- ggplot(rawdata_readLength,aes(x=read_length))+
    geom_histogram(binwidth = 1000)+theme_classic()+
    geom_vline(xintercept = 7000,color="blue")+
    lims(x=c(0,longest_read))

print(paste0("The total number of sequenced bases is ",max(rawdata_readLength$Cumsum_forward)))

genome_size=12302133
expected_genome_coverage=100
bases_needed <- genome_size*expected_genome_coverage

cumSumPlot_01 <- ggplot(rawdata_readLength,aes(x=read_length,y=Cumsum_reverse))+
    geom_point()+theme_classic()+
    geom_vline(xintercept = 7000,color="blue")+
    geom_hline(yintercept = bases_needed,color="red")+
    labs(y="Reverse Cummulative read sum")+
    lims(x=c(0,longest_read))+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank())

cumSumPlot_01+histoPlot_01+plot_layout(nrow=2,heights = c(1,2))

##FILTERED STEPS

filtered_readLength <- read_csv("C:/Users/Davin/Desktop/steps/003_raw_filter_clean/vis scripts/002004/20181119_1648_002-004_filtered_readLength.txt",  col_names = "read_length") %>% arrange(., read_length) %>% add_column(reads="filtered reads")
rawdata_readLength <- rawdata_readLength %>% add_column(reads="raw reads")
	
filtered_readLength$Cumsum_forward <- cumsum(filtered_readLength$read_length)
filtered_readLength$Cumsum_reverse <- max(filtered_readLength$Cumsum_forward) - cumsum(filtered_readLength$read_length)

print(paste0("The total number of sequenced bases is ",max(filtered_readLength$Cumsum_forward)))

reads_total <- rbind(rawdata_readLength,filtered_readLength)

histoPlot_01 <- ggplot(reads_total,aes(x=read_length,color=reads,fill=reads))+
    geom_histogram(binwidth = 1000,alpha=0.5, position="identity")+theme_classic()+
    geom_vline(xintercept = 7000,color="blue")+
    lims(x=c(0,longest_read))
	
cumSumPlot_01 <- ggplot(reads_total,aes(x=read_length,y=Cumsum_reverse,color=reads,fill=reads))+
    geom_point()+theme_classic()+
    geom_vline(xintercept = 7000,color="blue")+
    #geom_hline(yintercept = Half_sequencing_output,color="red")+
    labs(y="Reverse Cummulative read sum")+
    lims(x=c(0,longest_read))+
    theme(axis.title.x = element_blank(),axis.text.x = element_blank())

cumSumPlot_01+histoPlot_01+plot_layout(nrow=2,heights = c(1,2))
