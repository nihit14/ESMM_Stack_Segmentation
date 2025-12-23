library(ggplot2)
setwd("/jmsh/projects/researchers/bins14/AG_Israel/")

# read segmentaion bedfiles 
files <- c("results/Adipocytes-Z000000T7_MethylSeekR/Adipocytes-Z000000T7.MethylSeekR.segments.bed.gz",
           "results/Adipocytes_smooth_ESMM_4N_relabelled.bed.gz",
           "results/Adipocytes-Z000000T7_1bp_ESMM_4N_relabelled.bed.gz",
           "results/Adipocytes_200bp_smooth_ESMM_4N_relabelled.bed.gz")
                          

preprocess <- function(file_to_read,tag){
  df <- read.table(gzfile(file_to_read),header = FALSE)
  df <- df[,c(1:4)]
  df$length <- df$V3-df$V2
  df <- df[,c(4,5)]
  colnames(df) <- c("states","length")
  df$model <- tag
  return(df)
  
}

text_size=20
etext_size=22

# preparing dataframe 
mseekr <- preprocess(files[1],"mseekr")
mseekr$states[which(mseekr$states=="FMR")] <- "HMD"
esmm_smooth <- preprocess(files[2],"esmm_smooth")
esmm_raw <- preprocess(files[3],"esmm")
esmm_200bp_smooth <- preprocess(files[4],"esmm_smooth_200bp")
df <- rbind(mseekr,esmm_raw,esmm_smooth,esmm_200bp_smooth)

df$model <- factor(df$model,levels=c("mseekr","esmm_smooth","esmm_smooth_200bp","esmm"))


# State length of each model 
gplot <- ggplot(df, aes(x=states, y=length, fill=model)) + 
geom_boxplot() +  
  labs(x= "States", y = "Average segment length (bp)",fill="Model") + 
  theme(axis.text=element_text(size=etext_size), axis.title=element_text(size=text_size)) +  
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  ggtitle("States length") + xlim("HMD","PMD","LMR","UMR") + 
  theme_classic() + scale_fill_manual(values=c("esmm"="#FF7F0E",
                                               "esmm_smooth"="#3F7F93",
                                               "mseekr"="#2CA02C",
                                               "esmm_smooth_200bp"="#880808"))


gplot <- gplot + theme(axis.text = element_text(size = 20)) # changes axis labels

gplot <- gplot + theme(axis.title = element_text(size = 20)) # change axis titles

gplot <- gplot + theme(text = element_text(size = 20)) #+ guides(fill = "none")   # remove guide false to see legends

ggsave(filename = "Analysis/Adipocytes_avg_length.pdf",plot =  gplot,
       width = 8,height = 6)
ggsave(filename = "Analysis/Adipocytes_avg_length.png",plot =  gplot,
       dpi = 300,
       width = 8,height = 6)
# Segment percentage of hg38 comparison

counts <- read.table("/walter/projects/researchers/bins14/IHEC/counts/counts.bed",header = FALSE)
counts_length <- sum(counts$V3-counts$V2)
#options(scipen = 0)
options(scipen = 999)

heatmap <- aggregate(. ~ states + model, data = df, sum, na.rm = TRUE) 

heatmap$total <- (heatmap$length/counts_length)*100
      
hplot <- ggplot(heatmap, aes(x=states, y=total, fill=model)) + 
  geom_col(width=.5, position = "dodge") +  
  labs(x= "States", y = "Segment coverage %",fill="Model") + 
  theme(axis.text=element_text(size=etext_size), axis.title=element_text(size=text_size)) +  
  ggtitle("Segment coverage") + xlim("HMD","PMD","LMR","UMR") + 
  theme_classic() + scale_fill_manual(values=c("esmm"="#FF7F0E",
                                               "esmm_smooth"="#3F7F93",
                                               "mseekr"="#2CA02C",
                                               "esmm_smooth_200bp"="#880808"))

hplot <- hplot + theme(axis.text = element_text(size = 20)) # changes axis labels

hplot <- hplot + theme(axis.title = element_text(size = 20)) # change axis titles

hplot <- hplot + theme(text = element_text(size = 20)) #+ guides(fill = "none")   # remove guide false to see legends

ggsave(filename = "Analysis/Adipocytes_coverage.pdf",plot =  hplot,
       width = 8,height = 6)
ggsave(filename = "Analysis/Adipocytes_coverage.png",plot =  hplot,
       dpi = 300,
       width = 8,height = 6)

# Average methylation 
# read segmentaion bedfiles 
files <- c("Analysis/Adipocytes-Z000000T7_MethylSeekR.avg.methylation.tab",
           "Analysis/Adipocytes-Z000000T7_Adipocytes-Z000000T7_1bp_ESMM.avg.methylation.tab",
           "Analysis/Adipocytes-Z000000T7_Adipocytes-Z000000T7_1000bp_smooth_ESMM.avg.methylation.tab",
           "Analysis/Adipocytes-Z000000T7_Adipocytes-Z000000T7_200bp_smooth_ESMM.avg.methylation.tab")

preprocess_meth <- function(file_to_read,tag){
  df <- read.table(file_to_read,header = FALSE)
  colnames(df) <- c("states","meth")
  df$model <- tag
  return(df)
}

# preparing dataframe 
mseekr <- preprocess_meth(files[1],"mseekr")
mseekr$states[which(mseekr$states=="FMR")] <- "HMD"
esmm_raw <- preprocess_meth(files[2],"esmm")
esmm_smooth <- preprocess_meth(files[3],"esmm_smooth")
esmm_200bp_smooth <- preprocess_meth(files[4],"esmm_smooth_200bp")

df <- rbind(mseekr,esmm_raw,esmm_smooth,esmm_200bp_smooth)

df$model <- factor(df$model,levels=c("mseekr","esmm_smooth","esmm_smooth_200bp","esmm"))

iplot <- ggplot(df, aes(x=states, y=meth, fill=model)) + 
  geom_boxplot(width=.75) +  
  labs(x= "States", y = "Average methylation",fill="Model") + 
  theme(axis.text=element_text(size=etext_size), axis.title=element_text(size=text_size)) +  
  ggtitle("States methylation percentage") + xlim("HMD","PMD","LMR","UMR") + 
  theme_minimal() + scale_fill_manual(values=c("esmm"="#FF7F0E",
                                               "esmm_smooth"="#3F7F93",
                                               "mseekr"="#2CA02C",
                                               "esmm_smooth_200bp"="#880808"))

iplot <- iplot + theme(axis.text = element_text(size = 20)) # changes axis labels

iplot <- iplot + theme(axis.title = element_text(size = 20)) # change axis titles

iplot <- iplot + theme(text = element_text(size = 20)) #+ guides(fill = "none")   # remove guide false to see legends

ggsave(filename = "Analysis/Adipocytes_avg_meth.pdf",plot =  iplot,
       width = 8,height = 6)
ggsave(filename = "Analysis/Adipocytes_avg_meth.png",plot =  iplot,
       dpi = 300,
       width = 8,height = 6)




files <- c("/jmsh/projects/researchers/bins14/AG_Israel/Analysis/Adipocytes-Z000000T7_MethylSeekR.CpG_counts.tab",
           "/jmsh/projects/researchers/bins14/AG_Israel/Analysis/Adipocytes-Z000000T7_Adipocytes-Z000000T7_1bp_ESMM.CpG_counts.tab",
           "/jmsh/projects/researchers/bins14/AG_Israel/Analysis/Adipocytes-Z000000T7_Adipocytes-Z000000T7_1000bp_smooth_ESMM.CpG_counts.tab",
           "/jmsh/projects/researchers/bins14/AG_Israel/Analysis/Adipocytes-Z000000T7_Adipocytes-Z000000T7_200bp_smooth_ESMM.CpG_counts.tab")
 

preprocess_CpG <- function(file_to_read,tag){
  df <- read.table(file_to_read,header = FALSE)
  colnames(df) <- c("states","CpG","CpG_density")
  df$model <- tag
  return(df)
}

mseekr <- preprocess_CpG(files[1],"mseekr")
mseekr$states[which(mseekr$states=="FMR")] <- "HMD"
esmm_raw <- preprocess_CpG(files[2],"esmm")
esmm_smooth <- preprocess_CpG(files[3],"esmm_smooth")
esmm_200bp_smooth <- preprocess_CpG(files[4],"esmm_smooth_200bp")


df <- rbind(mseekr,esmm_raw,esmm_smooth,esmm_200bp_smooth)

df$model <- factor(df$model,levels=c("mseekr","esmm_smooth","esmm_smooth_200bp","esmm"))

jplot <- ggplot(df, aes(x=states, y=CpG, fill=model)) + 
  geom_boxplot(width=.75) +  
  labs(x= "States", y = "# CpGs",fill="Model") + 
  theme(axis.text=element_text(size=etext_size), axis.title=element_text(size=text_size)) +  
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  ggtitle("States CpG counts") + xlim("HMD","PMD","LMR","UMR") + 
  theme_minimal() + scale_fill_manual(values=c("esmm"="#FF7F0E",
                                               "esmm_smooth"="#3F7F93",
                                               "mseekr"="#2CA02C",
                                               "esmm_smooth_200bp"="#880808"))


kplot <- ggplot(df, aes(x=states, y=log(CpG_density), fill=model)) + 
  geom_boxplot(width=.75) +  
  labs(x= "States", y = "log(CpG Density)",fill="Model") + 
  theme(axis.text=element_text(size=etext_size), axis.title=element_text(size=text_size)) +  
  ggtitle("States CpG Density") + xlim("HMD","PMD","LMR","UMR") + 
  theme_minimal() + scale_fill_manual(values=c("esmm"="#FF7F0E",
                                               "esmm_smooth"="#3F7F93",
                                               "mseekr"="#2CA02C",
                                               "esmm_smooth_200bp"="#880808")) 

kplot <- kplot + theme(axis.text = element_text(size = 20)) # changes axis labels

kplot <- kplot + theme(axis.title = element_text(size = 20)) # change axis titles

kplot <- kplot + theme(text = element_text(size = 20)) #+ guides(fill = "none")   # remove guide false to see legends

ggsave(filename = "Analysis/Adipocytes_CpG_density.pdf",plot =  kplot,
       width = 8,height = 6)
ggsave(filename = "Analysis/Adipocytes_CpG_density.png",plot =  kplot,
       dpi = 300,
       width = 8,height = 6)
