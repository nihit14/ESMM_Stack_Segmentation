library(ComplexHeatmap)
library(ggplot2)
library(reshape2)

setwd("/jmsh/projects/researchers/bins14/AG_Israel/")

column_names <- read.table("/jmsh/external/nihit/Israeli_methylation_dataset/matrices/files_list_order.txt")
df_column_names <- sapply(column_names$V1, function(x) gsub("^[^_]+_|-Z000.*$", "", x))
df_column_names[grep("CNVS-NORM-cfDNA", df_column_names)] <- "CNVS-NORM-cfDNA"
df_column_names[grep("CNVS-NORM-WBC-", df_column_names)] <- "CNVS-NORM-WBC"



samples <- as.data.frame(table(df_column_names))
colnames(samples) <- c("cell_type","counts")

samples$tissue_type <- sapply(as.character(samples$cell_type), function(x) strsplit(x,"-",fixed = TRUE)[[1]][1])

tissues <- unique(samples$tissue_type)

colors_37 <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#969696", "#9ECAE1", "#A1D99B", "#BC80BD", "#FDB462",
  "#80B1D3", "#FB8072", "#B3DE69", "#FCCDE5", "#D9D9D9",
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FDC086", "#CCEBC5",
  "#FFED6F", "#E41A1C"
)

colors.hex <- setNames(colors_37[seq_along(tissues)], tissues)




p <- ggplot(samples, aes(x=counts,y = reorder(cell_type, counts),fill = tissue_type)) +
  geom_bar(stat="identity") +theme_minimal() + 
  scale_fill_manual(values =  c(colors.hex)) +  
  geom_text(aes(label = counts),hjust = -0.1,size = 3)  + 
  expand_limits(x = max(samples$counts) * 1.1) +
  xlab("Counts") + ylab("Samples") 
p <- p + theme(axis.text = element_text(size = 14,color = "black")) # changes axis labels

p <- p + theme(axis.title = element_text(size = 14,color = "black")) # change axis titles

p <- p + theme(text = element_text(size = 14,color = "black"))  
length(unique(samples$tissue_type))
