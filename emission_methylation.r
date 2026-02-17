library(ComplexHeatmap)

setwd("/jmsh/projects/researchers/bins14/AG_Israel/results_v2/stack_model_BB_207_samples_50N/")

file_emission <- "stack_model_BB_207_samples_50N/emissions_model-stack_model_BB_207_samples_50N_converted_final-model-stack_model_BB_207_samples_50N_50_N.txt"
emission <- read.table(file_emission,
                       header = TRUE,sep = "\t")

# correcting two column  names GSM5652340 -> GSM5652350
# and GSM5652240 -> GSM5652250 # grep("2350$",colnames(emission)) & grep("2250$",colnames(emission))
colnames(emission)[166] <- "GSM5652340"
colnames(emission)[66] <- "GSM5652240"


column_names <- read.csv("/jmsh/external/nihit/Israeli_methylation_dataset/matrices/metadata.csv")
column_names2 <- read.csv("../../metadata/41586_2022_5580_MOESM4_ESM.csv",header = T,skip = 2)

column_names$V1 <- column_names$Sample.name
column_names2$V1 <- column_names2$Sample.name

column_names$V1 <- sapply(column_names$file_name, function(x) {
  sub(".*?(Z.*)", "\\1", x)
})

column_names2$V1 <- sapply(column_names2$V1, function(x) {
  sub(".*?(Z.*)", "\\1", x)
})


columns_merged <- merge(column_names, column_names2, by = "V1", all.x = TRUE,sort = FALSE)
columns_merged[is.na(columns_merged)] <- "N/A"

# emission matrix columns mapping 
column_maps <- read.table("ids_maps.txt",header = FALSE)
column_maps$id <- column_maps$V2

columns_merged2 <- merge(column_maps, columns_merged, by = "id", all.x = TRUE,sort = FALSE)
columns_merged2[is.na(columns_merged2)] <- "N/A"
columns_merged2$Cell.type <- factor(columns_merged2$Cell.type,
                                    levels=c(unique(columns_merged2$Cell.type))) 
columns_merged2$Source.Tissue <- factor(columns_merged2$Source.Tissue,
                                    levels=c(unique(columns_merged2$Source.Tissue)))

# reorder 
#

df_column_names <- sapply(df_column_names$V1, function(x) gsub("^[^_]+_|-Z000.*$", "", x))
df_tissue_names <-  sapply(df_column_names, function(x) strsplit(x,"-",fixed = TRUE)[[1]][1])
table(df_tissue_names)

df_column_names[grep("CNVS-NORM-cfDNA", df_column_names)] <- "CNVS-NORM-cfDNA"
df_column_names[grep("CNVS-NORM-WBC-", df_column_names)] <- "CNVS-NORM-WBC"

df_column_names
# [1] "Aorta-Endothel"
emission_mat <- as.matrix(emission[,c(2:ncol(emission))])
rownames(emission_mat) <- c(paste("State",emission$States))
dim(emission_mat)

col_ha <- HeatmapAnnotation(
  Sex = columns_merged2$sex,
  Tissues = columns_merged2$Source.Tissue,
  Cell_type = columns_merged2$Cell.type,
  annotation_name_gp = gpar(fontsize = 20),
  annotation_legend_param = list(
    Cell_type = list(
      ncol = 2,
      title_gp  = gpar(fontsize = 20, fontface = "bold"),
      labels_gp = gpar(fontsize = 20)
    ),
    Sex = list(
      ncol = 2,
      title_gp  = gpar(fontsize = 20, fontface = "bold"),
      labels_gp = gpar(fontsize = 20)
    ),
    Tissues = list(
      ncol = 2,
      title_gp  = gpar(fontsize = 20, fontface = "bold"),
      labels_gp = gpar(fontsize = 20)
    )
  )
)



row_ha = rowAnnotation(States=c(rownames(emission_mat)),
                       annotation_name_gp = gpar(fontsize=20),
                       annotation_legend_param = list(
                         States=list(
                           title_gp=gpar(fontsize=20,fontface="bold"),
                           labels_gp = gpar(fontsize = 20))
                       ))

# row_ha = rowAnnotation(States=c(levels(rownames(emission_mat))),
#                        col = list(States = colors.hex),
#                        annotation_name_gp = gpar(fontsize=20),
#                        annotation_legend_param = list(
#                          States=list(
#                            title_gp=gpar(fontsize=20,fontface="bold"),
#                            labels_gp = gpar(fontsize = 20))
#                        ))
col_fun2 = colorRamp2(c(0, 0.2,0.5,0.75, 1), c("white","red","yellow", "blue", "green"))
col_fun = colorRamp2(c(0,0.5,1), c("blue","white","red"))
png("plots/stack_segmentation_chr12-normEmission_40N.clustered.png",
    width = 30,
    height = 20,res = 300,units = "in")

pdf("stack_model_BB_207_samples_50N/plots/stack_segmentation_chr1-normEmission_50N.pdf",
    width = 30,
    height = 20)
Heatmap(emission_mat,show_row_names = TRUE,show_column_names=FALSE,
             show_row_dend = FALSE,
             cluster_columns = FALSE ,cluster_rows = TRUE,
              top_annotation = col_ha,
             col = col_fun,
             cluster_row_slices = FALSE,
             column_title_gp = gpar(fontsize = 20),  
             column_names_gp = gpar(fontsize = 20),   
             row_title_gp = gpar(fontsize=20),
             row_names_gp = gpar(fontsize = 20),row_gap = unit(5, "mm"),
             heatmap_legend_param = list(
               title = "Min-Max",direction="horizontal",                     
               title_gp = gpar(fontsize = 20,fontface="bold"),       
               labels_gp = gpar(fontsize = 20)        
             ))
dev.off()

