#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(bsseq))
suppressPackageStartupMessages(library("BiocParallel"))

# Reading arguments
data=read.table(args[1],header=F)
smoothing=as.integer(args[2])
target<-paste(args[3])
cores=as.integer(args[4])
options(scipen = 100, digits = 4)

# save all columns as vectors
chromosome <- as.vector(as.matrix(data[1]))
position <- as.matrix(data[2])
end <- as.matrix(data[3])
coverage <-as.matrix(data[5]) # Since converted from the beta file using wgbstools that's why coverage is 5th column
methylation <-as.matrix(data[4]) # Since converted from the beta file using wgbstools that's why Methylated reads is 4th column


# colnames for current version has to be set to NULL
colnames(methylation) <- NULL
colnames(coverage) <- NULL

# create bsseq object
BS_methyl_smooth <- BSseq(M=methylation, Cov=coverage, pos= position, chr=chromosome) 

# conditional ns value based on smoothing window 

if (smoothing == 1000) {
    nCpGs=70
} else if (smoothing ==  200) {
    nCpGs=25
} else {
    nCpGs=50
}
print(paste("Smoothing setting used nCPGs:", nCpGs))
print(paste("Smoothing setting used smoothing window:", smoothing))



# perform smoothing
BS_methyl_after_smooth <- BSmooth(BS_methyl_smooth, ns=nCpGs,h=smoothing,BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))

# Storing the  smooth values
meth=getMeth(BS_methyl_after_smooth,type="smooth")
cov=getCoverage(BS_methyl_after_smooth,type="Cov")
smooth_meth=round(cov*meth,0) # estimate the methylated reads 

data_bed_file = data.frame(chromosome,position,end,cov,smooth_meth)
colnames(data_bed_file) <- c("chr","start","end","Cov","Meth")

# remove chr coordinates to fit to ESMM 
data_bed_file$chr <- sub("^chr", "", data_bed_file$chr)

# correct for methylation values if greater than coverage 
data_bed_file$Meth <- ifelse(
  data_bed_file$Meth > data_bed_file$Cov,
  data_bed_file$Cov,
  data_bed_file$Meth
)

write.table(data_bed_file,file=paste(target,"_",smoothing,"_smooth.txt",sep=""),quote=FALSE,sep="\t",col.names=T,row.names=F)


print("Smoothing ended")
rm(list=ls())
quit()