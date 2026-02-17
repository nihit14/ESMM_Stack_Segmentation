library(jsonlite)
args = commandArgs(trailingOnly=TRUE)

# arguments
file = args[1]
n=as.integer(args[2])
output=paste(args[3])
name=as.character(args[4])
model=as.character(args[5])

# extract data and filename
filename=tools::file_path_sans_ext(basename(args[1]))
json_data <- fromJSON(file,flatten=TRUE[[1]])


NBI_mean <- function(p,r){
  pn = as.numeric(p)
  rn = as.numeric(r)
  p = pn
  r = rn
  dist_mean = r * (1-p)/p
  return(dist_mean)
}

SI_mean <- function(mu,sigma,v){
  mn = as.numeric(mu)
  sn = as.numeric(sigma)
  vn = as.numeric(v)
  mu = mn
  sigma = sn
  v = vn
  alpha = sqrt((1 / (sigma*sigma)) + (2*mu / sigma))
  w = sqrt(mu*mu + alpha*alpha) - mu
  dist_mean = mu * besselK(w,v+1) / besselK(w,v)
  return(dist_mean)
}

BNB_mean <- function(r, alpha,beta){
  alphan = as.numeric(alpha)
  betan = as.numeric(beta)
  rn = as.numeric(r)
  alpha = alphan
  beta =betan
  r = rn
  dist_mean = r*beta/(alpha -1)
  return(dist_mean)
}
rows= length(json_data$marker)
number_of_histone_marks = rows

BB_mean <- function(alpha,beta){
  alphan <- as.numeric(alpha)
  betan <- as.numeric(beta)
  alpha <- alphan
  beta <- betan
  dist_mean <- alpha/(alpha+beta)
  return(dist_mean)
}

if (model== "Coverage"){
  rows= length(json_data$coverage_marker)
  number_of_histone_marks = rows
} else {
  rows= length(json_data$marker)
  number_of_histone_marks = rows
}

dist_mean = data.frame(matrix(nrow=rows,ncol=n))
for (h in seq_len(n)){
  df = as.data.frame(json_data$emission[[h]])
  for (i in seq_len(rows)){
    dist <- df$distribution[i]
    if (dist == "NBI") {
      dist_mean[i,h] <- NBI_mean(df$parameters.p[i], df$parameters.r[i])
    } else if (dist == "BNB") {
      dist_mean[i,h] <- BNB_mean(df$parameters.r[i], df$parameters.alpha[i], df$parameters.beta[i])
    } else if (dist == "SI") {
      dist_mean[i,h] <- SI_mean(df$parameters.mu[i], df$parameters.sigma[i], df$parameters.v[i])
    } else if (dist == "BB") {
      dist_mean[i,h] <- BB_mean(df$parameters.alpha[i], df$parameters.beta[i])
    }
  }
}


States=as.data.frame(c(seq_len(n)))
df_new = cbind(States,t(dist_mean))
colnames(df_new) = c("States",json_data$marker)
rownames(df_new) = NULL

write.table(df_new,file=paste(output,paste("emissions",name,filename,n,"N.txt",sep="_"),sep="/"),sep="\t",quote=FALSE,row.names=F) 