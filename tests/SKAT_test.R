## Script for SKAT-O implementation
## Environment clean and library load
rm(list=ls())

## load libraries
if(!require(stringr)){
  install.packages("stringr",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(stringr)
}

if(!require(stringi)){
  install.packages("stringi",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(stringi)
}

if(!require(tidyr)){
  install.packages("tidyr",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(tidyr)
}

if(!require(dplyr)){
  install.packages("dplyr",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(dplyr)
}

if(!require(data.table)){
  install.packages("data.table",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(data.table)
}
if(!require(parallel)){
  install.packages("parallel",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(parallel)
}
if(!require(SKAT)){
  install.packages("SKAT",repos = "https://mirrors.ebi.ac.uk/CRAN/")
  library(SKAT)
}

tryCatchAdv <- function(expr)
{
  
  # Initial settings
  V <- NA
  S <- "succeeded"
  M <- NA
  
  # Warning handler
  w.handler <- function(w){
    
    # Record information about warning
    S <<- "warning"
    M <<- w
    # <<- is used for assignment outside the function scope (i.e. in the external environment)
    # http://stackoverflow.com/questions/2628621/how-do-you-use-scoping-assignment-in-r
    
    # Execute expression again, suppressing warnings
    invokeRestart("muffleWarning")
    
  }
  
  # Error handler
  e.handler <- function(e){
    
    # Record information about error
    S <<- "error"
    M <<- e
    
    # Return NA as result
    return(NA)
    
  }
  
  # Try to execute the expression, use the above handlers
  V <- withCallingHandlers(tryCatch(expr, error=e.handler), warning=w.handler)
  
  # Return value
  list(value = V,
       status = S,
       message = M)
}


## Load data
load("association.RData")

## Set missing to 9
vcf_file[vcf_file == -9] <- 9

## Split by gene
split_region <- split(vcf_file,as.factor(vcf_file$GENE))

## Weights
allele_Freq <- vcf_file$AF
weights <- Get_Logistic_Weights_MAF(MAF = allele_Freq)

## Setting binary phenotype
binary <- af_unaf
binary[binary == T] <- 1
binary[binary == F] <- 0

## Generating matrix_list
gene_matrix_list <- mapply(function(x, i){
  mat <- as.data.frame(x)
  if(nrow(mat) > 1){
  col.n <- colnames(mat[c(18:ncol(mat))])
  row.n <- as.character(mat[,1])
  mat <- as.matrix(t(mat[c(18:ncol(mat))]))
  rownames(mat) <- col.n
  colnames(mat) <- row.n
  return(mat)
  }
}, split_region, names(split_region),SIMPLIFY = F)

## SKAT null model
obj <- SKAT_Null_Model(binary ~ PCAs$PC1, out_type="D")

## Performing binary SKAT
out <- mapply(function(x, i){
  tryCatchAdv({
  cat(paste("[ASSOCIATION TESTS][SKAT] Testing ",i," \n",sep = ""))
  SKATBinary(as.matrix(x), obj, method = "SKAT", kernel = "linear.weighted")
  })
}, gene_matrix_list, names(gene_matrix_list),SIMPLIFY = F)

# ## Generating matrix_list using exac weights
# gene_matrix_list <- mapply(function(x, i){
#   mat <- as.data.frame(x)
#   col.n <- colnames(mat[c(5,18:ncol(mat))])
#   row.n <- as.character(mat[,1])
#   mat <- as.matrix(t(mat[c(5,18:ncol(mat))]))
#   rownames(mat) <- col.n
#   colnames(mat) <- row.n
#   return(mat)
# }, split_region, names(split_region),SIMPLIFY = F)
# 
# ## SKAT null model
# obj <- SKAT_Null_Model(binary ~ PCAs$PC1, out_type="D")
# 
# ## Performing binary SKAT
# out <- mapply(function(x, i){
#   #tryCatch({
#   cat(paste("[ASSOCIATION TESTS][SKAT] Testing ",i," \n",sep = ""))
#   SKATBinary(as.matrix(x[-1,]), obj, method = "SKAT", kernel = "linear.weighted",weights = dbeta(x[1,],1,25))
#   #},error=function(e) NULL)
# }, gene_matrix_list, names(gene_matrix_list),SIMPLIFY = F)

# Error in SKATExactBin.SKATO_GetQParam(Z, res, idx, 0, pi1, pr$prob_k,  :
#                                         Error in ResampleSTAT!
#                                         Calls: mapply ... SKATExactBin.Moment -> SKATExactBin.SKATO_GetQParam
#                                       In addition: There were 50 or more warnings (use warnings() to see the first 50)
#                                       

skat_results <- do.call(rbind,mapply(function(x, i){
    if(length(x$value) == 1){
      data.frame(gene=i,pvalue=NA,warning=x$status,description=x$message$message)
    } else {
    unique(data.frame(gene=i,pvalue=x$value$p.value,warning=x$status,description=ifelse(!is.na(x$message),x$message$message,NA)))
    }
}, out, names(out),SIMPLIFY = F))

skat_results <- skat_results[order(skat_results$pvalue),]

write.table(skat_results,"skat_output.results",quote=F,row.names=F,col.names=T,sep="\t")

png(filename = "SKAT_QQ.png",width = 8,height = 8,units = "in",res = 600)
QQPlot_Adj(Pval = as.numeric(unlist(lapply(out,function(x) if(length(x$value) > 1 ){x$value$p.value}))),
           MAP = as.numeric(unlist(lapply(out,function(x) if(length(x$value) > 1 ){x$value$MAP}))))
dev.off()
