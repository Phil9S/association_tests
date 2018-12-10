rm(list = ls())
## Enable cmd line args
args = commandArgs(trailingOnly=TRUE)

## load libraries
cat("[ASSOCIATION TESTS][BEVIMED] Loading libraries\n")
suppressMessages(
  if(!require(stringr)){
    install.packages("stringr",repos = "https://mirrors.ebi.ac.uk/CRAN/")
    suppressMessages(library(stringr))
  }
)
suppressMessages(
  if(!require(stringi)){
    install.packages("stringi",repos = "https://mirrors.ebi.ac.uk/CRAN/")
    suppressMessages(library(stringi))
  }
)
suppressMessages(
  if(!require(BeviMed)){
    install.packages("BeviMed",repos = "https://mirrors.ebi.ac.uk/CRAN/")
    suppressMessages(library(BeviMed))
  }
)
suppressMessages(
  if(!require(tidyr)){
    install.packages("tidyr",repos = "https://mirrors.ebi.ac.uk/CRAN/")
    suppressMessages(library(tidyr))
  }
)
suppressMessages(
  if(!require(dplyr)){
    install.packages("dplyr",repos = "https://mirrors.ebi.ac.uk/CRAN/")
    suppressMessages(library(dplyr))
  }
)
suppressMessages(
  if(!require(data.table)){
    install.packages("data.table",repos = "https://mirrors.ebi.ac.uk/CRAN/")
    suppressMessages(library(data.table))
  }
)
suppressMessages(  
  if(!require(parallel)){
    install.packages("parallel",repos = "https://mirrors.ebi.ac.uk/CRAN/")
    suppressMessages(library(parallel))
  }
)
## Disable scientific notation
options(scipen=999)
## Load data
cat("[ASSOCIATION TESTS][BEVIMED] Loading association data\n")
load("association.RData")

if(args[1] == "PATHWAY"){
  ## Scan in data per line into vector
  x <- scan(args[2], what="", sep="\n")
  ## split into list of values and make first item the name
  y <- strsplit(x, "\t") 
  names(y) <- sapply(y, `[[`, 1)
  paths <- lapply(y, `[`, -1)
}

if(args[1] == "PATHWAY"){
  cat("[ASSOCIATION TESTS][BEVIMED] Performing PATHWAY association analysis\n")
  ## Performing matrix transformation and bevimed analysis - pathway
  results <- lapply(seq_along(paths),FUN=function(x){
                  vars <- vcf_file[vcf_file$GENE %in% as.character(unlist((paths[x]))),][c(1,18:ncol(vcf_file))]
                  rownames(vars) <- vars[,1]
                  vars <- as.matrix(t(vars[-1]))
                  c(list(path=names(paths)[[x]]),summary(bevimed(y = af_unaf,G = vars)))
                })
  
  results_table <- do.call(what=rbind, lapply(results, function(x) data.frame(
    Gene=x[["path"]],
    Prob_assoc=sum(x[["prob_association"]]),
    Prob_dominance=x[["prob_association"]]["dominant"]/sum(x[["prob_association"]]),
    Prob_recessive=x[["prob_association"]]["recessive"]/sum(x[["prob_association"]]),
    check.names=FALSE,
    stringsAsFactors=FALSE
  )))
  
  signif <- lapply(results,FUN = function(x){
    if(x[[1]] %in% results_table$Gene[results_table$Prob_assoc > 0.9])
    { return(x)
    } else 
    { return(NULL)
    }})
  ## Remove null entries in list
  signif <- signif[!sapply(signif,is.null)]
  
  ## Using new "significant calls" and pull more information regarding variants involved, genes and consequences
  signif_calls <- do.call(what=rbind, lapply(signif,function(x) data.frame(
    Pathway=x[['path']],
    Prob_gene_assoc=sum(x[["prob_association"]]),
    Variant=x[['variant_names']],
    Expected_explained=x[['models']][['dominant']]['expected_explained'],
    Explaining_variants=x[['models']][['dominant']]['explaining_variants'],
    Pathogenic_prob_dom=x[['models']][['dominant']][['conditional_prob_pathogenic']],
    Case_count=x[['models']][['dominant']][['variant_counts']][['TRUE']],
    Control_count=x[['models']][['dominant']][['variant_counts']][['FALSE']],
    check.names = F,stringsAsFactors = F
  )))
  
  if(is.null(signif_calls) == FALSE){
  signif_calls$Gene <- vcf_file$GENE[vcf_file$calcID %in% signif_calls$Variant]
  }
  
} else {
  cat("[ASSOCIATION TESTS][BEVIMED] Performing GENE association analysis\n")
  ## Performing matrix transformation and bevimed analysis - per gene
  results <- lapply(as.list(unique(vcf_file$GENE)), FUN=function(gene){
		  cat(paste("[ASSOCIATION TESTS][BEVIMED] Testing ",gene," \n",sep=""))
                  vars <- vcf_file[vcf_file$GENE == gene,][c(1,18:ncol(vcf_file))]
                  rownames(vars) <- vars[,1]
                  vars <- as.matrix(t(vars[-1]))
                  c(list(GENE=gene),summary(bevimed(y = af_unaf,G = vars)))
                })


  ## Results indexing to form output table of genes and their associations
  results_table <- do.call(what=rbind, lapply(results,function(x) data.frame(
                Gene=x[["GENE"]],
                Prob_assoc=sum(x[["prob_association"]]),
                Prob_dominance=x[["prob_association"]]["dominant"]/sum(x[["prob_association"]]),
                Prob_recessive=x[["prob_association"]]["recessive"]/sum(x[["prob_association"]]),
                check.names=FALSE,
                stringsAsFactors=FALSE
  )))
  
  ## Take results table entries where probable association exceeds 0.9 make new list
  signif <- lapply(results,FUN = function(x){
      if(x[[1]] %in% results_table$Gene[results_table$Prob_assoc > 0.9])
        { return(x)
      } else 
        { return(NULL)
      }})
  ## Remove null entries in list
  signif <- signif[!sapply(signif,is.null)]

  ## Using new "significant calls" and pull more information regarding variants involved, genes and consequences
  signif_calls <- do.call(what=rbind, lapply(signif,function(x) data.frame(
                    Gene=x[['GENE']],
                    Prob_gene_assoc=sum(x[["prob_association"]]),
                    Variant=x[['variant_names']],
                    Expected_explained=x[['models']][['dominant']]['expected_explained'],
                    Explaining_variants=x[['models']][['dominant']]['explaining_variants'],
                    Pathogenic_prob_dom=x[['models']][['dominant']][['conditional_prob_pathogenic']],
                    Case_count=x[['models']][['dominant']][['variant_counts']][['TRUE']],
                    Control_count=x[['models']][['dominant']][['variant_counts']][['FALSE']],
                    check.names = F,stringsAsFactors = F
                  )))
}
cat("[ASSOCIATION TESTS][BEVIMED] Formatting output for writing\n")
if(is.null(signif_calls) == FALSE){
  signif_calls$Consequence <- vcf_file$CONSEQ[vcf_file$calcID %in% signif_calls$Variant]
  signif_calls$ExAC_all <- vcf_file$AF[vcf_file$calcID %in% signif_calls$Variant]
  rownames(signif_calls) <- c()
}

## Write results to new output
cat("[ASSOCIATION TESTS][BEVIMED] Writing output\n")
write.table(results_table,file = "bevimed_summary.results",quote = F,sep = "\t",col.names = T,row.names = F)
if(is.null(signif_calls) == FALSE){
  write.table(signif_calls,file = "bevimed_pathogenic.results",quote = F,sep = "\t",row.names = F,col.names = T)
}
