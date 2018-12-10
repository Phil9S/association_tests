## Script counts the number of affected samples harbouring 1 or more variants in a gene
## Environment clean and library load
rm(list=ls())

## load libraries
cat("[ASSOCIATION TESTS][FISHERS] Loading libraries\n")
suppressMessages(
  if(!require(stringr)){
    install.packages("stringr",repos = "https://mirrors.ebi.ac.uk/CRAN/")
    library(stringr)
  }
)
suppressMessages(
  if(!require(stringi)){
    install.packages("stringi",repos = "https://mirrors.ebi.ac.uk/CRAN/")
    library(stringi)
  }
)
suppressMessages(
  if(!require(tidyr)){
    install.packages("tidyr",repos = "https://mirrors.ebi.ac.uk/CRAN/")
    library(tidyr)
  }
)
suppressMessages(
  if(!require(dplyr)){
    install.packages("dplyr",repos = "https://mirrors.ebi.ac.uk/CRAN/")
    library(dplyr)
  }
)
suppressMessages(
  if(!require(data.table)){
    install.packages("data.table",repos = "https://mirrors.ebi.ac.uk/CRAN/")
    library(data.table)
  }
)
suppressMessages(
  if(!require(parallel)){
    install.packages("parallel",repos = "https://mirrors.ebi.ac.uk/CRAN/")
    library(parallel)
  }
)
## Load data
cat("[ASSOCIATION TESTS][FISHERS] Loading association data\n")
load("association.RData")

## Fishers by site
fisher_variants <- vcf_file[c(1,18:ncol(vcf_file))]
fisher_variants[fisher_variants == -9] <- NA

cat("[ASSOCIATION TESTS][FISHERS] Generating counts matrix - Individual sites\n")
counts.matrix <- t(apply(as.matrix(fisher_variants[-1]),1,function(x) {
                    c(sum(x[af_unaf] == 1,na.rm = T) + sum(x[af_unaf] == 2,na.rm = T)*2,
                    sum(x[!af_unaf] == 1,na.rm = T) + sum(x[!af_unaf] == 2,na.rm = T)*2,
                    sum(x[af_unaf] == 0,na.rm = T)*2,
                    sum(x[!af_unaf] == 0,na.rm = T)*2)
                  }))
cat("[ASSOCIATION TESTS][FISHERS] Performing Fisher's exact test - Individual sites\n")
fisher.out <- apply(counts.matrix,1,function(x) fisher.test(matrix(x,nrow = 2),alternative = "g")[[1]])
cat("[ASSOCIATION TESTS][FISHERS] Performing multiple testing correction - Individual sites\n")
fdr.adj <- p.adjust(fisher.out,method = "fdr",length(fisher.out))
bonferroni.adj <- p.adjust(fisher.out,method = "bonferroni",length(fisher.out))

fisher_sites <- cbind(vcf_file$calcID,vcf_file$GENE,as.data.frame(counts.matrix),fisher.out,fdr.adj,bonferroni.adj)
names(fisher_sites) <- c("Site","Gene","Case_Mut_AC","Cont_Mut_AC","Case_Ref_AC","Cont_Ref_AC","pvalue","FDR_adj","Bonferroni_adj")
fisher_sites <- fisher_sites[order(fisher_sites$pvalue),]

rm(counts.matrix,fisher_variants,bonferroni.adj,fdr.adj,fisher.out)

## Filter by consequence - Truncations
cat("[ASSOCIATION TESTS][FISHERS] Generating counts matrix - By gene | By Truncation\n")
fisher_agg_trunc <- vcf_file[vcf_file$CONSEQ %in% c("stopgain","frameshift deletion","frameshift insertion"),][c(3,18:ncol(vcf_file))]
fisher_agg_trunc[fisher_agg_trunc == -9] <- NA

fisher_gene_agg <- as.data.frame(fisher_agg_trunc %>% group_by(GENE) %>% summarise_all(sum))

counts.matrix <- t(apply(as.matrix(fisher_gene_agg[-1]),1,function(x) {
  c(sum(x[af_unaf] > 0,na.rm = T),
    sum(x[!af_unaf] > 0,na.rm = T),
    sum(x[af_unaf] == 0,na.rm = T),
    sum(x[!af_unaf] == 0,na.rm = T))
}))
cat("[ASSOCIATION TESTS][FISHERS] Performing Fisher's exact test - By gene | By Truncation\n")
fisher.out <- apply(counts.matrix,1,function(x) fisher.test(matrix(x,nrow = 2),alternative = "g")[[1]])
cat("[ASSOCIATION TESTS][FISHERS] Performing multiple testing correction - By gene | By Truncation\n")
fdr.adj <- p.adjust(fisher.out,method = "fdr",length(fisher.out))
bonferroni.adj <- p.adjust(fisher.out,method = "bonferroni",length(fisher.out))

fisher_gene_trunc <- cbind(fisher_gene_agg$GENE,as.data.frame(counts.matrix),fisher.out,fdr.adj,bonferroni.adj)
names(fisher_gene_trunc) <- c("Gene","Case_Aff","Cont_Aff","Case_Unaff","Cont_Unaff","pvalue","FDR_adj","Bonferroni_adj")
fisher_gene_trunc <- fisher_gene_trunc[order(fisher_gene_trunc$pvalue),]

rm(counts.matrix,fisher_agg_trunc,fisher_gene_agg,bonferroni.adj,fdr.adj,fisher.out)

## Filter by consequence - nonsynonymous
cat("[ASSOCIATION TESTS][FISHERS] Generating counts matrix - By gene | By nonsynonymous\n")
fisher_agg_nsyn <- vcf_file[vcf_file$CONSEQ == "nonsynonymous_SNV",][c(3,18:ncol(vcf_file))]
fisher_agg_nsyn[fisher_agg_nsyn == -9] <- NA

fisher_gene_agg <- as.data.frame(fisher_agg_nsyn %>% group_by(GENE) %>% summarise_all(sum))

counts.matrix <- t(apply(as.matrix(fisher_gene_agg[-1]),1,function(x) {
  c(sum(x[af_unaf] > 0,na.rm = T),
    sum(x[!af_unaf] > 0,na.rm = T),
    sum(x[af_unaf] == 0,na.rm = T),
    sum(x[!af_unaf] == 0,na.rm = T))
}))

cat("[ASSOCIATION TESTS][FISHERS] Performing Fisher's exact test - By gene | By nonsynonymous\n")
fisher.out <- apply(counts.matrix,1,function(x) fisher.test(matrix(x,nrow = 2),alternative = "g")[[1]])
cat("[ASSOCIATION TESTS][FISHERS] Performing multiple testing correction - By gene | By nonsynonymous\n")
fdr.adj <- p.adjust(fisher.out,method = "fdr",length(fisher.out))
bonferroni.adj <- p.adjust(fisher.out,method = "bonferroni",length(fisher.out))

fisher_gene_nsyn <- cbind(fisher_gene_agg$GENE,as.data.frame(counts.matrix),fisher.out,fdr.adj,bonferroni.adj)
names(fisher_gene_nsyn) <- c("Gene","Case_Aff","Cont_Aff","Case_Unaff","Cont_Unaff","pvalue","FDR_adj","Bonferroni_adj")
fisher_gene_nsyn <- fisher_gene_nsyn[order(fisher_gene_nsyn$pvalue),]

rm(counts.matrix,fisher_agg_nsyn,fisher_gene_agg,bonferroni.adj,fdr.adj,fisher.out)

## Filter by gene
cat("[ASSOCIATION TESTS][FISHERS] Generating counts matrix - By gene \n")
fisher_agg_all <- vcf_file[c(3,18:ncol(vcf_file))]
fisher_agg_all[fisher_agg_all == -9] <- NA

fisher_gene_agg <- as.data.frame(fisher_agg_all %>% group_by(GENE) %>% summarise_all(sum))

counts.matrix <- t(apply(as.matrix(fisher_gene_agg[-1]),1,function(x) {
  c(sum(x[af_unaf] > 0,na.rm = T),
    sum(x[!af_unaf] > 0,na.rm = T),
    sum(x[af_unaf] == 0,na.rm = T),
    sum(x[!af_unaf] == 0,na.rm = T))
}))

cat("[ASSOCIATION TESTS][FISHERS] Performing Fisher's exact test - By gene\n")
fisher.out <- apply(counts.matrix,1,function(x) fisher.test(matrix(x,nrow = 2),alternative = "g")[[1]])
cat("[ASSOCIATION TESTS][FISHERS] Performing multiple testing correction - By gene\n")
fdr.adj <- p.adjust(fisher.out,method = "fdr",length(fisher.out))
bonferroni.adj <- p.adjust(fisher.out,method = "bonferroni",length(fisher.out))

fisher_gene_all<- cbind(fisher_gene_agg$GENE,as.data.frame(counts.matrix),fisher.out,fdr.adj,bonferroni.adj)
names(fisher_gene_all) <- c("Gene","Case_Aff","Cont_Aff","Case_Unaff","Cont_Unaff","pvalue","FDR_adj","Bonferroni_adj")
fisher_gene_all <- fisher_gene_all[order(fisher_gene_all$pvalue),]

rm(counts.matrix,fisher_agg_all,fisher_gene_agg,bonferroni.adj,fdr.adj,fisher.out)

cat("[ASSOCIATION TESTS][FISHERS] Writing output files\n")
write.table(fisher_sites,"fisher_sites_test.results",quote = F,col.names = T,row.names = F,sep = "\t")
write.table(fisher_gene_all,"fisher_gene_all_test.results",quote = F,col.names = T,row.names = F,sep = "\t")
write.table(fisher_gene_nsyn,"fisher_gene_nsyn_test.results",quote = F,col.names = T,row.names = F,sep = "\t")
write.table(fisher_gene_trunc,"fisher_gene_trunc_test.results",quote = F,col.names = T,row.names = F,sep = "\t")
