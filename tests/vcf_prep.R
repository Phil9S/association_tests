rm(list = ls())
## Enable cmd line args
args = commandArgs(trailingOnly=TRUE)

## load libraries
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

## Disable scientific notation
options(scipen=999)

print(args[1])
print(args[2])
print(args[3])
print(args[4])
print(args[5])
print(args[6])
print(args[7])
print(args[8])

## Load list of excluded samples
samples_rm <- c()

#"synonymous_SNV","nonframeshift_deletion","nonframeshift_insertion","stoploss","nonsynonymous_SNV"
CONSEQS <- c("nonframeshift_deletion","nonframeshift_insertion","frameshift_deletion","stopgain","frameshift_insertion","splicing","nonsynonymous_SNV","synonymous_SNV","stoploss")
CONSEQS <- CONSEQS[CONSEQS %in% as.character(unlist(strsplit(args[8],",")))]

## Load column header from vcf file 
col_headers <- fread(args[1],nrows = 1,header = F,sep="\t")
col_headers <- gsub(pattern = "#",replacement = "",col_headers)
## Load vcf data and append to column headers
vcf_file <- fread(input = args[1],skip = 1,stringsAsFactors = F,header = F,,sep="\t")
names(vcf_file) <- col_headers
vcf_file <- as.data.frame(vcf_file)

## Load case/control T/F data - Sort them into the same order as vcf (vcftools and annovar like to switch column orders) - assign T/F vector to af_unaf
samples <- read.table(args[2])
samples <- samples[match(names(vcf_file[10:ncol(vcf_file)]),samples$V1),]
samples <- samples[!samples$V1 %in% samples_rm,]
af_unaf <- as.logical(samples[,2])

## Converting genotypes
out <- mclapply(vcf_file[10:ncol(vcf_file)],mc.preschedule = T,mc.cores = 2, function(x){
  x <- stri_replace_all_regex(x,pattern = "^0/0.*",replacement = 0)
  x <- stri_replace_all_regex(x,pattern = "^0/1.*",replacement = 1)
  x <- stri_replace_all_regex(x,pattern = "^1/1.*",replacement = 2)
  x <- stri_replace_all_regex(x,pattern = "^\\./\\..*",replacement = -9)
})

vcf_file[10:ncol(vcf_file)] <- do.call(cbind.data.frame, out)

## Coerce all genotype columns as.numeric 
vcf_file[,10:ncol(vcf_file)] <- as.data.frame(sapply(vcf_file[,10:ncol(vcf_file)], function(f) as.numeric(as.character(f))),stringsAsFactors=F)

## Adding internal AFs
cases <- as.character(samples$V1[samples$V2 == TRUE])
controls <- as.character(samples$V1[samples$V2 == FALSE])

intAF <- data.frame(intAF_cases=as.numeric(seq_len(nrow(vcf_file))),
                    intAF_controls=as.numeric(seq_len(nrow(vcf_file))),
                    intAF_set=as.numeric(seq_len(nrow(vcf_file))))

intAF$intAF_cases <- signif(apply(vcf_file[,which(names(vcf_file) %in% cases)],1,function(x) sum(x != 0) / length(cases)),digits = 2)
intAF$intAF_controls <- signif(apply(vcf_file[,which(names(vcf_file) %in% controls)],1,function(x) sum(x != 0) / length(controls)),digits = 2)
intAF$intAF_set <- signif(apply(vcf_file[,which(names(vcf_file) %in% c(cases,controls))],1,function(x) sum(x != 0) / length(c(cases,controls))),digits = 2)

vcf_file <- cbind(intAF,vcf_file)
## Retrieving annotation information on Allele freq, region function, and mutation consequence
info_split <- as.data.frame(
  str_split(string = gsub(".*Func\\.refGene=(.*?);.*Gene\\.refGene=(.*?);.*ExonicFunc\\.refGene=(.*?);.*ExAC_ALL=(.*?);.*",
                          replacement = "\\1#\\2#\\3#\\4", 
                          x = vcf_file$INFO),
            pattern = "#",
            simplify = T),
  stringsAsFactors = F)

## Naming annoation information - replacing uncoded ; values added by annovar for some reason - Set unknown/novel AF to 0 and setting to numeric
names(info_split) <- c("FUNC","GENE","CONSEQ","AF")
info_split$FUNC <- gsub(pattern = "\\\\x3b",replacement = ";",info_split$FUNC)
info_split$GENE <- gsub(pattern = "\\\\x3b",replacement = ";",info_split$GENE)
info_split$AF[info_split$AF == "."] <- 0
info_split$AF <- as.numeric(info_split$AF)

## bind new columns to vcf data
vcf_file <- cbind(info_split,vcf_file)

## Adding chr:pos ids to rsID column and replacing "." missing value - making unique ID for each variant calcID - binding calcID to vcf data
vcf_file$ID[vcf_file$ID == "."] <- paste(vcf_file$CHROM[vcf_file$ID == "."],vcf_file$POS[vcf_file$ID == "."],sep = ":")
calcID <- paste(vcf_file$ID,vcf_file$REF,vcf_file$ALT,sep = "_")
vcf_file <- cbind(calcID,vcf_file)

## Clean up
rm(info_split,calcID,out)

## Add splicing anntation to CONSEQ field
vcf_file$CONSEQ[vcf_file$FUNC == "splicing" | vcf_file$FUNC == "exonic;splicing"] <- "splicing"

## PCA analysis
vcf_file_PCA <- vcf_file[vcf_file$AF > 0.05,18:ncol(vcf_file)]
vcf_file_PCA_t <- as.data.frame(t(vcf_file_PCA))
vcf_file_PCA_t <- vcf_file_PCA_t[,which(apply(vcf_file_PCA_t,2,var)!=0)]

PCA.out <- prcomp(vcf_file_PCA_t,center = T, scale. = T)
PCAs <- as.data.frame(PCA.out$x[,1:5])
rm(vcf_file_PCA,vcf_file_PCA_t,PCA.out)

colours <- ifelse(af_unaf == T,"blue","grey")
colours[which(abs(PCAs$PC1) > sd(PCAs$PC1)*3 & abs(PCAs$PC2) > sd(PCAs$PC2)*3)] <- "red"
colours[which(abs(PCAs$PC2) > sd(PCAs$PC2)*3 & abs(PCAs$PC3) > sd(PCAs$PC3)*3)] <- "red"

png(filename = "association_test_PCA.png",width = 16,height = 8,units = "in",res = 600)
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
plot(abs(PCAs$PC1),abs(PCAs$PC2),col=colours,xlab = "PC1",ylab = "PC2",main = "Association tests - PC1 ~ PC2",sub = "Red = Excluded")
abline(v=sd(PCAs$PC1)*3, col="red")
abline(h=sd(PCAs$PC2)*3, col="red")

plot(abs(PCAs$PC2),abs(PCAs$PC3),col=colours,xlab = "PC2",ylab = "PC3",main = "Association tests - PC2 ~ PC3",sub = "Red = Excluded")
abline(v=sd(PCAs$PC2)*3, col="red")
abline(h=sd(PCAs$PC3)*3, col="red")
dev.off()

PCA_excluded <- unique(c(rownames(PCAs[abs(PCAs$PC1) > sd(PCAs$PC1)*3 & abs(PCAs$PC2) > sd(PCAs$PC2)*3,]),
                  rownames(PCAs[abs(PCAs$PC2) > sd(PCAs$PC2)*3 & abs(PCAs$PC3) > sd(PCAs$PC3)*3,])))

writeLines(PCA_excluded,sep = "\n",con = "PCA_excluded.samples")
## Adjust Samples and data to remove samples
PCAs <- PCAs[!rownames(PCAs) %in% PCA_excluded,]
samples <- samples[!samples$V1 %in% PCA_excluded,]

af_unaf <- as.logical(samples[,2])

vcf_file <- vcf_file[!colnames(vcf_file) %in% PCA_excluded]

## Filtering vcf data on various traits - Rarity, Function, Consequence 
vcf_file <- vcf_file[vcf_file$CONSEQ %in% CONSEQS,]
vcf_file <- vcf_file[vcf_file$AF < args[3],] # 0.005
vcf_file <- vcf_file[vcf_file$intAF_set < args[4],] # 0.05
vcf_file <- vcf_file[vcf_file$intAF_cases < args[5],] # 0.2
vcf_file <- vcf_file[vcf_file$intAF_controls < args[6],] #0.2

save(vcf_file,PCAs,samples,af_unaf,file="association.RData")
