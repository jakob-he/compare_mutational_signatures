library(VariantAnnotation)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(signeR)
library(rtracklayer)


mut <- read.table(system.file("extdata","21_breast_cancers.mutations.txt",
                              package="signeR"), header=TRUE, check.names=FALSE)
opp <- read.table(system.file("extdata","21_breast_cancers.opportunity.txt",
                              package="signeR"))

set.seed(1)

setwd("~/Documents/Uni/Master/mutsig")

get_mut_signatures <- function(mut_counts){
  #load mutation counts and generate oppurtunity matrix
  #target_regions <- import(con="data/Exons_refseq_genes", format="bed")
  mut <- read.table(mut_counts, header=TRUE, check.names=FALSE, sep =',', row.names = 1)
  mut = mut[rowSums(mut!=0) > 0,]
  #opp <- genOpportunityFromGenome(BSgenome.Hsapiens.1000genomes.hs37d5,target_regions,nsamples=nrow(mut))
  
  #run EM
  signatures <- signeR(M=mut, Opport=NA, main_eval=100, EM_eval=50, EMit_lim=20, nlim=c(1,7))
  
  return(signatures)
} 

#Set K based on reference paper
#ks <- c(3,2,3,4,3,3,5,4,3,3,5,3,3,3,3,4,4,3,4,3,4,3,3,3,4,4,2,5,4,4)


#For all cancer classes
files <- list.files(path="data/Cmatrix/", full.names=TRUE, recursive=FALSE)

all_sigs <- list()
number_of_sigs <- c()
for (idx in seq_along(files)){
  signatures <- get_mut_signatures(files[idx])
  all_sigs[[length(all_sigs)+1]] <- signatures
  number_of_sigs <- c(number_of_sigs,signatures$Nsign)
}

for (idx in seq_along(all_sigs)){
  tiff(paste('figures/signeR/',strsplit(basename(files[idx]),"\\.")[[1]][1],".png",sep=""), units="in", width=8, height=5, res=200)
  print(SignPlot(all_sigs[[idx]]$SignExposures))
  dev.off()
}

#save(all_sigs,file = "data/R/signer_signatures.RData")
#save(number_of_sigs,file = "data/R/signer_signature_distribution.RData")
