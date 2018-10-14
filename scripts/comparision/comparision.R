library(signeR)
library(SomaticSignatures)
library(ggplot2)


calc_dist <- function(vector1,vector2){
  if(length(vector1) != length(vector2)){
    print("Arrays have to be of the same length!")
  }
  else{
    distances <- c()
    for (idx in seq_along(vector1)){
      distances <- c(distances, abs(vector1[idx]-vector2[idx]))
    }
    return(distances)
  }
}

#set working directory
setwd("~/Documents/Uni/Master/mutsig")

#load signatures data
load("data/R/somaticsignatures_signatures_rssk.RData")
som_sigs <- all_sigs
somsigs_dist <- c(4,5,3,2,3,2,3,2,4,5,2,4,6,6,3,3,3,2,5,2,3,2,3,4,5,2,4,3,3,3)

load("data/R/signer_signatures.RData")
signer_sigs <- all_sigs
load("data/R/signer_signature_distribution.RData")
signer_dist <- number_of_sigs

pmsig_dist <-  c(3,2,3,4,3,3,5,4,3,3,5,3,3,3,3,4,4,3,4,3,4,3,3,3,4,4,2,5,4,4)

#compare number of signatures
distances <- calc_dist(signer_dist,somsigs_dist)
