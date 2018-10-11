library(SomaticSignatures)
library(VariantAnnotation)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(ggplot2)
library(gtools)
library(gridExtra)
library(grid)

set.seed(1)

setwd("~/Documents/Uni/Master/mutsig")


get_mut_signatures <- function(vcf, n_sigs){
  #read VCF files as VRange object
  vrange <- readVcfAsVRanges(vcf)
  motifs = mutationContext(vrange, BSgenome.Hsapiens.1000genomes.hs37d5)
  
  #get mutation spectrum 
  mm = motifMatrix(motifs, group = "SID",normalize = TRUE)
  mm = mm[, colSums(mm != 0) > 0]
  mm = mm[rowSums(mm!=0) > 0,]

  #Decompose the motif spectrum 
  sigs_nmf = identifySignatures(mm, n_sigs, nmfDecomposition)
  
  return(sigs_nmf)
} 

compute_euc_distance <- function(array1,array2){
  euc_distance <- sqrt(sum((array1-array2) ^ 2))
  return(euc_distance)
}




#Set K based on reference paper
ks <- c(3,2,3,4,3,3,5,4,3,3,5,3,3,3,3,4,4,3,4,3,4,3,3,3,4,4,2,5,4,4)

#ALL/AML test
ALL <- "data/VCF/ALL.vcf"
AML <- "data/VCF/AML.vcf"
all_sigs <- get_mut_signatures(ALL,ks[1])
aml_sigs <- get_mut_signatures(AML,ks[2])

for(aml_idx in c(1:ks[2])){
  for (all_idx in c(1:ks[1])){
    print(compute_euc_distance(aml_sigs@signatures[,aml_idx],all_sigs@signatures[,all_idx]))
  }
}

#For all cancer classes
files <- list.files(path="data/VCF", full.names=TRUE, recursive=FALSE)

all_mutational_sigs <- list()
unique_sigs <- list()
cancer_classes <- list()
euc_distances <- c()
for (idx in c(1:length(files))){
  print(files[idx])
  f_sigs <- get_mut_signatures(files[idx],ks[idx])
  all_mutational_sigs[basename(files[idx])] <- f_sigs
  
  for (sig_idx in c(1:ks[idx])){
    merged <- FALSE
    for (uni_sig_idx in seq_along(unique_sigs)){
      dist <- compute_euc_distance(unique_sigs[[uni_sig_idx]],f_sigs@signatures[,sig_idx])
      #merge signatures if euc distance is below 7
      if(dist<7){
        unique_sigs[[uni_sig_idx]] <- rowMeans(cbind(unique_sigs[[uni_sig_idx]],f_sigs@signatures[,sig_idx]),na.rm = TRUE)
        cancer_classes[[uni_sig_idx]] <- c(cancer_classes[[uni_sig_idx]],paste(basename(files[idx])))
        merged <- TRUE
        break
      }
    }
    if (merged == FALSE){
      unique_sigs[[paste(basename(files[idx]),sig_idx,sep = "_")]] <- f_sigs@signatures[,sig_idx]
      cancer_classes[[paste(basename(files[idx]),sig_idx,sep = "_")]] <- c(paste(basename(files[idx])))
    }
    euc_distances <- c(euc_distances,lapply(unique_sigs,compute_euc_distance,array1 = f_sigs@signatures[,sig_idx]))
  }
}


#FIGURES

cancer_classes <- lapply(cancer_classes,unique)
plots <- list()
for(name in names(unique_sigs)){
  cancer_number <- strsplit(name,"_")
  cancer <- cancer_number[[1]][1]
  number <- cancer_number[[1]][2]
  all_mutational_sigs[cancer][[1]]@signatures[,paste('S',number,sep ="")] <- unique_sigs[name][[1]]
  plots[[name]] <- plotSignatures(all_mutational_sigs[cancer][[1]], percent = TRUE) +coord_cartesian(ylim = c(0, 20)) 
}

#save plots

for (idx in seq_along(plots)){
  tiff(paste('images/signature',idx,".tiff",sep=""), units="in", width=5, height=5, res=300)
  print(plots[[idx]] + ggtitle(paste("Signature",idx,sep=" ")))
  dev.off()
}

#save signature overview plot

ys <- c()
xs <- c()

cancers <- lapply(files,function(x){strsplit(basename(x),"\\.")[[1]][1]})
  
for (idx in seq_along(unique_sigs)){
  y <- idx
  for (cancer_class in cancer_classes[[idx]]){
    xs <- c(xs,which(strsplit(cancer_class,"\\.")[[1]][1] == cancers))
    ys <- c(ys,y)
  }
}

plotdata <- data.frame(cbind(xs,ys))

tiff('images/cancer_class_distribution.tiff', units="in", width=10, height=8, res=300)
ggplot(data = plotdata) + geom_point(aes(x = xs, y = ys), color = "darkgreen") + scale_x_continuous(breaks = c(1:30),labels = cancers) +
  scale_y_continuous(breaks = pretty(plotdata$ys, n = 27)) + labs(x = "Cancer classes", y = "Signatures") + theme(axis.text.x = element_text(angle=45, hjust =1))
dev.off()

