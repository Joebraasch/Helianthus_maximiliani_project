################################
# Site frequency spectrum - HM #
################################

###
# Aim
###

# This script aims to estimate the site frequency spectrum (SFS) for each HM 
# populations independently and then for populations grouped by seed source.

##############################
# Load data and dependencies #
##############################
source("vcf2sfs.R")
mygt <- vcf2gt("Data/SNP_FINAL_noLD.vcf", "popmap.txt")

source <- unique(mygt$popmap)
for(i in 1:length(source)){
  cat(paste("#-- Now processing:", source[i], "--#\n"))
  mysfs <- gt2sfs.raw(mygt, source[i])
  mysfsF <- fold.sfs(mysfs)
  assign(source[i], mysfsF[which(mysfsF!=0)])
}
sourceloc <- list(ExSitu, Native_remnant, Producer, Selected)
names(sourceloc) <- c("Ex situ", "Wild Contemporary", "Commercial", "Selected")

# Here you will observe twice the number of individuals in mysfs as we are working
# with a diploid species.


################
# Estimate SFS #
################

#-- Seed Source --#
yrange <- max(range(sourceloc[[4]]))
names(sourceloc)
col <- c("gold", "forestgreen", "sienna", "orange")
name <- c("ExSitu", "Native_remnant", "Producer", "Selected")

# Without selected #
sourceloc <- sourceloc[-4]
names(sourceloc)
col <- c("gold", "forestgreen", "sienna")

for(i in 1:length(sourceloc)){
  if(i==1){
    sp <- as.vector(sourceloc[[i]])
    xaxis <- (0:(length(sp)-1))/(length(which(mygt$popmap==name[i]))*2)
    if(length(sourceloc)==4) plot(sp~xaxis, col=col[i], xlab = "Minor allele frequency", ylab="Number of SNPs", type="l", ylim=c(0, yrange+50))
    if(length(sourceloc)==3) plot(sp~xaxis, col=col[i], xlab = "Minor allele frequency", ylab="Number of SNPs", type="l", ylim=c(0, 350))
    points(sp~xaxis, pch=i, col=col[i])
  }
  if(i!=1){
    sp <- as.vector(sourceloc[[i]])
    xaxis <- (0:(length(sp)-1))/(length(which(mygt$popmap==name[i]))*2)
    lines(sp~xaxis, col=col[i], pch=i)
    points(sp~xaxis, pch=i, col=col[i])
  }
}
legend("topright", pch = c(1:4), col=c("gold", "forestgreen", "sienna", "orange"), legend = names(sourceloc), bty = "n")
title("Folded Site Frequency Spectrums (SFS)")
