#----------------------------#
# Loading necessary packages #
#----------------------------#
library(adegenet)
library(hierfstat)
library(boot)
library(simpleboot)
library(pegas)
library(strataG)
library(ggplot2)
#library(radiator)
library(svMisc)

#-----------------------#
# clean_genind function #
#-----------------------#

# This function takes a sub-sampled genind object and cleans it to highlight putative private alleles

clean_genind <- function(data){
  df <- genind2df(data, sep = ",")
  gen <- df2genind(df[,-1], sep = ",",  ind.names = rownames(df), loc.names = colnames(df)[2:ncol(df)], 
                   pop = df$pop, ploidy = 2, type = "codom", ncode = NULL)
  return(gen)
}

#----------------------------------#
# Ho_Fst_down_sampling_v2 function #
#----------------------------------#

# This function down-samples SNPs based on observed heterozygocity and Fst estimates (bivariate approach)
# It also provide the opportunity to write some of commonly used population genetics file formats.

Ho_Fst_down_sampling_v2 <- function(data, genepop=T, cluster=T){
  cat("\n#-- Calculating Ho and Fst per locus --#\n")
  div <- basic.stats(data)
  bi <- div$perloc[,c(1,8)]
  cat("\n#-- Sub-sampling parameters --#\n")
  cat(paste("\nRange of Ho:\n"))
  print(range(bi[,1]))
  cat("\nRange of Fst:\n")
  print(range(bi[,2]))
  if(cluster==T){
    cat("\n")
    cat("What value would you like to use as minimum (Ho)? "); minH <- as.numeric(readLines("stdin",n=1))
    cat("What value would you like to use as maximum (Ho)? "); maxH <- as.numeric(readLines("stdin",n=1))
    cat("What interval length would you like to use (Ho)? "); binH <- as.numeric(readLines("stdin",n=1))
    cat("What value would you like to use as minimum (Fst)? "); minF <- as.numeric(readLines("stdin",n=1))
    cat("What value would you like to use as maximum (Fst)? "); maxF <- as.numeric(readLines("stdin",n=1))
    cat("What interval length would you like to use (Fst)? "); binF <- as.numeric(readLines("stdin",n=1))
  }
  if(cluster==F){
    cat("\n")
    minH <- as.numeric(readline("What value would you like to use as minimum (Ho)? "))
    maxH <- as.numeric(readline("What value would you like to use as maximum (Ho)? "))
    binH <- as.numeric(readline("What interval length would you like to use (Ho)? "))
    minF <- as.numeric(readline("What value would you like to use as minimum (Fst)? "))
    maxF <- as.numeric(readline("What value would you like to use as maximum (Fst)? "))
    binF <- as.numeric(readline("What interval length would you like to use (Fst)? "))
  }
  write(paste("Ho - minH: ", minH, ", maxH:", maxH, " and binH: ", binH,
              "\nFst - minF: ", minF, ", maxF:", maxF, " and binF: ", binF,sep = ""), file = "Ho_Fst_down_sampling.log", append = F)
  intervalH <- seq(minH, maxH, binH)
  intervalF <- seq(minF, maxF, binF)
  bins = list()
  c <- 1
  cat("\n#-- Sub-sampling of bins --#\n")
  for(i in 1:(length(intervalH)-1)){
    miH <- intervalH[i]
    maH <- intervalH[i+1]
    for(j in 1:(length(intervalF)-1)){
      miF <- intervalF[j]
      maF <- intervalF[j+1]
      sel <- which(bi[,1]>=miH & bi[,1]<maH & bi[,2]>=miF & bi[,2]<maF)
      sel.names <- rownames(bi[sel,])
      bins[[c]] <- sel.names
      names(bins)[c] <- paste("binH_", miH,"_", maH, "_binF_", miF, "_", maF, sep = "")
      c <- c+1
    }
  }
  len = NULL
  for(i in 1:length(bins)){
    len <- append(len, length(bins[[i]]))
  }
  tot <- sum(len)
  perc <- len/tot
  loci.sub = list()
  for(i in 1:length(bins)){
    temp <- bins[[i]]
    nsamp <- sample(temp, size = ceiling((length(temp)*perc[i])), replace = F)
    loci.sub[[i]] <- nsamp
    names(loci.sub)[i] <- names(bins)[i]
  }
  cat("\n#-- Convert genind into data frame --#\n")
  cat("[Converting missing genotype to -9,-9]\n")
  df <- genind2df(data, usepop = T, oneColPerAll = F, sep=",")
  for(j in 2:ncol(df)){
    sel <- which(df[,j]=="NA")
    df[sel,j] <- "-9,-9"
    sel <- which(is.na(df[,j])==T)
    df[sel,j] <- "-9,-9"
    progress(j, max.value = ncol(df))
  }
  loci.list = NULL
  for(i in 1:length(loci.sub)){
    temp <- loci.sub[[i]]
    if(length(temp)!=0){loci.list <- append(loci.list, temp)}
  }
  loci.list <- sort(loci.list)
  write.table(loci.list, file = "loci_list.txt", row.names = F, quote = F, sep = "\n", col.names = F)
  cat("\n\n[Creating a new down-sampled data frame]\n")
  for(i in 1:length(loci.list)){
    pos <- which(colnames(df)==loci.list[i])
    if(i==1){dfnew <- df[,c(1,pos)]}
    if(i!=1){dfnew <- cbind(dfnew, df[,pos])}
    progress(i, max.value = length(loci.list))
  }
  cat("\n\n#-- Create a new genind object with sub-sampled loci --#\n")
  indname <- rownames(dfnew)
  popname <- dfnew[,1]
  dfnew <- dfnew[,-1]
  locname <- loci.list
  genindobj <- df2genind(dfnew, sep = ",", ncode = NULL, ind.names = indname,
                         loc.names = locname, pop = popname, NA.char = "-9", ploidy = 2,
                         type = "codom")
  divsub <- basic.stats(genindobj)
  bisub <- divsub$perloc[,c(1,8)]
  write.table(bi, file = "Ho_Fst_ori.txt", sep="\t", quote = F)
  write.table(bisub, file = "Ho_Fst_down.txt", sep="\t", quote = F)
  if(genepop==T){
    cat("\n#-- Converting genind object in other formats --#\n")
    tidy <- tidy_genind(data = genindobj, keep.allele.names = F, gds = F, verbose = T)
  }
  if(genepop==T){write_genepop(data=tidy, genepop.header = "Subsampling", filename = "output_subsampling")}
  return(genindobj)
  cat("\n#-- DONE --#")
}

#-----------------------#
# get_vcf_info function #
#-----------------------#

# This function access the given VCF file and extract chromosomes ID as well as the position
# of variants on that chromosome.

get_vcf_info <- function(vcf){
  gt <- extract.gt(vcf)
  chrom <- getCHROM(vcf)
  pos <- getPOS(vcf)
  locpositions <- matrix(ncol=3, nrow=nrow(gt))
  colnames(locpositions) <- c("locname", "CHROM","POS")
  locpositions[,1] <- rownames(gt)
  locpositions[,2] <- chrom
  locpositions[,3] <- pos
  locpositions <- data.frame(locpositions, stringsAsFactors = F)
  return(locpositions)
}

#---------------------#
# HWE_filter function #
#---------------------#

# This function perform a HWE test per locus and population and filters out markers with a pvalue
# below a user define threshold

HWE_filter <- function(data, th=0.05, cutoff=0.25){
  cat("\n#-- Performing HWE test --#\n")
  sepop <- seppop(data)
  HWE=list()
  for(i in 1:length(sepop)){
    temp <- sepop[[i]]
    HW <- hw.test(temp, B=0)
    HWE[[i]] <- HW
    names(HWE)[i] <- names(sepop)[i]
  }
  PopPr <- matrix(ncol=length(HWE), nrow = nrow(HWE[[1]]))
  rownames(PopPr) <- rownames(HWE[[1]])
  colnames(PopPr) <- names(HWE)
  cat("\n#-- Identifying loci failing HWE test --#\n")
  for(i in 1:length(HWE)){
    temp <- HWE[[i]]
    for(j in 1:nrow(temp)){
      if(temp[j,3]<th){temp[j,3] <- 1} else {temp[j,3] <- 0}
    }
    HWE[[i]] <- temp
  }
  for(i in 1:length(HWE)){
    temp <- HWE[[i]]
    PopPr[,i] <- temp[,3]
  }
  PopPr <- data.frame(PopPr, stringsAsFactors = F)
  for(i in 1:nrow(PopPr)){
    PopPr$Prop[i] <- sum(PopPr[i, 1:length(HWE)])/length(HWE) 
    if(PopPr$Prop[i]>cutoff){PopPr$fail[i] <- 1} else {PopPr$fail[i] <- 0}
  }
  fail <- subset(PopPr, PopPr$fail==1)
  locname <- rownames(fail)
  df <- genind2df(data, sep = "/", usepop = T, oneColPerAll = F)
  for(j in 2:ncol(df)){
    sel <- which(df[,j]=="NA")
    df[sel,j] <- "-9/-9"
    sel <- which(is.na(df[,j])==T)
    df[sel,j] <- "-9/-9"
  }
  loci=NULL
  for(i in 1:length(locname)){
    loci <- append(loci, which(colnames(df)==locname[i]))
  }
  df <- df[,-loci]
  indname <- rownames(df)
  popname <- df$pop
  locname <- colnames(df)[2:ncol(df)]
  df <- df[,-1]
  cat("\n#-- Creating new genind object --#\n")
  genindobj<- df2genind(df, sep = "/", ncode = NULL, ind.names = indname, pop=popname,
                                 loc.names = locname, NA.char = "-9", ploidy = 2,
                                 type = "codom")
  cat(paste("\n#-- Kept", ncol(genindobj@tab)/2, "of", ncol(data@tab)/2, "loci --#\n"))
    return(genindobj)
}

#---------------------#
# Fis_filter function #
#---------------------#

# This function filters out loci that do not fall within the specified interval.

Fis_filter <- function(data, min=-0.5, max=0.5){
  cat("\n#-- Extracting Fis values for each loci --#\n")
  stats <- basic.stats(data)
  Fis <- stats$perloc$Fis
  keep <- which(Fis>=min & Fis<=max)
  temp <- genind2df(data, sep=",", usepop = T, oneColPerAll = F)
  cat("\n#-- Keeping loci with Fis values within the specified interval --#\n")
  for(j in 2:ncol(temp)){
    sel <- which(temp[,j]=="NA")
    temp[sel,j] <- "-9,-9"
    sel <- which(is.na(temp[,j])==T)
    temp[sel,j] <- "-9,-9"
  }
  keep <- keep+1
  indname <- rownames(temp)
  popname <- temp$pop
  temp <- temp[,keep]
  locname <- colnames(temp)
  cat("\n#-- Creating a new genind object --#\n")
  genindobj <- df2genind(temp, sep = ",", ncode = NULL, ind.names = indname,
                                loc.names = locname, pop = popname, NA.char = "-9", ploidy = 2,
                                type = "codom")
  cat(paste("\n#-- Kept", ncol(genindobj@tab)/2, "of", ncol(data@tab)/2, "loci --#\n"))
  return(genindobj)
}

#----------------#
# LD_Ne function #
#----------------#

# This function implements the Ne (calculated following Waples et al. 2016 that corrects for 
# physical linkage among loci) estimation  described in Braasch et al. 2018. It samples n loci 
# of the total number of loci and estimate Ne. This process is repeated over multiple iterations.

LD_Ne <- function(data, maf=0.05, n=20, it=30000, digits=0){
  Nesim = list()
  for(p in 1:nPop(data)){
    Nesim[[p]] <- numeric(length=it)
  }
  cat("\n#-- Converting genind object into data frame --#\n")
  df <- genind2df(data, sep = "/", usepop = T, oneColPerAll = F)
  for(j in 2:ncol(df)){
    sel <- which(df[,j]=="NA")
    df[sel,j] <- "-9/-9"
    sel <- which(is.na(df[,j])==T)
    df[sel,j] <- "-9/-9"
  }
  rand <- 2:ncol(df)
  cat(paste("\n#-- Estimating Ne for", it, "iterations --#\n\n"))
  for(i in 1:it){
    cat(paste("#-- Now running iteration:", i, "--#\n"))
    samp <- sample(rand, n, replace = F)
    samp <- sort(samp)
    temp <- cbind(df$pop, df[,samp])
    indname <- rownames(temp)
    popname <- temp[,1]
    locname <- colnames(temp)[2:ncol(temp)]
    temp <- temp[,-1]
    gen <- df2genind(temp, sep = "/", ncode = NULL, ind.names = indname,  pop=popname,
                     loc.names = locname, NA.char = "-9", ploidy = 2,
                     type = "codom")
    gtypes <- genind2gtypes(gen)
    eff <- ldNe(gtypes, maf.threshold = maf, by.strata = F, ci = 0.95)
    eff <- data.frame(eff, stringsAsFactors = F)
    for(h in 1:nrow(eff)){
      Nesim[[h]][i] <- round(eff$Ne[h], digits = digits)
    }
  }
  cat("\n#-- Finilizing the results --#\n")
  for(o in 1:nrow(eff)){
    names(Nesim)[o] <- rownames(eff)[o]
  }
  return(Nesim)
}

#----------------------------#
# genind2coancestry function #
#----------------------------#

# This function takes a genind object, split it per populations and write files that 
# can be used with the software Coancestry.

genind2coancestry <- function(data){
  data <- seppop(data)
  for(i in 1:length(data)){
    cat(paste("\n #-----Now processing: ", names(data)[i], "-----#\n", sep = ""))
    df <- genind2df(data[[i]], usepop = T, oneColPerAll = T)
    for(j in 2:ncol(df)){
      sel <- which(df[,j]==0)
      df[sel,j] <- 2
    }
    for(j in 2:ncol(df)){
      sel <- which(df[,j]=="NA")
      df[sel,j] <- 0
    }
    df$pop <- rownames(df)
    rownames(df) = NULL
    write.table(df, file = paste(names(data)[i], ".txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
  }
}

#-----------------------------#
# process_coancestry function #
#-----------------------------#

# This function takes an R object return by the function coancestry 
# (package 'related') and provide estimate of relatedness and bootstrap 95% CIs.

process_coancestry <- function(data, method, digits=3){
  df <- data$relatedness
  pos = NULL
  pos <- which(names(df)==method)
  df2 <- cbind(df[,2:3], df[,pos])
  mu <- mean(df2[,3])
  oneboot <- one.boot(df2[,3], mean, R=2000)
  ci <- boot.ci(oneboot, conf=0.95, type = "basic")
  res <- matrix(ncol = 3, nrow=1)
  colnames(res) <- c("mean_relatedness", "relatedness_low", "relatedness_high")
  res[,1] <- round(mu, digits = digits)
  res[,2] <- round(ci$basic[4], digits = digits)
  res[,3] <- round(ci$basic[5], digits = digits)
  res <- data.frame(res)
  return(res)
}