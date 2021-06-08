#################################################
# Discriminant Analysis of Principal Components #
#################################################
setwd("~/Research/Helianthus_work")
###
# Aim
###

# Perform a DAPC using SNPs calles in Helianthis maximilinani

install.packages("adegenet", dependencies = TRUE)
install.packages("adegraphics", dependencies = TRUE)


##########################
# Load packages and data #
##########################
library(adegenet)
library(adegraphics)
library(ggplot2)
load("SNP_FINAL_genind.RData")
strata <- read.table("Strata.txt", sep="\t", h=T, stringsAsFactors = T)

###################
# All populations #
###################
## DAPC cross validation
#Replace NAs for DAPC
X <- scaleGen(genindHM12bi, NA.method="mean")


#Do cross validation. Initial centered best fits around 200 PCs
set.seed(999)
#crossval <- xvalDapc( X,genindHM12bi$pop,  n.pca = 150:250, n.rep = 300)
# 159 PCAs had highest mean successful assignment (0.8714) and lowest MSE (0.1334)
crossval <- xvalDapc( X,genindHM12bi$pop,  n.pca = 100:290, n.rep = 40)
# 175 PCAs had highest mean successful assignment (0.8776316) and lowest MSE (0.1270265)

#-- Perform DAPC --#
dapc1 <- dapc(genindHM12bi, genindHM12bi$pop, perc.pca = 90, n.da = 18)
round(dapc1$eig/sum(dapc1$eig), digits=3)

#-- Colour by groups --#
Each_pop <- seppop(genindHM12bi)

invs.list.TLi <- rownames(Each_pop$TLi@tab)
invs.list.TLiM <- rownames(Each_pop$TLiM@tab)
invs.list <- c(invs.list.TLi, invs.list.TLiM)

pops.TLi<- rep("TLi", length(invs.list.TLi))
pops.TliM <- rep("TLiM", length(invs.list.TLiM))
pops.list <- c(pops.TLi, pops.TliM)

type.TLi <- rep("Selected", length(invs.list.TLi))
type.TLiM <- rep("Selected", length(invs.list.TLiM))
type.list <- c(type.TLi, type.TLiM)

temp.frame <- data.frame(invs.list, pops.list, type.list)
colnames(temp.frame) <- colnames(strata)

strata <- rbind(strata, temp.frame)


colour=NULL
for(i in 1:nrow(strata)){
  if(strata$POPULATIONS[i]=="A" | strata$POPULATIONS[i]=="B" | strata$POPULATIONS[i]=="C" | 
     strata$POPULATIONS[i]=="D" | strata$POPULATIONS[i]=="E" | strata$POPULATIONS[i]=="F") {strata$COL[i] <- "forestgreen"}
  if(strata$POPULATIONS[i]=="RG0426171" | strata$POPULATIONS[i]=="RMAX17" | strata$POPULATIONS[i]=="MSG" | 
     strata$POPULATIONS[i]=="AG" | strata$POPULATIONS[i]=="PrairieRest"){strata$COL[i] <- "darkcyan"}
  if(strata$POPULATIONS[i]=="PI586893" | strata$POPULATIONS[i]=="PI586891" | strata$POPULATIONS[i]=="PI586892" |
     strata$POPULATIONS[i]=="PI586904" | strata$POPULATIONS[i]=="PI650011" | strata$POPULATIONS[i]=="PI650010"){strata$COL[i] <- "gold"}
  if(strata$POPULATIONS[i]=="TLiM" | strata$POPULATIONS[i]=="TLi"){strata$COL[i] <- "orange"}
  #if(strata$POPULATIONS[i]=="AG" | strata$POPULATIONS[i]=="PrairieRest"){strata$COL[i] <- "blue"}
}





# Graph the inside of the figure
s.class(dapc1$ind.coord, fac = pop(genindHM12bi),col=alpha(strata$COL[which(duplicated(strata$POPULATIONS)==F)], 0.6),
        paxes.draw=F, plabels.cex=0, plegend.drawKey=F, pgrid.draw=F, ppoint.cex=1, porigin.lty=2,
        pellipses.lwd=1, ellipseSize = 3.5, starSize = 1, ppoints.pch=21, pellipses.border="black",
        ppoints.col=alpha("black", 0.6))

#Produce frame for the figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=F)
plot(dapc1$ind.coord[,1], dapc1$ind.coord[,2], col=alpha(strata$COL, 0.6), pch=20, 
     cex=2, xlab="LD 1 (32.2% var. explained)", ylab="LD 2 (18.1% var. explained)", font.lab=2, main="Discriminant Analysis of Principal Components\n Helinathus maximiliani")
points(dapc1$ind.coord[,1], dapc1$ind.coord[,2], cex=1.45, col=alpha("black", 0.15))
abline(h=0, lty=2)
abline(v=0, lty=2)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=T)
text(dapc1$grp.coord[,1], dapc1$grp.coord[,2], levels(dapc1$grp))
legend("topright", inset=c(-0.3,0), legend = c("Ex situ", "Contemporary", "Commercial", "Selected"), pch=20, lty=1, pt.cex=3, cex=0.95, col = c("gold", "forestgreen", "darkcyan", "orange"), bty="n")


#################################
# Ex situ and Contemporary only #
#################################

#-- Create new genind object with data associated with ex situ and contemporary populations only --#
genindHM12bi_excont <- genindHM12bi[-which(genindHM12bi@pop=="TLi" | 
                                             genindHM12bi@pop=="TLiM" | 
                                             genindHM12bi@pop=="RG0426171" | 
                                             genindHM12bi@pop=="RMAX17" | 
                                             genindHM12bi@pop=="MSG" | 
                                             genindHM12bi@pop=="AG" | 
                                             genindHM12bi@pop=="PrairieRest")]



#Replace NAs for DAPC
X <- scaleGen(genindHM12bi_excont , NA.method="mean")


#Do cross validation. Initial centered best fits around 200 PCs
set.seed(999)
crossval <- xvalDapc( X,genindHM12bi_excont $pop,  n.pca = 100:290, n.rep = 40)
# 108 PCAs had highest mean successful assignment (0.9906250) and lowest MSE (0.02185018)

#-- Perform DAPC --#
#dapc1 <- dapc(genindHM12bi_excont, genindHM12bi_excont@pop, perc.pca = 90, n.da = 11); dapc1 
dapc1 <- dapc(genindHM12bi_excont, genindHM12bi_excont@pop, npca = 108, n.da = 11); dapc1 

round(dapc1$eig/sum(dapc1$eig), digits=3)

#-- Colour by groups --#
strata <- read.table("Strata.txt", sep="\t", h=T, stringsAsFactors = T)
sel <- which(strata[,3]=="Producer" | strata[,2]=="TLiM" | strata[,2]=="TLiM")
strata <- strata[-sel,]
colour=NULL
for(i in 1:nrow(strata)){
  if(strata$POPULATIONS[i]=="A" | strata$POPULATIONS[i]=="B" | strata$POPULATIONS[i]=="C" | 
     strata$POPULATIONS[i]=="D" | strata$POPULATIONS[i]=="E" | strata$POPULATIONS[i]=="F") {strata$COL[i] <- "forestgreen"}
  if(strata$POPULATIONS[i]=="PI586893" | strata$POPULATIONS[i]=="PI586891" | strata$POPULATIONS[i]=="PI586892" |
     strata$POPULATIONS[i]=="PI586904" | strata$POPULATIONS[i]=="PI650011" | strata$POPULATIONS[i]=="PI650010"){strata$COL[i] <- "gold"}
}


#Produce graph for the figure
s.class(dapc1$ind.coord, fac = pop(genindHM12bi_excont),col=alpha(strata$COL[which(duplicated(strata$POPULATIONS)==F)], 0.6),
        paxes.draw=F, plabels.cex=0, plegend.drawKey=F, pgrid.draw=F, ppoint.cex=1, porigin.lty=2,
        pellipses.lwd=1, ellipseSize = 3.5, starSize = 1, ppoints.pch=21, pellipses.border="black",
        ppoints.col=alpha("black", 0.6))


#Produce frame for the figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=F)
plot(dapc1$ind.coord[,1], dapc1$ind.coord[,2], col=alpha(strata$COL, 0.6), pch=20, 
     cex=2, xlab="LD 1 (24.1% var. explained)", ylab="LD 2 (16.0% var. explained)", font.lab=2, main="Discriminant Analysis of Principal Components\n Helinathus maximiliani")
points(dapc1$ind.coord[,1], dapc1$ind.coord[,2], cex=1.45, col=alpha("black", 0.15))
abline(h=0, lty=2)
abline(v=0, lty=2)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=T)
text(dapc1$grp.coord[,1], dapc1$grp.coord[,2], levels(dapc1$grp))
legend("topright", inset=c(-0.3,0), legend = c("Ex situ", "Contemporary"), pch=20, lty=1, pt.cex=3, cex=0.95, col = c("gold", "forestgreen"), bty="n")



#### Black and White

strata <- read.table("Strata.txt", sep="\t", h=T, stringsAsFactors = T)
sel <- which(strata[,3]=="Producer" | strata[,2]=="TLiM" | strata[,2]=="TLiM")
strata <- strata[-sel,]
colour=NULL
for(i in 1:nrow(strata)){
  if(strata$POPULATIONS[i]=="A" | strata$POPULATIONS[i]=="B" | strata$POPULATIONS[i]=="C" | 
     strata$POPULATIONS[i]=="D" | strata$POPULATIONS[i]=="E" | strata$POPULATIONS[i]=="F") {strata$COL[i] <- "darkgrey"}
  if(strata$POPULATIONS[i]=="RG0426171" | strata$POPULATIONS[i]=="RMAX17" | strata$POPULATIONS[i]=="MSG" | 
     strata$POPULATIONS[i]=="AG" | strata$POPULATIONS[i]=="PrairieRest"){strata$COL[i] <- "darkgrey"}
  if(strata$POPULATIONS[i]=="PI586893" | strata$POPULATIONS[i]=="PI586891" | strata$POPULATIONS[i]=="PI586892" |
     strata$POPULATIONS[i]=="PI586904" | strata$POPULATIONS[i]=="PI650011" | strata$POPULATIONS[i]=="PI650010"){strata$COL[i] <- "darkgrey"}
  if(strata$POPULATIONS[i]=="TLiM" | strata$POPULATIONS[i]=="TLi"){strata$COL[i] <- "darkgrey"}
}

for(i in 1:nrow(strata)){
  if(strata$POPULATIONS[i]=="A" | strata$POPULATIONS[i]=="B" | strata$POPULATIONS[i]=="C" | 
     strata$POPULATIONS[i]=="D" | strata$POPULATIONS[i]=="E" | strata$POPULATIONS[i]=="F") {strata$POINT[i] <- 22}
  if(strata$POPULATIONS[i]=="RG0426171" | strata$POPULATIONS[i]=="RMAX17" | strata$POPULATIONS[i]=="MSG" | 
     strata$POPULATIONS[i]=="AG" | strata$POPULATIONS[i]=="PrairieRest"){strata$POINT[i] <- 24}
  if(strata$POPULATIONS[i]=="PI586893" | strata$POPULATIONS[i]=="PI586891" | strata$POPULATIONS[i]=="PI586892" |
     strata$POPULATIONS[i]=="PI586904" | strata$POPULATIONS[i]=="PI650011" | strata$POPULATIONS[i]=="PI650010"){strata$POINT[i] <- 21}
  if(strata$POPULATIONS[i]=="TLiM" | strata$POPULATIONS[i]=="TLi"){strata$POINT[i] <- 23}
}



par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=F)
plot(dapc1$ind.coord[,1], dapc1$ind.coord[,2], col=alpha(strata$COL, 0.5), pch=strata$POINT, 
     cex=2, xlab="LD 1 (32.3% var. explained)", ylab="LD 2 (15.3% var. explained)", font.lab=2, main="Discriminant Analysis of Principal Components\n Helinathus maximiliani")
points(dapc1$ind.coord[,1], dapc1$ind.coord[,2], cex=2, col=alpha("black", 0.1), pch=strata$POINT)
abline(h=0, lty=2)
abline(v=0, lty=2)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=T)
text(dapc1$grp.coord[,1], dapc1$grp.coord[,2], levels(dapc1$grp))
legend("topright", inset=c(-0.2,0), legend = c("Ex situ", "Native"), pch=c(21,22), lty=1, pt.cex=2, cex=0.95, col = c("black"), bty="n")







########################################
# PCA of all loci           
########################################


load("SNP_FINAL_genind.RData")
#HM_genepop <- read.genepop("Test_data_3pops_10loci.gen", ncode = 3)
HM_genepop <- genindHM12bi

SNP_Frame <- HM_genepop@tab

sum(is.na(SNP_Frame))   # 226294 missing loci

#Replace missing data with average estimates and then perform PCA
PCA_in <- tab(HM_genepop, freq=TRUE, NA.method="mean")
HM_PCA <- dudi.pca(PCA_in, center=TRUE, scale=TRUE, scannf = FALSE, nf= ncol(PCA_in))

#### Calculate variance explained by PC1 and PC2
HM_PCA$eig[1]/sum(HM_PCA$eig)    #PC1 4.0 %
HM_PCA$eig[2]/sum(HM_PCA$eig)    #PC2 2.2 %

s.class(HM_PCA$li, fac=pop(HM_genepop),col=alpha(strata$COL[which(duplicated(strata$POPULATIONS)==F)], 0.6),axesel=FALSE, cstar=0, cpoint=3)


s.class(HM_PCA$li, fac = pop(genindHM12bi),col=alpha(strata$COL[which(duplicated(strata$POPULATIONS)==F)], 0.6),
        paxes.draw=F, plabels.cex=1, plegend.drawKey=F, pgrid.draw=F, ppoint.cex=1, porigin.lty=2,
        xlab="LD 1 (4.0% var. explained)", ylab="LD 2 (2.2% var. explained)",
        xlim = c(-150, 50),
        pellipses.lwd=1, ellipseSize = 3.5, starSize = 1, ppoints.pch=21, pellipses.border="black")


## Calculate sum of eigenvalues to estimate variation explained by each axes and cummulative axes

All_eig <- sum(HM_PCA$eig)  #25874
Half_var <- All_eig/2       #12937

sum(HM_PCA$eig[1:110])   # 12961 - just over 50% of the variation
