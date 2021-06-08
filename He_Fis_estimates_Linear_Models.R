######## make linear models to test for differences btwn populations and sources


library(adegenet)
library(lme4)
library(multcomp)
library(ggplot2)
library(lmerTest)

#Load data
load("SNP_FINAL_genind.RData")
HM_genepop <- genindHM12bi

#Separate genind by population
Each_pop <- seppop(HM_genepop)

##############
#   Estimates of He and Fis
##############
#### NATIVE POPS
# Calculate for each population
A_summary <- summary(Each_pop$A)
A_He <- A_summary$Hexp
A_Fis <- (A_summary$Hexp - A_summary$Hobs)/A_summary$Hexp
A_list <- rep("A", length(A_He))


B_summary <- summary(Each_pop$B)
B_He <- B_summary$Hexp
B_Fis <- (B_summary$Hexp - B_summary$Hobs)/B_summary$Hexp
B_list <- rep("B", length(B_He))

C_summary <- summary(Each_pop$C)
C_He <- C_summary$Hexp
C_Fis <- (C_summary$Hexp - C_summary$Hobs)/C_summary$Hexp
C_list <- rep("C", length(C_He))

D_summary <- summary(Each_pop$D)
D_He <- D_summary$Hexp
D_Fis <- (D_summary$Hexp - D_summary$Hobs)/D_summary$Hexp
D_list <- rep("D", length(D_He))

E_summary <- summary(Each_pop$E)
E_He <- E_summary$Hexp
E_Fis <- (E_summary$Hexp - E_summary$Hobs)/E_summary$Hexp
E_list <- rep("E", length(E_He))

F_summary <- summary(Each_pop$F)
F_He <- F_summary$Hexp
F_Fis <- (F_summary$Hexp - F_summary$Hobs)/F_summary$Hexp
F_list <- rep("F", length(F_He))

#Combine estimates for each population
Nat_He <- c(A_He, B_He, C_He, D_He, E_He, F_He)
Nat_Fis <- c(A_Fis, B_Fis, C_Fis, D_Fis, E_Fis, F_Fis)
Nat_pops <- c(A_list, B_list, C_list, D_list, E_list, F_list)
Nat_list <- rep("Native", length(Nat_pops))

#Create data frame with He and Fis
Nat_DF <- data.frame(Nat_He, Nat_Fis, Nat_pops, Nat_list)
colnames(Nat_DF) <- c("He", "Fis", "Pop", "Source")

### Commercial
AG_summary <- summary(Each_pop$AG)
AG_He <- AG_summary$Hexp
AG_Fis <- (AG_summary$Hexp - AG_summary$Hobs)/AG_summary$Hexp
AG_list <- rep("AG", length(AG_He))

MSG_summary <- summary(Each_pop$MSG)
MSG_He <- MSG_summary$Hexp
MSG_Fis <- (MSG_summary$Hexp - MSG_summary$Hobs)/MSG_summary$Hexp
MSG_list <- rep("MSG", length(MSG_He))

PrairieRest_summary <- summary(Each_pop$PrairieRest)
PrairieRest_He <- PrairieRest_summary$Hexp
PrairieRest_Fis <- (PrairieRest_summary$Hexp - PrairieRest_summary$Hobs)/PrairieRest_summary$Hexp
PrairieRest_list <- rep("PrairieRest", length(PrairieRest_He))

RG_summary <- summary(Each_pop$RG0426171)
RG_He <- RG_summary$Hexp
RG_Fis <- (RG_summary$Hexp - RG_summary$Hobs)/RG_summary$Hexp
RG_list <- rep("RG", length(RG_He))

RMAX_summary <- summary(Each_pop$RMAX17)
RMAX_He <- RMAX_summary$Hexp
RMAX_Fis <- (RMAX_summary$Hexp - RMAX_summary$Hobs)/RMAX_summary$Hexp
RMAX_list <- rep("RMAX", length(RMAX_He))

Com_He <- c(AG_He, MSG_He, PrairieRest_He, RG_He, RMAX_He)
Com_Fis <- c(AG_Fis, MSG_Fis, PrairieRest_Fis, RG_Fis, RMAX_Fis)
Com_pops <- c(AG_list, MSG_list, PrairieRest_list, RG_list, RMAX_list)
Com_list <- rep("Commercial", length(Com_He))

Com_DF <- data.frame(Com_He, Com_Fis, Com_pops, Com_list)
colnames(Com_DF) <- c("He", "Fis", "Pop", "Source")

#Ex situ pops  

ESB_summary <- summary(Each_pop$PI586891)
ESB_He <- ESB_summary$Hexp
ESB_Fis <- (ESB_summary$Hexp - ESB_summary$Hobs)/ESB_summary$Hexp
ESB_list <- rep("ES-B", length(ESB_He))

ESC_summary <- summary(Each_pop$PI586892)
ESC_He <- ESC_summary$Hexp
ESC_Fis <- (ESC_summary$Hexp - ESC_summary$Hobs)/ESC_summary$Hexp
ESC_list <- rep("ESC", length(ESC_He))

ESA_summary <- summary(Each_pop$PI586893)
ESA_He <- ESA_summary$Hexp
ESA_Fis <- (ESA_summary$Hexp - ESA_summary$Hobs)/ESA_summary$Hexp
ESA_list <- rep("ESA", length(ESA_He))

ESD_summary <- summary(Each_pop$PI586904)
ESD_He <- ESD_summary$Hexp
ESD_Fis <- (ESD_summary$Hexp - ESD_summary$Hobs)/ESD_summary$Hexp
ESD_list <- rep("ESD", length(ESD_He))

ESF_summary <- summary(Each_pop$PI650010)
ESF_He <- ESF_summary$Hexp
ESF_Fis <- (ESF_summary$Hexp - ESF_summary$Hobs)/ESF_summary$Hexp
ESF_list <- rep("ESF", length(ESF_He))

ESE_summary <- summary(Each_pop$PI650011)
ESE_He <- ESE_summary$Hexp
ESE_Fis <- (ESE_summary$Hexp - ESE_summary$Hobs)/ESE_summary$Hexp
ESE_list <- rep("ESE", length(ESE_He))

ES_He <- c(ESB_He, ESC_He, ESA_He, ESD_He, ESF_He, ESE_He)
ES_Fis <- c(ESB_Fis, ESC_Fis, ESA_Fis, ESD_Fis, ESF_Fis, ESE_Fis)
ES_pops <- c(ESB_list, ESC_list, ESA_list, ESD_list, ESF_list, ESE_list)
ES_list <- rep("Ex_situ", length(ES_He))

ES_DF <- data.frame(ES_He, ES_Fis, ES_pops, ES_list)
colnames(ES_DF) <- c("He", "Fis", "Pop", "Source")

# Selected pops  
TLi_summary <- summary(Each_pop$TLi)
TLi_He <- TLi_summary$Hexp
TLi_Fis <- (TLi_summary$Hexp - TLi_summary$Hobs)/TLi_summary$Hexp
TLi_list <- rep("TLi", length(TLi_He))

TLiM_summary <- summary(Each_pop$TLiM)
TLiM_He <- TLiM_summary$Hexp
TLiM_Fis <- (TLiM_summary$Hexp - TLiM_summary$Hobs)/TLiM_summary$Hexp
TLiM_list <- rep("TLiM", length(TLiM_He))

Sel_He <- c(TLi_He, TLiM_He)
Sel_Fis <- c(TLi_Fis, TLiM_Fis)
Sel_pops <- c(TLi_list, TLiM_list)
Sel_list <- rep("Selected", length(Sel_pops))

Sel_DF <- data.frame(Sel_He, Sel_Fis, Sel_pops, Sel_list)
colnames(Sel_DF) <- c("He", "Fis", "Pop", "Source")


## Combine all seed source types into a single dataframe
He_DF <- rbind(Nat_DF, Com_DF, ES_DF, Sel_DF)


######  Linear model

#He.full.LM <- lmer(He ~ Source + (1|Pop), data = He_DF)
Fis.full.LM <- lmer(Fis ~ Source + (1|Pop), data = He_DF)
#print(He.full.LM)

print(c("THIS IS THE FULL MODEL SUMMARY"))
#print(summary(He.full.LM))
print(summary(Fis.full.LM))
summary(Fis.full.LM)

print(c("BREAK BREAK BREAK BREAK"))
print(c("BREAK BREAK BREAK BREAK"))
print(c("THIS IS THE ANOVA SUMMARY FOR THE FULL MODEL"))

#print(anova(He.full.LM))
#print(summary(anova(He.full.LM)))

print(anova(Fis.full.LM))
print(summary(anova(Fis.full.LM)))
anova(Fis.full.LM)
summary(anova(Fis.full.LM))

print(c("BREAK BREAK BREAK BREAK"))
print(c("BREAK BREAK BREAK BREAK"))
print(c("THIS IS THE MODEL WITHOUT SOURCE TYPE /n"))

#He.NoSource.LM <- lmer(He ~ (1|Pop), data = He_DF)
#print(summary(He.NoSource.LM))

Fis.NoSource.LM <- lmer(Fis ~ (1|Pop), data = He_DF)
print(summary(Fis.NoSource.LM))

print(c("BREAK BREAK BREAK BREAK"))
print(c("BREAK BREAK BREAK BREAK"))
print(c("THIS IS THE COMPARISON OF THE FULL MODEL AND NO SOURCE MODEL"))

#print(anova(He.NoSource.LM, He.full.LM))
print(anova(Fis.NoSource.LM, Fis.full.LM))

print(c("BREAK BREAK BREAK BREAK"))
print(c("BREAK BREAK BREAK BREAK"))
print(c("THIS IS THE POST HOC TEST FOR THE FULL MODEL"))



post.hoc <- summary(glht(Fis.full.LM , linfct = mcp(Source = "Tukey")))
print(post.hoc)
post.hoc





### Load GenDiv table
df <- readxl::read_xlsx("GenDiv.xlsx"); df <- data.frame(df)

###Simple linear models for Ne and Theta
df <- df[-which(df$Group %in% c("Selected")),] # remove selected populatoins because of low replication

Ne_LM <- aov(Ne ~ Group, data = df)
summary(Ne_LM)


Theta_LM <- aov(r ~ Group, data = df)
summary(Theta_LM)
TukeyHSD(Theta_LM)

