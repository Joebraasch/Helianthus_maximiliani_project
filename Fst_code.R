
###################################
### Fst additional calculations
###################################
#install.packages("geodist", dependencies = TRUE)

library(geodist)
library(lme4)
library(lmerTest)
library(multcomp)
library(lmtest)
library(ComplexHeatmap)
library(ggplot2)

#load files
Pop_info <- read.table(file = "PopulationList.csv", sep = ",", header = TRUE)
Fst_table <- read.table(file = "Ordered_WC_Fst.csv", sep = ",", header = TRUE)
colnames(Fst_table) <- gsub("\\.", "-", colnames(Fst_table))
Fst_table <- Fst_table [,-1]
popnames <- colnames(Fst_table)
popnames <- popnames[-1]
popnames <- c(popnames, "S-2")
rownames(Fst_table) <- popnames
Fst_table <- round(Fst_table, 4)


#### Make Fst Heatmap
Heatmap(as.matrix(Fst_table), cluster_rows = FALSE, cluster_columns = FALSE, 
        na_col = "white", rect_gp = gpar(col = "white", lwd = 2), 
        heatmap_legend_param = list(
          title = "Fst", legend_height = unit(6, "cm")  
        ))



###Load in the table with pairwise Fst, population information, and gps coordinates
Fst_dists <- read.table(file = "Fst_pairwise_for_LM.csv", sep = ",", header = TRUE)


## First test if intra and inter population comparisons are sig diff
Inter <- subset(Fst_dists, Inter_Intra == "Inter")
Fst_dists$Inter_Intra <- gsub("Inter", "Inter-population comparisons", Fst_dists$Inter_Intra)

Intra <- subset(Fst_dists, Inter_Intra == "Intra")
Fst_dists$Inter_Intra <- gsub("Intra", "Intra-population comparisons", Fst_dists$Inter_Intra)

Inter_Intra <- wilcox.test(Inter$Fst, Intra$Fst) ## P < 0.001, W= 3845

ggplot(Fst_dists, aes(Inter_Intra,Fst)) +geom_boxplot()+
  xlab("")+
  theme_bw()+                      ###  everything after here removes background      
  theme(axis.line = element_line(colour = "black"),
        axis.text.y = element_text(color = "black", size=12),
        axis.text.x = element_text(color = "black", size=12),
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


### Separate GPS coordinates into smaller data frames to be used in distance calculations
Pop1_Gps <- Fst_dists[, c('Pop1_Long', 'Pop1_Lat')]
Pop2_Gps <- Fst_dists[, c('Pop2_Long', 'Pop2_Lat')]


### Calculate Euclidean distances using gdistance package. Store result as an ordered vector
Fst_dists$HaverDist <- geodist(Pop1_Gps, Pop2_Gps, paired = TRUE, measure = "haversine")

### Rescale distances to be km rather than m
Fst_dists$HaverDist <- Fst_dists$HaverDist/1000


### Remove comparisions with selected individuals, since we have no distance data 
Fst_dists <- subset(Fst_dists, HaverDist != "NaN")


### Remove population ES-E
#Fst_dists <- subset(Fst_dists, Pop1 != "ES-E")
#Fst_dists <- subset(Fst_dists, Pop2 != "ES-E")


### Simple linear model to compare Fst and Haverdistances
lm.distonly <- lm(Fst ~ HaverDist, data = Fst_dists)   ## Not sig - P = 0.42
lm.simple <- lm(Fst ~ HaverDist+comparison, data = Fst_dists)  ## Dist sig (0.017)
lm.cmponly <- lm(Fst ~ comparison, data = Fst_dists)  ##


post.hoc <- summary(glht(lm.simple , linfct = mcp(comparison = "Tukey")))
# WW/CC, WEs/EsC, WW/EsC, WW/EsEs, WW/WC

# Mixxed effects models to compare
modl <- lmer (Fst ~ sqrt(HaverDist) + ( sqrt(HaverDist)|comparison), data = Fst_dists) ## Mixxed Model with individual slope and intercept for each Population
mod1 <- lmer (Fst ~ sqrt(HaverDist) + (1|comparison), data = Fst_dists) ## Mixxed Model with individual slope and intercept for each Population
mod2 <- lm (Fst ~ sqrt(HaverDist), data = Fst_dists) ## Mixxed Model with individual slope and intercept for each Population
mod3 <- lmer (Fst ~  (sqrt(HaverDist)|comparison), data = Fst_dists) ## Mixxed Model with individual slope and intercept for each Population
mod4 <- lmer (Fst ~  (1|comparison), data = Fst_dists)


## Compare models with liklihood ratio test
lrtest(modl,mod1,mod2,mod3, mod4)

mm_plot <- ggplot(Fst_dists, aes(x =HaverDist, y = Fst, colour = comparison)) +
 # facet_wrap(~mountainRange, nrow=2) +   # a panel for each mountain range
  geom_point(alpha = 0.8) +
  geom_line(data = cbind(Fst_dists, pred = predict(mod3)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
  xlab("Distance (Km)")+
  theme_bw()+   
  scale_color_brewer(palette="Dark2")+ ###  everything after here removes background      
  theme(axis.line = element_line(colour = "black"),
        axis.text.y = element_text(color = "black", size=12),
        axis.text.x = element_text(color = "black", size=12),
        axis.title.y = element_text(size=14),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


values <- coef(modl)$comparison   #get intercepts and coefficients from the mixxed model
values[,3] <- rownames(values)
colnames(values)<- c("Intercept", "Slope", "comparison")


ggplot(Fst_dists, aes(HaverDist,Fst)) +geom_point()+facet_wrap(~comparison, ncol=3)+
  geom_abline(aes(intercept=Intercept, slope=Slope), data=values, colour="blue")+
  theme_bw()+                      ###  everything after here removes background      
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
