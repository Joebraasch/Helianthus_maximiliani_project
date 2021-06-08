# Helianthus_maximiliani_project

This repository includes data files and R code used to conduct the analyses in the manuscript “Testing for evolutionary change in restoration: a genomic comparison between ex situ, native and commercial seed sources of Helianthus maximiliani”


Data files:

SNP_FINAL_noLD.vcf:  The vcf file of all SNPs after filtering used in the analysis

SNP_FINAL_genind.RData: A genind file of all SNPs used for analyses in R

PopulationList.csv: A list of all H. maximiliani populations linking their ID in the vcf file with the IDs presented in the paper. This file also includes the Lat and Long coordinates for populations when available

LocationDataForR.csv: A file with Lat-Long coordinates. Only includes populations for which we had location data

Genetic_Diversity_Estimats.xlsx: Spreadsheet with all mean expected heterozygosity, mean inbreeding coefficients, effective population size, coancestry coefficient, and Tajima's D for all populations

R Code:

vcf2sfs.R: Functions used to produce the site frecuency spectrum with the file SFS2.R

SFS2.R: Code that produces a site frequency spectrum. References functions in the file vcf2sfs.R

He_Fis_estimates_Linear_Models: Code that uses the genind file to estimate expected heterozygosity and inbreeding coefficients for each locus in each population. Then uses these data to perform linear mixxed models to compare seed source types

Functions.R: Functions used to produce the genind file, clean the data, estimate Ne, and estimate coancestry

Fst_code.R: Code used to estimate pairwise Fst and perform analyses associated with these estimates

DAPC_and_PCA.R: Code used to conduct PCA and DAPC analyses using the genind file
