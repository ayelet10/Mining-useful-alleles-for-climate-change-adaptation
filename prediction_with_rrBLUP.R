## Prediction with rrBLUP mixed.solve ###
## output to be used for the following 2 approaches :
##  - breeding scheme of prediction only approach
##  - qtl mapping approach 
start_time <- Sys.time()
print(start_time)
library(genomicSimulation)
library(vcfR)
library(rrBLUP)
##############
## settings ##
##############
args <- commandArgs(trailingOnly = TRUE) # TRUE
print(args)
seed_num = args[1]
# seed_num = 1231227
tested_trait <- "sal_opt"  #"temp_opt" #  
# input.file_vcf <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/Case1.10krows.twoindivs.vcf"
# input.file_vcf <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/Case1.10krows.tenindivs.vcf" # 10 inds
# input.file_vcf <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/vcf_files/%s_VCF_causal.vcf", seed_num) 
# input.file_vcf <- sprintf("C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/%s_VCF_causal.vcf", seed_num)

# input.file_vcf <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_plusneut_MAF01.recode2.vcf.gz", seed_num) 
input.file_vcf <- "/home/salmanay/projects/breeding_simulations/1231227/1231227.1000_individuals_recode_vcf_out_25Jan.vcf"
# 1231227.1000_individuals_recode_vcf_out_25Jan.vcf

# input.file_vcf <- "/home/salmanay/projects/simulated_data_laurson2021/tests_on_subseted_vcf_file/Case1.10krows.twoindivs.vcf"

# my_ind_data <-"C:/Users/ayelet/Documents/UCDavis/scripts/simulations/Case1_1442299973452_ind.txt"
# my_ind_data <- "/home/salmanay/projects/simulated_data_laurson2021/Case1_1442299973452_ind.txt"
my_ind_data <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_ind.txt", seed_num)
# my_ind_data <-"C:/Users/ayelet/Documents/UCDavis/scripts/simulations/data_from_Lotterhos_preprint_2022/1231227_ind.txt"

my_pop_data <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_popInfo.txt", seed_num) # opt0 - salinity; opt1 -temp
# my_pop_data  <-  "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/data_from_Lotterhos_preprint_2022/1231227_popInfo.txt"

# my_wd <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/"
# my_wd <- "/home/salmanay/projects/breeding_simulations/genomicSimulation/"
my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/",seed_num )

dir.create(file.path(my_wd)) # if dir exists, it just gives a warning
setwd(file.path(my_wd))
###############
## read data ##
###############

vcf <- read.vcfR(input.file_vcf, verbose = FALSE )
my_geno <- extract.gt(vcf) #, as.numeric = TRUE)
# my_geno[1:10,1:10]
my_geno_no_na <- na.omit (my_geno) # exclude rows with NA       >>> to understand: why there are NAs in the simulated data??
my_geno_no_na <- my_geno_no_na     #[1:100,95:105]
# print(my_geno_no_na[1:4, 1:10]) # to see that the data was read fine
num_sites              <- nrow(my_geno_no_na)
num_individuals_tested <- ncol(my_geno_no_na)# 1000 #1000 # 
print(c("number of individuals:" ,num_individuals_tested ))
print(c("number of sites:" , num_sites))
# colnames(my_geno_no_na)[1] <- paste0("name\t", colnames(my_geno_no_na)[1])  # add a string and \t in the beginning of the colnames to fit the genomicSimulation allele file input format
my_geno_no_na[1:num_sites,1:num_individuals_tested] <- gsub("|", "", my_geno_no_na[1:num_sites,1:num_individuals_tested], fixed = TRUE)
# my_geno_no_na <- gsub("|", "", my_geno_no_na, fixed = TRUE)
# my_geno_no_na[1:2,1:2] <- gsub("|", "", my_geno_no_na[1:2,1:2], fixed = TRUE)


## read phenotype data ## >> to replace with environment data
ind_file <- read.table(file = my_ind_data, header = TRUE, stringsAsFactors = FALSE)
# pop_file <- read.table(file = my_pop_data, header = TRUE, stringsAsFactors = FALSE)


my_y <- ind_file[1:num_individuals_tested,c(tested_trait)]
#########################
## genomic prediction  ## 
#########################
## rrBLUP
# y = # observations

my_W2 <- matrix (as.numeric(t(my_geno_no_na)), ncol = num_sites, nrow =num_individuals_tested) # genomic data
# indices_ordered <- order(1:num_individuals_tested, decreasing=T)
# print(indices_ordered)
# my_W2_ordered <- my_W2[c(indices_ordered),]
# print()
# # unique(my_W2)
row.names(my_W2) <- colnames(my_geno_no_na)
# RRBLUP2 <- mixed.solve(y =my_y, Z = my_W2) # y - environment; Z - genomic data
RRBLUP2 <- mixed.solve(y =my_y, K = A.mat(my_W2)) # y - environment; K - genomic data

# additive genetic variance
RRBLUP2$Vu
# residual variance
RRBLUP2$Ve
# intercept 
RRBLUP2$beta
# additive genetic values
head(RRBLUP2$u)
tail(RRBLUP2$u)
# genomic h2
RRBLUP2$Vu / (RRBLUP2$Vu + RRBLUP2$Ve)
# ratio of variance components 
RRBLUP2$Ve / RRBLUP2$Vu
###############################
# get the selected landraces ##
###############################
# order(RRBLUP2$u,decreasing=T)[1:50]
# predicted_lr_geno <- my_W2[order(RRBLUP2$u,decreasing=T)[1:50] ] ## get the 50 best genotype
indices_ordered <- order(RRBLUP2$u,decreasing=T)  #[1:50]
# print(RRBLUP2$u)
# print(indices_ordered)

print(nrow(my_W2))
print(length(indices_ordered))
my_W2_ordered <- my_W2[indices_ordered,]  #>>>> problem here ## get the 50 best genotypes (change range to the number of desired best genotypes, 50)

ordered_u <- RRBLUP2$u[indices_ordered]
## get the id of the predicted
predicted_lr <- rownames(my_W2_ordered) ## get the ids for best genotype
# predicted_lr <- cbind(predicted_lr, ordered_u)
predicted_lr <- na.omit(predicted_lr[1:50] ) ## get the ids for best genotype
####################################
# save selected landraces to file ##
####################################
# write.csv(ordered_u, file = paste0(my_wd,"u_mixed_solve.csv") ) # one columns : ordered_u
# write.csv(predicted_lr, file = paste0(my_wd,"ids_predicted_landraces_mixed_solve.csv") , row.names = FALSE) # one columns : predicted_lr
# # write.csv(my_W2_ordered, file = paste0(my_wd,"genotypes_predicted_landraces_mixed_solve.csv"))

# a vcf file for the qtl mapping
# write.vcf(vcf[, c(1,indices_ordered)], sprintf("%s_genotypes_predicted_landraces_mixed_solve.vcf.gz", seed_num ) )
# a txt file for the breeding scheme  -- not sure it is needed because I would like to have each individual in a separate file to read into the breeding simulation ?
# write.table(my_geno_no_na[, c(1,indices_ordered)], file = "rrBLUP_predicted_allele_file_as_txt.txt", sep = "\t", quote = FALSE)

print(Sys.time())
print("before writing genotypes to txt files")

## write genotype of each predicted individual to single txt file ##
for(id_p in 1:length(predicted_lr)){
  print(predicted_lr[id_p])
  id_predicted_individual <- predicted_lr[id_p]
  ind_geno <- as.matrix(my_geno_no_na[, predicted_lr[id_p]])
  ind_geno <- cbind(ind_geno, rownames(ind_geno))
  colnames(ind_geno) <- c(id_predicted_individual,"name")
  
  source("/home/salmanay/projects/breeding_simulations/genomicSimulation/scripts/1_GEA_genomicSimulation.R") # to use the function UpdateMarkerNamesandDuplicants
  ## update markers names to be as in muts_full format (remove serial number and add number to markers in duplicated positions)
  ind_geno <- UpdateMarkerNamesandDuplicants(as.data.frame(ind_geno))  ## needs to be tested in this script 
  
  write.table(ind_geno[,c(2,1)], file = sprintf("%s%s_allele_file_as_txt.txt",my_wd, id_predicted_individual ), sep = "\t", quote = FALSE, row.names  = FALSE)
}


end_time <- Sys.time()
print(c("time (min): ", (end_time-start_time)/60))
## call breeding scheme script ##
# system("Rscript ")
## call sql script ##
# system("Rscript ")
#####################################################
# ## two examples on how the prediction can be used:
# #predict marker effects
# ans <- mixed.solve(y,Z=M)  #By default K = I
# accuracy <- cor(u,ans$u)
# 
# #predict breeding values
# ans <- mixed.solve(y,K=A.mat(M))
# accuracy <- cor(g,ans$u)