## Prediction with rrBLUP mixed.solve ###
## output to be used for the following 2 approaches :
##  - breeding scheme of prediction-only approach
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
tested_trait <- "sal_opt"  #"temp_opt" #   
input.file_vcf <- "/home/salmanay/projects/breeding_simulations/1231227/1231227.1000_individuals_recode_vcf_out_25Jan.vcf"
my_ind_data <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_ind.txt", seed_num)
my_pop_data <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_popInfo.txt", seed_num) 
my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/",seed_num )

dir.create(file.path(my_wd)) # if dir exists, it just gives a warning
setwd(file.path(my_wd))
###############
## read data ##
###############
vcf <- read.vcfR(input.file_vcf, verbose = FALSE )
my_geno <- extract.gt(vcf) 
my_geno_no_na <- na.omit (my_geno) 
my_geno_no_na <- my_geno_no_na     
num_sites              <- nrow(my_geno_no_na)
num_individuals_tested <- ncol(my_geno_no_na) 
print(c("number of individuals:" ,num_individuals_tested ))
print(c("number of sites:" , num_sites))
my_geno_no_na[1:num_sites,1:num_individuals_tested] <- gsub("|", "", my_geno_no_na[1:num_sites,1:num_individuals_tested], fixed = TRUE)

## read phenotype data ##
ind_file <- read.table(file = my_ind_data, header = TRUE, stringsAsFactors = FALSE)
my_y <- ind_file[1:num_individuals_tested,c(tested_trait)]
#########################
## genomic prediction  ## 
#########################
## rrBLUP
my_W2 <- matrix (as.numeric(t(my_geno_no_na)), ncol = num_sites, nrow =num_individuals_tested) # genomic data
row.names(my_W2) <- colnames(my_geno_no_na)
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
indices_ordered <- order(RRBLUP2$u,decreasing=T)  #[1:50]
my_W2_ordered <- my_W2[indices_ordered,]  
ordered_u <- RRBLUP2$u[indices_ordered]
## get the id of the predicted
predicted_lr <- rownames(my_W2_ordered) 
predicted_lr <- na.omit(predicted_lr[1:50] ) 
####################################
# save selected landraces to file ##
####################################
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
