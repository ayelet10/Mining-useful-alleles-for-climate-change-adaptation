## Prediction with rrBLUP mixed.solve ###
## output to be used for the following 2 approaches :
##  - breeding scheme of prediction only approach
##  - qtl mapping approach 
start_time <- Sys.time()
print(start_time)
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
# input_data_dir <- "/home/salmanay/projects/Lotterhos_preprint_data/sim_output_20220428/"
# input.file_vcf <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/Case1.10krows.twoindivs.vcf"
# input.file_vcf <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/Case1.10krows.tenindivs.vcf" # 10 inds
# input.file_vcf <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/vcf_files/%s_VCF_causal.vcf", seed_num) 
# input.file_vcf <- sprintf("C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/%s_VCF_causal.vcf", seed_num)
# input.file_vcf <- "/home/salmanay/projects/simulated_data_laurson2021/tests_on_subseted_vcf_file/Case1.10krows.twoindivs.vcf"

# input.file_vcf <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_plusneut_MAF01.recode2.vcf.gz", seed_num) 
# input.file_vcf <- "/home/salmanay/projects/breeding_simulations/1231227/1231227.1000_individuals_recode_vcf_out_25Jan.vcf"
# input.file_vcf <- sprintf("/home/salmanay/projects/breeding_simulations/%s/%s_plusneut_MAF01.recode2_1000_random_samples.recode.vcf" ,seed_num, seed_num)
# input.file_vcf <- sprintf("/home/salmanay/projects/breeding_simulations/%s/run3_all_causals/%s_plusneut_MAF01.recode2_1000_random_samples.recode.vcf" ,seed_num, seed_num)
input.file_vcf <- sprintf("/home/salmanay/projects/breeding_simulations/%s/sorted_concat_full_vcf_with_causals.vcf" ,seed_num) ## run3
input.file_vcf_maf <- sprintf("/home/salmanay/projects/breeding_simulations/%s/sorted_concat_full_vcf_with_causals_maf_005.recode.vcf" ,seed_num) ## run3

# input.file_vcf <- sprintf("%s_plusneut_MAF01.recode2_1000_random_samples.recode.vcf" ,seed_num)

# my_ind_data <-"C:/Users/ayelet/Documents/UCDavis/scripts/simulations/Case1_1442299973452_ind.txt"
# my_ind_data <- "/home/salmanay/projects/simulated_data_laurson2021/Case1_1442299973452_ind.txt"
# my_ind_data <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_ind.txt", seed_num)
# my_ind_data <-"C:/Users/ayelet/Documents/UCDavis/scripts/simulations/data_from_Lotterhos_preprint_2022/1231227_ind.txt"
# my_ind_data <- sprintf("/home/salmanay/projects/breeding_simulations/%s/%s_ind_subset_1000.txt", seed_num, seed_num) ## subset of 1000 randomly selected individuals
my_ind_data <- sprintf("/home/salmanay/projects/breeding_simulations/%s/%s.ind_subset_1000_indID_sorted.txt", seed_num, seed_num) ## updated August 2023 sorted by individuals as in the vcf file ## subset of 1000 randomly selected individuals

# my_pop_data <- sprintf("%s/%s_popInfo.txt", input_data_dir, seed_num) # opt0 - salinity; opt1 -tempinput_data_dir
# my_pop_data <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_popInfo.txt", seed_num) # opt0 - salinity; opt1 -temp
# my_pop_data  <-  "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/data_from_Lotterhos_preprint_2022/1231227_popInfo.txt"

# my_wd <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/"
# my_wd <- "/home/salmanay/projects/breeding_simulations/genomicSimulation/"
# my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/",seed_num )
# my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/results_pred_qtl_updated_map_and_ind_files/",seed_num )
my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/run3_all_causals/",seed_num )


dir.create(file.path(my_wd)) # if dir exists, it just gives a warning
setwd(file.path(my_wd))
#################################################
## function from R script 1_GEA_genomicSimulation.R  ##
UpdateMarkerNamesandDuplicants <- function(geno_mat){ 
  ## add a column to save the original names ##
  geno_mat$original_names <- geno_mat[,"name"]
  ## column with marker names to change should be "name" ##
  ## remove suffix with serial numbers ##
  geno_names <- geno_mat[,"name"]
  geno_mat[,"name"] <-  gsub("_[^_]+$", "", geno_names)
  ## remove prefix with mutation type ("1_")
  geno_mat[,"name"] <-  gsub(".*_", "",  geno_mat[,"name"])
  ## check for duplicated markers, by name ##
  ## for the second duplicated site, add _2 to the name ##
  duplicated_loci <- geno_mat[duplicated(geno_mat[,"name"]),] 
  unique_dupicated_sites <- unique(duplicated_loci[,"name"])
  for(u in unique_dupicated_sites){
    curr_dupl_rows <- duplicated_loci[duplicated_loci[,"name"] %in% u, ]
    i =1
    for(n in 1:nrow(curr_dupl_rows)){
      i = i+1
      curr_dupl_rows[n,]$name <- sprintf("%s_%s", curr_dupl_rows[n,"name"] , i)
    }
    geno_mat[rownames(geno_mat) %in% rownames(curr_dupl_rows),] <- curr_dupl_rows 
  }
  write.table(geno_mat[,c("name","original_names")], file = "key_marker_names_vcf_to_muts_all.txt", sep = "\t", quote = FALSE, row.names  = FALSE)
  return(geno_mat)
}

###############
## read data ##
###############
vcf_maf <- read.vcfR(input.file_vcf_maf, verbose = FALSE )
my_geno <- extract.gt(vcf_maf) #, as.numeric = TRUE)
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
my_y <- ind_file[1:num_individuals_tested_full, c(tested_trait)]


### read vcf before maf 0.05 filter - to be used in writing the genotypes to file ####
vcf <- read.vcfR(input.file_vcf, verbose = FALSE )
my_geno_full <- extract.gt(vcf) #, as.numeric = TRUE)
# my_geno[1:10,1:10]
my_geno_full_no_na <- na.omit (my_geno_full) # exclude rows with NA       >>> to understand: why there are NAs in the simulated data??
my_geno_full_no_na <- my_geno_full_no_na     #[1:100,95:105]
print(my_geno_full_no_na[1:4, 1:10]) # to see that the data was read fine
num_sites_full              <- nrow(my_geno_full_no_na)
num_individuals_tested_full <- ncol(my_geno_full_no_na)# 1000 #1000 # 
print(c("number of individuals full:" ,num_individuals_tested_full ))
print(c("number of sites full(no maf filter):" , num_sites_full))
my_geno_full_no_na[1:num_sites_full,1:num_individuals_tested_full] <- gsub("|", "", my_geno_full_no_na[1:num_sites_full,1:num_individuals_tested_full], fixed = TRUE)


#########################
## genomic prediction  ## 
#########################
## rrBLUP
# y = # observations

my_W2 <- matrix(as.numeric(t(my_geno_no_na)), ncol = num_sites, nrow =num_individuals_tested) # genomic data
# indices_ordered <- order(1:num_individuals_tested, decreasing=T)
# print(indices_ordered)
# my_W2_ordered <- my_W2[c(indices_ordered),]
# print()
# # unique(my_W2)
row.names(my_W2) <- colnames(my_geno_no_na)
# RRBLUP2 <- mixed.solve(y =my_y, Z = my_W2) # y - environment; Z - genomic data
RRBLUP2 <- mixed.solve(y =my_y, K = A.mat(my_W2), method="REML") # y - environment; K - genomic data; default method is REML (restricted), ML (full) is optional.

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

# ordered_u <- RRBLUP2$u[indices_ordered]
## get the id of the predicted
predicted_lr <- rownames(my_W2_ordered) ## get the ids for best genotype
# predicted_lr <- cbind(predicted_lr, ordered_u)
predicted_lr <- na.omit(predicted_lr[1:50] ) ## get the ids for best genotype
####################################
# save selected landraces to file ##
####################################
# write.csv(ordered_u, file = paste0(my_wd,"u_mixed_solve.csv") ) # one columns : ordered_u
write.csv(predicted_lr, file = paste0(my_wd,"ids_predicted_landraces_mixed_solve.csv") , row.names = FALSE) # one columns : predicted_lr
# # write.csv(my_W2_ordered, file = paste0(my_wd,"genotypes_predicted_landraces_mixed_solve.csv"))
# predicted_lr<- read.csv(paste0(my_wd,"ids_predicted_landraces_mixed_solve.csv"))
print(Sys.time())
print("before writing genotypes to txt files")

## write genotype of each predicted individual to single txt file ## in addition, save all 50 into a one file ##
all_predicted_indiv <- c()
for(id_p in 1:length(predicted_lr)){
  print(predicted_lr[id_p])
  id_predicted_individual <- predicted_lr[id_p]
  ind_geno <- as.matrix(my_geno_full_no_na[, predicted_lr[id_p]]) ## the vcf without maf filter
  ind_geno <- cbind(ind_geno, rownames(ind_geno))
  colnames(ind_geno) <- c(id_predicted_individual,"name")
  
  ### added to treat the causal mutation names that have "_c" in addition to the "_2"## 
  ## first, update the names that are not the causal (without "_c" in the name)
  df_ind_geno <- as.data.frame(ind_geno) 
  df_ind_geno$original_names <- df_ind_geno$name
  df_ind_geno$index <- 1:nrow(df_ind_geno)
  test_ind_no_c <- UpdateMarkerNamesandDuplicants(df_ind_geno[!grepl("_c", df_ind_geno$name),])
  
  ## keep the causals, to be merged after fixing the names of non-causals ##
  with_c <-  df_ind_geno[grepl("_c", df_ind_geno$name),]
  
  ## then, merge the fixed non-causal with the causals, remove the "_c" from the names and sort by the index column that was added 
  after_fix <- rbind(test_ind_no_c, with_c[grepl("_c", with_c$name),]) 
  ## try not to remove the "_c"  >>> added it to the tester, map and effect files
  # after_fix$name <- gsub("_c", "", after_fix$name)
  sorted_after_fix <- after_fix[ order(after_fix$index),]
  ###
  
  # source("/home/salmanay/projects/breeding_simulations/genomicSimulation/scripts/1_GEA_genomicSimulation.R") # to use the function UpdateMarkerNamesandDuplicants
  ## update markers names to be as in muts_full format (remove serial number and add number to markers in duplicated positions)
#  ind_geno <- UpdateMarkerNamesandDuplicants(as.data.frame(ind_geno))  ## needs to be tested in this script 
  write.table(sorted_after_fix[,c(2,1)], file = sprintf("%s%s_allele_file_as_txt.txt",my_wd, id_predicted_individual ), sep = "\t", quote = FALSE, row.names  = FALSE)
  # write.table(ind_geno[,c(2,1)], file = sprintf("%s%s_allele_file_as_txt_fixed.txt",my_wd, id_predicted_individual ), sep = "\t", quote = FALSE, row.names  = FALSE)
  if(id_p >1 ){
      all_predicted_indiv <- cbind(all_predicted_indiv, sorted_after_fix[,c(1)])
  }else{ # first occurrence, take marker names column 
    all_predicted_indiv <- sorted_after_fix[,c(2,1)]
  }

}
colnames(all_predicted_indiv) <- c("name" , predicted_lr)
write.table(all_predicted_indiv, file = sprintf("%s%s_all_predicted_indiv_allele_file_as_txt.txt",my_wd, seed_num ), sep = "\t", quote = FALSE, row.names  = FALSE)


end_time <- Sys.time()
print(end_time)
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