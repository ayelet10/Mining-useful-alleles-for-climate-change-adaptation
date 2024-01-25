### 1_GEA_breeding_scheme_genomicSimulation #### using lfmm2 ### (there is a separate R script for running with GEMMA) ##
start_time <- Sys.time()
print(start_time)
library(LEA)
library(vcfR)
library(data.table)
##############
## settings ##
##############
args <- commandArgs(trailingOnly = TRUE) # TRUE
print(args)
seed_num = args[1] # 1231227
# seed_num =  1231227
# seed_num =  1231118
tested_trait <-  "sal_opt" # "temp_opt" #
print(tested_trait)
# input.file_vcf <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_plusneut_MAF01.recode2.vcf.gz", seed_num) # non variable sites error ?
# input.file_vcf <- "/home/salmanay/projects/breeding_simulations/1231227/1231227_hundred_over_hundred_recode_vcf_out.vcf"    # non variable sites error ?
#input.file_vcf <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/vcf_files/%s_VCF_causal.vcf", seed_num)  # is read fine
   # 76 sites kept after vcftools -mac 1 for non variable sites removal
# input.file_vcf <- "/home/salmanay/projects/breeding_simulations/1231227/1231227_hundred_over_hundred_recode_vcf_out_mac.vcf.recode.vcf"

# my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/",seed_num)
my_input_dir <- sprintf("/home/salmanay/projects/breeding_simulations/%s/run3_all_causals/",seed_num)

# my_wd <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/"

## file without invariant sites
# input.file_vcf_mac <- sprintf("%s%s.1000_individuals_recode_vcf_out_1_Feb_mac1.vcf.recode.vcf", my_wd, seed_num)
# input.file_vcf_mac <- sprintf("%s%s_plusneut_MAF01.1000_random_samples_mac.recode.vcf", my_wd, seed_num)
# input.file_vcf_LD_pruned <- sprintf("%s%s.plusneut_MAF01.recode2_1000_random_samples_pruned_200_25_04.vcf", my_wd, seed_num)
# input.file_vcf_LD_pruned <- sprintf("%s%s.plusneut_MAF01_1000_random_samples_mac.recode.vcf", my_wd, seed_num) ## and invariate removed 
# input.file_vcf_LD_pruned <- sprintf("%s%s.plusneut_MAF01_1000_random_samples_maf_005.recode.vcf", my_wd, seed_num) ## and invariate removed 
input.file_vcf_LD_pruned <- sprintf("%s%s.plusneut_MAF01.recode2_1000_random_samples_pruned_200_25_04_vcf4_2.vcf", my_input_dir, seed_num) ## maf 0.05 filtered AND LD pruned

## vcf file with all site (including invariant)
# input.file_vcf <- sprintf("%s%s.1000_individuals_recode_vcf_out_1_Feb.vcf.recode.vcf", my_wd, seed_num)
# input.file_vcf <- sprintf("%s%s_plusneut_MAF01.recode2_1000_random_samples.recode.vcf", my_wd, seed_num)

## vcf file after replacing the ids in ID column (and removing and getting the header back and tab-delimiting the file)
#input.file_vcf <- sprintf("%ssorted_concat_full_vcf_with_causals.vcf", my_wd) #### sprintf("%sout1_added_header2_tabs.vcf", my_wd)
input.file_vcf <- sprintf("/home/salmanay/projects/breeding_simulations/%s/sorted_concat_full_vcf_with_causals.vcf", seed_num) ### sprintf("%sout1_added_header2_tabs.vcf", my_input_dir)


### added Nov 2023 ### vcf filtered for maf 0.05 to give the lfmm2.test function because of an invariate sites error
input.file_vcf_maf <- sprintf("/home/salmanay/projects/breeding_simulations/%s/sorted_concat_full_vcf_with_causals_maf_005.recode.vcf", seed_num)


# my_ind_data <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_ind.txt", seed_num)
# my_ind_data <- sprintf("%s%s_ind_subset_1000.txt", my_wd, seed_num) ## subset of 1000 randomly selected individuals
my_ind_data <-  sprintf("/home/salmanay/projects/breeding_simulations/%s/%s.ind_subset_1000_indID_sorted.txt", seed_num, seed_num)

# pruned_vcf_file <- "/home/salmanay/projects/breeding_simulations/1231227/1231227_plusneut_MAF01.recode2_1000_random_samples.recode.vcf_pruned.vcf"
# pruned_vcf_file <- "/home/salmanay/projects/breeding_simulations/1231227/1231227_plusneut_MAF01.recode2_1000_random_samples_pruned_100k_500_07.vcf"

# input.file_vcf <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/1231227_hundred_over_hundred_recode_vcf_out_mac.vcf.recode.vcf"    # 76 sites kept after vcftools -mac 1 for non variable sites removal
# input.file_vcf <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/1231227.1000_individuals_recode_vcf_out_25Jan_out_mac.vcf.recode.vcf"    # 76 sites kept after vcftools -mac 1 for non variable sites removal

# my_wd <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/"

# my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/gea_with_updated_traits_file/",seed_num)
my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/run3_all_causals/gea_with_updated_traits_file/",seed_num)

ped_file_name <- "sorted_concat_full_vcf_with_causals_recoded_vcftools_run3"
ped_file_path <- sprintf("/home/salmanay/projects/breeding_simulations/%s/%s", seed_num, ped_file_name)

dir.create(file.path(my_wd)) # if dir exists, it just gives a warning
setwd(file.path(my_wd))


## run local ##
# # my_wd <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/"
# # input.file_vcf <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/1231227.1000_individuals_recode_vcf_out_1_Feb.vcf.recode.vcf"
# seed_num =  1231227
# # input.file_vcf <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/1231227_hundred_over_hundred_recode_vcf_out_mac.vcf.recode.vcf"
# # # input.file_vcf <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/1231227.1_10_individuals_recode_vcf_out_16_Feb.vcf.recode.vcf"
# # my_ind_data <- sprintf("C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/%s_ind.txt", seed_num)
## input.file_vcf <-sprintf("C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/1231227.plusneut_MAF01.recode2_1000_random_samples_pruned_200_25_04.vcf")
## input.file_vcf <-sprintf("C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/1231227_VCF_causal.vcf")
###############
##  outputs  ## 
###############
best_markers_outfile <- "gea_best_markers_before_FDR_outfile.csv"
best_individuals_outfile     <- "gea_best_individuals_outfile.csv"
# mapped_best_markers_and_inds <- "gea_mapped_best_markers_and_inds_outfile.csv" ## "gea_mapped_best_markers_and_inds_outfile_July9.csv"
mapped_best_markers_and_inds <- "gea_mapped_best_markers_and_inds_outfile_2.csv" ## starting Dec 6 2023 - with a fix of unique individuals taken per marker
pv2_outfile <- sprintf("%s.pv2.txt",seed_num)
###############
## functions ## 
###############
# update marker names # removes the "chromosome" name and index name, keeps only the marker name # 1_1534 will be 1534
UpdateMarkerName <- function (data_with_markers_names){
  updated_marker_names <- sub("*._", "",data_with_markers_names )
  # updated_marker_names <- sub("_.*" , "" ,updated_marker_names )
  return(updated_marker_names)
}
## mapping of markers and individuals (who are the individuals with a specific marker etc.)
MapIndividualWithMarker <- function(my_genotypes_to_check, markers_lst, num_individuals_selected){
  ## for each marker, get the individuals that contain it# save to matrix #
  mapped_inds_to_markers <- c()
  for(m in 1:length(markers_lst)){
    randomly_sampled_curr_ind <- c()
    curr_mapped_marker <- c()
    marker_curr <- markers_lst[m]
    print(marker_curr)
    marker_curr_column_Y <- my_genotypes_to_check[which(rownames(my_genotypes_to_check) %in% marker_curr),]
    
    # marker_curr_column_Y <- my_genotypes_to_check[which(rownames(my_genotypes_to_check) == marker_curr),]
    print(head(marker_curr_column_Y))
    curr_selected_ind <- names(marker_curr_column_Y[marker_curr_column_Y %in% c("01", "10") ]) ## takes only individuals heterozygous for the marker

    # print(curr_selected_ind)
    if(length(curr_selected_ind) > num_individuals_selected){
      print("found more heterozygous than wanted, need to select among them")
      randomly_sampled_curr_ind <- sample(curr_selected_ind, num_individuals_selected) ## select random  individuals with the marker 
      curr_mapped_marker <- c(marker_curr,randomly_sampled_curr_ind)
      print(curr_mapped_marker)
    }else{
      print("found less heterozygous individuals than wanted, take all of them") ## although I guess chances are very low to not have 10 heterozygotes in pop of 1000
      ## add NA to curr_selected_ind to reach length of num_individuals_selected
      num_na_to_fill <- num_individuals_selected - length(curr_selected_ind)
      curr_selected_ind <- c(curr_selected_ind, rep(NA,num_na_to_fill))
      curr_mapped_marker <- c(marker_curr,curr_selected_ind)
      print(curr_mapped_marker)
    }
    mapped_inds_to_markers <- rbind(mapped_inds_to_markers, curr_mapped_marker)
  }
  colnames(mapped_inds_to_markers) <- c("marker",paste0(rep("individual",num_individuals_selected), 1:num_individuals_selected) )
  return(mapped_inds_to_markers)
}

UpdateMarkerNamesandDuplicants <- function(geno_mat){ 
  ## add a column to save the original names ##
  geno_mat$original_names <- geno_mat[,"name"]
  ## column with marker names to change should be "name" ##
  ## remove suffix with serial numbers ##
  geno_names <- geno_mat[,"name"]
  geno_mat[,"name"] <-  gsub("_[^_]+$", "", geno_names)
  ## remove prefix with mutation type ("1_")
  geno_mat[,"name"] <-  gsub(".*_", "",  geno_mat[,"name"])
  # geno_mat[,"name"] <- gsub("_", "-" , geno_mat[,"name"] )
  
  ## check for duplicated markers, by name ##
  ## for the second duplicated site, add _2 to the name ##
  duplicated_loci <- geno_mat[duplicated(geno_mat[,"name"]),] 
  unique_dupicated_sites <- unique(duplicated_loci[,"name"])

  for(u in unique_dupicated_sites){
    curr_dupl_rows <- duplicated_loci[duplicated_loci[,"name"] %in% u, ]
    i =1
    for(n in 1:nrow(curr_dupl_rows)){
      i = i+1
      # print(sprintf("%s_%s", curr_dupl_rows[n,"name"] , i))
      curr_dupl_rows[n,]$name <- sprintf("%s_%s", curr_dupl_rows[n,"name"] , i)
    }
    geno_mat[rownames(geno_mat) %in% rownames(curr_dupl_rows),] <- curr_dupl_rows 
  }
  # write.table(geno_mat[,c("name","original_names")], file = "key_marker_names_vcf_to_muts_all.txt", sep = "\t", quote = FALSE, row.names  = FALSE)
  return(geno_mat)
}

#######################
## read genomic data ##
#######################
vcf <- read.vcfR(input.file_vcf_LD_pruned, verbose = FALSE)
# vcf <- read.vcfR(input.file_vcf_mac, verbose = FALSE )
# vcf <- read.vcfR(input.file_vcf, verbose = FALSE )
my_genotypes <- extract.gt(vcf, IDtoRowNames =TRUE)  
# my_genotypes <- extract.gt(vcf)#, IDtoRowNames =FALSE)  ### IDtoRowNames set to false results with no markers names
# print(head(my_genotypes))
# my_genotypes <- gsub("|", "", my_genotypes, fixed = TRUE)  ## | in the original vcf (phased)
my_genotypes <- gsub("/", "", my_genotypes, fixed = TRUE)    ## / in the LD pruned vcf
# print(nrow(my_genotypes))
my_genotypes_no_na <- na.omit(my_genotypes) # exclude rows with NA
# my_genotypes_no_na <- gsub("|", "", my_genotypes_no_na, fixed = TRUE)
# my_genotypes_no_na_both_alleles <- gsub("|", "", my_genotypes_no_na, fixed = TRUE) ## to be used for using individuals files
## replace zeros and ones with numbers ## to avoid discarding of the leading zero in following steps #
my_genotypes_no_na[my_genotypes_no_na == "00"] <- "0"
my_genotypes_no_na[my_genotypes_no_na == "01"] <- "1"
my_genotypes_no_na[my_genotypes_no_na == "10"] <- "1"
my_genotypes_no_na[my_genotypes_no_na == "11"] <- "2"
my_genotypes_t <- t(my_genotypes_no_na)
num_sites <- nrow(my_genotypes_no_na)  #  14.12.22 >>>>  100 #
num_individuals_tested <- ncol(my_genotypes_no_na)
print(c("num_individuals_tested", num_individuals_tested))
num_markers_to_save <- 50
p_value_threshold <-  0.05  # 0.9# 
num_individuals_selected <- 10  ## number of individuals that have the marker, to be outputted to file
wanted_num_markers <- 50
threshold_percentage_pc_explains <- 0.01 ## number of pc's to take that explains at least 1% of the variance
#################
## run lfmm2  ### 
#################
# matrix of individuals and loci, in the lfmm format #
# Y <- my_hap_t[1:4,c(102,105,112,117, 121, 151, 156, 158,  175)]
# Y <- my_hap_t
# Y <- my_genotypes_t[1:2,c(105,112,117, 121, 151, 156, 158,  175)]
# Y <- my_genotypes_t[1:num_individuals_tested,c(105,112, 151, 156,  175)]
# Y <- my_genotypes_t[1:num_individuals_tested,c(1,3, 4, 8,  9)]
Y <- my_genotypes_t[1:num_individuals_tested,1:num_sites]
# Y <- my_genotypes_t[,1:num_sites] ## takes all 1000 randomly selected individuals data from subseted file 

print(Y[1:10,1:10])
class(Y) <- "numeric" # convert to numeric for the prcomp function ##  pay attention: in this step 00 will be converted to 0
print(Y[1:10,1:10])
## read phenotypic data from SliM simulation output #
ind_data <- read.table(my_ind_data, header = TRUE)
# head(ind_data)
print(c("number of unique values for tested trait", length(unique(ind_data[tested_trait,]))))
# unique(ind_data ["opt1"] )
# ind_data[c(1000:1004), c("opt0","opt1")]
# env_data <- ind_data[c(9990:9993), c("opt0","opt1")] ## subset for the four first individuals
# env_data <- ind_data[c(9990:9991), c("opt0","opt1")] ## subset for the two first individuals
# env_data <- ind_data[1:num_individuals_tested, c("opt0","opt1")]
X <- ind_data[1:num_individuals_tested, c(tested_trait)]  #
#################
## read vcf file with all sites
vcf_all <- read.vcfR(input.file_vcf, verbose = FALSE )
my_genotypes_all <- extract.gt(vcf_all)
# print(head(my_genotypes_all))
my_genotypes_all <- gsub("|", "", my_genotypes_all, fixed = TRUE)
# head(my_genotypes_all)
num_sites_all <- nrow(my_genotypes_all)
my_genotypes_all_no_na <- na.omit(my_genotypes_all)
# class(my_genotypes_all) <- "character"
my_genotypes_all_no_na[my_genotypes_all_no_na == "00"] <- "0"
my_genotypes_all_no_na[my_genotypes_all_no_na == "01"] <- "1"
my_genotypes_all_no_na[my_genotypes_all_no_na == "10"] <- "1"
my_genotypes_all_no_na[my_genotypes_all_no_na == "11"] <- "2"
my_genotypes_all_t <- t(my_genotypes_all_no_na) 

Y_all <- my_genotypes_all_t[1:num_individuals_tested,1:num_sites_all]
class(Y_all) <- "numeric"       ## to use in the test function
print(Y_all[1:10,1:10])
################# read vcf filtered for maf 0.05 , without LD pruning ##################
vcf_all_maf_005 <- read.vcfR(input.file_vcf_maf, verbose = FALSE )
my_genotypes_maf <- extract.gt(vcf_all_maf_005)
# print(head(my_genotypes_all))
my_genotypes_maf <- gsub("|", "", my_genotypes_maf, fixed = TRUE)
# head(my_genotypes_all)
num_sites_maf <- nrow(my_genotypes_maf)
my_genotypes_maf_no_na <- na.omit(my_genotypes_maf)
# class(my_genotypes_all) <- "character"
my_genotypes_maf_no_na[my_genotypes_maf_no_na == "00"] <- "0"
my_genotypes_maf_no_na[my_genotypes_maf_no_na == "01"] <- "1"
my_genotypes_maf_no_na[my_genotypes_maf_no_na == "10"] <- "1"
my_genotypes_maf_no_na[my_genotypes_maf_no_na == "11"] <- "2"
my_genotypes_maf_t <- t(my_genotypes_maf_no_na) 

Y_maf <- my_genotypes_maf_t[1:num_individuals_tested,1:num_sites_maf]
class(Y_maf) <- "numeric"       ## to use in the test function
print(Y_maf[1:10,1:10])


#################
###################
## select best K ###
####################

write.lfmm(Y_all,"genotypes.lfmm") # save as file to give to pca function, that its input is needed for T-W

# write.lfmm(Y,"genotypes.lfmm") # save as file to give to pca function, that its input is needed for T-W
pc <- pca("genotypes.lfmm")
tw <- tracy.widom(pc)
# Plot the percentage of variance explained by each component.
plot(tw$percentage)
# Display the p-values for the Tracy-Widom tests.
# print(tw$pvalues)
# K <- length(tw$pvalues[tw$pvalues < p_value_threshold ])  ###  >>> ???? not sure this is correct for selecting K? number of components that the pval < 0.05 ?
K <- length(tw$percentage[tw$percentage > threshold_percentage_pc_explains])


###################
# remove pca Project
remove.pcaProject("genotypes.pcaProject")
print(c("k from tracy-widom test: ", K))
if(K < 1){
  print("k is less than 1")
}
##################
# mod2 <- lfmm2(input = Y[1:10,1:10], env = X[1:10], K = K) # estimation function - recommended to use data after  LD pruning

mod2 <- lfmm2(input = Y, env = X, K = K) # estimation function - recommended to use data after  LD pruning
# mod2 <- lfmm2(input = Y_all, env = X, K = K) # estimation function - recommended to use data after  LD pruning

# Computing P-values and plotting their minus log10 values
## with LD pruned data ##
# pv2 <- lfmm2.test(object = mod2, input = Y, env = X, linear = TRUE, genomic.control = TRUE) # test function - recommended to use all genotypes (before LD pruning)
#### with all genotype data (no LD pruning) ##

## changed Nov 2023 - merged vcf (full with causal) filtered for maf 0.05 (without LD pruning)
pv2 <- lfmm2.test(object = mod2, input = Y_maf, env = X, linear = TRUE, genomic.control = TRUE) # test function - recommended to use all genotypes (before LD pruning)

# pv2 <- lfmm2.test(object = mod2, input = Y_all, env = X, linear = TRUE, genomic.control = TRUE) # test function - recommended to use all genotypes (before LD pruning)
# pv2 <- lfmm2.test(object = mod2, input = Y_all[1:10,1:30], env = X[1:10], linear = TRUE, genomic.control = TRUE) # test function - recommended to use all genotypes (before LD pruning)
write.table(pv2, file = pv2_outfile)  ## save the file as is 
## plots ##
pdf(file = "lfmm2_pv_plot_gif_plot.pdf")
plot(-log10(pv2$pvalues), col = "black", cex = .4, pch = 19)
############## GIF ####
pv2$gif  ## genomic inflation factor. make sure it is around 1 (more than one is more liberal)
#
# # p-values histogram # http://varianceexplained.org/statistics/interpreting-pvalue-histogram/
# #Compute the GIF
# lambda = median(pv2$zscores^2)/0.456 # the gif calculation
# ## check histogram of p values to make sure it is flat with a peak near zero
# hist(pv2$pvalues, col = "red") # to check that the genomic control correction (in lfmm2.test) was fine

# lambda = 1 ##  to decide if needed to be change and to what value - to check histogram by eye
# adj.p.values = pchisq(pv2$zscores^2/lambda, df = 1, lower = FALSE)
# # #histogram of p-values
# hist(adj.p.values, col = "red")
# dev.off()
# compute adjusted p-values from the combined z-scores
# adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)
##########################
## get the best markers ##
##########################
######### FDR ##########
## update threshold until we get 50 markers or reached 0.5 threshold ##
indices_good_adj_p_values <- c()
fdr_pvalue_jump <- 0 ## for the first iteration to start with p_value_threshold of 0.05

## making a one sided test (so that the negative effect sizes will get p-value > 0.5 and will not be significant)
## (the default output is a result of a two-sided t-test)
two_sided_p = pv2$pvalues
z_scores = pv2$zscores
positive_effects = z_scores>0
one_sided_p = two_sided_p/2
one_sided_p[!positive_effects] = 1-one_sided_p[!positive_effects]
summary(one_sided_p)
adj_p_values = p.adjust(one_sided_p, method='BH')
summary(adj_p_values)
## correct p-values. Can use different methods "fdr", "BH" (they are very similar)
# adj_p_values <- p.adjust(pv2$pvalues, method = "BH", n = length(pv2$pvalues)) # changed to BH

## file with adj_p_values and marker names to give to --clump ## 
pvalues_pv2 <- pv2$pvalues
pv2_sub_names <- sub("Response ", "", rownames(pvalues_pv2))
names_and_adj_p_values <- cbind(pv2_sub_names, adj_p_values)
colnames(names_and_adj_p_values) <- c("SNP", "P" )
write.table(names_and_adj_p_values, file = "adj_p_values.csv", quote = FALSE, row.names = FALSE) ## needs to be a tab-delimited with columns: SNP  P  ### --clump-p2 0.01 (the default)
head(names_and_adj_p_values)
#####  take pvalues until reached 50 markers or threshold of 0.5 #####
while(length(indices_good_adj_p_values) < wanted_num_markers){
  print("number of markers in current iteration:")
  print(length(indices_good_adj_p_values))
  ## correct p-values. Can use different methods "fdr", "BH" (they are very similar)
  # adj_p_values <- p.adjust(pv2$pvalues, method = "BH", n = length(pv2$pvalues)) # changed to BH
  indices_good_adj_p_values <- which(adj_p_values < p_value_threshold)
  if(p_value_threshold < 0.49){
    fdr_pvalue_jump <- 0.05 #
    p_value_threshold <- p_value_threshold + fdr_pvalue_jump
    print(p_value_threshold)
  }else{
      break
  }
}

system(sprintf("plink --file %s --clump adj_p_values.csv --clump-p1 %s --clump-p2 %s --clump-kb 15" , ped_file_path, p_value_threshold, p_value_threshold*10) )## calls plink 1.07 # define the pvalue threshold
clumping_res <- read.table(sprintf("%splink.clumped", my_wd), skip = 1) ## the file result from the clumping
head(clumping_res)


print(c("p_value_threshold: ", p_value_threshold))
print(c("number of markers below threshold: ",  length(indices_good_adj_p_values)))
while((length(indices_good_adj_p_values) > 49) & (nrow(clumping_res) < 49) & (p_value_threshold < 0.49)){
  print("re-evaluating pvalue threshold")
  ## in case there were 50 markers but some were removed after the clumping, elevate the fdr pvalue to get more markers and do clumping again 
  ##  allow higher p-value to get more markers (because some of the markers were found as clustered by the clumping)
  p_value_threshold <- p_value_threshold + fdr_pvalue_jump
  indices_good_adj_p_values <- which(adj_p_values < p_value_threshold)
  ## re-run the clumping ##
#  system(sprintf("plink --file %s --clump adj_p_values.csv --clump-p1 %s --clump-p2 %s --clump-kb 15" ,ped_file_name, p_value_threshold, p_value_threshold*10) )## calls plink 1.07 # define the pvalue threshold
  system(sprintf("plink --file %s --clump adj_p_values.csv --clump-p1 %s --clump-p2 %s --clump-kb 15" ,ped_file_path, p_value_threshold, p_value_threshold*10) )## calls plink 1.07 # define the pvalue threshold
  
  clumping_res <- read.table(sprintf("%splink.clumped", my_wd), skip = 1) ## the file result from the clumping
}

plot(-log10(adj_p_values))
abline(h=-log10(p_value_threshold), col="blue")
text(1,0, sprintf("abline( h = %s)", round(-log10(p_value_threshold), digits = 4)), col = "gray60", adj = c(0, -.1))

# ## check histogram of p values to make sure it is flat with a peak near zero
hist(pv2$pvalues, col = "red") # to check that the genomic control correction (in lfmm2.test) was fine
hist(adj_p_values, col = "red") # to check that the genomic control correction (in lfmm2.test) was fine
dev.off()
print(c("p_value_threshold: ", p_value_threshold))
print(c("number of markers below threshold: ",  length(indices_good_adj_p_values)))
snps_independent <- clumping_res[order(clumping_res$V5),]$V3
print(snps_independent)

pv2_below_threshold_after_fdr <- pv2$pvalues[indices_good_adj_p_values, ]

best_markers_names <- names(sort(pv2_below_threshold_after_fdr))
best_markers_names <- sub("Response ", "", best_markers_names)
print(best_markers_names[1:10])
print(sprintf("best markers, pvalue < %s", p_value_threshold))


# write.csv(best_markers_names_before_fdr, file = paste0(my_wd, best_markers_outfile))
write.csv(best_markers_names, file = paste0(my_wd, "best_GEA_markers_after_fdr.csv"))
if(length(unique(best_markers_names))==1){
  if(is.na(unique(best_markers_names))){
    stop("no markers found.Stopping", call. = FALSE)
  }
}

sub_names <- sub("Response ", "",names(pv2_below_threshold_after_fdr))
names(pv2_below_threshold_after_fdr) <- sub_names
write.table(pv2_below_threshold_after_fdr, file = "pv2_below_threshold_after_fdr.csv", quote = FALSE) ## needs to be a tab delimited with columns: SNP  P
###  clumping (by pvalue threshold) ### >>> should use pv2 or adjusted pvalues ??

# rownames(pvalues_pv2) <- pv2_sub_names
# pvalues_pv2$SNPS <- c()
# pvalues_pv2$SNPS <- pv2_sub_names
# to_save_pvalues_pv2 <- cbind(pv2_sub_names,pvalues_pv2 )
# length(pv2_sub_names)
# nrow(to_save_pvalues_pv2)
# colnames(to_save_pvalues_pv2) <- c("SNP", "P" )
# write.table(to_save_pvalues_pv2, file = "pv2.csv", quote = FALSE, row.names = FALSE) ## needs to be a tab-delimited with columns: SNP  P  ### --clump-p2 0.01 (the default)
# system(sprintf("plink --file out1_added_header2_tabs_recoded_vcftools --clump pv2.csv --clump-kb 15" ) )## calls plink 1.07 # define the pvalue threshold
# system(sprintf("plink --file out1_added_header2_tabs_recoded_vcftools --clump pv2.csv --clump-p1 %s --clump-p2 %s --clump-kb 15" , p_value_threshold, p_value_threshold*10) )## calls plink 1.07 # define the pvalue threshold




# plink --file /home/salmanay/projects/breeding_simulations/1231213/out1_added_header2_tabs_recoded_vcftools --clump /home/salmanay/projects/breeding_simulations/1231213/pv2_below_threshold_after_fdr.csv

### # # # #
# print all sites p-values to file                >> only for testing. remove later
# write.csv(pv2, file = paste0(my_wd,"lfmm2_test_results.csv"))

## check to make sure that there are markers selected #
if(length(unique(snps_independent))==0){
  stop("no independent snps found. Stopping", call. = FALSE)
}
# print("markers found. Continue to breeding scheme")

independent_best_markers_names <- best_markers_names[best_markers_names %in% snps_independent]
print("best markers after clumping test:")
print(independent_best_markers_names)
### take best 50 markers ###
fifty_best_markers <- independent_best_markers_names[1:num_markers_to_save] ##
fifty_best_markers <- na.omit(fifty_best_markers)
print(head(fifty_best_markers))
##  only HETEROZYGOUS and not homozygous individuals are taken (for the breeding scheme to get segregation of the marker in the F1 pop)
## it is defined in the function MapIndividualWithMarker ##
# mapped_markers_and_inds <- MapIndividualWithMarker(my_genotypes, fifty_best_markers, num_individuals_selected) # best_markers_names
mapped_markers_and_inds <- MapIndividualWithMarker(my_genotypes_all, fifty_best_markers, num_individuals_selected) # best_markers_names

print(head(mapped_markers_and_inds))
write.csv(mapped_markers_and_inds, file = paste0(my_wd, mapped_best_markers_and_inds))
saved_individuals_lst <- c()
## save first individual per marker to file #
for(id_p in 1:nrow(mapped_markers_and_inds)){
  print(id_p)
  # print(mapped_markers_and_inds[id_p])
  predicted_marker <- mapped_markers_and_inds[id_p,"marker"] #
  ind_to_save <- mapped_markers_and_inds[id_p,"individual1"] #mapped_markers_and_inds[id_p,2] second column is the first individual with the marker
  if(ind_to_save %in% saved_individuals_lst){
    ## if individual was selected to represent previous marker, take the next individual
    ind_to_save <- mapped_markers_and_inds[id_p,"individual2"] #mapped_markers_and_inds[id_p,2] second column is the first individual with the marker
    mapped_markers_and_inds[id_p,"individual1"] <- ind_to_save ## replace the first individual with the second that was not already picked by other markers
  }
  if(ind_to_save %in% saved_individuals_lst){
    ## if individual was selected to represent previous marker, take the next individual
    ind_to_save <- mapped_markers_and_inds[id_p,"individual3"] #mapped_markers_and_inds[id_p,2] second column is the first individual with the marker
    mapped_markers_and_inds[id_p,"individual1"] <- ind_to_save ## replace the first individual with the third that was not already picked by other markers
  }

  saved_individuals_lst <- c(saved_individuals_lst,ind_to_save)
  print(c("predicted_marker:", predicted_marker))
  print(c("individual_with_the_marker_to_save:", ind_to_save))
#  ind_geno <- my_genotypes_all[,ind_to_save]  ## if not working - check to make sure the relevant individuals are in
  # ind_geno <- my_genotypes[,ind_to_save] ## works with this file

  ### change Nov 2023 to keep the "_c" in causal mutations names ###
  ind_geno <- as.matrix(my_genotypes_all[,ind_to_save]) ## the vcf without maf filter
  ind_geno <- cbind(ind_geno, rownames(ind_geno))
  colnames(ind_geno) <- c(ind_to_save,"name")
  ### added to treat the causal mutation names that have "_c" in addition to the "_2"##
  ## first, update the names that are not the causal (without "_c" in the name)
  df_ind_geno <- as.data.frame(ind_geno)
  df_ind_geno$original_names <- df_ind_geno$name ## df_ind_geno$ind_geno #
  df_ind_geno$index <- 1:nrow(df_ind_geno)
  test_ind_no_c <- UpdateMarkerNamesandDuplicants(df_ind_geno[!grepl("_c", df_ind_geno$name),])

  ## keep the causals, to be merged after fixing the names of non-causals ##
  with_c <-  df_ind_geno[grepl("_c", df_ind_geno$name),]

  ## then, merge the fixed non-causal with the causals and sort by the index column that was added
  after_fix <- rbind(test_ind_no_c, with_c[grepl("_c", with_c$name),])
  sorted_after_fix <- after_fix[ order(after_fix$index),]
  print(head(sorted_after_fix))

  ## create a genotype file in the format to be read by genomicSimulation ##
#  ind_geno <- cbind(ind_geno, names(ind_geno))
#  colnames(ind_geno) <- c(ind_to_save,"name")
  
  
  ## update markers names to be as in muts_full format (remove serial number, mut type number and add number to markers in duplicated positions)
#  ind_geno <- UpdateMarkerNamesandDuplicants(as.data.frame(ind_geno))
  ## save individual data to file ##
  # write.table(ind_geno[,c(2,1)], file = sprintf("%s%s_updated_markers_names_allele_file_as_txt2.txt",my_wd, ind_to_save), sep = "\t", quote = FALSE, row.names  = FALSE)
#  write.table(ind_geno[,c("name",ind_to_save)], file = sprintf("%s%s_updated_markers_names_allele_file_as_txt_GEA.txt", my_wd, ind_to_save), sep = "\t", quote = FALSE, row.names  = FALSE)
  print("before printing to file")
  write.table(sorted_after_fix[,c(2,1)], file = sprintf("%s%s_updated_markers_names_allele_file_as_txt_GEA.txt",my_wd, ind_to_save ), sep = "\t", quote = FALSE, row.names  = FALSE)
  ind_to_save <- c()
}

end_time <- Sys.time()
print(end_time)
print(c("time (min): ", (end_time-start_time)/60))

## relevant tutorial:
## https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html
## clumping with plink:
## https://zzz.bwh.harvard.edu/plink/clump.shtml
##############################################################