### 1_GEA_breeding_scheme_genomicSimulation #### 
start_time <- Sys.time()
print(start_time)
# library(lfmm)
library(LEA)
library(vcfR)
library(genomicSimulation)
library(data.table)
##############
## settings ##
##############
args <- commandArgs(trailingOnly = TRUE) # TRUE
print(args)
seed_num = args[1] # 1231227 #
tested_trait <-  "sal_opt" # "temp_opt" #
print(tested_trait)

input.file_vcf <- "/home/salmanay/projects/breeding_simulations/1231227/1231227.1000_individuals_recode_vcf_out_1_Feb.vcf.recode.vcf"

my_ind_data <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_ind.txt", seed_num)
my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/",seed_num )

dir.create(file.path(my_wd)) # if dir exists, it just gives a warning
setwd(file.path(my_wd))

###############
##  outputs  ## 
###############
best_markers_outfile <- "gea_best_markers_before_FDR_outfile.csv"
best_individuals_outfile     <- "gea_best_individuals_outfile.csv"
mapped_best_markers_and_inds <- "gea_mapped_best_markers_and_inds_outfile.csv"
p_value_threshold <-  0.05  
num_individuals_selected <- 10
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
MapIndividualWithMarker <- function (my_genotypes_to_check, markers_lst){
  ## for each marker, get the individuals that contain it# save to matrix #
  mapped_inds_to_markers <- c()
  for(m in 1:length(markers_lst)){
    randomly_sampled_curr_ind <- c()
    curr_mapped_marker <- c()
    marker_curr <- markers_lst[m]
    marker_curr_column_Y <- my_genotypes_to_check[which(rownames(my_genotypes_to_check) == marker_curr),]
    curr_selected_ind <- names(marker_curr_column_Y [marker_curr_column_Y %in% c("01", "10") ]) ## takes only heterozygous individuals for the marker
    if(length(curr_selected_ind) > num_individuals_selected){
      randomly_sampled_curr_ind <- sample(curr_selected_ind, num_individuals_selected) ## select random  individuals with the marker 
      curr_mapped_marker <- c(marker_curr,randomly_sampled_curr_ind)
    }else{
      ## add NA to curr_selected_ind to reach length of num_individuals_selected
       num_na_to_fill <- num_individuals_selected - length(curr_selected_ind)
       curr_selected_ind <- c(curr_selected_ind, rep(NA,num_na_to_fill))
    }
    curr_mapped_marker <- c(marker_curr,randomly_sampled_curr_ind)
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
      print(sprintf("%s_%s", curr_dupl_rows[n,"name"] , i))
      curr_dupl_rows[n,]$name <- sprintf("%s_%s", curr_dupl_rows[n,"name"] , i)
    }
    geno_mat[rownames(geno_mat) %in% rownames(curr_dupl_rows),] <- curr_dupl_rows 
  }
  write.table(geno_mat[,c("name","original_names")], file = "key_marker_names_vcf_to_muts_all.txt", sep = "\t", quote = FALSE, row.names  = FALSE)
  return(geno_mat)
}
#######################
## read genomic data ##
#######################
vcf <- read.vcfR(input.file_vcf_mac, verbose = FALSE )
my_genotypes <- extract.gt(vcf)
my_genotypes <- gsub("|", "", my_genotypes, fixed = TRUE)
print(nrow(my_genotypes))
my_genotypes_no_na <- na.omit (my_genotypes) # exclude rows with NA
my_genotypes_no_na <- gsub("|", "", my_genotypes_no_na, fixed = TRUE)
## replace zeros and ones with numbers ## to avoid discarding of the leading zero in following steps #
my_genotypes_no_na[my_genotypes_no_na == "00"] <- "0"
my_genotypes_no_na[my_genotypes_no_na == "01"] <- "1"
my_genotypes_no_na[my_genotypes_no_na == "10"] <- "1"
my_genotypes_no_na[my_genotypes_no_na == "11"] <- "2"
my_genotypes_t <- t(my_genotypes_no_na)
# impute()
num_sites <- nrow(my_genotypes_no_na)  
num_individuals_tested <- ncol(my_genotypes_no_na)
num_markers_to_save <- 50
#################
## run lfmm2  ###  
#################
# matrix of individuals and loci, in the lfmm format #
Y <- my_genotypes_t[1:num_individuals_tested,1:num_sites]

class(Y) <- "numeric" # convert to numeric for the prcomp function ##  pay attention: in this step 00 will be converted to 0

## read phenotypic data from SliM simulation output #
ind_data <- read.table(my_ind_data, header = TRUE)
print(c("number of unique values for tested trait", length(unique(ind_data [tested_trait] ))))
X <- ind_data[1:num_individuals_tested, c(tested_trait)]  #

###################
## select best K ###
####################
write.lfmm(Y,"genotypes.lfmm") # save as file to give to pca function, that its input is needed for T-W
pc <- pca("genotypes.lfmm")
tw <- tracy.widom(pc)
# Plot the percentage of variance explained by each component.
# plot(tw$percentage)
# Display the p-values for the Tracy-Widom tests.
print(tw$pvalues)
K <- length(tw$pvalues[tw$pvalues < p_value_threshold ])  
# remove pca Project
remove.pcaProject("genotypes.pcaProject")
print(c("k from tracy-widom test: ", K))
if(K < 1){
  print("k less than 1")
}
##################
mod2 <- lfmm2(input = Y, env = X, K = K)
# Computing P-values and plotting their minus log10 values
pv2 <- lfmm2.test(object = mod2, input = Y, env = X, linear = TRUE, genomic.control = TRUE)
# print(pv2)
pdf(file = "lfmm2_pv_plot_gif_plot.pdf")
plot(-log10(pv2$pvalues), col = "black", cex = .4, pch = 19)

############## GIF ####
# pv2$gif  ## genomic inflation factor. make sure it is around 1 (more than one is more liberal)
# ## check histogram of p values to make sure it is flat with a peak near zero
hist(pv2$pvalues, col = "red") # to check that the genomic control correction (in lfmm2.test) was fine
dev.off()
# compute adjusted p-values from the combined z-scores
# adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)

######### FDR ##########
## correct p-values. Can use different methods fdr, "BH" (they are very similar)
adj_p_values <- p.adjust(pv2$pvalues, method = "BH", n = length(pv2$pvalues)) # changed to BH
names_pv2 <- rownames(unlist(pv2$pvalues))
##########################
## get the best markers ##
##########################
indices_good_adj_p_values   <- which(adj_p_values < p_value_threshold)
pv2_below_threshold_after_fdr <- pv2$pvalues[indices_good_adj_p_values, ]
best_markers_names <- names(sort(pv2_below_threshold_after_fdr))

# # ## without the fdr correction:
pv2_below_threshold <- pv2$pvalues[which(unlist(pv2$pvalues) < p_value_threshold),]
best_markers_names_before_fdr  <- names(sort(pv2_below_threshold))
print(best_markers_names_before_fdr[1:10])

print("after FDR:")
print(best_markers_names[1:10])
best_markers_names <- sub("Response ", "", best_markers_names)
print(sprintf("best markers, pvalues < %s", p_value_threshold))
write.csv(best_markers_names_before_fdr, file = paste0(my_wd, best_markers_outfile))
write.csv(best_markers_names, file = paste0(my_wd, "best_GEA_markers_after_fdr.csv"))

## check to make sure that there are markers selected #
if(length (unique(best_markers_names))==1){
  if(is.na(unique(best_markers_names))){
   stop("no markers found.Stopping", call. = FALSE)
  }
}
print("markers found. Continue to breeding scheme")

### take best 50 markers ###
fifty_best_markers <- best_markers_names[1:num_markers_to_save] ##
fifty_best_markers <- na.omit(fifty_best_markers)

##  only HETEROZYGOUS and not homozygous individuals are taken (for the breeding scheme to get segregation of the marker in the F1 pop)
## it is defined in the function MapIndividualWithMarker ##
mapped_markers_and_inds <- MapIndividualWithMarker(my_genotypes, fifty_best_markers) # best_markers_names
write.csv(mapped_markers_and_inds, file = paste0(my_wd, mapped_best_markers_and_inds))

## read vcf file with all sites
vcf_all <- read.vcfR(input.file_vcf, verbose = FALSE )
my_genotypes_all <- extract.gt(vcf_all)
my_genotypes_all <- gsub("|", "", my_genotypes_all, fixed = TRUE)

## save first individual per marker to file #
for(id_p in 1:nrow(mapped_markers_and_inds)){
  # print(mapped_markers_and_inds[id_p])
  predicted_marker <- mapped_markers_and_inds[id_p,"marker"] #
  ind_to_save <- mapped_markers_and_inds[id_p,2] # second column is the first individual with the marker
  print(c("predicted_marker:", predicted_marker))
  ind_geno <- my_genotypes_all[,ind_to_save]
  print(head(ind_geno))
  ## create a genotype file in the format to be read by genomicSimulation #
  ind_geno <- cbind(ind_geno, names(ind_geno))
  colnames(ind_geno) <- c(ind_to_save,"name")
  
  ## update markers names to be as in muts_full format (remove serial number, mut type number and add number to markers in duplicated positions)
  ind_geno <- UpdateMarkerNamesandDuplicants(as.data.frame(ind_geno))

  # write.table(ind_geno[,c(2,1)], file = sprintf("%s%s_updated_markers_names_allele_file_as_txt2.txt",my_wd, ind_to_save), sep = "\t", quote = FALSE, row.names  = FALSE)
  write.table(ind_geno[,c("name",ind_to_save)], file = sprintf("%s%s_updated_markers_names_allele_file_as_txt.txt",my_wd, ind_to_save), sep = "\t", quote = FALSE, row.names  = FALSE)
  
  }
end_time <- Sys.time()
print(c("time (min): ", (end_time-start_time)/60))
