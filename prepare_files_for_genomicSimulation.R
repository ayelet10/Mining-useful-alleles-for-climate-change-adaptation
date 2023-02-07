### Rscript to prepare files for the breeding simulation with genomicSimulation R package ##
### genetic map file, mutations effect sizes file, tester genotype file as txt ## 

args <- commandArgs(trailingOnly = TRUE) # TRUE
print(args)
seed_num =args[1]  
# seed_num =1231227  #1231214 ##  #
my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/",seed_num )

dir.create(file.path(my_wd)) # if dir exists, it just gives a warning
setwd(file.path(my_wd))
library(vcfR)

#################
### functions ###
#################
ChangeDuplicatedPositions <- function(markers_mat) { 
  ## column with positions to change should be "pos" ##
  ## check for duplicated markers, by position
  ## per duplicated site, add 1/ number of duplicants of a site to get different adjacent bp sites ##
  duplicated_loci <- markers_mat[duplicated(markers_mat$pos),] 
  unique_dupicated_sites <- unique(duplicated_loci$pos)
  for(u in unique_dupicated_sites){
    curr_dupl_rows <- duplicated_loci[duplicated_loci$pos %in% u, ]
    # jump_between_dupl <- 1/(nrow(curr_dupl_rows)+ 1) ## added 1 to account for the first occurence # the jump changes according to the number of dupicants on the same position
    jump_between_dupl <- 0.00000001
    steps_to_add <- seq(jump_between_dupl,jump_between_dupl*nrow(curr_dupl_rows), by = jump_between_dupl)
    curr_dupl_rows$pos <- as.numeric(curr_dupl_rows$pos) + steps_to_add 
    
    markers_mat[rownames(markers_mat) %in% rownames(curr_dupl_rows),] <- curr_dupl_rows # replace duplicated rows with corrected positions rows
  }
  
  return(markers_mat)
}
##############################
## prepare genomic map file ## >> read the vcf and add unique chromosome number to every 50kb block. File header need to be: marker	chr	pos 
###################################################
# marker - the snp name                           #
# pos - the position on the "chromosome"          #
# chr - 1-20 linkage groups as chromosomes 1-20   #
###################################################
# full_vcf_file <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_plusneut_MAF01.recode2.vcf.gz", seed_num) # the tester is not dependent on the seed 
# full_vcf <- read.vcfR(full_vcf_file, verbose = FALSE )
## copy the first 3 columns to file for a genetic map creation  >> to try when farm is available
# vcf-subset -c name1,name2,name3 /home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/1231227_plusneut_MAF01.recode2.vcf.gz | tr "\t" "," > first_three_columns_vcf_out.csv 


# number of rows in vcf. 6 first rows are the header ### marker names are not read completly
# zcat /home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/1231227_plusneut_MAF01.recode2.vcf.gz | wc -l

# my_map <- read.csv(sprintf("%s.first_two_columns_recode_vcf_out.csv", seed_num) , skip = 5)
# my_map <- read.csv(sprintf("%s.first_two_columns_recode_vcf_out_15Jan.csv", seed_num) , skip = 5, sep = "\t")
############## 
# input.file_vcf <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_plusneut_MAF01.recode2.vcf.gz", seed_num) # non variable sites error ?
# input.file_vcf <- sprintf("/home/salmanay/projects/breeding_simulations/%s/%s_hundred_over_hundred_recode_vcf_out_mac.vcf.recode.vcf" , seed_num,seed_num)   # 76 sites kept after vcftools -mac 1 for non variable sites removal

# vcf <- read.vcfR(input.file_vcf, verbose = FALSE )
# my_genotypes <- extract.gt(vcf)
# print(head(my_genotypes[1:10,1:10]))
# print("row names:")
# print(rownames(my_genotypes[1:10,1:10]))
# marker_names <- rownames(my_genotypes)


###############################
## prepare effect sizes file ##
###############################
# output format: ## marker, allele, effect size ##
# # mutations_file <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_muts.txt", seed_num)
mutations_file <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_Rout_muts_full.txt", seed_num)
# mutations_file <- sprintf("C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/%s_Rout_muts_full_300_rows.txt", seed_num)
# # mutations_file <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/data_from_Lotterhos_preprint_2022/1231227_muts.txt"

seed_mutations <- read.csv(mutations_file, sep = " ")
# write.table( seed_mutations[1:300,], file = sprintf("%s_Rout_muts_full_300_rows.txt", seed_num), sep = " ")
# print(head(seed_mutations))
# print(unique(seed_mutations["mutname"]))

my_columns <- c("mutname", "mutSalEffect" )
with_mut <- seed_mutations[,my_columns]

with_mut <- cbind(with_mut , rep(1, nrow(with_mut))) # the allele 1
# # add rows for allele 0, (without the mutation), effect size set to 0
# without_mut <- cbind(seed_mutations[c("mutname" )] , rep(0, nrow(with_mut)),rep(0, nrow(with_mut)) )
with_mut[with_mut$mutSalEffect %in% NA,]$mutSalEffect <- 0
# add positions column # save only positions numbers from marker names#
print(head(with_mut))
with_mut$pos <- sub(".*-", "" , with_mut$mutname)
with_mut$pos <- sub("_.*" , "" ,with_mut$pos)
colnames(with_mut) <- c(my_columns, "allele", "pos")
with_mut[,"mutname"] <- gsub(".*-", "", with_mut[,"mutname"]) # update marker names format, remove "-" and digit before it 
print(length(with_mut$pos))
print(length(unique(with_mut$pos)))


# print(with_mut[!with_mut$mutname %in% 1,])
# colnames(without_mut) <- c(my_columns, "allele")
# all_muts <- rbind(with_mut, without_mut, deparse.level = 0  )
# write.table(with_mut[c("pos","allele","mutSalEffect")], file = sprintf("%s%s_effect_file_sal.txt",my_wd,seed_num ) , sep = " ", quote = FALSE,  col.names = FALSE, row.names = FALSE)

write.table(with_mut[c("mutname","allele","mutSalEffect")], file = sprintf("%s%s_effect_file_sal.txt",my_wd,seed_num ) , sep = " ", quote = FALSE,  col.names = FALSE, row.names = FALSE)
###############################################################################
## prepare tester genotype file for breeding schemes with genomicSimulation  ##
###############################################################################
number_sites <-  nrow(with_mut) 
rnames <-  with_mut$pos # with_mut$mutname # my_map$marker #  c(rep(1:10))# genes names from full vcf file pos #
cnames <-  "tester"
tester_geno <- matrix(rep("00" ,number_sites), ncol = 1, nrow = number_sites, dimnames=list(rnames,cnames))
tester_geno <- cbind(tester_geno, with_mut$mutname)   


tester_geno <- tester_geno[,c(2,1)]
colnames(tester_geno) <- c("name",cnames )
# tester_geno[,"name"] <- gsub("-", "_", tester_geno[,"name"]) # update marker names format, replace "-" with "_"
write.table(tester_geno, file = sprintf("%s_tester_allele_file_muts_full.txt", seed_num), sep = "\t", quote = FALSE, row.names = FALSE)

########################################################################
## prepare genomic map file #                                         ## 
## read the vcf and add unique chromosome number to every 50kb block.  #
## convert bp to centiMorgan by dividing in 1000                       #
## File header needs to be: marker	chr	pos                            #
########################################################################
## currently it holds duplicated sites  ## markers from SLiM that are at the same position but acted differently over time #
my_map <- with_mut[,c("allele", "pos", "mutname")]
colnames(my_map) <- c("chr", "pos","marker")
print("head my_map:")
print(head(my_map))
 ## the marker is the original position (before division to chromosomes)
my_map$chr <- sub("_.*", "", my_map$marker)

# ## add chromosome numbers - every 50kb is a new "chromosome" #
chr_start <- 1 # 1
chr_end   <- 50000 # 50000
chr_size  <- 50000

num_linkage_groups <- 20  # number of linkage groups (chromosomes)

print(length(my_map$pos))
print(length(unique(my_map$pos)))

my_map$pos <- as.numeric(my_map$pos)
my_map$pos_total <- as.numeric(my_map$pos)
# my_map$marker <- as.numeric(my_map$marker)
for(chr_num in c(1:num_linkage_groups)){ # loop to update the chr names and markers
  print(chr_num)
  my_map[(chr_start < my_map$pos & my_map$pos < chr_end), ]$chr <- chr_num #  updating the chr name
  
  # update the position by reducing 50kb, starting from chr2
  if(chr_num > 1){
      my_map[my_map$chr %in% chr_num, ]$pos <- (my_map[my_map$chr %in% chr_num, ]$pos_total - (chr_size*(chr_num-1)))
  }
  chr_start <- chr_start + chr_size
  chr_end   <- chr_end + chr_size
  print(chr_start)
  print(chr_end)
}
## change bp to centiMorgan # 0.001 between proximate base #
my_map$pos <- my_map$pos*0.001
## Per chromosome, Treat duplications, by position. Per duplicated site, add 1/ number of duplicants*1000000# 
for (chr in 1:num_linkage_groups){
  print(chr)
  curr_chr <- my_map[my_map$chr %in% chr,]
  my_map[my_map$chr %in% chr,] <- ChangeDuplicatedPositions(curr_chr)
}
write.table(my_map[,c("marker", "chr","pos" )], file= sprintf("%s_map_file_muts_full.txt", seed_num) , sep = "\t", quote = FALSE, row.names = FALSE)
#############################
## work with previous files #
#############################
# cnames <-  "allele"
# ## the mutations with effect are listed in %s_muts.txt file, all the other have 0 effect. ## not needed if using _Rout_muts_full. file (?)
# # all_markers_one <- matrix(rep("1" ,number_sites), ncol = 1, nrow = number_sites, dimnames=list(rnames))
# all_markers_one <- matrix(rep("1" ,number_sites), ncol = 1, nrow = number_sites, dimnames=list(rnames)) # alleles
# all_markers_one <- cbind(all_markers_one, my_map$marker)
# all_markers_one <- cbind(all_markers_one, 0) ## effects
# all_markers_one[all_markers_one[,2]%in% with_mut$mutID, 3] <- with_mut$mutSalEffect # in case there is a snp with effect size it will be assigned the effect size from the sim
# ##
# all_markers_zero <- matrix(rep("0" ,number_sites), ncol = 1, nrow = number_sites, dimnames=list(rnames))
# all_markers_zero <- cbind(all_markers_zero, my_map$marker)
# all_markers_zero <- cbind(all_markers_zero, 0) ## effects
# 
# all_markers <- rbind(all_markers_one,all_markers_zero )
# head(all_markers)
# write.table(all_markers[,c(2,1,3)], file = sprintf("%s%s_effect_file_all_markers.txt", my_wd, seed_num) , sep = " ", quote = FALSE,  col.names = FALSE, row.names = FALSE)
# 
############## general tester - not sure genomicSimulation can use it ###############
# number_sites <-  num_linkage_groups*chr_size #10 # to take from the full_vcf
# rnames <-  c(1:number_sites) #  c(rep(1:10))# genes names from full vcf file pos # 
# cnames <-  "tester"
# tester_geno <- matrix(rep("00" ,number_sites), ncol = 1, nrow = number_sites, dimnames=list(rnames,cnames))
# tester_geno <- cbind(tester_geno, rnames)
# head(tester_geno)
# tester_geno <- tester_geno[,c(2,1)]
# colnames(tester_geno) <- c("name",cnames )
# write.table(tester_geno, file = sprintf("tester_allele_file.txt", seed_num), sep = "\t", quote = FALSE, row.names = FALSE)

# ## check consistency of mutations between files ##
# eff_file <- read.table(sprintf("%s_effect_file.txt", seed_num ),sep = " " )
# eff_file <- read.table(sprintf("%s_effect_file_muts_only.txt", seed_num ),sep = " " )
# map_file <-read.table(sprintf("%s_map_file.txt", seed_num),sep = "\t" , header = TRUE)
# alleles <- read.table(sprintf("%s%s_updated_markers_names_allele_file_as_txt.txt",my_wd, pred_lr ),sep = "\t" , header = TRUE, colClasses=c('character','character'))
# 
# alleles[alleles$name %in% map_file$pos,] #
# eff_file[eff_file$V1 %in% map_file$pos,] # 0
# map_file[as.numeric(map_file$pos) %in% as.numeric(eff_file$V1) ,] # 0
# alleles[alleles$name %in% eff_file$V1,]