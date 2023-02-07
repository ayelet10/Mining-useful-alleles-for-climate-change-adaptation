### Rscript to prepare files for the breeding simulation with genomicSimulation R package ##
### genetic map file, mutations effect sizes file, tester genotype file as txt ## 

args <- commandArgs(trailingOnly = TRUE) # TRUE
print(args)
seed_num =args[1]  
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
    jump_between_dupl <- 0.00000001
    steps_to_add <- seq(jump_between_dupl,jump_between_dupl*nrow(curr_dupl_rows), by = jump_between_dupl)
    curr_dupl_rows$pos <- as.numeric(curr_dupl_rows$pos) + steps_to_add 
   # replace duplicated rows with corrected positions rows
    markers_mat[rownames(markers_mat) %in% rownames(curr_dupl_rows),] <- curr_dupl_rows 
  }
  
  return(markers_mat)
}

###############################
## prepare effect sizes file ##
###############################
# output format: ## marker, allele, effect size ##
mutations_file <- sprintf("/home/salmanay/projects/Lotterhos_preprint_data/for_JeffRI/%s_Rout_muts_full.txt", seed_num)
seed_mutations <- read.csv(mutations_file, sep = " ")
my_columns <- c("mutname", "mutSalEffect" )
with_mut <- seed_mutations[,my_columns]
with_mut <- cbind(with_mut , rep(1, nrow(with_mut))) # the allele 1
# # add rows for allele 0, (without the mutation), effect size set to 0
with_mut[with_mut$mutSalEffect %in% NA,]$mutSalEffect <- 0
# add positions column # save only positions numbers from marker names#
with_mut$pos <- sub(".*-", "" , with_mut$mutname)
with_mut$pos <- sub("_.*" , "" ,with_mut$pos)
colnames(with_mut) <- c(my_columns, "allele", "pos")
with_mut[,"mutname"] <- gsub(".*-", "", with_mut[,"mutname"]) # update marker names format, remove "-" and digit before it 
print(length(with_mut$pos))
print(length(unique(with_mut$pos)))

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
