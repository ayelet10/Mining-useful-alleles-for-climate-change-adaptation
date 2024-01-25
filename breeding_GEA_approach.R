### breeding scheme for GEA approach #1  ###
### breeding with genomicSimulation package ##

 
## read best_markers_outfile from the GEA script
## read selected individuals and tester genotype files 
start_time <- Sys.time()
print(start_time)
library(genomicSimulation)
library(data.table)


args <- commandArgs(trailingOnly = TRUE) # TRUE
print(args)
seed_num = args[1]
mapped_best_markers_and_inds_file = args[2] # gea_mapped_best_markers_and_inds_outfile.csv
my_wd = args[3]
bc1s4_out_file <- args[4]  ## sprintf("%s_GEMMA_bc1s4_fitness_outfile.txt", seed_num)

all_bc1s4_out_file <- sprintf("%s_all_bc1s4_out_file.txt") # added Dec 7
g0_out_file  <- sprintf("%s_g0_out_file.txt") # added Dec 7
f1_file <- sprintf("%s_f1_file.txt") # added Dec 7
f1_with_marker_file <- sprintf("%s_f1_with_marker_file.txt") # added Dec 7
f1_without_marker_file <- sprintf("%s_f1_without_marker_file.txt") # added Dec 7


# seed_num = 1231227
# seed_num =  1231689
### test - no f1 with marker #
# seed_num = 1231409
# seed_num =1231141
# seed_num =1231277
# seed_num=1233159
set.seed(seed_num) ## added 21 September 2023
###############
## functions ## 
###############
source("/home/salmanay/projects/breeding_simulations/genomicSimulation/scripts/functions_script_breeding.R")
#############

############
## inputs ## 
############
###  markers selected by GEA method; genotypes file
# best_markers_outfile <- "gea_best_markers_outfile.csv"
# mapped_best_markers_and_inds_file <- "gemma_gea_mapped_best_markers_and_inds_outfile.csv" ## change according to method (lfmm2 / gemma)
# my_wd                <- sprintf("/home/salmanay/projects/breeding_simulations/%s/gemma_with_updated_traits_file/",seed_num )## change according to method (lfmm2 / gemma)

####
# input_files_dir <- sprintf("/home/salmanay/projects/breeding_simulations/%s/",seed_num )
input_files_dir <- sprintf("/home/salmanay/projects/breeding_simulations/%s/run3_all_causals/",seed_num )

# mapped_best_markers_and_inds_file <- "gea_mapped_best_markers_and_inds_outfile_July9.csv"
# mapped_best_markers_and_inds_file <- "gemma_gea_mapped_best_markers_and_inds_outfile.csv"
# my_wd                <- sprintf("/home/salmanay/projects/breeding_simulations/%s/",seed_num )
# my_wd <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/"
# my_wd                <- sprintf("C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/test_specific_seed/%s/",seed_num )
# my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/run3_all_causals/gemma_with_updated_traits_file/", seed_num)
setwd(my_wd)
mapped_markers_and_inds <- read.csv(paste0(my_wd, mapped_best_markers_and_inds_file))

bc1s4_offsprings_num <- 100
f1_num_offsprings <- 100
number_of_offsprings_bc1S1_S2_S3 <- 20
# effect.file <- sprintf("%s%s_effect_file_all_markers_as_zero_effect_added_true_effects.txt", my_wd, seed_num) # temp file manually edited
effect.file <- sprintf("%s%s_effect_file_with_neutral_and_causal_sal_run3.txt", input_files_dir, seed_num)
# map.file <- sprintf("%s%s_map_file_muts_full_with_neutral.txt", my_wd, seed_num)
# map.file    <- sprintf("%s%s_map_file_muts_full_with_neutral_fixed.txt", input_files_dir, seed_num)## 15 June 2023. fix of chr names (it caused mis-reading of the wrong rows)
map.file    <- sprintf("%s%s_map_file_muts_full_with_neutral_fixed_run3.txt", input_files_dir, seed_num)## 15 June 2023. fix of chr names (it caused mis-reading of the wrong rows)

# read tester allele file
# tester_geno_file <- sprintf("%s%s_tester_allele_file_muts_full_with_neutral.txt", input_files_dir, seed_num)
tester_geno_file <- sprintf("%s%s_tester_allele_file_muts_full_with_neutral_22.txt", input_files_dir, seed_num)

# tester_geno_file <- sprintf("%s_tester_allele_file_muts_full.txt", seed_num)
# tester_geno_file <- sprintf("%s_tester_allele_file_10_snps.txt", seed_num)
elite_line_id <- "tester"

## for testing - files with the same number of markers as the individual txt file. 
# tester_geno_file <- sprintf("%s_tester_allele_file.txt", seed_num)
# effect.file <- sprintf("%s%s_effect_file_sal.txt", my_wd, seed_num)
# map.file <- sprintf("%s%s_map_file.txt", my_wd, seed_num)

## test seed 1231689 - file with "11" at 50005, 100004 ## 
##
# tester_geno_file <- sprintf("%s_tester_with_markers_at_chr_ends.txt",  seed_num)
# effect.file <- sprintf("1231409_effect_file_test.txt")
#####################
### delete out file if already exists ### 
#Check its existence
if (file.exists(bc1s4_out_file)) {
  #Delete file if it exists
  file.remove(bc1s4_out_file)
}
################
### breeding ###
################
## load the data into genomicSimulation ###
all_bc1S4_results <- c()   ## collect verified markers (markers that passed the t.test)
# tests with 10 bp 
# map.file  <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/my_data_to_qtl2/test_3_subset_10_sites/1231227_map_file_muts_full.txt"
# tester_geno_file <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/my_data_to_qtl2/test_3_subset_10_sites/1231227_tester_allele_file_muts_full.txt"
# effect.file  <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/my_data_to_qtl2/test_3_subset_10_sites/1231227_effect_file_sal.txt"
# curr_marker_original_name <- 95
# predicted_lr_txt_file <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/my_data_to_qtl2/test_3_subset_10_sites/tsk_488_updated_markers_names_allele_file_as_txt.txt"
# ###
# map.file  <- "/home/salmanay/projects/breeding_simulations/genomicSimulation/temp_tests/test_tsk_488_10_sites/1231227_map_file_muts_full.txt"
# tester_geno_file <- "/home/salmanay/projects/breeding_simulations/genomicSimulation/temp_tests/test_tsk_488_10_sites/1231227_tester_allele_file_muts_full.txt"
# effect.file  <- "/home/salmanay/projects/breeding_simulations/genomicSimulation/temp_tests/test_tsk_488_10_sites/1231227_effect_file_sal.txt"
# curr_marker_original_name <- 95
# predicted_lr_txt_file <- "/home/salmanay/projects/breeding_simulations/genomicSimulation/temp_tests/test_tsk_488_10_sites/tsk_488_updated_markers_names_allele_file_as_txt.txt"
# curr_ind_id <- "tsk_488"

for(l in 1:nrow(mapped_markers_and_inds)){ # loop over markers
  bc1s4_res <- c()
  all_bc1s4_res <- c()#all genotypes # added Dec 7
  g0_res <- c()   # added Dec 7
  f1_res <- c()   # added Dec 7
  f1_with_marker_res <- c() # added Dec 7
  f1_without_marker_res <- c() # added Dec 7
  
  fitness_homozygous_bc1S4 <- c()
  # per selected predicted individual: #
  curr_ind_id <- mapped_markers_and_inds[l,]$individual1 ## take only the first individual. column2 (there could be more than one)
  curr_marker_original_name <- mapped_markers_and_inds[l,]$marker
  print(c("curr_marker_original_name:", curr_marker_original_name))
 # curr_marker<- gsub("_c", "", curr_marker_original_name)  # Nov 12 - the marker names are in the format of pos_2_c (in case of two mutations in the same pos and it is causal)
  curr_marker <-  curr_marker_original_name  ## not sure if GetNameInMutsFullFormat step is needed anymore, looks like the marker names are in the 'muts_full" format e.g the format in the map file
  f1_geno_file <- sprintf("%s_%s_f1_genotypes_gea_breeding.txt", seed_num, curr_ind_id)
  print(c("individual:", curr_ind_id))
  print(c("marker:", curr_marker))
  curr_marker_ind_in_map <- findMarkerIndexInMap(map.file, curr_marker)
  print(c("curr_marker_ind_in_map" ,curr_marker_ind_in_map))
  # curr_ind_id <- "tsk_84"
  # curr_ind_id <- "tsk_488"
  
  predicted_lr_txt_file <- sprintf("%s%s_updated_markers_names_allele_file_as_txt_GEA.txt",my_wd, curr_ind_id)
  ## load tester data
  g0 <- genomicSimulation::load.data(tester_geno_file, map.file = map.file, effect.file = effect.file ) # load tester
  # g0 <- genomicSimulation::load.data(predicted_lr_txt_file, map.file = map.file, effect.file = effect.file ) # load tester
  # g0_geno_file <- sprintf("%s_%s_g0_genotypes_gea_breeding.txt", seed_num, curr_ind_id)
  # save.genotypes(filename = "g0_geno_file_test_2.txt", group = g0, type= "T")
  
  ## add predicted individual data ##
  g0 <- genomicSimulation::load.more.genotypes(predicted_lr_txt_file)
  g0_info <- data.frame(Index=see.group.data(g0,data.type ="X"), GEBV=see.group.data(g0,data.type = "BV"), G=see.group.data(g0,data.type ="G"))
  fitness_g0    <- fitness_function(g0_info$GEBV)     # added Dec 7
  g0_res <- cbind(fitness_g0, curr_ind_id, curr_marker) # added Dec 7
  
  print(sprintf("selected individual %s GEBV: %s", curr_ind_id, g0_info$GEBV))
  print(sprintf("selected individual %s calculated fitness: %s", curr_ind_id, fitness_function(g0_info$GEBV)))

# head(g0_info$G)
  print("g0 at marker:")
  g0_at_marker <-substr((g0_info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)
  print(g0_at_marker)
  if("00" %in% g0_at_marker){
    f1_marker_error_message <- sprintf("%s marker name. genotype at marker is 00", curr_marker)
    print(f1_marker_error_message)
    write(f1_marker_error_message, file = "marker_error_in_f1", append = TRUE)
    clear.simdata()
    break
  }
  
  
  f1 <- cross.combinations(first.parents = elite_line_id, second.parents = curr_ind_id, give.names = TRUE, name.prefix = "f1", save.genotype= FALSE, offspring = f1_num_offsprings)
  save.genotypes(filename = f1_geno_file, group = f1)
  # ### reading f1 from file ##
  # f1_geno_from_file <- read.table(f1_geno_file, sep = "\t", header = TRUE,  colClasses = "character",check.names = FALSE)
  # rownames(f1_geno_from_file) <- f1_geno_from_file[,1]
  # f1_geno_from_file <- f1_geno_from_file[,-1]
  # t_f1_geno_from_file <- t(f1_geno_from_file)
  # t_f1_geno_from_file[which(t_f1_geno_from_file[rownames(t_f1_geno_from_file) %in% curr_marker,] %in% c("01","10"))]
  # marker_t_f1_geno_from_file <- t_f1_geno_from_file[rownames(t_f1_geno_from_file) %in% curr_marker,]
  # names_individuals_with_marker <- names(marker_t_f1_geno_from_file[marker_t_f1_geno_from_file %in% c("01","10")][1:10])
  
  ## check if marker is significantly good ## compare individuals with the marker to the tester # 
  f1_info <- getGroupInfoAsDF(f1) #data.frame(Index=see.group.data(f1,data.type ="X"), GEBV=see.group.data(f1,data.type = "BV"), G=see.group.data(f1,data.type ="G"))
  fitness_f1    <- fitness_function(f1_info$GEBV)     # added Dec 7
  f1_res <- cbind(fitness_f1, curr_ind_id, curr_marker) # added Dec 7
  
  print("f1 at marker:")
  f1_at_marker<- substr((f1_info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)
  print(f1_at_marker)
  
  print(c("curr_marker_ind_in_map", curr_marker_ind_in_map))
  f1_inds_with_marker    <- getIndsWithMarkerFromInfo(f1_info, curr_marker_ind_in_map, "with_marker")    #f1_info[f1_info$G %in% c("11","01","10"),]
  print("summary f1_inds_with_marker:")
  summary(f1_inds_with_marker$GEBV)
  print(substr(f1_inds_with_marker$G[1:10], curr_marker_ind_in_map*2-20, curr_marker_ind_in_map*2+20))
  f1_inds_without_marker <- getIndsWithMarkerFromInfo(f1_info, curr_marker_ind_in_map, "without_marker") #f1_info[f1_info$G %in% c("00"),]
  print("summary f1_inds_without_marker:")
  summary(f1_inds_without_marker$GEBV)
  print(substr(f1_inds_without_marker$G[1:10], curr_marker_ind_in_map*2-20, curr_marker_ind_in_map*2+20))
  
  fitness_f1_with_marker    <- fitness_function(f1_inds_with_marker$GEBV)     # added Dec 7
  f1_with_marker_res <- cbind(fitness_f1_with_marker, curr_ind_id, curr_marker) # added Dec 7
  fitness_f1_without_marker    <- fitness_function(f1_inds_without_marker$GEBV)     # added Dec 7
  f1_without_marker_res <- cbind(fitness_f1_without_marker, curr_ind_id, curr_marker) # added Dec 7
  
  # print(f1_hetero_individuals)
  # print(nrow(f1_hetero_individuals))
  if(nrow(f1_inds_with_marker) < 1 ){
    f1_error_message <- sprintf("%s marker name. no individuals with the marker found in f1. Move to next marker", curr_marker)
    print(f1_error_message)
    write(f1_error_message, file = "error_in_f1", append = TRUE)

    clear.simdata()
    next
  }else{
    ## to add a check that both groups (with vs. without the marker) have the same size ? or use average ##
    # avg_fitt_inds_with_marker <- mean(fitt_inds_with_marker) / length(fitt_inds_with_marker)
    # avg_fitt_inds_without_marker <- mean(fitt_inds_without_marker) / length(fitt_inds_without_marker)
    if(nrow(f1_inds_with_marker) < 10 | nrow(f1_inds_without_marker) < 10){
      num_f1_inds_to_compare <- min(nrow(f1_inds_with_marker), nrow(f1_inds_without_marker))
    }else{
      num_f1_inds_to_compare <- 10
    }
    ## get the effect sizes of the selected individuals ### to make sure that we take marker that has a positive effect size
    effect_sizes_with_marker    <- (f1_inds_with_marker$GEBV)[1:num_f1_inds_to_compare]  #to take the 1-10 individuals
    effect_sizes_without_marker <- (f1_inds_without_marker$GEBV)[1:num_f1_inds_to_compare] #to take the 1-10 individuals
    print("effect_sizes_without_marker:")
    print(effect_sizes_without_marker)
    print("effect_sizes_with_marker:")
    print(effect_sizes_with_marker)
    print(c("num_f1_inds_to_compare:", num_f1_inds_to_compare))
    ## test to make sure that the mean effect_sizes_with_marker is higher that the mean effect_sizes_without_marker
    ## if it fails - do not continue with this marker ##
    if(mean(effect_sizes_with_marker) > mean(effect_sizes_without_marker)){
      print("mean effect_sizes_with_marker is higher that the mean effect_sizes_without_marker")
    }else{
      print("mean effect_sizes_with_marker is LOWER that the mean effect_sizes_without_marker !! skipping to next marker")
      next
    }
    ## compare fitness ##
    fitt_inds_with_marker    <- fitness_function(effect_sizes_with_marker) # the updated function for the new optimum environment. 0.82 is the new optimum
    fitt_inds_without_marker <- fitness_function(effect_sizes_without_marker)
    # tester_fitt              <- fitness_function(tester_effect_size)
    print(c("fitt_inds_with_marker:", fitt_inds_with_marker))

    ### compare fitness with vs. without the marker ##
    # fitt_t_test_res <-  t.test(fitt_inds_with_marker, fitt_inds_without_marker, var.equal=TRUE, alternative = "greater") # t.test(x,y,var.equal=TRUE)
    fitt_t_test_res <- try(t.test(fitt_inds_with_marker, fitt_inds_without_marker, var.equal=TRUE, alternative = "greater"))
    if(grepl(pattern = "Error", x = fitt_t_test_res[1])){
      message <- sprintf("%s. error in t-test. Skipping to next marker ", curr_marker)
      print(message)
      write(message, file = "error_in_t_test", append = TRUE)
      write(fitt_t_test_res, file = "error_in_t_test", append = TRUE)
      clear.simdata()
      next
    } else {
      fitt_t_test_res
    }

    # print(t_test_res) # effect sizes
    print(fitt_t_test_res) # fitness
    # if(exists(fitt_t_test_res)){

      if(fitt_t_test_res$p.value <= 0.05){
        print("found good marker")
        # curr_res <- c(curr_ind_id, max(f1_inds_with_marker$GEBV))
        curr_res <- c(curr_ind_id, max(fitt_inds_with_marker))
        # save one f1 descendant to new group (a heterozygous is currently taken)
        f1_hetero_individuals <- getIndsWithMarkerFromInfo(f1_info, curr_marker_ind_in_map, "hetero") ##
        name_f1_to_cross <- f1_hetero_individuals$Name[1]  ## take the first f1 heterozygote

        ### bc1. second parent is the F1 ### 1 offspring ###
        bc1 <- cross.combinations(first.parents = elite_line_id ,second.parents = name_f1_to_cross ,give.names = TRUE, name.prefix = "bc1")
        bc1_group_info <- getGroupInfoAsDF(bc1)
        bc1_geno_at_marker <- substr((bc1_group_info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)
        print(substr((bc1_group_info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2))

        while(!bc1_geno_at_marker %in% c("21", "12","02","20") ){ ##  until I get a heterozygote  ## changed Jan-25-24, was:  c("10", "01")
          delete.group(bc1)
          print("bc1 deleted")
          bc1 <- cross.combinations(first.parents = elite_line_id ,second.parents = name_f1_to_cross, give.names = TRUE, name.prefix = "bc1")
          bc1_group_info <- getGroupInfoAsDF(bc1)
          bc1_geno_at_marker <- substr((bc1_group_info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)
          print(bc1_geno_at_marker)
        }

        print("before bc1s1")
        bc1S1 <- ChooseIndividualByStateAndSelf(bc1_group_info, curr_marker_ind_in_map, number_of_offsprings_bc1S1_S2_S3, "hetero")

        print(see.existing.groups())
        bc1S1_group_info <- getGroupInfoAsDF(bc1S1)
        #  delete.group(bc1)
        bc1S2 <- ChooseIndividualByStateAndSelf(bc1S1_group_info, curr_marker_ind_in_map, number_of_offsprings_bc1S1_S2_S3, "hetero")
        print(see.existing.groups())
        
        bc1S2_group_info <- getGroupInfoAsDF(bc1S2)
        # #  delete.group(bc1S1)
        bc1S3 <- ChooseIndividualByStateAndSelf(bc1S2_group_info, curr_marker_ind_in_map, number_of_offsprings_bc1S1_S2_S3, "hetero")
        print(see.existing.groups())
        
        bc1S3_group_info <- getGroupInfoAsDF(bc1S3)
        # #  delete.group(bc1S2)
        bc1S4 <- ChooseIndividualByStateAndSelf(bc1S3_group_info, curr_marker_ind_in_map, bc1s4_offsprings_num, "hetero") ## 100 offsprings
        print(see.existing.groups())
        
        bc1S4_group_info <- getGroupInfoAsDF(bc1S4)
        fitness_bc1S4    <- fitness_function(bc1S4_group_info$GEBV)     # added Dec 7
        all_bc1s4_res <- cbind(fitness_bc1S4, curr_ind_id, curr_marker) # added Dec 7
        ## save info of bc1s4 homozygous to the marker ##
        bc1S4_homozygous_to_marker <- getIndsWithMarkerFromInfo(bc1S4_group_info, curr_marker_ind_in_map, "homo") ## >> no heterozygous found in this stage when I tested with tsk_488, marker 133127
        new_group_bc1S4_homozygous_to_marker <- make.group(bc1S4_homozygous_to_marker$Index)
        print(see.existing.groups())
        
        delete.group(bc1S4)
        homozygous_bc1S4_group_info <- getGroupInfoAsDF(new_group_bc1S4_homozygous_to_marker)
        fitness_homozygous_bc1S4    <- fitness_function(homozygous_bc1S4_group_info$GEBV)

        bc1s4_res <- cbind(fitness_homozygous_bc1S4, curr_ind_id, curr_marker)


      }else{
        message <- sprintf("f1 fitness comparison - the marker %s is not significant. Skipping to next marker ", curr_marker)
        print(message)
        clear.simdata()
        next
      }
    # }else{
    #   message <- sprintf("error in t-test. Skipping to next marker ", curr_marker)
    #   print(message)
    #   clear.simdata()
    #   next
    # }
  }
  write.table(bc1s4_res, file = bc1s4_out_file, append = TRUE, col.names = F, row.names = F)
  write.table(all_bc1s4_res, file = all_bc1s4_out_file, append = TRUE, col.names = F, row.names = F) # added Dec 7
  write.table(g0_res, file = g0_out_file, append = TRUE, col.names = F, row.names = F) # added Dec 7
  write.table(f1_res, file = f1_out_file, append = TRUE, col.names = F, row.names = F) # added Dec 7
  write.table(f1_with_marker_res, file = f1_with_marker_file, append = TRUE, col.names = F, row.names = F) # added Dec 7
  write.table(f1_without_marker_res, file = f1_without_marker_file, append = TRUE, col.names = F, row.names = F) # added Dec 7

  write.table(f2_homo_to_marker_res, file = f2_homo_to_marker_file, append = TRUE, col.names = F, row.names = F) # added Dec 7
  write.table(f2_hetero_to_marker_res, file = f2_hetero_to_marker_file, append = TRUE, col.names = F, row.names = F) # added Dec 7
  write.table(f2_without_marker_res, file = f2_without_marker_file, append = TRUE, col.names = F, row.names = F) # added Dec 7
  
  write.table(f2_homo_to_marker_res, file = f2_homo_to_marker_file, append = TRUE, col.names = F, row.names = F) # added Dec 7
  write.table(f2_hetero_to_marker_res, file = f2_hetero_to_marker_file, append = TRUE, col.names = F, row.names = F) # added Dec 7
  write.table(f2_without_marker_res, file = f2_without_marker_file, append = TRUE, col.names = F, row.names = F) # added Dec 7

  write.table(f2_homo_to_marker_res, file = f2_homo_to_marker_file, append = TRUE, col.names = F, row.names = F) # added Dec 7
  write.table(f2_hetero_to_marker_res, file = f2_hetero_to_marker_file, append = TRUE, col.names = F, row.names = F) # added Dec 7
  write.table(f2_without_marker_res, file = f2_without_marker_file, append = TRUE, col.names = F, row.names = F) # added Dec 7
  
  
  clear.simdata() ## needs to clear previous simulation data to run a new one
}