### breeding scheme for GEA approach #1  ###
### breeding with genomicSimulation package ##

## read data from vcf 
## read best_markers_outfile from the GEA script

start_time <- Sys.time()
print(start_time)
library(vcfR)
library(genomicSimulation)
library(data.table)


args <- commandArgs(trailingOnly = TRUE) 
print(args)
seed_num = args[1]
###############
## functions ## 
###############
findMarkerIndexInMap <- function(map_file, marker_name){
  ## get the index of a marker in the map file ##
  map_file <- map.file
  my_map <- read.table(map_file, header = TRUE)
  marker_ind <- which(my_map$marker %in% marker_name)
  return(marker_ind)
}

getMarkerEffectSize <- function (marker_name, effect_file_path){
  print(marker_name)
  effect_file <- read.csv(effect_file_path, sep = " ", header = FALSE)
  marker_effect_size <- effect_file[as.character(effect_file$V1) %in% marker_name,]$V3
  return(marker_effect_size)
}

getIndsWithMarkerFromInfo <- function(group_info_data, curr_marker_index, state){
  ## states : "homo", "hetero", "without_marker", "with_marker" ##
  if(state == "homo"){
    states_vec <- c("11")
  }else if(state == "hetero"){
    states_vec <- c("01","10")
  }else if (state == "without_marker"){
    states_vec <- c("00")
  }else if (state == "with_marker"){
    states_vec <- c("01","10", "11")
  }
  # get the marker genotypes per individuals in the group #
  group_info_data$G <- substr((group_info_data$G), curr_marker_index*2-1, curr_marker_index*2)
  group_inds_with_marker <- group_info_data[group_info_data$G %in% states_vec,]
  
  return(group_inds_with_marker)
}

SelfIndividualWithMarker <- function(breeding_group_info, number_of_offsprings, index_of_marker_in_map){
  ## from a group, take individual that has the marker (homo or hetero), save it to a new group and self. Return the offsprings group.
    # group_info <- data.frame(Index=see.group.data(breeding_group,data.type ="X"), GEBV=see.group.data(breeding_group,data.type = "BV"), G=see.group.data(breeding_group,data.type ="G"))
  with_marker_individuals <- getIndsWithMarkerFromInfo(breeding_group_info, index_of_marker_in_map, "with_marker")
  print(with_marker_individuals)
  new_group <- make.group(with_marker_individuals[1,1])   
  bc_group <- self.n.times(new_group, n = 1, offspring = number_of_offsprings) 
  return (bc_group)
}

GetNameInMutsFullFormat <- function(name_of_marker_from_vcf, path_to_names_key_file){
  # names_keys <- read.table("key_marker_names_vcf_to_muts_all.txt", header = TRUE )
  names_keys <-  read.table(path_to_names_key_file, header = TRUE )
  print(head(names_keys))
  name_of_marker_as_in_muts_all <- names_keys[names_keys$original_names %in% name_of_marker_from_vcf, ]$name
  return(name_of_marker_as_in_muts_all)
}
#############
## outputs ##
#############
best_individuals_outfile     <- "gea_best_individuals_outfile.csv"
mapped_best_markers_and_inds <- "gea_mapped_best_markers_and_inds_outfile.csv"
out_file                     <- "GEA_approach_good_markers_bc1S4_results.csv"
good_markers_out_file <- "GEA_approach_good_markers.csv"
############
## inputs ## 
############
###  markers selected by GEA method; genotypes file
best_markers_outfile <- "gea_best_markers_outfile.csv"
mapped_best_markers_and_inds_file <- "gea_mapped_best_markers_and_inds_outfile.csv"
my_wd                <- sprintf("/home/salmanay/projects/breeding_simulations/%s/",seed_num )
setwd(my_wd)
mapped_markers_and_inds <- read.csv(paste0(my_wd, mapped_best_markers_and_inds_file))
bc1s4_offsprings_num <- 100
f1_num_offsprings <- 100
number_of_offsprings_bc1S1_S2_S3 <- 20
effect.file <- sprintf("%s%s_effect_file_sal.txt", my_wd, seed_num)
map.file <- sprintf("%s%s_map_file_muts_full.txt", my_wd, seed_num)
# read tester allele file
tester_geno_file <- sprintf("%s_tester_allele_file_muts_full.txt", seed_num)
elite_line_id <- "tester"
#####################
## breeding scheme ##
#####################
## load the data into genomicSimulation ###
varified_good_markers <- c()
all_bc1S4_results <- c()   ## collect verified markers (markers that passed the t.test)
for (l in 1:nrow(mapped_markers_and_inds) ){ # loop per marker
  t_test_res = c()
  # per selected predicted individual: #
  curr_ind_id <- mapped_markers_and_inds[l,]$individual1 ## take only the first individual. column2 (there could be more than one)
  curr_marker_original_name <- mapped_markers_and_inds[l,]$marker
  path_to_names_keys_file <- paste0(my_wd, "key_marker_names_vcf_to_muts_all.txt")
  print(path_to_names_keys_file)
  curr_marker <- GetNameInMutsFullFormat(curr_marker_original_name, path_to_names_keys_file) # update names. keep position and consider duplicants.
 
  print(curr_ind_id)
  print(curr_marker)
  curr_marker_ind_in_map <- findMarkerIndexInMap(tester_geno_file, curr_marker)
  print(curr_marker_ind_in_map)
  predicted_lr_txt_file <- sprintf("%s%s_updated_markers_names_allele_file_as_txt.txt",my_wd, curr_ind_id )
  
  g0 <- genomicSimulation::load.data(tester_geno_file, map.file = map.file, effect.file = effect.file ) # load tester
    ## tester's effect size ##
  tester_effect_size <- see.group.data(g0, data.type = "BV")
  ## add predicted individual data ##
  g0 <- genomicSimulation::load.more.genotypes(predicted_lr_txt_file)
  f1 <- cross.combinations(first.parents = elite_line_id, second.parents = curr_ind_id, give.names = TRUE, name.prefix = "f1", save.genotype= TRUE, offspring = f1_num_offsprings)
  print(see.group.data(group = f1, data.type = "PED"))
  ## check if marker is significantly good ## compare individuals with the marker to the tester # 
  f1_info <- data.frame(Index=see.group.data(f1,data.type ="X"), GEBV=see.group.data(f1,data.type = "BV"), G=see.group.data(f1,data.type ="G"))
  print(curr_marker_ind_in_map)
  f1_info$G <- substr((f1_info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2) ##indices of the selected marker (the row index minus one, the row index)
  f1_inds_with_marker <- f1_info[f1_info$G %in% c("11","01","10"),]
  f1_inds_without_marker <- f1_info[f1_info$G %in% c("00"),]

  if(nrow(f1_inds_with_marker)< 1  ){
    print("no individuals with the marker found in f1. Move to next marker")
    next
  }
  ## get the effect sizes of the selected individuals ###  >> currently it takes the 1:10 individuals
  effect_sizes_with_marker <- (f1_inds_with_marker$GEBV)[1:10]
  effect_sizes_without_marker <- (f1_inds_without_marker$GEBV)[1:10]
  print("effect_sizes_with_marker:")
  print(effect_sizes_with_marker)
  ## compare fitness ##
  fitt_inds_with_marker <- exp(-1/2*(effect_sizes_with_marker-0.82)^2/.5) # the updated function for the new optimum environment. 0.82 is the new optimum
  fitt_inds_without_marker <- exp(-1/2*(effect_sizes_without_marker-0.82)^2/.5)
  tester_fitt <- exp(-1/2*(tester_effect_size-0.82)^2/.5)
  print("fitt_inds_with_marker:")
  print(fitt_inds_with_marker)
  ## to add a check that both groups (with vs. without the marker) have the same size ? ##
  
  ### compare fitness with vs. without the marker ##
  fitt_t_test_res <-  t.test(fitt_inds_with_marker,fitt_inds_without_marker,var.equal=TRUE, alternative = "greater") # t.test(x,y,var.equal=TRUE)
  # print(t_test_res) # effect sizes
  print(fitt_t_test_res) # fitness
  
  if(fitt_t_test_res$p.value <= 0.05){
    print("found good marker")
    curr_res <- c(curr_ind_id, max(f1_inds_with_marker$GEBV)) 
    
    # save one f1 descendant to new group      
    new_f1_group <- make.group(f1_inds_with_marker[1,1]) 
    
    # delete.group(f1) # important step otherwise the session terminates
    
    ### bc1. second parent is the F1 ###
    bc1 <- cross.combinations(first.parents = elite_line_id ,second.parents = new_f1_group) # 
    bc1_group_info <- data.frame(Index=see.group.data(bc1,data.type ="X"), GEBV=see.group.data(bc1,data.type = "BV"), G=see.group.data(bc1,data.type ="G"))
    bc1_geno_at_marker <- substr((bc1_group_info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)
    while(bc1_geno_at_marker == "00"){
      delete.group(bc1)
      bc1 <- cross.combinations(first.parents = elite_line_id ,second.parents = new_f1_group)
      bc1_group_info <- data.frame(Index=see.group.data(bc1,data.type ="X"), GEBV=see.group.data(bc1,data.type = "BV"), G=see.group.data(bc1,data.type ="G"))
      bc1_geno_at_marker <- substr((bc1_group_info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)
      print(bc1_geno_at_marker)
    }
    
    
    bc1S1 <- SelfIndividualWithMarker(bc1_group_info, curr_marker_ind_in_map, number_of_offsprings_bc1S1_S2_S3)
    #  delete.group(bc1)
    bc1S2 <- SelfIndividualWithMarker(bc1S1, curr_marker_ind_in_map, number_of_offsprings_bc1S1_S2_S3)
    # #  delete.group(bc1S1)
    bc1S3 <- SelfIndividualWithMarker(bc1S2, curr_marker_ind_in_map, number_of_offsprings_bc1S1_S2_S3)
    # #  delete.group(bc1S2)
    bc1S4 = self.n.times(bc1S3,n = 1, offspring = bc1s4_offsprings_num)
    #  delete.group(bc1S3)
    
    ## phenotype and genotype ##
    print(see.group.data(bc1S4, "D")) ## IDs
    print(see.group.data(bc1S4, "X")) ## Indexes 
    print(see.group.data(bc1S4, "PED") )## Indexes 
    
    
    ### confirm that the marker has significant effect ##
    # seen_genotype <- see.group.data(group = g0, data.type = "G")
    info <- data.frame(Index=see.group.data(bc1S4,"X"), GEBV=see.group.data(bc1S4,"BV"),G=see.group.data(bc1S4,"G") )
    delete.group(bc1S4)
    ###  select individuals with the marker / without the marker ##
    info$G <- substr((info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2 ) ##indices of the selected marker. >>> check here - curr_marker_index ?
    inds_with_marker <- info[info$G %in% c("11","01","10"),]
    ## get the effect sizes of the selected individuals ###  >> currently it takes the 1:10 individuals
    bc1S4_effect_sizes_with_marker <- (inds_with_marker$GEBV)[1:10]
    print("bc1S4_effect_sizes_with_marker:")
    print(bc1S4_effect_sizes_with_marker)
    
    ## compare fitness ## individuals with marker compared to the tester #
    bc1S4_fitt_inds_with_marker <- exp(-1/2*(bc1S4_effect_sizes_with_marker-0.82)^2/.5) # the updated function for the new optimum environment. 0.82 is the new optimum
    print("bc1S4_fitt_inds_with_marker:")
    print(bc1S4_fitt_inds_with_marker)
    avg_bc1S4_fitt_inds_with_marker <- mean(bc1S4_fitt_inds_with_marker)
    
    ## compare effect sizes with the marker vs. the tester's fitness ### >>> which test to use?
    bc1S4_fitt_t_test_res <-  t.test(avg_bc1S4_fitt_inds_with_marker, tester_fitt, var.equal=TRUE, alternative = "greater") # t.test(x,y,var.equal=TRUE)
    print(bc1S4_t_test_res) # effect sizes
    print(bc1S4_fitt_t_test_res) # fitness
    
    if(bc1S4_fitt_t_test_res$p.value <= 0.05) {
      print("found good marker")
      marker_effect_size <- getMarkerEffectSize(curr_marker, effect.file)
      curr_marker_and_effect <- c(curr_marker,marker_effect_size)
      varified_good_markers <- rbind(varified_good_markers,curr_marker_and_effect )
      bc1s4_res <- c(curr_ind_id, max(inds_with_marker$GEBV)) 
      ## to add - save also genotype data
      all_bc1S4_results <- rbind(all_bc1S4_results, bc1s4_res)
    }
    
  }else{
    message <- sprintf("f1 fitness comparison - the marker %s is not significant. Skipping to next marker ", curr_marker)
    print(message)
    clear.simdata()
    next
  }

  clear.simdata() ## needs to clear previous simulation data to run a new one
}
write.csv(all_bc1S4_results, file = out_file , append = TRUE)
write.csv(varified_good_markers, file = good_markers_out_file , append = TRUE)
