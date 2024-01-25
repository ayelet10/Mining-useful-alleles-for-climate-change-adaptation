### the breeding scheme for the prediction only approach ##
## based on the R script 3_prediction_only__genomicSimulation.R ##
start_time <- Sys.time()
print(start_time)
library(genomicSimulation)
library(rrBLUP)

args <- commandArgs(trailingOnly = TRUE) # TRUE
print(args)
seed_num = args[1]
## to run local ##
# seed_num = 1231227
# my_wd <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/"
# my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/results_pred_qtl_updated_map_and_ind_files/",seed_num )
my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/run3_all_causals/",seed_num )

# my_input_files_dir <- sprintf("/home/salmanay/projects/breeding_simulations/%s/",seed_num )
my_input_files_dir <- sprintf("/home/salmanay/projects/breeding_simulations/%s/run3_all_causals/",seed_num )

dir.create(file.path(my_wd)) # if dir exists, it just gives a warning
setwd(file.path(my_wd))
set.seed(seed_num) ## added 21 September 2023

fitness_function <- function(effect_size){
  fitness <- exp(-1/2*(effect_size-0.82)^2/.5) # the updated function for the new optimum environment. 0.82 is the new optimum 
  return(fitness)
}

SelectBestPredicted <- function(group_info, predicted_fitness,percent_to_select){
  ## select best 20% individuals based on the Genomic prediction model and predicted fitness#
  indices_ordered <- order(predicted_fitness, decreasing=T)
  my_group_info_ordered <- group_info[indices_ordered,]
  ## get the id of the predicted
  predicted_group <- my_group_info_ordered$Index
  predicted_best  <- predicted_group[1:(percent_to_select*length(predicted_group))] ## select top 20%
  return(predicted_best)
}

CalAccuracy <- function(group_info, group_name, pred_fitness_group){
  group_fitness      <- fitness_function(group_info$GEBV)
  fitness_accuracy   <- cor(pred_fitness_group, group_fitness) #, use = "complete")
  fitness_accuracy_res  <- c(group_name, fitness_accuracy)
  return(fitness_accuracy_res)
}
getGroupInfoAsDF <- function(group_name){
  df_group_info <- data.frame(Index=see.group.data(group_name,data.type ="X"), GEBV=see.group.data(group_name,data.type = "BV"), G=see.group.data(group_name,data.type ="G"))
  return(df_group_info)
}
#################
## input data  ##
#################
# effect.file <- sprintf("%s_effect_file.txt", seed_num )
# effect.file <- "effect_file_10_snps.txt"
# effect.file <- sprintf("%s%s_effect_file_all_markers.txt", my_wd, seed_num)
# effect.file <- sprintf("%s%s_effect_file_all_markers_as_zero_effect_added_true_effects.txt", my_wd, seed_num)
# effect.file <- sprintf("%s%s_effect_file_sal.txt", my_wd, seed_num)
# effect.file <- sprintf("%s%s_effect_file_with_neutral_sal.txt", my_input_files_dir, seed_num)
effect.file <- sprintf("%s%s_effect_file_with_neutral_and_causal_sal_run3.txt", my_input_files_dir, seed_num)


# map.file <- sprintf("%s_map_file.txt", seed_num)
# map.file <- sprintf("%s%s_map_file_muts_full.txt", my_wd, seed_num)
# map.file <- sprintf("%s%s_map_file_muts_full_with_neutral.txt", my_wd, seed_num)
# map.file    <- sprintf("%s%s_map_file_muts_full_with_neutral_fixed.txt", my_input_files_dir, seed_num)## 15 June 2023. fix of chr names (caused mis-reading of the wrong rows)
map.file    <- sprintf("%s%s_map_file_muts_full_with_neutral_fixed_run3.txt", my_input_files_dir, seed_num)## 15 June 2023. fix of chr names (caused mis-reading of the wrong rows)



# read tester allele file
tester_geno_file <- sprintf("%s%s_tester_allele_file_muts_full_with_neutral_22.txt",my_input_files_dir, seed_num) # Jan-11-24
# tester_geno_file <- sprintf("%s%s_tester_allele_file_muts_full_with_neutral.txt",my_input_files_dir, seed_num)
# tester_geno_file <- sprintf("%s_tester_allele_file_muts_full.txt", seed_num)

#### tests ##
# tester_geno <- read.table(sprintf("%s_tester_allele_file.txt", seed_num))
# tester_geno_file <- sprintf("%s_tester_allele_file.txt", seed_num)
# tester_geno_file <- sprintf("%s_tester_allele_file_10_snps.txt", seed_num) ## for testing the script
#############################
## ids predicted landraces ## 
#############################
## predicted ids from "prediction_with_rrBLUP.R" run
ids_predicted <- read.csv(file = paste0(my_wd,"ids_predicted_landraces_mixed_solve.csv") )
# ids_predicted <- read.csv(file = sprintf("%s%sids_predicted_landraces_mixed_solve.csv", my_wd, seed_num)) ## to add seed number to file name in prediction script
# ids_predicted <- read.csv("C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/ids_predicted_landraces_mixed_solve.csv") 

predicted_lr  <- ids_predicted$x
num_predicted_individuals <- length(predicted_lr)
elite_line_id <- "tester"
f1_num_inds   <-  100  #2 #100 ## each combination tester*predicted lr will give 100 offsprings
num_total_f1  <-   num_predicted_individuals*f1_num_inds   # 5000  ## 300 ##
percent_to_select_bc1_to_bc1s3 <- 0.2  # 20% percent progeny to select in bc1, bc1s1, bc1s2, bc1s3 
percent_to_select_bc1s4 <- 0.125  # 20% percent progeny to select in bc1, bc1s1, bc1s2, bc1s3

##########################
## outputs              ## 
##########################
fitness_file    <- sprintf("%s_fitness_prediction_simulation_bc1s4.txt",seed_num)
# best_bc1s4_fitness_file <- sprintf("%s_best_individual_fitness_prediction_simulation_bc1s4.txt",seed_num) ## added Dec 4 2023
accuracies_file <- sprintf("%s_genomic_selection_model_accuracies.txt", seed_num)
accuracies_all  <- c()
best_actual_fitness_all_steps_file <- sprintf("%s_genomic_selection_best_actual_fitness_all_steps.txt", seed_num)
best_actual_fitness_all_steps <- c()
bc1s4_geno_file <- sprintf("%s_genomic_selection_bc1s4_geno.txt", seed_num)
##########################
## recurrent selection  ## 
##########################
## cross all predicted individuals to elite and do recurrent selection ## the input file contains genotypes of all the predicted individuals ##
# pred_lr <- "tsk_488" #"tsk_9455"
# print(pred_lr)
# predicted_lr_txt_file <- sprintf("%s%s_allele_file_as_txt.txt",my_wd, pred_lr )
# predicted_lr_txt_file <- sprintf("%s%s_allele_file_10_snps.txt",my_wd, pred_lr )

# ReplaceMarkerNames(predicted_lr_txt_file) # marker names updated to match names in map file. example: "1_13_1" will be "13"
# predicted_lr_txt_file <- sprintf("%s%s_updated_markers_names_allele_file_as_txt.txt",my_wd, pred_lr )


predicted_lr_txt_file <- sprintf("%s%s_all_predicted_indiv_allele_file_as_txt.txt", my_wd, seed_num)

# predicted_lr_txt_file <- "/home/salmanay/projects/breeding_simulations/1231099/1231099_all_predicted_indiv_allele_file_as_txt_4_individuals.txt"
# predicted_lr_txt_file <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/my_data_to_qtl2/test_3_subset_10_sites/tsk_488_pred_sites_shared_with_tester_fake_3_inds.txt" ## fake data 3 individuals in file
# predicted_lr_txt_file <- "tsk_630_updated_markers_names_allele_file_as_txt.txt" ##

num_my_predicted <- length(predicted_lr)
g0 <- genomicSimulation::load.data(tester_geno_file, map.file = map.file, effect.file = effect.file) # load tester 
print(see.existing.groups())
# see.group.data(g0, data.type = 'X')
save.genome.model(file  = "g0_genome_model_breeding_scheme_prediction.txt")
genome_model  <- read.csv("g0_genome_model_breeding_scheme_prediction.txt", sep = "\t")
markers_names <- genome_model$name
# tsk_8143 <- "tsk_8143_updated_markers_names_allele_file_as_txt_GEA.txt"
# g0 <- load.more.genotypes(tsk_8143)

## add predicted individuals data
g0 <- load.more.genotypes(predicted_lr_txt_file)
print(see.existing.groups())
warnings()

# get actual fitness to record #
g0_group_info <- getGroupInfoAsDF(g0)
g0_fitness <- fitness_function(g0_group_info$GEBV)
g0_fitness <- c("g0", g0_fitness)
max_g0_fit <- c("g0", max(g0_fitness))
best_actual_fitness_all_steps <- rbind(best_actual_fitness_all_steps, max_g0_fit)

## tests ##
# f1 <- cross.combinations(first.parents = elite_line_id, second.parents = c("tsk_630"), give.names = TRUE, name.prefix = "f1",offspring = 1)
# f1 <- cross.combinations(first.parents = rep(elite_line_id, times = 3), second.parents = c("tsk_488", "tsk_488_a", "tsk_488_b"), give.names = TRUE, name.prefix = "f1",offspring = f1_num_inds)

# first cross: elite X predicted landrace #
f1 <- cross.combinations(first.parents = rep(elite_line_id, times = num_my_predicted), second.parents = predicted_lr, give.names = TRUE, name.prefix = "f1",offspring = f1_num_inds)
f1_indexes <- see.group.data(group = f1, data.type = "X")

see.existing.groups()
print(f1_indexes[1:10])
# get actual fitness to record #
f1_group_info <- getGroupInfoAsDF(f1)
f1_fitness <- fitness_function(f1_group_info$GEBV)
max_f1_fit <- c("f1", max(f1_fitness))
best_actual_fitness_all_steps <- rbind(best_actual_fitness_all_steps, max_f1_fit)

print("after f1 crosses")
bc1 <- cross.combinations(first.parents = rep(elite_line_id, num_total_f1), second.parents = f1_indexes)  ##  back cross the f1 with the elite

see.existing.groups()
bc1_group_info <- getGroupInfoAsDF(bc1)  # a function from the script "breeding_GEA_approach.R"
#######################################################
### Genomic model on bc1 to select best individuals ###
### (bc1 is genotyped and phenotyped)               ###
#######################################################
## select random 80% of bc1 to be the training set ##
train <- sample(1:num_total_f1, 0.8*num_total_f1)
test  <- setdiff(1:num_total_f1, train)

bc1_geno=see.group.data(bc1, data.type = "G")
bc1_geno_mat <- do.call(rbind,lapply(bc1_geno,function(bc1_geno) colSums(matrix(as.numeric(strsplit(bc1_geno,'')[[1]]),nr=2))))
colnames(bc1_geno_mat) <- markers_names
rownames(bc1_geno_mat) <- 1:nrow(bc1_geno_mat)

bc1_group_info$fitness <- fitness_function(bc1_group_info$GEBV)
fitness_test           <- bc1_group_info$fitness[test]
fitness_train          <- bc1_group_info$fitness[train]

## separate to test and train first, later convert the G to a matrix
test_bc1_group_G  <- bc1_geno_mat[rownames(bc1_geno_mat) %in% test, ]
train_bc1_group_G <- bc1_geno_mat[rownames(bc1_geno_mat) %in% train, ]

print("before rrblup")
print(Sys.time())
## run the mixed.solve ##
dim(train_bc1_group_G)
dim(fitness_train)
print(head(fitness_train))
RRBLUP2 <- mixed.solve(y = fitness_train, Z = train_bc1_group_G, method="REML") ###
print("after rrblup")
print(Sys.time())

## test prediction accuracy on the 20% bc1 test set ##
markers_estimations <- as.matrix(RRBLUP2$u)
print(head(markers_estimations))
pred_fitness_test   <- test_bc1_group_G %*% markers_estimations  ## %*% matrix multiplication
print(head(pred_fitness_test))

## predicted bc1 fitness 
pred_fitness_bc1   <- bc1_geno_mat %*% markers_estimations
print("after calculating fitness from marker estimations")
print(head(pred_fitness_bc1))
### save test accuracy results ###
fitness_accuracy_test         <- cor(pred_fitness_test, fitness_test) #, use = "complete")
test_fitness_accuracy_to_save <- c("test", fitness_accuracy_test)
accuracies_all                <- rbind(accuracies_all, test_fitness_accuracy_to_save)
### save bc1 accuracy results ###
fitness_accuracy_bc1         <- cor(pred_fitness_bc1, bc1_group_info$fitness) #, use = "complete")
bc1_fitness_accuracy_to_save <- c("bc1", fitness_accuracy_bc1)
accuracies_all               <- rbind(accuracies_all, bc1_fitness_accuracy_to_save)

## bc1 - select best bc1 individuals based on the GS model and self best individuals ##
indices_ordered  <- order(pred_fitness_bc1, decreasing=T)  #[1:50]
my_bc1_group_info_ordered <- bc1_group_info[indices_ordered,]
## get the id of the predicted
predicted_bc1 <- my_bc1_group_info_ordered$Index  ## get the ids for best genotype
print(c("num predicted_bc1", length(predicted_bc1)))
predicted_best_bc1 <- predicted_bc1[1:(percent_to_select_bc1_to_bc1s3*length(predicted_bc1))] ## select top 20%
print(c("num predicted_best_bc1", length(predicted_best_bc1)))

# get actual fitness to record #
max_bc1_fit <- c("bc1", max(bc1_group_info$fitness))
best_actual_fitness_all_steps <- rbind(best_actual_fitness_all_steps, max_bc1_fit)


best_bc1_group     <- make.group(as.integer(predicted_best_bc1)) ## create new group with best bc1
bc1s1              <- self.n.times(best_bc1_group, n = 1)        ## self to get bc1S1
print("after bc1S1")
see.existing.groups()
bc1s1_group_info   <- getGroupInfoAsDF(bc1s1)
# get actual fitness to record #
bc1s1_fitness <- fitness_function(bc1s1_group_info$GEBV)
max_bc1s1_fit <- c("bc1s1", max(bc1s1_fitness))
best_actual_fitness_all_steps <- rbind(best_actual_fitness_all_steps, max_bc1s1_fit)
## get predicted by the GS model
bc1s1_geno=see.group.data(bc1s1, data.type = "G")
bc1s1_geno_mat <- do.call(rbind,lapply(bc1s1_geno,function(bc1s1_geno) colSums(matrix(as.numeric(strsplit(bc1s1_geno,'')[[1]]),nr=2))))
colnames(bc1s1_geno_mat) <- markers_names
rownames(bc1s1_geno_mat) <- 1:nrow(bc1s1_geno_mat)

print(bc1s1_geno_mat[1:10,1:10])
pred_fitness_bc1s1 <-  bc1s1_geno_mat %*% markers_estimations  ## %*% matrix multiplication
predicted_bc1S1    <- SelectBestPredicted(bc1s1_group_info, pred_fitness_bc1s1,percent_to_select_bc1_to_bc1s3)
## bc1S1 accuracy
bc1S1_fitness_accuracy <- CalAccuracy(bc1s1_group_info, "bc1S1", pred_fitness_bc1s1)
accuracies_all         <- rbind(accuracies_all, bc1S1_fitness_accuracy)
print("after bc1 accuracy")
print(accuracies_all)

## bc1S1 - select best 20% bc1S1 based on the Genomic prediction model #
best_bc1s1_group   <- make.group(as.integer(predicted_bc1S1)) ## create new group with best bc1
bc1S2              <- self.n.times(best_bc1s1_group, n = 1) ## self to get bc1S2
bc1s2_group_info   <- getGroupInfoAsDF(bc1S2)  # a function from the script "GEA_breeding_scheme.R"
# get actual fitness to record #
bc1s2_fitness <- fitness_function(bc1s2_group_info$GEBV)
max_bc1s2_fit <- c("bc1s2", max(bc1s2_fitness))
best_actual_fitness_all_steps <- rbind(best_actual_fitness_all_steps, max_bc1s2_fit)

bc1s2_geno=see.group.data(bc1S2, data.type = "G")
bc1s2_geno_mat <- do.call(rbind,lapply(bc1s2_geno,function(bc1s2_geno) colSums(matrix(as.numeric(strsplit(bc1s2_geno,'')[[1]]),nr=2))))
colnames(bc1s2_geno_mat) <- markers_names
rownames(bc1s2_geno_mat) <- 1:nrow(bc1s2_geno_mat)

## get predicted by the GS model
pred_fitness_bc1s2 <- bc1s2_geno_mat %*% markers_estimations  ## %*% matrix multiplication
predicted_bc1s2    <- SelectBestPredicted(bc1s2_group_info, pred_fitness_bc1s2,percent_to_select_bc1_to_bc1s3)
## bc1S2 accuracy
bc1S2_fitness_accuracy <- CalAccuracy(bc1s2_group_info, "bc1S2", pred_fitness_bc1s2)
accuracies_all         <- rbind(accuracies_all, bc1S2_fitness_accuracy)
print("after bc1s1 accuracy")
print(accuracies_all)

## bc1S2 - select best 20% bc1S2 based on the Genomic prediction model #
best_bc1s2_group     <- make.group(as.integer(predicted_bc1s2)) ## create new group with best bc1
bc1S3                <- self.n.times(best_bc1s2_group, n = 1) ## self to get bc1S3
bc1S3_group_info     <- getGroupInfoAsDF(bc1S3)  # a function from the script "GEA_breeding_scheme.R"
# get actual fitness to record #
bc1s3_fitness <- fitness_function(bc1S3_group_info$GEBV)
max_bc1s3_fit <- c("bc1S3", max(bc1s3_fitness))
best_actual_fitness_all_steps <- rbind(best_actual_fitness_all_steps, max_bc1s3_fit)
## get predicted by the GS model
bc1S3_geno=see.group.data(bc1S3, data.type = "G")
bc1S3_geno_mat <- do.call(rbind,lapply(bc1S3_geno,function(bc1S3_geno) colSums(matrix(as.numeric(strsplit(bc1S3_geno,'')[[1]]),nr=2))))
colnames(bc1S3_geno_mat) <- markers_names
rownames(bc1S3_geno_mat) <- 1:nrow(bc1S3_geno_mat)

pred_fitness_bc1s3 <- bc1S3_geno_mat %*% markers_estimations  ## %*% matrix multiplication
predicted_bc1S3    <- SelectBestPredicted(bc1S3_group_info, pred_fitness_bc1s3,percent_to_select_bc1_to_bc1s3)
## bc1S3 accuracy ##
bc1S3_fitness_accuracy <- CalAccuracy(bc1S3_group_info, "bc1S3", pred_fitness_bc1s3)
accuracies_all         <- rbind(accuracies_all, bc1S3_fitness_accuracy)
print("after bcs31 accuracy")
print(accuracies_all)

## bc1S3 - select best 20% bc1S3 based on the Genomic prediction model #
best_bc1s3_group     <- make.group(as.integer(predicted_bc1S3)) ## create new group with best bc1
bc1S4                <- self.n.times(best_bc1s3_group, n = 1) ## self to get bc1S4
bc1S4_group_info     <- getGroupInfoAsDF(bc1S4)  # a function from the script "GEA_breeding_scheme.R"
## get predicted by the GS model
bc1S4_geno=see.group.data(bc1S4, data.type = "G")
bc1S4_geno_mat <- do.call(rbind,lapply(bc1S4_geno,function(bc1S4_geno) colSums(matrix(as.numeric(strsplit(bc1S4_geno,'')[[1]]),nr=2))))
colnames(bc1S4_geno_mat) <- markers_names
rownames(bc1S4_geno_mat) <- 1:nrow(bc1S4_geno_mat)
pred_fitness_bc1S4 <-  bc1S4_geno_mat %*% markers_estimations  ## %*% matrix multiplication
print(pred_fitness_bc1S4)
### bc1s4 accuracy
bc1S4_fitness_accuracy <- CalAccuracy(bc1S4_group_info, "bc1S4", pred_fitness_bc1S4)
accuracies_all         <- rbind(accuracies_all, bc1S4_fitness_accuracy)
print("after bc1s4 accuracy")
print(accuracies_all)
# save bc1s4 fitness to file # 
bc1S4_fitness <- fitness_function(bc1S4_group_info$GEBV)
write.table(bc1S4_fitness,file = fitness_file)
write.table(accuracies_all, file = accuracies_file)
### save BC1S4 individual that has highest predicted fitness ## added Dec 4 2023 ##
# predicted_best_bc1S4_individual    <- SelectBestPredicted(bc1S4_group_info, pred_fitness_bc1S4,percent_to_select_bc1s4)
# best_bc1s4_group     <- make.group(as.integer(predicted_best_bc1S4_individual)) ## create new group with best bc1
# best_bc1s4_group_info <- getGroupInfoAsDF(best_bc1s4_group)
# best_bc1S4_fitness <- fitness_function(best_bc1s4_group_info$GEBV)
# write.table(best_bc1S4_fitness, file = best_bc1s4_fitness_file)
# get actual fitness to record #
bc1s4_fitness <- fitness_function(bc1S4_group_info$GEBV)
max_bc1s4_fit <- c("bc1S4", max(bc1s4_fitness))
best_actual_fitness_all_steps <- rbind(best_actual_fitness_all_steps, max_bc1s4_fit)
write.table(best_actual_fitness_all_steps, file = best_actual_fitness_all_steps_file)
# save bc1s4 genotypes# 
save.genotypes(filename = bc1s4_geno_file, group = bc1S4)
print(see.existing.groups())
clear.simdata()

end_time <- Sys.time()
print(c("end time: ", end_time))
####################