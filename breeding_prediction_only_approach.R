### the breeding scheme for the prediction only approach ##
## based on the R script 3_prediction_only__genomicSimulation.R ##
start_time <- Sys.time()
print(start_time)
library(genomicSimulation)

args <- commandArgs(trailingOnly = TRUE) # TRUE
print(args)
seed_num = args[1]
# seed_num = 1231227 #args[1]

my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/",seed_num )
# my_wd <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/"
dir.create(file.path(my_wd)) # if dir exists, it just gives a warning
setwd(file.path(my_wd))
gebv_out_f_extention <-  "_prediction_simulation_bc1s4_gebv_3.txt"
percentage_inds =  100 # 20  percentage of decadents to take in bc1s3 and bc1s4
bc1_n <- 2 # 250
bc1s2_n <-  2 # 250
# bc1s3_n <- 1 # 50 
# allele_file <- "tsk_9455_allele_file_as_txt.txt"
## write allele file with marker names as in map file
ReplaceMarkerNames <- function (allele_file){
   my_alleles <- read.table(allele_file, sep = "\t" , header = TRUE ,colClasses = c("character","character"))
   # my_alleles <- read.table(predicted_lr_txt_file, sep = "\t", header = TRUE )
   my_alleles$name <- sub("*._", "",my_alleles$name )
   my_alleles$name <- sub("_.*" , "" ,my_alleles$name )
   write.table(my_alleles, sprintf("%s%s_updated_markers_names_allele_file_as_txt.txt", my_wd, pred_lr), sep = "\t", quote = FALSE, row.names  = FALSE)
}

## function - read gebv file and calculate fitness ##
CalFitnessFromGEBV <- function(gebv_file_name_extension, predicted_line){
  gebv_file_name <- paste0(predicted_line, gebv_file_name_extension )
  print(gebv_file_name)
  gebv_out_f <- read.table(gebv_file_name)
  y <- gebv_out_f$V2
  gebv_out_f$fitt <- exp(-1/2*(y-0.82)^2/.5) # the updated function for the new optimum environment. 0.82 is the new optimum 
  return(gebv_out_f)
}
#################
## input data  ##
#################
# effect.file <- sprintf("%s_effect_file.txt", seed_num )
# effect.file <- "effect_file_10_snps.txt"
# effect.file <- sprintf("%s%s_effect_file_all_markers.txt", my_wd, seed_num)
# effect.file <- sprintf("%s%s_effect_file_all_markers_as_zero_effect_added_true_effects.txt", my_wd, seed_num)
effect.file <- sprintf("%s%s_effect_file_sal.txt", my_wd, seed_num)

# map.file <- sprintf("%s_map_file.txt", seed_num)
map.file <- sprintf("%s%s_map_file_muts_full.txt", my_wd, seed_num)

# read tester allele file
tester_geno_file <- sprintf("%s_tester_allele_file_muts_full.txt", seed_num)
# tester_geno <- read.table(sprintf("%s_tester_allele_file.txt", seed_num))
# tester_geno_file <- sprintf("%s_tester_allele_file.txt", seed_num)
# tester_geno_file <- sprintf("%s_tester_allele_file_10_snps.txt", seed_num)

elite_line_id <- "tester"
#############################
## ids predicted landraces ## 
#############################
## predicted ids from "prediction_with_rrBLUP.R" run
ids_predicted <- read.csv(file = paste0(my_wd,"ids_predicted_landraces_mixed_solve.csv") )
predicted_lr <- ids_predicted$x

##########################
## recurrent selection  ## 
##########################
## for each selected landrace, cross to elite and do recurrent selection ##
for (pred_lr in predicted_lr[1]){
  # pred_lr <- "tsk_9455"
  print(pred_lr)
  # predicted_lr_txt_file <- sprintf("%s%s_allele_file_as_txt.txt",my_wd, pred_lr )
  # predicted_lr_txt_file <- sprintf("%s%s_allele_file_10_snps.txt",my_wd, pred_lr )
  
  # ReplaceMarkerNames(predicted_lr_txt_file) # marker names updated to match names in map file. example: "1_13_1" will be "13"
  predicted_lr_txt_file <- sprintf("%s%s_updated_markers_names_allele_file_as_txt.txt",my_wd, pred_lr )

  g0 <- genomicSimulation::load.data(tester_geno_file, map.file = map.file, effect.file = effect.file) # load tester 

  ## add predicted individual data
  g0 <- load.more.genotypes(predicted_lr_txt_file)
  # load.different.effects(effect.file)

  # first cross will be elite X predicted landrace - define the parents by their id
  # f1 <- cross.combinations( first.parents = elite_line_id ,second.parents = c("30001729"))
  f1 <- cross.combinations( first.parents = elite_line_id ,second.parents = pred_lr, give.names = TRUE, name.prefix = "f1")
  
  
  # see.group.data(group = g0, data.type = "D")
  # see.group.data(group = f1, data.type = "G")
  print(see.group.data(group = f1, data.type = "PED"))
  # print(see.group.data(group = f1, data.type = "GEBV"))
  
  # Delete groups we are not currently using, to free up some memory.
  # delete.group(g0) ## keep it for the 'elite' line

  bc1 <- cross.combinations(first.parents = elite_line_id, second.parents = "f13", offspring = bc1_n)  ##  back cross the f1 with the elite
  print(see.group.data(group = bc1, data.type = "PED"))
  # see.group.data(group = bc1, data.type = "G")
  
  ## bc1 <- cross.combinations(first.parents = elite_line_id, second.parents = c("30001729"))   # this works
#  delete.group(f1)
  bc1S1 = self.n.times(bc1, n = 1)
  # see.group.data(group = bc1S1, data.type = "G")
  print("bc1S1 PED:")
  print(see.group.data(bc1S1, "PED"))
#  delete.group(bc1)
  bc1S2 = self.n.times(bc1S1,n = bc1s2_n)

#  delete.group(bc1S1)
  # print("bc1S2 GEBV:")
  # print(see.group.data(bc1S2, "GEBV"))
  see.group.data(group = bc1S2, data.type = "X")
  
  # load.different.effects(effect.file)
  ##  I - select best individuals based on gebv and self best individuals ##
  bc1S2_selected = select.by.gebv(bc1S2, low.score.best = FALSE, percentage = percentage_inds) ## once individuals are selected from a group - it seems that the group is empty
  print("bc1S2_selected PED:")
  print(see.group.data(bc1S2_selected, data.type ="PED"))
  print(see.group.data(group = bc1S2_selected, data.type = "D"))  ## show ID's
  print(see.group.data(group = bc1S2_selected, data.type = "X"))  ## show indexes
  bc1S3 = self.n.times(bc1S2_selected,n = 1)

  ## II - select best individuals based on gebv and self best individuals ##
  bc1S3_selected = select.by.gebv(bc1S3, low.score.best = FALSE, percentage = percentage_inds)
  # print("bc1S3_selected GEBV:")
  # print(see.group.data(bc1S3_selected, "GEBV"))
  
  bc1S4 = self.n.times(bc1S3_selected,n = 1)
  # print("bc1S4 GEBV:")
  # print(see.group.data(bc1S4, "GEBV"))
  
  # save GEBVs to file
  save.GEBVs(filename = sprintf("f1_%s%s", pred_lr, gebv_out_f_extention ), group = f1) # _prediction_simulation_bc1s4_gebv_3.txt
  
  save.GEBVs(filename = sprintf("g0_%s%s", pred_lr, gebv_out_f_extention ), group = g0) # _prediction_simulation_bc1s4_gebv_3.txt
  
  save.GEBVs(filename = sprintf("bc1S3_selected_%s%s", pred_lr, gebv_out_f_extention ), group = bc1S3_selected) # _prediction_simulation_bc1s4_gebv_3.txt
  
  ## save bc1S4 GEBV file to further use #
  save.GEBVs(filename = sprintf("%s%s", pred_lr, gebv_out_f_extention ), group = bc1S4) # _prediction_simulation_bc1s4_gebv_3.txt
  
  gebv_fitt <- CalFitnessFromGEBV(gebv_out_f_extention, pred_lr)
  write.table(gebv_fitt, file = sprintf("fitt_%s%s",pred_lr, gebv_out_f_extention ))
  save.pedigrees("a.txt", bc1S4, type="P")

  print(see.existing.groups())
  clear.simdata() ## needs to clear previous simulation data to run a new one
  
}
end_time <- Sys.time()
print(c("time (min): ", (end_time-start_time)/60))

# my_map <- read.csv(map.file, sep = "\t")
# tester_alleles <- read.csv(tester_geno_file, sep = "\t")
# nrow(tester_alleles[tester_alleles$name %in% my_map$marker,] )
# nrow(my_map[my_map$marker  %in%  tester_alleles$name,])
# eff_file$V1 %in% my_map$marker # FALSE 
