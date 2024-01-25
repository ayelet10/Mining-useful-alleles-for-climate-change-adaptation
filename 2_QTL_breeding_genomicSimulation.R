#########2_QTL_breeding_genomicSimulation.R  ######
start_time <- Sys.time()
print(start_time)
library(genomicSimulation)
library(qtl2)
library(stringr) # to use with the function getPercentLandraceFromInfo
library(dplyr)
## https://kbroman.org/qtl2/assets/vignettes/user_guide.html


args <- commandArgs(trailingOnly = TRUE) # TRUE
print(args)
seed_num = args[1]
dir_outputs = args[2] # "_tester_as_2" # "" # "_random"  ## "_individuals_selected_by_sal_opt" # "_random_corrected_f2_index"## change per run 
ids_file_name = args[3]                # "ids_predicted_landraces_mixed_solve.csv" # "ids_randomly_selected_lr.csv" # "ids_selected_lr_by_sal_value.csv"
do_qtl_mapping <- args[4] #FALSE ## TRUE ## option to *not* run the qtl mapping and instead read the qtl peaks from file 
# 
# dir_outputs = "_tester_as_2" # "" # "_random"  ## "_individuals_selected_by_sal_opt" # "_random_corrected_f2_index"## change per run
# ids_file_name =  "ids_predicted_landraces_mixed_solve.csv" # "ids_randomly_selected_lr.csv" # "ids_selected_lr_by_sal_value.csv"
# do_qtl_mapping <- FALSE
# seed_num = 1231119 ##1231227 #1231409 #
set.seed(seed_num) ## added 21 September 2023

# my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/", seed_num)
my_wd <- sprintf("/home/salmanay/projects/breeding_simulations/%s/run3_all_causals/", seed_num)

# my_wd <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/"
# my_wd <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/my_data_to_qtl2/test_2_subset_100_sites/"
# my_wd <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/my_data_to_qtl2/test_3_subset_10_sites/"
# my_wd <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/my_data_to_qtl2/test_full_genome/"

# my_wd <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/test_specific_seed/1231409/"
# out_dir <- sprintf("/home/salmanay/projects/breeding_simulations/%s/test", seed_num)

# out_dir <- sprintf("/home/salmanay/projects/breeding_simulations/%s/results_pred_qtl_updated_map_and_ind_files/",seed_num )

out_dir <- sprintf("/home/salmanay/projects/breeding_simulations/%s/run3_all_causals/qtl%s/",seed_num, dir_outputs ) ## run3
# out_dir <- sprintf("/home/salmanay/projects/breeding_simulations/%s/run3_all_causals/qtl_random_individuals/",seed_num ) ## run3 - test with random individuals
# out_dir <- sprintf("/home/salmanay/projects/breeding_simulations/%s/run3_all_causals/qtl_individuals_selected_by_sal_opt/",seed_num ) ## run3 - test with individuals selected by env (sal_opt)



## file with list of selected individuals 
ids_predicted <- read.csv(file = paste0(my_wd,ids_file_name)) ## selected by the prediction script

# ids_predicted <- read.csv(file = paste0(my_wd,"ids_predicted_landraces_mixed_solve.csv")) ## selected by the prediction script
# ids_predicted <- read.csv( file = paste0(my_wd,"ids_randomly_selected_lr.csv"))           ## selected randomly                (select_individuals_and_save_to_file.R)
# ids_predicted <- read.csv(file = paste0(my_wd,"ids_selected_lr_by_sal_value.csv"))        ## selected by env value of sal_opt (select_individuals_and_save_to_file.R)

# out_dir <- sprintf("/home/salmanay/projects/breeding_simulations/%s/qtl_res_map_file_fixed/", seed_num) ## also: takes markers with positive effect size
dir.create(out_dir) ## all outputs will be saved here

# my_yaml <- "C:/Users/ayelet/Documents/UCDavis/scripts/simulations/genomicSimulation_tests/my_data_to_qtl2/test_1/test1.yaml"
# setwd(my_wd)
setwd(out_dir) ## the new qtl files will be saved here

# pred_tsk_488 <- read.table("tsk_488_updated_markers_names_allele_file_as_txt.txt", header = TRUE, sep = "\t", colClasses=c('character','character')) 
# head(pred_tsk_488$name)
# head(tester_geno)
# pred_sites_shared_with_tester <-  pred_tsk_488[pred_tsk_488$name %in% tester_geno$name,] # 324 rows not in tester file # 31332 are in tester file
# write.table(pred_sites_shared_with_tester, file= "tsk_488_pred_sites_shared_with_tester.txt", row.names = FALSE, quote = FALSE,sep = "\t")
source("/home/salmanay/projects/breeding_simulations/genomicSimulation/scripts/functions_script_breeding.R")
#################
## input data  ##
#################
# effect.file <- sprintf("%s%s_effect_file_with_neutral_sal.txt", my_wd, seed_num)
effect.file <- sprintf("%s%s_effect_file_with_neutral_and_causal_sal_run3.txt", my_wd, seed_num)

# map.file  <- sprintf("%s%s_map_file_muts_full_with_neutral.txt", my_wd, seed_num)
# map.file    <- sprintf("%s%s_map_file_muts_full_with_neutral_fixed.txt", my_wd, seed_num)## 15 June 2023. fix of chr names (caused mis-reading of the wrong rows)
map.file    <- sprintf("%s%s_map_file_muts_full_with_neutral_fixed_run3.txt", my_wd, seed_num)## 15 June 2023. fix of chr names (caused mis-reading of the wrong rows)

# read tester allele file
if(do_qtl_mapping==TRUE){
  print("reading tester 00")
  tester_geno_file <- sprintf("%s%s_tester_allele_file_muts_full_with_neutral.txt", my_wd, seed_num) ## tester as "00"
}else{   ## option to run only the breeding without the qtl mapping
  print("reading tester 22")
  tester_geno_file <- sprintf("%s%s_tester_allele_file_muts_full_with_neutral_22.txt", my_wd, seed_num) # Jan-11-24. tester as "22" 
}
# tester_geno_file <- sprintf("%s_tester_allele_file_muts_full.txt", seed_num)
# effect.file <- sprintf("%s%s_effect_file_sal.txt", my_wd, seed_num)
# map.file <- sprintf("%s_map_file_muts_full.txt", seed_num) # gmap
# map.file.mbp <- sprintf("%s_map_file_mbp.txt", seed_num) # pmap
# tester_geno_file <- sprintf("%s_tester_allele_file_muts_full.txt", seed_num)
elite_line_id        <- "tester"
num_f2_offsprings    <- 500 # 
bc1s4_offsprings_num <- 100 # 
number_of_offsprings_bc1S1_S2_S3 <- 20 #

################
##  out file  ##
################
# qtl_out_file <- sprintf("%s_all_qtls_bc1S4_homozygous_to_qtls.txt", seed_num) ## all qtls == one qtl per predicted individual #
qtl_out_file <- sprintf("%s%s%s_all_qtls_bc1S4_homozygous_to_qtls.txt", my_wd, seed_num, dir_outputs) ## all qtls == one qtl per predicted individual #
# no_f2_homozygotes_file <- sprintf("%s%s%s_no_f2_homozygotes_file.txt", my_wd, seed_num, dir_outputs)
#############################
## ids predicted landraces ## 
#############################
## chose 10 predicted individuals from the prediction step: prediction_with_rrBLUP.R script ##
predicted_lr <- ids_predicted$x[1:10]
all_bc1S4_results <- c()
# no_f2_homozygotes_lst <- c()
## per individual ##
## for each selected landrace, cross to elite and do recurrent selection ##
# pred_lr <- predicted_lr[10]
for (pred_lr in predicted_lr){
  all_gener_percent_landrace <- data.frame()
  homozygous_bc1S4_group_info <- c()
  peaks_found <- c()
  print(pred_lr)
  
  #################
  ## outputs     ##
  #################
  ## outputs per predicted individual ##
  qtl_control_file  <- sprintf("%s_%s_qtl_control_file.yaml", seed_num, pred_lr)
  zipped_qtl_control_file <- sprintf("%s_%s_qtl_control_file.zip", seed_num, pred_lr)
  f2_geno_file      <- sprintf("%s_%s_f2_genotypes.txt", seed_num, pred_lr)
  f2_pheno_file     <- sprintf("%s_%s_f2_GEBVs.txt", seed_num, pred_lr)
  f2_geno_file_csv  <- sprintf("%s_%s_f2_genotypes.csv", seed_num, pred_lr)
  f2_pheno_file_csv <- sprintf("%s_%s_f2_fitness.csv", seed_num, pred_lr)
  map_file_csv      <- sprintf("%s_%s_map_file_muts_full.csv",seed_num,pred_lr)  ## pred_lr added Jan-16-24
  pmap_file_csv     <- sprintf("%s_%s_map_file_mbp.csv",seed_num,pred_lr)         ## pred_lr added Jan-16-24
  
  # map_file_csv      <- sprintf("%s_map_file_muts_full.csv",seed_num)
  # pmap_file_csv     <- sprintf("%s_map_file_mbp.csv",seed_num)
  
  
  peaks_found_file  <- sprintf("%s_%s_peaks_found.csv", seed_num, pred_lr)
  peaks_effect_sizes <- sprintf("%s_%s_peaks_effect_sizes.csv", seed_num, pred_lr)
  f2_hetero_to_marker_geno_file <- sprintf("%s_%s_f2_hetero_to_marker_genotype.csv", seed_num, pred_lr)
  bc1_fitness_file_csv <- sprintf("%s_%s_bc1_hetero_to_marker_fitness.csv", seed_num, pred_lr)
  # bc1_hetero_to_marker_geno_file <- sprintf("%s_%s_bc1_hetero_to_marker_genotype.csv", seed_num, pred_lr)
  bc1_geno_file <- sprintf("%s_%s_bc1_genotype.csv", seed_num, pred_lr)
  # f2_geno_file_22      <- sprintf("%s_%s_f2_genotypes_22.txt", seed_num, pred_lr) ## genotypes with 2 allele from the tester
  f2_without_marker_geno_file <- sprintf("%s_%s_f2_wo_marker_genotypes.txt", seed_num, pred_lr) ## 
  f2_homo_to_marker_geno_file <- sprintf("%s_%s_f2_homo_to_marker_genotypes.txt", seed_num, pred_lr) ## 
  all_f2_hetero_to_marker_geno_file <- sprintf("%s_%s_all_f2_hetero_to_marker_genotype.csv", seed_num, pred_lr)
  
  f1_geno_file <- sprintf("%s_%s_f1_genotypes.txt", seed_num, pred_lr)
  bc1s1_homo_to_marker_geno_file <- sprintf("%s_%s_bc1s1_homo_to_marker_genotypes.txt", seed_num, pred_lr)
  bc1s1_hetero_to_marker_geno_file <- sprintf("%s_%s_bc1s1_hetero_to_marker_genotypes.txt", seed_num, pred_lr)
  bc1s1_without_marker_geno_file <- sprintf("%s_%s_bc1s1_wo_marker_genotypes.txt", seed_num, pred_lr)
  
  bc1s2_homo_to_marker_geno_file <- sprintf("%s_%s_bc1s2_homo_to_marker_genotypes.txt", seed_num, pred_lr)
  bc1s2_hetero_to_marker_geno_file <- sprintf("%s_%s_bc1s2_hetero_to_marker_genotypes.txt", seed_num, pred_lr)
  bc1s2_without_marker_geno_file <- sprintf("%s_%s_bc1s2_wo_marker_genotypes.txt", seed_num, pred_lr)
  
  bc1s3_homo_to_marker_geno_file <- sprintf("%s_%s_bc1s3_homo_to_marker_genotypes.txt", seed_num, pred_lr)
  bc1s3_hetero_to_marker_geno_file <- sprintf("%s_%s_bc1s3_hetero_to_marker_genotypes.txt", seed_num, pred_lr)
  bc1s3_without_marker_geno_file <- sprintf("%s_%s_bc1s3_wo_marker_genotypes.txt", seed_num, pred_lr)
  
  bc1s4_homo_to_marker_geno_file <- sprintf("%s_%s_bc1s4_homo_to_marker_genotypes.txt", seed_num, pred_lr)
  bc1s4_hetero_to_marker_geno_file <- sprintf("%s_%s_bc1s4_hetero_to_marker_genotypes.txt", seed_num, pred_lr)
  bc1s4_without_marker_geno_file <- sprintf("%s_%s_bc1s4_wo_marker_genotypes.txt", seed_num, pred_lr)
  
  
  
  # qtl_control_file  <- sprintf("%s%s_%s_qtl_control_file.yaml",out_dir, seed_num, pred_lr)
  # zipped_qtl_control_file <- sprintf("%s%s_%s_qtl_control_file.zip",out_dir, seed_num, pred_lr)
  # f2_geno_file      <- sprintf("%s%s_%s_f2_genotypes.txt",out_dir, seed_num, pred_lr)
  # f2_pheno_file     <- sprintf("%s%s_%s_f2_GEBVs.txt",out_dir, seed_num, pred_lr)
  # f2_geno_file_csv  <- sprintf("%s%s_%s_f2_genotypes.csv",out_dir, seed_num, pred_lr)
  # f2_pheno_file_csv <- sprintf("%s%s_%s_f2_fitness.csv",out_dir, seed_num, pred_lr)
  # map_file_csv      <- sprintf("%s%s_map_file_muts_full.csv",out_dir,seed_num)
  # pmap_file_csv     <- sprintf("%s%s_map_file_mbp.csv",out_dir,seed_num)
  # peaks_found_file  <- sprintf("%s%s_%s_peaks_found.csv",out_dir, seed_num, pred_lr)
  ##############
  ## breeding ##
  ##############
  # predicted_lr_txt_file <- sprintf("%s%s_updated_markers_names_allele_file_as_txt.txt",my_wd, pred_lr )
  predicted_lr_txt_file <- sprintf("%s%s_allele_file_as_txt.txt", my_wd, pred_lr)
  
  
  # predicted_lr_txt_file <- "tsk_488_pred_sites_shared_with_tester.txt"
  
  ## load tester 
  g0 <- genomicSimulation::load.data(tester_geno_file, map.file = map.file, effect.file = effect.file) 
  print(see.existing.groups())
  see.group.data(g0, data.type = "N")
  ## save tester's effect size ##
  # tester_effect_size <- see.group.data(g0, data.type = "BV")
  ## add predicted individual data
  g0 <- genomicSimulation::load.more.genotypes(predicted_lr_txt_file)
  print(see.existing.groups())
  see.group.data(g0, data.type = "N")
  g0_info <- getGroupInfoAsDF(g0) ## added Dec-7-2023
  fitness_g0    <- fitness_function(g0_info$GEBV) ## added Dec-7-2023
  write.table(fitness_g0, file = sprintf("%s_%s_g0_fitness.txt",seed_num, pred_lr)) ## added Dec-7-2023
  
  # first cross will be elite X predicted landrace - define the parents by their id
  f1 <- cross.combinations(first.parents = elite_line_id ,second.parents = pred_lr, give.names = TRUE, name.prefix = "f1")
  f1_info <- getGroupInfoAsDF(f1)
  # SaveGEBVandFitness(f1, "f1", seed_num, pred_lr) ## added Dec-7-2023 # does not work
  fitness_f1    <- fitness_function(f1_info$GEBV)
  write.table(fitness_f1, file = sprintf("%s_%s_f1_fitness.txt", seed_num, pred_lr))
  save.genotypes(filename = f1_geno_file, group = f1)
  
  ## create f2 population ##
  f2 <- self.n.times(f1, n=1, offspring = num_f2_offsprings)
  f2_info <- data.frame(Index=see.group.data(f2,data.type ="X"), GEBV=see.group.data(f2,data.type = "BV"), G=see.group.data(f2,data.type ="G"))
  
  if(do_qtl_mapping==TRUE){ ## option to run only the breeding without the qtl mapping
    print("doing qtl mapping")
    save.genotypes(filename = f2_geno_file, group = f2) ## after deleting g0 and f1, f2 is the only group that will be printed. ERROR: the tester and f1 are printed --- needs to be deleted
    # save.genotypes(filename = paste0("groupf2", f2_geno_file),group = f2) ## after deleting g0 and f1, f2 is the only group that will be printed. ERROR: the tester and f1 are printed --- needs to be deleted
    save.GEBVs(filename = f2_pheno_file, group = f2)
    
    f2_geno <- read.table(f2_geno_file, sep = "\t", header = TRUE, check.names = FALSE, colClasses = 'character')
    # print(head(f2_geno))
    f2_geno_filtered <- f2_geno[vapply(f2_geno, function(x) length(unique(x)) > 1, logical(1L))] ## filter out columns with the same genotype in all individuals
    # print(head(f2_geno_filtered))
    print("genotype length (ncol f2)")
    print(ncol(f2_geno_filtered))
    write.csv(f2_geno_filtered, file = f2_geno_file_csv, row.names = FALSE) 
    # write.csv(f2_geno, file = f2_geno_file_csv, row.names = FALSE) 
    
    f2_pheno <- read.table(f2_pheno_file, header = FALSE, check.names = FALSE)
    head(f2_pheno[1:2, 1:2])
    colnames(f2_pheno) <- c("id","gebv")
    f2_pheno$fitness   <- fitness_function(f2_pheno$gebv) ## calculate fitness from GEBV ( =effect sizes)
    write.table(f2_pheno[,c("id","fitness")], file = f2_pheno_file_csv, sep = "," , col.names = TRUE, row.names = FALSE)
    
    map_as_txt <- read.table(map.file, sep = "\t", header = TRUE, check.names = FALSE)
    # print(head(map_as_txt))
    ## filter map to match markers in f2_geno_filtered (after f2_geno was filtered)
    filtered_map_as_txt <- map_as_txt[map_as_txt$marker %in% colnames(f2_geno_filtered)[2:ncol(f2_geno_filtered)],]
    write.csv(filtered_map_as_txt ,file = map_file_csv, row.names = FALSE )
    # write.csv(map_as_txt ,file = map_file_csv, row.names = FALSE )
    ## physical map (Mbp)  >> make sure that it is not needed
    # p_map_as_txt <- read.table(map.file.mbp, sep = "\t", header = TRUE, check.names = FALSE)
    # head(p_map_as_txt)
    # write.csv(p_map_as_txt ,file = pmap_file_csv , row.names = FALSE)
    #####################################
    ## prepare files for qtl analysis ##
    ####################################
    geno_codes <-  c("1", "2", "2", "3") # c("1", "1", "1" , "1", "2", "2","2", "2", "3")
    names(geno_codes) <-   c("00", "10", "01", "11") # c("00", "22", "20", "02", "10", "01", "12", "21" ,"11") 
    write_control_file(qtl_control_file, crosstype = "f2", geno_file = f2_geno_file_csv, pheno_file = f2_pheno_file_csv, # pmap_file =pmap_file_csv,
                       gmap_file = map_file_csv, geno_transposed = FALSE, pheno_transposed = FALSE, overwrite = TRUE, alleles = c("0", "1"), geno_codes = geno_codes )

    ## run qtl analysis ##
    my_data_files <- read_cross2(qtl_control_file) # to read zipped files use:  zipped_qtl_control_file
    map <- insert_pseudomarkers(my_data_files$gmap, step=0)  # insert pseudomarkers into the genetic map >> remove. not needed for our data. step = 0 is no pseudomarkers
    pr  <- calc_genoprob(my_data_files, map, error_prob=0.002) # , cores=4)# calculate the QTL genotype probabilities   ## error_prob = ? should be zero, no genotyping errors in our simulated data ##
    out <- scan1(pr, my_data_files$pheno) #, Xcovar=Xcovar) # , cores=4)                                                ## Xcovar needed ?? ##
    pdf(file = sprintf("%s_qtl_results_2.pdf", pred_lr))
    par(mar=c(5.1, 4.1, 1.1, 1.1))
    ymx <- maxlod(out) # overall maximum LOD score
    plot(out, map, lodcolumn=1, col="slateblue", ylim=c(0, ymx*1.02))
    legend("topleft", lwd=2, col=c("slateblue", "violetred"), colnames(out), bg="gray90")
    dev.off()
    ## did we find advantageous qtls? ##
    # identify a set of LOD peaks that exceed some threshold
    peaks_found <- find_peaks(out, map, threshold=4, peakdrop=1.8, drop=1.5) #  ##  ?? ## peakdrop indicates the amount that the LOD curve must drop below the lowest of two adjacent peaks, threshold = 4 is restrictive
 #   peaks_found <- find_peaks(out, map, threshold=0, peakdrop=Inf) #  ##  Jan-24
    
    write.csv(peaks_found, file = peaks_found_file)
    # peaks_found <- read.table(peaks_found_file, sep= ",", header = T) ## for testing
    ## in case no marker is found - move to next individual ##
    if (nrow(peaks_found) == 0){
      print("no peaks were found in qtl mapping")
      clear.simdata()
      ## delete the two large file with genotypes ## 
      if (file.exists(f2_geno_file_csv)) {
        #Delete file if it exists
        file.remove(f2_geno_file_csv)
      }
      if (file.exists(f2_geno_file)) {
        #Delete file if it exists
        file.remove(f2_geno_file)
      }
      next
    }
    
  
    map_file <- read.table(map.file, header= T)
    effect_file <- read.table(effect.file) ## marker, allele, effect size
    colnames(effect_file) <- c("marker", "allele", "effect_size")
    
    peaks_found$effect_size_in_qtl_position <-NA
    peaks_found$sum_f1_effect_sizes_in_qtl_range <-NA
    peaks_found$qtl_marker_name <-NA
    peaks_found$start_marker <-NA
    peaks_found$end_marker <-NA
    peaks_found$marker_effect_size <- NA
    
    ### Nov 2023 fix to avoid the _c _2 in marker names so that it can be treated as a number ##
    peaks_found$ci_lo <- gsub("_.*", "", peaks_found$ci_lo )
    peaks_found$ci_hi <- gsub("_.*", "", peaks_found$ci_hi )
    ### end of fix ###
    
    ## get the markers in the range ##
    for(p in 1:nrow(peaks_found)){
      peaks_found[p,]$ci_lo
      peaks_found[p,]$ci_hi
      peaks_found[p,]$lod

      peaks_found[p,]$qtl_marker_name <- map_file[map_file$pos %in% peaks_found[p,]$pos &   map_file$chr %in% peaks_found[p,]$chr,]$marker
      peaks_found[p,]$start_marker    <- map_file[map_file$pos %in% peaks_found[p,]$ci_lo & map_file$chr %in% peaks_found[p,]$chr,]$marker
      peaks_found[p,]$end_marker      <- map_file[map_file$pos %in% peaks_found[p,]$ci_hi & map_file$chr %in% peaks_found[p,]$chr,]$marker
      
      qtl_position <- map_file[map_file$pos %in% peaks_found[p,]$pos & map_file$chr %in% peaks_found[p,]$chr,]
      qtl_effect_size <-  effect_file[effect_file$marker %in% qtl_position$marker,]$effect_size
      peaks_found[p,]$effect_size_in_qtl_position <- qtl_effect_size
      start_range <- as.numeric(rownames(map_file[map_file$pos == peaks_found[p,]$ci_lo & map_file$chr %in% peaks_found[p,]$chr,]))
      end_range   <- as.numeric(rownames(map_file[map_file$pos == peaks_found[p,]$ci_hi & map_file$chr %in% peaks_found[p,]$chr,]))
      markers_in_range <- map_file[start_range:end_range,]
      # changed Nov 2023    # markers_in_range <- map_file[map_file$pos >= peaks_found[p,]$ci_lo & map_file$pos <= peaks_found[p,]$ci_hi & map_file$chr %in% peaks_found[p,]$chr,]
      
      ## take the genotypes of markers in the range and calculate the effect size in the region ##
      curr_start_in_map <- findMarkerIndexInMap(map.file, peaks_found[p,]$start_marker)
      curr_end_in_map <- findMarkerIndexInMap(map.file, peaks_found[p,]$end_marker)
      f1_geno_at_qtl_range <- substr((f1_info$G), curr_start_in_map*2-1, curr_end_in_map*2)

      effect_sizes_replicated <-  rep(effect_file[effect_file$marker %in% markers_in_range$marker,]$effect_size, each=2)# replicate each effect size to match the genotype string
      peaks_found[p,]$sum_f1_effect_sizes_in_qtl_range  <- sum(as.numeric(unlist(strsplit(f1_geno_at_qtl_range, ""))) *effect_sizes_replicated)
      
      ## calculate coef, per marker(by the chromosome)
      ## get the name of the marker and calculate the effect size out of the coef table
      curr_chr <- peaks_found[p,]$chr
      curr_m_name <- filtered_map_as_txt[(filtered_map_as_txt$pos %in% peaks_found[p,]$pos) & (filtered_map_as_txt$chr %in% curr_chr),]$marker
      coef <- scan1coef(pr[,curr_chr], my_data_files$pheno)#, addcovar=covar)
      #     peaks_found$marker_effect_size <- coef[rownames(coef) %in% curr_m_name,"11"] - coef[rownames(coef) %in% curr_m_name,"00"] ## this was run but is incorrect- the effect of the last is given to all
      peaks_found[p,]$marker_effect_size <- (coef[rownames(coef) %in% curr_m_name,"11"] - coef[rownames(coef) %in% curr_m_name,"00"])/2 ## Jan-12 fix. divide in two to get the effect per single allele A, not AA
      # plot_coef(coef, map = map,scan1_output = out)
    }
    sorted_peaks_found <- peaks_found[order(peaks_found$lod),]
    write.csv(sorted_peaks_found, file = peaks_effect_sizes)
    
    sort(unique(effect_file$effect_size))
    summary(effect_file$effect_size)
    
    
    ### continue breeding without running qtl mapping (read qtl table from file)
    } else{ ## option to skip the qtl mapping, read the qtl results from file and continue with the breeding
    peaks_found <- read.csv(peaks_effect_sizes)
    ### make sure to take qtls with positive effect size ###
    peaks_found_positive_effect <- peaks_found[peaks_found$marker_effect_size > 0,]
    if(!nrow(peaks_found_positive_effect) >0){
      clear.simdata()
      next
    }
    
    map_file <- read.table(file = map.file, header = T)

    ### make sure to take qtls with positive effect size ###
    peaks_found_positive_effect <- peaks_found[peaks_found$marker_effect_size > 0,]
    if(!nrow(peaks_found_positive_effect) >0){
      print("no peaks with positive effect size")
      clear.simdata()
      ## delete the two large files with genotypes ## 
      if (file.exists(f2_geno_file_csv)) {
        #Delete file if it exists
        file.remove(f2_geno_file_csv)
      }
      if (file.exists(f2_geno_file)) {
        #Delete file if it exists
        file.remove(f2_geno_file)
      }
      next
    }
    
    ## get the best qtl (highest LOD score) ##
    pos_best_qtl <- peaks_found_positive_effect[peaks_found_positive_effect$lod == max(peaks_found_positive_effect$lod),]$pos
    chr_best_qtl <- peaks_found_positive_effect[peaks_found_positive_effect$lod == max(peaks_found_positive_effect$lod),]$chr
    print(c("pos_best_qtl:", pos_best_qtl))
    print(c("chr_best_qtl:", chr_best_qtl))
    
    ## the best qtl ## 
    curr_marker_name <- map_file[(map_file$pos %in% pos_best_qtl) & (map_file$chr %in% chr_best_qtl),]$marker ## Jan-16-24: not sure why the filtered map is needed and why the full map is not good for this

    print(c("curr_marker_name:", curr_marker_name))
    curr_marker_ind_in_map <- findMarkerIndexInMap(map.file, curr_marker_name) 
    print(c("curr_marker_ind_in_map", curr_marker_ind_in_map))
    
    ##  F2 ##
    FitnessFromInfO_DivideToGeno(f2_info, curr_marker_ind_in_map, "f2", seed_num, pred_lr)
    f2_mean_perc_landrace <- getPercentLandraceFromInfo(f2_info, curr_marker_ind_in_map, "f2") ## per genotype at marker
    all_gener_percent_landrace <- dplyr::bind_rows(all_gener_percent_landrace, as.data.frame(f2_mean_perc_landrace))
    
    ## get the f2 individuals that are homozygous to the best qtl  ##  
    all_f2_homo_to_qtl <- getIndsWithMarkerFromInfo(f2_info, curr_marker_index = curr_marker_ind_in_map, state = "homo")
    if(!length(all_f2_homo_to_qtl$Index) > 0){
      print("no f2 homozygotes to marker found")
      write(sprintf("no individuals homozygotes to the best qtl %s in f2",curr_marker_name), file = sprintf("qtl_no_homozygotes_to_marker_%s_%s",curr_marker_name, pred_lr) ) ## added Jan-17-24
      next
    }
    all_f2_homo_to_qtl_group <- make.group(all_f2_homo_to_qtl$Index) ## added Dec-4-2023
    all_f2_homo_to_qtl_group_info <- getGroupInfoAsDF(all_f2_homo_to_qtl_group) ## Jan-22-24
    save.genotypes(filename = f2_homo_to_marker_geno_file, group = all_f2_homo_to_qtl_group)## Jan-25-24
    
    
    ## take the one f2 individual homozygous to the marker and with **highest fitness** ##
    #  f2_homo_to_qtl <- all_f2_homo_to_qtl[all_f2_homo_to_qtl$GEBV %in% max(all_f2_homo_to_qtl$GEBV),][1,] ## Dec 7 2023- max GEBV does not correspond to max fitness...
    # f2_homo_to_qtl <- all_f2_homo_to_qtl[all_f2_homo_to_qtl$fitness %in% max(all_f2_homo_to_qtl$fitness),][1,] ## Dec 8
    # f2_homo_to_qtl_index <- f2_homo_fitness_and_indices[f2_homo_fitness_and_indices$fitness %in% max(f2_homo_fitness_and_indices$fitness),]$index
    
    # f2_homo_to_qtl <- getIndsWithMarkerFromInfo(f2_info, curr_marker_index = curr_marker_ind_in_map, state = "homo")[1,]  ## take the first appearance ## function from script "functions_script_breeding.R"
    
    ### get the heterozygotes to the qtl ###
    all_f2_hetero_to_qtl <- getIndsWithMarkerFromInfo(f2_info, curr_marker_index = curr_marker_ind_in_map, state = "hetero")
    all_f2_hetero_to_qtl_group <- make.group(all_f2_hetero_to_qtl$Index) ## added Dec-4-2023
    save.genotypes(filename = all_f2_hetero_to_marker_geno_file, group = all_f2_hetero_to_qtl_group)  
    
    ## select random hetero f2 to continue the breeding with  ## Jan-10-24
    f2_hetero_to_qtl<- all_f2_hetero_to_qtl[1,] # randomly take the first occurrence
    
    if(is.na(f2_hetero_to_qtl$Index)){
      # while(nrow(f2_homo_to_qtl) < 1){
      write(sprintf("no individuals heterozygotes to the best qtl %s in f2",curr_marker_name), file = sprintf("qtl_no_heterozygotes_to_marker_%s",curr_marker_name) )
      print(sprintf("no individuals heterozygotes to the best qtl %s in f2",curr_marker_name)) ## this could point on a problem ##
      clear.simdata()
      ## delete the two large file with genotypes ## 
      if (file.exists(f2_geno_file_csv)) {
        #Delete file if it exists
        file.remove(f2_geno_file_csv)
      }
      if (file.exists(f2_geno_file)) {
        #Delete file if it exists
        file.remove(f2_geno_file)
      }
      next
    }

    ## continue breeding ##  bc1-bc1s4 should be similar to GEA breeding scheme stages #
    
    print(c("index of f2 hetero to qtl to continue in breeding scheme:", f2_hetero_to_qtl$Index))
    f2_hetero_to_qtl_group <- make.group(as.integer(f2_hetero_to_qtl$Index)) ## the selected f2 hetero individual to continue the breeding with
    f2_hetero_to_qtl_group_info <- getGroupInfoAsDF(f2_hetero_to_qtl_group)
    substr((f2_hetero_to_qtl_group_info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)

    ## save genotype of the selected f2 heterozygous to the marker ##
    save.genotypes(filename = f2_hetero_to_marker_geno_file, group = f2_hetero_to_qtl_group)  
    SaveGEBVandFitness(f2_hetero_to_qtl_group, "random_f2_hetero_to_qtl", seed_num, pred_lr) ## added Dec-4-2023 # changed to hetero Jan-10-24
    see.existing.groups()
    ####

    ##  selected F2 to continue with ##
    selected_hetero_f2_mean_perc_landrace <- getPercentLandraceFromInfo(f2_hetero_to_qtl_group_info, curr_marker_ind_in_map, "selected_hetero_f2")
    all_gener_percent_landrace <- bind_rows(all_gener_percent_landrace, as.data.frame(selected_hetero_f2_mean_perc_landrace))
    
    f2_without_the_qtl <- getIndsWithMarkerFromInfo(f2_info, curr_marker_index = curr_marker_ind_in_map, state = "without_marker")
    if(length(f2_without_the_qtl$Index) > 0){
      f2_without_the_qtl_group <- make.group(f2_without_the_qtl$Index) ## added Dec-4-2023
      save.genotypes(filename = f2_without_marker_geno_file, group = f2_without_the_qtl_group)
    }
    
    ####
    ## make bc1 ##
    #  bc1 <- cross.combinations(first.parents = elite_line_id ,second.parents = new_f2_group)
    # bc1 <- cross.combinations(first.parents = elite_line_id ,second.parents = f2_homo_to_qtl$Index)
    # bc1 <- cross.combinations(first.parents = elite_line_id ,second.parents = f2_homo_to_qtl_group_info$Index) ## fix 13.6.23 - seem to work well. 
    bc1 <- cross.combinations(first.parents = elite_line_id ,second.parents = f2_hetero_to_qtl_group_info$Index) ## Jan-10-24 - take random hetero f2 to continue 
    
    ## make sure bc1 is hetero and self it##
    bc1_group_info     <- data.frame(Index=see.group.data(bc1,data.type = "X"), GEBV=see.group.data(bc1,data.type = "BV"), G=see.group.data(bc1,data.type ="G"))
    bc1_geno_at_marker <- substr((bc1_group_info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)
    if(!bc1_geno_at_marker %in% c("21", "12","02","20")){  ## c("10", "01", "21", "12")
      print("bc1 not heterozygous.") ## this should not happen because parent is a homozygote.>> Jan-10-24 - the f2 parent is now hetero so it could give bc1 without the marker. need to make sure bc1 is hetero
      print(bc1_geno_at_marker)
      
      while(!bc1_geno_at_marker %in% c("21", "12","02","20") ){  ## c("10", "01", "21", "12"),##  until I get a heterozygote ## added  Jan-10-24 
        delete.group(bc1)
        print("bc1 deleted")
        bc1 <- cross.combinations(first.parents = elite_line_id ,second.parents = f2_hetero_to_qtl_group_info$Index, give.names = TRUE, name.prefix = "bc1")
        bc1_group_info <- getGroupInfoAsDF(bc1)
        bc1_geno_at_marker <- substr((bc1_group_info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)
        print(bc1_geno_at_marker)
      }
    } 

    ##  bc1  ##
    FitnessFromInfO_DivideToGeno(bc1_group_info, curr_marker_ind_in_map, "bc1", seed_num, pred_lr)
    bc1_mean_perc_landrace <- getPercentLandraceFromInfo(bc1_group_info, curr_marker_ind_in_map, "bc1")
    all_gener_percent_landrace <- bind_rows(all_gener_percent_landrace, as.data.frame(bc1_mean_perc_landrace))
    save.genotypes(filename = bc1_geno_file, group = bc1)  
    
    print("after bc1")
    see.existing.groups()
    
    ## bc1s1 ##
    bc1S1 <- ChooseIndividualByStateAndSelf(bc1_group_info, curr_marker_ind_in_map, number_of_offsprings_bc1S1_S2_S3, "hetero")
    bc1S1_group_info <- getGroupInfoAsDF(bc1S1)
    substr((bc1S1_group_info$G),curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)[1:10]
    FitnessFromInfO_DivideToGeno(bc1S1_group_info, curr_marker_ind_in_map, "bc1s1", seed_num, pred_lr)
    bc1S1_mean_perc_landrace <- getPercentLandraceFromInfo(bc1S1_group_info, curr_marker_ind_in_map, "bc1S1")
    all_gener_percent_landrace <- bind_rows(all_gener_percent_landrace, as.data.frame(bc1S1_mean_perc_landrace))
    
    bc1S1_hetero_to_qtl <- getIndsWithMarkerFromInfo(bc1S1_group_info, curr_marker_index = curr_marker_ind_in_map, state = "hetero")## added Dec-7-2023
    bc1S1_hetero_to_qtl_group <- make.group(bc1S1_hetero_to_qtl$Index) ## added Dec-7-2023
    save.genotypes(filename = bc1s1_hetero_to_marker_geno_file, group = bc1S1_hetero_to_qtl_group)
 
    bc1S1_homo_to_qtl <- getIndsWithMarkerFromInfo(bc1S1_group_info, curr_marker_index = curr_marker_ind_in_map, state = "homo")## added Dec-7-2023
    if(length(bc1S1_homo_to_qtl$Index) > 0){
      bc1S1_homo_to_qtl_group <- make.group(bc1S1_homo_to_qtl$Index) ## added Dec-7-2023
      save.genotypes(filename = bc1s1_homo_to_marker_geno_file, group = bc1S1_homo_to_qtl_group)
    }else{
      print("no bc1S1_homo_to_qtl")
    }

    bc1S1_without_qtl <- getIndsWithMarkerFromInfo(bc1S1_group_info, curr_marker_index = curr_marker_ind_in_map, state = "without_marker")## added Dec-7-2023
    if(length(bc1S1_without_qtl$Index) > 0){
      bc1S1_without_qtl_group <- make.group(bc1S1_without_qtl$Index) ## added Dec-7-2023
      save.genotypes(filename = bc1s1_without_marker_geno_file, group = bc1S1_without_qtl_group)
    }else{
      print("no bc1S1_without_qtl_group")
    }
  
    ### BC1S2 ####
    bc1S2 <- ChooseIndividualByStateAndSelf(bc1S1_group_info, curr_marker_ind_in_map, number_of_offsprings_bc1S1_S2_S3, "hetero")
    bc1S2_group_info <- getGroupInfoAsDF(bc1S2)
    substr((bc1S2_group_info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)[1:10]
    FitnessFromInfO_DivideToGeno(bc1S2_group_info, curr_marker_ind_in_map, "bc1s2", seed_num, pred_lr)
    bc1S2_mean_perc_landrace <- getPercentLandraceFromInfo(bc1S2_group_info, curr_marker_ind_in_map, "bc1s2")
    all_gener_percent_landrace <- bind_rows(all_gener_percent_landrace, as.data.frame(bc1S2_mean_perc_landrace))
    
    bc1S2_homo_to_qtl <- getIndsWithMarkerFromInfo(bc1S2_group_info, curr_marker_index = curr_marker_ind_in_map, state = "homo")## added Dec-7-2023
    if(length(bc1S2_homo_to_qtl$Index) > 0){
      bc1S2_homo_to_qtl_group <- make.group(bc1S2_homo_to_qtl$Index) ## added Dec-7-2023
      save.genotypes(filename = bc1s2_homo_to_marker_geno_file, group = bc1S2_homo_to_qtl_group)
    }else{
      print("no bc1S2_homo_to_qtl_group")
    }
    #
    bc1S2_hetero_to_qtl <- getIndsWithMarkerFromInfo(bc1S2_group_info, curr_marker_index = curr_marker_ind_in_map, state = "hetero")## added Dec-7-2023
    if(length(bc1S2_hetero_to_qtl$Index) > 0){
      bc1S2_hetero_to_qtl_group <- make.group(bc1S2_hetero_to_qtl$Index) ## added Dec-7-2023
      save.genotypes(filename = bc1s2_hetero_to_marker_geno_file, group = bc1S2_hetero_to_qtl_group)
    }else{
      print("no bc1S2_hetero_to_qtl_group")
    }
    #
    bc1S2_wo_qtl <- getIndsWithMarkerFromInfo(bc1S2_group_info, curr_marker_index = curr_marker_ind_in_map, state = "without_marker")## added Dec-7-2023
    if(length(bc1S2_wo_qtl$Index) > 0){
      bc1S2_wo_qtl_group <- make.group(bc1S2_wo_qtl$Index) ## added Dec-7-2023
      save.genotypes(filename = bc1s2_without_marker_geno_file, group = bc1S2_wo_qtl_group)
    }else{
      print("no bc1S2_wo_qtl_group")
    }
    
    ### bc1s3 ###
    bc1S3 <- ChooseIndividualByStateAndSelf(bc1S2_group_info, curr_marker_ind_in_map, number_of_offsprings_bc1S1_S2_S3, "hetero")
    bc1S3_group_info <- getGroupInfoAsDF(bc1S3)
    substr((bc1S3_group_info$G),curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)[1:10]
    FitnessFromInfO_DivideToGeno(bc1S3_group_info, curr_marker_ind_in_map, "bc1s3", seed_num, pred_lr)
    bc1S3_mean_perc_landrace <- getPercentLandraceFromInfo(bc1S3_group_info, curr_marker_ind_in_map, "bc1s3")
    all_gener_percent_landrace <- bind_rows(all_gener_percent_landrace, as.data.frame(bc1S3_mean_perc_landrace))
    
    bc1S3_homo_to_qtl <- getIndsWithMarkerFromInfo(bc1S3_group_info, curr_marker_index = curr_marker_ind_in_map, state = "homo")## added Dec-7-2023
    if(length(bc1S3_homo_to_qtl$Index) > 0){
      bc1S3_homo_to_qtl_group <- make.group(bc1S3_homo_to_qtl$Index) ## added Dec-7-2023
      save.genotypes(filename = bc1s3_homo_to_marker_geno_file, group = bc1S3_homo_to_qtl_group)
    }else{
      print("no bc1S3_homo_to_qtl_group")
    }
    #
    bc1S3_hetero_to_qtl <- getIndsWithMarkerFromInfo(bc1S3_group_info, curr_marker_index = curr_marker_ind_in_map, state = "hetero")## added Dec-7-2023
    if(length(bc1S3_hetero_to_qtl$Index) > 0){
      bc1S3_hetero_to_qtl_group <- make.group(bc1S3_hetero_to_qtl$Index) ## added Dec-7-2023
      save.genotypes(filename = bc1s3_hetero_to_marker_geno_file, group = bc1S3_hetero_to_qtl_group)
    }else{
      print("no bc1S3_hetero_to_qtl_group")
    }
    #
    bc1S3_wo_qtl <- getIndsWithMarkerFromInfo(bc1S3_group_info, curr_marker_index = curr_marker_ind_in_map, state = "without_marker")## added Dec-7-2023
    if(length(bc1S3_wo_qtl$Index) > 0){
      bc1S3_wo_qtl_group <- make.group(bc1S3_wo_qtl$Index) ## added Dec-7-2023
      save.genotypes(filename = bc1s3_without_marker_geno_file, group = bc1S3_wo_qtl_group)
    }else{
      print("no bc1S3_wo_qtl_group")
    }
    
    ### bc1s4 ###
    bc1S4 <- ChooseIndividualByStateAndSelf(bc1S3_group_info, curr_marker_ind_in_map, bc1s4_offsprings_num, "hetero") ## 100 offsprings 
    bc1S4_group_info <- getGroupInfoAsDF(bc1S4)
    substr((bc1S4_group_info$G),curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)[1:5]
    FitnessFromInfO_DivideToGeno(bc1S4_group_info, curr_marker_ind_in_map, "bc1s4", seed_num, pred_lr)
    bc1S4_mean_perc_landrace <- getPercentLandraceFromInfo(bc1S4_group_info, curr_marker_ind_in_map, "bc1s4")
    all_gener_percent_landrace <- bind_rows(all_gener_percent_landrace, as.data.frame(bc1S4_mean_perc_landrace))
    

    ## save info of bc1s4 homozygous to the marker ##
    bc1S4_homozygous_to_marker <- getIndsWithMarkerFromInfo(bc1S4_group_info, curr_marker_ind_in_map, "homo") ##
    if(length(bc1S4_homozygous_to_marker$Index) > 0){
      new_group_bc1S4_homozygous_to_marker <- make.group(bc1S4_homozygous_to_marker$Index) 
      save.genotypes(filename = bc1s4_homo_to_marker_geno_file, group = new_group_bc1S4_homozygous_to_marker)
    }else{
      print("no new_group_bc1S4_homozygous_to_marker")
    }
    
    homozygous_bc1S4_group_info <- getGroupInfoAsDF(new_group_bc1S4_homozygous_to_marker)
    substr((homozygous_bc1S4_group_info$G), curr_marker_ind_in_map*2-1, curr_marker_ind_in_map*2)
    fitness_homozygous_bc1S4    <- fitness_function(homozygous_bc1S4_group_info$GEBV)
    print(see.existing.groups())
    bc1s4_res <- cbind(fitness_homozygous_bc1S4, pred_lr, curr_marker_name)
    ## append QTL analysis results for all predicted individuals results file
    write.table(bc1s4_res, file = qtl_out_file, sep = ",", append = T, col.names = F, row.names = F)

    
    bc1S4_hetero_to_qtl <- getIndsWithMarkerFromInfo(bc1S4_group_info, curr_marker_index = curr_marker_ind_in_map, state = "hetero")## added Dec-7-2023
    if(length(bc1S4_hetero_to_qtl$Index) > 0){
      bc1S4_hetero_to_qtl_group <- make.group(bc1S4_hetero_to_qtl$Index) ## added Dec-7-2023
      save.genotypes(filename = bc1s4_hetero_to_marker_geno_file, group = bc1S4_hetero_to_qtl_group)
    }else{
      print("no bc1S4_hetero_to_qtl_group")
    }
    #
    bc1S4_wo_qtl <- getIndsWithMarkerFromInfo(bc1S4_group_info, curr_marker_index = curr_marker_ind_in_map, state = "without_marker")## added Dec-7-2023
    if(length(bc1S4_wo_qtl$Index) > 0){
      bc1S4_wo_qtl_group <- make.group(bc1S4_wo_qtl$Index) ## added Dec-7-2023
      save.genotypes(filename = bc1s4_without_marker_geno_file, group = bc1S4_hetero_to_qtl_group)
    }else{
      print("no bc1S4_wo_qtl_group")
    }
    #######################
    clear.simdata()

  }## end of else - remove after re-running all seeds with the fix of taking only markers with positive effect size
  clear.simdata()
  
  # colnames(all_gener_percent_landrace) <- c("generation", "mean_percent_landrace", "std")
  write.csv(all_gener_percent_landrace, file = sprintf("mean_percent_landrace_%s_%s.csv",seed_num, pred_lr ))
}
## QTL analysis results for all predicted individuals
# write.table(all_bc1S4_results, file = qtl_out_file , sep = ",")

end_time <- Sys.time()
print(end_time)
print(c("time (min): ", (end_time-start_time)/60))

########################