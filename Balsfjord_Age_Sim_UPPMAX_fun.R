#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se

simulate_linkage_breakdown_GW <- function(sim_linkage_prof, n_gen = 100){
  f1 <- Vectorize(FUN = function(x,y){min(which(y > x))}, vectorize.args = "x")
  gen_map_max <- max(sim_linkage_prof$interval_dist)
  rec_vec <- runif(n=n_gen, max = max(c(gen_map_max, 100)))
  #Core simulation block
  rec_df <- data.frame(rec_pos = 1:n_gen, lower_bp_tmp = NA, upper_bp_tmp = NA, size = NA, stringsAsFactors = F)
  rec_df$upper_bp_tmp[1] <- sum(sim_linkage_prof$interval_size)
  rec_df$lower_bp_tmp[1] <- 0
  rec_df$rec_pos <- NA 
  rec_df$rec_event <- (rec_vec <= gen_map_max)
  rec_int_vec <- f1(x = rec_vec[rec_df$rec_event], y = sim_linkage_prof$interval_dist) 
  rec_df$rec_interval[rec_df$rec_event] <- rec_int_vec 
  rec_df$rec_event_shift[rec_df$rec_event] <- sample.int(1e6, sum(rec_df$rec_event), replace = T)/1e6
  rec_df$int_start[rec_df$rec_event] <- sim_linkage_prof[rec_int_vec,"BIN_START"]
  rec_df$int_size[rec_df$rec_event] <- sim_linkage_prof[rec_int_vec,"interval_size"]
  rec_df$rec_pos <- rec_df$int_start + round(rec_df$int_size*rec_df$rec_event_shift)
  return(invisible(rec_df))
}

age_per_pos <- function(focus_pos, focus_size, rec_df){
  rec_df$before_focus <- rec_df$rec_pos <= focus_pos
  rec_df$after_focus <- rec_df$rec_pos > focus_pos
  
  rec_df$upper_bp_tmp[which(rec_df$after_focus & rec_df$rec_event)] <- cummin(rec_df$rec_pos[which(rec_df$after_focus & rec_df$rec_event)])
  rec_df$upper_bp <- rep(rec_df$upper_bp_tmp[!is.na(rec_df$upper_bp_tmp)], times = diff(c(which(!is.na(rec_df$upper_bp_tmp)),(length(rec_df$upper_bp_tmp)+1))))
  rec_df$lower_bp_tmp[which(rec_df$before_focus & rec_df$rec_event)] <- cummax(rec_df$rec_pos[which(rec_df$before_focus & rec_df$rec_event)])
  rec_df$lower_bp <- rep(rec_df$lower_bp_tmp[!is.na(rec_df$lower_bp_tmp)], times = diff(c(which(!is.na(rec_df$lower_bp_tmp)),(length(rec_df$lower_bp_tmp)+1))))
  rec_df$size <- rec_df$upper_bp - rec_df$lower_bp
  age_est <- which.min(abs(focus_size - rec_df$size))
  if(all(rec_df$size > focus_size)) age_est <- 2*length(rec_df$size)
  return(age_est)
}

age_pos_vec <- Vectorize("age_per_pos",vectorize.args = c("focus_pos", "focus_size"))

chr_age_estimation <- function(chr, chr_target_GR, gen_map, gen_prof, size_df, n_generations = 2000, n_runs = 100){
  reg_linkage_max <- max(gen_map[gen_map$chr == sub("chr", "", chr),]$all_map)
  phys_size <- size_df[size_df$name == chr,"size"]
  #Evening out the rates across each "step" of the map
  #Adding refinement by using recombination profile scaled by length of map per chromosome
  reg_prof <- gen_prof[gen_prof$CHR ==  sub("chr", "", chr),]
  reg_prof[,"interval_dist"] <- cumsum((reg_prof$RR/sum(reg_prof$RR, na.rm = T)) * reg_linkage_max)
  reg_prof$interval_size <- diff(c(reg_prof$BIN_START, phys_size))
  
  chr_pos_vec <- chr_target_GR@ranges@start + chr_target_GR@ranges@width/2
  chr_size_vec <- chr_target_GR@ranges@width
  age_est_df <- data.frame(hap = names(chr_target_GR), gr = paste(chr_target_GR), stringsAsFactors = F)
  for (i in 1:n_runs){
    chr_rec_df <- simulate_linkage_breakdown_GW(sim_linkage_prof = reg_prof, n_gen = n_generations)
    chr_age_vec <-  age_pos_vec(chr_pos_vec, chr_size_vec, chr_rec_df)
    age_est_df[,paste("age", i, sep = "_")] <- chr_age_vec
  }
  age_est_df$median_age <- apply(age_est_df[, grepl("age_", names(age_est_df))],1, "median")
  return(invisible(age_est_df))
}