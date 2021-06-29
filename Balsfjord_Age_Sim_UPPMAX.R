#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se

argv <- commandArgs(trailingOnly = TRUE)
chr <- argv[1]

save(Ch_v2_linkage_map, Ch_v2_recombination_profile, Ch_v2.0.2_sizes, intro_20k_lists, file = "~/Projects/Herring/data/Balsfjord/age_sim_ref_data.RData")
source("./Balsfjord_Age_Sim_UPPMAX_fun.R")
load("./age_sim_ref_data.RData")

target_intro_GR <- unlist(intro_20k_lists$atl_red)
#atl_GW_results <- list()
pac_target_intro_GR <- unlist(intro_20k_lists$pac_red)

print(paste0("Starting chr: ", chr))
flush.console()
chr_GR <- target_intro_GR[target_intro_GR@seqnames == chr]
pac_chr_GR <- pac_target_intro_GR[pac_target_intro_GR@seqnames == chr]
atl_tmp <- chr_age_estimation(chr = chr, chr_target_GR = chr_GR, gen_map = Ch_v2_linkage_map, gen_prof = Ch_v2_recombination_profile, size_df = Ch_v2.0.2_sizes, n_runs = 1000, n_generations = 1e6)
assign(paste("atl_age", chr, sep = "_"),atl_tmp )
pac_tmp <-  chr_age_estimation(chr = chr, chr_target_GR = pac_chr_GR, gen_map = Ch_v2_linkage_map, gen_prof = Ch_v2_recombination_profile, size_df = Ch_v2.0.2_sizes, n_runs = 1000, n_generations = 1e6)
assign(paste("pac_age", chr, sep = "_"),pac_tmp )

save(list=c(paste("pac_age", chr, sep = "_"), paste("atl_age", chr, sep = "_")), file = paste0("./", chr, "_age_est_UPPMAX.RData"))
