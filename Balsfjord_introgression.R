#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se


#Example workflow:
require(stringdist)

#Update for current genome (v2.0.2)
load("~/Projects/Herring/data/v2.0.2_genotypes/herring_79.RData")
load("~/Projects/Herring/data/v2.0.2_genotypes/h79_sw.RData")
hws6_dist_v2.0.2_df <- introgression_contrast(target="HWS6", ref_1 = "A[MF][0-9]", ref_2 = "HWS1", w_df=h79_w, sample_list=herring_79$sample_list, geno=herring_79$geno)
save(hws6_dist_v2.0.2_df, file = "~/Projects/Herring/data/v2.0.2_genotypes/hws6_introgression_df_v2.0.2.Rdata")
hws6_v2.0.2_intro <- introgression_plot(hws6_dist_v2.0.2_df, sample_list=herring_79$sample_list, snp_numbers=h79_s$n_snps, snp_cutoff = 500)
#Sanity check of intorgression
which(h79_w$scaffold == "chr6" & h79_w$start == 18500001)
#[1] 1813
plot_top_regions_2(top_pos = 1813, win_df = h79_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/Balsfjord/", snp_col = "blue")
#Checks out
which(h79_w$scaffold == "chr21" & h79_w$start == 13700001)
#[1] 6155
plot_top_regions_2(top_pos = 6155, win_df = h79_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/Balsfjord/HapDist_v2.0.2/", snp_col = "blue")

#Z-score type threshold
assoc_t <- 1/10^(mean(log10(hws6_v2.0.2_intro$dist$HWS61_Balsfjord_Atlantic_1_ratio), na.rm  =T) - qnorm(1 -(0.05/length(hws6_v2.0.2_intro$dist$HWS61_Balsfjord_Atlantic_1_ratio))) * sd(log10(hws6_v2.0.2_intro$dist$HWS61_Balsfjord_Atlantic_1_ratio), na.rm  =T))
hws6_v2.0.2_intro <- introgression_plot(hws6_dist_v2.0.2_df, sample_list=herring_79$sample_list, snp_numbers=h79_s$n_snps, pdf_file="~/Projects/Herring/doc/Balsfjord/HWS6_v2.0.2_intro.pdf", snp_cutoff = 500, assoc_tresh = assoc_t)






#Smaller windows
require(HaploDistScan)
h79_20k_w <- construct_window_annotation(herring_79$geno, interval = 2e4)
hws6_dist_v2.0.2_20k_df <- introgression_contrast(target="HWS6", ref_1 = "A[MF][0-9]", ref_2 = "HWS1", w_df=h79_20k_w, sample_list=herring_79$sample_list, geno=herring_79$geno)
save(hws6_dist_v2.0.2_20k_df, file = "~/Projects/Herring/data/v2.0.2_genotypes/hws6_introgression_df_v2.0.2_20k.Rdata")
diff_SNP_count_vec <- rowSums(hws6_dist_v2.0.2_20k_df[,6:21])/length(6:21)*2
hws6_v2.0.2_20k_intro <- introgression_plot_2(hws6_dist_v2.0.2_20k_df, sample_list=herring_79$sample_list, snp_numbers=diff_SNP_count_vec, pdf_file="~/Projects/Herring/doc/Balsfjord/HWS6_v2.0.2_20k_intro_recalc5.pdf", min_diff = 20, snp_cutoff = 50, assoc_tresh = 8, assoc_dir = "down")

#White Sea sample
hws3_dist_v2.0.2_20k_df <- introgression_contrast(target="HWS3", ref_1 = "A[MF][0-9]", ref_2 = "HWS1", w_df=h79_20k_w, sample_list=herring_79$sample_list, geno=herring_79$geno)
save(hws3_dist_v2.0.2_20k_df, file = "~/Projects/Herring/data/Balsfjord/hws3_introgression_df_v2.0.2_20k.Rdata")
hws3_diff_SNP_count_vec <- rowSums(hws3_dist_v2.0.2_20k_df[,6:21])/length(6:21)*2
hws3_v2.0.2_20k_intro <- introgression_plot_2(hws3_dist_v2.0.2_20k_df, sample_list=herring_79$sample_list, snp_numbers=hws3_diff_SNP_count_vec, pdf_file="~/Projects/Herring/doc/Balsfjord/HWS3_v2.0.2_20k_intro.pdf", min_diff = 20, snp_cutoff = 50, assoc_tresh = 8, assoc_dir = "down")

#Pechora Sea sample
hws2_dist_v2.0.2_20k_df <- introgression_contrast(target="HWS2", ref_1 = "A[MF][0-9]", ref_2 = "HWS1", w_df=h79_20k_w, sample_list=herring_79$sample_list, geno=herring_79$geno)
save(hws2_dist_v2.0.2_20k_df, file = "~/Projects/Herring/data/Balsfjord/hws2_introgression_df_v2.0.2_20k.Rdata")
hws2_diff_SNP_count_vec <- rowSums(hws2_dist_v2.0.2_20k_df[,6:21])/length(6:21)*2
hws2_v2.0.2_20k_intro <- introgression_plot_2(hws2_dist_v2.0.2_20k_df, sample_list=herring_79$sample_list, snp_numbers=hws2_diff_SNP_count_vec, pdf_file="~/Projects/Herring/doc/Balsfjord/HWS2_v2.0.2_20k_intro.pdf", min_diff = 20, snp_cutoff = 50, assoc_tresh = 8, assoc_dir = "down")

#KB Spring sample
hws4_dist_v2.0.2_20k_df <- introgression_contrast(target="HWS4", ref_1 = "A[MF][0-9]", ref_2 = "HWS1", w_df=h79_20k_w, sample_list=herring_79$sample_list, geno=herring_79$geno)
save(hws4_dist_v2.0.2_20k_df, file = "~/Projects/Herring/data/Balsfjord/hws4_introgression_df_v2.0.2_20k.Rdata")
hws4_diff_SNP_count_vec <- rowSums(hws4_dist_v2.0.2_20k_df[,6:21])/length(6:21)*2
hws4_v2.0.2_20k_intro <- introgression_plot_2(hws4_dist_v2.0.2_20k_df, sample_list=herring_79$sample_list, snp_numbers=hws4_diff_SNP_count_vec, pdf_file="~/Projects/Herring/doc/Balsfjord/HWS4_v2.0.2_20k_intro.pdf", min_diff = 20, snp_cutoff = 50, assoc_tresh = 8, assoc_dir = "down")

#KB Autumn sample
hws5_dist_v2.0.2_20k_df <- introgression_contrast(target="HWS5", ref_1 = "A[MF][0-9]", ref_2 = "HWS1", w_df=h79_20k_w, sample_list=herring_79$sample_list, geno=herring_79$geno)
save(hws5_dist_v2.0.2_20k_df, file = "~/Projects/Herring/data/Balsfjord/hws5_introgression_df_v2.0.2_20k.Rdata")
hws5_diff_SNP_count_vec <- rowSums(hws5_dist_v2.0.2_20k_df[,6:21])/length(6:21)*2
hws5_v2.0.2_20k_intro <- introgression_plot_2(hws5_dist_v2.0.2_20k_df, sample_list=herring_79$sample_list, snp_numbers=hws5_diff_SNP_count_vec, pdf_file="~/Projects/Herring/doc/Balsfjord/HWS5_v2.0.2_20k_intro.pdf", min_diff = 20, snp_cutoff = 50, assoc_tresh = 8, assoc_dir = "down")


ratio_cols <- grep("ratio", names(hws6_v2.0.2_20k_intro$dist))
#Comparing with recombination rate data
Ch_v2_recombination_profile <- read.table("~/Projects/Herring/data/HiC_assemblies/Version_2_release/Baltic.LDhat.100K.txt", stringsAsFactors=F, header = T)



#Examining the 20k version
require(GenomicRanges)
intro_20k_lists <- complile_intro_lists(intr_obj = hws6_v2.0.2_20k_intro, win_df = h79_20k_w, fuse_thresh = 5e4, assoc_t = 8)
summary(unlist(intro_20k_lists$atl_red)@ranges@width)
hist(unlist(intro_20k_lists$atl_red)@ranges@width, breaks = seq(from = 0, to = 3.5e6, by = 2e4))

summary(unlist(intro_20k_lists$pac_red)@ranges@width)
hist(unlist(intro_20k_lists$pac_red)@ranges@width, breaks = seq(from = 0, to = 3.5e6, by = 2e4))

atl_intro_cov_20k <- calc_cov(intro_20k_lists$atl_red)
pdf(width = 15, height = 5,  file = "~/Projects/Herring/doc/Balsfjord/HWS6_v2.0.2_intro_20k_cov.pdf")
plot(x = atl_intro_cov_20k$global_start, y = atl_intro_cov_20k$cov, type = "n", main = "Introgression \"Coverage\"", xlab = "", ylab = "Number of haplotypes") 
segments(x0 = atl_intro_cov_20k$global_start, x1 = atl_intro_cov_20k$global_end, y0 = atl_intro_cov_20k$cov, lwd = 4, col = atl_intro_cov_20k$col)
dev.off()

pac_intro_cov_20k <- calc_cov(intro_20k_lists$pac_red)
typical_dist <- mean(hws6_v2.0.2_20k_intro$hap$ref_1_mean/hws6_v2.0.2_20k_intro$hap$ref_2_mean)
intro_haplotype_plot(intro_obj = hws6_v2.0.2_20k_intro, rec_profile = Ch_v2_recombination_profile, 
                     p_cov = pac_intro_cov_20k, a_cov = atl_intro_cov_20k, 
                     p_thresh = typical_dist*assoc_t,  a_thresh = typical_dist/assoc_t,
                     pdf_file = "~/Projects/Herring/doc/Balsfjord/HWS6_v2.0.2__20k_haps.pdf")

#Recurring Atlantic regions
recurring_intro_df <- atl_intro_cov_20k[atl_intro_cov_20k$cov == 8,c(2:4)]
names(recurring_intro_df)[1] <- "seqnames"
recurring_intro_GR <- GRanges(recurring_intro_df)
rec_gtf <- cluhar_v2.0.2_gtf[unique(findOverlaps(recurring_intro_GR, cluhar_v2.0.2_gtf)@to)]
rec_gtf_flank <- cluhar_v2.0.2_gtf[unique(findOverlaps(recurring_intro_GR, cluhar_v2.0.2_gtf, maxgap = 5e4)@to)]
rec_genes <- rec_gtf_flank[rec_gtf_flank$type == "gene"]
rec_genes$intro_type = "Atlantic"
rec_genes$location_type = "flanking_50kb"
rec_genes$location_type[rec_genes$gene_id %in% rec_gtf$gene_id] <- "internal"
rec_exons <- rec_gtf[rec_gtf$type == "exon"]
rec_exon_footprint <- reduce(rec_exons, ignore.strand = T)
#Exon content of introgressed regions
sum(rec_exon_footprint@ranges@width)/sum(recurring_intro_GR@ranges@width)
#[1] 0.08214032

genomic_exon_footprint <- reduce(cluhar_v2.0.2_gtf[cluhar_v2.0.2_gtf$type == "exon"], ignore.strand = T)
#Exon content of the chromosomes
sum(genomic_exon_footprint@ranges@width)/sum(Ch_v2.0.2@ranges@width[1:26])
#[1] 0.0890693
exp_n_exons <- round(length(genomic_exon_footprint)*(sum(recurring_intro_GR@ranges@width)/sum(Ch_v2.0.2_sizes$size[1:26])))
chisq.test(x = matrix(data = c(length(rec_exon_footprint), exp_n_exons), ncol = 2))

#Recurring Pacific regions
recurring_intro_pac_df <- pac_intro_cov_20k[pac_intro_cov_20k$cov == 8,c(2:4)]
names(recurring_intro_pac_df)[1] <- "seqnames"
recurring_intro_pac_GR <- GRanges(recurring_intro_pac_df)
rec_pac_gtf <- cluhar_v2.0.2_gtf[unique(findOverlaps(recurring_intro_pac_GR, cluhar_v2.0.2_gtf)@to)]
rec_pac_gtf_flank <- cluhar_v2.0.2_gtf[unique(findOverlaps(recurring_intro_pac_GR, cluhar_v2.0.2_gtf, maxgap = 5e4)@to)]
rec_pac_genes <- rec_pac_gtf_flank[rec_pac_gtf_flank$type == "gene"]
rec_pac_genes$intro_type = "Pacific"
rec_pac_genes$location_type = "flanking_50kb"
rec_pac_genes$location_type[rec_pac_genes$gene_id %in% rec_pac_gtf$gene_id] <- "internal"
rec_pac_exons <- rec_pac_gtf[rec_pac_gtf$type == "exon"]
rec_pac_exon_footprint <- reduce(rec_pac_exons, ignore.strand = T)
exp_n_exons_pac <- round(length(genomic_exon_footprint)*(sum(recurring_intro_pac_GR@ranges@width)/sum(Ch_v2.0.2_sizes$size[1:26])))

#Exon content of introgressed regions
sum(rec_pac_exon_footprint@ranges@width)/sum(recurring_intro_pac_GR@ranges@width)
#[1] 0.1574758

chisq.test(x = matrix(data = c(length(rec_pac_exon_footprint), exp_n_exons_pac), ncol = 2))

save(rec_pac_genes, rec_genes, file = "~/Projects/Herring/data/Balsfjord/genes_in_recurring_intro_regions.Rdata")




#Writing a list with all detected introgresissions
atl_20k_intro_unlist <- unlist(intro_20k_lists$atl_red)
atl_20k_intro_unlist$intro_type = "Atlantic"
atl_20k_intro_unlist$haplo_id = names(atl_20k_intro_unlist)
pac_20k_intro_unlist <- unlist(intro_20k_lists$pac_red)
pac_20k_intro_unlist$intro_type = "Pacific"
pac_20k_intro_unlist$haplo_id = names(pac_20k_intro_unlist)


write.table(x = c(atl_20k_intro_unlist, pac_20k_intro_unlist), quote = F, col.names = T, row.names = F, sep = "\t", file = "~/Projects/Herring/data/Balsfjord/intro_20k_all.txt")

#Plotting recurring windows
atl_intro_cov_20k_no_red <- calc_cov(intro_20k_lists$atl)
h79_20k_w_GR <- GRanges(seqnames = h79_20k_w$scaffold, ranges = IRanges(start = h79_20k_w$start, end =  h79_20k_w$stop))
rec_win_df <- atl_intro_cov_20k_no_red[atl_intro_cov_20k_no_red$cov >= 8,]
rec_win_GR <- GRanges(seqnames = rec_win_df$group_name, ranges = IRanges(start = rec_win_df$start, end =  rec_win_df$end))
#rec_win_hits <- findOverlaps(h79_20k_w_GR, rec_win_GR)
rec_win_hits <- findOverlaps(h79_20k_w_GR, recurring_intro_GR) #Using reduced version
h79_20k_w_GR[rec_win_hits@from]
plot_top_regions_2(top_pos = rec_win_hits@from, win_df = h79_20k_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/Balsfjord/HapDist_v2.0.2/Recurring_windows/", snp_col = "blue")

#attempting to switch clade order
weight_vec <- rep(10, length(herring_79$sample_list))
weight_vec[grep("Pacific", herring_79$sample_list)] <- 0.1
weight_vec[grep("HWS", herring_79$sample_list)] <- 5
weight_vec[grep("HWS14_Japan_SeaofJapan_2|HWS13_Japan_SeaofJapan_2", herring_79$sample_list)] <- 50
weight_vec <- -weight_vec
#rec_win_HapDist <- plot_top_regions_3(top_pos = rec_win_hits@from, win_df = h79_20k_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/Balsfjord/HapDist_v2.0.2/Reorder/", snp_col = "blue", dend_weight_vec = weight_vec)
rec_win_HapDist <- plot_top_regions_3(top_pos = rec_win_hits@from, win_df = h79_20k_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/Balsfjord/HapDist_v2.0.2/Rec_win_div/", snp_col = "blue", dend_weight_vec = weight_vec)

#Using the distance matrix to get within-group diveristy
Atl_haps <- grep("Atlantic|Baltic", herring_79$sample_list, value = T)
Atl_haps <- grep("Balsfjord", Atl_haps, value = T, invert = T)
HWS_haps <- grep("HWS", herring_79$sample_list, value = T)
HWS_haps <- grep("Balsfjord", HWS_haps, value = T, invert = T)

atl_div_vec <- numeric()
HWS_div_vec <- numeric()
dxy_vec <- numeric()
for(reg_loc in names(rec_win_HapDist$dist_mat)){
  atl_div_vec[reg_loc] <- mean(rec_win_HapDist$dist_mat[[reg_loc]][Atl_haps, Atl_haps], na.rm = T)/2e4
  HWS_div_vec[reg_loc] <- mean(rec_win_HapDist$dist_mat[[reg_loc]][HWS_haps, HWS_haps], na.rm = T)/2e4
  dxy_vec[reg_loc] <- mean(rec_win_HapDist$dist_mat[[reg_loc]][Atl_haps, HWS_haps], na.rm = T)/2e4
}


pac_intro_cov_20k_no_red <- calc_cov(intro_20k_lists$pac)
pac_rec_win_df <- pac_intro_cov_20k_no_red[pac_intro_cov_20k_no_red$cov >= 8,]
pac_rec_win_GR <- GRanges(seqnames = pac_rec_win_df$group_name, ranges = IRanges(start = pac_rec_win_df$start, end =  pac_rec_win_df$end))
#pac_rec_win_hits <- findOverlaps(h79_20k_w_GR, pac_rec_win_GR)
pac_rec_win_hits <- findOverlaps(h79_20k_w_GR, recurring_intro_pac_GR) #Using reduced version
h79_20k_w_GR[pac_rec_win_hits@from]
plot_top_regions_2(top_pos = pac_rec_win_hits@from, win_df = h79_20k_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/Balsfjord/HapDist_v2.0.2/Recurring_windows_pac/", snp_col = "blue")
#Pac_rec_win_HapDist <- plot_top_regions_3(top_pos = 5849, win_df = h79_20k_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/Balsfjord/HapDist_v2.0.2/Reorder/", snp_col = "blue", dend_weight_vec = weight_vec)
pac_rec_win_HapDist <- plot_top_regions_3(top_pos = pac_rec_win_hits@from, win_df = h79_20k_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/Balsfjord/HapDist_v2.0.2/Rec_win_div/", snp_col = "blue", dend_weight_vec = weight_vec)

atl_div_vec_pac <- numeric()
HWS_div_vec_pac <- numeric()
dxy_vec_pac <- numeric()
for(reg_loc in names(pac_rec_win_HapDist$dist_mat)){
  atl_div_vec_pac[reg_loc] <- mean(pac_rec_win_HapDist$dist_mat[[reg_loc]][Atl_haps, Atl_haps], na.rm = T)/2e4
  HWS_div_vec_pac[reg_loc] <- mean(pac_rec_win_HapDist$dist_mat[[reg_loc]][HWS_haps, HWS_haps], na.rm = T)/2e4
  dxy_vec_pac[reg_loc] <- mean(pac_rec_win_HapDist$dist_mat[[reg_loc]][Atl_haps, HWS_haps], na.rm = T)/2e4
}

bg_win_HapDist <- plot_top_regions_3(top_pos = sample(dim(h79_20k_w)[1], 1000), win_df = h79_20k_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/Balsfjord/HapDist_v2.0.2/Rec_win_div/background/", snp_col = "blue", dend_weight_vec = weight_vec)
atl_div_vec_bg <- numeric()
HWS_div_vec_bg <- numeric()
dxy_vec_bg <- numeric()
for(reg_loc in names(bg_win_HapDist$dist_mat)){
  win_size <- (as.numeric(sub(".+_([0-9]+)_to_([0-9]+)_kb", "\\2", reg_loc)) - as.numeric(sub(".+_([0-9]+)_to_([0-9]+)_kb", "\\1", reg_loc))) * 1e3 # Adjusting for partial windows
  atl_div_vec_bg[reg_loc] <- mean(bg_win_HapDist$dist_mat[[reg_loc]][Atl_haps, Atl_haps], na.rm = T)/win_size
  HWS_div_vec_bg[reg_loc] <- mean(bg_win_HapDist$dist_mat[[reg_loc]][HWS_haps, HWS_haps], na.rm = T)/win_size
  dxy_vec_bg[reg_loc] <- mean(bg_win_HapDist$dist_mat[[reg_loc]][Atl_haps, HWS_haps], na.rm = T)/win_size
}

box_plot_df <- data.frame(win_type = rep(c("atl","pac", "bg"), times = c(length(atl_div_vec), length(atl_div_vec_pac), length(atl_div_vec_bg))), atl_div = c(atl_div_vec, atl_div_vec_pac, atl_div_vec_bg), HWS_div = c(HWS_div_vec, HWS_div_vec_pac, HWS_div_vec_bg), dxy = c(dxy_vec, dxy_vec_pac, dxy_vec_bg))
box_plot_df$win_type <- factor(box_plot_df$win_type, levels = c("atl", "pac", "bg"))
rownames(box_plot_df) <- c(paste0(names(atl_div_vec), "_atl"), paste0(names(atl_div_vec_pac),"_pac"), paste0(names(atl_div_vec_bg), "_bg"))


#Recombination rate in recurring windows
#Ch_v2_recombination_profile from elsewhere
Ch_v2_recomb_prof_GR <- GRanges(seqnames = paste0("chr", Ch_v2_recombination_profile$CHR), ranges = IRanges(star = Ch_v2_recombination_profile$BIN_START, width = 1e5))
Ch_v2_recomb_prof_GR$RR <- Ch_v2_recombination_profile$RR
atl_rec_v_recomb <- findOverlaps(Ch_v2_recomb_prof_GR, recurring_intro_GR)
pac_rec_v_recomb <- findOverlaps(Ch_v2_recomb_prof_GR, recurring_intro_pac_GR)
Ch_v2_recomb_prof_GR$recurring <- "none"
Ch_v2_recomb_prof_GR$recurring[atl_rec_v_recomb@from] <- "atl"
Ch_v2_recomb_prof_GR$recurring[pac_rec_v_recomb@from] <- "pac"
Ch_v2_recomb_prof_GR$recurring <- factor(Ch_v2_recomb_prof_GR$recurring, levels = c("atl", "pac", "none"))
aggregate(Ch_v2_recomb_prof_GR$RR ~ Ch_v2_recomb_prof_GR$recurring, FUN = "summary")

#Summary statistic boxplots
pdf(file = "~/Projects/Herring/doc/Balsfjord/Recurring_intro_region_boxplots.pdf")
boxplot(box_plot_df$atl_div ~ box_plot_df$win_type, main = "Nucleotide diversity among Atlantic haplotypes", xlab = "Region type", ylab = "Pi")
boxplot(box_plot_df$HWS_div ~ box_plot_df$win_type, main = "Nucleotide diversity among Arctic/Japanese haplotypes", xlab = "Region type", ylab = "Pi")
boxplot(box_plot_df$dxy ~ box_plot_df$win_type, main = "Dxy between Atlantic- vs Arctic/Japanese haplotypes", xlab = "Region type", ylab = "Dxy")
boxplot(Ch_v2_recomb_prof_GR$RR ~ Ch_v2_recomb_prof_GR$recurring, main = "Recombination rate", xlab = "Region type", ylab = "RR")
dev.off()

#tests
t.test(box_plot_df$atl_div ~ box_plot_df$win_type, subset =box_plot_df$win_type != "pac" )
t.test(box_plot_df$atl_div ~ box_plot_df$win_type, subset =box_plot_df$win_type != "pac" )$p.value
t.test(box_plot_df$atl_div ~ box_plot_df$win_type, subset =box_plot_df$win_type != "atl" )$p.value

t.test(box_plot_df$HWS_div ~ box_plot_df$win_type, subset =box_plot_df$win_type != "pac" )
t.test(box_plot_df$HWS_div ~ box_plot_df$win_type, subset =box_plot_df$win_type != "atl" )

t.test(box_plot_df$dxy ~ box_plot_df$win_type, subset =box_plot_df$win_type != "pac" )
t.test(box_plot_df$dxy ~ box_plot_df$win_type, subset =box_plot_df$win_type != "atl" )

t.test(Ch_v2_recomb_prof_GR$RR ~ Ch_v2_recomb_prof_GR$recurring, subset = Ch_v2_recomb_prof_GR$recurring != "pac" )
t.test(Ch_v2_recomb_prof_GR$RR ~ Ch_v2_recomb_prof_GR$recurring, subset = Ch_v2_recomb_prof_GR$recurring != "atl" )



exp_pac_count <- mean(hws6_v2.0.2_20k_intro$hap$ref_2_mean)
exp_atl_count <- mean(hws6_v2.0.2_20k_intro$hap$ref_1_mean)
mean_diverese_SNP_count <- exp_pac_count + exp_atl_count #Average number of variant SNPs (between references) per window
#Hyopthetical ratio = 0.18; 15.5 % SNPs diverges from Atl, 85.5% from Pac
hyp_count_atl <- mean_diverese_SNP_count * 0.155
hyp_count_pac <- mean_diverese_SNP_count * 0.855
fisher.test(x = matrix(round(c(hyp_count_atl, hyp_count_pac, exp_atl_count, exp_pac_count)), nrow = 2, byrow = F))


#Looking at decay around recurring signals - Chr4 19 Mb
which(h79_w$scaffold == "chr4" & h79_w$start== 19e6+1) #1179
plot_top_regions_2(top_pos = c(1179, 1180), win_df = h79_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/Balsfjord/HapDist_v2.0.2/Recurring_windows/", snp_col = "blue")

#Second sub-peak
plot_top_regions_2(top_pos = c(1184,1185), win_df = h79_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/Balsfjord/HapDist_v2.0.2/Recurring_windows/", snp_col = "blue")
#The pattern is not the same, likely another origin
#Overlapping gene
glu_met_r <- cluhar_v2.0.2_gtf[cluhar_v2.0.2_gtf$gene_id == "ENSCHAG00000006184" & cluhar_v2.0.2_gtf$type == "gene"]

#pdf(file = "~/Projects/Herring/doc/Balsfjord/HWS6_chr4_19Mb.pdf", width = 10)
pdf(file = "~/Projects/Herring/doc/Balsfjord/HWS6_chr4_19Mb_reorder.pdf", width = 10)
plot(x = 1, y = 1, type = "n", main = "Introgression \"Coverage\" & haplotypes", xlab = "", ylab = "Number of haplotypes", xlim = c(18e6, 20e6), ylim = c(-9,9), axes =F )
axis(2, at = c(0,4,8))
axis(1)
segments(x0 = atl_intro_cov_20k$start[atl_intro_cov_20k$group_name == "chr4"], x1 = atl_intro_cov_20k$end[atl_intro_cov_20k$group_name == "chr4"], y0 = atl_intro_cov_20k$cov[atl_intro_cov_20k$group_name == "chr4"], lwd = 4, col = "grey70")
segments(x0 = glu_met_r@ranges@start, x1 = glu_met_r@ranges@start + glu_met_r@ranges@width, y0 = 0, lwd = 4, col = "darkorchid")
text(x = glu_met_r@ranges@start + glu_met_r@ranges@width/2, y = 0, labels = "GRM4", cex = 0.8, pos = 3)

#points(x = Ch_v2_recombination_profile$BIN_START[Ch_v2_recombination_profile$CHR == 4], y = (Ch_v2_recombination_profile$RR[Ch_v2_recombination_profile$CHR == 4]/4), col = "darkorchid")
#hap_col_vec <- rep(c("firebrick1", "darkorchid1"), times = 4)
#hap_order_vec <- c(1,4,8,6,3,5,2,7) #Haplotypes grouped by likely IBD, and dcreasing decay
hap_order_vec <- c(1:8) #Haplotypes grouped by individual
hap_col_vec <- c("black","red", "blue")[c(1,2,1,1,1,1,3,1)]

for(hap in hap_order_vec){
  target_segments <- which(intro_20k_lists$atl_red[[hap]]@seqnames=="chr4")
  x0_vec <- intro_20k_lists$atl_red[[hap]]@ranges@start[target_segments]
  x1_vec <- intro_20k_lists$atl_red[[hap]]@ranges@width[target_segments] +  x0_vec
  y_pos <- -which(hap_order_vec == hap)
  segments(x0 = x0_vec, x1 =x1_vec, y0 = y_pos, lwd = 4, col = hap_col_vec[hap])
  text(x = 18e6, y = y_pos, labels = names(intro_20k_lists$atl_red)[hap], pos = 4, cex = 0.8, col = hap_col_vec[hap])
}

#Approximating the orignal haplotype
atl_intro_cov_20k[atl_intro_cov_20k$group_name == "chr4" & atl_intro_cov_20k$end > 18e6 & atl_intro_cov_20k$start < 20e6,]
abline(v = c(18680000, 19420001))
dev.off()

Ch_v2_linkage_map <- read.table("~/Projects/Herring/data/HiC_assemblies/Version_2_release/20190412_v2.0.2_linkage_map.txt", header = T, stringsAsFactors=F)
r_avg <- ((59.179 - 52.863)/(23279630-9351242))/100 #Average recombination rate across the region in M/base (from map)
r_reg <- r_avg*(19420001 -18680000) #Chance of recombinaion in the region per generation ~ 0.5%; matches 6cM/12Mb
1/r_reg #Generations between recombinations, on average (~300)





###
pdf("~/Projects/Herring/doc/Balsfjord/age_estimate/GW_size_rank.pdf")
atl_widths <- unlist(intro_20k_lists$atl_red)@ranges@width
pac_widths <- unlist(intro_20k_lists$pac_red)@ranges@width
y_vec <- jitter(c(atl_widths,pac_widths), amount = 5000)
col_vec <- rep(c("red", "blue"), times = c(length(atl_widths), length(pac_widths)))
order_vec <- order(y_vec, decreasing = T)
plot(y = y_vec[order_vec] , x = 1:length(y_vec), pch = 20, col = col_vec[order_vec], cex = 2, xlab = "Introgression Size Index", ylab = "Region Size")
text(x = 5e3, y = 1.5e6, labels = paste0("P(wilcox) = ", signif(wilcox.test(atl_widths, pac_widths)$p.value, digits = 2)))
legend(x= "topright", legend = c("Pacific introgressions", "Atlantic introgressions"), pch = 20, col = c("blue", "red"), pt.cex = 2)
dev.off()
###

#Examining Atl v Pac Fst profiles
PvA_low_Fst <- read.table("~/Projects/Herring/data/Balsfjord/low_fst_atl_pac.bed", sep = "\t")
names(PvA_low_Fst) <- c("seqnames", "start", "end", "width")
PvA_low_Fst_GR <- GRanges(PvA_low_Fst)
findOverlaps(recurring_intro_pac_GR, PvA_low_Fst_GR)
findOverlaps(recurring_intro_GR, PvA_low_Fst_GR)

PvA_Fst <- read.csv("~/Projects/Herring/data/Balsfjord/atl_pac_bal_combined_diversity.csv", sep = ",")
PvA_Fst_alternatives <- read.csv("~/Projects/Herring/data/Balsfjord/df.alt.pac.bal.pbs.csv", sep = ",")
PvA_GR <- GRanges(PvA_Fst$chr, ranges = IRanges(start = PvA_Fst$BIN_START, end = PvA_Fst$BIN_END))
high_Fst <- which(PvA_Fst$AtlPac_Fst > 0.6)
PvA_high_Fst_GR <- GRanges(PvA_Fst$chr, ranges = IRanges(start = PvA_Fst$BIN_START, end = PvA_Fst$BIN_END))[high_Fst]
Atl_v_high_Fst <- findOverlaps(recurring_intro_GR, PvA_high_Fst_GR)
PvA_high_Fst_df <- as.data.frame(PvA_high_Fst_GR, stringsAsFactors = F)
PvA_high_Fst_df[,"global_start"] <- PvA_high_Fst_df[,"start"] + Ch_v2.0.2_sizes[match(PvA_high_Fst_df[, "seqnames"], Ch_v2.0.2_sizes[,"name"]), "offset"]
PvA_high_Fst_df[,"global_end"] <-  PvA_high_Fst_df[,"end"] + Ch_v2.0.2_sizes[match(PvA_high_Fst_df[, "seqnames"], Ch_v2.0.2_sizes[,"name"]), "offset"]
PvA_high_Fst_df[,"atl_rec_overlap"] <- F
PvA_high_Fst_df[Atl_v_high_Fst@to,"atl_rec_overlap"] <- T



#Whole genome
plot(x = atl_intro_cov_20k$global_start, y = atl_intro_cov_20k$cov, type = "n", main = "Introgression \"Coverage\"", xlab = "", ylab = "Number of haplotypes", ylim = c(0,10)) 
segments(x0 = atl_intro_cov_20k$global_start, x1 = atl_intro_cov_20k$global_end, y0 = atl_intro_cov_20k$cov, lwd = 4, col = atl_intro_cov_20k$col)
target_blocks <- !PvA_high_Fst_df$atl_rec_overlap
segments(x0 = PvA_high_Fst_df$global_start[target_blocks], x1 = PvA_high_Fst_df$global_end[target_blocks], y0 = 9, lwd = 4, col = "blue")
target_blocks <- PvA_high_Fst_df$atl_rec_overlap
segments(x0 = PvA_high_Fst_df$global_start[target_blocks ], x1 = PvA_high_Fst_df$global_end[target_blocks], y0 = 9.5, lwd = 4, col = "red")

#
# target_chr <- "chr23"
# plot(x = atl_intro_cov_20k$start[atl_intro_cov_20k$group_name == target_chr], y = atl_intro_cov_20k$cov[atl_intro_cov_20k$group_name == target_chr], type = "n", main = paste("Introgression \"Coverage\";", target_chr ), xlab = "", ylab = "Number of haplotypes", ylim = c(0,10)) 
# segments(x0 = atl_intro_cov_20k$start[atl_intro_cov_20k$group_name == target_chr], x1 = atl_intro_cov_20k$end[atl_intro_cov_20k$group_name == target_chr], y0 = atl_intro_cov_20k$cov[atl_intro_cov_20k$group_name == target_chr], lwd = 4, col = atl_intro_cov_20k$col[atl_intro_cov_20k$group_name == target_chr])
# target_blocks <- PvA_high_Fst_df$seqnames == target_chr & !PvA_high_Fst_df$atl_rec_overlap
# segments(x0 = PvA_high_Fst_df$start[target_blocks ], x1 = PvA_high_Fst_df$end[target_blocks], y0 = 9, lwd = 4, col = "blue")
# target_blocks <- PvA_high_Fst_df$seqnames == target_chr & PvA_high_Fst_df$atl_rec_overlap
# segments(x0 = PvA_high_Fst_df$start[target_blocks ], x1 = PvA_high_Fst_df$end[target_blocks], y0 = 9.5, lwd = 4, col = "red")

#Chromosomes with overlap
fst_v_intro_plot("chr4", cov_df = atl_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_GR)
fst_v_intro_plot("chr6", cov_df = atl_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_GR)
fst_v_intro_plot("chr8", cov_df = atl_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_GR)
fst_v_intro_plot("chr13", cov_df = atl_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_GR)
fst_v_intro_plot("chr19", cov_df = atl_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_GR)
fst_v_intro_plot("chr23", cov_df = atl_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_GR)

#Chromosomes without overlap
fst_v_intro_plot("chr17", cov_df = atl_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_GR)
fst_v_intro_plot("chr18", cov_df = atl_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_GR)
fst_v_intro_plot("chr21", cov_df = atl_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_GR)

#Pacific version
unique(recurring_intro_pac_df$seqnames)
findOverlaps(recurring_intro_pac_GR, PvA_high_Fst_GR)
#Chromosomes without overlap
fst_v_intro_plot("chr1", cov_df = pac_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_pac_GR)
fst_v_intro_plot("chr3", cov_df = pac_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_pac_GR)
fst_v_intro_plot("chr4", cov_df = pac_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_pac_GR)
fst_v_intro_plot("chr6", cov_df = pac_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_pac_GR)
fst_v_intro_plot("chr7", cov_df = pac_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_pac_GR)
fst_v_intro_plot("chr8", cov_df = pac_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_pac_GR)
fst_v_intro_plot("chr10", cov_df = pac_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_pac_GR)
fst_v_intro_plot("chr12", cov_df = pac_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_pac_GR)
fst_v_intro_plot("chr15", cov_df = pac_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_pac_GR)
fst_v_intro_plot("chr18", cov_df = pac_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR = recurring_intro_pac_GR)

#Combined, all chromosomes
pdf("~/Projects/Herring/doc/Balsfjord/Fst_v_introgression.pdf", width = 12, height = 5)
for(chr in paste0("chr", 1:26)){
  fst_v_intro_plot2(target_chr = chr, cov_df1 = atl_intro_cov_20k, cov_df2 = pac_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR1 = recurring_intro_GR , target_reg_GR2 = recurring_intro_pac_GR)
}
dev.off()

#Less stringent
intro_7up_pac_df <- pac_intro_cov_20k[pac_intro_cov_20k$cov >= 7,c(2:4)]
names(intro_7up_pac_df )[1] <- "seqnames"
intro_7up_pac_GR <- GRanges(intro_7up_pac_df)
intro_7up_atl_df <- atl_intro_cov_20k[atl_intro_cov_20k$cov >= 7,c(2:4)]
names(intro_7up_atl_df)[1] <- "seqnames"
intro_7up_atl_GR <- GRanges(intro_7up_atl_df)

pdf("~/Projects/Herring/doc/Balsfjord/Fst_v_introgression_7up.pdf", width = 12, height = 5)
for(chr in paste0("chr", 1:26)){
  fst_v_intro_plot2(target_chr = chr, cov_df1 = atl_intro_cov_20k, cov_df2 = pac_intro_cov_20k, full_Fst = PvA_Fst, target_reg_GR1 = intro_7up_atl_GR , target_reg_GR2 = intro_7up_pac_GR)
}
dev.off()

#Correlation
pac_intro_cov_20k_GR <- GRanges(seqnames = pac_intro_cov_20k$group_name, ranges = IRanges(start = pac_intro_cov_20k$start, end = pac_intro_cov_20k$end))
PvA_v_Pac_cov <- findOverlaps(PvA_GR, pac_intro_cov_20k_GR, select = "arbitrary")
PvA_Fst$pac_intro_cov <- pac_intro_cov_20k$cov[PvA_v_Pac_cov]
atl_intro_cov_20k_GR <- GRanges(seqnames = atl_intro_cov_20k$group_name, ranges = IRanges(start = atl_intro_cov_20k$start, end = atl_intro_cov_20k$end))
PvA_v_Atl_cov <- findOverlaps(PvA_GR, atl_intro_cov_20k_GR, select = "arbitrary")
PvA_Fst$atl_intro_cov <- atl_intro_cov_20k$cov[PvA_v_Atl_cov]
png("~/Projects/Herring/doc/Balsfjord/HWS6v1_Fst_t_test.png", height = 1000, width = 1500)
plot(y = PvA_Fst$AtlPac_Fst, x = PvA_Fst$atl_intro_cov, pch = 20, col = "darkorange", ylab = "Fst", xlab = "Number of introgressed haplotypes", xlim = c(0,9)) 
points(y = PvA_Fst$AtlPac_Fst, x = PvA_Fst$pac_intro_cov+0.1, pch = 20, col = "olivedrab")
reg_coef <- lm(PvA_Fst$AtlPac_Fst[PvA_Fst$pac_intro_cov >= 1]~PvA_Fst$pac_intro_cov[PvA_Fst$pac_intro_cov >= 1])$coef
lines(x = c(0.8,8.5), y = c(0.8,8.5)*reg_coef[2]+reg_coef[1],col = "olivedrab",lwd = 2)
reg_coef <- lm(PvA_Fst$AtlPac_Fst[PvA_Fst$atl_intro_cov >= 1]~PvA_Fst$atl_intro_cov[PvA_Fst$atl_intro_cov >= 1])$coef
lines(x = c(0.8,8.5), y = c(0.8,8.5)*reg_coef[2]+reg_coef[1], col = "darkorange", lwd = 2)
for(cov_val in 1:8){
  t_test <- t.test(x= PvA_Fst$AtlPac_Fst[PvA_Fst$atl_intro_cov == cov_val], y = PvA_Fst$AtlPac_Fst[PvA_Fst$pac_intro_cov == cov_val])
  text(x = cov_val, y = 0.95, labels = paste("p =", signif(t_test$p.value, digits = 2)), cex = 0.6)
  text(x = cov_val, y = 1, labels = paste("est =", signif(diff(t_test$estimate), digits = 2)), cex = 0.6)
}
abline(h = 0.6, col = "grey40")
dev.off()

boxplot(PvA_Fst$AtlPac_Fst~PvA_Fst$atl_intro_cov, at = (0:8)-0.15, boxwex = 0.2, xlab = "Introgression haplotypes", ylab = "Fst (Pac v Atl)", axes = F, col = "darkorange", notch = T)
boxplot(PvA_Fst$AtlPac_Fst~PvA_Fst$pac_intro_cov, at = (0:8) + 0.15, boxwex = 0.2, add = T, axes = F, col = "olivedrab", notch = T) 
axis(1, at = 0:8)
axis(2)
#Percentage of the bins that have coverage 7 or 8
sum(PvA_Fst$atl_intro_cov >= 7, na.rm = T)/dim(PvA_Fst)[1]*100 
#[1] 0.3131698
#Percentage of the high (>= 0.6) Fst-bins that also have coverage 7 or 8 
sum(PvA_Fst$atl_intro_cov >= 7 & PvA_Fst$AtlPac_Fst >= 0.6, na.rm = T)/sum(PvA_Fst$AtlPac_Fst >= 0.6, na.rm = T)*100
#[1] 10.88957
#M-value
log2(10.88957/0.3131698)
intro_Mval_df <- data.frame(cov = 0:8, stringsAsFactors = F)
for(i in 0:8){
  intro_Mval_df[intro_Mval_df$cov == i,"Atl_size_frac"] <- sum(PvA_Fst$atl_intro_cov == i, na.rm = T)/dim(PvA_Fst)[1]
  intro_Mval_df[intro_Mval_df$cov == i,"Atl_high_Fst_frac"] <-sum(PvA_Fst$atl_intro_cov == i & PvA_Fst$AtlPac_Fst >= 0.6, na.rm = T)/sum(PvA_Fst$AtlPac_Fst >= 0.6, na.rm = T)
  intro_Mval_df[intro_Mval_df$cov == i,"Pac_size_frac"] <- sum(PvA_Fst$pac_intro_cov == i, na.rm = T)/dim(PvA_Fst)[1]
  intro_Mval_df[intro_Mval_df$cov == i,"Pac_high_Fst_frac"] <-sum(PvA_Fst$pac_intro_cov == i & PvA_Fst$AtlPac_Fst >= 0.6, na.rm = T)/sum(PvA_Fst$AtlPac_Fst >= 0.6, na.rm = T)
}
intro_Mval_df[,"Atl_ratio"] <- intro_Mval_df$Atl_high_Fst_frac/intro_Mval_df$Atl_size_frac
intro_Mval_df[,"Atl_M_val"] <- log2(intro_Mval_df[,"Atl_ratio"])
intro_Mval_df[,"Pac_ratio"] <- intro_Mval_df$Pac_high_Fst_frac/intro_Mval_df$Pac_size_frac
intro_Mval_df[,"Pac_M_val"] <- log2(intro_Mval_df[,"Pac_ratio"])

#Binomial test vs the genomic average
p_vec <- numeric()
for(i in 1:8){
  n_tries <- sum(PvA_Fst$atl_intro_cov == i, na.rm = T)
  p_succ <- sum(PvA_Fst$AtlPac_Fst >= 0.6, na.rm = T)/sum(!is.na(PvA_Fst$AtlPac_Fst))
  n_obs <- sum(PvA_Fst$atl_intro_cov == i & PvA_Fst$AtlPac_Fst >= 0.6, na.rm = T)
  p_vec[i] <- binom.test(x = n_obs, n = n_tries, p = p_succ)$p.value
}

pdf("~/Projects/Herring/doc/Balsfjord/Fst_v_introgression_M_val.pdf")
boxplot(PvA_Fst$AtlPac_Fst~PvA_Fst$atl_intro_cov, at = (0:8)-0.15, boxwex = 0.2, xlab = "Introgression haplotype count", ylab = "", axes = F, col = "darkorange", notch = T, ylim = c(0,1.7))
boxplot(PvA_Fst$AtlPac_Fst~PvA_Fst$pac_intro_cov, at = (0:8) + 0.15, boxwex = 0.2, add = T, axes = F, col = "olivedrab", notch = T) 
axis(1, at = 0:8)
axis(2, at = c(0,0.5,1))
lines(x = 0:8, y = (intro_Mval_df$Atl_M_val)/20+1.2, lwd = 2, col = "darkorange") 
lines(x = 0:8, y = (intro_Mval_df$Pac_M_val)/20+1.2, lwd = 2, col = "olivedrab") 
axis(2, at = c(1.2,1.3,1.4,1.5,1.6)-0.1, labels = c(-2,0,2,4,6), las = 2)
lines(x = c(-1, 8), y = c(1.2, 1.2))
mtext(text = c("Fst (Pac v Atl)", "M value (Fst >0.6)"), 2, at = c(0.5, 1.3), line = 2)
legend(x = "topleft", legend = c("Atlantic introgressions", "Pacific introgressions"), fill = c("darkorange","olivedrab"))
dev.off()

#Same, for Pacific "introgressions"
sum(PvA_Fst$pac_intro_cov >= 7, na.rm = T)/dim(PvA_Fst)[1]*100
#[1] 0.2079242
#Percentage of the high (>= 0.6) Fst-bins that also have coverage 7 or 8 
sum(PvA_Fst$pac_intro_cov >= 7 & PvA_Fst$AtlPac_Fst >= 0.6, na.rm = T)/sum(PvA_Fst$AtlPac_Fst >= 0.6, na.rm = T)*100
#[1] 0

plot(x = 0:8, y = intro_Mval_df$Atl_M_val, lwd = 3, col = "darkorange", xlab = "Introgression haplotype count", ylab = "Representation of >0.6 Fst windows (M value)", type = "l") 
lines(x = 0:8, y = intro_Mval_df$Pac_M_val, lwd = 3, col = "olivedrab") 




#Rough counts of callable fraction of the genome
hist(PvA_Fst$AtlPac_Fst[PvA_Fst$atl_intro_cov >= 1 | PvA_Fst$pac_intro_cov >= 1])
#Fst of >0.1 seems sufficient to call an introgression
#callable_size <- sum(PvA_Fst$AtlPac_Fst >= 0.10, na.rm = T) * 5000
callable_size <- sum(PvA_Fst$AtlPac_Fst >= 0.15, na.rm = T) * 5000
hap_intro_lengths <- unlist(lapply(intro_20k_lists$atl_red, function(x){sum(x@ranges@width)}))
hap_intro_frac <- max(hap_intro_lengths)/callable_size
hap_intro_lengths_pac <- unlist(lapply(intro_20k_lists$pac_red, function(x){sum(x@ranges@width)}))
hap_intro_frac_pac <- max(hap_intro_lengths_pac)/callable_size

#In bases
exp_8_overlap <- hap_intro_frac^8*callable_size
exp_7_overlap <- (choose(8,7)*hap_intro_frac^7*(1-hap_intro_frac)^1 + choose(8,8)*hap_intro_frac^8*(1-hap_intro_frac)^0)*callable_size
obs_8_overlap <- sum(atl_intro_cov_20k$width[atl_intro_cov_20k$cov >= 8])
obs_7_overlap <- sum(atl_intro_cov_20k$width[atl_intro_cov_20k$cov >= 7])
#M-value
log2(obs_8_overlap/exp_8_overlap)
log2(obs_7_overlap/exp_7_overlap)

exp_overlap_vec <- numeric(9)
obs_overlap_vec <- numeric(9)
exp_overlap_vec_pac <- numeric(9)
obs_overlap_vec_pac <- numeric(9)
for(i in 0:8){
  #exp_overlap_vec[i] <- dbinom(i,8, hap_intro_frac)*callable_size
  #exp_overlap_vec[i+1] <- dbinom(i,8, mean(hap_intro_lengths/callable_size))*callable_size
  exp_overlap_vec[i+1] <- dbinom(i,8, max(hap_intro_lengths/callable_size))*callable_size
  obs_overlap_vec[i+1] <- sum(atl_intro_cov_20k$width[atl_intro_cov_20k$cov == i])
  #exp_overlap_vec_pac[i+1] <- dbinom(i,8, mean(hap_intro_lengths_pac/callable_size))*callable_size
  exp_overlap_vec_pac[i+1] <- dbinom(i,8, max(hap_intro_lengths_pac/callable_size))*callable_size
  obs_overlap_vec_pac[i+1] <- sum(pac_intro_cov_20k$width[pac_intro_cov_20k$cov == i])
} 
obs_overlap_vec[1] <- callable_size - sum(obs_overlap_vec[2:9])
obs_overlap_vec_pac[1] <- callable_size - sum(obs_overlap_vec_pac[2:9])

pdf("~/Projects/Herring/doc/Balsfjord/intro_M_values.pdf")
plot(x = 0:8, y = log2(obs_overlap_vec/exp_overlap_vec), xlab = "Introgression coverage", pch = 16, ylab = "M value (Observed over Expected)", ylim = c(-3,25), cex = 1.5)
points(x = 0:8, y = log2(obs_overlap_vec_pac/exp_overlap_vec_pac), pch = 16, col  = "red", cex = 1.5)
legend(x = "topleft", legend = c("Atlantic", "Pacific"), pch = 16, col = c("black", "red"), cex = 1.5)
dev.off()

## Whole-genome trees of relevant samples - TBD
pool_freq <- read.table("~/Projects/Herring/data/v2.0.2_genotypes/60.Neff.freq.gz", stringsAsFactors = F, sep = "\t", header = T)
save(pool_freq, file = "~/Projects/Herring/data/v2.0.2_genotypes/pool_freq_60Neff.RData")
f3 <- Vectorize(FUN = function(X,Y, data_mat){sum(abs(data_mat[,X] - data_mat[,Y]), na.rm = T)/sum(!is.na(data_mat[,X]) & !is.na(data_mat[,Y]))}, vectorize.args = c("X", "Y"))
tree_SNPs <- sample.int(n = dim(pool_freq)[1], size = 1e5, replace = F)
pool_raw_dist <- outer(X = 3:dim(pool_freq)[2], Y = 3:dim(pool_freq)[2], FUN = "f3", data_mat = pool_freq[tree_SNPs,])
rownames(pool_raw_dist) <- names(pool_freq)[3:dim(pool_freq)[2]]
colnames(pool_raw_dist) <- names(pool_freq)[3:dim(pool_freq)[2]]
save(pool_raw_dist, file = "~/Projects/Herring/data/v2.0.2_genotypes/pool_freq_60Neff_dist.RData")

#Trees by APE
require(Biostrings)
require(ape)
#HapDNA <- DNAStringSet(hap_vec)
#HapDNA_bin <- as.DNAbin(HapDNA)
pool_dist <- as.dist(pool_raw_dist)
tree <- bionj(pool_dist)
tree$tip.label <- paste0("___", colnames(pool_raw_dist), "___")

#Original color scheme
#tip_col_vec <- rep("darkorange1", length(tree$tip.label))
#tip_col_vec[grep("Baltic", tree$tip.label)] <- "darkorange4"
#tip_col_vec[grep("Landvik", tree$tip.label)] <- "darkorange1"

#tip_col_vec[grep("HWS", tree$tip.label)] <- "olivedrab4"
#tip_col_vec[grep("PB8", tree$tip.label)] <- "olivedrab2"
#tip_col_vec[grep("HWS6", tree$tip.label)] <- "darkolivegreen"

#Different color scheme
tip_col_vec <- rep("darkorange1", length(tree$tip.label))
tip_col_vec[grep("Baltic", tree$tip.label)] <- "darkorange4"
tip_col_vec[grep("Landvik", tree$tip.label)] <- "darkorange1"

tip_col_vec[grep("HWS", tree$tip.label)] <- "olivedrab4"
tip_col_vec[grep("HWS1", tree$tip.label)] <- "steelblue2"
tip_col_vec[grep("PB8", tree$tip.label)] <- "olivedrab3"
tip_col_vec[grep("HWS6", tree$tip.label)] <- "darkorchid"

plot_tree <- tree 
plot_tree$tip.label <- rep("\U2022", length(tree$tip.label))


png(file = "~/Projects/Herring/doc/Balsfjord/pool_tree_dots.png", height = 1000, width = 1000)
plot.phylo(plot_tree, type="unrooted", lab4ut = "axial", tip.color = tip_col_vec, cex = 3, edge.width = 2.5)
dev.off()

#tip_color_vec <- rep("black", times = length(sample_list))
#tip_color_vec[grep("HWS", sample_list)] <- "darkorchid4"
#tip_color_vec[grep("Pacific", sample_list)] <- "blue"
#tip_color_vec[grep("North-Sea", sample_list)] <- "red"
####


#Sample locations
Pool_loc <- read.csv("~/Projects/Herring/data/genomic_data/NorthSea/GENSINC_pool.csv", stringsAsFactors = F)
pool_plot_vec <- c(29, 30, 31, 32, 47)

require(maptools)
require(raster)
require(rgdal)

wrld_map <- readShapePoly("~/Projects/Herring/data/SNP_chip/SEC16B/10m-admin-0-countries/10m_admin_0_countries.shp")
wrld_map <- readOGR("~/Projects/Herring/data/SNP_chip/SEC16B/10m-admin-0-countries/10m_admin_0_countries.shp")

#cr <- colorRamp(c("blue", "red"))
pdf(file = "~/Projects/Herring/doc/Balsfjord/sample_map.pdf", height = 5, width = 10)
plot(x = 15, y = 57, type = "n", xlim = c(10,140), ylim = c(35,72), xlab = "Longitude", ylab = "Latitude", main = "")
plot(wrld_map, add = T, col = "grey90", border = "grey70")
points(x= Pool_loc$lon[pool_plot_vec], y = Pool_loc$lat[pool_plot_vec], col = "olivedrab", cex = 1.5, pch = 16)
#text(x= pool_lon, y = pool_lat, labels= paste("'A' freq: ", round(pool_means[pool_reorder], digits= 2), sep  =""),pos = 3, cex = 0.5)
#text(x= pool_lon, y = pool_lat, labels= pool_names[pool_reorder], pos = 1, cex = 0.5)
dev.off()

#Fst of HWS6 vs other HWS
HWS6v1_fst <- read.table("~/Projects/Herring/data/Balsfjord/fst/HWS1_v_HWS6_50k.windowed.weir.fst", stringsAsFactors = F, header = T, sep = "\t")
HWS6v1_fst[,"col"] <- c("grey30", "grey70")[(as.integer(sub("chr|unplaced_scaffold", "", HWS6v1_fst$CHR)) %% 2) + 1]
HWS6v1_fst[,"cumulative_pos"] <- HWS6v1_fst$BIN_START + Ch_v2.0.2_sizes[match(HWS6v1_fst$CHR, Ch_v2.0.2_sizes$name),"offset"]
hist(HWS6v1_fst$N_VARIANTS)
win_filter <- HWS6v1_fst$N_VARIANTS > 100 & HWS6v1_fst$N_VARIANTS < 1100
plot(x = HWS6v1_fst$cumulative_pos[win_filter], y = HWS6v1_fst$WEIGHTED_FST[win_filter], col = HWS6v1_fst$col[win_filter], pch = 16, cex = 0.5)
plot(x = HWS6v1_fst$BIN_START[win_filter & HWS6v1_fst$CHROM == "chr6"], y = HWS6v1_fst$WEIGHTED_FST[win_filter & HWS6v1_fst$CHROM == "chr6"], col = HWS6v1_fst$col[win_filter & HWS6v1_fst$CHROM == "chr6"], pch = 16, cex = 0.5)

HWS6v1_daf <- pool_freq[,1:2]
HWS6v1_daf$HWS6v1_daf <- abs(pool_freq$HWS6_Balsfjord_Atlantic - pool_freq$HWS1_Japan_SeaOfJapan)
HWS6v1_daf$HWS6v1to5_daf <- abs(pool_freq$HWS6_Balsfjord_Atlantic - rowMeans(pool_freq[, grep("HWS[12345]", names(pool_freq))]))

HWS6v1_daf[,"col"] <- c("grey30", "grey70")[(as.integer(sub("chr|unplaced_scaffold", "", HWS6v1_daf$CHROM)) %% 2) + 1]
HWS6v1_daf[,"cumulative_pos"] <- HWS6v1_daf$POS + Ch_v2.0.2_sizes[match(HWS6v1_daf$CHR, Ch_v2.0.2_sizes$name),"offset"]
reliable_pos <- paste(HWS6v1_daf$CHROM, HWS6v1_daf$POS, sep  = "_") %in% paste(S_v_N$CHR, S_v_N$POS , sep  = "_")
HWS6v1_daf <- HWS6v1_daf[reliable_pos,]


png("~/Projects/Herring/doc/Balsfjord/HWS6v1to5_GW_daf.png", height = 700, width = 1500)
plot(x = HWS6v1_daf$cumulative_pos, y = HWS6v1_daf$HWS6v1to5_daf, col = HWS6v1_daf$col, pch = 16, cex = 0.2, main = "DAF")
pos_filter <- !is.na(HWS6v1_daf$HWS6v1to5_daf)
y_vec <- filter(x = HWS6v1_daf$HWS6v1to5_daf[pos_filter], filter = rep(1/100, 100))
points(x = HWS6v1_daf$cumulative_pos[pos_filter], y = y_vec, col = "red", pch = 16, cex = 0.4)
dev.off()

png("~/Projects/Herring/doc/Balsfjord/HWS6v1to5_chr6_daf.png", height = 700, width = 1500)
plot(x = HWS6v1_daf$POS[HWS6v1_daf$CHROM == "chr6"], y = HWS6v1_daf$HWS6v1to5_daf[HWS6v1_daf$CHROM == "chr6"], col = HWS6v1_daf$col[HWS6v1_daf$CHROM == "chr6"], pch = 16, cex = 0.2, main = "Chr 6; DAF")
pos_filter <- HWS6v1_daf$CHROM == "chr6" & !is.na(HWS6v1_daf$HWS6v1to5_daf)
y_vec <- filter(x = HWS6v1_daf$HWS6v1to5_daf[pos_filter], filter = rep(1/100, 100))
points(x = HWS6v1_daf$POS[pos_filter], y = y_vec, col = "red", pch = 16, cex = 0.4)
abline(v = c(14100001,14320000))
dev.off()

png("~/Projects/Herring/doc/Balsfjord/HWS6v1to5_chr4_daf.png", height = 700, width = 1500)
plot(x = HWS6v1_daf$POS[HWS6v1_daf$CHROM == "chr4"], y = HWS6v1_daf$HWS6v1to5_daf[HWS6v1_daf$CHROM == "chr4"], col = HWS6v1_daf$col[HWS6v1_daf$CHROM == "chr4"], pch = 16, cex = 0.2, main = "Chr 4; DAF")
pos_filter <- HWS6v1_daf$CHROM == "chr4" & !is.na(HWS6v1_daf$HWS6v1to5_daf)
y_vec <- filter(x = HWS6v1_daf$HWS6v1to5_daf[pos_filter], filter = rep(1/100, 100))
points(x = HWS6v1_daf$POS[pos_filter], y = y_vec, col = "red", pch = 16, cex = 0.4)
abline(v = c(19020001,19080000))
dev.off()

pdf("~/Projects/Herring/doc/Balsfjord/HWS6v1_Fst.pdf", width = 10)
plot(x = HWS6v1_fst$cumulative_pos[win_filter], y = HWS6v1_fst$WEIGHTED_FST[win_filter], col = HWS6v1_fst$col[win_filter], pch = 16, cex = 0.5, main = "Windowed Fst")
plot(x = HWS6v1_fst$BIN_START[win_filter & HWS6v1_fst$CHROM == "chr6"], y = HWS6v1_fst$WEIGHTED_FST[win_filter & HWS6v1_fst$CHROM == "chr6"], col = HWS6v1_fst$col[win_filter & HWS6v1_fst$CHROM == "chr6"], pch = 16, cex = 0.5, main = "Chr 6; Fst")
abline(v = c(14100001,14320000))
plot(x = HWS6v1_fst$BIN_START[win_filter & HWS6v1_fst$CHROM == "chr4"], y = HWS6v1_fst$WEIGHTED_FST[win_filter & HWS6v1_fst$CHROM == "chr4"], col = HWS6v1_fst$col[win_filter & HWS6v1_fst$CHROM == "chr4"], pch = 16, cex = 0.5, main = "Chr 4; Fst")
abline(v = c(19020001,19080000))
dev.off()

#Heatmaps of recurring regions
#balsfjord_pool_vec <- c(grep("Pacific|HWS", names(pool_freq)), grep("Pacific|HWS", names(pool_freq), invert = T))
balsfjord_pool_vec <- c(grep("Pacific|HWS", names(pool_freq), invert = T), grep("HWS", names(pool_freq)), grep("Pacific", names(pool_freq)))

balsfjord_pool_vec <- balsfjord_pool_vec[!(balsfjord_pool_vec %in% 1:2)]
#correcting miss-spelled names
plot_pool_names <- names(pool_freq)[balsfjord_pool_vec]
plot_pool_names <- sub("Pacific_Pacific", "Strait_of_Georgia_Pacific", plot_pool_names)
plot_pool_names <- sub("Galve", "Gävle", plot_pool_names)
plot_pool_names <- sub("HastKar", "Hästskär", plot_pool_names)
plot_pool_names <- sub("Traslovslage", "Träslövsläge", plot_pool_names)
plot_pool_names <- sub("LandvikS17_Norway_Baltic_Spring", "S17_Landvik_Atlantic_Spring", plot_pool_names)
plot_pool_names <- sub("Rugen", "Rügen", plot_pool_names)
plot_pool_names <- sub("Riga", "Gulf_of_Riga", plot_pool_names)
plot_pool_names <- sub("Lindas", "Lindås", plot_pool_names)
plot_pool_names <- sub("Ringkobing", "Ringkøbing", plot_pool_names)
plot_pool_names <- sub("NSSH", "Norway", plot_pool_names)
plot_pool_names <- sub("DalFB", "Fortune_Bay", plot_pool_names)
plot_pool_names <- sub("DalBoB", "Bonavista_Bay", plot_pool_names)
plot_pool_names <- sub("DalNsF", "Northumberland_Strait", plot_pool_names)
plot_pool_names <- sub("DalGeB", "German_Banks", plot_pool_names)
plot_pool_names <- sub("DalInB", "Inner_Baie_Des_Chaleurs", plot_pool_names)
plot_pool_names <- sub("DalNsS", "Northumberland_Strait", plot_pool_names)
plot_pool_names <- sub("HGS20_CapeWrath_Atlantic_Spring", "HGS20_Isle_of_Skye_Atlantic_Spring", plot_pool_names)
plot_pool_names <- sub("TysklandS18_Germany_Baltic", "Ariadnegrund_Baltic", plot_pool_names)
plot_pool_names <- sub("PB9_Kattegat_Atlantic_Spring", "PB9_Kattegat (Björköfjorden)_Atlantic_Spring", plot_pool_names)
plot_pool_names <- sub("HGS8_KattegatNorth_Atlantic_Spring", "HGS8_NorthKattegat_Atlantic_Spring", plot_pool_names)
plot_pool_names <- sub("PB10_Skagerrak_Atlantic_Spring", "PB10_Skagerrak (Brofjorden)_Atlantic_Spring", plot_pool_names)
plot_pool_names <- sub("PB2_Iceland_Atlantic_Spring", "PB2_Iceland_Atlantic_Spring", plot_pool_names)
plot_pool_names <- sub("HGS23_Clyde_Atlantic_Spring", "HGS23_Ballantrae (Clyde)_Atlantic_Spring", plot_pool_names)
plot_pool_names <- sub("HGS21_Hebrides_Atlantic_Mixed", "HGS21_West_of_Hebrides_Atlantic_Mixed", plot_pool_names)
plot_pool_names <- sub("HGS17_IsleOfMan_IrishSea_Autumn", "HGS17_Douglas_Bank (Isle_of_Man)_IrishSea_Autumn", plot_pool_names)
#plot_pool_names <- sub("HGS10_Downs_EnglishChannel_Winter", "HGS10_Downs_Atlantic_Winter", plot_pool_names)
plot_pool_names <- sub("N_NorthSea_Atlantic_Autumn", "N_North Sea_Atlantic_Autumn", plot_pool_names)
plot_pool_names <- sub("HWS4_KandalakshaBay_WhiteSea", "HWS4_KandalakshaBay_WhiteSea_Spring", plot_pool_names)
plot_pool_names <- sub("HWS5_KandalakshaBay_WhiteSea", "HWS5_KandalakshaBay_WhiteSea_Summer", plot_pool_names)
plot_pool_names <- sub("Q_Norway_Atlantic_Atlantic_Spring", "Q_Norway_Atlantic_Spring", plot_pool_names)
plot_pool_names <- sub("HWS3_WhiteSea_WhiteSea", "HWS3_White Sea_WhiteSea", plot_pool_names)

plot_pool_names <- gsub("_(WhiteSea|Barents|Atlantic|Baltic|Pacific|Irish|NorthSea|SeaOfJapan|English|Spring|Summer|Autumn|Winter)", ", \\1", plot_pool_names)

plot_pool_names <- gsub("([a-z])([A-Z])", "\\1_\\2", plot_pool_names)

plot_pool_names <- sub("[A-Z0-9]+_", "",plot_pool_names)
plot_pool_names <- gsub("_", " ", plot_pool_names)

pool_name_df <- data.frame(raw = names(pool_freq)[balsfjord_pool_vec], polished = plot_pool_names, stringsAsFactors = F)
save(pool_name_df, file = "~/Projects/Herring/data/Balsfjord/pool_names.RData")

#atl_pdf <- paste0("~/Projects/Herring/doc/Balsfjord/heatmaps/Supp_fig_S6/Atlantic_rec_reg_flip_HM.pdf")
atl_pdf <- paste0("~/Projects/Herring/doc/Balsfjord/heatmaps/Supp_fig_S6/Atlantic_rec_reg_flip_rename_HM.pdf")
pdf(atl_pdf, width = 10, height = 10)
for(i in 1:length(recurring_intro_GR)){
  tmp_freq <- pool_freq[pool_freq$CHROM == as.character(recurring_intro_GR@seqnames[i]) & pool_freq$POS > recurring_intro_GR@ranges@start[i] & pool_freq$POS < recurring_intro_GR@ranges@start[i] + recurring_intro_GR@ranges@width[i],]
  rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
  #tmp_freq$NonSyn <- Baltic_EastAtl_ann_df[match(rownames(tmp_freq), Baltic_EastAtl_ann_df$SNP_id), "NonSyn"]
  #tmp_pdf <- paste0("~/Projects/Herring/doc/Balsfjord/heatmaps/", recurring_intro_GR[i], "_flip_HM.pdf")
  #tmp_pdf <- gsub("[:-]", "_", tmp_pdf)
  #pool_freq_hm(tmp_freq, pdf_file = tmp_pdf, pool_order_vec = balsfjord_pool_vec)
  tmp_reg_name <- paste0(recurring_intro_GR[i]@seqnames, ": ", round(recurring_intro_GR[i]@ranges@start/1e6, digits = 2), " to ",round((recurring_intro_GR[i]@ranges@start + recurring_intro_GR[i]@ranges@width)/1e6, digits = 2), " Mb")
  pool_freq_hm_v2(tmp_freq, pool_order_vec = balsfjord_pool_vec, reg_name = tmp_reg_name, row_lab_vec = plot_pool_names, margins = c(12, 15))
}
dev.off()

#pac_pdf <- paste0("~/Projects/Herring/doc/Balsfjord/heatmaps/Supp_fig_S6/Pacific_rec_reg_flip_HM.pdf")
pac_pdf <- paste0("~/Projects/Herring/doc/Balsfjord/heatmaps/Supp_fig_S6/Pacific_rec_reg_flip_rename_HM.pdf")
pdf(pac_pdf, height = 10, width = 10)
for(i in 1:length(recurring_intro_pac_GR)){
  tmp_freq <- pool_freq[pool_freq$CHROM == as.character(recurring_intro_pac_GR@seqnames[i]) & pool_freq$POS > recurring_intro_pac_GR@ranges@start[i] & pool_freq$POS < recurring_intro_pac_GR@ranges@start[i] + recurring_intro_pac_GR@ranges@width[i],]
  rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
  #tmp_freq$NonSyn <- Baltic_EastAtl_ann_df[match(rownames(tmp_freq), Baltic_EastAtl_ann_df$SNP_id), "NonSyn"]
  #tmp_pdf <- paste0("~/Projects/Herring/doc/Balsfjord/heatmaps/", recurring_intro_pac_GR[i], "_pac_flip_HM.pdf")
  #tmp_pdf <- gsub("[:-]", "_", tmp_pdf)
  #pool_freq_hm(tmp_freq, pdf_file = tmp_pdf, pool_order_vec = balsfjord_pool_vec)
  tmp_reg_name <- paste0(recurring_intro_pac_GR[i]@seqnames, ": ", round(recurring_intro_pac_GR[i]@ranges@start/1e6, digits = 2), " to ",round((recurring_intro_pac_GR[i]@ranges@start + recurring_intro_pac_GR[i]@ranges@width)/1e6, digits = 2), " Mb")
  pool_freq_hm_v2(tmp_freq, pool_order_vec = balsfjord_pool_vec, reg_name = tmp_reg_name,  row_lab_vec = plot_pool_names, margins = c(12, 15))
}
dev.off()

#Table of recurring regions
recurring_intro_GR$type <- "Atlantic"
recurring_intro_pac_GR$type <- "Pacific"
write.table(x = as.data.frame(c(recurring_intro_GR, recurring_intro_pac_GR)), file = "~/Projects/Herring/doc/Balsfjord/Draft/Supplematary_Table_S2.txt", quote = F, sep = "\t", col.names = T, row.names = F)



#Other groups
HWS2_intro_20k_lists <- complile_intro_lists(intr_obj = hws2_v2.0.2_20k_intro, win_df = h79_20k_w, fuse_thresh = 5e4, assoc_t = 8)
summary(unlist(HWS2_intro_20k_lists$atl_red)@ranges@width)
sum(unlist(HWS2_intro_20k_lists$atl_red)@ranges@width)/1e6

HWS3_intro_20k_lists <- complile_intro_lists(intr_obj = hws3_v2.0.2_20k_intro, win_df = h79_20k_w, fuse_thresh = 5e4, assoc_t = 8)
summary(unlist(HWS3_intro_20k_lists$atl_red)@ranges@width)
sum(unlist(HWS3_intro_20k_lists$atl_red)@ranges@width)/1e6

HWS4_intro_20k_lists <- complile_intro_lists(intr_obj = hws4_v2.0.2_20k_intro, win_df = h79_20k_w, fuse_thresh = 5e4, assoc_t = 8)
summary(unlist(HWS4_intro_20k_lists$atl_red)@ranges@width)
sum(unlist(HWS4_intro_20k_lists$atl_red)@ranges@width)/1e6

HWS5_intro_20k_lists <- complile_intro_lists(intr_obj = hws5_v2.0.2_20k_intro, win_df = h79_20k_w, fuse_thresh = 5e4, assoc_t = 8)
summary(unlist(HWS5_intro_20k_lists$atl_red)@ranges@width)
sum(unlist(HWS5_intro_20k_lists$atl_red)@ranges@width)/1e6

HWS_tot <- c(intro_20k_lists$atl_red, HWS2_intro_20k_lists$atl_red, HWS3_intro_20k_lists$atl_red, HWS4_intro_20k_lists$atl_red,HWS5_intro_20k_lists$atl_red)
HWS_tot_Atl_cov <- calc_cov(HWS_tot)

pdf("~/Projects/Herring/doc/Balsfjord/Total_intro_cov.pdf", width = 12)
plot(x = HWS_tot_Atl_cov$global_start, y = HWS_tot_Atl_cov$cov, type = "n", main = "Introgression \"Coverage\"", xlab = "", ylab = "Number of haplotypes", ylim = c(0,40)) 
segments(x0 = HWS_tot_Atl_cov$global_start, x1 = HWS_tot_Atl_cov$global_end, y0 = HWS_tot_Atl_cov$cov, lwd = 4, col = HWS_tot_Atl_cov$col)
for(i in 1:26){
  t_chr <- paste0("chr", i)
  plot(x = HWS_tot_Atl_cov$start[HWS_tot_Atl_cov$group_name ==t_chr], y = HWS_tot_Atl_cov$cov[HWS_tot_Atl_cov$group_name == t_chr], type = "n", main = paste0(t_chr, "; Total Introgression \"Coverage\""), xlab = "", ylab = "Number of haplotypes", ylim = c(0,40)) 
  segments(x0 = HWS_tot_Atl_cov$start[HWS_tot_Atl_cov$group_name ==t_chr], x1 = HWS_tot_Atl_cov$end[HWS_tot_Atl_cov$group_name == t_chr], y0 = HWS_tot_Atl_cov$cov[HWS_tot_Atl_cov$group_name == t_chr], lwd = 4, col = HWS_tot_Atl_cov$col[HWS_tot_Atl_cov$group_name == t_chr])
}
dev.off()

h79_20k_w_GR <- GRanges(seqnames = h79_20k_w$scaffold, ranges = IRanges(start = h79_20k_w$start, end =  h79_20k_w$stop))
tot_rec_win_df <- HWS_tot_Atl_cov[HWS_tot_Atl_cov$cov >= 25,]
tot_rec_win_GR <- GRanges(seqnames = tot_rec_win_df$group_name, ranges = IRanges(start = tot_rec_win_df$start, end =  tot_rec_win_df$end))
tot_rec_win_hits <- findOverlaps(h79_20k_w_GR, tot_rec_win_GR)
h79_20k_w_GR[tot_rec_win_hits@from]
plot_top_regions_2(top_pos = tot_rec_win_hits@from, win_df = h79_20k_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/Balsfjord/HapDist_v2.0.2/Total_recurring_windows/", snp_col = "blue")



#Comparing age with frequency
atl_age_GR <- GRanges(names(atl_age_vec), age = atl_age_vec) #, rec = F
#atl_age_vs_rec <- findOverlaps(atl_age_GR, recurring_intro_GR, type = "within")
#atl_age_GR$rec[unique(atl_age_vs_rec@from)] <- T
Tot_cov_GR <- GRanges(seqnames = HWS_tot_Atl_cov$group_name, ranges = IRanges(start = HWS_tot_Atl_cov$start, end = HWS_tot_Atl_cov$end), cov = HWS_tot_Atl_cov$cov, rec = F, rec_block = NA)
Tot_cov_vs_rec <- findOverlaps(Tot_cov_GR, recurring_intro_GR)
Tot_cov_GR$rec[Tot_cov_vs_rec@from] <- T
Tot_cov_GR$rec_block[Tot_cov_vs_rec@from] <- Tot_cov_vs_rec@to

age_vs_cov <- findOverlaps(atl_age_GR, Tot_cov_GR)
age_vs_cov_GR <- atl_age_GR[age_vs_cov@from]
age_vs_cov_GR$cov_block <- age_vs_cov@to
age_vs_cov_GR$cov <- Tot_cov_GR$cov[age_vs_cov@to]
age_vs_cov_GR$rec <- Tot_cov_GR$rec[age_vs_cov@to]
age_vs_cov_GR$rec_block <- Tot_cov_GR$rec_block[age_vs_cov@to]

plot(x = age_vs_cov_GR$cov, y = log10(age_vs_cov_GR$age), pch = 20, cex = 0.4, col = "grey50")
col_vec <- hcl.colors(max(age_vs_cov_GR$rec_block, na.rm = T))
points(pch = 20, x = age_vs_cov_GR$cov[age_vs_cov_GR$rec], y = log10(age_vs_cov_GR$age[age_vs_cov_GR$rec]), col = col_vec[age_vs_cov_GR$rec_block[age_vs_cov_GR$rec]])

High_tot_cov_GR <- reduce(Tot_cov_GR[Tot_cov_GR$cov > 25], min.gapwidth = 1e3)
for(i in 1:length(High_tot_cov_GR)){
  tmp_freq <- pool_freq[pool_freq$CHROM == as.character(High_tot_cov_GR@seqnames[i]) & pool_freq$POS > High_tot_cov_GR@ranges@start[i] & pool_freq$POS < High_tot_cov_GR@ranges@start[i] + High_tot_cov_GR@ranges@width[i],]
  rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
  tmp_freq$NonSyn <- Baltic_EastAtl_ann_df[match(rownames(tmp_freq), Baltic_EastAtl_ann_df$SNP_id), "NonSyn"]
  tmp_pdf <- paste0("~/Projects/Herring/doc/Balsfjord/heatmaps/", High_tot_cov_GR[i], "_tot_HM.pdf")
  tmp_pdf <- gsub("[:-]", "_", tmp_pdf)
  pool_freq_hm(tmp_freq, pdf_file = tmp_pdf, pool_order_vec = balsfjord_pool_vec)
}



#Chr17, 14.7 Mb (i = 5 above)
#Chr17_14.7_freq <- tmp_freq
missense_GR <- GRanges(seqnames = Chr17_14.7_freq$CHROM[which(Chr17_14.7_freq$NonSyn)], ranges = IRanges(start = Chr17_14.7_freq$POS[which(Chr17_14.7_freq$NonSyn)], end = Chr17_14.7_freq$POS[which(Chr17_14.7_freq$NonSyn)]))
ann_hits <- findOverlaps(missense_GR, cluhar_v2.0.2_gtf)
Chr17_14.7_ann <- cluhar_v2.0.2_gtf[unique(ann_hits@to)]
Chr17_14.7_ann[Chr17_14.7_ann$type == "gene"]
findOverlaps(missense_GR, Chr17_14.7_ann[Chr17_14.7_ann$type == "gene" & Chr17_14.7_ann$gene_id == "ENSCHAG00000021745"])
findOverlaps(missense_GR, Chr17_14.7_ann[Chr17_14.7_ann$type == "gene" & Chr17_14.7_ann$gene_id == "ENSCHAG00000021773"])
findOverlaps(missense_GR, Chr17_14.7_ann[Chr17_14.7_ann$type == "gene" & Chr17_14.7_ann$gene_id == "ENSCHAG00000021832"])

Chr4_19.05Mb_ann <- SNP_eff_intersect(recurring_intro_GR[1], non_ann_df = pool_freq, ann_df = Baltic_EastAtl_ann_df)
pool_freq_hm(Chr4_19.05Mb_ann, pdf_file = "~/Projects/Herring/doc/Balsfjord/heatmaps/Chr4_19.05Mb_ann.pdf", pool_order_vec = balsfjord_pool_vec)
pool_freq_hm(Chr4_19.05Mb_ann[which(Chr4_19.05Mb_ann$NonSyn),], pdf_file = "~/Projects/Herring/doc/Balsfjord/heatmaps/Chr4_19.05Mb_NonSyn.pdf", pool_order_vec = balsfjord_pool_vec)



#Comparing with the Sprat
gene21745_tr1_GR <- Chr17_14.7_ann[Chr17_14.7_ann$type == "exon" & Chr17_14.7_ann$transcript_id == "ENSCHAT00000050871"]
strand(gene21745_tr1_GR) <- "*"
gene21745_tr1 <- DNAStringSet(reverseComplement(unlist(Ch_v2.0.2[gene21745_tr1_GR])))
names(gene21745_tr1) <- "ENSCHAT00000050871"
writeXStringSet(x = gene21745_tr1, filepath = "~/Projects/Herring/data/Balsfjord/gene21745_tr1_cDNA.fa")
#blastn -task blastn -query ./gene21745_tr1_cDNA.fa -outfmt 6 -culling_limit 5 -db ~/Projects/Sprat/data/assemblies/IPA/Sprat_IPA -out gene21745_tr1_cDNA_v_Sprat_IPA.blastout
gene21745_tr1_blastout <- read.table("~/Projects/Herring/data/Balsfjord/gene21745_tr1_cDNA_v_Sprat_IPA.blastout", sep = "\t", stringsAsFactors = F)
hit_filter <- gene21745_tr1_blastout[,2] == "ctg.003110F"
plot(x = range(gene21745_tr1_blastout[hit_filter,7:8]), y = range(gene21745_tr1_blastout[hit_filter,9:10]), type = "n")
segments(x0 =gene21745_tr1_blastout[hit_filter,7], x1 = gene21745_tr1_blastout[hit_filter,8], y0 = gene21745_tr1_blastout[hit_filter,9], y1 = gene21745_tr1_blastout[hit_filter,10])
gene21745_IPA_pri_GR <- reduce(GRanges(seqnames = "ctg.003110F", ranges = IRanges(start = gene21745_tr1_blastout[hit_filter,9], end = gene21745_tr1_blastout[hit_filter,10])))
gene21745_IPA_pri <- DNAStringSet(unlist(IPA_pri[gene21745_IPA_pri_GR]))
names(gene21745_IPA_pri) <- "gene21745_IPA_pri"
writeXStringSet(x = gene21745_IPA_pri, filepath = "~/Projects/Herring/data/Balsfjord/gene21745_IPA_pri.fa")

for(i in length(gene21745_tr1_GR):1){
  if (i == length(gene21745_tr1_GR)){
    gene21745_VCF <- reverseComplement(herring_vcf_extract(gene21745_tr1_GR[i]))
    name_vec <- names(gene21745_VCF)
  }
  if (i != length(gene21745_tr1_GR)){
    tmp_VCF <- reverseComplement(herring_vcf_extract(gene21745_tr1_GR[i]))
    gene21745_VCF <- DNAStringSet(paste0(as.character(gene21745_VCF), as.character(tmp_VCF)))
  }
  names(gene21745_VCF) <-  name_vec
}
writeXStringSet(x = gene21745_VCF, filepath = "~/Projects/Herring/data/Balsfjord/gene21745_from_VCF.fa")
gene21745_IPA_pri_DNAbin <- read.FASTA("~/Projects/Herring/data/Balsfjord/gene21745_IPA_pri.fa")
gene21745_VCF_DNAbin <- read.FASTA("~/Projects/Herring/data/Balsfjord/gene21745_from_VCF.fa")
chr17_14.69_Mb <- c(gene21745_VCF_DNAbin, gene21745_IPA_pri_DNAbin)
chr17_14.69_clustalo <- clustalomega(x = chr17_14.69_Mb)
save(chr17_14.69_clustalo, file = "~/Projects/Sprat/data/tree/chr17_14.69_clustalo.RData")
aln_basic_plots(clustal_obj = chr17_14.69_clustalo, pdf_file = "~/Projects/Sprat/doc/chr17_14.69_Mb_trees.pdf", rooted_version = F)
dist.dna(chr17_14.69_clustalo)
chr17_14.69_diff <- chr17_14.69_clustalo[,which(chr17_14.69_clustalo[1,] != chr17_14.69_clustalo[8,])]
aln_basic_plots(clustal_obj = chr17_14.69_diff, pdf_file = "~/Projects/Sprat/doc/chr17_14.69_Mb_trees.pdf", rooted_version = F)
chr17_14.69_diff_tree <- bionj(dist.aa(chr17_14.69_diff)/length(chr17_14.69_diff[1,]))
plot.phylo(chr17_14.69_diff_tree, type = "unrooted", lab4ut = "axial")

#construct version from missense annotation
miss_ann_raw <- Baltic_EastAtl_ann_df[match(rownames(Chr17_14.7_freq)[which(Chr17_14.7_freq$NonSyn)], Baltic_EastAtl_ann_df$SNP_id),c("CHR", "POS", "raw")]
strsplit_vec <- Vectorize(function(x, y){unlist(strsplit(x, split = "\t"))[y]})
miss_ann_raw$ref_base <- strsplit_vec(miss_ann_raw$raw, 3)
miss_ann_raw$alt_base <- strsplit_vec(miss_ann_raw$raw, 4)
tmp_region_seq <- Ch_v2.0.2["chr17"]
for(i in 1:length(missense_GR)) subseq(tmp_region_seq, start = missense_GR@ranges@start[i], width = 1) <- miss_ann_raw$alt_base[i]
alt_version_cDNA <- DNAStringSet(reverseComplement(unlist(tmp_region_seq[gene21745_tr1_GR])))
names(alt_version_cDNA) <- "Ch_missense"
writeXStringSet(x = alt_version_cDNA, filepath = "~/Projects/Herring/data/Balsfjord/gene21745_Ch_missense.fa")
Missense_DNAbin <- read.FASTA("~/Projects/Herring/data/Balsfjord/gene21745_Ch_missense.fa")
chr17_14.69_Mb <- c(gene21745_VCF_DNAbin, gene21745_IPA_pri_DNAbin, Missense_DNAbin)
chr17_14.69_clustalo <- clustalomega(x = chr17_14.69_Mb)
chr17_14.69_diff <- chr17_14.69_clustalo[,which(chr17_14.69_clustalo[1,] != chr17_14.69_clustalo[13,] & chr17_14.69_clustalo[13,] == chr17_14.69_clustalo[8,])]
checkAlignment(chr17_14.69_diff[c(1,2,3,12,13,5,8,11),])
chr17_14.69_diff_tree <- bionj(dist.aa(chr17_14.69_diff)/length(chr17_14.69_diff[1,]))
plot.phylo(chr17_14.69_diff_tree, type = "unrooted", lab4ut = "axial")


wu_fj_clustalo <- transcript_extraction(transcript_id = "ENSCHAT00000051055", ann_obj = cluhar_v2.0.2_gtf, sprat_target = "ctg.001436F", sprat_dir = "neg", missense_seq = tmp_region_seq)
#wu_fj_diff <- wu_fj_clustalo[,which(wu_fj_clustalo[1,] != wu_fj_clustalo[8,])]
wu_fj_diff <- wu_fj_clustalo[,which(wu_fj_clustalo[3,] != wu_fj_clustalo[13,] & wu_fj_clustalo[13,] == wu_fj_clustalo[8,])]#Missense positions that differ between NSSH and Sea of Japan
wu_fj_diff_tree <- bionj(dist.aa(wu_fj_diff)/length(wu_fj_diff[1,]))

png("~/Projects/Herring/doc/Balsfjord/wu_fj_full_aln.png", width = 1200, height = 600)
image(wu_fj_clustalo)
dev.off()
png("~/Projects/Herring/doc/Balsfjord/wu_fj_missense_aln.png", width = 1200, height = 600)
image(wu_fj_diff[c(1,2,3,12,13,5,8,11),])
dev.off()
pdf("~/Projects/Herring/doc/Balsfjord/wu_fj_missense_tree.pdf")
plot.phylo(wu_fj_diff_tree, type = "unrooted", lab4ut = "axial")
dev.off()

gene21745_clustalo <- transcript_extraction(transcript_id = "ENSCHAT00000050871", ann_obj = cluhar_v2.0.2_gtf, sprat_target = "ctg.003110F", sprat_dir = "pos", missense_seq = tmp_region_seq)
gene21745_diff <- gene21745_clustalo[,which(gene21745_clustalo[3,] != gene21745_clustalo[13,] & gene21745_clustalo[13,] == gene21745_clustalo[8,])]#Missense positions that differ between NSSH and Sea of Japan
checkAlignment(gene21745_diff)

#Mitochondria-linked genes
Mito_GO_gene_ID <- read.table("~/Projects/Herring/data/Balsfjord/mitochondria/Ch_v2.0.2_genes_with_mtGO-0005739_summary.txt", sep = "\t", header = F, stringsAsFactors = F)
Mito_GO_genes <- cluhar_v2.0.2_gtf[cluhar_v2.0.2_gtf$gene_id %in% Mito_GO_gene_ID$V1 & cluhar_v2.0.2_gtf$type == "gene"]
#atl_intro_cov_20k_GR$cov <- atl_intro_cov_20k$cov #Not scaled by size - biased
PvA_GR$atl_intro_cov <- PvA_Fst$atl_intro_cov
PvA_GR$pac_intro_cov <- PvA_Fst$pac_intro_cov
Mito_GO_gene_hits <- findOverlaps(Mito_GO_genes, PvA_GR)

summary(PvA_GR$atl_intro_cov[unique(Mito_GO_gene_hits@to)])
hist(PvA_GR$atl_intro_cov[unique(Mito_GO_gene_hits@to)])
summary(PvA_GR$atl_intro_cov[-unique(Mito_GO_gene_hits@to)])
hist(PvA_GR$atl_intro_cov[-unique(Mito_GO_gene_hits@to)])
wilcox.test(PvA_GR$atl_intro_cov[unique(Mito_GO_gene_hits@to)], PvA_GR$atl_intro_cov[-unique(Mito_GO_gene_hits@to)])

summary(PvA_GR$pac_intro_cov[unique(Mito_GO_gene_hits@to)])
hist(PvA_GR$pac_intro_cov[unique(Mito_GO_gene_hits@to)])
summary(PvA_GR$pac_intro_cov[-unique(Mito_GO_gene_hits@to)])
hist(PvA_GR$pac_intro_cov[-unique(Mito_GO_gene_hits@to)])
wilcox.test(PvA_GR$pac_intro_cov[unique(Mito_GO_gene_hits@to)], PvA_GR$pac_intro_cov[-unique(Mito_GO_gene_hits@to)])



#Support functions
pool_freq_hm_v2 <- function(freq_df, pdf_file = NULL, pool_order_vec, hm_col = c("blue4", "gold1"), margins = c(12,12), reg_name = NULL, row_lab_vec = NULL, ...){
  if(!is.null(pdf_file)) {
    pdf(file = pdf_file, height = 10, width = 10)
  }
  par(xpd = NA)
  if(is.null(row_lab_vec)){
    row_lab_vec <- sub("[A-Z0-9]+_", "", colnames(freq_df)[pool_order_vec])
  }
  
  target_snp_id <- rownames(freq_df)
  #pool_order_vec <- c(17, 12, 20, 6, 7, 8, 18, 5, 4, 21, 19, 10, 9, 3, 2, 11, 14, 15, 16, 13, 1)
  cr <- colorRamp(hm_col)
  crp = colorRampPalette(colors = hm_col)(500)
  col_colors <- rep("black", dim(freq_df)[1])
  col_colors[freq_df$NonSyn] <- "darkorchid"
  heatmap(t(as.matrix(freq_df[,pool_order_vec])), scale = "none", Colv = NA, Rowv = NA, labCol = target_snp_id, labRow = row_lab_vec, col = crp, margins = margins, ...) #ColSideColors = col_colors
  scale_x_vec <- seq(from = par("usr")[1], to = par("usr")[2]*0.5, length.out = 505)
  rect(ybottom = par("usr")[3] - (par("usr")[4]-par("usr")[3])*0.06, xleft = scale_x_vec[-(501:505)], ytop = par("usr")[3] - (par("usr")[4] -par("usr")[3])*0.04, xright =  scale_x_vec[-(1:5)], col = rgb(cr((1:500)/500), maxColorValue=255), border = NA)
  text(x=scale_x_vec[c(10,253,505)], y = par("usr")[3] - (par("usr")[4]-par("usr")[3])*0.08, labels = paste0(c(0, 50, 100), "%"))
  if(!is.null(reg_name)) {
    text(x=scale_x_vec[253], y = par("usr")[3] - (par("usr")[4]-par("usr")[3])*0.02, labels = reg_name, cex = 1.3)
  }
  if(!is.null(pdf_file)) {
    dev.off()
  }
  #return(par("usr"))
}


SNP_eff_intersect <- function(target_GR, non_ann_df, ann_df){
  tmp_freq <- non_ann_df[non_ann_df$CHROM == as.character(target_GR@seqnames) & non_ann_df$POS > target_GR@ranges@start & non_ann_df$POS < target_GR@ranges@start + target_GR@ranges@width,]
  rownames(tmp_freq) <- paste(tmp_freq$CHROM, tmp_freq$POS, sep = "_")
  tmp_freq$NonSyn <- ann_df[match(rownames(tmp_freq), ann_df$SNP_id), "NonSyn"]
  return(tmp_freq)
}

#HaploDist-type plots with clade re-ordering
plot_top_regions_3 <- function(top_pos, win_df, geno, sample_list = NULL, plot_dir = "", snp_col = "orange", snp_lwd = 0.1, snp_alpha = 50, dend_weight_vec = NULL){
  require(stringdist)
  require(diptest)
  hclust_list <- list()
  dist_list <- list()
  for (pos in top_pos){
    win <- win_df[pos,]
    scaff_geno <- geno[geno[,1] ==  win[1,1],]
    target_markers <- scaff_geno[,2] >= win[1,2] & scaff_geno[,2] < win[1,3]
    target_geno <- scaff_geno[target_markers,]
    n_haps <- nchar(target_geno[1,3])
    hap_vec <- array()
    for(i in 1:n_haps){
      hap_vec[i] <- paste(substr(target_geno[,3], i, i), collapse = "")
    }
    hap_dist_mat <- outer(hap_vec, hap_vec, FUN = "stringdist", method = "hamming")
    diag(hap_dist_mat) <- NA
    hap_dist_data <- array(data = hap_dist_mat[!is.na(hap_dist_mat)])
    
    if (is.null(sample_list)) sample_list <- 1:n_haps
    rownames(hap_dist_mat) <- sample_list
    colnames(hap_dist_mat) <- sample_list
    win_loc <- paste(win[1], "_", round(win[2]/1000), "_to_", round(win[3]/1000), "_kb", sep = "")
    pdf(paste(plot_dir, win_loc, ".pdf", sep = ""), width = 14, height = 10, paper = "a4r")
    
    #Distance distribution histogram
    hist(hap_dist_data, main = win_loc, xlab = "Edit distance")
    
    #Dendrogram
    hap_hc <- hclust(as.dist(hap_dist_mat))
    plot(hap_hc, xlab ="", ann = F, cex = min(c(1, 60/length(sample_list))))
    title(main = win_loc)
    
    #Heatmap
    #block_hm <- heatmap(x= hap_dist_mat, scale = "none", margins = c(4,4), Rowv = as.dendrogram(hap_hc), Colv = as.dendrogram(hap_hc), cexRow = min(c(1, 35/length(sample_list))), labCol = "")
    reorder_dendro <- reorder(as.dendrogram(hap_hc), wts = dend_weight_vec)
    block_hm <- heatmap(x= hap_dist_mat, scale = "none", margins = c(4,4), Rowv = reorder_dendro, Colv = reorder_dendro, cexRow = min(c(1, 35/length(sample_list))), labCol = "")
    
    #Haplotype visualisation
    ext_start <- max(pos-5, min(which(win_df[,1] == win[1,1])))
    ext_stop <- min(pos+5, max(which(win_df[,1] == win[1,1])))
    ext_target_m <- scaff_geno[,2] >= win_df[ext_start,2] & scaff_geno[,2] < win_df[ext_stop,3]
    ext_target_g <- scaff_geno[ext_target_m,]
    for(i in 1:n_haps){
      hap_vec[i] <- paste(substr(ext_target_g[,3], i, i), collapse = "")
    }
    plot(1:n_haps, type = "n", xlim = range(ext_target_g[,2]), axes = F, ann = F)
    axis(2, labels = sample_list[block_hm$rowInd], at = 1:n_haps, las = 2, cex.axis= min(c(0.6, 20/n_haps)), col=1:2, pos = range(ext_target_g[,2])[1] - 0.0175*diff(range(ext_target_g[,2])))
    axis(1)
    title(main = win_loc)
    rect(range(ext_target_g[,2])[1]- 0.0025*diff(range(ext_target_g[,2])),0,range(ext_target_g[,2])[2] + 0.0025*diff(range(ext_target_g[,2])), n_haps+1, col = "grey90")
    rect(range(target_geno[,2])[1],0.1,range(target_geno[,2])[2], n_haps+0.9, col = "grey80", border = NA)
    #rect(range(ext_target_g[,2])[1] - 0.015*diff(range(ext_target_g[,2])), which(block_hm$rowInd == 1)-0.4, range(ext_target_g[,2])[1] - 0.005*diff(range(ext_target_g[,2])), which(block_hm$rowInd == 1) +0.4, col = "orange")
    #ref_col <- factor(unlist(strsplit(hap_vec[1], "")), levels = c("G", "C", "T", "A"))
    allele_df <- data.frame("G" = 1:dim(ext_target_g)[1], "C" = 0, "T" = 0, "A" = 0, "major" = NA, stringsAsFactors = F)
    allele_df[,"G"] <- nchar(gsub("[^G]", "",ext_target_g[,3]))
    allele_df[,"C"] <- nchar(gsub("[^C]", "",ext_target_g[,3]))
    allele_df[,"T"] <- nchar(gsub("[^T]", "",ext_target_g[,3]))
    allele_df[,"A"] <- nchar(gsub("[^A]", "",ext_target_g[,3]))			
    allele_max_vec <- pmax(allele_df[,"G"], allele_df[,"C"], allele_df[,"T"], allele_df[,"A"])
    allele_df[allele_df[,"G"] == allele_max_vec,"major"] <- "G"
    allele_df[allele_df[,"C"] == allele_max_vec,"major"] <- "C"
    allele_df[allele_df[,"T"] == allele_max_vec,"major"] <- "T"
    allele_df[allele_df[,"A"] == allele_max_vec,"major"] <- "A"
    ref_col <- factor(allele_df[,"major"], levels = c("G", "C", "T", "A"))
    for(i in 1:n_haps){
      tmp_col <- factor(unlist(strsplit(hap_vec[i], "")), levels = c("G", "C", "T", "A"))
      y_vec <- rep(which(block_hm$rowInd == i), length(tmp_col))[tmp_col != ref_col]
      x_vec <- ext_target_g[,2][tmp_col != ref_col]
      #rect(x_vec, y_vec-0.4 , x_vec + 0.001*diff(range(ext_target_g[,2])), y_vec+0.4, col = rgb(t(col2rgb("orange")), alpha = 100, max= 255), border = NA)
      segments(x0 = x_vec, y0 = y_vec-0.4, y1 = y_vec+0.4, col = rgb(t(col2rgb(snp_col)),alpha = min(c(snp_alpha,100*(10000/length(tmp_col)))) , max= 255), lwd = snp_lwd)
      #alpha = min(c(50,100*(1000/length(x_vec))))
    }
    segments(x0 = range(target_geno[,2]), y0 = 0.05, y1 = n_haps + 0.95)
    dev.off()
    hclust_list[[win_loc]] <- hap_hc
    dist_list[[win_loc]] <- hap_dist_mat
  }
  return(invisible(list(hclust = hclust_list, dist_mat = dist_list)))
}



transcript_extraction <- function(transcript_id, ann_obj, sprat_target, sprat_dir, missense_seq = NULL){
  tr_GR <- ann_obj[ann_obj$type == "exon" & ann_obj$transcript_id == transcript_id]
  tr_GR <- tr_GR[order(tr_GR@ranges@start)]
  if(as.character(tr_GR@strand[1]) == "+") tr_dir <- "pos"
  if(as.character(tr_GR@strand[1]) == "-") tr_dir <- "neg"
  strand(tr_GR) <- "*"
  if(tr_dir == "neg") tr1 <- DNAStringSet(reverseComplement(unlist(Ch_v2.0.2[tr_GR])))
  if(tr_dir == "pos") tr1 <- DNAStringSet(unlist(Ch_v2.0.2[tr_GR]))
  names(tr1) <- transcript_id
  writeXStringSet(x = tr1, filepath = "~/Projects/Herring/data/Balsfjord/tmp_files/tr1_cDNA.fa")
  system_cmd <- "~/Software/ncbi-blast-2.2.29+/bin/blastn -task blastn -query ~/Projects/Herring/data/Balsfjord/tmp_files/tr1_cDNA.fa -outfmt 6 -culling_limit 5 -db ~/Projects/Sprat/data/assemblies/IPA/Sprat_IPA -out ~/Projects/Herring/data/Balsfjord/tmp_files/tr1_cDNA_v_Sprat_IPA.blastout"
  system(system_cmd)
  tr1_blastout <- read.table("~/Projects/Herring/data/Balsfjord/tmp_files/tr1_cDNA_v_Sprat_IPA.blastout", sep = "\t", stringsAsFactors = F)
  hit_filter <- tr1_blastout[,2] == sprat_target
  if(sprat_dir == "pos") {
    IPA_pri_GR <- reduce(GRanges(seqnames = sprat_target, ranges = IRanges(start = tr1_blastout[hit_filter,9], end = tr1_blastout[hit_filter,10])))
    IPA_pri_tr <- DNAStringSet(unlist(IPA_pri[IPA_pri_GR]))
  }
  if(sprat_dir == "neg"){
    IPA_pri_GR <- reduce(GRanges(seqnames = sprat_target, ranges = IRanges(start = tr1_blastout[hit_filter,10], end = tr1_blastout[hit_filter,9])))
    IPA_pri_tr <- reverseComplement(DNAStringSet(unlist(IPA_pri[IPA_pri_GR])))
  }
  
  names(IPA_pri_tr) <- "IPA_pri_tr"
  writeXStringSet(x = IPA_pri_tr, filepath = "~/Projects/Herring/data/Balsfjord/tmp_files/IPA_pri_tr.fa")
  
  if(tr_dir == "neg"){
    for(i in length(tr_GR):1){
      if (i == length(tr_GR)){
        tr_VCF <- reverseComplement(herring_vcf_extract(tr_GR[i]))
        name_vec <- names(tr_VCF)
      }
      if (i != length(tr_GR)){
        tmp_VCF <- reverseComplement(herring_vcf_extract(tr_GR[i]))
        tr_VCF <- DNAStringSet(paste0(as.character(tr_VCF), as.character(tmp_VCF)))
      }
      names(tr_VCF) <-  name_vec
    }
  }
  if(tr_dir == "pos"){
    for(i in 1:length(tr_GR)){
      if (i == 1){
        tr_VCF <- herring_vcf_extract(tr_GR[i])
        name_vec <- names(tr_VCF)
      }
      if (i != 1){
        tmp_VCF <- herring_vcf_extract(tr_GR[i])
        tr_VCF <- DNAStringSet(paste0(as.character(tr_VCF), as.character(tmp_VCF)))
      }
      names(tr_VCF) <-  name_vec
    }
  }
  
  writeXStringSet(x = tr_VCF, filepath = "~/Projects/Herring/data/Balsfjord/tmp_files/tr_from_VCF.fa")
  
  IPA_pri_DNAbin <- read.FASTA("~/Projects/Herring/data/Balsfjord/tmp_files/IPA_pri_tr.fa")
  VCF_DNAbin <- read.FASTA("~/Projects/Herring/data/Balsfjord/tmp_files/tr_from_VCF.fa")
  tr_combined <- c(VCF_DNAbin, IPA_pri_DNAbin)
  if(!is.null(missense_seq)){
    if(tr_dir == "pos") alt_version_cDNA <- DNAStringSet(unlist(missense_seq[tr_GR]))
    if(tr_dir == "neg") alt_version_cDNA <- DNAStringSet(reverseComplement(unlist(missense_seq[tr_GR])))
    names(alt_version_cDNA) <- "Ch_missense"
    writeXStringSet(x = alt_version_cDNA, filepath = "~/Projects/Herring/data/Balsfjord/tmp_files/missense.fa")
    missense_DNAbin <- read.FASTA("~/Projects/Herring/data/Balsfjord/tmp_files/missense.fa")
    tr_combined <- c(VCF_DNAbin, IPA_pri_DNAbin, missense_DNAbin)
  }
  
  tr_clustalo <- clustalomega(x = tr_combined)
  return(tr_clustalo)
}

aln_basic_plots <- function(clustal_obj, pdf_file, rooted_version = T){
  #old.par <- par(no.readonly = T)
  png(filename = sub("pdf", "png", pdf_file), width = 1000, height = 1000)
  checkAlignment(clustal_obj)
  dev.off()
  #par(old.par)
  pdf(pdf_file)
  SvH_tree <- bionj(dist.dna(clustal_obj))
  plot.phylo(SvH_tree, type = "unrooted", lab4ut = "axial")
  if(rooted_version){
    SvH_tree_rooted <- root(SvH_tree, outgroup = grep("ctg", SvH_tree$tip.label, value = T), resolve.root = T)
    plot.phylo(SvH_tree_rooted)
  }
  dev.off()
}

herring_vcf_extract <- function(herring_GR, herring_samples = c("BM15_HastKar_Baltic_Spring", "Gavle54_Gavle_Baltic_Autumn", "NSSH34_Norway_Atlantic_Spring", "Z14_IsleofMan_Atlantic_Autumn", "HWS22_PechoraSea_BarentsSea", "HWS42_KandalakshaBay_WhiteSea_Spring", "HWS53_KandalakshaBay_WhiteSea_Summer", "HWS12_Japan_SeaofJapan", "HWS32_WhiteSea_WhiteSea", "HWS63_Balsfjord_Atlantic", "Pacific3_Vancouver_Pacific")){
  require(Biostrings)
  require(GenomicRanges)
  require(ape)
  old_wd <-  getwd()
  setwd("~/Projects/Sprat/data/tree/")
  if(file.exists("./tmp_H_extract.fasta")){
    file.remove("./tmp_H_extract.fasta")
    print("Cleared old temp file!")
  } 
  for(herring_sample in herring_samples){
    system_cmd <- paste0("export PATH=$PATH:/Users/mapet205/Software/local/bin; samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta ", herring_GR, " | bcftools consensus -H 1 --sample ", herring_sample, " ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz >> ./tmp_H_extract.fasta")
    system(system_cmd)
  }
  
  raw_seq_list <- readDNAStringSet("./tmp_H_extract.fasta")
  names(raw_seq_list)[1:length(herring_samples)] <- herring_samples
  setwd(old_wd)
  return(raw_seq_list)
}



calc_cov <- function(GR_list, size_df = Ch_v2.0.2_sizes){
  cov_obj <- coverage(unlist(GR_list))
  cov_df <- as.data.frame(ranges(cov_obj))
  cov_df[,"cov"] <- unlist(runValue(cov_obj))
  cov_df[,"global_start"] <- cov_df[,"start"] + size_df[match(cov_df[, "group_name"], size_df[,"name"]), "offset"]
  cov_df[,"global_end"] <-  cov_df[,"end"] + size_df[match(cov_df[, "group_name"], size_df[,"name"]), "offset"]
  cov_df[,"col"] <-  size_df[match(cov_df[, "group_name"], size_df[,"name"]), "col"]
  return(cov_df)
}

cov_and_hap_plot <- function(region_GR, cov_df = atl_intro_cov_20k, intro_list = intro_20k_lists$atl_red, hap_order_vec = NULL, pdf_file = NULL){
  introgression_pos <- region_GR@ranges@start + region_GR@ranges@width/2
  focal_GR <- GRanges(seqnames = region_GR@seqnames, ranges = IRanges(start = introgression_pos, end = introgression_pos))
  target_intros <- GRanges()
  for(i in 1:length(intro_list)){
    target_intros[names(intro_list)[i]] <- intro_list[[i]][findOverlaps(focal_GR, intro_list[[i]])@to]
  }
  
  x_lim <- c(min(target_intros@ranges@start) - 5e5, max(target_intros@ranges@start + target_intros@ranges@width) + 5e5)
  target_chr <- as.character(region_GR@seqnames)
  if(!is.null(pdf_file)) pdf(pdf_file)
  plot(x = 1, y = 1, type = "n", main = "Introgression \"Coverage\" & haplotypes", xlab = "", ylab = "Number of haplotypes", xlim = x_lim, ylim = c(-9,9))   
  rect(xleft =  region_GR@ranges@start, xright =  region_GR@ranges@start +  region_GR@ranges@width, ytop = 11, ybottom = -11, col = "grey90", border = NA)
  segments(x0 = cov_df$start[cov_df$group_name == target_chr], x1 = cov_df$end[cov_df$group_name == target_chr], y0 = cov_df$cov[cov_df$group_name == target_chr], lwd = 4, col = "grey70")
  #hap_order_vec <- c(1,4,8,6,3,5,2,7) #Haplotypes grouped by likely IBD, and dcreasing decay
  #hap_col_vec <- c("black","red", "blue")[c(1,2,1,1,1,1,3,1)]
  hap_col_vec <- rep(c("red", "blue"), times = ceiling(length(intro_list)/2))
  if(is.null(hap_order_vec))hap_order_vec <- 1:length(intro_list)
  for(hap in hap_order_vec){
    target_segments <- which(intro_list[[hap]]@seqnames==target_chr)
    x0_vec <- intro_list[[hap]]@ranges@start[target_segments]
    x1_vec <- intro_list[[hap]]@ranges@width[target_segments] +  x0_vec
    y_pos <- -which(hap_order_vec == hap)
    segments(x0 = x0_vec, x1 =x1_vec, y0 = y_pos, lwd = 4, col = hap_col_vec[hap])
    text(x = x_lim[1], y = y_pos+0.3, labels = names(intro_20k_lists$atl_red)[hap], pos = 4, cex = 0.6, col = hap_col_vec[hap])
  }
  if(!is.null(pdf_file)) dev.off()
}

fst_v_intro_plot <- function(target_chr, cov_df, full_Fst, target_reg_GR, Fst_cutoff = 0.6, size_df = Ch_v2.0.2_sizes, flank_ext = -1){
  #target_chr <- "chr23"
  high_Fst <- which(full_Fst$AtlPac_Fst > Fst_cutoff)
  high_Fst_GR <- GRanges(full_Fst$chr, ranges = IRanges(start = full_Fst$BIN_START, end = full_Fst$BIN_END))[high_Fst]
  Target_v_high_Fst <- findOverlaps(target_reg_GR, high_Fst_GR, maxgap = flank_ext)
  high_Fst_df <- as.data.frame(high_Fst_GR, stringsAsFactors = F)
  high_Fst_df[,"global_start"] <- high_Fst_df[,"start"] + size_df[match(high_Fst_df[, "seqnames"], size_df[,"name"]), "offset"]
  high_Fst_df[,"global_end"] <-  high_Fst_df[,"end"] + size_df[match(high_Fst_df[, "seqnames"], size_df[,"name"]), "offset"]
  high_Fst_df[,"atl_rec_overlap"] <- F
  high_Fst_df[Target_v_high_Fst@to,"atl_rec_overlap"] <- T
  
  plot(x = cov_df$start[cov_df$group_name == target_chr], y = cov_df$cov[cov_df$group_name == target_chr], type = "n", main = paste("Introgression \"Coverage\";", target_chr ), xlab = "", ylab = "Number of haplotypes", ylim = c(0,11)) 
  segments(x0 = cov_df$start[cov_df$group_name == target_chr], x1 = cov_df$end[cov_df$group_name == target_chr], y0 = cov_df$cov[cov_df$group_name == target_chr], lwd = 3, col = cov_df$col[cov_df$group_name == target_chr])
  target_blocks <- high_Fst_df$seqnames == target_chr & !high_Fst_df$atl_rec_overlap
  if(sum(target_blocks > 0)) segments(x0 = high_Fst_df$start[target_blocks ], x1 = high_Fst_df$end[target_blocks], y0 = 10.5, lwd = 3, col = "blue")
  target_blocks <- high_Fst_df$seqnames == target_chr & high_Fst_df$atl_rec_overlap
  if(sum(target_blocks > 0)) segments(x0 = high_Fst_df$start[target_blocks ], x1 = high_Fst_df$end[target_blocks], y0 = 10.7, lwd = 3, col = "red")
  segments(x0 = full_Fst$BIN_START[full_Fst$chr == target_chr], x1 = full_Fst$BIN_END[full_Fst$chr == target_chr], y0 = full_Fst$AtlPac_Fst[full_Fst$chr == target_chr] + 8.5, lwd = 3, col = "darkorchid")
}

fst_v_intro_plot2 <- function(target_chr, cov_df1, cov_df2, full_Fst, target_reg_GR1, target_reg_GR2, Fst_cutoff = 0.6, size_df = Ch_v2.0.2_sizes, flank_ext = -1){
  #target_chr <- "chr23"
  high_Fst <- which(full_Fst$AtlPac_Fst > Fst_cutoff)
  high_Fst_GR <- GRanges(full_Fst$chr, ranges = IRanges(start = full_Fst$BIN_START, end = full_Fst$BIN_END))[high_Fst]
  T1_v_high_Fst <- findOverlaps(target_reg_GR1, high_Fst_GR, maxgap = flank_ext)
  T2_v_high_Fst <- findOverlaps(target_reg_GR2, high_Fst_GR, maxgap = flank_ext)
  high_Fst_df <- as.data.frame(high_Fst_GR, stringsAsFactors = F)
  #high_Fst_df[,"global_start"] <- high_Fst_df[,"start"] + size_df[match(high_Fst_df[, "seqnames"], size_df[,"name"]), "offset"]
  #high_Fst_df[,"global_end"] <-  high_Fst_df[,"end"] + size_df[match(high_Fst_df[, "seqnames"], size_df[,"name"]), "offset"]
  high_Fst_df[,"T1_rec_overlap"] <- F
  high_Fst_df[T1_v_high_Fst@to,"T1_rec_overlap"] <- T
  high_Fst_df[,"T2_rec_overlap"] <- F
  high_Fst_df[T2_v_high_Fst@to,"T2_rec_overlap"] <- T
  
  plot(x = 1, y = 1, type = "n", main = paste("Introgression vs Fst;", target_chr ), xlab = "", ylab = "Number of haplotypes", ylim = c(0,2.5), xlim = c(0,size_df[match(target_chr, size_df[,"name"]), "size"]), axes = F)
  axis(1)
  abline(h = c(1.05, 1.25), col = "grey50")
  legend(x = "topleft", legend = c("Fst profile", "High Fst" , "High Fst and Atl", "High Fst and Pac"), lwd = 3, col = c("darkorchid", "steelblue", "darkorange3", "olivedrab3"), cex = 0.8)
  legend(x = "topright", legend = c("Atlantic introgressions", "Pacific introgressions"), lwd = 3, col = c("darkorange2", "olivedrab2"), cex = 0.8)
  segments(x0 = cov_df1$start[cov_df1$group_name == target_chr], x1 = cov_df1$end[cov_df1$group_name == target_chr], y0 = cov_df1$cov[cov_df1$group_name == target_chr]/8, lwd = 3, col = "darkorange2")
  segments(x0 = cov_df2$start[cov_df2$group_name == target_chr], x1 = cov_df2$end[cov_df2$group_name == target_chr], y0 = (cov_df2$cov[cov_df2$group_name == target_chr]/8)+0.02, lwd = 3, col = "olivedrab2")
  target_blocks <- high_Fst_df$seqnames == target_chr & !high_Fst_df$T1_rec_overlap & !high_Fst_df$T2_rec_overlap
  if(sum(target_blocks > 0)) segments(x0 = high_Fst_df$start[target_blocks], x1 = high_Fst_df$end[target_blocks], y0 = 1.15, lwd = 3, col = "steelblue")
  target_blocks <- high_Fst_df$seqnames == target_chr & high_Fst_df$T1_rec_overlap
  if(sum(target_blocks > 0)) segments(x0 = high_Fst_df$start[target_blocks ], x1 = high_Fst_df$end[target_blocks], y0 = 1.20, lwd = 4, col = "darkorange3")
  target_blocks <- high_Fst_df$seqnames == target_chr & high_Fst_df$T2_rec_overlap
  if(sum(target_blocks > 0)) segments(x0 = high_Fst_df$start[target_blocks ], x1 = high_Fst_df$end[target_blocks], y0 = 1.10, lwd = 4, col = "olivedrab3")
  segments(x0 = full_Fst$BIN_START[full_Fst$chr == target_chr], x1 = full_Fst$BIN_END[full_Fst$chr == target_chr], y0 = full_Fst$AtlPac_Fst[full_Fst$chr == target_chr] + 1.3, lwd = 3, col = "darkorchid")
  abline(h = Fst_cutoff + 1.3, col = "grey40")
}


complile_intro_lists <- function(intr_obj, win_df,  fuse_thresh = 2e5, assoc_t = 5){
  ratio_cols <- grep("ratio", names(intr_obj$dist))
  #GRanges-based quantifications
  a_GR_list <- GRangesList()
  a_GR_red_list <- GRangesList()
  p_GR_list <- GRangesList()
  p_GR_red_list <- GRangesList()
 
  for (i in 1:length(ratio_cols)) {
    hap_name <- sub("_ratio", "", names(intr_obj$dist)[ratio_cols[i]])
    dist_vec <- intr_obj$dist[, ratio_cols[i]]
    a_win <- which(dist_vec <= median(dist_vec, na.rm = T) / assoc_t)
    a_GR <-
      GRanges(
        seqnames = win_df$scaffold[a_win],
        ranges =  IRanges(start = win_df$start[a_win], end = win_df$stop[a_win])
      )
    p_win <- which(dist_vec >= median(dist_vec, na.rm = T) * assoc_t)
    p_GR <-
      GRanges(
        seqnames = win_df$scaffold[p_win],
        ranges =  IRanges(start = win_df$start[p_win], end = win_df$stop[p_win])
      )
    a_GR_list[[hap_name]] <- a_GR
    a_GR_red_list[[hap_name]] <- reduce(a_GR, min.gapwidth = fuse_thresh)
    p_GR_list[[hap_name]] <- p_GR
    p_GR_red_list[[hap_name]] <- reduce(p_GR, min.gapwidth = fuse_thresh)
  }
  return(list(atl = a_GR_list, atl_red = a_GR_red_list, pac = p_GR_list, pac_red = p_GR_red_list))
}

intro_haplotype_plot <- function(intro_obj, rec_profile, a_cov, p_cov, pdf_file = "~/Projects/Herring/doc/Balsfjord/haps.pdf", p_thresh, a_thresh){
  pdf(width = 10, file = pdf_file)
  for(chr in 1:26){
    cur_chr <- paste0("chr", chr)
    for(i in 1:length(ratio_cols)){
      y_vec <- log10(intro_obj$dist[intro_obj$dist$scaffold == cur_chr,ratio_cols[i]])
      main_title = paste(names(intro_obj$dist)[ratio_cols[i]], "; Chr ", chr, sep = "") 
      plot(x = intro_obj$dist$start[intro_obj$dist$scaffold == cur_chr], y = y_vec, ylim = c(-4,7.5), col = i, lwd = 2, main = main_title, xlab = "Position", ylab = "log10(Distance Ratio)") 
      abline(h = c(0, log10(p_thresh), log10(a_thresh), mean(y_vec, na.rm = T)), col = c("black", "grey70", "grey70",i), lty = 1:3)
      points(x = rec_profile$BIN_START[rec_profile$CHR==chr], y = (rec_profile$RR[rec_profile$CHR==chr]/4)-4, col = "darkorchid")
      segments(x0 = p_cov$start[p_cov$group_name == cur_chr], x1 = p_cov$end[p_cov$group_name == cur_chr], y0 = (p_cov$cov[p_cov$group_name == cur_chr])/8+3.5, lwd = 2, col = "darkorange")
      segments(x0 = a_cov$start[a_cov$group_name == cur_chr], x1 = a_cov$end[a_cov$group_name == cur_chr], y0 = (a_cov$cov[a_cov$group_name == cur_chr])/8+4.6, lwd = 2, col = "darkgreen")
      
      legend(x = "topleft", pch = c(1,1,NA, NA), lty = c(NA, NA, 1, 1), col = c(i, "darkorchid", "darkorange", "darkgreen"), legend = c("Distance ratio", "Recombination rate (Not to scale)", "Pacific introgressions", "Atlantic introgressions"))
    }
  }
  dev.off()
}

introgression_contrast <- function(target = "HWS-6", ref_1 = "A[MF]", ref_2 = "HWS-1", w_df, sample_list, geno){
  #target_haps <- grep("HWS-6", h67_clean_samples)
  #ref_atl_haps <- grep("A[MF]", h67_clean_samples)
  #ref_pac_haps <- grep("Paci", h67_clean_samples)
  #ref_pac_haps <- grep("HWS-1", h67_clean_samples)
  target_haps <- grep(target, sample_list)
  ref_1_haps <- grep(ref_1, sample_list)
  ref_2_haps <- grep(ref_2, sample_list)
  
  #hws6_dist_df <- h67_w
  dist_df <- w_df
  
  for (w_idx in 1:dim(dist_df)[1]){
    win <- dist_df[w_idx,]
    scaff_geno <- geno[geno[,1] ==  win[1,1],]
    target_markers <- scaff_geno[,2] >= win[1,2] & scaff_geno[,2] < win[1,3]
    target_geno <- scaff_geno[target_markers,]
    if(dim(target_geno)[1] > 0){
      n_haps <- nchar(target_geno[1,3])
      hap_vec <- array()
      for(i in 1:n_haps){
        hap_vec[i] <- paste(substr(target_geno[,3], i, i), collapse = "")
      }
      hap_dist_mat <- outer(hap_vec, hap_vec, FUN = "stringdist", method = "hamming")
      diag(hap_dist_mat) <- NA
    } else{
      print(paste("Win", w_idx, "empty; skipping."))
    }
    for(th in target_haps){
      #hws6_dist_df[w_idx, paste(h67_clean_samples[th],"_atl", sep = "")] <- mean(hap_dist_mat[th, ref_atl_haps])
      #hws6_dist_df[w_idx, paste(h67_clean_samples[th],"_pac", sep = "")] <- mean(hap_dist_mat[th, ref_pac_haps])
      dist_df[w_idx, paste(sample_list[th],"_1", sep = "")] <- min(hap_dist_mat[th, ref_1_haps])
      dist_df[w_idx, paste(sample_list[th],"_2", sep = "")] <- min(hap_dist_mat[th, ref_2_haps])
    }
  }
  return(dist_df)
}

introgression_plot <- function(dist_df, sample_list,  snp_numbers, assoc_tresh = 100, snp_cutoff = 200, min_diff = 100, pdf_file = "", eps = 1e-3){
  target <- sub("[0-9][._].+", "", names(dist_df)[6])
  if (pdf_file != "") pdf(file = pdf_file, width = 15, height = 5)
  target_haps <- grep(target, sample_list)
  plot(x= 1, y=1, xlim = c(0,max(dist_df[,4])), ylim = c(log10(eps)-3, -log10(eps)+3), type = "n", main = paste(target,"x vs Atlantic & Sea of Japan", sep = ""), xlab = "Cumulative position", ylab = "Log10 of distance ratio")
  hap_ratio_df <- data.frame(hap = sample_list[target_haps], stringsAsFactors= F)
  #assoc_tresh <- 100
  #snp_cutoff <- 200
  #min_diff <- 100
  
  for(bf_hap in 1:dim(hap_ratio_df)[1]){
    
    bf_ind <- ceiling(bf_hap/2)
    hap_no <- 2 - bf_hap%%2
    hap_ratio_df[bf_hap, "ref_1_mean"] <- mean(dist_df[,2*(bf_hap-1)+6])
    hap_ratio_df[bf_hap, "ref_2_mean"] <- mean(dist_df[,2*(bf_hap-1)+7])
    mean_ratio <- hap_ratio_df[bf_hap, "ref_1_mean"]/hap_ratio_df[bf_hap, "ref_2_mean"]
    
    rv <- (dist_df[,2*(bf_hap-1)+6]+eps)/(dist_df[,2*(bf_hap-1)+7]+eps)
    rv[snp_numbers <= snp_cutoff|(dist_df[,2*(bf_hap-1)+6] < min_diff & dist_df[,2*(bf_hap-1)+7] < min_diff)] <- NA #1
    
    dist_df[,paste(hap_ratio_df[bf_hap, "hap"], "_ratio", sep = "")] <- rv
    assoc_filter <- which(rv > assoc_tresh * mean_ratio | rv < mean_ratio/assoc_tresh) #which(rv > assoc_tresh | rv < 1/assoc_tresh)
    hap_ratio_df[bf_hap, "pac"] <- sum(rv > mean_ratio*assoc_tresh, na.rm = T)
    hap_ratio_df[bf_hap, "atl"] <- sum(rv < mean_ratio/assoc_tresh, na.rm = T)
    bg_col_vec <- c("grey20", "grey80")[match(as.integer(as.factor(dist_df[, 1])), unique(as.integer(as.factor(dist_df[, 1]))))%%2 + 1]
    points(y = log10(dist_df[-assoc_filter,paste(hap_ratio_df[bf_hap, "hap"], "_ratio", sep = "")]), x = dist_df[-assoc_filter,4], col = bg_col_vec[-assoc_filter], pch = 16, cex=0.2)
    points(y = log10(dist_df[assoc_filter,paste(hap_ratio_df[bf_hap, "hap"], "_ratio", sep = "")]), x = dist_df[assoc_filter,4], col = rgb(t(col2rgb(bf_ind+1)), alpha = 120, maxColorValue = 255), pch = 15+hap_no, cex = 0.8)
  }
  if (pdf_file != "")dev.off()
  return(list(dist = dist_df, hap= hap_ratio_df))
}

introgression_plot_2 <- function(dist_df, sample_list,  snp_numbers, assoc_tresh = 100, snp_cutoff = 200, min_diff = 100, pdf_file = "", eps = 1e-3, assoc_dir = "both"){
  target <- sub("[0-9][._].+", "", names(dist_df)[6])
  if (pdf_file != "") pdf(file = pdf_file, width = 15, height = 5)
  target_haps <- grep(target, sample_list)
  plot(x= 1, y=1, xlim = c(0,max(dist_df[,4])), ylim = c(log10(eps)-1, -log10(eps)+1), type = "n", main = paste(target,"x vs Atlantic & Sea of Japan", sep = ""), xlab = "Cumulative position", ylab = "Log10 of distance ratio")
  hap_ratio_df <- data.frame(hap = sample_list[target_haps], stringsAsFactors= F)
  #assoc_tresh <- 100
  #snp_cutoff <- 200
  #min_diff <- 100
  
  for(bf_hap in 1:dim(hap_ratio_df)[1]){
    bf_ind <- ceiling(bf_hap/2)
    hap_no <- 2 - bf_hap%%2
    hap_ratio_df[bf_hap, "ref_1_mean"] <- mean(dist_df[,2*(bf_hap-1)+6])
    hap_ratio_df[bf_hap, "ref_2_mean"] <- mean(dist_df[,2*(bf_hap-1)+7])
    mean_ratio <- hap_ratio_df[bf_hap, "ref_1_mean"]/hap_ratio_df[bf_hap, "ref_2_mean"]
    
    rv <- (dist_df[,2*(bf_hap-1)+6]+eps)/(dist_df[,2*(bf_hap-1)+7]+eps)
    rv[snp_numbers <= snp_cutoff|(dist_df[,2*(bf_hap-1)+6] < min_diff & dist_df[,2*(bf_hap-1)+7] < min_diff)] <- NA #1
    
    dist_df[,paste(hap_ratio_df[bf_hap, "hap"], "_ratio", sep = "")] <- rv
    if(assoc_dir == "both"){
      assoc_filter <- which(rv > assoc_tresh * mean_ratio | rv < mean_ratio/assoc_tresh) #which(rv > assoc_tresh | rv < 1/assoc_tresh)
    }
    if(assoc_dir == "up"){
      assoc_filter <- which(rv > assoc_tresh * mean_ratio ) #which(rv > assoc_tresh | rv < 1/assoc_tresh)
    }
    if(assoc_dir == "down") {
      assoc_filter <- which(rv < mean_ratio/assoc_tresh) #which(rv > assoc_tresh | rv < 1/assoc_tresh)
    }
    hap_ratio_df[bf_hap, "pac"] <- sum(rv > mean_ratio*assoc_tresh, na.rm = T)
    hap_ratio_df[bf_hap, "atl"] <- sum(rv < mean_ratio/assoc_tresh, na.rm = T)
    bg_col_vec <- c("grey20", "grey80")[match(as.integer(as.factor(dist_df[, 1])), unique(as.integer(as.factor(dist_df[, 1]))))%%2 + 1]
    points(y = log10(dist_df[-assoc_filter,paste(hap_ratio_df[bf_hap, "hap"], "_ratio", sep = "")]), x = dist_df[-assoc_filter,4], col = bg_col_vec[-assoc_filter], pch = 16, cex=0.2)
    points(y = log10(dist_df[assoc_filter,paste(hap_ratio_df[bf_hap, "hap"], "_ratio", sep = "")]), x = dist_df[assoc_filter,4], col = rgb(t(col2rgb(bf_ind+1)), alpha = 120, maxColorValue = 255), pch = 15+hap_no, cex = 0.8)
  }
  if (pdf_file != "")dev.off()
  return(list(dist = dist_df, hap= hap_ratio_df))
}


introgression_corr_plot <- function(intro_obj, sample_list, assoc_tresh = 100, pdf_file = "", snp_set = "signif", eps = 1e-3){
  
  dist_df <- intro_obj$dist
  ratio_mean <- sum(intro_obj$hap[,"ref_1_mean"])/sum(intro_obj$hap[,"ref_2_mean"])
  #Finding windows below 1e-5
  ratio_cols <- grep("ratio", names(dist_df))
  dist_df[is.na(dist_df)] <- 1
  min_idx <- matrix(ncol=2, data = c(1:dim(dist_df)[1],max.col(-dist_df[,ratio_cols])))
  ref_1_wins <- which(dist_df[,ratio_cols][min_idx] < ratio_mean/assoc_tresh)
  #plot_top_regions_2(atl_wins, h67_w, herring_67$geno, h67_clean_samples, plot_dir="~/Projects/Herring/doc/Hap_dist_scan/67_inds/Full_set/Balsfjord_introgression/atl/", snp_col="blue", snp_lwd = 0.25, snp_alpha = 120)
  
  #Finding windows above 1e5
  max_idx <- matrix(ncol=2, data = c(1:dim(dist_df)[1],max.col(dist_df[,ratio_cols])))
  ref_2_wins<- which(dist_df[,ratio_cols][max_idx] > ratio_mean*assoc_tresh)
  #plot_top_regions_2(pac_wins, h67_w, herring_67$geno, h67_clean_samples, plot_dir="~/Projects/Herring/doc/Hap_dist_scan/67_inds/Full_set/Balsfjord_introgression/pac/", snp_col="blue", snp_lwd = 0.25, snp_alpha = 120)
  signif_win <- unique(c(ref_1_wins, ref_2_wins))
  plot_title <- "Introgression windows"	
  if (snp_set == "all"){
    signif_win <- 1:dim(dist_df)[1]
    plot_title <- "All windows"	
  }
  #Halpotype correlation plot
  if(pdf_file != "") pdf(file = pdf_file, width=10, height=10)
  
  
  plot(log10(dist_df[signif_win,22]), log10(dist_df[signif_win,24]), xlim = c(log10(eps)-3, -log10(eps)+3), ylim = c(log10(eps)-3, -log10(eps)+3), type = "n", xlab = "Haplotype 1", ylab = "Haplotype 2", main = plot_title)
  #plot(log10(hws6_dist_df[signif_win,22]), log10(hws6_dist_df[signif_win,24]), xlim = c(-3, 3), ylim = c(-3, 3), type = "n", xlab = "Haplotype 1", ylab = "Haplotype 2", main = "All windows")
  box_cutoff_up <- log10(ratio_mean*assoc_tresh)
  box_cutoff_down <- log10(ratio_mean/assoc_tresh)
  rect(box_cutoff_down,box_cutoff_down,box_cutoff_up,box_cutoff_up, col = "grey90")
  rect(c(-10, box_cutoff_up),c(-10,box_cutoff_up),c(box_cutoff_down, 10),c(box_cutoff_down, 10), col = c("steelblue", "lightskyblue"))
  rect(c(-10, box_cutoff_up),c(box_cutoff_up,-10),c(box_cutoff_down, 10),c(10, box_cutoff_down), col = "salmon")
  rect(c(-10,box_cutoff_up,box_cutoff_down, box_cutoff_down),c(box_cutoff_down,box_cutoff_down,box_cutoff_up,-10),c(box_cutoff_down,10,box_cutoff_up,box_cutoff_up),c(box_cutoff_up,box_cutoff_up,10,box_cutoff_down), col = c("mediumpurple3", "mediumpurple1")[c(1,2,2,1)])
  corr_data_collection <- matrix(nrow=3, ncol = 3, data = 0)
  
  ratio_offset <- min(ratio_cols)-1
  n_ratio_cols <- length(ratio_cols)
  for(i in 1:(n_ratio_cols-1)){
    for (j in (i+1):n_ratio_cols){
      x_vec <- log10(dist_df[signif_win,ratio_offset +i])
      y_vec <- log10(dist_df[signif_win,ratio_offset +j])
      points(x = x_vec, y = y_vec, pch = 16, cex = 0.5, col = rgb(0,0,0, 0.5))
      #points(log10(hws6_dist_df[,21+i]), log10(hws6_dist_df[,21+j]), pch = 16, cex = 0.5, col = rgb(0,0,0, 0.2))
      corr_data_collection[1,1] <- corr_data_collection[1,1] + sum(x_vec < box_cutoff_down & y_vec < box_cutoff_down)
      corr_data_collection[1,2] <- corr_data_collection[1,2] + sum(x_vec < box_cutoff_down & y_vec >= box_cutoff_down & y_vec <= box_cutoff_up)
      corr_data_collection[1,3] <- corr_data_collection[1,3] + sum(x_vec < box_cutoff_down & y_vec > box_cutoff_up)
      corr_data_collection[2,1] <- corr_data_collection[2,1] + sum(x_vec >= box_cutoff_down & x_vec <= box_cutoff_up & y_vec < box_cutoff_down)
      corr_data_collection[2,2] <- corr_data_collection[2,2] + sum(x_vec >= box_cutoff_down & x_vec <= box_cutoff_up & y_vec >= box_cutoff_down & y_vec <= box_cutoff_up)
      corr_data_collection[2,3] <- corr_data_collection[2,3] + sum(x_vec >= box_cutoff_down & x_vec <= box_cutoff_up & y_vec > box_cutoff_up)
      corr_data_collection[3,1] <- corr_data_collection[3,1] + sum(x_vec > box_cutoff_up & y_vec < box_cutoff_down)
      corr_data_collection[3,2] <- corr_data_collection[3,2] + sum(x_vec > box_cutoff_up & y_vec >= box_cutoff_down & y_vec <= box_cutoff_up)
      corr_data_collection[3,3] <- corr_data_collection[3,3] + sum(x_vec > box_cutoff_up & y_vec > box_cutoff_up)
    }
  }
  corr_labels <- paste(c("Atl/Atl", "Atl/n.s.", "Atl/Pac", "n.s./Atl", "", "n.s./Pac", "Pac/Atl", "Pac/n.s.", "Pac/Pac"), "\nn = ", corr_data_collection)
  corr_labels[5] <- ""
  text(x = rep(c(box_cutoff_down+log10(eps), 0, box_cutoff_up-log10(eps)), times = 3), y = rep(c(box_cutoff_down+log10(eps), 0, box_cutoff_up-log10(eps)), each = 3), labels = corr_labels)
  if(pdf_file != "") dev.off()
  return(list(point_dist=corr_data_collection, ref_1 = ref_1_wins, ref_2 = ref_2_wins))
}


