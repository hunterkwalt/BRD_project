#### Re-estimate virus abundances using counts instead of depth
### Hunter K. Walt 
## Final version: 27 June 2025

library(data.table)
library(tidyverse)
library(tidyHeatmap)
library(grid)
library(ggfortify)
library(cowplot)
library(GGally)
library(ggpubr)
library(combinat)
library(patchwork)
library(gggenes)
library(ComplexHeatmap)
library(UpSetR)
library(ComplexUpset)
library(ggbreak)
library(edgeR)
library(data.table)

#working directory
wd <- "/mnt/md0/brd_project"
setwd(wd)

#read in necessary files
cts_data <- "/mnt/md0/brd_project/virus_counts/full_read_counts_no_polyA.txt"
vir_to_scaffold <- "virus_to_scaf.txt"
il_samps <- "virus_coverage/depth_stats/illumina_mapped.txt"
np_samps <- "virus_coverage/depth_stats/nanopore_mapped.txt"
il_mapped <- "virus_coverage/depth_stats/illumina_nopolya_mapped.txt"
locs <- "SRA_locs.txt"
symps <- "SRA_symptoms.txt"
libs <- "library_sizes/library_sizes_all_samples.tsv"
ctl_samples <- "control_samples_to_remove.txt"
il <- "virus_coverage/depth_stats/illumina_depth_nopolyA.tsv"
np <- "virus_coverage/depth_stats/nanopore_depth_nopolyA.tsv"

#read
cts <- fread(cts_data, header = F)
ils <- read.table(file = il_samps, header = F)
ilm <- read.table(file = il_mapped, header = F)
nps <- read.table(file = np_samps, header = F)
vps <- read.table(file = vir_to_scaffold, header = F)
loc <- read.table(file = locs, header = F)
symp <- read.table(file = symps, header = F)
lib <- read.table(file = libs, header = T)
ctl <- read.table(file = ctl_samples, header = F)
#read in depth datasets
ild <- fread(file = il, sep = "\t")
npd <- fread(file = np, sep = "\t")

#add colnames
names(cts) <- c("SRA", "Virus", "Length", "Counts","Unmapped_counts")
names(ils) <- "SRA"
names(nps) <- "SRA"
names(vps) <- c("Virus_name", "Virus")
names(loc) <- c("SRA", "Location")
names(symp) <- c("SRA", "Disease_State")#remove unmapped counts rows
names(ctl) <- "SRA"
cts_v <- subset(cts, Virus != "Unmapped")
cts_v$Concordant <- cts_v$Counts - cts_v$Unmapped_counts
cts_f <- cts_v[,-c(3,4,5)]
names(cts_f)[3] <- "Counts"

#make counts matrix
cts_w <- data.frame(cts_f %>% pivot_wider(id_cols = "Virus", names_from = "SRA", values_from = Counts))
ctsv <- merge(cts_w, vps, by = "Virus")
rownames(ctsv) <- ctsv$Virus_name
cts_w <- ctsv[,-c(1,ncol(ctsv))]

#make length matrix
len <- data.frame(cts_v %>% pivot_wider(id_cols = "Virus", names_from = "SRA", values_from = Length))
rownames(len) <- len$Virus
len <- len[,-c(1)]
length_df <- data.frame(Virus = rownames(len), Length = unique(len$SRR11836760))
length_df <- merge(length_df, vps, by = "Virus")
length_df_final <- length_df[,c(3,2)] 
names(length_df_final) <- c("Virus", "Length")

#get library sizes
cts$Counts[which(cts$Unmapped_counts > 0)] <- cts$Unmapped_counts[which(cts$Unmapped_counts > 0)]
lib_size <- cts %>% group_by(SRA) %>% summarise(LibSize = sum(Counts))
lib_size <- lib_size[-which(lib_size$SRA %in% ctl$SRA), ]


#filter from count data
cts_w <- cts_w[,-which(names(cts_w) %in% ctl$SRA)]


## Normalize for library size and genome length

# Creating a DGEList object for use in edgeR.
y <- DGEList(counts = cts_w, lib.size = lib_size$LibSize, genes = length_df_final)
y <- calcNormFactors(y, method = "TMMwsp")


# Get RPKMs. TMM normalized RPKM values can be statistically compared.
rpkms <- edgeR::rpkm(y, log = F)
#rpkms <- edgeR::rpkm(y, gene.length = length_df_final$Length)

write.csv(data.frame(rpkms), file = "viral_RPKMs.tsv", quote = F, sep = '\t', row.names = T)


### rpkms Make heatmap
f_rpkms <- rpkms[,which(colSums(rpkms) > 1)] #cutoff if sample doesnt have at least an RPKM of one. Filters samples with no virus to make heatmap more visually pleasing
#RPKM values less than one are turned to zero
f_rpkms[f_rpkms < 1] <- 0
#log2 transform
f_rpkms_log <- log(f_rpkms + 1, 2)
symp_sub <- symp[which(symp$SRA %in% colnames(f_rpkms_log)),]
loc_df <- loc[which(loc$SRA %in% colnames(f_rpkms_log)),]
loc_pal <- c('Ireland'='palegreen3', 'Australia'='plum4',
             "Canada"='orange3', "USA:Kansas" = "steelblue")
trpkms_log <- t(f_rpkms_log) 
#match rows for labeling
loc_match <- loc_df[match(rownames(trpkms_log), loc_df$SRA), ] #match rows for labeling
symp_match <- symp_sub[match(rownames(trpkms_log), symp_sub$SRA), ]

# make heatmap
new_hm_rpkm <- Heatmap(trpkms_log, show_row_names = F, row_split = symp_match$Disease_State, col = c("steelblue", "white", "firebrick1"),
                      left_annotation = rowAnnotation(Location = loc_match$Location, col = list(Location = loc_pal)), name = "Log2 RPKM") %>% draw() %>% grid.grabExpr()

ggsave(plot = new_hm_rpkm, filename = "/mnt/md0/brd_project/figures/new_figures/RPKM_heatmap.pdf", width = 6, height = 6, dpi = 300)



#test for higher healthy rpkm in all viruses
virus_test <- data.frame(rpkms)
virus_test$Virus <- row.names(virus_test)
virus_long <- virus_test %>% pivot_longer(!Virus, names_to = "SRA", values_to = "rpkm")
virus_sub <- subset(virus_long, Virus != "BVDV3")
#virus_sub <- virus_long
virus_symp <- merge(virus_sub, symp, by = "SRA")
virus_final <- subset(virus_symp, rpkm >= 1)
virus_final$Log2rpkm <- log(virus_final$rpkm, 2)

#make boxplot
all_viruses_rpkm <- ggboxplot(virus_final, x = "Disease_State", y = "Log2rpkm", color = "Disease_State", palette = c("#4682B4", "firebrick1"), add = "jitter") +
  theme(legend.position = "right") + xlab("Disease State") + ylim(0,20) + facet_wrap(~Virus, ncol =2, scales = "free") + 
  stat_compare_means(aes(label = paste0("p =", ..p.format..)), label.x = 1.4, label.y = 17)

ggsave(plot = all_viruses_rpkm, width = 7.5, height = 9, 
       filename = "/mnt/md0/brd_project/figures/new_figures/rpkm_comparisons_plot.pdf", dpi = 300)


## correlations between viruses based on RPKM
cors <- ggpairs(data.frame(t(f_rpkms_log)), diag = list(continuous = "blankDiag"), ggplot2::aes(alpha = 0.5))

ggsave(plot = cors, width = 11, height = 9, 
       filename = "/mnt/md0/brd_project/figures/new_figures/correlation_matrix.pdf", dpi = 300)





######################## UPSET PLOT  ###########################
# make new dataframe for summary table. If the RPKM is greater than or equal to 1, virus is present.

rpkms[rpkms >= 1] <- TRUE
rpkms[rpkms < 1] <- FALSE
rpkms_t <-  t(rpkms)
write.csv(rpkms_t[which(rowSums(rpkms_t) == 0),1], file = 'no_virus_rpkm.csv')
write.csv(rpkms, file = 'virus_presence_absence_rpkm.csv')

m1 <- make_comb_mat(rpkms_t)
ComplexHeatmap::comb_name(m1)
#m2 <- make_comb_mat(fd_mat, mode = "union")

virus_sum <- UpSet(m1, comb_order = order(comb_size(m1)), pt_size = unit(4, "mm"),
      top_annotation = upset_top_annotation(m1, add_numbers = TRUE)) %>% draw() %>% grid.grabExpr()



ggsave(plot = virus_sum, width = 10.5, height = 4, filename = "/mnt/md0/brd_project/figures/new_figures/virus_summary_figure_rpkm.pdf")
ggsave(plot = virus_fig, width = 10, height = 10, filename = "../figures/new_figures/virus_summary_figure_filtered.pdf")


################ VIRUS DEPTH FIGURE ################
#make new colnames
names(ild) <- c("Virus","Position", ilm$V1)
names(npd) <- c("Virus","Position", nps$V1)
names(ctl) <- "SRA"

#merge dataframes
fd <- cbind(npd, ild[,3:ncol(ild)])

#change names from scaffolds to viruses
for(x in 1:length(vps$V1)){
  fd$Virus[which(fd$Virus %in% vps[x,2])] <- vps[x,1]
}


#make data long for analysis
fd_long <- fd %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")
fd_long <- fd_long[-which(fd_long$SRA %in% ctl$SRA), ]

#subset by virus
bcov <- subset(fd_long, Virus == "BCoV")
brav <- subset(fd_long, Virus == "BrAv")
brbv <- subset(fd_long, Virus == "BrBv")
beve <- subset(fd_long, Virus == "BEVE")
bhv <- subset(fd_long, Virus == "BHV1")
bnv <- subset(fd_long, Virus == "BNV")
brsv <- subset(fd_long, Virus == "BRSV")
bpiv3 <- subset(fd_long, Virus == "BPIV3")
bvdv1 <- subset(fd_long, Virus == "BVDV1")
bvdv2 <- subset(fd_long, Virus == "BVDV2")
bvdv3 <- subset(fd_long, Virus == "BVDV3")


#remove samples that do not have the specific virus
tst <- bcov %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["BCoV",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
bcov <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")

tst <- brav %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["BrAv",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
brav <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")

tst <- brbv %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["BrBv",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
brbv <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")

tst <- beve %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["BEVE",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
beve <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")

tst <- bhv %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["BHV1",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
bhv <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")

tst <- bnv %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["BNV",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
bnv <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")

tst <- brsv %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["BRSV",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
brsv <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")

tst <- bpiv3 %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["BPIV3",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
bpiv3 <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")

tst <- bvdv1 %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["BVDV1",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
bvdv1 <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")

tst <- bvdv2 %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["BVDV2",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
bvdv2 <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")

tst <- bvdv3 %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["BVDV3",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
bvdv3 <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")

# add disease state
names(symp) <- c("SRA", "Disease_State")
brav_ds <- merge(brav, symp, by = "SRA")
brbv_ds <- merge(brbv, symp, by = "SRA")
bcov_ds <- merge(bcov, symp, by = "SRA")
beve_ds <- merge(beve, symp, by = "SRA")
bhb_ds <- merge(bhv, symp, by = "SRA")
bnv_ds <- merge(bnv, symp, by = "SRA")
brsv_ds <- merge(brsv, symp, by = "SRA")
bpiv_ds <- merge(bpiv3, symp, by = "SRA")
bvdv1_ds <- merge(bvdv1, symp, by = "SRA")
bvdv2_ds <- merge(bvdv2, symp, by = "SRA")
bvdv3_ds <- merge(bvdv3, symp, by = "SRA")


#get mean coverage healthy vs brd
brav_mean <- brav_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))
brbv_mean <- brbv_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))
beve_mean <- beve_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))
bnv_mean <- bnv_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))
bhv_mean <- bhb_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))
brsv_mean <- brsv_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))
bpiv_mean <- bpiv_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))
bvdv1_mean <- bvdv1_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))
bvdv2_mean <- bvdv2_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))
bvdv3_mean <- bvdv3_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))
bcov_mean <- bcov_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))


# Generate plots
options(scipen = 999)
bcov_plot <- ggplot(bcov_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  #geom_ribbon(aes(x=Position, ymax = high, ymin = low, color = Disease_State), alpha = 0.3) +
  ylab("Average Depth") +
  ggtitle("BCoV Coverage")+
  #ylim(2,12) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11))


brav_plot <- ggplot(brav_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Average Depth") +
  ggtitle("BrAv Coverage")+
  #ylim(2,12) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11))

brbv_plot <- ggplot(brbv_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Average Depth") +
  ggtitle("BrBv Coverage")+
  #ylim(0,300) +
  theme_classic() + #xlim(0,7500) + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11))

beve_plot <- ggplot(beve_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Average Depth") +
  ggtitle("BEVE Coverage")+
  #ylim(2,12) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11))


bvdv1_plot <- ggplot(bvdv1_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Average Depth") +
  ggtitle("BVDV1 Coverage")+
  #ylim(2,12) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11))


bvdv2_plot <- ggplot(bvdv2_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Average Depth") +
  ggtitle("BVDV2 Coverage")+
  #ylim(2,12) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), legend.position = "none")

bvdv3_plot <- ggplot(bvdv3_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Average Depth") +
  ggtitle("BVDV3 Coverage")+
  #ylim(2,12) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), legend.position = "none", text = element_text(size = 11))


bhv_plot <- ggplot(bhv_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Average Depth") +
  ggtitle("BHV1 Coverage")+
  #ylim(2,12) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11))

brsv_plot <- ggplot(brsv_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Average Depth") +
  ggtitle("BRSV Coverage")+
  #ylim(2,12) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11))

bpiv_plot <- ggplot(bpiv_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Average Depth") +
  ggtitle("BPIV3 Coverage")+
  #ylim(2,12) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11))

bnv_plot <- ggplot(bnv_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Average Depth") +
  ggtitle("BNV Coverage")+
  #ylim(2,12) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11))

cov_fig <- brsv_plot +  bpiv_plot + bhv_plot + bcov_plot + brav_plot + brbv_plot + bvdv1_plot + bvdv2_plot + bvdv3_plot + beve_plot + bnv_plot + 
  patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(tag_levels = "A")

cov_fig
#cov_fig <- brsv_plot +  bpiv_plot + bhv_plot + bcov_plot + brav_plot + brbv_plot + bvdv1_plot + bvdv2_plot  + beve_plot + bnv_plot + 
#  patchwork::plot_layout(guides = "collect", ncol = 2) + patchwork::plot_annotation(tag_levels = "A")


ggsave(filename = "figures/new_figures/coverage_plots_rpkm_filtered_samples.pdf", plot = cov_fig, height = 6, width = 12)


### VIRAL STRUCTURES ###
str_vir <- "/mnt/md0/brd_project/virus_coords.tsv"
str <- read.table(str_vir, sep = "\t")
names(str) <- c("Virus", "Chr", "Start", "End", "Strand", "Gene")
str$Strand <- gsub("\\-", FALSE, str$Strand)
str$Strand <- gsub("\\+", TRUE, str$Strand)

bcstr <- subset(str, Virus == "BcoV")
bevstr <- subset(str, Virus == "BEVE")
bhvstr <- subset(str, Virus == "BHV1")
bnvstr <- subset(str, Virus == "BNV")
bpivstr <- subset(str, Virus == "BPIV3")
bravstr <- subset(str, Virus == "BrAv")
brbvstr <- subset(str, Virus == "BrBv")
brsvstr <- subset(str, Virus == "BRSV")
bvdv1str <- subset(str, Virus == "BVDV1")
bvdv2str <- subset(str, Virus == "BVDV2")
bvdv3str <- subset(str, Virus == "BVDV3")

cov_plot <- ggplot(bcstr, aes(xmin = Start, xmax = End, y = Virus, 
                              forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_manual(values = "darkgrey") +
  theme_genes() + 
  theme(legend.position = "none")

beve_plot <- ggplot(bevstr, aes(xmin = Start, xmax = End, y = Virus, 
                                forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() + 
  theme(legend.position = "none")

bhv_plot <- ggplot(bhvstr, aes(xmin = Start, xmax = End, y = Virus, 
                               forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() + 
  theme(legend.position = "none")

bnv_plot <- ggplot(bnvstr, aes(xmin = Start, xmax = End, y = Virus, 
                               forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Reds") +
  theme_genes() + 
  theme(legend.position = "none")

bpiv_plot <- ggplot(bpivstr, aes(xmin = Start, xmax = End, y = Virus,  
                                 forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() + 
  theme(legend.position = "none")

brav_plot <- ggplot(bravstr, aes(xmin = Start, xmax = End, y = Virus,
                                 forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() + 
  theme(legend.position = "none")

brbv_plot <- ggplot(brbvstr, aes(xmin = Start, xmax = End, y = Virus, 
                                 forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() + 
  theme(legend.position = "none")

brsv_plot <- ggplot(brsvstr, aes(xmin = Start, xmax = End, y = Virus, 
                                 forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() + 
  theme(legend.position = "none")

bvdv1_plot <- ggplot(bvdv1str, aes(xmin = Start, xmax = End, y = Virus, 
                                   forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() + 
  theme(legend.position = "none")

bvdv2_plot <- ggplot(bvdv2str, aes(xmin = Start, xmax = End, y = Virus,
                                   forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() + 
  theme(legend.position = "none")

bvdv3_plot <- ggplot(bvdv3str, aes(xmin = Start, xmax = End, y = Virus, 
                                   forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() + 
  theme(legend.position = "none")

str_fig <- brsv_plot +  bpiv_plot + bhv_plot + cov_plot + brav_plot + brbv_plot + bvdv1_plot + bvdv2_plot + bvdv3_plot + beve_plot + bnv_plot + 
  patchwork::plot_layout(guides = "collect") + patchwork::plot_annotation(tag_levels = "A")
# ggsave(filename = "../figures/structure_plots.pdf", plot = str_fig, height = 10, width = 9)


svg( "figures/new_figures/structure_plots.svg", height = 6, width = 12)
str_fig
dev.off()
