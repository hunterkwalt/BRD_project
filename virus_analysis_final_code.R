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
library(tximport)
library(tximportData)
library(scales)
library(rstatix)
library(tidygraph)
library(ggraph)
library(latex2exp)

#working directory
wd <- "/mnt/md0/brd_project"
setwd(wd)

#read in necessary files
cts_data <- "/mnt/md0/brd_project/virus_counts/full_final_read_counts.tsv"
vir_to_scaffold <- "virus_to_scaf.txt"
il_samps <- "virus_coverage/depth_stats/illumina_mapped.txt"
np_samps <- "virus_coverage/depth_stats/nanopore_mapped.txt"
il_mapped <- "virus_coverage/depth_stats/illumina_nopolya_mapped.txt"
locs <- "SRA_locs.txt"
symps <- "SRA_symptoms.txt"
libs <- "library_sizes/library_sizes_all_samples.tsv"
ctl_samples <- "control_samples_to_remove.txt"
il <- "/mnt/md0/brd_project/virus_coverage/depth_stats/illumina_final_depth.tsv"
np <- "/mnt/md0/brd_project/virus_coverage/depth_stats/nanopore_final_depth.tsv"


#sequencing type data
nano <- "/mnt/md0/brd_project/virus_coverage/depth_stats/nanopore_mapped_final.txt"
illum <- "/mnt/md0/brd_project/virus_coverage/depth_stats/illumina_mapped_final.txt"
nan <- read.table(nano, header = F)
ill <- read.table(illum, header = F)
colnames(nan) <- "SRA"
colnames(ill) <- "SRA"
nan$Seqtype <- "Nanopore"
ill$Seqtype <- "Illumina"
st <- rbind(nan, ill)

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

# metadata
int <- merge(st, symp, by.x = "SRA", by.y = "V1")
md <- merge(int,loc, by.x = "SRA", by.y = "V1")
names(md) <- c("SRA", "Sequencing_type", "Disease_State", "Location")
final <- md[-which(md$SRA %in% ctl$V1), ]
write.csv(md, quote = F, file = "/mnt/md0/brd_project/Supplementary_Table1.csv")

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
ctsv <- ctsv[,-c(1)]
ctsv <- aggregate(.~Virus_name,data=ctsv,FUN=sum)
rownames(ctsv) <- ctsv$Virus_name
cts_w <- ctsv[,-c(1)]

#make length matrix
len <- data.frame(cts_v %>% pivot_wider(id_cols = "Virus", names_from = "SRA", values_from = Length))
rownames(len) <- len$Virus
len <- len[,-c(1)]
length_df <- data.frame(Virus = rownames(len), Length = unique(len$SRR11836760))
length_df <- merge(length_df, vps, by = "Virus")
length_df_final <- length_df[,c(3,2)] 
names(length_df_final) <- c("Virus", "Length")
length_df_final <- aggregate(.~Virus,data=length_df_final,FUN=sum)

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
rpkms <- rpkms[-which(rownames(rpkms) == "BRV"), ]

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
virus_symp$Log2rpkm <- log(virus_symp$rpkm +1, 2)

###### make Figure 2 #########


all_viruses_rpkm <- ggboxplot(virus_symp, x = "Disease_State", y = "Log2rpkm", color = "Disease_State", palette = c("#4682B4", "firebrick1"), add = "jitter") +
  theme(legend.position = "right", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) + ylim(0,19) + labs(color = "Disease State") +
  facet_wrap(~Virus, ncol =2, scales = "free") + stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") +
  stat_compare_means(aes(label = paste0("p =", ..p.format..)), label.x = 1.4, label.y = 16, method = "wilcox.test", p.adjust.method = "bonferroni")

ggsave(plot = all_viruses_rpkm, width = 7.5, height = 9, 
       filename = "/mnt/md0/brd_project/figures/new_figures/final_figures/rpkm_comparisons_plot.pdf", dpi = 300)


######### Full viral abundance ############
va <- virus_symp %>% group_by(SRA, Disease_State) %>% summarise(Viral_abundance = sum(Log2rpkm))
vast <- merge(va, st, by = "SRA")
vastl <- merge(vast, loc_df, by = "SRA")

#compare viral abundance by disease state
viral_abundance <- ggboxplot(vastl, x = "Disease_State", y = "Viral_abundance", color = "Disease_State", palette = c("#4682B4", "firebrick1"), add = "jitter") +
  theme(legend.position = "right", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
        plot.title = element_text(size = 11), text = element_text(size = 11)) + 
  ggtitle("Viral Abundance by Disease State") +
  stat_compare_means(aes(label = paste0("p =", ..p.format.., "**"), ), method = "wilcox.test", label.x = 1.3, label.y = 37) + ylab("Total Viral Abundance") +
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") + labs(color = "Disease State")

#compare viral abundance by location
mycomps <- list(c("Canada", "Australia"), c("Ireland", "Australia"), c("Ireland", "Canada"))
viral_abundanceloc <- ggboxplot(vastl, x = "Location", y = "Viral_abundance", color = "Location", palette = loc_pal, add = "jitter") +
  theme(legend.position = "right", axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title = element_blank(), plot.title = element_text(size = 11), text = element_text(size = 11)) + 
  ggtitle("Viral Abundance by Location") +
  stat_compare_means(comparisons = mycomps, method = "wilcox.test") + ylab("Total Viral Abundance") +
  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") + labs(color = "Location")


viral_abundance
viral_abundanceloc


######################## ONEHOT Encode ###########################
# make new dataframe for summary table. If the RPKM is greater than or equal to 1, virus is present.

rpkms[rpkms >= 1] <- TRUE
rpkms[rpkms < 1] <- FALSE
rpkms_t <-  t(rpkms)
onehot <- rpkms_t %>% data.frame() %>% mutate("SRA" = rownames(rpkms_t)) %>% pivot_longer(!SRA, values_to = "Presence/Absence")
names(onehot)[2] <- "Virus"
merge_counts <- merge(data.frame(rpkms_t), md, by.x = "row.names", by.y = "SRA")

bhv_bvdv <- merge_counts[,c(4, 11, 17)][which(rowSums(merge_counts[,c(4, 11)]) == 2), ] 
names(merge_counts)[1] <- "SRA"

#load in bacterial data
bac_pa <- "/mnt/md0/brd_project/bacteria_presence_absence.csv"
bac_data <- read.csv(bac_pa)
all_counts <- merge(merge_counts, bac_data, by = "SRA")
#format to long
bac_long <-bac_data %>% pivot_longer(!SRA, names_to = "Virus", values_to = "Presence/Absence")

#add bacteria data to virus data
virus_bac <- rbind(onehot, bac_long)
healthy_samps <- merge_counts$SRA[which(merge_counts$Disease_State == "Healthy")]
brd_samps <- merge_counts$SRA[which(merge_counts$Disease_State == "BRD")]

#split between healthy and diseased.
healthy_onehot <- virus_bac[which(virus_bac$SRA %in% healthy_samps), ]
brd_onehot <- virus_bac[which(virus_bac$SRA %in% brd_samps), ]

write.csv(rpkms_t[which(rowSums(rpkms_t) == 0),1], file = 'no_virus_rpkm.csv')
write.csv(rpkms, file = 'virus_presence_absence_rpkm.csv')
write.csv(onehot, file = 'virus_presence_absence_onehot.csv')
write.csv(virus_bac, file = '/mnt/md0/brd_project/co_occur_with_bacteria_onehot.csv')
write.csv(healthy_onehot, file = '/mnt/md0/brd_project/healthy_onehot.csv')
write.csv(brd_onehot, file = '/mnt/md0/brd_project/brd_onehot.csv')





########## Get viral diversity (i.e. number of viral taxa in each sample) ##########
vir_md <- merge(rpkms_t, loc, by.x = "row.names", by.y = "SRA")
vir_mdf <- merge(vir_md, symp, by.x = "Row.names", by.y = "SRA")
viral_counts <- vir_mdf %>% group_by(Location, Disease_State) %>% summarise(across(where(is.numeric), list(sum = sum)))
names(viral_counts) <- gsub("_sum", "", names(viral_counts))

obs <- data.frame(Location = vir_mdf$Location,
                  Disease_State = vir_mdf$Disease_State,
                  rowSums(x = vir_mdf[,-c(1,16,17)]))
names(obs)[3] <- "Observed_Diversity"

#plot by disease state
virdiv <- ggboxplot(obs, x = "Disease_State", y = "Observed_Diversity", color = "Disease_State", palette = c("#4682B4", "firebrick1"), add = "point",
                    add.params = list(position = position_jitter(width =  0.22, height = 0.01))) +
  theme(legend.position = "right", plot.title = element_text(size = 11), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        text = element_text(size = 11)) +
  ggtitle("Viral Diversity by Disease State") +  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") +
  stat_compare_means(aes(label = paste0("p =", ..p.format..), ), method = "wilcox.test", label.x = 1.3) + ylab("# Viruses") + labs(color = "Disease State")


#plot by location
mycomps <- list(c("Canada", "Australia"), c("Ireland", "Australia"), c("Ireland", "Canada"))
virdivloc <- ggboxplot(obs, x = "Location", y = "Observed_Diversity", color = "Location", palette = loc_pal, add = "point",
                       add.params = list(position = position_jitter(width =  0.22, height = 0.01))) +
  theme(legend.position = "right", plot.title = element_text(size = 11), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title = element_blank(), text = element_text(size = 11)) + scale_y_continuous(breaks = c(0,2,4,6,8)) +
  ggtitle("Viral Diversity by Location") +  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") +
  stat_compare_means(comparisons = mycomps, method = "wilcox.test", label.y = c(7,8,9)) + ylab("# Viruses") + labs(color = "Location")


#make full plot for figure 1
design  <- "AABC
            AADE"
hmw <- wrap_plots(new_hm_rpkm, viral_abundance, viral_abundanceloc, virdiv, virdivloc, design = design, guides = "collect") + plot_annotation(tag_levels = "A")   
hmw

ggsave(plot = hmw, filename = "/mnt/md0/brd_project/figures/new_figures/final_figures/viral_figure_final_withloc.pdf", width = 13.5, height = 7.5)



################ VIRUS DEPTH FIGURE (Figure 3) ################
#make new colnames
names(ild) <- c("Virus","Position", ilm$V1)
names(npd) <- c("Virus","Position", nps$SRA)
names(ctl) <- "SRA"

#merge dataframes
fd <- cbind(npd, ild[,3:ncol(ild)])

vseg <- "/mnt/md0/brd_project/virus_to_scaf_idv_segments.txt"
vps <- read.table(vseg)
names(vps) <- c("Virus_name", "Contig")

#change names from scaffolds to viruses
for(x in 1:length(vps$Virus_name)){
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
idv <- subset(fd_long, Virus == "IDV")
bpv <- subset(fd_long, Virus == "BPV")
bhkv <- subset(fd_long, Virus == "BHKV")
idv <- fd_long[grep("IDV",fd_long$Virus),]
idv <- idv |>
  separate_wider_delim(Virus, delim = "_", names = c("Virus", "Segment"))
idv <- idv[order(idv$SRA),]
idv_positions <- seq(1,nrow(idv))
idv$Position <- idv_positions

#get segements to read straight through instead of starting over
segtst <- subset(idv, SRA == "SRR11836760")
strait <- rep(seq(1,nrow(segtst)), length(unique(idv$SRA)))
idv$Position <- strait

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

tst <- idv %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["IDV",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
idv <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")

tst <- bpv %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["BPV",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
bpv <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")

tst <- bhkv %>% pivot_wider(values_from = Coverage, names_from = SRA)
vir_oc <- names(which(rpkms["BHKV",] == 1))
tst3 <- tst[,c("Virus","Position",vir_oc)]
bhkv <- tst3 %>% pivot_longer(!c(Virus,Position), names_to = c("SRA"), values_to = "Coverage")


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
idv_ds <- merge(idv, symp, by = "SRA")
bpv_ds <- merge(bpv, symp, by = "SRA")
bhkv_ds <- merge(bhkv, symp, by = "SRA")


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
bhkv_mean <- bhkv_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))
bpv_mean <- bpv_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))
idv_mean <- idv_ds %>% group_by(Disease_State, Position) %>%
  summarise_at(vars(Coverage), list(Mean_Coverage = mean))

# Generate plots
options(scipen = 999)
bcov_plot <- ggplot(bcov_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  #geom_ribbon(aes(x=Position, ymax = high, ymin = low, color = Disease_State), alpha = 0.3) +
  ylab("Avg. Depth") +
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3)) + # kb
  ggtitle("BCoV Coverage")+
  #ylim(2,12) +
  theme_classic() + labs(color = "Disease State") + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), legend.title = element_text(size = 14), legend.text = element_text(size = 12))


brav_plot <- ggplot(brav_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Avg. Depth") +
  ggtitle("BrAv Coverage")+
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3)) + # kb
  #ylim(2,12) +
  theme_classic() + labs(color = "Disease State") + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), axis.title.y = element_blank(), legend.title = element_text(size = 14), legend.text = element_text(size = 12))

brbv_plot <- ggplot(brbv_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Avg. Depth") +
  ggtitle("BrBv Coverage")+
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3)) + # kb
  #ylim(0,300) +
  theme_classic() + labs(color = "Disease State") +  #xlim(0,7500) + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), axis.title.y = element_blank(), legend.title = element_text(size = 14), legend.text = element_text(size = 12))

beve_plot <- ggplot(beve_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Avg. Depth") +
  ggtitle("BEVE Coverage")+
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3)) + # kb
  #ylim(2,12) +
  theme_classic() + labs(color = "Disease State") + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), legend.title = element_text(size = 14), legend.text = element_text(size = 12))


bvdv1_plot <- ggplot(bvdv1_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Avg. Depth") +
  ggtitle("BVDV1 Coverage")+
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3)) + # kb
  #ylim(2,12) +
  theme_classic() + labs(color = "Disease State") + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), legend.title = element_text(size = 14), legend.text = element_text(size = 12))


bvdv2_plot <- ggplot(bvdv2_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line() +
  ylab("Avg. Depth") +
  ggtitle("BVDV2 Coverage")+
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3)) + # kb
  #ylim(2,12) +
  theme_classic() + labs(color = "Disease State") + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), legend.position = "none", axis.title.y = element_blank())

bvdv3_plot <- ggplot(bvdv3_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Avg. Depth") +
  ggtitle("BVDV3 Coverage")+
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3)) + # kb
  #ylim(2,12) +
  theme_classic() + labs(color = "Disease State") + 
  theme(axis.title.x = element_blank(), legend.position = "none", text = element_text(size = 11), axis.title.y = element_blank(), legend.title = element_text(size = 14), legend.text = element_text(size = 12))


bhv_plot <- ggplot(bhv_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Avg. Depth") +
  ggtitle("BHV1 Coverage")+
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3)) + # kb
  #ylim(2,12) +
  theme_classic() + labs(color = "Disease State") + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), axis.title.y = element_blank(), legend.title = element_text(size = 14), legend.text = element_text(size = 12))

brsv_plot <- ggplot(brsv_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Avg. Depth") +
  ggtitle("BRSV Coverage")+
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3)) + # kb
  #ylim(2,12) +
  theme_classic() + labs(color = "Disease State") + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), legend.title = element_text(size = 14), legend.text = element_text(size = 12))

bpiv_plot <- ggplot(bpiv_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Avg. Depth") +
  ggtitle("BPIV3 Coverage")+
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3)) + # kb
  #ylim(2,12) +
  theme_classic() + labs(color = "Disease State") + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), axis.title.y = element_blank(), legend.title = element_text(size = 14), legend.text = element_text(size = 12))

bnv_plot <- ggplot(bnv_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Avg. Depth") +
  ggtitle("BNV Coverage")+
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3)) + # kb
  #ylim(2,12) +
  theme_classic() + labs(color = "Disease State") + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), axis.title.y = element_blank(), legend.title = element_text(size = 14), legend.text = element_text(size = 12))

bhkv_plot <- ggplot(bhkv_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("") +
  ggtitle("BHKV Coverage")+
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3)) + # kb
  #ylim(2,12) +
  theme_classic() + labs(color = "Disease State") + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), legend.title = element_text(size = 14), legend.text = element_text(size = 12))


bpv_plot <- ggplot(bpv_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State)) + 
  geom_line()+
  ylab("Avg. Depth") +
  ggtitle("BPV Coverage")+
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3)) + # kb
  #ylim(2,12) +
  theme_classic() + labs(color = "Disease State") + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), legend.title = element_text(size = 14), legend.text = element_text(size = 12))


idv_plot <- ggplot(idv_mean, aes(x=Position, y=Mean_Coverage, color = Disease_State), legend.title = element_text(size = 14), legend.text = element_text(size = 12)) + 
  geom_line()+
  ylab("Avg. Depth") +
  ggtitle("IDV Coverage")+
  scale_x_continuous(labels = label_number(suffix = "Kb", scale = 1e-3, accuracy = 1)) + # kb
  #ylim(2,12) +
  theme_classic() + labs(color = "Disease State") + 
  theme(axis.title.x = element_blank(), text = element_text(size = 11), axis.title.y = element_blank(), legend.title = element_text(size = 14), legend.text = element_text(size = 12))

# design <- "
# 123
# 456
# 789
# 101112
# 131415"

cov_fig <- brsv_plot +  bpiv_plot + bhv_plot + bcov_plot + brav_plot + brbv_plot + bvdv1_plot + bvdv2_plot + bvdv3_plot + beve_plot + bnv_plot + bhkv_plot + bpv_plot + idv_plot +
  guide_area() + patchwork::plot_layout(guides = "collect", ncol = 3) + patchwork::plot_annotation(tag_levels = "A")

cov_fig
#cov_fig <- brsv_plot +  bpiv_plot + bhv_plot + bcov_plot + brav_plot + brbv_plot + bvdv1_plot + bvdv2_plot  + beve_plot + bnv_plot + 
#  patchwork::plot_layout(guides = "collect", ncol = 2) + patchwork::plot_annotation(tag_levels = "A")


ggsave(filename = "figures/new_figures/coverage_plots_rpkm_filtered_samples_added.pdf", plot = cov_fig, height = 7, width = 8)


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
idvstr <- subset(str, Virus == "IDV")
bhkvstr <- subset(str, Virus == "BHKV")
bpvstr <- subset(str, Virus == "BPV")


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

idv_plot <- ggplot(idvstr, aes(xmin = Start, xmax = End, y = Virus, 
                               forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() + 
  theme(legend.position = "none")

bhkv_plot <- ggplot(bhkvstr, aes(xmin = Start, xmax = End, y = Virus, 
                                 forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() + 
  theme(legend.position = "none")

bpv_plot <- ggplot(bpvstr, aes(xmin = Start, xmax = End, y = Virus, 
                               forward = Strand)) +
  geom_gene_arrow(fill = "skyblue", alpha= 0.9) +
  facet_wrap(~Virus, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  theme_genes() + 
  theme(legend.position = "none")

str_fig <- brsv_plot +  bpiv_plot + bhv_plot + cov_plot + brav_plot + brbv_plot + bvdv1_plot + bvdv2_plot + bvdv3_plot + beve_plot + bnv_plot +  bhkv_plot + bpv_plot + idv_plot +
  guide_area() + patchwork::plot_layout(guides = "collect", ncol = 3) + patchwork::plot_annotation(tag_levels = "A")
# ggsave(filename = "../figures/structure_plots.pdf", plot = str_fig, height = 10, width = 9)


svg( "figures/new_figures/structure_plots.svg", height = 7, width = 8)
str_fig
dev.off()




###### ZHANG'S Metric - figure 5 This requires the .csv file output from calculating zhang's metric ########
z<- "/mnt/md0/brd_project/association/zhangs_w_bact.csv"
zm <- read.csv(z)
names(zm) <- c("Virus1", "Virus2", "Zhangs_A_to_B","Zhangs_B_to_A", 
               "SupportA", "SupportB", "Union") 

get_lower_tri<-function(cormat1){
  cormat1[upper.tri(cormat1)] <- NA
  return(cormat1)}


# zml <- zm %>% pivot_longer(!Virus1, names_to = "Virus2")
# zml$value[which(zml$value == "Void")] <- NA
# zml$value <- as.numeric(zml$value)
zmt<-zm
#zmt <- subset(zm, SupportA > 0.01 & SupportB > 0.01)
zmt[grepl("_", zmt)]
coal <- subset(zmt, Virus1 == "BHV1_BPIV3_BRSV" | Virus2 == "BHV1_BPIV3_BRSV")

zmtf <- filter(zmt, !grepl("_",Virus1))
zmtf <- filter(zmt, !grepl("_",Virus2))
# zmtf$Virus1 <- factor(zmtf$Virus1, levels = c("BRSV", "BPIV3", "BHV1", "BEVE", "BrAv", "BcoV", "BrBv", "BNV", "BVDV1", "BVDV2", "BVDV3"))
# zmtf$Virus2 <- factor(zmtf$Virus2, levels = c("BRSV", "BPIV3", "BHV1", "BEVE", "BrAv", "BcoV", "BrBv", "BNV", "BVDV1", "BVDV2", "BVDV3"))


melted_cormat1 <- melt(zmtf, na.rm = TRUE)
melted_cormat1 <- subset(melted_cormat1, variable == "Zhangs_A_to_B")
melted_cormat2 <- melted_cormat1[which(melted_cormat1$Virus1 != melted_cormat1$Virus2),]
ttt <- melted_cormat2 %>% pivot_wider(names_from = Virus1, values_from = value)
tttt <- as.matrix(ttt[,3:18])
rownames(tttt) <- ttt$Virus2
ab <- melt(tttt)
names(ab) <- c("Virus1", "Virus2", "value")


melted_cormat1 <- melt(zmtf, na.rm = TRUE)
melted_cormat1 <- subset(melted_cormat1, variable == "Zhangs_B_to_A")
melted_cormat2 <- melted_cormat1[which(melted_cormat1$Virus1 != melted_cormat1$Virus2),]
ttt <- melted_cormat2 %>% pivot_wider(names_from = Virus1, values_from = value)
tttt <- as.matrix(ttt[,3:18])
rownames(tttt) <- ttt$Virus2
full_t <- rcompanion::fullPTable(tttt)
t(full_t)
ba <- melt(tttt)
names(ba) <- c("Virus2", "Virus1", "value")
test <- na.omit(rbind(ab,ba[,c(2,1,3)]))
test$Virus1 <- factor(test$Virus1, levels = c("BRSV", "BPIV3", "BHV1", "BEVE", "BrAv", "BCoV", "BrBv", "BNV", "BVDV1", "IDV", "BHKV", "BPV", "Mannheimia", "Pasteurella", "Histophilus", "Mycoplasma", "Moraxella"))
test$Virus2 <- factor(test$Virus2, levels = c("BRSV", "BPIV3", "BHV1", "BEVE", "BrAv", "BCoV", "BrBv", "BNV", "BVDV1", "IDV", "BHKV", "BPV", "Mannheimia", "Pasteurella", "Histophilus", "Mycoplasma", "Moraxella"))

zhang_table1 <- na.omit(test) %>%
  ggplot(aes(Virus2, Virus1, fill = value)) +
  geom_hline(yintercept = 1:50, color = "darkgray", linetype = "dotted") +
  #geom_vline(xintercept = 1:50, color = "lightgray", alpha = 0.7, linetype = "dotted") +
  geom_abline(color = "darkgray", linetype = "dotted") +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "steelblue1", high = "red", mid = "white", 
                       na.value = "grey", midpoint = 0, limit = c(-1,1.000000001), 
                       space = "Lab", name = NULL) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  geom_text(aes(Virus2, Virus1, label = round(value, 2)), color = "black", size = 2.5) +
  labs(y = "Microbe A", x = "Microbe B", title = TeX("Zhang's Association Metric: Microbe A$\\rightarrow$B")) +
  theme_bw() +
  #geom_vline(xintercept = 6.5, color = "blue", linetype = "dashed") +
  #geom_vline(xintercept = 12.5, color = "blue", linetype = "dashed") +
  #geom_vline(xintercept = 19.5, color = "blue", linetype = "dashed") +
  #geom_hline(yintercept = 6.5, color = "blue", linetype = "dashed") +
  #geom_hline(yintercept = 12.5, color = "blue", linetype = "dashed") +
  #geom_hline(yintercept = 19.5, color = "blue", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1, 
                                   face = "italic", color = "black"),
        axis.text.y = element_text(size = 11, face = "italic", 
                                   color = "black"),
        strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        legend.text = element_text(size = 13),
        plot.title = element_text(hjust = 0.5, size = 13),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right") 
#coord_fixed()

zhang_table2 <- na.omit(test) %>%
  ggplot(aes(Virus2, Virus1, fill = value)) +
  #geom_hline(yintercept = 1:50, color = "lightgray", alpha = 0.7, linetype = "dotted") +
  geom_vline(xintercept = 1:50, color = "darkgray", linetype = "dotted") +
  geom_abline(color = "darkgray", linetype = "dotted") +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "steelblue1", high = "red", mid = "white", 
                       na.value = "grey", midpoint = 0, limit = c(-1,1.000000001), 
                       space = "Lab", name = NULL) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  geom_text(aes(Virus2, Virus1, label = round(value, 2)), color = "black", size = 2.5) +
  labs(y = "Microbe A", x = "Microbe B", title = TeX("Zhang's Association Metric: Microbe B$\\rightarrow$A")) +
  theme_bw() +
  #geom_vline(xintercept = 6.5, color = "blue", linetype = "dashed") +
  #geom_vline(xintercept = 12.5, color = "blue", linetype = "dashed") +
  #geom_vline(xintercept = 19.5, color = "blue", linetype = "dashed") +
  #geom_hline(yintercept = 6.5, color = "blue", linetype = "dashed") +
  #geom_hline(yintercept = 12.5, color = "blue", linetype = "dashed") +
  #geom_hline(yintercept = 19.5, color = "blue", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1, 
                                   face = "italic", color = "black"),
        axis.text.y = element_text(size = 11, face = "italic", 
                                   color = "black"),
        strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        legend.text = element_text(size = 13),
        plot.title = element_text(hjust = 0.5, size = 13),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right") 
#coord_fixed()


t2 <- zhang_table1 + zhang_table2  + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")


ggsave("/mnt/md0/brd_project/figures/new_figures/final_figures/zhang_association_table_with_microbes.pdf", 
       plot = t2,
       height = 5, width = 13.5)


