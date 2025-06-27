### Analysis of BRD-associated microbiome results from Kraken/Bracken
## Hunter K. Walt 6/11/2025

library(tidyverse)
library(phyloseq)
library(vegan)
library(RColorBrewer)
library(ggplot2)
library(stats)
library(rstatix)
library(ggpubr)

#working directory
wd <- "/mnt/md0/brd_project/bacteria"
setwd(wd)

#read in metadata
locs <- "/mnt/md0/brd_project/SRA_locs.txt"
symps <- "/mnt/md0/brd_project/SRA_symptoms.txt"
bm <- "BRD_bracken_genera.biom" 
ctl_samples <- "/mnt/md0/brd_project/control_samples_to_remove.txt"
novir <- "/mnt/md0/brd_project/no_virus_rpkm.csv"


################ VIRUS SAMPLES ########################

#import biom file from kraken-biom
data <-import_biom(bm, parseFunction=parse_taxonomy_default)
data@tax_table@.Data <- substring(data@tax_table@.Data, 4)
#rename columns to species level
colnames(data@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#View(data@tax_table@.Data)

#get rid of contaminant human reads (Only for genus and species)
#I added a taxa to this that is likely a false positive from Kraken
data <- subset_taxa(data, (!Genus %in% "Homo"))
data <- subset_taxa(data, (!Genus %in% "Callorhinchus"))
data <- subset_taxa(data, (!Genus %in% "Bos"))


#load in sample metadata
loc_df <- read.table(file = locs, header = F, sep = "\t")
symp_df <- read.table(file = symps, header = F, sep = "\t")
ctl_df <- read.table(file = ctl_samples, header = F, sep = "\t")
colnames(loc_df) <- c("SRA", "Location")
colnames(symp_df) <- c("SRA", "Disease_State")
md <- merge(loc_df, symp_df, by="SRA")
row.names(md) <- md$SRA
nv <- read.csv(novir)
novirus <- nv$X



#change names of samples from files
#change this based on what taxonomic level you imported from Bracken. e.g., if you imported genus level data,
##change "_._report_bracken_families" to "_._report_bracken_genuses" 
colnames(data@otu_table@.Data) <- gsub("\\.k2_bracken_genuses","", colnames(data@otu_table@.Data))

#create sample data for phyloseq object
samdf <- md[which(md$SRA %in% colnames(otu_table(data))), ]

#make sample names the rows
data@sam_data<-sample_data(samdf)
virus <- colnames(otu_table(data))[which(!colnames(otu_table(data)) %in% novirus)]
#filter out samples with less than 5 OTUs detected.
data.alpha_f <- estimate_richness(data)
#data <- prune_samples(row.names(data.alpha_f[which(data.alpha_f$Observed >= 10),]) , data)
data <- prune_samples(colnames(otu_table(data))[which(!colnames(otu_table(data)) %in% ctl_df$V1)], data)
data <- prune_samples(virus, data)
data.alpha <- estimate_richness(data)

#data <- prune_taxa(taxa_sums(data) != 0, data) 

######### alpha diveristy ###########
brd_pal <- c("Healthy" = "#4682B4", "BRD" = "firebrick1")

shannon <- data.alpha[,"Shannon"]
shannon_div <- data.frame(row.names = rownames(data.alpha),shannon, sample_data(data)$Disease_State, sample_data(data)$Location)
colnames(shannon_div) <- c("Shannon", "Disease_State","Location")

#plot boxes
b <- ggplot(shannon_div, aes(x=Disease_State,y=Shannon, fill=Disease_State)) + geom_boxplot() 
alpha_div <- b + ggtitle("Alpha Diversity") + ylab("Shannon Diversity Index")+ scale_fill_manual(values = brd_pal) + theme_classic2()
alpha_div

#add stats
pairwise <- shannon_div %>% wilcox_test(Shannon ~ Disease_State) %>% adjust_pvalue(method = "bonferroni") %>% add_significance()
pairwise <- pairwise %>% add_xy_position(x = "Disease_State")

b_l <- ggboxplot(shannon_div, x="Disease_State",y="Shannon", color="Disease_State", add = "jitter", palette = brd_pal)
alpha_div_l <- b_l + ggtitle("Alpha Diversity") + ylab("Shannon Diversity") +
  theme(plot.title =  element_text(size = 14), axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),  legend.position = "right") +
  stat_pvalue_manual(pairwise, label = "p.adj.signif", tip.length = 0.01, y.position = 5)
alpha_div_l

#ggsave('/mnt/md0/brd_project/figures/new_figures/bac_alpha_div.pdf', plot = alpha_div_l, width = 6, height = 4, dpi = 300)


hist(data.alpha$Shannon)
shapiro.test(data.alpha$Shannon)

######### BETA DIVERSITY ##########

ps.prop <- transform_sample_counts(data, function(otu) otu/sum(otu))

ldf <- psmelt(ps.prop)
gtone <- ldf %>% 
  group_by(Genus) %>% 
  summarize(Mean = mean(Abundance *100)) %>%
  arrange(-Mean) 

lgt_sub <- subset(gtone, Mean >= 1)

gtone <- ldf %>% 
  group_by(Family) %>% 
  summarize(Mean = mean(Abundance *100)) %>%
  arrange(-Mean) 
fgt_sub <- subset(gtone, Mean >= 1)

### make PCoA of larvae and frass individually.
ord.nmds.bray <- ordinate(ps.prop, method="PCoA", distance="bray")
loc_pal <- c('Ireland'='palegreen3', 'Australia'='plum4',
             "Canada"='orange3', "USA:Kansas" = "steelblue")
#plot_ordination(larvae, ord.nmds.bray, color="Diet", title="Bray PCoA") #+ geom_point(aes(size = 1)) + scale_size_continuous(guide = "none")
pcoa <- plot_ordination(ps.prop, ord.nmds.bray, title="Beta Diversity", color = "Location", shape = "Disease_State") + 
  geom_point(size = 3) + 
  scale_size_continuous(guide = "none") +
  #xlab("PCo1 (65.2%)")+ ylab("PCo2 (15.3%)") + theme(text = element_text(size = 12), legend.position = "none")   + 
  scale_color_manual(values = loc_pal) + 
  scale_shape_manual(values = c(4,1)) +
  #stat_ellipse(aes(color=Location, linetype = Disease_State))  + 
  theme_classic()
pcoa

# braydist <- distance(ps.prop, "bray")
# sampledf <-  data.frame(sample_data(ps.prop))
# ps.prop@sam_data<-sample_data(sampledf)
# 
# beta <- betadisper(braydist, sampledf$Disease_State)
# permutest(beta)
# 
# adonis2(braydist ~ Location * Disease_State, data = sampledf, permutations = 1000, by = "terms")
# 
# div <- alpha_div_l + pcoa
#ggsave('/mnt/md0/brd_project/figures/new_figures/bac_beta_div.pdf', plot = pcoa, width = 6, height = 4, dpi = 300)

########## Relative abundance ########
for_gg <- psmelt(ps.prop)
for_gg$Genus[which(!for_gg$Genus %in% lgt_sub$Genus)] <- "Other (<1%)"
#change factor levels
for_gg$Genus <- factor(for_gg$Genus, levels = c(unique(for_gg$Genus)[-c(10)], "Other (<1%)"))
#species
#for_gg$Sample <- gsub("_report_bracken_species","", for_gg$Sample)

#genus
for_gg$Sample <- gsub("_report_bracken_genuses","", for_gg$Sample)

#make pct level bar chart
ncols <- length(unique(for_gg$Genus))
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(16,"Paired"))(ncols)
mycolors[length(mycolors)] <- "darkgrey"
pct_larvae <- ggplot(data = for_gg, 
                     mapping = aes_string(x = "Sample",y = "Abundance"))  + 
  theme(legend.position = "right", 
        legend.text = element_text(face = "italic", size = 11),
        text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="fill") +
  facet_wrap(~Disease_State, scales="free_x") + 
  scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors)+ ggtitle("Relative Abundance")

div <- (alpha_div_l  + pcoa) / pct_larvae
div
ggsave("/mnt/md0/brd_project/figures/new_figures/bacterial_diversity_virus_samples.pdf", plot = div, width = 12, height = 10, dpi = 300)



######## Wilcoxon Method ##########3
library(rstatix)
gen_glom <- tax_glom(ps.prop,  taxrank = "Genus")
gen_df <- psmelt(gen_glom)


brd_nv <- subset(gen_df, Genus == "Mannheimia" | Genus == "Pasteurella" | Genus == "Mycoplasma" | Genus == "Histophilus" | Genus == "Trueperella" | Genus == "Moraxella")
#brd_nv <- brd[c(which(brd$Sample %in% novirus)),]

#stats 
stat.test <- brd_nv %>%
  group_by(Genus) %>%
  t_test(Abundance ~ Disease_State) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
#write.csv(stat.test, file = "/mnt/md0/baird_microbiome_data/tables/fungal_stats/fam_abun_stats.csv", row.names = F)
#add_coordinates
stat.test <- stat.test %>%
  add_xy_position(x = "Disease_State", dodge = 1)

#stat.test$y.position <- .209

bxp <- ggboxplot(brd_nv, x = "Disease_State", y = "Abundance", color = "Disease_State", facet.by = "Genus",scales = "free",
                 palette = brd_pal, add = "jitter", size = 0.9, title = "BRD Bacteria Abundance - Virus Samples") + 
  theme(legend.position = "right", axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 14)) + 
  ylab("Relative Abundance")

lsig_bxp <- bxp + 
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0.01,
    hide.ns = F
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
lsig_bxp

ggsave("/mnt/md0/brd_project/figures/new_figures/BRD_bacteria_virus_samples.pdf", plot = lsig_bxp, width = 8, height = 6, dpi = 300)




################### NO VIRUS SAMPLES ####################################


#import biom file from kraken-biom
data <-import_biom(bm, parseFunction=parse_taxonomy_default)
data@tax_table@.Data <- substring(data@tax_table@.Data, 4)
#rename columns to species level
colnames(data@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#View(data@tax_table@.Data)

#get rid of contaminant human reads (Only for genus and species)
#I added a taxa to this that is likely a false positive from Kraken
data <- subset_taxa(data, (!Genus %in% "Homo"))
data <- subset_taxa(data, (!Genus %in% "Callorhinchus"))
data <- subset_taxa(data, (!Genus %in% "Bos"))


#load in sample metadata
loc_df <- read.table(file = locs, header = F, sep = "\t")
symp_df <- read.table(file = symps, header = F, sep = "\t")
ctl_df <- read.table(file = ctl_samples, header = F, sep = "\t")
colnames(loc_df) <- c("SRA", "Location")
colnames(symp_df) <- c("SRA", "Disease_State")
md <- merge(loc_df, symp_df, by="SRA")
row.names(md) <- md$SRA
nv <- read.csv(novir)
novirus <- nv$X



#change names of samples from files
#change this based on what taxonomic level you imported from Bracken. e.g., if you imported genus level data,
##change "_._report_bracken_families" to "_._report_bracken_genuses" 
colnames(data@otu_table@.Data) <- gsub("\\.k2_bracken_genuses","", colnames(data@otu_table@.Data))

#create sample data for phyloseq object
samdf <- md[which(md$SRA %in% colnames(otu_table(data))), ]

#make sample names the rows
data@sam_data<-sample_data(samdf)
virus <- colnames(otu_table(data))[which(colnames(otu_table(data)) %in% novirus)]
#filter out samples with less than 5 OTUs detected.
data.alpha_f <- estimate_richness(data)
#data <- prune_samples(row.names(data.alpha_f[which(data.alpha_f$Observed >= 10),]) , data)
data <- prune_samples(colnames(otu_table(data))[which(!colnames(otu_table(data)) %in% ctl_df$V1)], data)
data <- prune_samples(virus, data)
data.alpha <- estimate_richness(data)

#data <- prune_taxa(taxa_sums(data) != 0, data) 

######### alpha diveristy ###########
brd_pal <- c("Healthy" = "#4682B4", "BRD" = "firebrick1")

shannon <- data.alpha[,"Shannon"]
shannon_div <- data.frame(row.names = rownames(data.alpha),shannon, sample_data(data)$Disease_State, sample_data(data)$Location)
colnames(shannon_div) <- c("Shannon", "Disease_State","Location")

#plot boxes
b <- ggplot(shannon_div, aes(x=Disease_State,y=Shannon, fill=Disease_State)) + geom_boxplot() 
alpha_div <- b + ggtitle("Alpha Diversity") + ylab("Shannon Diversity Index")+ scale_fill_manual(values = brd_pal) + theme_classic2()
alpha_div

#add stats
pairwise <- shannon_div %>% wilcox_test(Shannon ~ Disease_State) %>% adjust_pvalue(method = "bonferroni") %>% add_significance()
pairwise <- pairwise %>% add_xy_position(x = "Disease_State")

b_l <- ggboxplot(shannon_div, x="Disease_State",y="Shannon", color="Disease_State", add = "jitter", palette = brd_pal)
alpha_div_l <- b_l + ggtitle("Alpha Diversity") + ylab("Shannon Diversity") +
  theme(plot.title =  element_text(size = 14), axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),  legend.position = "right") +
  stat_pvalue_manual(pairwise, label = "p.adj.signif", tip.length = 0.01, y.position = 5)
alpha_div_l

#ggsave('/mnt/md0/brd_project/figures/new_figures/bac_alpha_div.pdf', plot = alpha_div_l, width = 6, height = 4, dpi = 300)


hist(data.alpha$Shannon)
shapiro.test(data.alpha$Shannon)

######### BETA DIVERSITY ##########

ps.prop <- transform_sample_counts(data, function(otu) otu/sum(otu))

ldf <- psmelt(ps.prop)
gtone <- ldf %>% 
  group_by(Genus) %>% 
  summarize(Mean = mean(Abundance *100)) %>%
  arrange(-Mean) 

lgt_sub <- subset(gtone, Mean >= 1)

gtone <- ldf %>% 
  group_by(Family) %>% 
  summarize(Mean = mean(Abundance *100)) %>%
  arrange(-Mean) 
fgt_sub <- subset(gtone, Mean >= 1)

### make PCoA of larvae and frass individually.
ord.nmds.bray <- ordinate(ps.prop, method="PCoA", distance="bray")
loc_pal <- c('Ireland'='palegreen3', 'Australia'='plum4',
             "Canada"='orange3', "USA:Kansas" = "steelblue")
#plot_ordination(larvae, ord.nmds.bray, color="Diet", title="Bray PCoA") #+ geom_point(aes(size = 1)) + scale_size_continuous(guide = "none")
pcoa <- plot_ordination(ps.prop, ord.nmds.bray, title="Beta Diversity", color = "Location", shape = "Disease_State") + 
  geom_point(size = 3) + 
  scale_size_continuous(guide = "none") +
  #xlab("PCo1 (65.2%)")+ ylab("PCo2 (15.3%)") + theme(text = element_text(size = 12), legend.position = "none")   + 
  scale_color_manual(values = loc_pal) + 
  scale_shape_manual(values = c(4,1)) +
  #stat_ellipse(aes(color=Location, linetype = Disease_State))  + 
  theme_classic()
pcoa

# braydist <- distance(ps.prop, "bray")
# sampledf <-  data.frame(sample_data(ps.prop))
# ps.prop@sam_data<-sample_data(sampledf)
# 
# beta <- betadisper(braydist, sampledf$Disease_State)
# permutest(beta)
# 
# adonis2(braydist ~ Location * Disease_State, data = sampledf, permutations = 1000, by = "terms")
# 
# div <- alpha_div_l + pcoa
#ggsave('/mnt/md0/brd_project/figures/new_figures/bac_beta_div.pdf', plot = pcoa, width = 6, height = 4, dpi = 300)

########## Relative abundance ########
for_gg <- psmelt(ps.prop)
for_gg$Genus[which(!for_gg$Genus %in% lgt_sub$Genus)] <- "Other (<1%)"
#change factor levels
for_gg$Genus <- factor(for_gg$Genus, levels = c(unique(for_gg$Genus)[-c(3)], "Other (<1%)"))
#species
#for_gg$Sample <- gsub("_report_bracken_species","", for_gg$Sample)

#genus
for_gg$Sample <- gsub("_report_bracken_genuses","", for_gg$Sample)

#make pct level bar chart
ncols <- length(unique(for_gg$Genus))
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(16,"Paired"))(ncols)
mycolors[length(mycolors)] <- "darkgrey"
pct_larvae <- ggplot(data = for_gg, 
                     mapping = aes_string(x = "Sample",y = "Abundance"))  + 
  theme(legend.position = "right", 
        legend.text = element_text(face = "italic", size = 11),
        text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="fill") +
  facet_wrap(~Disease_State, scales="free_x") + 
  scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors)+ ggtitle("Relative Abundance")

div <- (alpha_div_l  + pcoa) / pct_larvae
div
ggsave("/mnt/md0/brd_project/figures/new_figures/bacterial_diversity_no_virus_samples.pdf", plot = div, width = 12, height = 10, dpi = 300)


########### Diff Abun #########
### with help from https://oxfordcms.github.io/OCMS-blog/bioinformatics/Example-16S-rRNA-Analysis
# data@otu_table <- otu_table(data) + 1
# #ps.prop <- transform_sample_counts(data, function(otu) otu/sum(otu))
# hbdds <- phyloseq_to_deseq2(data, ~Disease_State + Location)
# 
# #hbdds <- estimateSizeFactors(hbdds, type = "iterate")
# hbdds <- DESeq(hbdds)
# 
# #keep <- rowSums(counts(hbdds) >= 10) >= 24
# #hbdds <- hbdds[keep,]
# res <- results(hbdds)
# 
# #transform counts
# nt <- normTransform(hbdds)
# ntdf <- data.frame(assay(nt))
# 
# 
# res2 <- data.frame(res@listData)
# rownames(res2) <- rownames(res)
# 
# #get diferentially abundant taxa
# res.diff <- res2[res2$padj < 0.05 & !(is.na(res2$padj)),]
# 
# # Get number different
# ndiff <- nrow(res.diff)
# df <- data.frame("Number DE" = ndiff)
# knitr::kable(df, caption="number of genes differentially expressed at p < 0.05")
# 
# 
# # Get matrix with diff ASV's
# ntdf.diff <- ntdf[rownames(res.diff),]
# ntdf.diff <- ntdf.diff[order(rownames(ntdf.diff)),]
# 
# #get max annotation for taxa
# sig_dif_taxa <- data.frame(tax_table(data)[which(rownames(tax_table(data)) %in% rownames(ntdf.diff)),] )
# sig_dif_taxa <- sig_dif_taxa[order(rownames(sig_dif_taxa)),]
# #sig_dif_taxa$max_annotation <- sig_dif_taxa[cbind(1:nrow(sig_dif_taxa), max.col(!is.na(sig_dif_taxa), ties.method = 'last'))]
# #sig_dif_taxa <- sig_dif_taxa[-which(duplicated(sig_dif_taxa$Genus)),] # uncultured bacteria
# #row.names(sig_dif_taxa) <- sig_dif_taxa$Genus
# 
# #ntdf.diff$tax <- sig_dif_taxa$Phylum
# 
# df <- as.data.frame(colData(hbdds)[,c(2,3)])
# colnames(ntdf.diff) <- rownames(df)
# 
# 
# hm <- pheatmap(as.matrix(ntdf.diff), cluster_rows=T, show_rownames=T, show_colnames = F,
#                cluster_cols=T, annotation_col=df, labels_row = sig_dif_taxa$Genus, fontsize = 9)
# 

######## Wilcoxon Method ##########3
library(rstatix)
gen_glom <- tax_glom(ps.prop,  taxrank = "Genus")
gen_df <- psmelt(gen_glom)


brd_nv <- subset(gen_df, Genus == "Mannheimia" | Genus == "Pasteurella" | Genus == "Mycoplasma" | Genus == "Histophilus" | Genus == "Trueperella" | Genus == "Moraxella")
#brd_nv <- brd[c(which(brd$Sample %in% novirus)),]

#stats 
stat.test <- brd_nv %>%
  group_by(Genus) %>%
  t_test(Abundance ~ Disease_State) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
#write.csv(stat.test, file = "/mnt/md0/baird_microbiome_data/tables/fungal_stats/fam_abun_stats.csv", row.names = F)
#add_coordinates
stat.test <- stat.test %>%
  add_xy_position(x = "Disease_State", dodge = 1)

#stat.test$y.position <- .209

bxp <- ggboxplot(brd_nv, x = "Disease_State", y = "Abundance", color = "Disease_State", facet.by = "Genus",scales = "free",
                 palette = brd_pal, add = "jitter", size = 0.9, title = "BRD Bacteria Abundance - No Virus Samples") + 
  theme(legend.position = "right", axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 14)) + 
  ylab("Relative Abundance")

lsig_bxp <- bxp + 
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0.01,
    hide.ns = F
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
lsig_bxp

ggsave("/mnt/md0/brd_project/figures/new_figures/BRD_bacteria_no_virus_samples.pdf", plot = lsig_bxp, width = 8, height = 6, dpi = 300)




############### ALL SAMPLES #######################
#import biom file from kraken-biom
data <-import_biom(bm, parseFunction=parse_taxonomy_default)
data@tax_table@.Data <- substring(data@tax_table@.Data, 4)
#rename columns to species level
colnames(data@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#View(data@tax_table@.Data)

#get rid of contaminant human reads (Only for genus and species)
#I added a taxa to this that is likely a false positive from Kraken
data <- subset_taxa(data, (!Genus %in% "Homo"))
data <- subset_taxa(data, (!Genus %in% "Callorhinchus"))
data <- subset_taxa(data, (!Genus %in% "Bos"))


#load in sample metadata
loc_df <- read.table(file = locs, header = F, sep = "\t")
symp_df <- read.table(file = symps, header = F, sep = "\t")
ctl_df <- read.table(file = ctl_samples, header = F, sep = "\t")
colnames(loc_df) <- c("SRA", "Location")
colnames(symp_df) <- c("SRA", "Disease_State")
md <- merge(loc_df, symp_df, by="SRA")
row.names(md) <- md$SRA
nv <- read.csv(novir)
novirus <- nv$X



#change names of samples from files
#change this based on what taxonomic level you imported from Bracken. e.g., if you imported genus level data,
##change "_._report_bracken_families" to "_._report_bracken_genuses" 
colnames(data@otu_table@.Data) <- gsub("\\.k2_bracken_genuses","", colnames(data@otu_table@.Data))

#create sample data for phyloseq object
samdf <- md[which(md$SRA %in% colnames(otu_table(data))), ]

#make sample names the rows
data@sam_data<-sample_data(samdf)
virus <- colnames(otu_table(data))[which(colnames(otu_table(data)) %in% novirus)]
#filter out samples with less than 5 OTUs detected.
data.alpha_f <- estimate_richness(data)
#data <- prune_samples(row.names(data.alpha_f[which(data.alpha_f$Observed >= 10),]) , data)
data <- prune_samples(colnames(otu_table(data))[which(!colnames(otu_table(data)) %in% ctl_df$V1)], data)
#data <- prune_samples(virus, data)
data.alpha <- estimate_richness(data)

#data <- prune_taxa(taxa_sums(data) != 0, data) 

######### alpha diveristy ###########
brd_pal <- c("Healthy" = "#4682B4", "BRD" = "firebrick1")

shannon <- data.alpha[,"Shannon"]
shannon_div <- data.frame(row.names = rownames(data.alpha),shannon, sample_data(data)$Disease_State, sample_data(data)$Location)
colnames(shannon_div) <- c("Shannon", "Disease_State","Location")

#plot boxes
b <- ggplot(shannon_div, aes(x=Disease_State,y=Shannon, fill=Disease_State)) + geom_boxplot() 
alpha_div <- b + ggtitle("Alpha Diversity") + ylab("Shannon Diversity Index")+ scale_fill_manual(values = brd_pal) + theme_classic2()
alpha_div

#add stats
pairwise <- shannon_div %>% wilcox_test(Shannon ~ Disease_State) %>% adjust_pvalue(method = "bonferroni") %>% add_significance()
pairwise <- pairwise %>% add_xy_position(x = "Disease_State")

b_l <- ggboxplot(shannon_div, x="Disease_State",y="Shannon", color="Disease_State", add = "jitter", palette = brd_pal)
alpha_div_l <- b_l + ggtitle("Alpha Diversity") + ylab("Shannon Diversity") +
  theme(plot.title =  element_text(size = 14), axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),  legend.position = "right") +
  stat_pvalue_manual(pairwise, label = "p.adj.signif", tip.length = 0.01, y.position = 5)
alpha_div_l

#ggsave('/mnt/md0/brd_project/figures/new_figures/bac_alpha_div.pdf', plot = alpha_div_l, width = 6, height = 4, dpi = 300)


hist(data.alpha$Shannon)
shapiro.test(data.alpha$Shannon)

######### BETA DIVERSITY ##########

ps.prop <- transform_sample_counts(data, function(otu) otu/sum(otu))

ldf <- psmelt(ps.prop)
gtone <- ldf %>% 
  group_by(Genus) %>% 
  summarize(Mean = mean(Abundance *100)) %>%
  arrange(-Mean) 

lgt_sub <- subset(gtone, Mean >= 1)

gtone <- ldf %>% 
  group_by(Family) %>% 
  summarize(Mean = mean(Abundance *100)) %>%
  arrange(-Mean) 
fgt_sub <- subset(gtone, Mean >= 1)

### make PCoA of larvae and frass individually.
ord.nmds.bray <- ordinate(ps.prop, method="PCoA", distance="bray")
loc_pal <- c('Ireland'='palegreen3', 'Australia'='plum4',
             "Canada"='orange3', "USA:Kansas" = "steelblue")
#plot_ordination(larvae, ord.nmds.bray, color="Diet", title="Bray PCoA") #+ geom_point(aes(size = 1)) + scale_size_continuous(guide = "none")
pcoa <- plot_ordination(ps.prop, ord.nmds.bray, title="Beta Diversity", color = "Location", shape = "Disease_State") + 
  geom_point(size = 3) + 
  scale_size_continuous(guide = "none") +
  #xlab("PCo1 (65.2%)")+ ylab("PCo2 (15.3%)") + theme(text = element_text(size = 12), legend.position = "none")   + 
  scale_color_manual(values = loc_pal) + 
  scale_shape_manual(values = c(4,1)) +
  #stat_ellipse(aes(color=Location, linetype = Disease_State))  + 
  theme_classic()
pcoa

# braydist <- distance(ps.prop, "bray")
# sampledf <-  data.frame(sample_data(ps.prop))
# ps.prop@sam_data<-sample_data(sampledf)
# 
# beta <- betadisper(braydist, sampledf$Disease_State)
# permutest(beta)
# 
# adonis2(braydist ~ Location * Disease_State, data = sampledf, permutations = 1000, by = "terms")
# 
# div <- alpha_div_l + pcoa
#ggsave('/mnt/md0/brd_project/figures/new_figures/bac_beta_div.pdf', plot = pcoa, width = 6, height = 4, dpi = 300)

########## Relative abundance ########
for_gg <- psmelt(ps.prop)
for_gg$Genus[which(!for_gg$Genus %in% lgt_sub$Genus)] <- "Other (<1%)"
#change factor levels
for_gg$Genus <- factor(for_gg$Genus, levels = c(unique(for_gg$Genus)[-c(4)], "Other (<1%)"))
#species
#for_gg$Sample <- gsub("_report_bracken_species","", for_gg$Sample)

#genus
for_gg$Sample <- gsub("_report_bracken_genuses","", for_gg$Sample)

#make pct level bar chart
ncols <- length(unique(for_gg$Genus))
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(16,"Paired"))(ncols)
mycolors[length(mycolors)] <- "darkgrey"
pct_larvae <- ggplot(data = for_gg, 
                     mapping = aes_string(x = "Sample",y = "Abundance"))  + 
  theme(legend.position = "right", 
        legend.text = element_text(face = "italic", size = 11),
        text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="fill") +
  facet_wrap(~Disease_State, scales="free_x") + 
  scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors)+ ggtitle("Relative Abundance")

div <- (alpha_div_l  + pcoa) / pct_larvae
div
ggsave("/mnt/md0/brd_project/figures/new_figures/bacterial_diversity_all_samples.pdf", plot = div, width = 12, height = 10, dpi = 300)


########### Diff Abun #########
### with help from https://oxfordcms.github.io/OCMS-blog/bioinformatics/Example-16S-rRNA-Analysis
# data@otu_table <- otu_table(data) + 1
# #ps.prop <- transform_sample_counts(data, function(otu) otu/sum(otu))
# hbdds <- phyloseq_to_deseq2(data, ~Disease_State + Location)
# 
# #hbdds <- estimateSizeFactors(hbdds, type = "iterate")
# hbdds <- DESeq(hbdds)
# 
# #keep <- rowSums(counts(hbdds) >= 10) >= 24
# #hbdds <- hbdds[keep,]
# res <- results(hbdds)
# 
# #transform counts
# nt <- normTransform(hbdds)
# ntdf <- data.frame(assay(nt))
# 
# 
# res2 <- data.frame(res@listData)
# rownames(res2) <- rownames(res)
# 
# #get diferentially abundant taxa
# res.diff <- res2[res2$padj < 0.05 & !(is.na(res2$padj)),]
# 
# # Get number different
# ndiff <- nrow(res.diff)
# df <- data.frame("Number DE" = ndiff)
# knitr::kable(df, caption="number of genes differentially expressed at p < 0.05")
# 
# 
# # Get matrix with diff ASV's
# ntdf.diff <- ntdf[rownames(res.diff),]
# ntdf.diff <- ntdf.diff[order(rownames(ntdf.diff)),]
# 
# #get max annotation for taxa
# sig_dif_taxa <- data.frame(tax_table(data)[which(rownames(tax_table(data)) %in% rownames(ntdf.diff)),] )
# sig_dif_taxa <- sig_dif_taxa[order(rownames(sig_dif_taxa)),]
# #sig_dif_taxa$max_annotation <- sig_dif_taxa[cbind(1:nrow(sig_dif_taxa), max.col(!is.na(sig_dif_taxa), ties.method = 'last'))]
# #sig_dif_taxa <- sig_dif_taxa[-which(duplicated(sig_dif_taxa$Genus)),] # uncultured bacteria
# #row.names(sig_dif_taxa) <- sig_dif_taxa$Genus
# 
# #ntdf.diff$tax <- sig_dif_taxa$Phylum
# 
# df <- as.data.frame(colData(hbdds)[,c(2,3)])
# colnames(ntdf.diff) <- rownames(df)
# 
# 
# hm <- pheatmap(as.matrix(ntdf.diff), cluster_rows=T, show_rownames=T, show_colnames = F,
#                cluster_cols=T, annotation_col=df, labels_row = sig_dif_taxa$Genus, fontsize = 9)
# 

######## Wilcoxon Method ##########3
library(rstatix)
gen_glom <- tax_glom(ps.prop,  taxrank = "Genus")
gen_df <- psmelt(gen_glom)


brd_nv <- subset(gen_df, Genus == "Mannheimia" | Genus == "Pasteurella" | Genus == "Mycoplasma" | Genus == "Histophilus" | Genus == "Trueperella" | Genus == "Moraxella")
#brd_nv <- brd[c(which(brd$Sample %in% novirus)),]

#stats 
stat.test <- brd_nv %>%
  group_by(Genus) %>%
  t_test(Abundance ~ Disease_State) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
#write.csv(stat.test, file = "/mnt/md0/baird_microbiome_data/tables/fungal_stats/fam_abun_stats.csv", row.names = F)
#add_coordinates
stat.test <- stat.test %>%
  add_xy_position(x = "Disease_State", dodge = 1)

#stat.test$y.position <- .209

bxp <- ggboxplot(brd_nv, x = "Disease_State", y = "Abundance", color = "Disease_State", facet.by = "Genus",scales = "free",
                 palette = brd_pal, add = "jitter", size = 0.9, title = "BRD Bacteria Abundance - No Virus Samples") + 
  theme(legend.position = "right", axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 14)) + 
  ylab("Relative Abundance")

lsig_bxp <- bxp + 
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0.01,
    hide.ns = F
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1)))
lsig_bxp

ggsave("/mnt/md0/brd_project/figures/new_figures/BRD_bacteria_all_samples.pdf", plot = lsig_bxp, width = 8, height = 6, dpi = 300)
