### Analysis of BRD-associated microbiome results from Kraken/Bracken
## Hunter K. Walt 7/14/2025

library(tidyverse)
library(phyloseq)
library(vegan)
library(RColorBrewer)
library(ggplot2)
library(stats)
library(rstatix)
library(ggpubr)
library(DESeq2)
library(patchwork)

#working directory
wd <- "/mnt/md0/brd_project/bacteria/standard_kraken/"
setwd(wd)

#read in metadata
locs <- "/mnt/md0/brd_project/SRA_locs.txt"
symps <- "/mnt/md0/brd_project/SRA_symptoms.txt"
bm <- "BRD_bracken_genera_REVISED.biom" 
ctl_samples <- "/mnt/md0/brd_project/control_samples_to_remove.txt"

#sequencing type
nano <- "/mnt/md0/brd_project/virus_coverage/depth_stats/nanopore_mapped_final.txt"
illum <- "/mnt/md0/brd_project/virus_coverage/depth_stats/illumina_mapped_final.txt"
nan <- read.table(nano, header = F)
ill <- read.table(illum, header = F)
colnames(nan) <- "SRA"
colnames(ill) <- "SRA"
nan$Seqtype <- "Nanopore"
ill$Seqtype <- "Illumina"
st <- rbind(nan, ill)




############### ALL SAMPLES #######################
#import biom file from kraken-biom
data <-import_biom(bm, parseFunction=parse_taxonomy_default)
data@tax_table@.Data <- substring(data@tax_table@.Data, 4)
#rename columns to species level
colnames(data@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#View(data@tax_table@.Data)

#get rid of contaminant reads 
#I added a taxa to this that is likely a false positive from Kraken
data <- subset_taxa(data, (!Phylum %in% "Chordata"))
data <- subset_taxa(data, (!Phylum %in% grep("viricota", unique(data.frame(tax_table(data))$Phylum), value = T))) #get rid of viruses
phy <- tax_glom(data, taxrank = "Phylum")

# psm <- psmelt(physeq = phy)
# psm$Sample <- gsub("\\.k2_bracken_genuses", "", psm$Sample) 
# sum_cts <- psm %>% 
#   group_by(Sample) %>% 
#   summarize(total_microbial_counts = sum(Abundance))

#load in sample metadata
loc_df <- read.table(file = locs, header = F, sep = "\t")
symp_df <- read.table(file = symps, header = F, sep = "\t")
ctl_df <- read.table(file = ctl_samples, header = F, sep = "\t")
colnames(loc_df) <- c("SRA", "Location")
colnames(symp_df) <- c("SRA", "Disease_State")
md <- merge(loc_df, symp_df, by="SRA")
row.names(md) <- md$SRA



#change names of samples from files
#change this based on what taxonomic level you imported from Bracken. e.g., if you imported genus level data,
##change "_._report_bracken_families" to "_._report_bracken_genuses" 
#colnames(data@otu_table@.Data) <- gsub("\\.k2_bracken_species","", colnames(data@otu_table@.Data))
colnames(data@otu_table@.Data) <- gsub("\\.k2_bracken_genuses","", colnames(data@otu_table@.Data))
#create sample data for phyloseq object
samdf <- md[which(md$SRA %in% colnames(otu_table(data))), ]
samdf <- merge(samdf, st, by = "SRA")
row.names(samdf) <- samdf$SRA
#make sample names the rows
data@sam_data<-sample_data(samdf)
data <- prune_samples(colnames(otu_table(data))[which(!colnames(otu_table(data)) %in% ctl_df$V1)], data)
#data <- subset_samples(data, Seqtype != "Nanopore" | Location != "Ireland")
data.alpha <- estimate_richness(data, measures = "Shannon")
data@otu_table <- otu_table(data) + 1

#normalize data
dds = phyloseq_to_deseq2(data, ~Location + Disease_State)
# Make a copy of your phyloseq object, which you will then modify with VST values
physeqvsd = data
dds <- DESeq(dds)
vsd = getVarianceStabilizedData(dds)
mddf <- samdf[which(samdf$SRA %in% rownames(t(vsd))),]
new <- limma::removeBatchEffect(vsd, batch = mddf$Seqtype, covariates = model.matrix(~Location, mddf)[,-1])
cts <- counts(dds, normalized = T)
otu_table(physeqvsd) <- otu_table(vsd, taxa_are_rows = TRUE)
#distance(physeqvsd, ...)
physeq2 = transform_sample_counts(physeqvsd, round)

#normalization introduced negative numbers, so transforming the data so that the lowest numbers are zero.
physeq2@otu_table <- otu_table(physeq2) + abs(min(otu_table(physeq2)))


data.alpha <- estimate_richness(physeq2, measures = "Shannon")

######### alpha diveristy ###########
brd_pal <- c("Healthy" = "#4682B4", "BRD" = "firebrick1")

shannon <- data.alpha[,"Shannon"]
shannon_div <- data.frame(row.names = rownames(data.alpha),shannon, sample_data(data)$Disease_State, sample_data(data)$Location)
colnames(shannon_div) <- c("Shannon", "Disease_State","Location")
shannon_div <- merge(shannon_div, st, by.x = "row.names", by.y = "SRA")
write.csv(shannon_div, file = "shannon_diversity_REVISED.csv") #write out data

shannon_model <- lm(Shannon ~ Disease_State + Location + Seqtype, data = shannon_div)
shannon_pvalue <- broom::tidy(summary(shannon_model))

#plot boxes
b <- ggplot(shannon_div, aes(x=Disease_State,y=Shannon, fill=Disease_State)) + geom_boxplot() 
alpha_div <- b + ggtitle("Alpha Diversity") + ylab("Shannon Diversity Index")+ scale_fill_manual(values = brd_pal) + theme_classic2()
#alpha_div

#add stats
pairwise <- subset(shannon_div, Seqtype == "Illumina" & Location != "Australia") %>% wilcox_test(Shannon ~ Disease_State) %>% adjust_pvalue(method = "bonferroni") %>% add_significance()
pairwise <- pairwise %>% add_xy_position(x = "Disease_State")



b_l <- ggboxplot(subset(shannon_div, Location != "USA:Kansas"), x="Disease_State" ,y="Shannon",
                 color="Disease_State", add = "jitter", palette = brd_pal)
alpha_div_a <- b_l + ggtitle("Alpha Diversity") + ylab("Shannon Diversity") +
  xlab("Disease State") +
  theme(plot.title =  element_text(size = 12), 
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(), 
        legend.position = "none") + facet_wrap(~ Location + Seqtype) +  stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") +
  #stat_pvalue_manual(pairwise, label = "p.adj.signif", tip.length = 0.01, y.position = 7.8) 
  stat_compare_means(label =  "p.signif", label.x = 1.5, label.y = 7.74) + labs(color = "Disease State", fill = "Disease State")
alpha_div_a 

######### BETA DIVERSITY ##########

#ps.prop <- transform_sample_counts(physeqvsd, function(otu) otu/sum(otu))

#(Seqtype != "Nanopore" | Location != "Ireland")

ps.prop <- physeq2


# ldf <- psmelt(ps.prop)
# gtone <- ldf %>% 
#   group_by(Genus) %>% 
#   summarize(Mean = mean(Abundance *100)) %>%
#   arrange(-Mean) 
# 
# lgt_sub <- subset(gtone, Mean >= 0.1)
# 
# gtone <- ldf %>% 
#   group_by(Family) %>% 
#   summarize(Mean = mean(Abundance *100)) %>%
#   arrange(-Mean) 
# fgt_sub <- subset(gtone, Mean >= 1)

### make PCoA of larvae and frass individually.
ord.nmds.bray <- ordinate(ps.prop, method="PCoA", distance="bray")
loc_pal <- c('Ireland'='palegreen3', 'Australia'='plum4',
             "Canada"='orange3', "USA:Kansas" = "steelblue")
#plot_ordination(larvae, ord.nmds.bray, color="Diet", title="Bray PCoA") #+ geom_point(aes(size = 1)) + scale_size_continuous(guide = "none")
pcoa <- plot_ordination(ps.prop, ord.nmds.bray, title="Beta Diversity PCoA", color = "Location", shape = "Disease_State") + 
  geom_point(size = 3) + 
  scale_size_continuous(guide = "none") +
  xlab("PCo1 (59.4%)")+ ylab("PCo2 (25.2%)") + theme(text = element_text(size = 12), legend.position = "none")   + 
  scale_color_manual(values = loc_pal) + 
  scale_shape_manual(values = c(4,1)) +
  #stat_ellipse(aes(linetype = Disease_State))  + 
  theme_classic() +  labs(shape = "Disease State") +
  annotate("text", label = "Disease State: p = 0.002\nLocation: p < 0.001\nSequencing Platform: p < 0.001", x = 0.1, y = 0.2, size = 3)
pcoa

braydist <- phyloseq::distance(ps.prop, "bray")
#subtst <- prune_samples(dat$SRA, ps.prop)
sampledf <-  data.frame(sample_data(ps.prop))
#braysub <- phyloseq::distance(subtst, "bray")

ps.prop@sam_data<-sample_data(sampledf)

beta <- betadisper(braydist, sampledf$Disease_State)
permutest(beta)

adonis2(braydist ~ Disease_State + Location + Seqtype, data = sampledf, permutations = 1000, by = "terms")


######## Wilcoxon Method Diff Abundance ######## ##3
library(rstatix)
gen_glom <- tax_glom(ps.prop,  taxrank = "Genus")
gen_df <- psmelt(gen_glom)

mp <- lmerTest::lmer(Abundance ~ Disease_State + (1|Location) + (1|Seqtype), data = subset(gen_df, Genus == "Pasteurella"))
summary(mp)

mp <- lmerTest::lmer(Abundance ~ Disease_State + (1|Location) + (1|Seqtype), data = subset(gen_df, Genus == "Histophilus"))
summary(mp)

mp <- lmerTest::lmer(Abundance ~ Disease_State + (1|Location) + (1|Seqtype), data = subset(gen_df, Genus == "Mycoplasma"))
summary(mp)

mp <- lmerTest::lmer(Abundance ~ Disease_State + (1|Location) + (1|Seqtype), data = subset(gen_df, Genus == "Moraxella"))
summary(mp)

mp <- lmerTest::lmer(Abundance ~ Disease_State + (1|Location) + (1|Seqtype), data = subset(gen_df, Genus == "Mannheimia"))
summary(mp)

#subset BRD taxa
brd_nv <- subset(gen_df, Genus == "Mannheimia" | Genus == "Pasteurella" | Genus == "Mycoplasmopsis" | Genus == "Histophilus" | Genus == "Moraxella")
brd_nv$Disease_State <- factor(brd_nv$Disease_State, levels = c("Healthy", "BRD"))

qqnorm(brd_nv$Abundance, pch = 1, frame = FALSE)
qqline(brd_nv$Abundance, col = "steelblue", lwd = 2)

#stats 
stat.test <- subset(brd_nv, Location != "USA:Kansas") %>%
  group_by(Genus) %>%
  wilcox_test(Abundance ~ Disease_State) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
#write.csv(stat.test, file = "/mnt/md0/baird_microbiome_data/tables/fungal_stats/fam_abun_stats.csv", row.names = F)
#add_coordinates
stat.test <- stat.test %>%
  add_xy_position(x = "Disease_State", dodge = 1)


bxp <- ggdotplot(subset(brd_nv, Location != "USA:Kansas" & Abundance > 4), x = "Disease_State", y = "Abundance", color = "Disease_State", facet.by = "Genus",scales = "free",
                 palette = brd_pal, add = "jitter", title = "BRD Bacteria Abundance", size = 0.5, fill = "Disease_State") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        #axis.text.x = element_blank(), 
        #axis.ticks.x = element_blank(), 
        text = element_text(size = 12)) + 
  ylab("Normalized Abundance")  + stat_summary(fun = mean, geom = "point", shape = 15, size = 3, color = "black") #+ stat_compare_means()

bxp


mycomps <- list(c("Canada", "Australia"), c("Ireland", "Australia"), c("Ireland", "Canada"))
lsig_bxp <- bxp + 
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif", tip.length = 0.01, label.size = 4,
    hide.ns = F
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) + labs(color = "Disease State") 
lsig_bxp 

#control facel labels
bxp_f <- facet(lsig_bxp, nrow = 1, facet.by = "Genus", panel.labs.font = list(face = "italic"))

#make microbiome figure
alt_microbiome_fig <- (alpha_div_a | pcoa) / bxp_f + plot_annotation(tag_levels = "A")
alt_microbiome_fig

ggsave("/mnt/md0/brd_project/figures/new_figures/final_figures/BRD_final_bacteria_figurefinal_REVISED.pdf", plot = alt_microbiome_fig, width = 11, height = 7.3, dpi = 300)


### make co-occurence
PA <- ifelse(brd_nv$Abundance > 0, 1, 0)
brd_nv$Present <- PA
brd_sub <- brd_nv[,c(4,13,14)]
brd_bac_pa <- brd_sub |> pivot_wider(names_from = Genus, values_from = Present)
write.csv( brd_bac_pa, file = "/mnt/md0/brd_project/bacteria_presence_absence_REVISED.csv", quote = F, row.names = F)



############### DIFFERENTIAL ABUNDANCE #############



#ancombc2 seqtype is a random effect
out <- ANCOMBC::ancombc2(
  data = data,
  tax_level = "Genus",
  fix_formula = "Disease_State + Location + Seqtype",
  #rand_formula = NULL,
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  alpha = 0.05
  #global = FALSE,
  #pairwise = FALSE, 
  #dunnet = FALSE, 
  #trend = FALSE, 
)

#get result object
res <- out$res
res_disease <- data.frame(
  Taxon = res$taxon,
  #beta = res$beta_Disease_StateHealthy,
  se = res$se_Disease_StateHealthy,
  p = res$p_Disease_StateHealthy,
  q = res$q_Disease_StateHealthy,
  diff_abn = res$diff_Disease_StateHealthy
)


res <- out$res

df <- data.frame(
  Taxon = res$taxon,
  beta = res$lfc_Disease_StateHealthy,
  q = res$q_Disease_StateHealthy
)

ggplot(df, aes(beta, -log10(q))) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme_classic()

sig <- subset(df, q < 0.05 & abs(beta) > 1 ) %>% mutate(pos = beta >= 0)

dafull <- ggplot(sig,
       aes(x = beta*-1,
           y = reorder(Taxon, beta),
           color = pos)) +
  geom_point(size = 3) +
  geom_segment(aes(xend=0, yend = Taxon), size = 1.2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values = c("firebrick1", "#4682B4" ), 
                     labels = c("BRD", "Healthy")) +
  labs(color = "Differential Abundance", 
       title = "Differentially Abundant Genera") +
  ylab("Taxon") +
  xlab("Log2 Fold Change") +
  theme_classic()  


#stratify analysis
data_ill <- subset_samples(data, Seqtype == "Illumina")
data_nan <- subset_samples(data, Location == "Canada")





#ancombc2 seqtype is a random effect
out_ill <- ANCOMBC::ancombc2(
  data = data_ill,
  tax_level = "Genus",
  fix_formula = "Disease_State + Location",
  #rand_formula = NULL,
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  alpha = 0.05
  #global = FALSE,
  #pairwise = FALSE, 
  #dunnet = FALSE, 
  #trend = FALSE, 
)


res_ill <- out_ill$res

df_ill <- data.frame(
  Taxon = res_ill$taxon,
  beta = res_ill$lfc_Disease_StateHealthy,
  q = res_ill$q_Disease_StateHealthy
)

sig_ill <- subset(df_ill, q < 0.05 & abs(beta) > 1 )  %>% mutate(pos = beta >= 0)


illumina_fig <- ggplot(sig_ill,
       aes(x = beta*-1,
           y = reorder(Taxon, beta),
           color = pos)) +
  geom_point(size = 2) +
  geom_segment(aes(xend=0, yend = Taxon), size = 1) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values = c("firebrick1", "#4682B4" ), 
                     labels = c("BRD", "Healthy")) +
  labs(color = "Differential Abundance", 
       title = "Differentially Abundant Genera - Illumina Only") +
  ylab("Taxon") +
  xlab("Log2 Fold Change") +
  theme_classic()  



#ancombc2 seqtype is a random effect
out_nan <- ANCOMBC::ancombc2(
  data = data_nan,
  tax_level = "Genus",
  fix_formula = "Disease_State",
  #rand_formula = NULL,
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  alpha = 0.05
  #global = FALSE,
  #pairwise = FALSE, 
  #dunnet = FALSE, 
  #trend = FALSE, 
)


res_nan <- out_nan$res

df_nan <- data.frame(
  Taxon = res_nan$taxon,
  beta = res_nan$lfc_Disease_StateHealthy,
  q = res_nan$q_Disease_StateHealthy
)

sig_nan <- subset(df_nan, q < 0.05 & abs(beta) > 1 )


ggplot(sig_nan,
       aes(x = beta,
           y = reorder(Taxon, beta))) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic()


des <- 
"
AB
AC
"

revised_figure <- alpha_div_a + pcoa + dafull + plot_layout(design = des) + plot_annotation(tag_levels = "A")

ggsave(plot = revised_figure, filename = "/mnt/md0/brd_project/figures/revised_figures/REVISED_MICROBIOME.pdf", width = 11, height = 6.5)
ggsave(plot = illumina_fig, filename = "/mnt/md0/brd_project/figures/revised_figures/Illumina_DA.pdf", width = 9, height = 5)


write.csv(gen_df, file = "bacteria_abundances_REVISED.csv")
