#cow immune script
library(tximport)
library(DESeq2)
library(readr)
library("tximportData")
library(GenomicFeatures)
library(ggpubr)
library(ggrepel)
library(org.Bt.eg.db)
organism = org.Bt.eg.db
#library(pathview)
library(edgeR)
library(microbiome)
#library(tidyplots)
library(cowplot)
library(clusterProfiler)
library(rrvgo)

#disease state metadata
symps <- "/mnt/md0/brd_project/SRA_symptoms.txt"
symp <- read.table(file = symps, header = F)
names(symp) <- c("SRA", "Disease_state")

#location data
locs <- "/mnt/md0/brd_project/SRA_locs.txt"
loc <- read.table(file = locs, header = F)
names(loc) <- c("SRA", "Location")
loc$Location <- gsub("USA:Kansas", "USA.Kansas", loc$Location)

#sequencing type data
nano <- "/mnt/md0/brd_project/virus_coverage/depth_stats/nanopore_mapped.txt"
illum <- "/mnt/md0/brd_project/virus_coverage/depth_stats/illumina_mapped.txt"
nan <- read.table(nano, header = F)
ill <- read.table(illum, header = F)
colnames(nan) <- "SRA"
colnames(ill) <- "SRA"
nan$Seqtype <- "Nanopore"
ill$Seqtype <- "Illumina"
st <- rbind(nan, ill)
ols <- c("SRR23246845", "SRR23246847", "SRR23246844", "SRR23246846", "SRR1776513")

#samples with no virus
novir <- "/mnt/md0/brd_project/no_virus_rpkm.csv"
ctl_samples <- "/mnt/md0/brd_project/control_samples_to_remove.txt"
nv <- read.csv(novir)
names(nv) <- c("SRA", "NotImportant")
ctl <- read.table(file = ctl_samples, header = F)

#path to quant files
directory = "/mnt/md0/brd_project/cow_immune/salmon_long_filtered_250"

#quantfiles list: these are just a list of the sample names
CI <- read.table(file = "/mnt/md0/brd_project/cow_immune/samps_with_nanopore.txt", header = F)
#CI <- data.frame(CI[-which(CI$V1 %in% nv$SRA),]) # FOR VIRUSES
#CI <- data.frame(CI[which(CI$V1 %in% nv$SRA),]) # FOR NO VIRUSES
colnames(CI) <- "V1"
CI <- data.frame(CI[-which(CI$V1 %in% c(ctl$V1, ols)),])
colnames(CI) <- "V1"

#make file paths s R can find the files
CI_files <- file.path(directory, paste0(CI$V1,"_quant"), "quant.sf")
names(CI_files) <- CI$V1
#CI_files <- CI_files[-c(167,200)] # ONLY FOR NO VIRUS SAMPLES
#CI_files <- CI_files[can_sras] # ONLY FOR NO VIRUS SAMPLES

#annotation file to make txdb
gff = "/mnt/md0/brd_project/cow_immune/GCF_002263795.3_ARS-UCD2.0_genomic.gff.gz"

#to save txdb
txdb.filename = "/mnt/md0/brd_project/cow_immune/cow_annotation.sqlite"

#make and save txdb
txdb <- txdbmaker::makeTxDbFromGFF(file = gff)
saveDb(txdb, txdb.filename)

#maketx2gene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene <- tx2gene[-c(1),]

#I had to remove the following SRA datasets because salmon failed for them. They were genomic data, anyway.
#rm -r SRR1776541_quant/ SRR1776564_quant/ SRR1776524_quant/ SRR1776531_quant/ SRR1776533_quant/ SRR1776537_quant/

#make object TI
txi_CI <- tximport(CI_files, type = "salmon", tx2gene = tx2gene)
names(txi_CI)

head(txi_CI$counts)
txi_CI$counts <- txi_CI$counts + 1 
which(colSums(txi_CI$counts) < 500)

#import CI to DEseq
CI_table <- symp[which(symp$SRA %in% colnames(txi_CI$counts)),]
#CI_table[c(grep("Illumina", CI_table$Seqtype)), ]

rownames(CI_table) <- colnames(txi_CI$counts)
names(CI_table) <- c("SRA", "Disease_State")
CI_table <- merge(CI_table, loc, by="SRA")
#CI_table <- subset(CI_table, select = c("Disease_State", "Location"))
CI_table <- merge(CI_table, st, by= "SRA")
#CI_table <- subset(CI_table, select = c("Disease_State"))

CIdds <- DESeqDataSetFromTximport(txi_CI, CI_table, ~Disease_State + Location + Seqtype)
#CIdds <- DESeqDataSetFromTximport(txi_CI, CI_table, ~Disease_State)
#can_sras <- CI_table$SRA[grep("Canada", CI_table$Location)] # FOR NO VIRUS ONLY
#CIdds <- CIdds[,can_sras]

keep <- rowSums(counts(CIdds) >= 10) >= 20
#keep <- rowSums(counts(CIdds) >= 10) >= 10
#lowreads <- which(colSums(counts(CIdds)) < 500)
CIdds <- CIdds[keep,]
#CIdds <- CIdds[-c(lowreads),]

#Make DEseq
CIdds <- DESeq(CIdds, BPPARAM = 42)

CIresults <- results(CIdds, alpha = 0.05, contrast = c("Disease_State", "BRD", "Healthy"))

vsd <- vst(CIdds, blind = F)


### make pCA
#vsd <- varianceStabilizingTransformation(CIdds, blind = F)
pca <- plotPCA(vsd,ntop=10000, intgroup=c("Location", "Disease_State"), pcsToUse = c(1,2)) + geom_point(alpha = .7) +
  theme_cowplot() + ggtitle("Cow Gene Expression PCA - No Batch Correction") + scale_color_brewer(palette = "Paired") #stat_ellipse() + theme(plot.title = element_text(face = "plain"))  
  #scale_color_manual(values = c("BRD" = "firebrick1", "Healthy"= "#4682B4"))
pca


### Remove batch Effects
mat <- vsd
m <- assay(vsd)
#mm <- model.matrix(~Disease_State + Location + Seqtype, colData(vsd))
new <- limma::removeBatchEffect(m, batch = CI_table$Seqtype, covariates = model.matrix(~Location, CI_table)[,-1])
assay(mat) <- new
assay(mat)
batch <- plotPCA(mat, ntop=10000, intgroup=c("Location", "Disease_State"))+ geom_point(alpha = .7) +
  theme_cowplot() + ggtitle("Cow Gene Expression PCA - Batch Corrected") + scale_color_brewer(palette = "Paired")

#make plot for supp info.
corrected <- pca +  batch + patchwork::plot_layout(ncol = 2, guides = "collect") + patchwork::plot_annotation(tag_levels = "A") 

ggsave("/mnt/md0/brd_project/figures/revised_figures/Batch_correction_PCA.pdf", plot = corrected, width = 14, height = 6)
ggsave("/mnt/md0/brd_project/figures/revised_figures/Batch_correction_PCA.png", plot = corrected, width = 14, height = 6, dpi = 300)


summary(CIresults)

#subset significant genes
sig_res <- subset(CIresults, padj <= 0.05)
#write.csv(sig_res, file = "/mnt/md0/brd_project/cow_immune/significant_DE_genes.csv", quote = F)
sig_up <- subset(sig_res, log2FoldChange > 0)
sig_down <- subset(sig_res, log2FoldChange < 0)

#make new column for diff expression
CIresults$diffexpressed <- NA
CIresults$diffexpressed[which((rownames(CIresults) %in% rownames(sig_up)))] <- "Upregulated"
CIresults$diffexpressed[which((rownames(CIresults) %in% rownames(sig_down)))] <- "Downregulated"
CIresults$diffexpressed <- CIresults$diffexpressed %>% tidyr::replace_na("Not Significant")
CIresults$diffexpressed <- factor(CIresults$diffexpressed, levels = c("Upregulated", "Downregulated", "Not Significant"))

#make volcano plot
#colors
col_scheme = c("Upregulated" = "firebrick1", "Downregulated"= "#4682B4", "Not Significant" = "grey")
#make plot
vol_plot <- ggplot(CIresults, aes(x=log2FoldChange, y=-log(padj,10))) +
  geom_point(aes(color = diffexpressed),size = 2.5, alpha = 0.6) + 
  scale_color_manual(values = col_scheme) +
  ggtitle("BRD vs. Healthy") + ylab("-log10padj") + xlab("Log2FoldChange") +
  theme_cowplot() +theme(plot.title = element_text(face = "plain")) +
  guides(color=guide_legend(title="Diff. Expression")) +
  geom_hline(yintercept = -log(0.05,10), linetype = "dashed", color = "grey")  
#call plot
vol_plot


#get genes for gene set enrichment analysis
gene_list <- CIresults$log2FoldChange
names(gene_list) <- rownames(CIresults)
gene_list_sorted <- sort(gene_list, decreasing = T)

gse <- gseGO(geneList=gene_list_sorted, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE,
             eps = 0,
             OrgDb = org.Bt.eg.db, 
             pAdjustMethod = "BH")

require(DOSE)
gse@result$Description <- gsub("biological process involved in ", "", gse@result$Description)
gse@result$Description <- gsub(" between organisms", "", gse@result$Description)
GO_DEGS <- clusterProfiler::dotplot(gse, showCategory=20, title = "Enriched GO Terms - DEGs", font.size = 9) 
#emapplot(gse, showCategory=10)
#cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 50)
#ridgeplot(gse) + labs(x = "enrichment distribution")

#### KEGG

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-clusterProfiler::bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = CIresults[rownames(CIresults) %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)


kegg_organism = "bta"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")

t15_pathways <- clusterProfiler::dotplot(kk2, showCategory = 15, title = "Enriched Pathways - DEGs", 
                                         font.size = 9.5)
#png(filename = "/mnt/md0/brd_project/cow_immune/top15_enriched_pathways.png", width = 10, height = 8, res = 300, units = "in")
t15_pathways
#dev.off()


#Make GO figure
GO_immune <- gse %>% filter(p.adjust < 0.05) %>%
  ggplot(showCategory = 20,
         aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(size = Count)) +
  #scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(1, 3)) +
  #scale_color_manual(values = "green") +
  theme_minimal() + 
  theme(axis.text = element_text(size = 7.5)) +
  xlab("Gene Ratio") +
  ylab(NULL) + 
  ggtitle("GO Enrichment - DEGs")


#make KEGG Figure
KEGG_immune <- kk2 %>% filter(p.adjust < 0.05) %>%
  ggplot(showCategory = 20,
         aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(size = Count)) +
  #scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(1, 3)) +
  #scale_color_manual(values = "green") +
  theme_minimal() + 
  theme(axis.text = element_text(size = 9)) +
  xlab("Gene Ratio") +
  ylab(NULL) + 
  ggtitle("KEGG Enrichment - Pathways")

# GO_immune
# KEGG_immune

des <- 
"
AB
AC
"

###make host transcriptome final figure
#im_fig <- vol_plot + pca + GO_immune  + KEGG_immune + patchwork::plot_annotation(tag_levels = "A")
im_fig <- vol_plot | (GO_immune / KEGG_immune) + patchwork::plot_annotation(tag_levels = "A") #+ plot_layout(design = des)
#im_fig <- vol_plot + GO_immune + KEGG_immune + patchwork::plot_annotation(tag_levels = "A") + plot_layout(design = des, wid)

im_fig
#ggsave(plot = im_fig, filename =  "/mnt/md0/brd_project/figures/revised_figures/immune_fig_allSamples_REVISED.pdf", width = 11, height = 6, units = "in", dpi=300)
#ggsave(plot = im_fig, filename =  "/mnt/md0/brd_project/figures/revised_figures/immune_fig_allSamples_REVISED.png", width = 11, height = 6, units = "in", dpi=300)

godf <- data.frame(gse@result)
keggdf <- data.frame(kk2@result)


deglist <- list("DEGs" = data.frame(sig_res),
    "GO Enrichment" = godf,
     "KEGG Enrichment" = keggdf)
#save supplementary information
openxlsx::write.xlsx(deglist, "/mnt/md0/brd_project/SupplementaryTable_DEGs.xlsx")

############## WGCNA #############

#get batch-corrected counts and filter them based on variance
expr <- t(new)
gene_var <- colVars(expr)
keep <- gene_var > quantile(gene_var, 0.75)
expr_filt <- expr[, keep]

library(WGCNA)

sampleTree <- hclust(
  dist(expr_filt),
  method = "average"
)

plot(sampleTree)

powers <- c(1:50)

sft <- pickSoftThreshold(
  expr_filt,
  powerVector = powers,
  networkType = "signed"
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.80, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 26
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(expr_filt,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)


# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)


# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(expr_filt, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")


# pick out a few modules of interest here
modules_of_interest = c("grey", "turquoise", "blue", "yellow", "brown", "green")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id


# Get normalized expression for those genes
expr_filt[1:5,1:10]
texp <- t(expr_filt)
subexpr = texp[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")


#### CHECK FOR GO TERM ENRICHMENT IN THE WGCNA GROUPS
greygenes <- subset(submod, colors == "grey")$gene_id
greygo <- enrichGO(greygenes, 
             ont ="BP", 
             keyType = "SYMBOL", 
             #nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             #verbose = TRUE,
             #eps = 0,
             OrgDb = org.Bt.eg.db, 
             pAdjustMethod = "BH")


turqgenes <- subset(submod, colors == "turquoise")$gene_id 
turqgo <- enrichGO(turqgenes, 
                   ont ="BP", 
                   keyType = "SYMBOL", 
                   #nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   #verbose = TRUE,
                   #eps = 0,
                   OrgDb = org.Bt.eg.db, 
                   pAdjustMethod = "BH")

bluegenes <- subset(submod, colors == "blue")$gene_id 
bluego <- enrichGO(bluegenes, 
                   ont ="BP", 
                   keyType = "SYMBOL", 
                   #nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   #verbose = TRUE,
                   #eps = 0,
                   OrgDb = org.Bt.eg.db, 
                   pAdjustMethod = "BH")


yellowgenes <- subset(submod, colors == "yellow")$gene_id 
yellowgo <- enrichGO(yellowgenes, 
                   ont ="BP", 
                   keyType = "SYMBOL", 
                   #nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   #verbose = TRUE,
                   #eps = 0,
                   OrgDb = org.Bt.eg.db, 
                   pAdjustMethod = "BH")

browngenes <- subset(submod, colors == "brown")$gene_id 
browngo <- enrichGO(browngenes, 
                   ont ="BP", 
                   keyType = "SYMBOL", 
                   #nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   #verbose = TRUE,
                   #eps = 0,
                   OrgDb = org.Bt.eg.db, 
                   pAdjustMethod = "BH")

greengenes <- subset(submod, colors == "green")$gene_id 
greengo <- enrichGO(greengenes, 
                   ont ="BP", 
                   keyType = "SYMBOL", 
                   #nPerm = 10000, 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   #verbose = TRUE,
                   #eps = 0,
                   OrgDb = org.Bt.eg.db, 
                   pAdjustMethod = "BH")

View(greengo@result) # no sig
View(bluego@result) # no sig - immune
View(yellowgo@result) # protein glycosylation
View(browngo@result) #cellular respiration
View(turqgo@result) # Immune
View(yellowgo@result) # Viral immunity

# save as supplemental data
Aen <- data.frame(turqgo@result)
Ben <- data.frame(browngo@result)
Cen <- data.frame(bluego@result)
Den <- data.frame(greengo@result)
Een <- data.frame(yellowgo@result)
Uen <- data.frame(yellowgo@result)

A_ids <- data.frame("GeneName" = mapIds(org.Bt.eg.db, turqgenes, 'GENENAME', 'SYMBOL'))
B_ids <- data.frame("GeneName" = mapIds(org.Bt.eg.db, browngenes, 'GENENAME', 'SYMBOL'))
C_ids <- data.frame("GeneName" = mapIds(org.Bt.eg.db, bluegenes, 'GENENAME', 'SYMBOL'))
D_ids <- data.frame("GeneName" = mapIds(org.Bt.eg.db, greengenes, 'GENENAME', 'SYMBOL'))
E_ids <- data.frame("GeneName" = mapIds(org.Bt.eg.db, yellowgenes, 'GENENAME', 'SYMBOL'))
Un_ids <- data.frame("GeneName" = mapIds(org.Bt.eg.db, greygenes, 'GENENAME', 'SYMBOL'))

datalist <- list("ModuleA" = A_ids,
                 "ModuleB" = B_ids,
                 "ModuleC" = C_ids,
                 "ModuleD" = D_ids,
                 "ModuleE" = E_ids,
                 "Unclassified" = Un_ids,
                 "ModuleA Enrich" = Aen,
                 "ModuleB Enrich" = Ben,
                 "ModuleC Enrich" = Cen,
                 "ModuleD Enrich" = Den,
                 "ModuleE Enrich" = Een,
                 "Unclassified Enrich" = Uen)

openxlsx::write.xlsx(datalist, file = "/mnt/md0/brd_project/Supplementary_table_WGCNA.xlsx", rowNames = T)

## Turqoise group is enriched with viral GO terms.

#Correlate to other viruses/viral abundance
viral_rpkm <- "/mnt/md0/brd_project/viral_RPKMs_REVISED.tsv"
bac_abun <- "/mnt/md0/brd_project/bacteria/standard_kraken/bacteria_abundances_REVISED.csv"
rpkms <- read.csv(viral_rpkm)
ba <- read.csv(bac_abun)
bal <- ba %>% dplyr::select(SRA, Abundance, Genus) %>% pivot_wider(id_cols = SRA, values_from = "Abundance", names_from = "Genus")
data.table::setDT(rpkms)
rpkmt <- data.table::transpose(rpkms, make.names = "X", keep.names = "SRA")
rpkmt$Viral_abundance <- rowSums(rpkmt[,2:ncol(rpkmt)])
rpkmt <- rpkmt[-which(rpkmt$SRA %in% ols),]
alphadiv <- "/mnt/md0/brd_project/bacteria/standard_kraken/shannon_diversity_REVISED.csv"
adiv <- read.csv(alphadiv)
adiv <- adiv[,2:3]
names(adiv)[1] <- "SRA"

dat <- data.frame(
  SRA = rownames(MEs0),
  MEturquoise = MEs0$MEturquoise,
  MEgrey = MEs0$MEgrey,
  MEgreen = MEs0$MEgreen,
  MEbrown = MEs0$MEbrown,
  MEblue = MEs0$MEblue,
  MEyellow = MEs0$MEyellow
) |>
  dplyr::inner_join(rpkmt, by = "SRA") %>%
  dplyr::inner_join(CI_table, by = "SRA") %>%
  dplyr::inner_join(bal, by = "SRA") %>%
  dplyr::inner_join(adiv, by = "SRA")
dat$Disease_State <- factor(dat$Disease_State, levels = c("BRD", "Healthy"))

dat <- subset(dat, Location != "USA.Kansas")

### Which modules are correlated to disease state?
ME_results <- lapply(
  c(colnames(dat)[2:7]),
  function(v) {

    fite <- lm(
     dat[,v] ~ Disease_State +
        Seqtype +
        Location,
      data = dat
    )
    # fitv <- residuals(lm(
    #   dat[[v]] ~ Disease_State +
    #     Seqtype +
    #     Location,
    #   data = dat
    # ))
    print(summary(fite))
    data.frame(
      ME = v,
      p = broom::tidy(summary(fite))$p.value[[2]],
      beta = broom::tidy(summary(fite))$estimate[[2]]
    )
  }
)

ME_results <- bind_rows(ME_results)

ME_results$FDR <-
  p.adjust(ME_results$p, "BH")

lm_brown <- lm(
  MEbrown ~ Disease_State +
    Seqtype +
    Location,
  data = dat
)
summary(lm_brown)

lm_grey <- lm(
  MEgrey ~ Disease_State +
    Seqtype +
    Location,
  data = dat
)
summary(lm_grey)

lm_yellow <- lm(
  MEyellow ~ Disease_State +
    Seqtype +
    Location,
  data = dat
)
summary(lm_yellow)

lm_turqoise <- lm(
  MEturquoise ~ Disease_State +
    Shannon +
    IDV +
    BHV1 +
    Mycoplasmopsis +
    Pasteurella +
    Histophilus +
    Mannheimia +
    Moraxella +
    Seqtype +
    Location,
  data = dat
)
bturq <- broom::tidy(summary(lm_turqoise))


# lm_grey <- lm(
#   MEgrey ~ Disease_State +
#     Seqtype +
#     Location,
#   data = dat
# )
# summary(lm_grey)


lm_grey <- lm(
  MEgrey ~ Disease_State +
    Shannon +
    IDV +
    BHV1 +
    Mycoplasmopsis +
    Pasteurella +
    Histophilus +
    Mannheimia +
    Moraxella +
    Seqtype +
    Location,
  data = dat
)
bgray <- broom::tidy(summary(lm_grey))

lm_yellow <- lm(
  MEyellow ~ Disease_State +
    Shannon +
    IDV +
    BHV1 +
    Mycoplasmopsis +
    Pasteurella +
    Histophilus +
    Mannheimia +
    Moraxella +
    Seqtype +
    Location,
  data = dat
)
byel <- broom::tidy(summary(lm_yellow))


lm_brown <- lm(
  MEbrown ~ Disease_State +
    Shannon +
    IDV +
    BHV1 +
    Mycoplasmopsis +
    Pasteurella +
    Histophilus +
    Mannheimia +
    Moraxella + 
    Seqtype +
    Location,
  data = dat
)
bbrow <- broom::tidy(summary(lm_brown))



lm_blue <- lm(
  MEblue ~ Disease_State +
    Shannon +
    IDV +
    BHV1 +
    Mycoplasmopsis +
    Pasteurella +
    Histophilus +
    Mannheimia +
    Moraxella + 
    Seqtype +
    Location,
  data = dat
)
bblue <- broom::tidy(summary(lm_blue))


lm_green <- lm(
  MEgreen ~ Disease_State +
    Shannon +
    IDV +
    BHV1 +
    Mycoplasmopsis +
    Pasteurella +
    Histophilus +
    Mannheimia +
    Moraxella + 
    Seqtype +
    Location,
  data = dat
)
bgre <- broom::tidy(summary(lm_green))

ripA <- relaimpo::calc.relimp(lm_turqoise,type = "lmg")
ripC <-relaimpo::calc.relimp(lm_blue,type = "lmg")
ripB <-relaimpo::calc.relimp(lm_brown,type = "lmg")
ripE <-relaimpo::calc.relimp(lm_yellow,type = "lmg")
ripD <-relaimpo::calc.relimp(lm_green,type = "lmg")
ripU <-relaimpo::calc.relimp(lm_grey,type = "lmg")

riplist <- list(ModuleA = ripA@lmg, ModuleB = ripB@lmg,
            ModuleC = ripC@lmg, ModuleD = ripD@lmg,
            ModuleE = ripE@lmg, Unclassified = ripU@lmg)

riplist <- list(ripA@lmg, ripB@lmg,
                ripC@lmg, ripD@lmg,
                ripE@lmg, ripU@lmg)

ripdf <- bind_rows(riplist, .id = "id")

blis <- list(ModuleA = data.frame(bturq), ModuleB = data.frame(bbrow), 
             ModuleC = data.frame(bblue), ModuleD = data.frame(bgre), ModuleE = data.frame(byel), Unclassified = data.frame(bgray))

for (x in 1:length(blis)){
  blis[[x]]$FDR <- p.adjust(blis[[x]][,5], method = "BH")
}

blis
openxlsx::write.xlsx(blis, file = "/mnt/md0/brd_project/supplementary_table_Table3Stats.xlsx")

########## Raw correlations
virus_results <- lapply(
  c(colnames(rpkmt)[2:16]),
  function(v) {
    
    ct <- cor.test(
      dat$MEturquoise,
      dat[,v],
      method = "spearman",
      exact = FALSE
    )
    
    data.frame(
      Virus = v,
      p = ct$p.value,
      rho = ct$estimate[[1]]
    )
  }
)

virus_results <- bind_rows(virus_results)

virus_results$FDR <-
  p.adjust(virus_results$p, "BH")

brbv_lm <- lm(
  MEbrown ~ BrBv +
    Disease_State +
    Seqtype +
    Location,
  data = dat
)

summary(brbv_lm)

beve_lm <- lm(
  MEblue ~ IDV +
    Disease_State +
    Seqtype +
    Location,
  data = dat
)
summary(beve_lm)


idv_lm <- lm(
  MEblue ~ IDV +
    Disease_State +
    Seqtype +
    Location,
  data = dat
)
summary(idv_lm)

va_lm <- lm(
  MEturquoise ~ BEVE +
    Disease_State +
    Seqtype +
    Location,
  data = dat
)
summary(va_lm)


bcov_lm <- lm(
  MEgrey ~ BCoV +
    Disease_State +
    Seqtype +
    Location,
  data = dat
)
summary(bcov_lm)

bcov_lmy <- lm(
  MEyellow ~ BPV +
    Disease_State +
    Seqtype +
    Location,
  data = dat
)
summary(bcov_lmy)


#correlate to alpha diversity
adcor <- cor.test(dat$MEblue, dat$Shannon, method = "spearman")
adcor

ad_lm <- lm(MEyellow ~ Shannon +
                Disease_State,
                #Seqtype +
                #Location,
              data = subset(dat, Location == "Australia"))

ad_lm <- lm(MEturquoise ~ Shannon +
              Disease_State,
              #Seqtype +
              #Location,
            data = subset(dat, Seqtype == "Nanopore"))
summary(ad_lm)

# bacdat <- na.omit(dat)
# 
# bacteria_results <- lapply(
#   c(names(bacdat)[21:ncol(bacdat)]),
#   function(v) {
#     
#     fite <- residuals(lm(
#       MEturquoise ~ bacdat$Disease_State +
#         bacdat$Seqtype +
#         bacdat$Location,
#       data = bacdat
#     ))
#     fitv <- residuals(lm(
#       bacdat[[v]] ~ Disease_State +
#         bacdat$Seqtype +
#         bacdat$Location,
#       data = bacdat
#     ))
#     
#     ct <- cor.test(
#       fite,
#       fitv,
#       method = "spearman"
#     )
#     
#     data.frame(
#       bacteria = v,
#       p = ct$p.value,
#       rho = ct$estimate[[1]]
#     )
#   }
# )
# 
# bacteria_results <- bind_rows(bacteria_results)
# 
# bacteria_results$FDR <-
#   p.adjust(bacteria_results$p, "BH")



