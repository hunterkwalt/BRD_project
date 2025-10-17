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
library(pathview)
library(edgeR)
library(microbiome)
library(tidyplots)
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
CI <- data.frame(CI[-which(CI$V1 %in% ctl$V1),])
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
txdb <- makeTxDbFromGFF(file = gff)
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

vsd <- vst(CIdds)
#vsd <- varianceStabilizingTransformation(CIdds, blind = F)
pca <- plotPCA(vsd,ntop=1000, intgroup=c("Disease_State"), pcsToUse = c(1,2)) + geom_point(alpha = .7) +
  theme_cowplot() + ggtitle("Cow Gene Expression PCA") + stat_ellipse() + theme(plot.title = element_text(face = "plain")) + 
  scale_color_manual(values = c("BRD" = "firebrick1", "Healthy"= "#4682B4"))
pca


mat <- vsd
m <- assay(vsd)
#mm <- model.matrix(~Disease_State + Location + Seqtype, colData(vsd))
new <- limma::removeBatchEffect(m, batch1 = CI_table$Seqtype, batch2 = CI_table$Location)
assay(mat) <- new
plotPCA(mat, ntop=1000, intgroup=c("Location", "Disease_State"))

#ggsave("/mnt/md0/brd_project/figures/new_figures/expression_PCA.pdf", plot = pca)

#library("pheatmap")
# select <- order(rowMeans(counts(CIdds,normalized=TRUE)),
#                 decreasing=TRUE)[1:50]
# df <- as.data.frame(colData(CIdds)[,c("Disease_State")])
# hm = pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=FALSE, show_colnames = F,
#          cluster_cols=T, annotation_col=df)
# 
# #ggsave("heatmap_DE.pdf", plot = hm)
# 
# #png("../figures/deheatmp.png", width = 6, height = 8.5, units = "in", res = 300)
# hm
# #dev.off()
# #summarize results
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




# igs <- read.table("../cow_immune/immune_genes.txt", header = F)
# 
# idxs <- c()
# for (i in 1:length(igs$V1)){
#   ids <- grep(igs$V1[i], row.names(sig_res))
#   idxs <- c(idxs, ids)
# }
# 
# immune <- data.frame(sig_res[idxs,])
#write.csv(immune, "../cow_immune/immune_subset.csv", quote = F)

#grepl(igs$V1, row.names(sig_res))

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
ids<-bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
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

im_fig <- vol_plot + pca + GO_immune  + KEGG_immune + patchwork::plot_annotation(tag_levels = "A")
im_fig <- vol_plot | (GO_immune / KEGG_immune) + patchwork::plot_annotation(tag_levels = "A") #+ plot_layout(widths = c(2,1))
im_fig
ggsave(plot = im_fig, filename =  "/mnt/md0/brd_project/figures/new_figures/immune_fig_allSamples_filt250.pdf", width = 11, height = 6, units = "in", dpi=300)

