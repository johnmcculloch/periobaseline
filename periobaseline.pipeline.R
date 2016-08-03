#####Pipeline for alpha and beta analysis
#####By John McCulloch
#####Version 0.9
#####October 2015

####1. Load all necessary dependencies
library("phyloseq")
library("metagenomeSeq")
library("ggplot2")
library("plyr")
library("DESeq2")

####2. Prepare the data
###2.1 Import data from Mothur
mothur <- import_mothur(mothur_shared_file = "periobaseline.an.shared", mothur_constaxonomy_file = "periobaseline.an.cons.taxonomy", parseFunction = parse_taxonomy_default)

###2.2 create phyloseq_taxtable
periobaseline.phyloseq.taxtable <- cbind.data.frame(tax_table(mothur))
periobaseline.phyloseq.taxtable$Rank1<-as.character(periobaseline.phyloseq.taxtable$Rank1)
periobaseline.phyloseq.taxtable$Rank2<-as.character(periobaseline.phyloseq.taxtable$Rank2)
periobaseline.phyloseq.taxtable$Rank3<-as.character(periobaseline.phyloseq.taxtable$Rank3)
periobaseline.phyloseq.taxtable$Rank4<-as.character(periobaseline.phyloseq.taxtable$Rank4)
periobaseline.phyloseq.taxtable$Rank5<-as.character(periobaseline.phyloseq.taxtable$Rank5)
periobaseline.phyloseq.taxtable$Rank6<-as.character(periobaseline.phyloseq.taxtable$Rank6)
periobaseline.phyloseq.taxtable$Rank7<-as.character(periobaseline.phyloseq.taxtable$Rank7)
colnames(periobaseline.phyloseq.taxtable)<-c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

###2.3 create OTU table
periobaseline.otutable <- otu_table(mothur)
otutable<-cbind.data.frame(periobaseline.otutable[,1:ncol(periobaseline.otutable)])

###2.4 create Pheno table
periobaseline.pheno <- read.table(file="periobaseline.pheno.tsv", header=T)
rownames(periobaseline.pheno) <- periobaseline.pheno$Sample
periobaseline.pheno$Sample <- NULL

###2.5 Create a non-normalised class object in phyloseq (periobaseline_physeq)

OT <- as.matrix(otutable)
OTUphyloseq = otu_table(OT, taxa_are_rows = TRUE)

TT <- as.matrix(periobaseline.phyloseq.taxtable)
TAXphyloseq = tax_table(TT)

PTphyloseq = sample_data(periobaseline.pheno)

periobaseline_physeq = phyloseq(OTUphyloseq, TAXphyloseq, PTphyloseq)


###2.6 Create a non-normalised class object in metagenomeseq (periobaseline_mgseq)
phenotypeData = AnnotatedDataFrame(periobaseline.pheno)
OTUdata = AnnotatedDataFrame(periobaseline.phyloseq.taxtable)
periobaseline_mgseq = newMRexperiment(otutable, phenoData=phenotypeData, featureData=OTUdata)

###2.7 Create a normalised class object in metagenomeseq (periobaseline_mgseq_norm)

p = cumNormStatFast(periobaseline_mgseq) 
periobaseline_mgseq_norm = cumNorm(periobaseline_mgseq, p = p)
periobaseline.normalized.matrix <- MRcounts(periobaseline_mgseq_norm, norm = TRUE)

###2.8 Create a normalised class object in phyloseq (periobaseline_physeq_norm)
OTnorm <- as.matrix(periobaseline.normalized.matrix)
OTUphyloseqnorm = otu_table(OTnorm, taxa_are_rows = TRUE)

TT <- as.matrix(periobaseline.phyloseq.taxtable)
TAXphyloseq = tax_table(TT)
PTphyloseq = sample_data(periobaseline.pheno)

periobaseline_physeq_norm = phyloseq(OTUphyloseqnorm, TAXphyloseq, PTphyloseq)

###2.9 Create non-normalised class object in metagenomeseq with descriptive OTU names (periobaseline_mgseq_ren)

OTU_names <- rownames(periobaseline.phyloseq.taxtable)
OTU_names <- paste(OTU_names, periobaseline.phyloseq.taxtable$Genus, periobaseline.phyloseq.taxtable$Species, sep="-")
periobaseline.mgseq.taxtable <- periobaseline.phyloseq.taxtable
rownames(periobaseline.mgseq.taxtable) <- OTU_names

otutable_mgseq<-cbind.data.frame(otutable[,1:ncol(otutable)])
rownames(otutable_mgseq) <- OTU_names

#Create MRexperiment in metagenomeseq
OTUdata_ren = AnnotatedDataFrame(periobaseline.mgseq.taxtable)
periobaseline_mgseq_ren = newMRexperiment(otutable_mgseq, phenoData=phenotypeData, featureData=OTUdata_ren)

####3. Calculate alpha diversity measures
#Use periobaseline_physeq
#set theme and colours
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <- function(palname = pal, ...) {
    scale_colour_brewer(palette = palname, ...)
}
scale_fill_discrete <- function(palname = pal, ...) {
    scale_fill_brewer(palette = palname, ...)
}
###3.1 Calculate richness measures per patient
patient.obs <- plot_richness(periobaseline_physeq, x="Patient", color="abs_status", measures = c("Observed"))
patient.chao <- plot_richness(periobaseline_physeq, x="Patient", color="abs_status", measures = c("Chao1"))
patient.invsimpson <- plot_richness(periobaseline_physeq, x="Patient", color="abs_status", measures = c("InvSimpson"))
patient.shannon <- plot_richness(periobaseline_physeq, x="Patient", color="abs_status", measures = c("Shannon"))


#save results
ggsave(patient.obs, file="alpha_patient_observed.pdf", units="cm", width=29, height=21)
ggsave(patient.chao, file="alpha_patient_chao.pdf", units="cm", width=29, height=21)
ggsave(patient.invsimpson, file="alpha_patient_invsimpson.pdf", units="cm", width=29, height=21)
ggsave(patient.shannon, file="alpha_patient_shannon.pdf", units="cm", width=29, height=21)

###3.2 Calculate richness measures by type of sample
periobaseline_physeq.abs_status = merge_samples(periobaseline_physeq, "abs_status")
# repair variables that were damaged during merge (coerced to numeric)
sample_data(periobaseline_physeq.abs_status)$abs_status <- factor(sample_names(periobaseline_physeq.abs_status))

alpha_sampletype =  plot_richness(periobaseline_physeq.abs_status, color = "abs_status", measures = c("Chao1", "InvSimpson"))
alpha_sampletype + geom_point(size = 5, alpha = 0.7)
ggsave(alpha_sampletype, file="alpha_sampletype.pdf", units="cm", width=29, height=21)

####4. Draw up abundance barplots
#use periobaseline_physeq_norm
##transform to relabund
periobaseline_physeq_norm_relabund = transform_sample_counts(periobaseline_physeq_norm, function(x) x / sum(x))
periobaseline_physeq_norm_filter = filter_taxa(periobaseline_physeq_norm_relabund, function(x) mean(x) > 1e-5, TRUE)

##plot abund
title.phyla = "Abundance of Phyla per type of site"
abund.phyla = plot_bar(periobaseline_physeq_norm_filter, "Phylum", "Abundance", title=title.phyla, facet_grid="abs_status~.")
ggsave(abund.phyla, file="periobaseline_phylum_relabund.pdf", units="cm", width=29, height=21)

title.genus = "Abundance of Genera per type of site"
abund.genus = plot_bar(periobaseline_physeq_norm_filter, "Genus", "Abundance", title=title.genus, facet_grid="abs_status~.")
ggsave(abund.genus, file="periobaseline_genus_relabund.pdf", units="cm", width=29, height=21)

####5. Ordination plots
##use periobaseline_physeq_norm
###5.1 ordination plot of all individuals in the same plot
physeq.ord.bray.nmds <- ordinate(periobaseline_physeq_norm, "NMDS", "bray")
plot.physeq.ord.bray.nmds = plot_ordination(periobaseline_physeq_norm, physeq.ord.bray.nmds, type = "samples", color = "abs_status", shape = "Health")
ggtitle("samples")
ggsave(plot.physeq.ord.bray.nmds, file="periobaseline_ord_nmds_bray.pdf", units="cm", width=29, height=21)

physeq.ord.jaccard.nmds <- ordinate(periobaseline_physeq_norm, "NMDS", "jaccard")
plot.physeq.ord.jaccard.nmds = plot_ordination(periobaseline_physeq_norm, physeq.ord.jaccard.nmds, type = "samples", color = "abs_status", shape = "Health")
ggtitle("samples")
ggsave(plot.physeq.ord.jaccard.nmds, file="periobaseline_ord_nmds_jaccard.pdf", units="cm", width=29, height=21)

physeq.ord.bray.mds <- ordinate(periobaseline_physeq_norm, "MDS", "bray")
plot.physeq.ord.bray.mds = plot_ordination(periobaseline_physeq_norm, physeq.ord.bray.mds, type = "samples", color = "abs_status", shape = "Health")
ggtitle("Bray MDS")
ggsave(plot.physeq.ord.bray.mds, file="periobaseline_ord_pcoa_bray.pdf", units="cm", width=29, height=21)

physeq.ord.jaccard.mds <- ordinate(periobaseline_physeq_norm, "MDS", "jaccard")
plot.physeq.ord.jaccard.mds = plot_ordination(periobaseline_physeq_norm, physeq.ord.jaccard.mds, type = "samples", color = "abs_status", shape = "Health")
ggtitle("samples")
ggsave(plot.physeq.ord.jaccard.mds, file="periobaseline_ord_pcoa_jaccard.pdf", units="cm", width=29, height=21)

###5.2 Ordination of each individual separately
##First, separate each patient
## This is a bad idea

####6. Plots
###6.1 Abundance heatmap using metagenomeseq
###Use periobaseline_mgseq_ren

# plotMRheatmap of top 30 OTUs with highest variance
pdf(file="periobaseline_mrheat_health.pdf")
par(oma=c(5,4,4,5)+0.1)
trials = pData(periobaseline_mgseq_ren)$abs_status
heatmapColColors = brewer.pal(12,"Set3")[as.integer(factor(trials))]
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
plotMRheatmap(obj = periobaseline_mgseq_ren, n = 30, norm = TRUE, log = TRUE, fun = sd, cexRow = 0.4, cexCol = 0.3, trace = "none", col = heatmapCols, ColSideColors = heatmapColColors)
cores <- c("#8DD3C7", "#FB8072", "#BEBADA", "#FFFFB3")
legend(x=0.005, y=0.90, legend=c("Control","Shallow","Intermediate","Deep"), fill=cores, border=TRUE, bty="n", y.intersp = 0.7, cex=0.35)
dev.off()

# plotMRheatmap of top 107 OTUs with highest variance
pdf(file="periobaseline_mrheat_health_107.pdf")
par(oma=c(5,4,4,5)+0.1)
trials = pData(periobaseline_mgseq_ren)$abs_status
heatmapColColors = brewer.pal(12,"Set3")[as.integer(factor(trials))]
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
plotMRheatmap(obj = periobaseline_mgseq_ren, n = 107, norm = TRUE, log = TRUE, fun = sd, cexRow = 0.4, cexCol = 0.3, trace = "none", dendrogram = "none", col = heatmapCols, ColSideColors = heatmapColColors)
cores <- c("#8DD3C7", "#FB8072", "#BEBADA", "#FFFFB3")
legend(x=0.005, y=0.90, legend=c("Control","Shallow","Intermediate","Deep"), fill=cores, border=TRUE, bty="n", y.intersp = 0.7, cex=0.35)
dev.off()



###6.2 Correlation plot
pdf(file="periobaseline_correlation_plot.pdf")
par(mar=c(5,4,4,5)+0.1)
plotCorr(obj = periobaseline_mgseq_ren, n = 30, norm = TRUE, log = TRUE, fun = cor, dendrogram="none", cexRow = 0.25, cexCol = 0.25, trace = "none", col = heatmapCols)
dev.off()


####7. PA using Deseq2
###Use periobaseline_physeq_norm

library("DESeq2")
#Get foldchange table
diagdds = phyloseq_to_deseq2(periobaseline_physeq_norm, ~ Health)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(periobaseline_physeq_norm)[rownames(sigtab), ], "matrix"))
head(sigtab)

#Export sigtab as tsv
write.table(sigtab, file="periobaseline_significant_foldchange.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

#Plot results
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {scale_fill_brewer(palette = palname, ...)}

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
foldchangenorm<-ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#save
ggsave(foldchangenorm, file="periobaseline_significant_foldchange_plot.pdf", units="cm", width=29, height=21)


### TO DO: regenerate PA testing using only patients.
nondentists <- subset_samples(periobaseline_physeq_norm, Profession=="nondentist")
periobaseline_physeq_onlypatients <- subset_samples(nondentists, Patient!="4Ju")

library("DESeq2")
#Get foldchange table
diagdds_OP = phyloseq_to_deseq2(periobaseline_physeq_onlypatients, ~ Health)
diagdds_OP = DESeq(diagdds_OP, test="Wald", fitType="parametric")
res_OP = results(diagdds_OP, cooksCutoff = FALSE)
alpha = 0.05
sigtab_OP = res_OP[which(res_OP$padj < alpha), ]
sigtab_OP = cbind(as(sigtab_OP, "data.frame"), as(tax_table(periobaseline_physeq_onlypatients)[rownames(sigtab_OP), ], "matrix"))
head(sigtab_OP)

#Export sigtab as tsv
write.table(sigtab_OP, file="periobaseline_significant_foldchange_onlypatients.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

#Plot results
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {scale_fill_brewer(palette = palname, ...)}

# Genus order
x = tapply(sigtab_OP$log2FoldChange, sigtab_OP$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_OP$Genus = factor(as.character(sigtab_OP$Genus), levels=names(x))
foldchangenorm_OP<-ggplot(sigtab_OP, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#save
ggsave(foldchangenorm_OP, file="periobaseline_significant_foldchange_onlypatients_plot.pdf", units="cm", width=29, height=21)







#### Create an MR table for each metagenome with rows for each phylum or Genus. Rows as metagenomes and columns as taxa, preferentially.


###Agglomerate OTUs into phyla
periobaseline_physeq_norm_phylum = tax_glom(periobaseline_physeq_norm, "Phylum")
#Extract data frame containing OTU counts
periobaseline_physeq_norm_phylum_otutable <-  as.data.frame(otu_table(periobaseline_physeq_norm_phylum))
#Rename rows
rownames(periobaseline_physeq_norm_phylum_otutable) <- as.data.frame(tax_table(periobaseline_physeq_norm_phylum))$Phylum
#Round to entire sequences
periobaseline_physeq_norm_phylum_otutable <- round(periobaseline_physeq_norm_phylum_otutable, 0)
#Transpose
periobaseline_physeq_norm_phylum_otutable <- t(periobaseline_physeq_norm_phylum_otutable)
#Export
write.table(periobaseline_physeq_norm_phylum_otutable, file="periobaseline_physeq_norm_phylum_abundances.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

###Agglomerate OTUs into genera
periobaseline_physeq_norm_genera = tax_glom(periobaseline_physeq_norm, "Genus")
#Extract data frame containing OTU counts
periobaseline_physeq_norm_genera_otutable <-  as.data.frame(otu_table(periobaseline_physeq_norm_genera))
#Rename rows
rownames(periobaseline_physeq_norm_genera_otutable) <- make.unique(as.character(as.data.frame(tax_table(periobaseline_physeq_norm_genera))$Genus))
#Round to entire sequences
periobaseline_physeq_norm_genera_otutable <- round(periobaseline_physeq_norm_genera_otutable, 0)
#Transpose
periobaseline_physeq_norm_genera_otutable <- t(periobaseline_physeq_norm_genera_otutable)
#Export
write.table(periobaseline_physeq_norm_genera_otutable, file="periobaseline_physeq_norm_genera_abundances.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

###Export normalised OTU count with Genus and Species names. Use periobaseline_mgseq_ren.

periobaseline_physeq_norm_species_otutable <- MRcounts(periobaseline_mgseq_ren, norm=TRUE)
#Round to entire sequences
periobaseline_physeq_norm_species_otutable <- round(periobaseline_physeq_norm_species_otutable, 0)
#Transpose
periobaseline_physeq_norm_species_otutable <- t(periobaseline_physeq_norm_species_otutable)
#Export
write.table(periobaseline_physeq_norm_species_otutable, file="periobaseline_physeq_norm_species_abundances.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")


####Extra
# plotMRheatmap of top 100 OTUs with highest variance
pdf(file="periobaseline_mrheat_100OTUs_health.pdf")
par(oma=c(5,4,4,5)+0.1)
trials = pData(periobaseline_mgseq_ren)$abs_status
heatmapColColors = brewer.pal(12,"Set3")[as.integer(factor(trials))]
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
plotMRheatmap(obj = periobaseline_mgseq_ren, n = 100, norm = TRUE, log = TRUE, fun = sd, cexRow = 0.4, cexCol = 0.3, trace = "none", col = heatmapCols, ColSideColors = heatmapColColors)
cores <- c("#8DD3C7", "#FB8072", "#BEBADA", "#FFFFB3")
legend(x=0.005, y=0.90, legend=c("Control","Shallow","Intermediate","Deep"), fill=cores, border=TRUE, bty="n", y.intersp = 0.7, cex=0.35)
dev.off()


###Tree
pdf(file="periobaseline_dissimilarity_clustering.pdf")
BCdist=distance(periobaseline_physeq_norm,"bray")
hclustavg<-hclust(BCdist,method= "average")
plot(hclustavg, hang=-1, main="AverageLinkage", cex=0.35)
#rect.hclust(hclustavg,k =3)

Jacdist=distance(periobaseline_physeq_norm,"jaccard")
hclustavg<-hclust(Jacdist,method= "average")
plot(hclustavg, hang=-1, main="AverageLinkage", cex=0.35)
#rect.hclust(hclustavg,k =3)
dev.off()
