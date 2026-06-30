# Load packages ##
library(readr)
library(dplyr)
library(ggplot2)
library(limma)
library(edgeR)
library(tidyverse)
library(eulerr)
############################################################################################
## Load in the data ##
############################################################################################

## Load in the data ##

## Load in the data ##
path = "/Users/aciernia/Sync/collaborations/Jessica_bolton/Bolton\ Collaboration/2025-10-12-PVNSeq"
setwd(path)
counts <- read.delim("R747_raw_counts.tsv", header = TRUE, row.names = 1)
head(counts)
tail(counts)

#remove tail (last 5 lines)
counts <- head(counts, -5)

tail(counts)




Samples <- colnames(counts)[2:ncol(counts)]
Samples<-as.data.frame(Samples)
head(Samples)

#matrix gene names and ensembl IDs (column bind)
genenames = cbind(counts$gene_name,rownames(counts))
genenames = as.data.frame(genenames)
colnames(genenames) = c("genesymbol","EnsemblID")
head(genenames)


############################################################################################
## add gene length
############################################################################################

# BiocManager::install(c("GenomicFeatures","TxDb.Mmusculus.UCSC.mm10.knownGene","org.Mm.eg.db"))
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(GenomicRanges)
library(dplyr)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

ex_by_gene <- exonsBy(txdb, by = "gene")                 # names: Entrez IDs
gene_len_bp <- sapply(reduce(ex_by_gene), function(gr) sum(width(gr)))

len_df <- tibble(ENTREZID = names(gene_len_bp),
                 length_bp = as.integer(gene_len_bp))

# Map Entrez -> Ensembl + Symbol
map_df <- AnnotationDbi::select(org.Mm.eg.db,
                                keys = len_df$ENTREZID,
                                keytype = "ENTREZID",
                                columns = c("ENSEMBL","SYMBOL")) %>%
  distinct(ENTREZID, .keep_all = TRUE)

lengths_df <- len_df %>%
  left_join(map_df, by = "ENTREZID") %>%
  dplyr::select(ENSEMBL, SYMBOL, ENTREZID, length_bp)

lengths_df$length_kb = lengths_df$length_bp/1000

symbol_length <- lengths_df %>% dplyr::select(ENSEMBL,length_kb ) %>% distinct()

#add to genenames
genenames <- merge(genenames, symbol_length, by.x="EnsemblID",by.y="ENSEMBL")

#reorder to match orginal
idx <- match(rownames(counts), genenames$EnsemblID)

genenames_aligned <- genenames[idx, , drop = FALSE]


############################################################################################
## metadata
############################################################################################

# derive sample information # derive sample infocountsrmation 
metadata <- read.csv("Metadata.csv")
head(metadata)


#group interaction, set levels
metadata$Group = factor(metadata$Group, levels = c("CTL_F","CTL_M","ELA_F","ELA_M"))
Group <- metadata$Group
Group
metadata$Maternal_Care = factor(metadata$Maternal_Care, levels = c("CTL","ELA"))
Maternal_care <- metadata$Maternal_Care
Maternal_care 
metadata$Sex = factor(metadata$Sex, levels = c("F","M"))
Sex <- metadata$Sex
Sex
metadata$Isolation = factor(metadata$Isolation, levels = c("Normal","Yellow"))
Isolation <- metadata$Isolation
Isolation

#make counts matrix
counts.matrix <- counts[2:ncol(counts)]

#match metadata
metadata$Sample <- gsub("-","\\.",metadata$Sample)
counts.matrix<- counts.matrix[,metadata$Sample]
colnames(counts.matrix)

#make matrix
counts.matrix <- as.matrix(counts.matrix)
head(counts.matrix)




#Filtering lowly expressed genes #
d0 <- DGEList(counts.matrix, group = metadata$Group)
d0$genes = genenames_aligned
head(d0)
dim(d0)
#55573    36

############################################################################################
##getting rid of low expressed genes (cutoff for the number of samples that we want above 1,getting rid of things that are less than one group)
############################################################################################
cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) 
#15069    36

colnames(d) <- gsub("\\.", "-", colnames(d))

############################################################################################
#Making a filtering plot
############################################################################################
library(RColorBrewer)
nsamples <- ncol(d)
colourCount = nsamples
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
fill=getPalette(colourCount)
#plot:
pdf('FilteringCPM_plots_cutoff1-norm.pdf')
par(mfrow=c(1,2))
#prefilter:
lcpm <- cpm(d0, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2,
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)

for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}
#log-CPM of zero threshold (equivalent to a CPM value of 1) used in the filtering ste
lcpm <- cpm(d, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2, 
     main="", xlab="")

title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}
#legend("topright", Samples, text.col=fill, bty="n")
dev.off()


############################################################################################
##Plot library size##
############################################################################################
pdf('LibrarySizes.pdf',w=30,h=8)
barplot(d$samples$lib.size,names=colnames(d),las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()

#Plot Log counts#
# Log2 counts per million from filtered data (unnormalised)
logcounts <- cpm(d, log = TRUE)
library(tidyverse)

# Log2 counts per million from filtered data (unnormalised)
logcounts <- cpm(d, log = TRUE)

# Convert to long format for ggplot
logcounts_long <- as.data.frame(logcounts) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "LogCPM")


metadata$Sample <- gsub("\\.","-",metadata$Sample )

# Join metadata by matching Sample IDs
logcounts_long <- merge(logcounts_long, metadata, by="Sample")

# Check if join worked
head(logcounts_long)

# Median logCPM across all samples
global_median <- median(logcounts_long$LogCPM)

# ---- Standard ggplot (all samples together) ----
pdf("LogCPM-d_ggplot.pdf", width = 14, height = 8)
ggplot(logcounts_long, aes(x = Sample, y = LogCPM, fill = Group)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_hline(yintercept = global_median, color = "blue", linetype = "dashed") +
  labs(title = "Boxplots of logCPMs (unnormalised, filtered)",
       x = "Samples", y = "Log2 counts per million")
dev.off()

# ---- Facetted version by Group ----
pdf("LogCPM-d_byGroup.pdf", width = 14, height = 8)
ggplot(logcounts_long, aes(x = Sample, y = LogCPM, fill = Group)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~Group, scales = "free_x") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_hline(yintercept = global_median, color = "blue", linetype = "dashed") +
  labs(title = "Boxplots of logCPMs (unnormalised, filtered)", 
       x = "Samples", y = "Log2 counts per million")
dev.off()


# Save facetted boxplot by Isolation
pdf("LogCPM-d_byIsolation.pdf", width = 14, height = 8)
ggplot(logcounts_long, aes(x = Sample, y = LogCPM, fill = Isolation)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~Isolation, scales = "free_x") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_hline(yintercept = global_median, color = "blue", linetype = "dashed") +
  labs(title = "Boxplots of logCPMs (unnormalised, filtered) by Isolation",
       x = "Samples", y = "Log2 counts per million")
dev.off()

# Order samples by Percent_Unique
logcounts_long$Sample <- factor(logcounts_long$Sample, 
                                levels = metadata$Sample[order(metadata$Percent_Unique)])
global_median <- median(logcounts_long$LogCPM)


############################################################################################
####Design model matrix####
############################################################################################
mm <- model.matrix(~0 + Group + Isolation)
##Run Calculation Normalization
DGE=calcNormFactors(d,method =c("TMM")) 

pdf('Voom_with_isolation.pdf',w=6,h=4)
v=voom(DGE,mm,plot=T)
dev.off()

# Get log2 counts per million normalized
logcounts2 <- cpm(DGE,log=TRUE)
# Check distributions of samples using boxplots
pdf('NormalizedLogCPM1.pdf',w=30,h=10)
boxplot(logcounts2, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (normalised)")
dev.off()


############################################################################################
##MDS Plots
############################################################################################

#Make MDS Plots#
pdf('MDSplots.pdf')
par(mfrow=c(2,2))
#Multidimensional scaling (MDS) plot * colour by Group
col.cell <- c("pink","cyan","lavender","skyblue")[Group]
data.frame(Sex,col.cell)
plotMDS(DGE, col=col.cell,dim.plot = c(1,2), cex=0.4)
legend("topleft",
       fill= c("pink","cyan","lavender","skyblue"),
       legend=levels(Group),
       cex = 0.5)
title("Group")

#Multidimensional scaling (MDS) plot * colour by sex
col.cell <- c("red","blue")[Sex]
data.frame(Sex,col.cell)
# MDS with group colouring
plotMDS(DGE, col=col.cell,dim.plot = c(1,2), cex=0.4)
legend("topleft",
       fill= c("red","blue"),
       legend=levels(Sex),
       cex = 0.5)
title("Sex")

#Multidimensional scaling (MDS) plot * colour by Isolation
col.cell <- c("grey","yellow")[Isolation]
data.frame(Isolation,col.cell)
# MDS with group colouring
plotMDS(DGE, col=col.cell,dim.plot = c(1,2), cex=0.4)
legend("topleft",
       fill= c("grey", "yellow"),
       legend=levels(Isolation),
       cex = 0.5)
title("Isolation")

#Multidimensional scaling (MDS) plot * colour by maternal care
col.cell <- c("purple", "orange")[Maternal_care]
data.frame(Maternal_care,col.cell)
# MDS with group colouring
plotMDS(DGE, col=col.cell, dim.plot = c(1,2), cex=0.4)
legend("topleft",
       fill= c("purple", "orange"),
       legend=levels(Maternal_care),
       cex = 0.5)
title("Maternal Care")

dev.off()


#Make MDS Plots for dim3/4#
pdf('MDSplots_Dim3,4.pdf')
par(mfrow=c(2,2))
#Multidimensional scaling (MDS) plot * colour by Group
col.cell <- c("pink","cyan","lavender","skyblue")[Group]
data.frame(Sex,col.cell)
plotMDS(DGE, col=col.cell,dim.plot = c(3,4), cex=0.4)
legend("topleft",
       fill= c("pink","cyan","lavender","skyblue"),
       legend=levels(Group),
       cex = 0.5)
title("Group")

#Multidimensional scaling (MDS) plot * colour by sex
col.cell <- c("red","blue")[Sex]
data.frame(Sex,col.cell)
# MDS with group colouring
plotMDS(DGE, col=col.cell,dim.plot = c(3,4), cex=0.4)
legend("topleft",
       fill= c("red","blue"),
       legend=levels(Sex),
       cex = 0.5)
title("Sex")

#Multidimensional scaling (MDS) plot * colour by Isolation
col.cell <- c("grey","yellow")[Isolation]
data.frame(Isolation,col.cell)
# MDS with group colouring
plotMDS(DGE, col=col.cell,dim.plot = c(3,4), cex=0.4)
legend("topleft",
       fill= c("grey", "yellow"),
       legend=levels(Isolation),
       cex = 0.5)
title("Isolation")

#Multidimensional scaling (MDS) plot * colour by maternal care
col.cell <- c("purple", "orange")[Maternal_care]
data.frame(Maternal_care,col.cell)
# MDS with group colouring
plotMDS(DGE, col=col.cell, dim.plot = c(3,4), cex=0.4)
legend("topleft",
       fill= c("purple", "orange"),
       legend=levels(Maternal_care),
       cex = 0.5)
title("Maternal Care")

dev.off()


############################################################################################
### make contrast matrix
############################################################################################
mm2 <- model.matrix(~0 + Group + Isolation, data = metadata)
# colnames(mm2)

colnames(mm2) 

contr.matrix <- makeContrasts(
  ###Compare within sex CTL to ELA###
  CTLvsELA_M = GroupELA_M - GroupCTL_M,
  CTLvsELA_F = GroupELA_F - GroupCTL_F,
  
  ###Comapre within group sex effects##
  CTL_MvsF = GroupCTL_M - GroupCTL_F,
  ELA_MvsF = GroupELA_M - GroupELA_F,
  
  ### Combined effects over sexes for ELA and CNT
  CTLvsELA_both = (GroupELA_M + GroupELA_F)/2 - (GroupCTL_M + GroupCTL_F)/2,
  
  ### Isolation effect ###
  IsolationYellow_vs_Normal = IsolationYellow,
  
  levels = mm2
  )

contr.matrix

##Run Calculation Normalization

pdf('Voom.pdf',w=6,h=4)
v=voom(DGE,mm2,plot=T)
dev.off()

fit <- lmFit(v, mm2)


tmp <- contrasts.fit(fit, contr.matrix)
tmp <- eBayes(tmp)
summary(decideTests(tmp))

pdf('PlotSA_VoomTrend.pdf',w=6,h=4)
plotSA(tmp, main="Final model: Mean variance trend")
dev.off()

dt <- decideTests(tmp)
summary(dt)
write.csv(summary(dt), file = "DEGcounts.csv")

for (contrast in colnames(contr.matrix)) {
  top <- topTable(tmp, coef = contrast, number = Inf, adjust.method = "BH")
  write.csv(top, paste0("DEGs_", contrast, ".csv"), row.names = TRUE)
}


############################################################################################
###GLIMMA interactive plot building###
############################################################################################
#BiocManager::install("Glimma")
#http://bioconductor.org/packages/release/bioc/vignettes/Glimma/inst/doc/Glimma.pdf
library(Glimma)
#writes out html file: ** create MDS glimma from DGE_filtered***
glMDSPlot(DGE, groups=Group)

for (COEF in 1:6) {
  glMDPlot(tmp, counts=DGE,transform=FALSE,anno=v$genes,
           coef=COEF, status=dt, main=colnames(tmp)[COEF],
           groups=Group, folder="glimma_results", launch=FALSE, html = paste("MD-Plot",colnames(contr.matrix)[COEF]))}



# Create output folder if it doesn't exist
if (!dir.exists("glimma_results2")) {
  dir.create("glimma_results2")
}

# Loop through contrasts for interactive MD plots, grouped by Isolation
for (COEF in 1:ncol(contr.matrix)) {
  glMDPlot(tmp,
           counts = v$E,       # voom adjusted CPM values
           transform = FALSE,           # if TRUE, log2 counts are used
           anno = v$genes,              # gene annotation
           coef = COEF,
           status = dt,                 # significant DEGs from decideTests
           main = colnames(contr.matrix)[COEF],
           groups = Isolation,         # <-- group by Isolation
           folder = "glimma_results",
           launch = FALSE,
           html = paste0("MD-Plot_", colnames(contr.matrix)[COEF], "_Isolation.html"))
}





############################################################################################
#####Volcano plots for each comparison
############################################################################################

library(calibrate)
library(dplyr)

####Make contrasts####
comparisons=(coef(tmp))
comparisons=colnames(comparisons)
comp_out <- as.data.frame(rownames(v$E))
names(comp_out) <- c("GeneID")
nrowkeep <- nrow(comp_out)

SumTableOut <- NULL

for(i in 1:length(comparisons)){
  #comparison name
  comp=comparisons[i]
  print(comp)
  #make comparisons 
  
  topT=topTreat(tmp,coef=i,number=nrowkeep,adjust.method="BH")
  print(nrow(topT[(topT$adj.P.Val<0.05),]))
  
  #LogFC values:https://support.bioconductor.org/p/82478/
  topT$direction <- c("none")
  topT$direction[which(topT$logFC > 0)] = c("Increase")
  topT$direction[which(topT$logFC < 0)] = c("Decrease")
  
  topT$significance <- c("nonDE")
  topT$significance[which(topT$adj.P.Val <0.05)] <- c("DE")
  
  #summary counts table based on Ensemble Gene ID counts:
  SumTable <- table(topT$significance,topT$direction)
  SumTable <- as.data.frame(SumTable)
  SumTable$comparison <- paste(comp)
  SumTableOut <- rbind(SumTable,SumTableOut)
  
  
  
  #gene gene names and expression levels
  topT2 <- topT
  topT2$comparison <- paste(comp)
  write.csv(topT2,file = paste(comp,"_DEgenes.csv"))
  
  #get master file:
  colnames(topT)[3:ncol(topT)] <- paste(colnames(topT)[3:ncol(topT)],comp)
  comp_out <- merge(comp_out,topT, by.x = "GeneID" , by.y= "EnsemblID")
  
  #data for plot with gene names:
  genenames <- topT2 %>% dplyr::select(adj.P.Val,logFC,genesymbol) %>% distinct()
  
  #names for plots
  plotname <- gsub("\\."," ",comp)
  plotname <- gsub("vs"," vs ",plotname)
  
  #volcano plot
  pdf(file = paste(comp,"_Volcano.pdf", sep=""), wi = 9, he = 6, useDingbats=F)
  
  with(genenames, plot(logFC, -log10(adj.P.Val), pch=20,col="gray", main=paste(plotname,"\nVolcano plot", sep=" "), ylab =c("-log10(adj.pvalue)"),xlab =c("Log Fold Change") ))
  
  #color points red when sig and log2 FC > 1 and blue if log2 FC < -1 
  with(subset(genenames, logFC < -1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  with(subset(genenames, logFC > 1 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  
  #add lines
  abline(h = -log10(.05), col = c("black"), lty = 2, lwd = 1)
  abline(v = c(-1,1), col = "black", lty = 2, lwd = 1)
  
  #Label points with the textxy function from the calibrate plot
  library(calibrate)
  with(subset(genenames, adj.P.Val<0.05 & abs(logFC)>1), textxy(logFC, -log10(adj.P.Val), labs=genesymbol, cex=.5))
  
  dev.off()
  
}

write.csv(SumTableOut,"SummaryTableDEgenes.csv")


############################################################################################
#####Average CPM for each condition#####
############################################################################################
###get log2CPM counts from voom and put in dataframe:
library(plotrix)
#average log2 CPM and sem
countdf <- as.data.frame(v$E)
countdf$Ensembl <- rownames(v$E)
#add gene names
genenames$EnsemblID <- rownames(genenames)
genenames <- genenames %>% dplyr::select(-adj.P.Val,-logFC)
DF <- merge(countdf,genenames, by.x ="Ensembl",by.y="EnsemblID")
#write as csv
write.csv(DF,file="log2CPMvalues.csv")
head(DF)

#summarize 
countdf2 <- DF %>%  group_by(Ensembl,genesymbol) %>% 
  tidyr::gather(sample,log2CPM, 2:37)
countdf2 <- as.data.frame(countdf2)

#add in metadata:
metadata$Sample <- gsub("\\.","-",metadata$Sample)
countdf3 <-merge(countdf2,metadata, by.x = "sample", by.y = "Sample")
head(countdf3)

#Write CSV for Gene Summary
GeneSummary <- countdf3 %>% group_by(Ensembl,genesymbol,Group) %>% 
  summarize(meanlog2CPM = mean(log2CPM)) %>%
  ungroup()  %>%
  tidyr::spread(Group ,meanlog2CPM)


GeneSummary <- as.data.frame(GeneSummary)
head(GeneSummary)
write.csv(GeneSummary, file = "AverageLog2CPM.csv")

library(dplyr)

library(dplyr)

# 1) Repair duplicate names in comp_out
comp_out2 <- as.data.frame(comp_out)
names(comp_out2) <- make.unique(names(comp_out2), sep = "__dup")

# 2) Keep only what you need from GeneSummary
gene_lookup <- GeneSummary %>%
  distinct(Ensembl, .keep_all = TRUE) %>%
  dplyr::select(Ensembl, genesymbol)

# 3) Join with clear suffixes
mout <- comp_out2 %>%
  left_join(gene_lookup, by = c("GeneID" = "Ensembl"), suffix = c(".comp", ".sum"))

# 4) Collapse to a single genesymbol column if comp_out already had one
if ("genesymbol.comp" %in% names(mout)) {
  mout <- mout %>%
    mutate(genesymbol = coalesce(genesymbol.comp, genesymbol)) %>%
    dplyr::select(-ends_with(".comp"), -ends_with(".sum"))
}


write.csv(mout,"AllDEG_AllConditions_Avelog2CPM.csv")


#rpkm
rpkm_mat   <- edgeR::rpkm(DGE, gene.length = v$genes$length_kb, log = FALSE, normalized.lib.sizes = TRUE)
rpkm_mat   <- as.data.frame(rpkm_mat)
rpkm_mat$Ensembl <- rownames(rpkm_mat)

rpkm_wide <- rpkm_mat %>% group_by(Ensembl) %>% gather(sample, RPKM,1:36)

#add genesymboles
DF <- merge(rpkm_wide, countdf3, by=c("sample","Ensembl"))

write.csv(DF,"RPKM_Data.csv")






############################################################################################
#####heatmap for BOTH degs #####
############################################################################################

# 3) wide matrix: one column per sample, rows = genes
rpkm_wide <- rpkm_long %>%
  pivot_wider(
    names_from  = sample,
    values_from = RPKM,
    values_fill = NA_real_        # or 0 if you prefer zeros
  )

#read in modified full list of DEGS
mout2 <- read.csv("AllDEG_AllConditions_Avelog2CPM.csv")

#extract DEGS for BOTH
ELAvsCTL_both <- mout %>%
  dplyr::select(`adj.P.Val CTLvsELA_both`, genesymbol, GeneID) %>%
  filter(`adj.P.Val CTLvsELA_both` < 0.05)

#filter matrix
rpkm_mat_filt <- rpkm_mat[rownames(rpkm_mat) %in% ELAvsCTL_both$genesymbol, , drop = FALSE]


library(dplyr)
library(pheatmap)

## ---- Column annotations (aligned to rpkm_mat_filt) ----
sample_meta <- DF %>%
  dplyr::select(sample, Maternal_Care, Sex, Isolation) %>%
  dplyr::distinct(sample, .keep_all = TRUE) %>%
  dplyr::rename(Sample = sample) %>%
  dplyr::filter(Sample %in% colnames(rpkm_mat_filt)) %>%
  dplyr::arrange(match(Sample, colnames(rpkm_mat_filt)))

col_annot <- sample_meta %>%
  dplyr::select(Maternal_Care, Sex, Isolation) %>%
  as.data.frame()
rownames(col_annot) <- sample_meta$Sample

## ---- Colors for annotations (as given) ----
# maternal care
var1 <- c("blue","orange")
names(var1) <- unique(col_annot$Maternal_Care)

# isolation
var2 <- c("grey","yellow")
names(var2) <- unique(col_annot$Isolation)

# sex
var3 <- c("purple","blue")
names(var3) <- unique(col_annot$Sex)

anno_colors <- list(
  `Maternal_Care` = var1,
  `Sex`           = var3,
  Isolation       = var2
)

## ---- Optional: row-scale (z-score) for visibility ----
log2_mat_scaled <- log2(rpkm_mat_filt+1)
mat_scaled <- t(scale(t(log2_mat_scaled)))  # NAs for constant rows are fine

## ---- Heatmap ----
pheatmap(
  mat_scaled,
  show_rownames   = FALSE,
  show_colnames   = TRUE,
  cluster_rows    = TRUE,
  cluster_cols    = T,
  fontsize        = 9,
  annotation_col  = col_annot,
  annotation_colors = anno_colors,
  main = "CTRLvsELA BOTH DEGs (log2RPKM row-zscore)"
)

## ---- Optional: save to PDF ----
pdf("CTRLvsELA BOTH DEGs (log2RPKM row-zscore).pdf", width = 9, height = 11)
pheatmap(
  mat_scaled,
  show_rownames   = FALSE,
  show_colnames   = TRUE,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  fontsize        = 9,
  annotation_col  = col_annot,
  annotation_colors = anno_colors,
  main = "CTRLvsELA BOTH DEGs (log2RPKM row-zscore)"
)
dev.off()

#order columns

# --- (1) Ensure annotation rows line up with matrix columns ---
stopifnot(all(colnames(rpkm_mat_filt) %in% rownames(col_annot)))
col_annot <- col_annot[colnames(rpkm_mat_filt), , drop = FALSE]

# --- (2) Optional: set explicit factor orders (edit to your desired order) ---
col_annot$Maternal_Care <- factor(col_annot$Maternal_Care, levels = c("CTL","ELA"))
col_annot$Sex          <- factor(col_annot$Sex,          levels = c("F","M"))
col_annot$Isolation    <- factor(col_annot$Isolation,    levels = c("Normal","Yellow"))

# --- (3) Order columns by Maternal_Care, then Sex, then Isolation (tie-breaker = sample name) ---
ord <- do.call(order, list(
  col_annot$Maternal_Care,
  col_annot$Sex,
  col_annot$Isolation,
  colnames(rpkm_mat_filt)
))

rpkm_mat_ord <- rpkm_mat_filt[, ord, drop = FALSE]
col_annot_ord <- col_annot[ord, , drop = FALSE]

# --- (4) Use log2 first, then row-wise z-score ---
rpkm_log2   <- log2(rpkm_mat_ord + 1)
mat_scaled  <- t(scale(t(rpkm_log2)))   # NAs for constant rows are fine

# --- (5) Heatmap: NO column clustering, keep your sorted order ---
pheatmap(
  mat_scaled,
  show_rownames    = F,
  show_colnames    = TRUE,
  cluster_rows     = TRUE,     # set FALSE if you also want to keep row order
  cluster_cols     = FALSE,    # <-- no column clustering
  fontsize         = 9,
  annotation_col   = col_annot_ord,
  annotation_colors = anno_colors,   # uses your existing color list
  main = "log2(RPKM+1), row-scaled — columns sorted by Maternal_Care → Sex → Isolation"
)

# Optional: save to PDF
pdf("CTRLvsELA BOTH DEGs (log2RPKM row-zscore)_sorted.pdf", width = 9, height = 11)
pheatmap(
  mat_scaled,
  show_rownames    = F,
  show_colnames    = TRUE,
  cluster_rows     = TRUE,     # set FALSE if you also want to keep row order
  cluster_cols     = FALSE,    # <-- no column clustering
  fontsize         = 9,
  annotation_col   = col_annot_ord,
  annotation_colors = anno_colors,   # uses your existing color list
  main = "log2(RPKM+1), row-scaled — columns sorted by Maternal_Care → Sex → Isolation"
)
dev.off()




save.image("Image.RData")
