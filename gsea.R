#######################################################################
# FINAL, WORKING GSEA PIPELINE
#######################################################################

library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)

set.seed(123)

setwd("/Users/aciernia/Sync/collaborations/Jessica_bolton/Bolton\ Collaboration/2025-10-12-PVNSeq")
#######################################################################
# 1. LOAD MSIGDB GENE SETS (SYMBOL-BASED)
#######################################################################

msig <- msigdbr(species = "Mus musculus")
msigdbr_collections()


HallmarkCurated <- msigdbr(species = "mouse", collection = "H")
HallmarkCurated <- HallmarkCurated %>% dplyr::select(gs_name, gene_symbol)

GO_sets <- msigdbr(species = "mouse", collection = "C5")
GO_sets <- GO_sets %>% filter(gs_collection_name != "Human Phenotype Ontology")
GO_sets <- GO_sets %>% dplyr::select(gs_name, gene_symbol)

Immune_sets <- msigdbr(species = "mouse", collection = "C7")
Immune_sets <-Immune_sets %>% filter(gs_collection_name != "HIPC Vaccine Response")
Immune_sets <-Immune_sets %>% dplyr::select(gs_name, gene_symbol)

Pathways <- msigdbr(species = "mouse", collection = "C2")
Pathways <- Pathways %>% filter(gs_collection_name != "Chemical and Genetic Perturbations")
Pathways <- Pathways %>% dplyr::select(gs_name, gene_symbol)

#######################################################################
# 2. MAP ENSEMBL → SYMBOL (WITHOUT SHRINKING GENE UNIVERSE)
#######################################################################

mout <- read.csv("AllDEG_AllConditions_Avelog2CPM.csv")


#######################################################################
# 3. PREPARE RANKED GENE LISTS
#######################################################################

geneLists <- mout %>% dplyr::select(contains("logFC"))
geneLists$SYMBOL <- mout$genesymbol.x

# Restrict to biologically meaningful contrasts
conditionnames <- grep(
  "CTLvsELA",
  colnames(geneLists),
  value = TRUE
)


####################################################
#run enrichments for each list x DB
####################################################
#loop through each LPS condition vs PBS log2FC and run GSEA enrichment for different MsigDBs
motifoutput <- NULL
enrich_Immuneoutput <- NULL
enrich_HCoutput <- NULL
enrich_GOout <- NULL

masterout <- NULL
for (i in conditionnames) {
  
  test1 <- geneLists[,i]
  names(test1) <- geneLists$SYMBOL
  
  # Sort
  test1 <- sort(test1, decreasing = TRUE)
  
  #remove duplicates
  test1 <- test1[!duplicated(names(test1))]
  
  # Sort
  test1 <- sort(test1, decreasing = TRUE)
  
  
  #pathway enrichment:
  enrich_motif <- GSEA(
    geneList  = test1,
    TERM2GENE = Pathways,
    pvalueCutoff = 1,
    verbose = FALSE)
  
  enrich_motif <- as.data.frame(enrich_motif )
  #if no enrichments, then skip:
  if(nrow(enrich_motif)>0) {
    enrich_motif$comparision <- as.character(colnames(geneLists[i]))
    enrich_motif$msigDB <- c("Pathways c2")
    motifoutput <- rbind(motifoutput,enrich_motif)}
  
  #Immune enrichment:
  enrich_Immune  <- GSEA(
    geneList  = test1,
    TERM2GENE = Immune_sets,
    pvalueCutoff = 1,
    verbose = FALSE)
  
  enrich_Immune <- as.data.frame(enrich_Immune )
  #if no enrichments, then skip:
  if(nrow(enrich_Immune)>0) {
    enrich_Immune$comparision <- as.character(colnames(geneLists[i]))
    enrich_Immune$msigDB <- c("immune gene set C7")
    enrich_Immuneoutput <- rbind(enrich_Immuneoutput,enrich_Immune)}
  
  #Immune enrichment:
  enrich_hallmark <- GSEA(
    geneList  = test1,
    TERM2GENE = HallmarkCurated,
    pvalueCutoff = 1,
    verbose = FALSE)
  
  enrich_HC <- as.data.frame(enrich_hallmark)
  #if no enrichments, then skip:
  if(nrow(enrich_HC)>0) {
    enrich_HC$comparision <- as.character(colnames(geneLists[i]))
    enrich_HC$msigDB <- c("Hallmark and Currated gene sets")
    enrich_HCoutput <- rbind(enrich_HCoutput,enrich_HC)}
  
  #GO enrichment:
  enrich_GO  <- GSEA(
    geneList  = test1,
    TERM2GENE = GO_sets,
    pvalueCutoff = 1,
    verbose = FALSE)
  
  enrich_GO <- as.data.frame( enrich_GO)
  #if no enrichments, then skip:
  if(nrow( enrich_GO)>0) {
    enrich_GO$comparision <- as.character(colnames(geneLists[i]))
    enrich_GO$msigDB <- c("GO terms")
    enrich_GOout <- rbind( enrich_GOout, enrich_GO)}
  
}

#add back descriptions
HallmarkCurated <- msigdbr(species = "mouse", collection = "H")
HallmarkCurated <- HallmarkCurated %>% dplyr::select(gs_name,gs_description) %>% unique()

GO_sets <- msigdbr(species = "mouse", collection = "C5")
GO_sets <- GO_sets %>% filter(gs_collection_name != "Human Phenotype Ontology")
GO_sets <- GO_sets %>% dplyr::select(gs_name, gs_description) %>% unique()

Immune_sets <- msigdbr(species = "mouse", collection = "C7")
Immune_sets <-Immune_sets %>% filter(gs_collection_name != "HIPC Vaccine Response")
Immune_sets <-Immune_sets %>% dplyr::select(gs_name, gs_description) %>% unique()

Pathways <- msigdbr(species = "mouse", collection = "C2")
Pathways <- Pathways %>% filter(gs_collection_name != "Chemical and Genetic Perturbations")
Pathways <- Pathways %>% dplyr::select(gs_name,gs_description) %>% unique()

pathway_output <- merge(Pathways, motifoutput, by.x="gs_name", by.y= "ID")
Immune_output <- merge(Immune_sets,enrich_Immuneoutput, by.x="gs_name", by.y= "ID")
GO_output <- merge(GO_sets,enrich_GOout, by.x="gs_name", by.y= "ID")
Hallmark_output <- merge(HallmarkCurated,enrich_HCoutput, by.x="gs_name", by.y= "ID")


#write out full files
write.csv(pathway_output,"Allpathway_GSEAenrichments.csv")
write.csv(Immune_output,"AllImmune_GSEAenrichments.csv")
write.csv(GO_output,"AllHallmark_GSEAenrichments.csv")
write.csv(Hallmark_output,"AllGO_GSEAenrichments.csv")

#combine
GSEA_all <- rbind(pathway_output,Immune_output, GO_output, Hallmark_output )

write.csv(GSEA_all,"AllGSEA_GSEAenrichments.csv")

#######################################################################

#NES = Normalized Enrichment Score
#In GSEA (Gene Set Enrichment Analysis), the NES is the primary effect‑size statistic. It tells you how strongly a pathway (gene set) is enriched, while making results comparable across pathways and analyses.
#ES (Enrichment Score): Walks down the ranked gene list and measures whether genes from a pathway cluster at the top or bottom.
#Normalization: The ES is divided by the mean ES from many permutations (or fgsea’s adaptive null) for that same pathway.
#NES → coordinated gene set behavior
#######################################################################

#clean and filter
library(dplyr)

gsea_clean <- GSEA_all %>%
  filter(
    !is.na(NES),
    !is.na(p.adjust),
    p.adjust < 0.05        # drop totally uninformative results
  ) %>%
  mutate(
    neglog10_padj = -log10(p.adjust),
    direction = ifelse(NES > 0, "Positive NES", "Negative NES")
  )

#Pick top 10 enrichments per comparison

top_gsea <- gsea_clean %>%
  group_by(comparision) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  ungroup()



# Second plot option: PBS diet effects
p <- top_gsea %>%
  ggplot(aes(
    x = comparision,
    y = gs_description,
    size = neglog10_padj,
    fill = NES
  )) +
  
  geom_point(
    shape = 21,        # filled circles only
    color = "black",   # black outline
    alpha = 0.85
  ) +
  
  scale_size_continuous(
    name = expression(-log[10]~adj~p),
    range = c(2, 9)
  ) +
  
  scale_fill_gradient2(
    name = "NES",
    low = "#2166AC",     # blue = negative NES
    mid = "white",
    high = "#B2182B",    # red = positive NES
    midpoint = 0
  ) +
  
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  ) +
  
  labs(
    title = "GSEA: CTRL vs ELA",
    x = "Comparison",
    y = "Gene set"
  )

ggsave(
  "GSEA_Dotplot_MSigDB.pdf",
  p,
  width = 20,
  height = 10
)


#pull genes for heatmap

############################################################
## GSEA-leading-edge heatmap 
## Includes row annotation for NES direction
############################################################

## Required libraries
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

############################################################
## 1. Load GSEA results (top NES pathways)
############################################################



############################################################
## 2. Extract leading‑edge genes from core_enrichment
############################################################

# Each core_enrichment entry is "GENE1/GENE2/GENE3/..."
genes_for_heatmap <- gsea_top$core_enrichment |>
  strsplit("/") |>
  unlist() |>
  unique()

message("Number of unique leading-edge genes: ",
        length(genes_for_heatmap))

############################################################
## 3. Subset expression matrix to these genes
############################################################


DF <- read.csv("RPKM_Data.csv")

expr <- DF %>% dplyr::select(genesymbol, Sex, Maternal_Care, sample, RPKM)

############################################################
## 3. Subset long-format expression DF
############################################################

# Keep only leading-edge genes and relevant samples (example: females)
DF_sub <- expr |>
  filter(genesymbol %in% genes_for_heatmap)

DF_sub <- DF_sub %>% arrange(Maternal_Care,Sex)

############################################################
## 4. Build gene × sample expression matrix
############################################################
library(tidyverse)

# Using  as expression values
expr_mat <- DF_sub |>
  dplyr::select(sample, genesymbol, RPKM) |>
  distinct() |>
  pivot_wider(
    names_from  = sample,
    values_from = RPKM
  ) |>
  as.data.frame()

# Set gene symbols as rownames
rownames(expr_mat) <- expr_mat$genesymbol
expr_mat$genesymbol <- NULL

# Convert to numeric matrix
expr_mat <- as.matrix(expr_mat)

############################################################
## 5. Row-wise Z-score (recommended)
############################################################
log2expr_mat <- log2(expr_mat+1)
expr_mat_z <- t(scale(t(log2expr_mat)))

############################################################
## 6. Column annotations from DF
############################################################
library(ComplexHeatmap)
# Build metadata per sample (unique rows only)
metadata <- DF |>
  dplyr::select(sample, Maternal_Care,Sex) |>
  distinct()

# Align metadata to matrix columns
metadata <- metadata[
  match(colnames(expr_mat_z), metadata$sample),
]

rownames(metadata) <- metadata$sample

ha_col <- HeatmapAnnotation(
  Maternal_Care = metadata$Maternal_Care,
  Sex = metadata$Sex,
  col = list(
    Maternal_Care = c(
      CTL    = "#56B4E9",
      ELA = "#E69F00"
    ),
    Sex = c(
      F = "#CC79A7",
      M = "#0072B2"
    )
  ),
  annotation_name_side = "left"
)

############################################################
## 7. Row annotation: NES direction
############################################################

gene_direction <- gsea_top |>
  dplyr::select(Direction, core_enrichment) |>
  separate_rows(core_enrichment, sep = "/") |>
  distinct(core_enrichment, .keep_all = TRUE) |>
  dplyr::rename(genesymbol = core_enrichment)

row_direction <- gene_direction$Direction[
  match(rownames(expr_mat_z), gene_direction$genesymbol)
]

ha_row <- rowAnnotation(
  `NES Direction` = row_direction,
  col = list(
    `NES Direction` = c(
      "Up in ELA"   = "#B2182B",
      "Down in ELA" = "#2166AC"
    )
  )
)

############################################################
## 8. Draw heatmap
############################################################
min(expr_mat_z)
max(expr_mat_z)

library(circlize)

h <- Heatmap(
  expr_mat_z,
  name = "Z-score",
  
  top_annotation  = ha_col,
  left_annotation = ha_row,
  
  show_row_names = TRUE,
  show_column_names = FALSE,
  
  cluster_rows = TRUE,
  cluster_columns = F,
  
  col = colorRamp2(
    c(-6, 0, 4),
    c("#2166AC", "white", "#B2182B")
  ),
  
  row_title = "Leading-edge genes\nTop GSEA pathways\nELA vs CTRL",
  column_title = "Samples"
)

pdf('Heatmpa_GSEA_ELAvCTRL_MSIGDB.pdf', w = 10, h = 20)
h
dev.off()





