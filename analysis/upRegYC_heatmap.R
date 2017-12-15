library(dplyr)
library(heatmaply)
library(htmlwidgets)
library(RColorBrewer)
library(tibble)

### Read DEA results
sigDEA <- read.table("results/OldvsYoungCorms.sigDE.marchRun.txt", header = TRUE)

# Find ID of those that are up-regulated in YC
ycUp <- subset(sigDEA, sigDEA$logFC <= -3)
dim(ycUp)
# [1] 2135    6

# Shorten contig name for easy display
rownames(ycUp) <-
  gsub("allCormLibraries_Trinity_21Dec_TRINITY_", "", rownames(ycUp))


### Read CPM counts
cpmCounts <- read.table("results/salmon-marchRun.normalized_cpm.lowExpFiltered.txt")

# Shorten contig name for easy display
rownames(cpmCounts) <-
  gsub("allCormLibraries_Trinity_21Dec_TRINITY_", "", rownames(cpmCounts))

# Subset CPM for dataset of interest
sigDE_cpmCounts <- cpmCounts[rownames(ycUp), ]
dim(sigDE_cpmCounts)
# [1] 2135    4


### Read BLAST annotations
annots <- read.delim("data/OldvsYoungCorms.sigDE.marchRun.blastxRefSeqPlant85.txt", 
                    header = TRUE, stringsAsFactors = FALSE)
str(annots)
# 'data.frame':	17954 obs. of  16 variables:

# flattenedAnnots <- annots %>% select (qseqid, salltitles) %>% group_by(qseqid) %>% 
#   summarise(annots = paste(salltitles, collapse=" ")) %>% 
#   as.data.frame() %>% column_to_rownames(var = "qseqid")

# Only retain the top hit from BLAST result
flattenedAnnots <- annots %>% select (qseqid, salltitles) %>% group_by(qseqid) %>%
  summarise(annots = salltitles[1]) %>%
  as.data.frame() %>% column_to_rownames(var = "qseqid")
str(flattenedAnnots)
# 'data.frame':	3563 obs. of  1 variable:

# Shorten contig name for easy display
rownames(flattenedAnnots) <-
  gsub("allCormLibraries_Trinity_21Dec_TRINITY_", "", rownames(flattenedAnnots))


# Subset CPM for dataset of interest
sigDE_cpmCounts_annots <- 
  merge(sigDE_cpmCounts, flattenedAnnots, by.x = 0, by.y = 0, all.x = TRUE)
sigDE_cpmCounts_annots <- column_to_rownames(sigDE_cpmCounts_annots, "Row.names")
str(sigDE_cpmCounts_annots)
# 'data.frame':	2135 obs. of  5 variables:

cellnote <- data.frame(row.names = row.names(sigDE_cpmCounts_annots), 
                       "OC1" = sigDE_cpmCounts_annots[,5],
                       "OC3" = sigDE_cpmCounts_annots[,5],
                       "YC1" = sigDE_cpmCounts_annots[,5],
                       "YC3" = sigDE_cpmCounts_annots[,5]
                       )

rDist <- as.dist(1 - cor(t(sigDE_cpmCounts_annots[,c(1:4)]), method = "pearson"))
rHclust <- hclust(rDist, "complete")
rDendrogram <- as.dendrogram(rHclust)


heatmaply(
  main = "Heatmap for DE contigs up-regulated in Old Corms with logFC >= 3",
  file = "results/upRegYC_heatmap.14Dec.html",
  sigDE_cpmCounts_annots[, c(1:4)],
  colors = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
  col_text_angle = 90,
  Rowv = rDendrogram,
  Colv = NULL,
  scale = "row",
  label_names = c("contig", "library", "CPM"),
  srtCol = 0,
  showticklabels = c(T, F),
  custom_hovertext = as.matrix(cellnote)
)
