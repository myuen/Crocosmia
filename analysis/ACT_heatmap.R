library(heatmaply)
library(htmlwidgets)
library(RColorBrewer)

cpmCounts <- read.table("results/salmon-marchRun.normalized_cpm.lowExpFiltered.txt")

# Shorten contig name for easy display
rownames(cpmCounts) <-
  gsub("allCormLibraries_Trinity_21Dec_TRINITY_", "", rownames(cpmCounts))

# Manual curated genes list
acts <- scan("data/putative_ACT.contigId.txt", what = "string")
# Read 61 items
acts <- gsub("allCormLibraries_Trinity_21Dec_TRINITY_", "", acts)

acts <- subset(cpmCounts, rownames(cpmCounts) %in% acts)
str(acts)
# 'data.frame':	51 obs. of  4 variables:

# Create content for mouse-over bubble
cellnote <- matrix(data = rep(rownames(acts), 4), 
                   nrow = length(rownames(acts)), 
                   ncol = length(colnames(acts)))
row.names(cellnote) <- rownames(acts)
colnames(cellnote) <- colnames(acts)

# Calculate distance matrix
rDist <- as.dist(1 - cor(t(acts), method = "pearson"))
rHclust <- hclust(rDist, "complete")
rDendrogram <- as.dendrogram(rHclust)


heatmaply(
  main = "Acyltransferase Heatmap",
  file = "results/act_heatmap.14Dec.html",
  as.matrix(acts[, c(1:4)]),
  colors = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
  col_text_angle = 90,
  draw_cellnote = FALSE,
  Rowv = rDendrogram,
  Colv = NULL,
  scale = "row",
  label_names = c("contig", "library", "CPM"),
  srtCol = 0,
  custom_hovertext = cellnote
)
