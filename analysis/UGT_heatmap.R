library(heatmaply)
library(htmlwidgets)
library(RColorBrewer)

cpmCounts <- read.table("results/salmon-marchRun.normalized_cpm.lowExpFiltered.txt")

# Shorten contig name for easy display
rownames(cpmCounts) <-
  gsub("allCormLibraries_Trinity_21Dec_TRINITY_", "", rownames(cpmCounts))

# Manual curated genes list
ugts <- scan("data/UGT_list_from_transcriptome.txt", what = "string")
# Read 159 items
ugts <- gsub("allCormLibraries_Trinity_21Dec_TRINITY_", "", ugts)

ugts <- subset(cpmCounts, rownames(cpmCounts) %in% ugts)
str(ugts)
# 'data.frame':	147 obs. of  4 variables:

# Create content for mouse-over bubble
cellnote <- matrix(data = rep(rownames(ugts), 4), 
                   nrow = length(rownames(ugts)), 
                   ncol = length(colnames(ugts)))
row.names(cellnote) <- rownames(ugts)
colnames(cellnote) <- colnames(ugts)

# Calculate distance matrix
rDist <- as.dist(1 - cor(t(ugts), method = "pearson"))
rHclust <- hclust(rDist, "complete")
rDendrogram <- as.dendrogram(rHclust)


heatmaply(
  main = "UDP-Glucuronosyltransferase Heatmap",
  file = "results/ugt_heatmap.14Dec.html",
  as.matrix(ugts[, c(1:4)]),
  colors = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
  col_text_angle = 90,
  draw_cellnote = FALSE,
  Rowv = rDendrogram,
  Colv = NULL,
  scale = "row",
  srtCol = 0,
  label_names = c("contig", "library", "CPM"),
  showticklabels = c(T, F),
  custom_hovertext = cellnote
)
