library(heatmaply)
library(htmlwidgets)
library(RColorBrewer)

cpmCounts <- read.table("results/salmon-marchRun.normalized_cpm.lowExpFiltered.txt")

# Shorten contig name for easy display
rownames(cpmCounts) <-
  gsub("allCormLibraries_Trinity_21Dec_TRINITY_", "", rownames(cpmCounts))

# Manual curated genes list for flavanoid pathway
manAnn <- read.delim("data/putative_flavanoid_pathway.contigId.txt", header = FALSE, 
                     row.names = 1, stringsAsFactors = FALSE)
colnames(manAnn) <- c("annot")
rownames(manAnn) <- gsub("allCormLibraries_Trinity_21Dec_TRINITY_", "", rownames(manAnn))

manAnn <- merge(manAnn, cpmCounts, by.x = 0, by.y = 0, sort = FALSE)
row.names(manAnn) <- manAnn$Row.names
str(manAnn)
# 'data.frame':	11 obs. of  6 variables:

# manAnnCNote <- matrix(rep(manAnn$annot, 4), nrow = 11, ncol = 4)
# row.names(manAnnCNote) <- row.names(manAnn)
# colnames(manAnnCNote) <- colnames(manAnn)[3:6]

heatmaply(
  main = "Flavonoid Pathway Expression Heatmap",
  # file = "results/flavonoidPwy_heatmap.14Dec.html",
  # file = "results/flavonoidPwy_heatmap.14Dec.pdf",
  as.matrix(manAnn[, c(3:6)]),
  draw_cellnote = FALSE,
  # cellnote = manAnnCNote,
  # cellnote_color = "white",
  # cellnote_textposition = "middle center",
  colors = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
  col_text_angle = 90,
  dendrogram = "none",
  scale = "row",
  # label_names = c("contig", "library", "CPM"),
  srtCol = 0
)
