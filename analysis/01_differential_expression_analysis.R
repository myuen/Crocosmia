library(dplyr)
library(edgeR)
library(stringr)

# There are in total of 4 libraries from 2 corm tissues, young 
# and old, with 2 biological replicates each.

# Read Salmopn output file names
salmonFiles <- list.files("./data/salmon-march-results/", full.names = TRUE)

readFile = function(x){
  content <- read.table(x, header = TRUE, colClasses = c("character", rep("numeric", 4)))
  libname <- str_extract(x, "[A-Z]{2}[0-9]{1}")
  # content <- content[,c("Name", "TPM")]
  content <- content[,c("Name", "NumReads")]
  colnames(content)[2] <- libname
  return(content)}

# Read each output file, extract appropriate columns and merge into single data frame
salmonRes <- lapply(salmonFiles, readFile)
salmonResDf <- bind_cols(salmonRes)
row.names(salmonResDf) <- salmonResDf$Name
salmonResDf <- salmonResDf[, c("OC1", "OC3", "YC1", "YC3")]
str(salmonResDf)
# 'data.frame':	56147 obs. of  4 variables:


## Differential Expression Analysis ###
# p-value cutoff
pCutoff <- 0.08

# log fold change cutoff
lfcCutoff <- 2


## Create library descriptions to group libraries
libDesc <- data.frame("libNames" = colnames(salmonResDf))
libDesc$libNames <- as.character(libDesc$libNames)
libDesc$age <- 
  as.factor(c(rep("old", 2), rep("young", 2)))
libDesc$age <- relevel(libDesc$age, ref = "young")
#   libNames   age
# 1      OC1   old
# 2      OC3   old
# 3      YC1 young
# 4      YC3 young


y <- DGEList(counts = salmonResDf, group = libDesc$age)
dim(y)
# [1] 56147     4

write.table(cpm(y), "results/salmon-marchRun.non-normalized_cpm.noLowExpFilter.txt",
            row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# Keep only genes with at least 1 count-per-million reads (cpm) in at least 2 samples
y <- y[(rowSums(cpm(y) >= 1) >= 2), ]
dim(y) 
# [1] 36471     4

write.table(cpm(y), "results/salmon-marchRun.non-normalized_cpm.lowExpFilter.txt",
            row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# Reset depth
y$samples$lib.size <- colSums(y$counts)

# Normalize by Depth
y <- calcNormFactors(y)

write.table(cpm(y), "results/salmon-marchRun.normalized_cpm.lowExpFiltered.txt",
            row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)


design <- model.matrix(~ libDesc$age)
colnames(design) <- gsub("libDesc\\$age", "", colnames(design))


# Voom normalization
v <- voom(y, design, plot = TRUE)

fit <- lmFit(v, design)

fit2 <- eBayes(fit)

summary(decideTests(fit2, adjust.method = "fdr", p.value = pCutoff, lfc = lfcCutoff))
#    (Intercept)   old
# -1         122  3208
# 0        15193 32596
# 1        21156   667

allResult <- topTable(fit2, number = Inf, adjust.method = "fdr", sort.by = "none")

write.table(allResult, "results/OldvsYoungCorms.allResults.marchRun.txt", 
            row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

write.table(
  subset(allResult, allResult$adj.P.Val <= pCutoff & abs(allResult$logFC) >= lfcCutoff), 
  "results/OldvsYoungCorms.sigDE.marchRun.txt", 
  row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
