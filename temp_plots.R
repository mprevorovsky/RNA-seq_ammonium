library(gplots)

counts <- read.csv('./coverage_RNA-seq/counts.TPM-normalized.csv', row.names = 1, header = TRUE)

pdf('counts_norm.pdf', height = 7, width = 7)
heatmap.2(as.matrix(counts[, -1]),
          dendrogram = 'row', 
          trace = 'none',
          scale = 'row', 
          density.info = 'none',
          keysize = 1,
          col = colorRampPalette(c("blue", 'black', 'yellow'))(20),
          cex.axis = 0.5,
          margins = c(8, 5),
          breaks = seq(-2, 2, 0.2))
dev.off()

deg_dir <- './DEG_RNA-seq/'

wt <- read.csv(paste0(deg_dir, 'DESeq2results_WT-YES+AMM_vs_WT-YES.csv'), header = TRUE, row.names = 1)
cbf11_wt <- read.csv(paste0(deg_dir, 'DESeq2results_cbf11-YES_vs_WT-YES.csv'), header = TRUE, row.names = 1)
cbf11 <- read.csv(paste0(deg_dir, 'DESeq2results_cbf11-YES+AMM_vs_cbf11-YES.csv'), header = TRUE, row.names = 1)
cbf11_wt_a <- read.csv(paste0(deg_dir, 'DESeq2results_cbf11-YES+AMM_vs_WT-YES+AMM.csv'), header = TRUE, row.names = 1)

fce <- cbind(wt[, c('gene_name', 'log2FoldChange')], cbf11_wt$log2FoldChange, 
             cbf11$log2FoldChange, cbf11_wt_a$log2FoldChange)
fce <- fce[complete.cases(fce), ]
fce <- fce[which(rowSums(abs(fce[, -1])) >= 2), ]
colnames(fce) <- c('gene_name', 'wt', 'cbf11_wt', 'cbf11', 'cbf11_wt_a')

pdf('log2fce.pdf', height = 7, width = 7)
par(cex = 0.5)
heatmap.2(as.matrix(fce[, -1]),
          Colv = FALSE,
          dendrogram = 'row', 
          trace = 'none',
          scale = 'none', 
          density.info = 'none',
          keysize = 1,
          col = colorRampPalette(c("blue", 'black', 'yellow'))(10),
          margins = c(8, 5),
          breaks = seq(-5, 5, 1))
dev.off()
