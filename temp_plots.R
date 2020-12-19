library(gplots)

deg_dir <- './DEG_RNA-seq/'
image_dir <- './images/'

counts_norm <- read.csv('./coverage_RNA-seq/counts.TPM-normalized.csv', row.names = 1, header = TRUE)
counts_norm[counts_norm == 0] <- NA
counts_norm_rel <- cbind(counts_norm$WT_YES.AMM_1 / counts_norm$WT_YES_1,
                         counts_norm$WT_YES.AMM_2 / counts_norm$WT_YES_2,
                         counts_norm$WT_YES.AMM_3 / counts_norm$WT_YES_3,
                         counts_norm$cbf11_YES_1 / counts_norm$WT_YES_1,
                         counts_norm$cbf11_YES_2 / counts_norm$WT_YES_2,
                         counts_norm$cbf11_YES_3 / counts_norm$WT_YES_3,
                         counts_norm$cbf11_YES.AMM_1 / counts_norm$cbf11_YES_1,
                         counts_norm$cbf11_YES.AMM_2 / counts_norm$cbf11_YES_2,
                         counts_norm$cbf11_YES.AMM_3 / counts_norm$cbf11_YES_3,
                         counts_norm$cbf11_YES.AMM_1 / counts_norm$WT_YES.AMM_1,
                         counts_norm$cbf11_YES.AMM_2 / counts_norm$WT_YES.AMM_2,
                         counts_norm$cbf11_YES.AMM_3 / counts_norm$WT_YES.AMM_3)
counts_norm_rel <- counts_norm_rel[complete.cases(counts_norm_rel), ]

pdf(paste0(image_dir, 'counts_norm_rel_heatmap.pdf'), height = 7, width = 7)
heatmap.2(log2(as.matrix(counts_norm_rel)),
          Colv = FALSE,
          dendrogram = 'row', 
          trace = 'none',
          scale = 'none', 
          density.info = 'none',
          keysize = 1,
          col = colorRampPalette(c("blue", 'black', 'yellow'))(20),
          cex.axis = 0.5,
          margins = c(8, 5),
          breaks = seq(-2, 2, 0.2))
dev.off()


wt <- read.csv(paste0(deg_dir, 'DESeq2results_WT-YES+AMM_vs_WT-YES.csv'), header = TRUE, row.names = 1)
cbf11_wt <- read.csv(paste0(deg_dir, 'DESeq2results_cbf11-YES_vs_WT-YES.csv'), header = TRUE, row.names = 1)
cbf11 <- read.csv(paste0(deg_dir, 'DESeq2results_cbf11-YES+AMM_vs_cbf11-YES.csv'), header = TRUE, row.names = 1)
cbf11_wt_a <- read.csv(paste0(deg_dir, 'DESeq2results_cbf11-YES+AMM_vs_WT-YES+AMM.csv'), header = TRUE, row.names = 1)

fce <- cbind(wt[, c('gene_name', 'log2FoldChange')], cbf11_wt$log2FoldChange, 
             cbf11$log2FoldChange, cbf11_wt_a$log2FoldChange)
fce <- fce[complete.cases(fce), ]
fce <- fce[which(rowSums(abs(fce[, -1])) >= 2), ]
colnames(fce) <- c('gene_name', 'wt', 'cbf11_wt', 'cbf11', 'cbf11_wt_a')

pdf(paste0(image_dir, 'log2fce_heatmap.pdf'), height = 7, width = 7)
heatmap.2(as.matrix(fce[, -1]),
          Colv = FALSE,
          dendrogram = 'row', 
          trace = 'none',
          scale = 'none', 
          density.info = 'none',
          keysize = 1,
          col = colorRampPalette(c("blue", 'black', 'yellow'))(10),
          margins = c(8, 5),
          breaks = seq(-5, 5, 1), 
          labRow = '',
          labCol = c('WT+A / WT', 'cbf11 / WT', 'cbf11+A / cbf11', 'cbf11+A / WT+A'),
          cexCol = 1)
dev.off()


sig_WTA_WT <- read.csv(paste0(deg_dir, 'DESeq2results_WT-YES+AMM_vs_WT-YES.SIG.csv'), 
                       row.names = 1, stringsAsFactors = FALSE, header = TRUE)
sig_WTA_WT <- sig_WTA_WT[which(rownames(sig_WTA_WT) %in% rownames(sig_11_WT)), ]
sig_WTA_WT <- sig_WTA_WT[order(rownames(sig_WTA_WT)), ]

sig_11_WT <- read.csv(paste0(deg_dir, 'DESeq2results_cbf11-YES_vs_WT-YES.SIG.csv'), 
                       row.names = 1, stringsAsFactors = FALSE, header = TRUE)
sig_11_WT <- sig_11_WT[which(rownames(sig_11_WT) %in% rownames(sig_WTA_WT)), ]
sig_11_WT <- sig_11_WT[order(rownames(sig_11_WT)), ]

pdf(paste0(image_dir, 'log2fce_scatter_sig.pdf'), height = 7, width = 7)
plot(sig_WTA_WT$log2FoldChange, sig_11_WT$log2FoldChange,
     xlim = range(c(sig_11_WT$log2FoldChange, sig_WTA_WT$log2FoldChange)),
     ylim = range(c(sig_11_WT$log2FoldChange, sig_WTA_WT$log2FoldChange)),
     col = rgb(0, 0, 0, 0.3),
     pch = 20,
     xlab = 'WT+A / WT (log2FCE)',
     ylab = 'cbf11 / WT (log2FCE)')
abline(v = 0, h = 0, a = 0, b = 1)
abline(lm(sig_11_WT$log2FoldChange ~ sig_WTA_WT$log2FoldChange), col = 'blue')
points(sig_WTA_WT[which(sig_WTA_WT$gene_name == 'lsd90'), 'log2FoldChange'], 
       sig_11_WT[which(sig_11_WT$gene_name == 'lsd90'), 'log2FoldChange'], 
       pch = 20, col = 'red')
dev.off()