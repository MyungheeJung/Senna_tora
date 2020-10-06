library(pheatmap); library(gplots); library(ggplot2); library(colorspace); library(wesanderson)

tpm = delim("Senna_tora_CHS-L.txt", header=T, row.names='gene')
bk = unique(c(seq(0, 0, length=10), 0, seq(0, 13, length=100)))
colors = colorRampPalette(c("white", "red"))(100)
pheatmap(tpm, 
         color = colors,
         breaks = bk,
         cluster_cols=F,
         cluster_rows=F,
         fontsize_col=9,
         fontsize_row=7,
         cellwidth = 260/ncol(tpm), 
         cellheight = 200/nrow(tpm), 
         border_color = "grey",
         show_rownames=T,
         gaps_col = c(1,2,3,4,5,6),
         labels_row = rownames(tpm),
         main=paste(c(nrow(tpm), " CHS-L genes"),collapse=""))