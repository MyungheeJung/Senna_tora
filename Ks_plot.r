library(ggplot2); library(ggrepel); library(gridExtra)
args = commandArgs(trailingOnly=TRUE)

library(limma)
Cats = c()
KaKs = data.frame()
for (i in args){
  Data = data.frame(category = strsplit2(i, '\\.')[1,1], read.table(i, sep = '', header = T, quote = '', fill = T))
  KaKs = rbind(KaKs, Data)
  Label = strsplit2(i, '\\/')
  Cats = c(Cats, strsplit2(Label[1, ncol(Label)], '\\.')[1,1])
}

KaKs$category = factor(KaKs$category, levels = Cats)
NumOrg = length(Cats)


Ks = ggplot(KaKs) + 
  geom_density(aes(Ks, fill = factor(category)), alpha = 0.3, colour = NA) + 
  scale_fill_manual(values = rainbow(NumOrg)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 4)) + 
  ggtitle('Ks') +  
  theme_bw(base_size = 6) +
  xlab('Ks') +
  theme(legend.position = 'bottom',legend.title = element_blank(),plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 20))

pdf(paste(c('./Plot_Ks.', Cats, '.pdf'), collapse = ''))
grid.arrange(Ks)
dev.off()