
OTU <- read.table("D://Zhenlai_rewrite/Correlation_test/OTU_Phylum.txt", sep="\t", header=T, row.names=1)

OTU = t(OTU)

OTU <- scale(OTU)

OTU[1:18,17:18]


cor_spearman <- cor(OTU[1:143,1:16],OTU[1:143,17:18], method = 'spearman')
cor_spearman


install.packages("corrplot") 
library(corrplot)


col1=colorRampPalette(colors =c("darkgreen","white","darkred"),space="Lab")

corrplot(cor_spearman, method = 'square', addCoef.col = 'black', col = col1(21), number.cex = 0.7, tl.cex = 0.8,
         tl.col="black", tl.srt = 45,tl.offset=0.7)




cor_spearman <- cor(OTU, method = 'spearman')

cor_spearman[cor_spearman >= 1] <- 0

col1=colorRampPalette(colors =c("darkgreen","white","darkred"),space="Lab")


res1 <- cor.mtest(cor_spearman, conf.level = .95)


corrplot(cor_spearman, method = 'number', number.cex = 0.7, diag = FALSE, 
         tl.pos="lt", tl.cex=0.8, tl.col="black",tl.srt = 45,tl.offset=0.7,
         pch.cex = 1.2, pch.col="black", col = col1(10), p.mat=res1$p, order = "AOE")

corrplot(cor_spearman, add = TRUE, type = 'upper', method = 'circle', diag = TRUE, tl.pos = 'n', tl.cex=1, tl.col="black", 
         cl.pos = 'n',col = col1(10), p.mat = res1$p, insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .8, pch.col = "black", order = "AOE")

 #输出，例如
write.table(cor_pearson, 'cor_pearson.txt', sep = '\t', col.names = NA, quote = FALSE)