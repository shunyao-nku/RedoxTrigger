
data(iris)
aa<-iris

data <- vegdist(aa[,1:4], method = "bray")

pcoa<- pcoa(data, correction = "none", rn = NULL)




library("ggplot2")

library("vegan")


phylum <- read.table("D://Zhenlai_rewrite/PCA/filtered_ASV-table_sandcore_water_PCoA.txt", sep="\t", header=T, row.names=1)
phylum<-t(phylum)

phylum<- decostand(phylum, method = 'hellinger')

bray_dis <- vegdist(phylum, method = 'bray')

#输出距离矩阵
write.table(as.matrix(bray_dis), 'D://Zhenlai_rewrite/PCA/bray_distance.txt', sep = '\t', col.names = NA, quote = FALSE)

pcoa <- cmdscale(bray_dis, k = (nrow(phylum) - 1), eig = TRUE)

#各 PCoA 轴的特征值
pcoa_eig <- pcoa$eig
#先评估下负特征值（末尾几个轴）
barplot(pcoa_eig)

pcoa_exp <- pcoa$eig/sum(pcoa$eig)

#样方排序坐标
site <- pcoa$point
#或
site <- scores(pcoa)
#输出
write.table(site, 'pcoa_site.txt', sep = '\t', col.names = NA, quote = FALSE)


group <- read.table("D://Zhenlai_rewrite/PCA/Group.txt", sep="\t", header=T, row.names=1)

pcoa_eig<-data.frame(scores(pcoa)[,1:2])

pcoa_eig$sample<-rownames(pcoa_eig)

pcoa_eig<-merge(pcoa_eig,group,by='sample')

pcoa_eig1<-round(100*pcoa$eig[1:2]/sum(pcoa$eig),2)

p <- ggplot(pcoa_eig, aes(Dim1, Dim2)) +
  
  geom_point(aes(color = group)) +
  
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
  
  scale_color_manual(values = c('red3','#B2DE8A','green3','#9CCCFA',"#8200C9","#FF7B00","#B05929"))+
  
  scale_fill_manual(values = c('red3','#B2DE8A','green3','#9CCCFA',"#8200C9","#FF7B00","#B05929")) +
  
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'))+
  
  labs(x = paste('PCoA1: ', pcoa_eig1[1], '%'), y = paste('PCoA2: ',  pcoa_eig1[2], '%')) +
  
  scale_shape_manual(values = c(17,16,16,16,17,17,19.18)) +
  
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  
  geom_hline(yintercept = 0, color = 'gray', size = 0.5)

  p
  
  ggsave('D://Zhenlai_rewrite/PCA/pcoa.pdf', p, width = 6, height = 4.7)
  