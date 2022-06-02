#top10丰度phylum或class的冲击组成图 (输入对应格式数据后，代码无需修改，就能出图)

otu<-read.table("D://个人资料/博一/镇赉岩心样品/Analysis/feature-table.txt",sep="\t", header=T, row.names=1)

design<-read.csv("D://个人资料/博一/镇赉岩心样品/Analysis/test_design.csv")

tax<-read.csv("D://test100/test_tax.csv",row.names = 1)

tax<-read.table("D://个人资料/博一/镇赉岩心样品/Analysis/taxonomy.txt",sep="\t", header=T, row.names=1)

otu_tax<-cbind(otu,tax)

#链接：https://pan.baidu.com/s/1KpV3jhlcEf8DakZFKjsxmg 提取码：qxxd

########### 对top 10 phylum基因进行统计

phyla<-aggregate(otu_tax[,1:24],by=list(otu_tax$Genus),FUN=sum)

rownames(phyla)<-phyla[,1]

mean<-phyla[,-1]

a<-as.numeric(length(rownames(mean)))

#筛选top 10丰度物种

mean<-mean[order(rowSums(mean),decreasing=T),]

top_phyla<-rbind(colSums(mean[14:a,]),mean[13:1,])

rownames(top_phyla)[1]<-"Others"

bb<-cbind(t(top_phyla),design)

#求处理直接丰度的平均值

taxonomy<-aggregate(bb[,1:14],by=list(bb$Treatment),FUN=mean)

rownames(taxonomy)<-taxonomy[,1]

tax<-as.data.frame(taxonomy)

library(ggprism)

library(ggplot2)

library(reshape2)

dat <- melt(tax, id = 'Group.1')

dat$Group.1 <-as.factor(dat$Group.1)



color<-c("gray34","#4169B2","#B1A4C0","#479E9B","#BB2BA0","#DDA0DD",
         
         "#BC8F8F","#DD5F60","#FFDAB9","#B4EEB4","#8FBC8F")

color<-c("#0055AA","#7FD2FF","#007ED3","#C40003","#00C19B","#EAC862","#B2DF8A","#FFACAA","#FF9D1E","#C3EF00","#CAB2D6","#894FC6","skyblue","red")

library(ggalluvial)

ggplot(data=dat,aes(x =dat$Group.1, y = value, fill =variable)) +
         
  geom_bar(position = "fill",stat = "identity",width = 0.6)+
  
  ylab("Relative abundance")+
  
  scale_fill_manual(values=c(color))+
  
  scale_y_reverse(expand = c(0,0),labels  = c("1","0.75","0.50","0.25","0"))+
  
  theme(axis.text.x = element_text(color = "black",size = 8,angle = 90))+
  
  theme(axis.text.y = element_text(color = "black",size = 8))+
  
  theme(legend.position = "right",legend.text = element_text(size = 7),
        
        panel.grid =element_blank())+scale_x_discrete(name='Well Number' #x轴坐标名称
                                                      
        )+
  
  guides(fill=guide_legend(title="Genus",color="black",reverse=TRUE))+theme_bw()+
  
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  
  theme(axis.text=element_text(colour='black',size=9))


