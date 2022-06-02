# top(n) abundance phylum or genus stack plot

library(ggalluvial)

# input feature table
otu<-read.table("D://Zhenlai_rewrite/meta/Result2/ASV_resampled_total.txt",sep="\t", header=T, row.names=1)

# input group information
design<-read.csv("D://Zhenlai_rewrite/meta/Result2/test_design.csv")

# input taxonomy annotation information
tax<-read.table("D://Zhenlai_rewrite/meta/Result2/filtered_taxonomy_sandcore_water_oil.txt",sep="\t", header=T, row.names=1)

# combine otu and taxonomy table
otu_tax<-cbind(otu,tax)

# Count abundance by phylum name
phyla<-aggregate(otu_tax[,1:171],by=list(otu_tax$Phylum),FUN=sum)

rownames(phyla)<-phyla[,1]

mean<-phyla[,-1]

a<-as.numeric(length(rownames(mean)))

# Screen Top13 abundance phyla and assigned other phyla as "others"
mean<-mean[order(rowSums(mean),decreasing=T),]

top_phyla<-rbind(colSums(mean[14:a,]),mean[13:1,])

rownames(top_phyla)[1]<-"Others"

bb<-cbind(t(top_phyla),design)


# Calculate the average abundance
taxonomy<-aggregate(bb[,1:14],by=list(bb$Treatment),FUN=mean)

rownames(taxonomy)<-taxonomy[,1]

tax<-as.data.frame(taxonomy)

library(ggprism)

library(ggplot2)

library(reshape2)

dat <- melt(tax, id = 'Group.1')

dat$Group.1 <-as.factor(dat$Group.1)


color<-c("gray34","#997D52","#4169B2","#B1A4C0","#479E9B","#BB2BA0","#DDA0DD",
         
         "#BC8F8F","#DD5F60","#F6906C","#FFDAB9","#A6CEE3","#B4EEB4","#8FBC8F")


# Draw stack graph
ggplot(data=dat,aes(x =dat$Group.1, y = value, fill =variable)) +
         
  geom_bar(position = "fill",stat = "identity",width = 0.6)+
  
  ylab("Relative abundance")+
  
  scale_fill_manual(values=c(color))+
  
  scale_y_reverse(expand = c(0,0),labels  = c("1","0.75","0.50","0.25","0"))+
  
  theme(axis.text.x = element_text(color = "black",size = 8,angle = 90))+
  
  theme(axis.text.y = element_text(color = "black",size = 8))+
  
  theme(legend.position = "right",legend.text = element_text(size = 7),
        
        panel.grid =element_blank())+scale_x_discrete(name='Habitat' # x-axis name
                                                      
        )+
  
  guides(fill=guide_legend(title="Genus",color="black",reverse=TRUE))+theme_bw()+
  
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  
  theme(axis.text=element_text(colour='black',size=9))

