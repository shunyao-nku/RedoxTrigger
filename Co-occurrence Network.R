# Co-occurrence network, topological properties and their corresponding random networks

library(igraph)
library(psych)
library(reshape2)

# read environmental factor table, feature table, taxonomy annotation table 
Ev <- read.table("D://Zhenlai/Network2/Ev.txt", sep="\t", header=T, row.names=1)
OTU <- read.table("D://Zhenlai/Subnetwork/PW8910/ASV_PW8910.txt", sep="\t", header=T, row.names=1)
tax <- read.table("D://Zhenlai/Subnetwork/PW567/filtered_taxonomy_sandcore_water_oil.txt", sep="\t", header=T)
names(tax)[1] <- "Id"

# If ASV features are too many, it is feasible to conduct an abundance filtration before calculating correlations
abundance=0.01
OTU <- OTU[,colSums(OTU)/sum(OTU)>=(abundance/100)]

# set confidential level
r.cutoff=0.7
p.cutoff=0.01

# Ev-OTU spearman's correlation calculation 
Ev=t(Ev)
OTU=t(OTU)
occor=corr.test(OTU, Ev,
                use="pairwise",
                method="spearman", # 可选pearson/kendall
                adjust="fdr",
                alpha=0.05)

# OTU-OTU spearman's correlation calculation 
OTU=t(OTU)
occor=corr.test(OTU,
     use="pairwise",
     method="spearman",
     adjust="fdr",
     alpha=0.05)

####数据量小时可以用psych包corr.test求相关性矩阵，数据量大时，可应用WGCNA中corAndPvalue,但p值需要借助其他函数矫正
occor<-corr.test(otu,use="pairwise",method="spearman",adjust="fdr",alpha=.05)
occor.r<-occor$r ###取相关性矩阵R值
occor.p<-occor$p ###取相关性矩阵p值
###确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
occor.r[occor.p>0.05|abs(occor.r)<0.6]<-0 

###构建igraph对象
igraph<-graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
igraph
###NOTE:可以设置weighted=NULL,但是此时要注意此函数只能识别相互作用矩阵内正整数，所以应用前请确保矩阵正确。
###可以按下面命令转换数据
#occor.r[occor.r!=0]<- 1
#igraph<-graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=NULL,diag=FALSE)

###是否去掉孤立顶点，根据自己实验而定
# remove isolated nodes，即去掉和所有otu均无相关性的otu 可省略，前期矩阵已处理过
bad.vs<-V(igraph)[degree(igraph) == 0]
igraph<-delete.vertices(igraph, bad.vs)
igraph

###将igraph weight属性赋值到igraph.weight
igraph.weight<-E(igraph)$weight

###做图前去掉igraph的weight权重，因为做图时某些layout会受到其影响
E(igraph)$weight<-NA
###简单出图
#设定随机种子数，后续出图都从同一随机种子数出发，保证前后出图形状相对应
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))


##如果构建网络时，weighted=NULL,此步骤不能统计
sum(igraph.weight>0)###number of postive correlation
sum(igraph.weight<0)###number of negative correlation

#set edge color，postive correlation 设定为red, negative correlation设定为blue
E.color<-igraph.weight
E.color[E.color>0]<-"red"
E.color[E.color<0]<-"blue"
E(igraph)$color <-as.character(E.color)


###改变edge颜色后出图
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))


##可以设定edge的宽度，例如将相关系数与edge width关联
#set edge width
E(igraph)$width<-abs(igraph.weight)*4

###改变edge宽度后出图
set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))


###另外可以设置vertices size, vertices color来表征更多维度的数据
####注意otu_pro.txt文件为我随机产生的数据，因此网络图可能不会产生特定的模式或规律。
otu_pro<-read.table("otu_pro.txt",head=T,row.names=1)
#set vertices size
igraph.size<-otu_pro[V(igraph)$name,]
igraph.size1<-log((igraph.size$abundance)*100)
V(igraph)$size<-igraph.size1

#set vertices color
igraph.col<-otu_pro[V(igraph)$name,]
levels(igraph.col$phylum)
levels(igraph.col$phylum)<-c("green","deeppink","deepskyblue","yellow","brown","pink","gray","cyan","peachpuff")
V(igraph)$color <- as.character(igraph.col$phylum)

set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))


####改变layout,layout有很多，具体查看igraph官方帮助文档。
set.seed(123)
plot(igraph,main="Co-occurrence network",layout=layout_with_kk,vertex.frame.color=NA,vertex.label=NA,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

set.seed(123)
plot(igraph,main="Co-occurrence network",layout=layout.fruchterman.reingold,vertex.frame.color=NA,vertex.label=NA,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))


###network property
##The size of the graph (number of edges)
num.edges<-length(E(igraph))##length(curve_multiple(igraph))
##Order (number of vertices) of a graph
num.vertices<-length(V(igraph))##length(diversity(igraph, weights = NULL, vids = V(igraph)))
###connectance
connectance<-edge_density(igraph,loops=FALSE)##同 graph.density;loops如果为TRUE,允许自身环（self loops即A--A或B--B）的存在
##average degree
average.degree<-mean(igraph::degree(igraph))##或者为2M/N,其中M 和N 分别表示网络的边数和节点数。
##average path length 平均路径长度
average.path.length<-average.path.length(igraph)  #同mean_distance(igraph) ##mean_distance calculates the average path length in a graph
##diameter
diameter<-diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
edge.connectivity<-edge_connectivity(igraph)
####clustering coefficient
clustering.coefficient<-transitivity(igraph)
no.clusters<-no.clusters(igraph)
centralization.betweenness<-centralization.betweenness(igraph)$centralization
centralization.degree<-centralization.degree(igraph)$centralization
###modularity
fc<-cluster_fast_greedy(igraph,weights =NULL)#####cluster_walktrap cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
modularity<-modularity(igraph,membership(fc))


####按照模块为节点配色
comps <- membership(fc)
colbar <- rainbow(max(comps))
V(igraph)$color <- colbar[comps] 

set.seed(123)
plot(igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))





#CK随机网络的构建
#构建与100个点359个边的随机网络，点和边的量，手动改成自己构建的网络图中点和边的量，方便比较
g <- erdos.renyi.game(832, 6591, "gnm") 
degree_distribution = degree_distribution(g)
degree_distribution
write.table(degree_distribution, 'degree_sandcore.txt', sep = '\t', col.names = NA, quote = FALSE)


#计算随机生成1000次指定网络节点数与网络边数随机网络的Average path length(APL)平均值和标准差
apl<-mean(replicate(1000, average.path.length(erdos.renyi.game(281, 923, "gnm"))))
apl

apl.sd<-sd(replicate(1000, average.path.length(erdos.renyi.game(281, 923, "gnm"))))
apl.sd

#计算随机生成10000次指定网络节点数与网络边数随机网络的聚集系数(Clustering coefficient)平均值和标准差
cc<-mean(replicate(1000, transitivity(erdos.renyi.game(281, 923, "gnm"))))
cc
cc.sd<-sd(replicate(1000, transitivity(erdos.renyi.game(281, 923, "gnm"))))
cc.sd

#计算随机生成1000次指定网络节点数与网络边数随机网络的modularity(M)平均值和标准差

f <- function() {
  g <- erdos.renyi.game(132, 163, "gnm")
  fc <- cluster_fast_greedy(g, weights=NULL)
  modularity(g, membership(fc))
}
mo<-mean(replicate(1000, f()))
mo
mo.sd<-sd(replicate(1000, f()))
mo.sd

#计算随机生成1000次指定网络节点数与网络边数随机网络的直径平均值和标准差
gd<-mean(replicate(1000, diameter(erdos.renyi.game(281, 923, "gnm"))))
gd
gd.sd<-sd(replicate(1000, diameter(erdos.renyi.game(281, 923, "gnm"))))
gd.sd

#计算随机生成1000次指定网络节点数与网络边数随机网络的平均度average degree平均值和标准差
ad<-mean(replicate(1000, igraph::degree(erdos.renyi.game(281, 923, "gnm"))))
ad
ad.sd<-sd(replicate(1000, igraph::degree(erdos.renyi.game(281, 923, "gnm"))))
ad.sd

#计算随机生成1000次指定网络节点数与网络边数随机网络的平均度graph density平均值和标准差
gd<-mean(replicate(1000, graph.density(erdos.renyi.game(281, 923, "gnm", loops=FALSE))))
gd


