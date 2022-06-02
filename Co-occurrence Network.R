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

# Acquire correaltion coefficients r and p values
r_matrix=occor$r
p_matrix=occor$p

r_matrix[p_matrix>p.cutoff|abs(r_matrix)<r.cutoff]=0
p_value=melt(p_matrix)
r_value=melt(r_matrix)
r_value=cbind(r_value, p_value$value)
r_value=subset(r_value, r_value[,3]!=0)
r_value=na.omit(r_value)
abs=abs(r_value$value)

# Add positive or negative correaltion information
linktype=r_value$value
linktype[linktype>0]=1
linktype[linktype<0]=-1
r_value=cbind(r_value, abs, linktype)
names(r_value) <- c("Source","Target","r_value","p_value", "abs_value", "linktype")
names(p_value) <- c("Source","Target","p_value")

# Output edge data used for co-occurrence network construction in Gephi or Cytoscape
write.csv(r_value,file="D://Zhenlai_rewrite/Subnetwork/PW8910/Edges.csv", row.names=FALSE)

node_OTU <- as.data.frame(as.data.frame(r_value[,1])[!duplicated(as.data.frame(r_value[,1])), ])
node_Ev <- as.data.frame(as.data.frame(r_value[,2])[!duplicated(as.data.frame(r_value[,2])), ])
names(node_OTU)="Id"
names(node_Ev)="Id"
list <- rbind(node_Ev, node_OTU)
list <- rbind(node_OTU)
list=subset(tax,Id %in% list$Id)
list$Label <- list$Id

# Output node data used for co-occurrence network construction in Gephi or Cytoscape
write.csv(list,file="D://Zhenlai_rewrite/Subnetwork/PW8910/Nodes.csv", row.names=FALSE)



###igraph
igraph<-graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
igraph
# remove isolated nodes
bad.vs<-V(igraph)[degree(igraph) == 0]
igraph<-delete.vertices(igraph, bad.vs)
igraph




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


