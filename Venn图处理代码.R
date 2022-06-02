
setwd("D://Zhenlai_rewrite/Venn/Results")

##读取数据，OTU 丰度表和注释信息
otu <- read.delim('0.001_group_ASV.txt', row.names = 1, stringsAsFactors = FALSE)
tax <- read.delim('taxonomy.txt', row.names = 1, stringsAsFactors = FALSE)

##转换为边列表，即 OTU 和样本的对应关系
#结果中，所有非 0 的值代表该 OTU 在该样本中存在丰度
edge <- otu
edge$OTU <- rownames(edge)
edge <- reshape2::melt(edge, id = 'OTU')

#删除不存在的边
#即去除 0 丰度的值，该 OTU 不在该样本中存在
edge <- subset(edge, value != 0)

#修改列名称（以便后续软件识别），其中权重可表示为 OTU 在样本中的丰度
names(edge) <- c('source', 'target', 'weight')

#添加一列“shared name”，以“->”连接 source 和 target
#便于 cytoscape 读取识别，并防止读取后的名称错乱
edge$'shared name' <- paste(edge$source, edge$target, sep = '->')

#输出边列表，后续可导入至 cytoscape 用于构建网络
write.table(edge, 'edge.txt', sep = '\t', quote = FALSE, row.names = FALSE)

##获取节点属性列表
#获得 OTU 在各组中的分布状态（Venn 分布）
otu[otu>0] <- 1
otu$venn <- apply(otu, 1, function(x) paste(x, collapse = '-'))
otu <- otu[unique(edge$source), ]

#OTU 属性列表
node_otu <- data.frame(
  'shared name' = rownames(otu),  #OTU 名称，以“shared name”命名列名称，便于 cytoscape 读取识别，并防止读取后的名称错乱
  group1 = tax[unique(edge$source),'Phylum'],  #用于后续按指定分组赋值颜色，这里定义为 OTU 所属的门分类水平
  group2 = 'otu',  #用于后续按指定分组定义形状，这里统一分组为“otu”
  group3 = otu$venn,  #用于后续按指定分组调整聚群，按 OTU 在各组中的分布状态定义
  stringsAsFactors = FALSE, 
  check.names = FALSE
)

#样本属性列表
otu <- otu[-ncol(otu)]

node_sample <- data.frame(
  'shared name' = names(otu),  #样本名称，以“shared name”命名列名称，便于 cytoscape 读取识别，并防止读取后的名称错乱
  group1 = names(otu),  #用于后续按指定分组赋值颜色，指定样本分组
  group2 = names(otu),  #用于后续按指定分组定义形状，指定样本分组
  group3 = names(otu),  #用于后续按指定分组调整聚群，指定样本分组
  stringsAsFactors = FALSE, 
  check.names = FALSE
)

#二者合并构建节点属性列表后输出，后续可导入至 cytoscape 用于调整节点可视化
node <- rbind(node_otu, node_sample)
write.table(node, 'node.txt', sep = '\t', quote = FALSE, row.names = FALSE)
