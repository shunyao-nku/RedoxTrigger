# Establish the correaltion beteween groups of samples to generate input file for Cytoscape
setwd("D://Zhenlai/Venn/Results")

## read data, ASV feature table and annotation information
otu <- read.delim('0.001_group_ASV.txt', row.names = 1, stringsAsFactors = FALSE)
tax <- read.delim('taxonomy.txt', row.names = 1, stringsAsFactors = FALSE)

##Convert to edge list (OTU versus samples), non-zero values indicate this ASV exists in the sample
edge <- otu
edge$OTU <- rownames(edge)
edge <- reshape2::melt(edge, id = 'OTU')

# delete edges that do not exist
edge <- subset(edge, value != 0)

# revise the name list
names(edge) <- c('source', 'target', 'weight')

# Add a column “shared name” showing conenctions of source and target
edge$'shared name' <- paste(edge$source, edge$target, sep = '->')

# output edge list used for subsequnet constructing Cytoscape Venn network
write.table(edge, 'edge.txt', sep = '\t', quote = FALSE, row.names = FALSE)

# output node property list
otu[otu>0] <- 1
otu$venn <- apply(otu, 1, function(x) paste(x, collapse = '-'))
otu <- otu[unique(edge$source), ]

# OTU property list
node_otu <- data.frame(
  'shared name' = rownames(otu),  
  group1 = tax[unique(edge$source),'Phylum'], 
  group2 = 'otu', 
  group3 = otu$venn,  
  stringsAsFactors = FALSE, 
  check.names = FALSE
)

# sample property list
otu <- otu[-ncol(otu)]

node_sample <- data.frame(
  'shared name' = names(otu),
  group1 = names(otu),  # to assign group colour
  group2 = names(otu),  # to assign group shape
  group3 = names(otu),  # to group
  stringsAsFactors = FALSE, 
  check.names = FALSE
)


