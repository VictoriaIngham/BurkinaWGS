library(tidyverse)
library(ggtree)
library(treeio)

library(stringr)


tree = read.nhx('/Users/vickyingham/Seafile/WGS/Analysis of vcf/RaXML of 3L/Galaxy27-[Best-scoring_ML_Tree].nhx')
tree@phylo$tip.label = str_extract(tree@phylo$tip.label, "[^_]+")

ggtree(tree) + 
  theme_tree2() + 
  geom_nodepoint(aes(alpha=0.4,color = "#009E73"),size=2) +
  geom_tippoint(aes(color = "#56B4E9"),alpha=0.9,pch=18,size=3) +
  theme(legend.position="none") +
  geom_tiplab(size=3)

ggtree(tree) + 
  theme_tree2() + 
  geom_tiplab(size=2)

