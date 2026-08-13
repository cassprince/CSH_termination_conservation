# Load packages and set working directory.
library(tidyverse)
library(treeio)
library(ggtree)
library(phangorn)
library(ggnewscale)
library(castor)
library(ggbreak)
library(ggprism)

setwd("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//CSH")

# Define function to format tables for the tree.
table_for_tree = function(phylum, gene, accessions){
  table = as.data.frame(prop.table(table(phylum, gene), margin = 1)*100)
  table_new = data.frame(cbind(gene = table$Freq[table$gene == "TRUE"]))
  rownames(table_new) = accessions
  return(table_new)
}

# Upload dataframe containing NCBI Assembly and Nucleotide accession numbers, assembly statistics, taxonomy, and gene presence.(Table S1).
df_GCF_nc_compressed = read.csv("df_csh_081126.csv") 

# Separate into one Nucleotide accession number per row.
df_GCF_nc = separate_rows(df_GCF_nc_compressed, nuccore, sep = ";")

# Make a smaller dataframe with only data from phyla with >10 genomes.
dfShort = df_GCF_nc_compressed %>%
  group_by(phylum) %>%
  filter(n() > 10)

n_vals = dfShort %>%
  group_by(phylum) %>%
  summarize(n = n())

random = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Daniel//random_genomes_dt_tree_051826.csv") %>%
  mutate(X = assembly) %>%
  filter(phylum != "Cyanobacteriota" & phylum != "Bacteroidota" & phylum != "Planctomycetota" & phylum != "Verrucomicrobiota")

random = column_to_rownames(random, 'X')

####### ---- Collapsed tree ---- #######

# Upload tree.
tree = read.newick("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Daniel//phyla_051826.tre")

# Fix GCF accession labels and drop two leaves that are duplicated phyla.
tree$tip.label = sub("_([A-Z]|_).*", "", tree$tip.label)
tree_trim = drop.tip(tree, tip = c("GCF_041929205.1", "GCF_046000445.1", "GCF_038400315.1", "GCF_033878955.1"))

# Midpoint root.
tree_mid = midpoint(tree_trim)

# Format tables for gene presence to append to the tree.
table_prfA = table_for_tree(dfShort$phylum, dfShort$prfA, random$assembly)
table_prfB = table_for_tree(dfShort$phylum, dfShort$prfB, random$assembly)
table_prfC = table_for_tree(dfShort$phylum, dfShort$prfC, random$assembly)
table_prmC = table_for_tree(dfShort$phylum, dfShort$prmC, random$assembly)
table_prfH = table_for_tree(dfShort$phylum, dfShort$prfH, random$assembly)
table_selB = table_for_tree(dfShort$phylum, dfShort$selB, random$assembly)
table_pylS = table_for_tree(dfShort$phylum, dfShort$pylS, random$assembly)

genes = data.frame(cbind(RF1 = table_prfA$gene, RF2 = table_prfB$gene, RF3 = table_prfC$gene, PrmC = table_prmC$gene, PrfH = table_prfH$gene, SelB = table_selB$gene, PylS = table_pylS$gene))
rownames(genes) = rownames(random)


# Plot the tree and heatmap.
p = ggtree(tree_mid) %<+% random + 
  xlim(0, 10) + 
  geom_tiplab(aes(label=phylum), size = 4) + 
  geom_nodepoint(aes(fill = as.numeric(label)*100), size = 2, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage") + 
  new_scale_fill()

gheatmap(p, genes, offset = 5.5, width=2, font.size=3.5, colnames = TRUE, color=NA, colnames_angle = 45, colnames_offset_y = -0.2, colnames_offset_x = -0.15) + 
  scale_fill_viridis_c(option="B", direction = -1, name="Percent\nwith gene")

ggsave("CSH_tree_072126.png", width = 8.5, height = 6.4, dpi = 600, units = "in")



#### Full SNP tree

df_all_hits = df_GCF_nc_compressed %>%
  select(assembly:pylS, -nuccore) %>%
  `rownames<-`(df_GCF_nc_compressed$assembly) %>%
  select(-assembly)

df_RFs = data.frame(prfA = df_GCF_nc_compressed$prfA, 
                    prfB = df_GCF_nc_compressed$prfB, 
                    prfC = df_GCF_nc_compressed$prfC, 
                    prmC = df_GCF_nc_compressed$prmC, 
                    row.names = df_GCF_nc_compressed$assembly)

df_prfA = data.frame(prfA = df_GCF_nc_compressed$prfA,
                    row.names = df_GCF_nc_compressed$assembly)
df_prfB = data.frame(prfB = df_GCF_nc_compressed$prfB,
                     row.names = df_GCF_nc_compressed$assembly)
df_prfC = data.frame(prfC = df_GCF_nc_compressed$prfC,
                     row.names = df_GCF_nc_compressed$assembly)
df_prmC = data.frame(prmC = df_GCF_nc_compressed$prmC,
                     row.names = df_GCF_nc_compressed$assembly)

df_prfH = data.frame(prfH = df_GCF_nc_compressed$prfH,
                     row.names = df_GCF_nc_compressed$assembly)

df_selB = data.frame(selB = df_GCF_nc_compressed$selB,
                     row.names = df_GCF_nc_compressed$assembly)
df_pylS = data.frame(pylS = df_GCF_nc_compressed$pylS,
                     row.names = df_GCF_nc_compressed$assembly)

df_phylum = data.frame(df_GCF_nc_compressed$phylum, row.names = df_GCF_nc_compressed$assembly)

# Upload tree.
tree = read.newick("tree_out_j10_081225.tre")

# Fix GCF accession labels and drop two leaves that are duplicated phyla.
tree$tip.label = sub("_([A-Z]|_).*", "", tree$tip.label)


tree$tip.label = str_extract(tree$tip.label, "[^_]*_[^_]*")



p = ggtree(tree, layout="circular", linewidth=0.2)

p2 = gheatmap(p, df_prfA,  width=0.07, offset = 0, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("FALSE" = "white", "TRUE" = "dodgerblue")) + 
  new_scale_fill()

p3 = gheatmap(p2, df_prfB,  width=0.07, offset = 0.18, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("FALSE" = "white", "TRUE" = "springgreen2")) + 
  new_scale_fill()

p4 = gheatmap(p3, df_prfC,  width=0.07, offset = 0.36, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("FALSE" = "white", "TRUE" = "red3")) + 
  new_scale_fill()

p5 = gheatmap(p4, df_prmC,  width=0.07, offset = 0.54, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("FALSE" = "white", "TRUE" = "yellow3")) + 
  new_scale_fill()


p6 = gheatmap(p5, df_prfH,  width=0.07, offset = 0.72, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("FALSE" = "white", "TRUE" = "pink")) + 
  new_scale_fill()

p7 = gheatmap(p6, df_selB,  width=0.07, offset = 0.96, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("FALSE" = "white", "TRUE" = "orange2")) + 
  new_scale_fill()

p8 = gheatmap(p7, df_pylS,  width=0.07, offset = 1.14, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("FALSE" = "white", "TRUE" = "purple2")) + 
  new_scale_fill()



gheatmap(p8, df_phylum, offset=1.6, width=0.1, font.size=1, colnames = FALSE, color=NA) +
  theme(text=element_text(size=10)) +
  scale_fill_discrete(name = "Phylum/group", na.value = "white")
