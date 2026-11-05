# Load packages and set working directory.
library(tidyverse)
library(treeio)
library(ggtree)
library(phangorn)
library(ggnewscale)
library(castor)
library(ggbreak)

setwd("C:/Users/princecr/Documents/GitHub/CSH_termination_conservation/data")

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

random = read.csv("random_genomes_dt_tree_051826.csv") %>%
  mutate(X = assembly) %>%
  filter(phylum != "Cyanobacteriota" & phylum != "Bacteroidota" & phylum != "Planctomycetota" & phylum != "Verrucomicrobiota")

random = column_to_rownames(random, 'X')

####### ---- Collapsed tree ---- #######

# Upload tree.
tree = read.newick("phyla_051826.tre")

# Fix GCF accession labels and drop two leaves that are duplicated phyla.
tree$tip.label = sub("_([A-Z]|_).*", "", tree$tip.label)
tree_trim = drop.tip(tree, tip = c("GCF_041929205.1", "GCF_046000445.1", "GCF_038400315.1", "GCF_033878955.1"))

# Midpoint root.
tree_mid = midpoint(tree_trim)

# Format tables for gene presence to append to the tree.
genes = dfShort %>%
  group_by(phylum) %>%
  summarize(n = n(), RF1 = sum(prfA)/n*100, RF2 = sum(prfB)/n*100, RF3 = sum(prfC)/n*100, PrmC = sum(prmC)/n*100, PrfH = sum(prfH)/n*100) %>%
  left_join(random %>% 
              mutate(phylum = gsub("Actinobacteriota", "Actinomycetota", phylum)) %>%
              select(phylum, assembly), by = join_by(phylum)) %>%
  filter(phylum != "Chlorobiota") %>%
  relocate(assembly) %>%
  as.data.frame()


rownames(genes) = genes$assembly


# Plot the tree and heatmap.
p = ggtree(tree_mid) %<+% genes + 
  xlim(0, 10) + 
  geom_tiplab(aes(label=phylum), size = 4) + 
  geom_nodepoint(aes(fill = as.numeric(label)*100), size = 2, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage", limits = c(0,100)) + 
  new_scale_fill()

gheatmap(p, genes[c("RF1", "RF2", "RF3", "PrmC")], offset = 5.2, width=1.2, font.size=3.5, colnames = TRUE, color=NA, colnames_angle = 45, colnames_offset_y = -0.2, colnames_offset_x = -0.15) + 
  scale_fill_viridis_c(option="B", direction = -1, name="Percent\nwith protein") +
  theme(legend.position = c(0.8,0.5))

ggsave("../figures/CSH_tree_110526.png", width = 7.5, height = 5.6, dpi = 600, units = "in")



gheatmap(p, genes[c("RF1", "RF2", "RF3", "PrmC", "PrfH")], offset = 5.2, width=1.6, font.size=3.5, colnames = TRUE, color=NA, colnames_angle = 45, colnames_offset_y = -0.2, colnames_offset_x = -0.15) + 
  scale_fill_viridis_c(option="B", direction = -1, name="Percent\nwith protein")  +
  theme(legend.position = c(0.82,0.5))

ggsave("../figures/CSH_tree_prfh_110526.png", width = 7.5, height =5.6, dpi = 600, units = "in")


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

tree_trim = drop.tip(tree, tip = c("GCF_007993755.1", "GCF_003258315.1", "GCF_007993775.1", "GCF_003260125.1", "GCF_004005565.1", "GCF_007997015.1", "GCF_016097415.1", "GCF_900128535.1", "GCF_918320885.1"))

p = ggtree(tree_trim, layout="circular", linewidth=0.15)

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
  new_scale_fill() +
  theme(legend.position = "none")


p6 = gheatmap(p5, df_prfH,  width=0.07, offset = 0.72, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("FALSE" = "white", "TRUE" = "orange2")) + 
  new_scale_fill() +
  theme(legend.position = "none")

p7 = gheatmap(p6, df_selB,  width=0.07, offset = 0.96, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("FALSE" = "white", "TRUE" = "orange2")) + 
  new_scale_fill()

p8 = gheatmap(p7, df_pylS,  width=0.07, offset = 1.14, colnames = FALSE, color=NA) +
  scale_fill_manual(values=c("FALSE" = "white", "TRUE" = "purple2")) + 
  new_scale_fill()


ggsave("../figures/CSH_fullphylo_prfH_110526.png", p6, width = 7.5, height = 6, dpi = 600, units = "in")

ggsave("../figures/CSH_fullphylo_noprfH_110526.png", p5, width = 7.5, height = 6, dpi = 600, units = "in")

p_phylum = gheatmap(p6, df_phylum, offset=1.6, width=0.1, font.size=1, colnames = FALSE, color=NA) +
  theme(text=element_text(size=10)) +
  scale_fill_discrete(name = "Phylum/group", na.value = "white") 

ggsave("../figures/CSH_fullphylo_phylogeny_110526.png", p_phylum, width = 10, height = 6, dpi = 600, units = "in")

