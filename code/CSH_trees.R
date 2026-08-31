library(ggtree)
library(treeio)
library(tidyverse)
library(castor)
library(phangorn)
library(janitor)

setwd("C:/Users/cassp/Documents/GitHub/CSH_termination_conservation/data")


read_phmmer = function(file) {
  df = data.frame(readLines(file)) %>%
    slice(4:n()) %>%
    separate(readLines.file., into = c(as.character(1:34)), sep = " +") %>%
    select(-2,-3,-4) %>%
    replace(is.na(.), " ") %>% #Replace NAs with spaces.
    unite("description", 16:31, sep = " ") %>%
    mutate_if(is.character, str_trim) %>% #Trim spaces off description.
    filter(`1` != "#") %>%
    mutate(species = str_extract(description, "(?<=\\[).*(?=\\])")) %>% #Extract species name from between [] in description. 
    `colnames<-`(c("accession", "full_eval", "full_score", "full_bias", "dom_eval", "dom_score", "dom_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description", "species")) %>% #Rename columns.
    mutate_at(vars(full_eval, full_score, full_bias, dom_eval, dom_score, dom_bias), as.numeric)
  return(df)
}


tree = read.newick("prok_rf_tree_31.tree")

df = read_phmmer("prok_rf_prot_HMMER_1e5.csv")

df_tree = df %>%
  mutate(description = sub("isoform.*\\[", "[", description))%>%
  distinct(description, .keep_all = TRUE) %>%
  select(accession, description) %>%
  mutate(species = str_extract(description, "(?<=\\[).*(?=\\])")) %>%
  mutate(prot = str_extract(description, ".*(?=\\[)")) %>%
  mutate(location = str_extract(description, "mitochondrial")) %>%
  mutate(location = str_extract(description, "chloroplastic|mitochondrial"))


rownames(df_tree) = df_tree$accession

tree_one = get_subtree_with_tips(tree, only_tips = df_one$accession)$subtree

tree_mid = midpoint(tree_one)

#ggtree(tree_one, size=0.2) %<+% df_one + xlim(0,10) +
  #geom_tiplab(aes(label = description), size = 2)

#ggsave("prf_rect_tree.png", width = 9.5, height = 7)

p = ggtree(tree_mid, size=0.2) %<+% df_tree + xlim(0,10) +
  geom_tiplab(aes(label = species), size = 2)

gheatmap(p, df_tree %>% select(location),offset = 0, width=0.4, font.size=1.5, colnames = FALSE, color=NA)

ggsave("C:/Users/cassp/Documents/GitHub/CSH_termination_conservation/figures/prf_spp_mid_tree.png", width = 9.5, height = 7)





etree = read.newick("erf_tree_31.tree")

edf = read_phmmer("erf_prot_HMMER_1e5.csv")

edf_one = edf %>%
  mutate(description = sub("isoform.*\\[", "[", description)) %>%
  mutate(description = sub("Eukaryotic", "eukaryotic", description)) %>%
  mutate(species = str_extract(description, "(?<=\\[).*(?=\\])")) %>%
  distinct(description, .keep_all = TRUE) %>%
  select(accession, description, species)

etree_one = get_subtree_with_tips(etree, only_tips = edf_one$accession)$subtree

etree_mid = midpoint(etree_one)

ggtree(etree_one, size=0.2) %<+% edf_one + xlim(0,10) +
  geom_tiplab(aes(label = description), size = 2)

#ggsave("erf_rect_tree.png", width = 7, height = 4.5)

ggtree(etree_mid, size=0.2) %<+% edf_one + xlim(0,10) +
  geom_tiplab(aes(label = species), size = 2)

ggsave("C:/Users/cassp/Documents/GitHub/CSH_termination_conservation/figures/erf_spp_mid_tree.png", width = 7, height = 4.5)



################
ggtree(tree_one, size=0.2, layout = "equal_angle") %<+% df_one +
  geom_tiplab(aes(label = description), size = 2)


ggtree(tree_one, size=0.2, layout = "equal_angle")

ggsave("rf_equalangle_tree.png", width = 7, height = 5)


df_tree = df %>%
  select(accession, description) %>%
  mutate(short_description = sub("peptide chain release factor", "RF", description)) %>%
  mutate(short_description = sub("eukaryotic translation release factor", "eRF", short_description)) %>%
  mutate(short_description = sub("eukaryotic release factor", "eRF", short_description)) %>%
  mutate(short_description = sub("mitochondrial translation release factor", "mtRF", short_description))  %>%
  mutate(short_description = sub("RF 1, mitochondrial", "mtRF", short_description))




