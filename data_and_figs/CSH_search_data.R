library(tidyverse)
library(janitor)
library(jsonlite)

setwd("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//CSH")

# Define function to read in HMMER files.
read_nhmmer = function(file){
  df = read.table(file, skip = 2, sep = "", colClasses = c("character", rep("NULL", 3), rep("character", 11), rep("NULL", 8)), fill = TRUE, row.names = NULL)
  df = df[,1:12]
  df = replace(df, df=='', NA)
  colnames(df) = c("accessions", "hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sqlen", "strand", "eval", "score", "bias")
  df = drop_na(df)
  df = df %>% mutate_at(c("hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sqlen", "eval", "score", "bias"), as.numeric)
  
  return(df)
}


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

# Upload assembly metadata for NCBI assemblies and select relevant data columns.

lines = readLines("metadata_081126.json")
lines = lapply(lines, fromJSON)
lines = lapply(lines, unlist)
df_metadata = bind_rows(lines)

df_metadata = df_metadata %>%
  select(accession, 
         assembly_name = assembly_info.assembly_name,
         assembly_level = assembly_info.assembly_level, 
         organism_name =  organism.organism_name, 
         taxid = organism.tax_id,
         n50 = assembly_stats.contig_n50, 
         gc_content = assembly_stats.gc_percent, 
         number_of_contigs = assembly_stats.number_of_contigs, 
         checkm_completeness = checkm_info.completeness, 
         checkm_contamination = checkm_info.contamination, 
         gene_counts.total = annotation_info.stats.gene_counts.total, 
         gene_counts.non_coding = annotation_info.stats.gene_counts.non_coding, 
         gene_counts.protein_coding = annotation_info.stats.gene_counts.protein_coding, 
         gene_counts.pseudogene = annotation_info.stats.gene_counts.pseudogene)

# Upload lineages acquired from NCBI taxdump and taxonkit based on each assembly's taxid.

lines = readLines("taxonomies_081126.json")
lines = lapply(lines, fromJSON)
lines = lapply(lines, unlist)
df_taxonomy = bind_rows(lines) 

df_taxonomy = df_taxonomy %>%
  select(taxid = taxonomy.tax_id, 
    kingdom = taxonomy.classification.kingdom.name,
    phylum = taxonomy.classification.phylum.name,
    class = taxonomy.classification.class.name,
    order = taxonomy.classification.order.name,
    family = taxonomy.classification.family.name,
    genus = taxonomy.classification.genus.name,
    species = taxonomy.classification.species.name,
    scientific_name = taxonomy.current_scientific_name.name
  )

df_GCF_tax_select = left_join(df_metadata, df_taxonomy, by = join_by(taxid)) 

# Upload NCBI RefSeq Assembly database accessions (beginning with "GCF_") and corresponding NCBI Nucleotide (nuccore) database accessions (beginning with "NC_" or "NZ_".

df_GCF_nc = read.csv("GCF_accessions_bact_clean_072126.csv")

# Join the NCBI accessions with their metadata.
df_GCF_nc_tax = left_join(df_GCF_nc, df_GCF_tax_select, by = join_by(assembly == accession)) %>%
  mutate(checkm_completeness = as.numeric(checkm_completeness)) %>%
  mutate(checkm_contamination = as.numeric(checkm_contamination)) 


# Import protein search data

df_prfA_gcf = read_tsv("prfA_GCF_WP_accs.tsv", col_names = c("assembly", "prot_acc"))
df_prfB_gcf = read_tsv("prfB_GCF_WP_accs.tsv", col_names = c("assembly", "prot_acc"))
df_prfC_gcf = read_tsv("prfC_GCF_WP_accs.tsv", col_names = c("assembly", "prot_acc"))
df_prfH_gcf = read_tsv("prfH_GCF_WP_accs.tsv", col_names = c("assembly", "prot_acc"))
df_prmC_gcf = read_tsv("prmC_GCF_WP_accs.tsv", col_names = c("assembly", "prot_acc"))
df_pylS_gcf = read_tsv("pylS_GCF_WP_accs.tsv", col_names = c("assembly", "prot_acc"))
df_selB_gcf = read_tsv("selB_GCF_WP_accs.tsv", col_names = c("assembly", "prot_acc"))

df = data.frame(cbind(df_GCF_nc_tax, prfA = (df_GCF_nc_tax$assembly %in% df_prfA_gcf$assembly), prfB = (df_GCF_nc_tax$assembly %in% df_prfB_gcf$assembly), prfC = (df_GCF_nc_tax$assembly %in% df_prfC_gcf$assembly), prfH = (df_GCF_nc_tax$assembly %in% df_prfH_gcf$assembly), prmC = (df_GCF_nc_tax$assembly %in% df_prmC_gcf$assembly), selB = (df_GCF_nc_tax$assembly %in% df_selB_gcf$assembly), pylS = (df_GCF_nc_tax$assembly %in% df_pylS_gcf$assembly)
))



# Multiple contigs, scaffolds, or chromosomes can make up an assembly. Identify gene hits across all contigs/scaffolds/chromosomes in each assembly. 
df_presence = df %>%
  group_by(assembly) %>%
  summarise(across(prfA:pylS, ~sum(.))) %>%
  mutate_if(is.numeric, ~1 * (. != 0))

bools = ifelse(df_presence[-1] == 1,"TRUE","FALSE")
bools = data.frame(cbind(assembly = df_presence$assembly, bools))

# Save the names of all Nucleotide accession numbers that were searched in each assembly.
names =  df %>%
  distinct(nuccore, .keep_all = TRUE) %>%
  group_by(assembly) %>%
  summarize(nuccore=paste(nuccore, collapse=";"))


# Prepare the final dataframe with Assembly and Nucleotide accession numbers, assembly statistics, taxonomy, and gene presence.
# Remove archaeal genomes. Filter based on CheckM completeness and contamination. Remove any duplicate genomes.
df_distinct = inner_join(names, bools, by = "assembly") %>%
  left_join(df[,1:21], by = "assembly", multiple = "any") %>%
  select(-X, -nuccore.y) %>%
  rename(nuccore = nuccore.x) %>% 
  filter(checkm_completeness > 80) %>%
  filter(checkm_contamination < 10) %>%
  distinct(assembly, .keep_all = TRUE) %>%
  drop_na(phylum) %>%
  mutate(across(c(prfA:pylS), as.logical))


# Rename taxonomy to be consistent with Coleman et al 2021? 
df_distinct$phylum[df_distinct$phylum == "delta/epsilon subdivisions"] = "Pseudomonadota"
df_distinct$phylum[df_distinct$phylum == "Pseudomonadota"] = "Proteobacteria"
df_distinct$phylum[df_distinct$phylum == "Terrabacteria group"] = df_distinct$class[df_distinct$phylum == "Terrabacteria group"]
df_distinct$phylum[df_distinct$phylum == "Bacillota"] = "Firmicutes"
#df_distinct$phylum[df_distinct$phylum == "Actinomycetota"] = "Actinobacteriota"
df_distinct$phylum[df_distinct$phylum == "Abditibacteriota"] = "Armatimonadota"
df_distinct$phylum[df_distinct$phylum == "Aquificota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Campylobacterota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Deferribacterota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Thermodesulfobacteriota"] = "Desulfuromonadota + Desulfobacterota"

df_distinct$phylum[df_distinct$phylum == "Calditrichota"] = "FCB group"
df_distinct$phylum[df_distinct$phylum == "Fibrobacterota"] = "FCB group"
df_distinct$phylum[df_distinct$phylum == "Bacteroidota"] = "FCB group"

df_distinct$phylum[df_distinct$phylum == "Verrucomicrobiota"] = "PVC group"
df_distinct$phylum[df_distinct$phylum == "Planctomycetota"] = "PVC group"
df_distinct$phylum[df_distinct$phylum == "Chlamydiota"] = "PVC group"

df_distinct$phylum[df_distinct$phylum == "Cyanobacteriota"] = "Cyanobacteriota/Melainabacteria group"

df_distinct$phylum[df_distinct$phylum == "Thermomicrobiota"] = "Chloroflexota"
df_distinct$phylum[df_distinct$phylum == "Proteobacteria"] = df_distinct$class[df_distinct$phylum == "Proteobacteria"]



df_distinct %>%
  summarize(n = n(), prfA = sum(prfA), prfB = sum(prfB), prfC = sum(prfC), prfH = sum(prfH), prmC = sum(prmC), pylS = sum(pylS), selB = sum(selB))

summary_table = df_distinct %>%
  group_by(phylum) %>%
  summarize(n = n(), prfA = sum(prfA), prfB = sum(prfB), prfC = sum(prfC), prfH = sum(prfH), prmC = sum(prmC), pylS = sum(pylS), selB = sum(selB))

freaks = df_distinct %>%
  filter(pylS == TRUE & selB == TRUE)

df_distinct %>%
  filter(pylS == TRUE & selB == FALSE)

# Write Table S#.
write.csv(df_distinct, file = "df_csh_081126.csv")


# ------ Save datafile with all accessions and taxonomies and metadata!

# df_distinct = left_join(names, df_GCF_nc_tax, by = "assembly", multiple = "first") %>%
#   select(-X, -nuccore.y, -domain.x, -domain.y, -organism) %>%
#   rename(nuccore = nuccore.x) %>%
#   drop_na(phylum)
# write.csv(df_distinct, file = "ref_genomes_metadata_072126.csv")


df_prfA = read_phmmer("prfA_prot_HMMER_5e2.csv")
df_prfB = read_phmmer("prfB_prot_HMMER_5e2.csv")
df_prfC = read_phmmer("prfC_prot_HMMER_5e2.csv")
df_prfH = read_phmmer("prfH_prot_HMMER_5e2.csv")
df_prmC = read_phmmer("prmC_prot_HMMER_5e2.csv")
df_selB = read_phmmer("selB_prot_HMMER_5e2.csv")
df_pylS = read_phmmer("pylS_prot_HMMER_5e2.csv")


df_prfA_filt = df_prfA %>%
  filter(full_score > 350)

df_prfB_filt = df_prfB %>%
  filter(full_score > 350)

df_prfC_filt = df_prfC %>%
  filter(full_score > 500)

df_prfH_filt = df_prfH %>%
  filter(full_score > 115) # getting lower scoring hits (~125) named "release factor-like proteins"

df_prmC_filt = df_prmC %>%
  filter(full_score > 200)

df_selB_filt = df_selB %>%
  filter(full_score > 220) # getting hits w score ~220 labelled "SelB containing proteins"... think more about the "dom" and other parts of phmmer output. what do these mean? are these relevant to me?

df_pylS_filt = df_pylS %>%
  filter(full_score > 150)



df_distinct %>%
  summarize(n = n(), prfB = sum(prfB))

no_prfA = df_distinct %>%
  filter(prfA == FALSE)

##### prfa workzone
prfA = read_nhmmer("prfA_HMMER_5e2.csv") 

ggplot(prfA, aes(x = score)) +
  geom_histogram()

prfA_gcf = prfA %>%
  left_join(df_GCF_nc, by = join_by(accessions == nuccore)) %>%
  select(-X) %>%
  distinct(assembly, .keep_all = TRUE) %>%
  mutate(alilen = hmmto - hmmfrom)

ggplot(prfA_gcf, aes(x = score)) +
  geom_histogram()

ggplot(prfA_gcf, aes(x = alilen)) +
  geom_histogram()

df_rf1 = data.frame(cbind(df_GCF_nc_tax, prfA_prot = (df_GCF_nc_tax$assembly %in% df_prfA_gcf$assembly), prfA_nuc = (df_GCF_nc_tax$assembly %in% prfA_gcf$assembly)
))


df_presence_rf1 = df_rf1 %>%
  group_by(assembly) %>%
  summarise(across(prfA_prot:prfA_nuc, ~sum(.))) %>%
  mutate_if(is.numeric, ~1 * (. != 0))

bools_rf1 = ifelse(df_presence_rf1[-1] == 1,"TRUE","FALSE")
bools_rf1 = data.frame(cbind(assembly = df_presence_rf1$assembly, bools_rf1))
names_rf1 =  df_rf1 %>%
  distinct(nuccore, .keep_all = TRUE) %>%
  group_by(assembly) %>%
  summarize(nuccore=paste(nuccore, collapse=";"))


# Prepare the final dataframe with Assembly and Nucleotide accession numbers, assembly statistics, taxonomy, and gene presence.
# Remove archaeal genomes. Filter based on CheckM completeness and contamination. Remove any duplicate genomes.
df_distinct_rf1 = inner_join(names_rf1, bools_rf1, by = "assembly") %>%
  left_join(df_rf1[,1:21], by = "assembly", multiple = "any") %>%
  select(-X, -nuccore.y) %>%
  rename(nuccore = nuccore.x) %>% 
  filter(checkm_completeness > 80) %>%
  filter(checkm_contamination < 10) %>%
  distinct(assembly, .keep_all = TRUE) %>%
  drop_na(phylum) %>%
  mutate(across(c(prfA_prot:prfA_nuc), as.logical))



no_prfa_1 = df_distinct_rf1 %>%
  filter(prfA_prot == FALSE)

no_prfa_2 = df_distinct_rf1 %>%
  filter(prfA_nuc == FALSE)



test = no_prfa_1 %>%
  left_join(prfA_gcf, by = join_by(assembly))


