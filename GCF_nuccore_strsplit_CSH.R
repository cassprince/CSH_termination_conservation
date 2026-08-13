library(tidyverse)

df = read.table("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//CSH//GCF_accessions_bact_072126.txt", sep = "\r", quote = "") %>%
  separate(col = V1, into = c("assembly", "nuccore"), sep = ">")
df$assembly = str_extract(df$assembly, "[^_]*_[^_]*")
df$assembly = str_match(df$assembly, "/(.*?)/")[,2]
df$nuccore = sub(" .*", "", df$nuccore)

write.csv(df, file = "C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//CSH//GCF_accessions_bact_clean_072126.csv")
