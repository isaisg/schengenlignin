library(ohchibi)



df <- read.table(file = "../rawdata/rnaseq_star_count.tsv",
                 header = F,
                 comment.char = "",quote = "")

Tab <- acast(data = df,formula = V2~V1,value.var = "V3")
#Create Tab
Tab <- Tab %>% t

rownames(Tab) <- rownames(Tab) %>%
  gsub(pattern = "_.*",replacement = "") %>%
  gsub(pattern = "-",replacement = "\\.")


Map <- read.table(file = "../rawdata/metadata.csv",header = T,sep = "\t")
Map$Name <- Map$Sample_Name %>% gsub(pattern = "_",replacement = "\\.")


Map <- match(rownames(Tab),Map$Name) %>%
  Map[.,]
Map$Id <- Map$Name %>% make.unique()

rownames(Tab) <- Map$Id
Tab <- Tab[,which(colSums(Tab) !=0)]


rownames(Map) <- Map$Id

Map$Genotype <- Map$Name %>% gsub(pattern = "\\..*",replacement = "") 
Map$Genotype[169:192] <- "ucc1.38ucc2GK"
Map$Genotype <- Map$Genotype %>% factor %>% relevel(ref = "WT")



Dat <- create_dataset(Tab = Tab %>% t,Map = Map)
saveRDS(object = Dat,file = "../cleandata/dat.star.RDS")
rm(list=ls())


