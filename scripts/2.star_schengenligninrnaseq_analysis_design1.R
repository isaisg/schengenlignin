library(ohchibi)
library(DESeq2)
library(clusterProfiler)
library(org.At.tair.db)
library(extrafont)
loadfonts(device = "pdf")
library(paletteer)
library(scales)
library(ggtree)

Dat <- readRDS(file = "../cleandata/dat.star.RDS")

Dat$Map$Genotype %>% unique

### Subset dataset for specific contrasts ### 
#Create dds object 
dds <- DESeqDataSetFromMatrix(countData = Dat$Tab,
                              colData = Dat$Map,design = ~ Rep + Genotype)

dds <- DESeq2::collapseReplicates(object = dds,dds$Name)

dds <- DESeq(object = dds)


#Create object to plot
mvst <- vst(object = dds,blind = F)

#Remove the batch effect from the vst matrix
mat <- assay(mvst)

design0 <- model.matrix(~ 0 + Genotype, colData(dds))
mat <- limma::removeBatchEffect(mat, mvst$Rep,design = design0)

Tab_z <- mat %>% t %>% scale (center = T,scale = F)


Mapa <- Dat$Map[,c(1,2,3,4,6)] %>% unique
rownames(Mapa) <- Mapa$Name

Dat_z <- create_dataset(Tab = Tab_z %>% t,Map = Mapa)
Dat_z_all <- Dat_z 

######## Compute contrasts of each mutant versus Col-0 ####
mgenos <- dds$Genotype %>% levels %>% 
  grep(pattern = "WT",invert = T,value = T) %>%
  grep(pattern = "ucc",invert = T,value = T)

Res <- NULL
for(genotype in mgenos){
  mcontrast <- c("Genotype",genotype,"WT")
  df_temp <- results(object = dds,contrast = mcontrast) %>% as.data.frame
  df_temp$Gene <- rownames(df_temp)
  df_temp$Genotype <- genotype
  rownames(df_temp) <- NULL
  Res <-  rbind(Res,df_temp)
}

Res$padjG <- Res$pvalue %>% p.adjust(method = "fdr")

######## Save this two structures #########
general_list <- list(Dat_z = Dat_z_all,
                     Res_DESeq2 = Res)


saveRDS(object = general_list,file = "../cleandata/design1_dat_deseq2.star.RDS")

## Remove the ucc genotypes ####
Dat_z <- Dat_z %>% subset.Dataset(Genotype != "ucc1") %>%
  subset.Dataset(Genotype != "ucc1.38ucc2GK",drop = T,clean = T)

mpca <- oh.pca(Tab = Dat_z$Tab  %>% t,Map = Dat_z$Map,
               center = F,scale = F,id_var = "Name")

p <- chibi.pca(list_ohpca = mpca,col_val = "Genotype",size = 10,size_panel_border = 2)


oh.save.pdf(p = p,outname = "pca_design1.star.pdf",outdir =  "../figures/",width = 10,height = 8)

mgenes <- Res %>% subset(padjG < 0.05) %$% Gene %>% as.character %>% unique 

###### Draw a heatmap with the chosen genes ##########
melted <- Dat_z$Tab %>% as.matrix %>% melt
colnames(melted) <- c("Gene","Name","value")
melted$Genotype <- melted$Name %>% gsub(pattern = "\\..*",replacement = "")

Tab_av <- acast(data = melted,formula = Gene~Genotype,
                fun.aggregate = mean,value.var  = "value")
Tab_av_sub <- which(rownames(Tab_av) %in% mgenes) %>%
  Tab_av[.,]

nrow(Tab_av_sub)
length(mgenes)

## Cluster the genes and genotypes ###
dist_genes <- as.dist(1-cor(Tab_av_sub %>% t))
mclust_genes <- hclust(d = dist_genes,method = "ward.D2") 
mclust_genes %>% plot
order_genes <- mclust_genes$order %>% mclust_genes$labels[.]

df_clust_genes <- mclust_genes %>% cutree(k = 7) %>%
  data.frame(Gene = names(.), ClusterGene = paste0("C",.),row.names = NULL)

df_clust_genes <- df_clust_genes[,-1]

mclust_genotype <- dist(Tab_av_sub %>% t) %>%
  hclust(method = "ward.D2")
order_genotypes <- mclust_genotype$order %>% mclust_genotype$labels[.]

mclust_genotype %>% plot

df_clust_genotype <- mclust_genotype %>% cutree(k = 2) %>%
  data.frame(Genotype = names(.), ClusterGenotype = paste0("C",.),row.names = NULL)
df_clust_genotype <- df_clust_genotype[,-1]


melted_sub <- Tab_av_sub %>% melt
colnames(melted_sub) <- c("Gene","Genotype","value")

melted_sub <- merge(melted_sub,df_clust_genes, by = "Gene")
melted_sub <- merge(melted_sub,df_clust_genotype, by = "Genotype")

melted_sub$Gene <- melted_sub$Gene %>% 
  factor(levels = order_genes)
melted_sub$Genotype <- melted_sub$Genotype %>%
  factor(levels = order_genotypes)

order_groups_genes <- with(melted_sub,order(Gene)) %>%
  melted_sub$ClusterGene[.] %>% as.character %>%
  unique

melted_sub$ClusterGene <- melted_sub$ClusterGene %>%
  factor(levels = order_groups_genes %>% rev)


order_groups_genotypes <- with(melted_sub,order(Gene)) %>%
  melted_sub$ClusterGenotype[.] %>% as.character %>%
  unique

melted_sub$ClusterGenotype <- melted_sub$ClusterGenotype %>%
  factor(levels = order_groups_genotypes %>% rev)



melted_sub$value %>% sort %>% plot

p_heatmap <- ggplot(data = melted_sub,mapping = aes(Genotype,Gene)) + 
  geom_raster(aes(fill = value),color = "#00000000") + 
  #geom_tile(aes(color = Significance),fill = '#00000000',
  #          size = 1,width = 0.95,height = 0.95) + 
  facet_grid(ClusterGene~ClusterGenotype,scales = "free",space = "free") +
  theme_ohchibi() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_55_98_c37",
                         limits = c(-1,1),oob = squish,name = "z-score") +
  scale_color_manual(values = "black",na.value = "#00000000")+
  theme(
    #axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15,angle = 45,vjust = 0.5,hjust = 0.5),
    strip.background.y = element_blank(),
    strip.text.y = element_text(size = 15,family = "Arial",face = "bold",angle = 0),
    strip.background.x = element_blank(),
    strip.text.x = element_blank(),
    panel.border = element_rect(size = 1),
    axis.title.y = element_blank(),
    panel.spacing.y  = unit(0.2, "lines")
  )  

### Tree part ###
tree <- mclust_genotype %>% as.phylo
p_tree_geno  <- ggtree(tree,ladderize = F,size = 0.7) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.09,0.09)) +
  coord_flip()  + scale_x_reverse()

tree <- mclust_genes %>% as.phylo
p_tree_genes <- ggtree(tree,ladderize = F,size = 0.05) + 
  #The expand parth is fundamental to make the tree fit the composition with the other plots
  scale_y_continuous(expand = c(0.001,0.001))

p_blank <- ggplot() +theme_void()

###### Arrange figures ##########

composition <- egg::ggarrange(p_blank,p_tree_geno,
                              p_tree_genes,p_heatmap,
                              ncol = 2,nrow = 2,byrow = T,
                              heights = c(0.05,1),
                              widths = c(0.2,1),padding = unit(20,"line"))
oh.save.pdf(p = composition,outname = "composition_heatmap_design1.star.pdf",
            outdir = "../figures/",width = 8,height = 12)


### Create table ###
df_end <- merge(melted_sub,Res[,c(2,5,7,8,9)], by = c("Gene","Genotype"),all.x = TRUE) 
df_end <- df_end[,-c(5)]
colnames(df_end) <- c("Gene","Genotype","Expression","ClusterGeneAll","log2FoldChange_vs_WT","pvalue_vs_WT","padj_vs_WT")

#Append the gene information
### Mapped to gene symbol ######
x <- org.At.tairSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
df_gene_symbol <-lapply(X = xx,FUN = function(x){paste(x,collapse = "/")})  %>%
  do.call(rbind,.)
df_gene_symbol <- data.frame(Gene = rownames(df_gene_symbol),
                             Symbol = df_gene_symbol[,1],row.names = NULL)

df_end$GeneSymbol <- match(df_end$Gene,df_gene_symbol$Gene) %>%
  df_gene_symbol$Symbol[.]
df_end <- df_end[,c(1,2,4,3,5:8)]
df_end$Gene <- df_end$Gene %>% factor(levels = order_genes %>% rev)
df_end <- with(df_end,order(Gene)) %>%
  df_end[.,] 


## Gene ontology analysis
mlist <- list(
  C4 = df_clust_genes %>% subset(ClusterGene == "C4") %$% Gene %>% as.character,
  C6 = df_clust_genes %>% subset(ClusterGene == "C6") %$% Gene %>% as.character,
  C5 = df_clust_genes %>% subset(ClusterGene == "C5") %$% Gene %>% as.character,
  C1 = df_clust_genes %>% subset(ClusterGene == "C1") %$% Gene %>% as.character,
  C3 = df_clust_genes %>% subset(ClusterGene == "C3") %$% Gene %>% as.character,
  C7 = df_clust_genes %>% subset(ClusterGene == "C7") %$% Gene %>% as.character,
  C2 = df_clust_genes %>% subset(ClusterGene == "C2") %$% Gene %>% as.character
  
)

cg <- compareCluster(geneCluster=mlist,
                     fun="enrichGO",
                     keyType       = "TAIR",
                     OrgDb         = org.At.tair.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.1)

p <- dotplot(cg, showCategory=20, includeAll=TRUE)+
  scale_color_paletteer_c( "viridis::plasma",na.value = "#BFBFBF") 

oh.save.pdf(p = p,outname = "rnaseq_go_enrichments_top20.star.pdf",outdir = "../figures/",width = 10,height = 14)

### Shot high level enrichment ###
df_high <- gofilter(x = cg,level = 4) 
p <- dotplot(df_high, showCategory=20, includeAll=TRUE)+
  scale_color_paletteer_c( "viridis::plasma",na.value = "#BFBFBF") 
oh.save.pdf(p = p,outname = "rnaseq_go_enrichments_level4.star.pdf",outdir = "../figures/",width = 10,height = 14)


### Remove redundancy ####

bp2 <- simplify(cg, cutoff=0.7, by="p.adjust", select_fun=min)
p <- dotplot(bp2, showCategory=20, includeAll=TRUE)+
  scale_color_paletteer_c( "viridis::plasma",na.value = "#BFBFBF") 

oh.save.pdf(p = p,outname = "rnaseq_go_enrichments_top20.removedredundant.star.pdf",outdir = "../figures/",width = 10,height = 14)


write.table(x = df_end,file = "../cleandata/res_modelvswt_and_expression_genotypes_design1.star.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)


rm(list=ls())
dev.off()
