library(ohchibi)
library(DESeq2)
library(clusterProfiler)
library(org.At.tair.db)
library(extrafont)
loadfonts(device = "pdf")
library(paletteer)
library(scales)
library(ggtree)

Res <- readRDS(file = "../cleandata/design1_dat_deseq2.star.RDS")
res <- Res$Res_DESeq2
Dat_z <- Res$Dat_z
Dat_z <- Dat_z %>% subset.Dataset(Genotype != "ucc1") %>%
  subset.Dataset(Genotype != "ucc1.38ucc2GK",drop = T,clean = T)

df_genes <- read.table(file = "../cleandata/res_modelvswt_and_expression_genotypes_design1.star.csv",header = T,
           comment.char = "",quote = "",sep = ",")

##Append the family
df_fam <- read.table(file = "../rawdata/gene_families_sep_29_09_update.txt",
                     header = T,sep = "\t",quote = "",comment.char = "",fill = NA)
df_fam$Gene <- df_fam$Genomic_Locus_Tag %>% gsub(pattern = " ",replacement = "") %>%
  toupper()


peroxidase_locus <- df_fam$Gene_Family %>% grep(pattern = "peroxidase",value = F) %>%
  df_fam[.,] %$% Gene %>% as.character




##### Peroxidase #########
#Determine which genes are at least significantly different from WT in 1 genotype
mgenes <- which(df_genes$Gene %in% peroxidase_locus) %>%
  df_genes[.,] %>%
  subset(ClusterGeneAll == "C2" | ClusterGeneAll == "C1") %$% Gene %>%
  as.character %>% unique
mgenes <-  which(res$Gene %in% mgenes) %>%
  res[.,] %>%
  subset(padjG < 0.1) %>%
  droplevels %$% Gene %>% as.character %>%
   unique
res <- which(res$Gene %in% mgenes) %>%
  res[.,] %>% droplevels


###### Draw a heatmap with the chosen genes ##########
melted <- Dat_z$Tab %>% as.matrix %>% melt
colnames(melted) <- c("Gene","Name","value")
melted$Genotype <- melted$Name %>% gsub(pattern = "\\..*",replacement = "")

Tab_av <- acast(data = melted,formula = Gene~Genotype,
                fun.aggregate = mean,value.var  = "value")

Tab_av_sub <- which(rownames(Tab_av) %in% mgenes) %>%
  Tab_av[.,]


#Append the significance vs WT information
df_border <- res[,c("Genotype","Gene","padjG")]
df_border$Border <- NA
df_border$Border[which(df_border$padjG < 0.1)] <- "q < 0.1"
colnames(df_border)[1:2] <- c("IdCols","IdRows")
res_heatmap <- chibi.heatmap(Tab = Tab_av_sub,hclust_method_rows = "ward.D",df_border = df_border,
              hclust_method_cols = "ward.D",range_fill_heatmap = 2,
              axis_ticks_row = T,size_axis_text_row = 10,k_cols = 2,k_rows = 3,
              size_strip_text_row = 0,panel_spacing = 0.2)
p_pero <- res_heatmap$p_heatmap
melted_sub <- res_heatmap$melted

palette_heatmap = "pals::kovesi.diverging_bwr_55_98_c37"
range_fill_heatmap = 2
size_strip_text_row = 10
size_strip_text_col = 0
size_axis_text_col = 10
size_axis_text_row = 10
axis_ticks_row = F
size_legend_text = 10
size_legend_title = 12
size_border_tile = 0.5
width_border_tile = 0.85
height_border_tile = 0.85
palette_border = c("black")
size_dendrogram = 0.3
panel_border_heatmap = 0.3
panel_spacing = 0.2
font_family = "Helvetica"
font_face = "plain"
legend_position = "right"

melted_sub$IdCols <- melted_sub$IdCols %>%
  factor(levels = c("WT","sgn3","sgn3cif","esb1sgn3","sm","myb36","esb1","Wtcif"))
melted_sub$ClusterCols <- melted_sub$ClusterCols %>% factor(levels = c("CC2","CC1"))

#Append the common name for the locus
#Tair symbol
x <- org.At.tairSYMBOL
# Get the tair identifiers that are mapped to a gene symbol
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
xx <- lapply(X = xx,FUN = function(x)paste(x,collapse = "|"))
df_xx <- data.frame(Gene = names(xx),Symbol = xx %>% unlist %>% as.character) 
melted_sub$Symbol <- match(melted_sub$IdRows,df_xx$Gene) %>%
  df_xx$Symbol[.]
melted_sub$Nom <- paste0(melted_sub$IdRows,"(",melted_sub$Symbol,")") %>%
  gsub(pattern = "\\(NA\\)",replacement = "")
order_noms <- with(melted_sub,order(IdRows)) %>%
  melted_sub$Nom[.] %>% unique
melted_sub$Nom <- melted_sub$Nom %>% factor(levels = order_noms)

p_pero_end <- ggplot(data = melted_sub,mapping = aes(IdCols,Nom)) +
  geom_raster(aes(fill = value),color = "#00000000") +
  geom_tile(aes(color = Border),fill = '#00000000',
            size = size_border_tile,width = width_border_tile,height = height_border_tile) +
  facet_grid(ClusterRows~ClusterCols,scales = "free",space = "free") +
  theme_ohchibi(font_family = font_family
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand  = c(0,0)) +
  scale_fill_paletteer_c(palette_heatmap,
                         limits = c(-range_fill_heatmap,range_fill_heatmap),oob = squish,name = "z-score") +
  scale_color_manual(values = palette_border,na.value = "#00000000") +
  theme(
    axis.text.y = element_text(family = font_family,face = font_face,size =size_axis_text_row),
    axis.title.x = element_blank(),
    axis.text.x = element_text(family = font_family,face = font_face,size =size_axis_text_col ,angle = 45,vjust = 1,hjust = 1),
    strip.background.y = element_blank(),
    strip.text.y = element_text(size = size_strip_text_row,family = font_family,face = font_face,angle = 0),
    strip.background.x = element_blank(),
    strip.text.x = element_text(size = size_strip_text_col,family = font_family,face = font_face),
    panel.border = element_rect(size = panel_border_heatmap),
    axis.title.y = element_blank(),
    panel.spacing.y  = unit(panel_spacing, "lines"),
    panel.spacing.x  = unit(panel_spacing, "lines"),
    legend.position = legend_position ,
    legend.text = element_text(family = font_family,face = font_face,size = size_legend_text),
    legend.title = element_text(family = font_family,face = font_face,size = size_legend_title)
  )
    



##### Laccases ######
Res <- readRDS(file = "../cleandata/design1_dat_deseq2.star.RDS")
res <- Res$Res_DESeq2
Dat_z <- Res$Dat_z
Dat_z <- Dat_z %>% subset.Dataset(Genotype != "ucc1") %>%
  subset.Dataset(Genotype != "ucc1.38ucc2GK",drop = T,clean = T)

df_lac <- data.frame(Gene = c("AT1G18140","AT2G29130","AT2G30210","AT2G38080",
                              "AT2G40370","AT2G46570","AT3G09220","AT5G01040",
                              "AT5G01050","AT5G01190","AT5G03260","AT5G05390",
                              "AT5G07130","AT5G09360","AT5G48100","AT5G58910",
                              "AT5G60020"),
                     Name = c("LAC1","LAC2","LAC3","LAC4",
                              "LAC5","LAC6","LAC7","LAC8",
                              "LAC9","LAC10","LAC11","LAC12",
                              "LAC13","LAC14","LAC15","LAC16",
                              "LAC17"),
                     Type = "Laccase")

laccase_locus <- df_lac$Gene

mgenes <- which(df_genes$Gene %in% laccase_locus) %>%
  df_genes[.,] %>%
  subset(ClusterGeneAll == "C2" | ClusterGeneAll == "C1") %$% Gene %>%
  as.character %>% unique
mgenes <-  which(res$Gene %in% mgenes) %>%
  res[.,] %>%
  subset(padjG < 0.1) %>%
  droplevels %$% Gene %>% as.character %>%
  unique
res <- which(res$Gene %in% mgenes) %>%
  res[.,] %>% droplevels

###### Draw a heatmap with the chosen genes ##########
melted <- Dat_z$Tab %>% as.matrix %>% melt
colnames(melted) <- c("Gene","Name","value")
melted$Genotype <- melted$Name %>% gsub(pattern = "\\..*",replacement = "")

Tab_av <- acast(data = melted,formula = Gene~Genotype,
                fun.aggregate = mean,value.var  = "value")

Tab_av_sub <- which(rownames(Tab_av) %in% mgenes) %>%
  Tab_av[.,]

#Learn how to put asymetrical palettes
#Append the significance vs WT information
df_border <- res[,c("Genotype","Gene","padjG")]
df_border$Border <- NA
df_border$Border[which(df_border$padjG < 0.1)] <- "q < 0.1"
colnames(df_border)[1:2] <- c("IdCols","IdRows")

res_heatmap <- chibi.heatmap(Tab = Tab_av_sub,hclust_method_rows = "ward.D",df_border = df_border,
                             hclust_method_cols = "ward.D",range_fill_heatmap = 2,
                             axis_ticks_row = T,size_axis_text_row = 10,k_cols = 2,k_rows = 2,
                             size_strip_text_row = 0,panel_spacing = 0.2)
p_lac <- res_heatmap$p_heatmap

melted_sub <- res_heatmap$melted

melted_sub$IdCols <- melted_sub$IdCols %>%
  factor(levels = c("WT","sgn3","sgn3cif","esb1sgn3","sm","myb36","esb1","Wtcif"))
melted_sub$ClusterCols <- melted_sub$ClusterCols %>% factor(levels = c("CC2","CC1"))

#Append the common name for the locus
#Tair symbol
x <- org.At.tairSYMBOL
# Get the tair identifiers that are mapped to a gene symbol
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
xx <- lapply(X = xx,FUN = function(x)paste(x,collapse = "|"))
df_xx <- data.frame(Gene = names(xx),Symbol = xx %>% unlist %>% as.character) 
melted_sub$Symbol <- match(melted_sub$IdRows,df_xx$Gene) %>%
  df_xx$Symbol[.]
melted_sub$Nom <- paste0(melted_sub$IdRows,"(",melted_sub$Symbol,")") %>%
  gsub(pattern = "\\(NA\\)",replacement = "")
order_noms <- with(melted_sub,order(IdRows)) %>%
  melted_sub$Nom[.] %>% unique
melted_sub$Nom <- melted_sub$Nom %>% factor(levels = order_noms)

p_lac_end <- ggplot(data = melted_sub,mapping = aes(IdCols,Nom)) +
  geom_raster(aes(fill = value),color = "#00000000") +
  geom_tile(aes(color = Border),fill = '#00000000',
            size = size_border_tile,width = width_border_tile,height = height_border_tile) +
  facet_grid(ClusterRows~ClusterCols,scales = "free",space = "free") +
  theme_ohchibi(font_family = font_family
  ) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand  = c(0,0)) +
  scale_fill_paletteer_c(palette_heatmap,
                         limits = c(-range_fill_heatmap,range_fill_heatmap),oob = squish,name = "z-score") +
  scale_color_manual(values = palette_border,na.value = "#00000000") +
  theme(
    axis.text.y = element_text(family = font_family,face = font_face,size =size_axis_text_row),
    axis.title.x = element_blank(),
    axis.text.x = element_text(family = font_family,face = font_face,size =size_axis_text_col ,angle = 45,vjust = 1,hjust = 1),
    strip.background.y = element_blank(),
    strip.text.y = element_text(size = size_strip_text_row,family = font_family,face = font_face,angle = 0),
    strip.background.x = element_blank(),
    strip.text.x = element_text(size = size_strip_text_col,family = font_family,face = font_face),
    panel.border = element_rect(size = panel_border_heatmap),
    axis.title.y = element_blank(),
    panel.spacing.y  = unit(panel_spacing, "lines"),
    panel.spacing.x  = unit(panel_spacing, "lines"),
    legend.position = legend_position ,
    legend.text = element_text(family = font_family,face = font_face,size = size_legend_text),
    legend.title = element_text(family = font_family,face = font_face,size = size_legend_title)
  )




#Save plots ###
oh.save.pdf(p = p_pero_end,outname = "heatmap_peroxidases.pdf",
            outdir = "../figures/",width = 10,height = 12)
oh.save.pdf(p = p_lac_end,outname = "heatmap_laccases.pdf",
            outdir = "../figures/",width = 10,height = 12)
rm(list=ls())
dev.off()
gc()
