#Set working directory
setwd("/Users/wrq796/OneDrive - Victoria University/TWAS meta analysis of age/")

#List tbl files
files <- list.files(pattern=".tbl",
           recursive = TRUE)

#Bacon process
library(bacon)
library(tidyverse)
bacon_process <- function(file = NULL)
{
  folder <- str_extract(file,
                        pattern = ".*(?=/)")
  dataset <- str_extract(file,
                         pattern = "(?<=/).*(?=\\.)")
  f <- read_tsv(file,
                show_col_types = FALSE,
                lazy = FALSE)

  if(max(f$EFFECTSIZE)>1)
      f <- f %>%
      mutate_at(vars(EFFECTSIZE,SE),
                function(x){x/100})

  bc <- bacon(NULL,
              f$EFFECTSIZE,
              f$SE)
  
  tiff(paste0('./',folder,'/QQ-plot_',dataset,'.tiff'),
       width =4,
       height = 2.5,
       units = 'in',
       res=200)
  print(plot(bc,
             type="qq"))
  #print(c(inflation(bc),bias(bc)))
  dev.off()
  
  f$EFFECTSIZE_CORR <- es(bc)[,1]
  f$SE_CORR <- se(bc)[,1]
  f$PVALUE_CORR <- pval(bc)[,1]
  
  #Multiply ES and SE by 100 because METAL rounds very small numbers
  f <- f %>%
    mutate_at(vars(EFFECTSIZE_CORR,SE_CORR),
              function(x){x*100})
    
  write.table(f,
              file=paste0('./',folder,'/',dataset,".tbl"),
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE,
              sep="\t")
}

sapply(files,
       bacon_process)

#Write code for METAL
#List tbl files
files <- list.files(pattern=".tbl",
                    recursive = TRUE)
file <- "METAL_commands.txt"
write.table(paste("SCHEME STDERR",
                  "MARKER GENE",
                  "EFFECT EFFECTSIZE_CORR",
                  "SEPARATOR TAB",
                  "STDERR SE_CORR",
                  "PVALUE PVALUE_CORR",
                  sep = "\r"),
            file = file,
            quote=F,
            row.names=F,
            col.names=F)
#Add datasets to analyse
add_datasets_METAl <- function(d = NULL,
                               f = NULL)
{
  write(paste("PROCESS",
              d,
              sep = "\t"),
        file=f,
        append=TRUE)
}
sapply(files,
       add_datasets_METAl,
       f = file)

#Add last line of code
write("ANALYZE HETEROGENEITY",
      file=file,
      append=TRUE)

#Load METAL results
library(tidyverse)
meta_res <- read_tsv("METAL_muscle_age.TBL")%>%
    dplyr::rename(`Gene name` = MarkerName)
max(meta_res$HetDf)

#Calculate nb of genes included in analysis for each nb of studies
nb_genes <- c()
for (i in 1:max(meta_res$HetDf))
{
  nb_genes <- c(nb_genes,nrow(meta_res %>% dplyr::filter(HetDf>=(i-1))))
}

cumsum <- tibble(`Number of included studies` = 1:max(meta_res$HetDf),
                 `Number of genes present in all studies` = nb_genes)

ggplot(cumsum,
       mapping = aes(x = `Number of included studies`,
                     y = `Number of genes present in all studies`))+
  geom_point()+
  geom_line()+
  theme_classic()

#View results including at least 15 studies
meta_res_robust <- meta_res %>% dplyr::filter(HetDf>=14)
padj <- p.adjust(meta_res_robust$`P-value`,method="BH")
meta_res_robust <- meta_res_robust %>%
  mutate(`Adjusted p-value` = padj)
meta_res_robust <- meta_res_robust %>%
  mutate(sig = ifelse(`Adjusted p-value`<0.005,
                      "Sig",
                      "Not Sig"))
nrow(meta_res_robust) #16,657

#####Arrange
sigdir <- meta_res_robust$sig
sigdir[meta_res_robust$sig=="Sig"&
         meta_res_robust$Effect<0]="Under"
sigdir[meta_res_robust$sig=="Sig"&
         meta_res_robust$Effect>0]="Over"
meta_res_robust <- meta_res_robust %>%
  mutate(sigdir = sigdir)
meta_res_robust <- meta_res_robust %>%
  arrange(desc(sigdir))

#Select relevant columns only and rename them
meta_res_robust <- meta_res_robust %>%
  dplyr::select(`Gene name`,
                Effect,
                StdErr,
                `P-value`,
                `Adjusted p-value`,
                HetDf,
                HetISq,
                HetPVal,
                sig,
                sigdir)%>%
    dplyr::rename(`logFC per year of age` = Effect,
                  SE = StdErr,
                  FDR = `Adjusted p-value`,
                  `Number of studies` = HetDf,
                  `Heterogeneity index (I2)` = HetISq,
                  `Heterogeneity p-value` = HetPVal,
                  Significance = sig,
                  Direction = sigdir)

nDMGs <- nrow(meta_res_robust %>% dplyr::filter(Significance=="Sig")) #how many DMPs = 1,961 DMPs ~ 12% of the transcriptome
nrow(meta_res_robust %>% dplyr::filter(Direction=="Under"))/nDMGs #% hypo DMPs = 51% hypo

write.table(meta_res_robust,
            file="METAL_muscle_age.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#Run fixed effects meta-analysis with metafor to compare results with METAL
library(metafor)
library(readxl)
dataset_summary <- read_excel('/Users/wrq796/OneDrive - Victoria University/mRNA datasets.xlsx',
                          na = c("","NA"),
                          n_max = 20)%>%
  dplyr::rename(prop_males = `% males`)
meta_res_robust <- read_delim("METAL_muscle_age.txt",
                              lazy = FALSE)

#Create empty list
L.names <- meta_res_robust$`Gene name`
L <- vector("list", length(L.names))
names(L) <- L.names

#Initiate list
tib <- tibble(`Gene name` = meta_res_robust$`Gene name`,
              Dataset = rep("Meta-analysis",nrow(meta_res_robust)),
              ES = meta_res_robust$`logFC per year of age`,
              SE = meta_res_robust$`SE`,
              PVAL = meta_res_robust$`P-value`,
              FDR = meta_res_robust$FDR,
              Type = rep("Meta-analysis",nrow(meta_res_robust)))
tib <- tib %>%
  mutate_if(is.numeric,
            signif,
            digits = 2)

#Initiate the list from meta-analysis
L <- split(tib, seq(nrow(tib))) #splits each row of the tibble into a list component
names(L) <- meta_res_robust$`Gene name`

#List tbl files
files <- list.files(pattern=".tbl",
                    recursive = TRUE)

#Add each study to each component of the list
for (s in files)
{
    print(s)
  file <- read_tsv(s,
                   lazy = F)
  dataset <- str_extract(s,
                         pattern = "(?<=/).*(?=\\.)")
  
  #Add FDR to the table
  file <- file %>%
    mutate(FDR = p.adjust(PVALUE_CORR))
  
  #Add Effect 
  #Create tibble
  summary <- tibble(`Gene name` = file$GENE,
                    Dataset = rep(dataset,nrow(file)),
                    ES = file$EFFECTSIZE_CORR,
                    SE = file$SE_CORR,
                    PVAL = file$PVALUE_CORR,
                    FDR = file$FDR,
                    Type = rep("Individual study",nrow(file)))
  summary <- summary %>%
    mutate_if(is.numeric,
              signif,
              digits = 2)
  
  summary <- summary %>%
    dplyr::filter(`Gene name` %in% names(L))
  subL <- split(summary, seq(nrow(summary)))
  names(subL) = summary$`Gene name`
  
  #Merge the pieces of the two lists that are in common
  L2 <- L[names(subL)]
  L3 <- L[setdiff(names(L),names(subL))]
  listinter <- Map(bind_rows,
                   L2,
                   subL)
  L <- c(listinter,
         L3)
}

saveRDS(L,
        file = "METAL_muscle_age.rds")

#Set working directory
setwd("/Users/wrq796/OneDrive - Victoria University/TWAS meta analysis of age/")
library(tidyverse)
L <- read_rds("METAL_muscle_age.rds")
library(metafor)
#Remove the NA from the list (I don't know how/why it got there)
L <- L[-which(sapply(L,is.null)==T)]

meta_RE <- function(tib = tib)
{
  tib <- tib %>%
    dplyr::filter(Type != "Meta-analysis")
  try({
    model <- rma(yi = ES,
                 sei = SE,
                 data=tib,
                 method = "EB",
                 control=list(stepadj=0.5,
                              maxiter=10000,
                              threshold=1e-8))
    return(c(as.numeric(model$b),
             model$se,
             model$pval))
  })
}

#Run on all genes
models <- t(sapply(L,
                 meta_RE))
gene <- rownames(models)
models <- as_tibble(models)
colnames(models) <- c("ES_RE","SE_RE","pval_RE")
models <- models %>%
  mutate(Gene = gene,
         FDR = p.adjust(pval_RE))
write.table(models,
            file="metafor_muscle_age_RE.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")


#Pathway enrichment
setwd("/Users/wrq796/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
#Load RE meta-analysis results from metafor for age
age <- read_delim('/Users/wrq796/OneDrive - Victoria University/TWAS meta analysis of age/metafor_muscle_age_RE.txt')
age <- age %>%
    dplyr::rename(`log2FC per year of age` = ES_RE,
                  SE = SE_RE,
                  `P-value age` = pval_RE,
                  `FDR age` = FDR)%>%
    arrange(`P-value age`)

library(clusterProfiler)
library(msigdbr)
cgp <- msigdbr(species = "Homo sapiens",
               category = "C2",
               subcategory = "CGP")%>% 
    dplyr::select(gs_name,
                  human_gene_symbol)

cp <- msigdbr(species = "Homo sapiens",
               category = "C2",
               subcategory = "CP")%>% 
    dplyr::select(gs_name,
                  human_gene_symbol)

go <- msigdbr(species = "Homo sapiens",
              category = "C5")%>%
    dplyr::filter(str_detect(gs_subcat,
                            "GO"))%>%
    dplyr::select(gs_name,
                  human_gene_symbol)
    
hpo <- msigdbr(species = "Homo sapiens",
              category = "C5",
              subcategory = "HPO")%>% 
    dplyr::select(gs_name,
                  human_gene_symbol)

c8 <- msigdbr(species = "Homo sapiens",
              category = "C8")%>%
    dplyr::filter(str_detect(gs_name,
                             "RUBENSTEIN"))%>%
    dplyr::select(gs_name,
                  human_gene_symbol)
#add genes differentially expressed in fibre types to the list of marker genes
library(readxl)
fibre_marker_genes <- read_excel('/Users/wrq796/OneDrive - Victoria University/Genes whose expression differs between type I and type II fibers in Rubenstein 2020.xlsx',
                                 skip = 1) %>%
  dplyr::filter(padj < 0.005) %>%
  dplyr::rename(human_gene_symbol = Gene) %>%
  dplyr::select(human_gene_symbol)

fibre_marker_genes <- fibre_marker_genes %>%
  mutate(gs_name = "RUBENSTEIN_SKELETAL_MUSCLE_FIBRE_TYPE_SPECIFIC_GENES")

c8 <- bind_rows(c8,
                fibre_marker_genes)

DNAm <- read_tsv('/Users/wrq796/OneDrive - Victoria University/Tables & Figures/DNAm.txt') %>%
  pull(`Annotated gene(s)`)
DNAm <- tibble(gs_name = "DNAm",
               human_gene_symbol = unique(unlist(strsplit(DNAm,
                 ";"))))%>%
  drop_na()
  
em <- enricher(age %>%
                 dplyr::filter(`FDR age` < 0.005) %>%
                   pull(Gene),
               universe = age %>%
                   pull(Gene),
               pvalueCutoff = 1,
               qvalueCutoff = 1,
               TERM2GENE=c8)
View(em@result)

write.table(em@result,
            file="ORA_meta-analysis RE_C8.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#Graph without info about genes inside circles
setwd("/Users/wrq796/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
go_enrichment <- read_tsv("ORA_meta-analysis RE_GO.txt")%>%
    mutate(Description = str_extract(ID,
                                     "(?<=_)(.+)"))%>%
    mutate(Description = str_replace_all(Description,
                                         "_",
                                         " "))%>%
    mutate(Description = str_to_sentence(Description),
           prop = 100*Count/as.numeric(str_extract(GeneRatio,
                                         "(?<=/)(.+)")))

library(ggrepel)
library(pals)
p <- ggplot(data = go_enrichment,
            mapping = aes(x = prop,
                          y = -log10(pvalue),
                          label = Description,
                          col = p.adjust<0.005
            ))+
    #ylim(0,12)+
    geom_point()+
    scale_color_manual(values = c("Black","Red"))+
    geom_label_repel(data = go_enrichment %>% filter(p.adjust<0.005),
                     size = 2,
                     box.padding = unit(0.45, "lines"),
                     point.padding = unit(0.45, "lines"),
                     max.overlaps = 20)+
    labs(x = "% of genes in GO term that are\ndifferentially expressed",
         y = "Significance (-log10(pvalue))")+
    theme_bw()+
    theme(legend.position = "none")

setwd("/Users/wrq796/OneDrive - Victoria University/Tables & Figures")
tiff('ORA_meta-analysis RE_GO.tiff',
     width =100,
     height = 100,
     units = 'mm',
     res=300)
p
dev.off()

#Try network
setwd("/Users/wrq796/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
mRNA <- read_tsv('mRNA.txt')

log2FC <- mRNA$`log2FC per year of age`
names(log2FC) <- mRNA$Gene

em@result <- em@result %>%
    mutate(Description = str_extract(ID,
                                     "(?<=_)(.+)"))%>%
    mutate(Description = str_replace_all(Description,
                                         "_",
                                         " "))%>%
    mutate(Description = str_to_sentence(Description))%>%
   mutate(Description =  stringr::str_wrap(Description, 20))

setwd("/Users/wrq796/OneDrive - Victoria University/Tables & Figures")
tiff('ORA_meta-analysis RE_GO_network.tiff',
     width =180,
     height = 170,
     units = 'mm',
     res=300)
cnetplot(em,
         showCategory = 6,
         foldChange=log2FC,
         cex_label_gene	= 0.6,
         cex_label_category	= 0.75)
dev.off()

#Graph without info about genes inside circles
setwd("/Users/wrq796/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
hpo_enrichment <- read_tsv("ORA_meta-analysis RE_HPO.txt")%>%
    mutate(Description = str_extract(ID,
                                     "(?<=_)(.+)"))%>%
    mutate(Description = str_replace_all(Description,
                                         "_",
                                         " "))%>%
    mutate(Description = str_to_sentence(Description),
           prop = 100*Count/as.numeric(str_extract(GeneRatio,
                                                   "(?<=/)(.+)")))

library(ggrepel)
library(pals)
p <- ggplot(data = hpo_enrichment,
            mapping = aes(x = prop,
                          y = -log10(pvalue),
                          label = Description,
                          col = p.adjust<0.005
            ))+
    #ylim(0,12)+
    geom_point()+
    scale_color_manual(values = c("Black","Red"))+
    geom_label_repel(data = hpo_enrichment %>% filter(p.adjust<0.05),
                     size = 2,
                     box.padding = unit(0.45, "lines"),
                     point.padding = unit(0.45, "lines"),
                     max.overlaps = 20)+
    labs(x = "% of genes in HPO term that are\ndifferentially expressed",
         y = "Significance (-log10(pvalue))")+
    theme_bw()+
    theme(legend.position = "none")

setwd("/Users/wrq796/OneDrive - Victoria University/Tables & Figures")
tiff('ORA_meta-analysis RE_HPO.tiff',
     width =100,
     height = 100,
     units = 'mm',
     res=300)
p
dev.off()

#Graph without info about genes inside circles
setwd("/Users/wrq796/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)

em <- enricher(age %>%
                   dplyr::filter(`FDR age` < 0.005) %>%
                   pull(Gene),
               universe = age %>%
                   pull(Gene),
               pvalueCutoff = 1,
               qvalueCutoff = 1,
               TERM2GENE=hpo)

#Try network
setwd("/Users/wrq796/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
mRNA <- read_tsv('mRNA.txt')

log2FC <- mRNA$`log2FC per year of age`
names(log2FC) <- mRNA$Gene

em@result <- em@result %>%
    mutate(Description = str_extract(ID,
                                     "(?<=_)(.+)"))%>%
    mutate(Description = str_replace_all(Description,
                                         "_",
                                         " "))%>%
    mutate(Description = str_to_sentence(Description))%>%
    mutate(Description =  stringr::str_wrap(Description, 20))

setwd("/Users/wrq796/OneDrive - Victoria University/Tables & Figures")
tiff('ORA_meta-analysis RE_HPO_circ.tiff',
     width =80,
     height = 80,
     units = 'mm',
     res=300)
cnetplot(em,
         showCategory = 2,
         circular = TRUE,
         colorEdge = TRUE,
         foldChange=log2FC,
         cex_label_gene	= 0.5,
         cex_label_category	= 0.75)+
    theme(legend.position="none") 
dev.off()

#Results for cell type
setwd("/Users/wrq796/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
#Add significance as -log10(pval)
em@result <- em@result %>%
  mutate(`Significance (-log10((p-value))` = -log10(pvalue))%>%
  mutate(Description = str_extract(ID,
                                   "(?<=_)(.+)"))%>%
  mutate(Description = str_replace_all(Description,
                                       "_",
                                       " "))%>%
  mutate(Description = str_replace_all(Description,
                                       "SKELETAL MUSCLE ",
                                       ""))%>%
  mutate(Description = str_to_sentence(Description),
         prop = 100*Count/as.numeric(str_extract(GeneRatio,
                                                 "(?<=/)(.+)")))

#Load logFC for smooth muscle cell marker genes
setwd("/Users/wrq796/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
mRNA <- read_tsv('mRNA.txt')

log2FC <- mRNA$`log2FC per year of age`
names(log2FC) <- mRNA$Gene

setwd("/Users/wrq796/OneDrive - Victoria University/Tables & Figures")
tiff('ORA_meta-analysis RE_C8.tiff',
     width =100,
     height = 150,
     units = 'mm',
     res=600)
dotplot(em,
        color = "pvalue",
        x = "`Significance (-log10((p-value))`",
        showCategory=30)+
  xlim(0,4) +
  geom_vline(xintercept = -log10(0.00045),
             lty=2)+
  theme(legend.position="none") 
dev.off()

#Add heatmap to show the increase/decrease of the marker genes
genes_to_extract <- em@result %>%
  slice(1)%>%
  pull(geneID)
genes_to_extract <- unlist(strsplit(genes_to_extract,
                                    "/"))

dat_heat <- tibble(Gene = genes_to_extract,
                   log2FC = log2FC[genes_to_extract]) %>%
  arrange(log2FC)
genes <- dat_heat$Gene
dat_heat <- dat_heat[,2]
rownames(dat_heat) <- genes

library("pheatmap")
library(grid)
library(RColorBrewer)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

tiff("Heatmap ORA smooth muscle cell markers.tiff",
     width =50,
     height = 100,
     units = 'mm',
     res=300)
pheatmap(dat_heat,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(100),
         #cutree_rows = 2,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = TRUE)
dev.off()
