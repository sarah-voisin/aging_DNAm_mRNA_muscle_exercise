#Set working directory
setwd("C:/Users/e5103562/OneDrive - Victoria University/EWAS meta analysis of age/Muscle")

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
                lazy = FALSE)
  
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
    mutate_at(vars(EFFECTSIZE:SE_CORR),
              function(x){x*100})
    
  write.table(f,
              file=paste0('./',folder,'/',dataset,".tbl"),
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE,
              sep="\t")
}

sapply(files[10],
       bacon_process)

#Write code for METAL
#List tbl files
files <- list.files(pattern=".tbl",
                    recursive = TRUE)
file <- "METAL_commands.txt"
write.table(paste("SCHEME STDERR",
                  "MARKER CPG",
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

#Load METAL results and annotate
#Load annotation
annotation <- read.delim("Annotation.txt")

library(tidyverse)
meta_res <- read_tsv("METAL_muscle_age.TBL")
#Change name to probeID
meta_res <- dplyr::rename(meta_res,
                          probeID = MarkerName)

#Add CHR and position to meta analysis results
meta_res_tot <- inner_join(x = meta_res,
                           y = annotation,
                           by = "probeID")
max(meta_res_tot$HetDf)

#Calculate nb of CpGs included in analysis for each nb of studies
nb_CpGs <- c()
for (i in 1:max(meta_res_tot$HetDf))
{
  nb_CpGs <- c(nb_CpGs,nrow(meta_res_tot %>% filter(HetDf>=(i-1))))
}

cumsum <- tibble(`Number of included studies` = 1:max(meta_res_tot$HetDf),
                 `Number of CpGs present in all studies` = nb_CpGs)

ggplot(cumsum,
       mapping = aes(x = `Number of included studies`,
                     y = `Number of CpGs present in all studies`))+
  geom_point()+
  geom_line()+
  theme_classic()

#View results including at least 10 studies
meta_res_robust <- meta_res_tot %>% filter(HetDf>=9)
padj <- p.adjust(meta_res_robust$`P-value`,method="BH")
meta_res_robust <- meta_res_robust %>%
  mutate(`Adjusted p-value` = padj)
meta_res_robust <- meta_res_robust %>%
  mutate(sig = ifelse(`Adjusted p-value`<0.005,
                      "Sig",
                      "Not Sig"))
nrow(meta_res_robust) #595,541

#####Arrange
meta_res_robust <- meta_res_robust %>%
  arrange(desc(`P-value`))

#Select relevant columns only and rename them
meta_res_robust <- meta_res_robust %>%
  dplyr::select(probeID,
                CpG_chrm,
                CpG_beg,
                CGIposition,
                GeneHancer_interaction,
                genesUniq,
                Effect,
                StdErr,
                `P-value`,
                `Adjusted p-value`,
                HetDf,
                HetISq,
                HetPVal)
meta_res_robust <- meta_res_robust %>%
  mutate(genesUniq = replace_na(genesUniq,""),
         GeneHancer_interaction = replace_na(GeneHancer_interaction,""))

#Replace CpG island names with better words
meta_res_robust$CGIposition[grep("Shore",meta_res_robust$CGIposition)]="Shore"
meta_res_robust$CGIposition[grep("Shelf",meta_res_robust$CGIposition)]="Shelf"
meta_res_robust$CGIposition[is.na(meta_res_robust$CGIposition)]="Open sea"

#Change HetDf by number of included studies
meta_res_robust$HetDf <- meta_res_robust$HetDf +1
#Rename columns
colnames(meta_res_robust) = c("CpG",
                              "Chromosome",
                              "Position (hg38)",
                              "CpG island position",
                              "GeneHancer interaction",
                              "Annotated gene(s)",
                              "Effect size",
                              "SE",
                              "P-value",
                              "FDR",
                              "Number of studies",
                              "Heterogeneity index (I2)",
                              "Heterogeneity p-value")

write.table(meta_res_robust,
            file="METAL_muscle_age.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#Run random effects meta-analysis with metafor
library(metafor)
library(readxl)
dataset_summary <- read_excel('C:/Users/e5103562/OneDrive - Victoria University/DNAm datasets.xlsx',
                          na = c("","NA"))%>%
  dplyr::rename(prop_males = `% Male`)
meta_res_robust <- read_delim("METAL_muscle_age.txt",
                              lazy = FALSE)

#Create empty list
L.names <- meta_res_robust$CpG
L <- vector("list", length(L.names))
names(L) <- L.names

#Initiate list
tib <- tibble(CpG = meta_res_robust$CpG,
              Dataset = rep("Meta-analysis",nrow(meta_res_robust)),
              ES = meta_res_robust$`Effect size`,
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
names(L) <- meta_res_robust$CpG

#List tbl files
files <- list.files(pattern=".tbl",
                    recursive = TRUE)

#Add each study to each component of the list
for (s in files)
{
  file <- read_tsv(s,
                   lazy = F)
  dataset <- str_extract(s,
                         pattern = "(?<=/).*(?=\\.)")
  
  #Add FDR to the table
  file <- file %>%
    mutate(FDR = p.adjust(PVALUE_CORR))
  
  #Add Effect 
  #Create tibble
  summary <- tibble(CpG = file$CPG,
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
    filter(CpG %in% names(L))
  subL <- split(summary, seq(nrow(summary)))
  names(subL) = summary$CpG
  
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
setwd("C:/Users/e5103562/OneDrive - Victoria University/EWAS meta analysis of age/Muscle")
library(tidyverse)
L <- read_rds("METAL_muscle_age.rds")
library(metafor)
library(readxl)
dataset_summary <- read_excel('C:/Users/e5103562/OneDrive - Victoria University/DNAm datasets.xlsx',
                              na = c("","NA"))%>%
    dplyr::rename(prop_males = `% Male`,
                  Dataset = `Dataset ID`)%>%
    filter(!is.na(prop_males))%>%
    select(Dataset,prop_males)

meta_RE <- function(tib = tib)
{
  tib <- tib %>%
    filter(Type != "Meta-analysis")
  tib <- left_join(tib,
                   dataset_summary,
                   by = "Dataset")
  try({
    model <- rma(yi = ES,
                 sei = SE,
                 data=tib,
                 method = "EB",
                 mods = ~ prop_males,
                 control=list(stepadj=0.5,
                              maxiter=10000,
                              threshold=1e-8))
    return(cbind(as.numeric(model$b),
             model$se,
             model$pval))
  })
}

#Run on all CpGs
models <- lapply(L,
                 meta_RE)
CpG <- rownames(models)
save(models,
     file = "EWAS of age RE with sex as moderator.Rd")
#Load results
load("EWAS of age RE with sex as moderator.Rd")

#Extract sex effect
ES_RE <- sapply(models,
                FUN = function(m){return(m[2,1])})
SE_RE <- sapply(models,
                FUN = function(m){return(m[2,2])})
pval_RE <- sapply(models,
                FUN = function(m){return(m[2,3])})
CpG = names(models)

#models <- as_tibble(models)
#colnames(models) <- c("ES_RE","SE_RE","pval_RE")
models <- tibble(CpG = CpG, 
                 ES_RE = ES_RE,
                 SE_RE = SE_RE,
                 pval_RE = pval_RE) %>%
  mutate(FDR = p.adjust(pval_RE))
write.table(models,
            file="metafor_muscle_age_sex_RE.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")


#Pathway enrichment --> try 
setwd("C:/Users/e5103562/OneDrive - Victoria University/EWAS meta analysis of age/Muscle")
library(tidyverse)
#Load METAL results that have been formatted to annotate each CpG
to_write <- read_delim('METAL_muscle_age.txt')
to_order <- tibble(CpG = to_write$CpG)

#Load RE meta-analysis results from metafor for age
age <- read_delim('metafor_muscle_age_sex_RE.txt')
age <- left_join(to_order,
                 age)
to_write <- to_write %>%
  mutate(`Effect size` = age$ES_RE,
         SE = age$SE_RE,
         `P-value` = age$pval_RE,
         FDR = age$FDR)%>%
  dplyr::rename(`Additional % DNAm change per year of age in males` = `Effect size`,
                `P-value sex-age` = `P-value`,
                `FDR sex-age` = FDR)%>%
  arrange(`P-value sex-age`)%>%
  mutate(Significance = ifelse(`FDR sex-age`<0.005,
                               "Sig",
                               "Not Sig"))
sigdir <- to_write$Significance
sigdir[to_write$Significance=="Sig"&
         to_write$`Additional % DNAm change per year of age in males`<0]="Hypo"
sigdir[to_write$Significance=="Sig"&
         to_write$`Additional % DNAm change per year of age in males`>0]="Hyper"
to_write$Direction <- sigdir

#Replace NA with no name
to_write$`Annotated gene(s)`[is.na(to_write$`Annotated gene(s)`)]=""

library(missMethyl)
#OUR own annotation for each CpG
RSanno <- DataFrame(chr = to_write$Chromosome,
                    pos = to_write$`Position (hg38)`,
                    Name = to_write$CpG,
                    UCSC_RefGene_Name = to_write$`Annotated gene(s)`,
                    UCSC_RefGene_Group = to_write$`Annotated gene(s)`,
                    row.names = to_write$CpG)

#Download MSigDB file
if (!file.exists("C:/Users/e5103562/OneDrive - Victoria University/MSigDB.zip"))
{
    setwd("C:/Users/e5103562/OneDrive - Victoria University/")
    download.file("http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2022.1.Hs/msigdb_v2022.1.Hs_files_to_download_locally.zip",
                  destfile="MSigDB.zip")
    unzip("MSigDB.zip")
    file.rename(from = "msigdb_v2022.1.Hs_files_to_download_locally",
                to = "MSigDB")
}

#Load gene sets
library(mitch)
setwd("C:/Users/e5103562/OneDrive - Victoria University/")
genesets <- list.files(pattern = ".gmt",
                       recursive = TRUE)
cgp <- genesets[which(str_detect(genesets,
                                      fixed("cgp"))==TRUE)]
cp <- genesets[which(str_detect(genesets,
                                fixed("cp"))==TRUE)]
go <- genesets[which(str_detect(genesets,
                                fixed(".go."))==TRUE)]
hpo <- genesets[which(str_detect(genesets,
                                fixed("hpo"))==TRUE)]
c8 <- genesets[which(str_detect(genesets,
                                fixed("c8"))==TRUE)]

cgp <- gmt_import(cgp[which(str_detect(cgp,
                                       "entrez")==TRUE)])
cp <- gmt_import(cp[which(str_detect(cp,
                                     "biocarta|kegg|pid|reactome|wiki|symbols",
                                     negate = TRUE)==TRUE)])
go <- gmt_import(go[which(str_detect(go,
                                     ".bp|.cc|.mf|symbols",
                                     negate = TRUE)==TRUE)])
hpo <- gmt_import(hpo[which(str_detect(hpo,
                                       "entrez")==TRUE)])
c8 <- gmt_import(c8[which(str_detect(c8,
                                       "entrez")==TRUE)])
c8 <- c8[which(str_detect(names(c8),
                          "SKELETAL_MUSCLE")==TRUE)]
c8 <- c8[which(str_detect(names(c8),
                          "FETAL",
                          negate=TRUE)==TRUE)]
#add genes differentially expressed in fibre types to the list of marker genes
library(readxl)
fibre_marker_genes <- read_excel('Genes whose expression differs between type I and type II fibers in Rubenstein 2020.xlsx',
                                 skip = 1) %>%
  filter(padj < 0.005)
library(org.Hs.eg.db)
library(annotate)
fibre_marker_genes <- fibre_marker_genes %>%
  mutate(Entrez = mapIds(org.Hs.eg.db,
                         keys = Gene,
                         column = "ENTREZID",
                         keytype = "SYMBOL"))

c8$RUBENSTEIN_SKELETAL_MUSCLE_FIBRE_TYPE_SPECIFIC_GENES <-fibre_marker_genes$Entrez

library(org.Hs.eg.db)
library(annotate)
mRNA <- read_tsv('C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures/mRNA.txt')%>%
  mutate(Entrez = mapIds(org.Hs.eg.db,
                         keys = Gene,
                         column = "ENTREZID",
                         keytype = "SYMBOL"))

####Enrichment for each gene set
go_enrichment <- gsameth(#sig.cpg = to_write %>% filter(`FDR age`<0.005) %>%pull(CpG),
    sig.cpg = to_write %>% dplyr::slice_head(n = 200) %>%pull(CpG),           
    all.cpg = to_write$CpG,
                collection = go,
                array.type = "EPIC",
                plot.bias = FALSE,
                prior.prob = TRUE,
                anno = RSanno,
                equiv.cpg = TRUE,
                fract.counts = TRUE,
                sig.genes = TRUE) 
go_enrichment$ID <- rownames(go_enrichment)
go_enrichment <- go_enrichment %>%
    arrange(P.DE)
View(go_enrichment)
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
write.table(c8_enrichment,
            file="GSEA_meta-analysis RE_C8.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#Graph without info about genes inside circles
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
go_enrichment <- read_tsv("GSEA_meta-analysis RE_GO.txt")%>%
    mutate(Description = str_extract(ID,
                            "(?<=_)(.+)"))%>%
    mutate(Description = str_replace_all(Description,
                            "_",
                            " "))%>%
    mutate(Description = str_to_sentence(Description))

library(ggrepel)
library(pals)
p <- ggplot(data = go_enrichment,
            mapping = aes(x = 100*DE/N,
                          y = -log10(P.DE),
                          label = Description,
                          col = FDR<0.005
            ))+
    ylim(0,12)+
    geom_point()+
    scale_color_manual(values = c("Black","Red"))+
    geom_label_repel(data = go_enrichment %>% filter(FDR<0.02),
                     size = 2,
                     box.padding = unit(0.45, "lines"),
                     point.padding = unit(0.45, "lines"),
                     max.overlaps = 20)+
    labs(x = "% of genes in GO term that are\ndifferentially methylated",
         y = "Significance (-log10(pvalue))")+
    theme_bw()+
    theme(legend.position = "none")

setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
tiff('GSEA_meta-analysis RE_GO.tiff',
     width =100,
     height = 100,
     units = 'mm',
     res=600)
p
dev.off()

#HPO
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
hpo_enrichment <- read_tsv("GSEA_meta-analysis RE_HPO.txt")%>%
    mutate(Description = str_extract(ID,
                                     "(?<=_)(.+)"))%>%
    mutate(Description = str_replace_all(Description,
                                         "_",
                                         " "))%>%
    mutate(Description = str_to_sentence(Description))

library(ggrepel)
library(pals)
p <- ggplot(data = hpo_enrichment,
            mapping = aes(x = 100*DE/N,
                          y = -log10(P.DE),
                          label = Description,
                          col = FDR<0.005
            ))+
    #ylim(0,12)+
    geom_point()+
    scale_color_manual(values = c("Black","Red"))+
    geom_label_repel(data = hpo_enrichment %>% filter(FDR<0.05),
                     size = 2,
                     box.padding = unit(0.45, "lines"),
                     point.padding = unit(0.45, "lines"),
                     max.overlaps = 20)+
    labs(x = "% of genes in HPO term that are\ndifferentially methylated",
         y = "Significance (-log10(pvalue))")+
    theme_bw()+
    theme(legend.position = "none")

setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
tiff('GSEA_meta-analysis RE_HPO.tiff',
     width =100,
     height = 100,
     units = 'mm',
     res=600)
p
dev.off()

#Graph heatmap and circular plot
library(enrichplot)
library(mitch)
setwd("C:/Users/e5103562/OneDrive - Victoria University/")
genesets <- list.files(pattern = ".gmt",
                       recursive = TRUE)
hpo <- genesets[which(str_detect(genesets,
                                 fixed("hpo"))==TRUE)]
hpo <- gmt_import(hpo[which(str_detect(hpo,
                                       "symbol")==TRUE)])

#Create an enrichResult object for the code to work
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
hpo_enrichment <- read_tsv("GSEA_meta-analysis RE_HPO.txt")%>%
    mutate(Description = str_extract(ID,
                                     "(?<=_)(.+)"))%>%
    mutate(Description = str_replace_all(Description,
                                         "_",
                                         " "))%>%
    mutate(Description = str_to_sentence(Description))

hpo_enrichment <- hpo_enrichment %>%
    mutate(SigGenesInSet = str_replace_all(SigGenesInSet,
                                           ",",
                                           "/"),
           qvalue=FDR,
           GeneRatio=paste(DE,N,sep="/"))%>%
    dplyr::rename(pvalue=P.DE,
                  geneID=SigGenesInSet,
                  `p.adjust`=FDR,
                  Count=DE)%>%
    select(-c(N))
rownames(hpo_enrichment) <- hpo_enrichment$ID

library(qdapTools)
x <- new("enrichResult",
         result         = hpo_enrichment,
         pvalueCutoff   = 1,
         pAdjustMethod  = "BH",
         qvalueCutoff   = 0.005,
         gene           = unique(unlist(strsplit(hpo_enrichment$geneID,
                                                 ","))),
         universe       = list2df(hpo)$X1,
         geneSets       = hpo,
         organism       = "Homo sapiens",
         keytype        = "SYMBOL",
         ontology       = "HPO",
         readable       = FALSE
)


setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
tiff('GSEA_meta-analysis RE_HPO_circ.tiff',
     width =90,
     height = 130,
     units = 'mm',
     res=300)
cnetplot(x,
         showCategory=2,
         circular = TRUE,
         colorEdge = TRUE,
         cex_label_gene	= 0.5,
         show.legend = FALSE,
         cex_label_category	= 0.75) +
    theme(legend.position="none") 
dev.off()


#Cell type proportions as a dotplot
library(enrichplot)
library(mitch)
setwd("C:/Users/e5103562/OneDrive - Victoria University/")
genesets <- list.files(pattern = ".gmt",
                       recursive = TRUE)
c8 <- genesets[which(str_detect(genesets,
                                fixed("c8"))==TRUE)]
c8 <- gmt_import(c8[which(str_detect(c8,
                                     "entrez")==TRUE)])
c8 <- c8[which(str_detect(names(c8),
                          "SKELETAL_MUSCLE")==TRUE)]
c8 <- c8[which(str_detect(names(c8),
                          "FETAL",
                          negate=TRUE)==TRUE)]

#Create an enrichResult object for the code to work
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
c8_enrichment <- read_tsv("GSEA_meta-analysis RE_C8.txt")%>%
    mutate(Description = str_extract(ID,
                                     "(?<=_)(.+)"))%>%
    mutate(Description = str_replace_all(Description,
                                         "_",
                                         " "))%>%
    mutate(Description = str_replace_all(Description,
                                         "SKELETAL MUSCLE ",
                                         ""))%>%
    mutate(Description = str_to_sentence(Description))

c8_enrichment <- c8_enrichment %>%
    mutate(SigGenesInSet = str_replace_all(SigGenesInSet,
                                           ",",
                                           "/"),
           `Significance (-log10((p-value))` = -log10(P.DE),
           P.DE=1,
           qvalue=FDR,
           GeneRatio=paste(DE,N,sep="/"))%>%
    dplyr::rename(pvalue=P.DE,
                  geneID=SigGenesInSet,
                  `p.adjust`=FDR,
                  Count=DE)%>%
    dplyr::select(-c(N))
rownames(c8_enrichment) <- c8_enrichment$ID

library(qdapTools)
library(clusterProfiler)
x <- new("enrichResult",
         result         = c8_enrichment,
         pvalueCutoff   = 1,
         pAdjustMethod  = "BH",
         qvalueCutoff   = 1,
         gene           = unique(unlist(strsplit(c8_enrichment$geneID,
                                                 ","))),
         universe       = list2df(c8)$X1,
         geneSets       = c8,
         organism       = "Homo sapiens",
         keytype        = "SYMBOL",
         ontology       = "c8",
         readable       = FALSE
)

setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
tiff('GSEA_meta-analysis RE_C8.tiff',
     width =100,
     height = 150,
     units = 'mm',
     res=600)
dotplot(x,
        color = "pvalue",
        x = "`Significance (-log10((p-value))`",
        showCategory=30)+
    xlim(0,4) +
    geom_vline(xintercept = -log10(0.00041),
               lty=2)+
    theme(legend.position="none") 
dev.off()

