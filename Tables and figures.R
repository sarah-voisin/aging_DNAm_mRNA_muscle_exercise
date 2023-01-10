###############################################
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
#Load METAL results that have been formatted to annotate each CpG
to_write <- read_delim('C:/Users/e5103562/OneDrive - Victoria University/EWAS meta analysis of age/Muscle/METAL_muscle_age.txt')
to_order <- tibble(CpG = to_write$CpG)

#Load RE meta-analysis results from metafor for age
age <- read_delim('C:/Users/e5103562/OneDrive - Victoria University/EWAS meta analysis of age/Muscle/metafor_muscle_age_RE.txt')
age <- left_join(to_order,
         age)
to_write <- to_write %>%
 mutate(`Effect size` = age$ES_RE,
     SE = age$SE_RE,
     `P-value` = age$pval_RE,
     FDR = age$FDR)%>%
 dplyr::rename(`% DNAm change per year of age` = `Effect size`,
        `P-value age` = `P-value`,
        `FDR age` = FDR)%>%
 arrange(`P-value age`)

#Volcano plot for EWAS RE meta-analysis of age
to_write <- to_write %>%
    mutate(Significance = ifelse(`FDR age`<0.005,
                                 "Sig",
                                 "Not Sig"))
sigdir <- to_write$Significance
sigdir[to_write$Significance=="Sig"&
           to_write$`% DNAm change per year of age`<0]="Hypo"
sigdir[to_write$Significance=="Sig"&
           to_write$`% DNAm change per year of age`>0]="Hyper"
to_write$Direction <- sigdir
tiff('Volcano plot_EWAS RE meta-analysis of age.tiff',
   width =100,
   height = 80,
   units = 'mm',
   res=300)
ggplot() +
 geom_point(data = to_write,
       aes(`% DNAm change per year of age`,
         -log10(`P-value age`),
         col=Direction),
       size=0.5)+
 scale_color_manual(values=c("red","blue","black"))+
 labs(y="-log10(p-value)")+
 theme_classic()+
  theme(legend.position = "none")
dev.off()

DMPs <- subset(to_write,
               `FDR age`<0.005)

#Show histogram of effect sizes for DMPs
tiff('Distribution of effect size in DMPs_EWAS RE meta-analysis of age.tiff',
   width =100,
   height = 65,
   units = 'mm',
   res=300)
ggplot(data = DMPs,
    aes(x=`% DNAm change per year of age`,
      fill=Direction)) +
 geom_histogram(colour = "black")+
 scale_fill_manual(values=c("red","blue"))+
 lims(x=c(min(to_write$`% DNAm change per year of age`),
      max(to_write$`% DNAm change per year of age`)))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

#Forest plot for top hypo and hyper DMP
L <- read_rds("C:/Users/e5103562/OneDrive - Victoria University/EWAS meta analysis of age/Muscle/METAL_muscle_age.rds")
library(readxl)
dataset_summary <- read_excel('C:/Users/e5103562/OneDrive - Victoria University/DNAm datasets.xlsx',
                              na = c("","NA"))%>%
    dplyr::rename(prop_males = `% Male`,
                  Dataset = `Dataset ID`)%>%
    filter(!is.na(prop_males))%>%
    select(Dataset,prop_males,`n after preprocessing`)
library(metafor)
meta_RE <- function(tib = NULL)
{
 tib <- tib %>%
  filter(Type != "Meta-analysis")
 try({
  model <- rma(yi = ES,
         sei = SE,
         data=tib,
         mods = ~ prop_males,
         method = "EB",
         control=list(stepadj=0.5,
               maxiter=10000,
               threshold=1e-8))
  return(model)
 })
}

top_hypo <- DMPs %>%
 filter(`% DNAm change per year of age` < 0) %>%
 dplyr::slice(1)%>%
 pull(CpG)
tib <- left_join(L[[top_hypo]],
         dataset_summary)%>%
 arrange(-`n after preprocessing`)

top_hyper <- DMPs %>%
 filter(`% DNAm change per year of age` > 0) %>%
 dplyr::slice(1)%>%
 pull(CpG)

CpG <- "cg22833204"

tib <- left_join(L[[CpG]],
         dataset_summary)%>%
 arrange(-prop_males)

tiff('C:/Users/e5103562/OneDrive - Victoria University/Forest plot_sex moderator STAT1.tiff',
   width =150,
   height = 180,
   units = 'mm',
   res=300)
forest(meta_RE(tib),
    header="Dataset ID",
    xlab = "% DNAm change per year of age",
    digits=c(2L,4L),
    slab = tib %>%
     filter(Dataset != "Meta-analysis")%>%
     pull(Dataset)
)
dev.off()

#Plot ES as a function of % males
library(ggrepel)
tiff('C:/Users/e5103562/OneDrive - Victoria University/ES vs prop males_sex moderator STAT1.tiff',
     width =150,
     height = 100,
     units = 'mm',
     res=300)
ggplot(tib,
       mapping = aes(x=prop_males,
                     y=ES,
                     col = prop_males))+
    geom_point(size=3)+
    geom_smooth(method="lm")+
    labs(x = "% of males in dataset",
         y = "% DNAm change per year of age")+
    geom_text_repel(aes(label = Dataset),
                    col = "black",
                     #size = 2,
                     #box.padding = unit(0.45, "lines"),
                     #point.padding = unit(0.45, "lines"),
                     max.overlaps = 30)+
    theme_bw()

dev.off()

#Only keep age-related DMPs
to_write <- DMPs

#Load RE meta-analysis results from metafor for VO2max
VO2max <- read_delim('C:/Users/e5103562/OneDrive - Victoria University/EWAS meta analysis of VO2max/metafor_muscle_VO2max_RE.txt')%>%
 dplyr::select(-FDR)%>%
 dplyr::rename(`% DNAm change per unit of VO2max` = ES_RE,
        `SE VO2max` = SE_RE,
        `P-value VO2max` = pval_RE)
to_write <- left_join(to_write,
          VO2max)%>%
 mutate(`FDR VO2max` = p.adjust(`P-value VO2max`))

#QQplot to show inflation of p-values
library(ggrepel)
ci = 0.95
log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

toplot <- to_write %>%
  dplyr::filter(!is.na(`P-value VO2max`))%>%
    arrange(`P-value VO2max`)
n <- nrow(toplot)

toplot <- toplot %>%
  mutate(Consistency = ifelse(sign(`% DNAm change per year of age`)==sign(`% DNAm change per unit of VO2max`),
                "Aging",
                "Rejunevating"),
         observed = -log10(`P-value VO2max`),
         expected = -log10(ppoints(n)),
         clower  = -log10(qbeta(p = (1 - ci) / 2,
                                  shape1 = 1:n,
                                  shape2 = n:1)),
         cupper  = -log10(qbeta(p = (1 + ci) / 2,
                                  shape1 = 1:n,
                                  shape2 = n:1)))

tiff('QQplot_VO2max_DNAm.tiff',
   width =150,
   height = 150,
   units = 'mm',
   res=300)
ggplot(toplot) +
    geom_ribbon(
        mapping = aes(x = expected,
                      ymin = clower,
                      ymax = cupper),
        alpha = 0.1
    ) +
    geom_point(aes(expected,
                   observed,
                   col = Consistency),
               #shape = 2,
               size = 2) +
    scale_color_manual(values = c("Black","Green"))+
    geom_abline(intercept = 0,
                slope = 1,
                alpha = 0.5) +
    geom_label_repel(data = toplot %>% dplyr::filter(`FDR VO2max` < 0.005),
                     aes(expected,
                         observed,
                         label = `Annotated gene(s)`),
                     #size = 2,
                     #box.padding = unit(0.45, "lines"),
                     #point.padding = unit(0.45, "lines"),
                     max.overlaps = 30
                     )+
    ylim(0,25)+
    xlab(log10Pe) +
    ylab(log10Po)+
    theme_bw(base_size = 24) +
    theme(axis.ticks = element_line(size = 0.5),
          panel.grid = element_blank(),
          legend.position = "none")
dev.off()

#Graph same QQ plot, but separately for hypo and hyper DMPs, to see if CRF preferably reverses hypo/hyper DNAm with age
toplot <- to_write %>%
  dplyr::filter(!is.na(`P-value VO2max`) &
                  Direction=="Hypo")%>%
  arrange(`P-value VO2max`)
n <- nrow(toplot) #2,296 hypoDMPs

toplot <- toplot %>%
  mutate(Consistency = ifelse(sign(`% DNAm change per year of age`)==sign(`% DNAm change per unit of VO2max`),
                              "Aging",
                              "Rejunevating"),
         observed = -log10(`P-value VO2max`),
         expected = -log10(ppoints(n)),
         clower  = -log10(qbeta(p = (1 - ci) / 2,
                                shape1 = 1:n,
                                shape2 = n:1)),
         cupper  = -log10(qbeta(p = (1 + ci) / 2,
                                shape1 = 1:n,
                                shape2 = n:1)))

tiff('QQplot_VO2max_DNAm_hypoDMPs.tiff',
     width =150,
     height = 150,
     units = 'mm',
     res=300)
ggplot(toplot) +
  geom_ribbon(
    mapping = aes(x = expected,
                  ymin = clower,
                  ymax = cupper),
    alpha = 0.1
  ) +
  geom_point(aes(expected,
                 observed,
                 col = Consistency),
             #shape = 2,
             size = 2) +
  scale_color_manual(values = c("Black","Green"))+
  geom_abline(intercept = 0,
              slope = 1,
              alpha = 0.5) +
  geom_label_repel(data = toplot %>% dplyr::filter(`FDR VO2max` < 0.005),
                   aes(expected,
                       observed,
                       label = `Annotated gene(s)`),
                   #size = 2,
                   #box.padding = unit(0.45, "lines"),
                   #point.padding = unit(0.45, "lines"),
                   max.overlaps = 30
  )+
  ylim(0,25)+
  xlim(0,3.5)+
  xlab(log10Pe) +
  ylab(log10Po)+
  theme_bw(base_size = 24) +
  theme(axis.ticks = element_line(size = 0.5),
        panel.grid = element_blank(),
        legend.position = "none")
dev.off()

#Load RE meta-analysis results from metafor for training
training <- read_delim('C:/Users/e5103562/OneDrive - Victoria University/EWAS meta analysis of training/metafor_muscle_training_RE.txt')%>%
 dplyr::select(-FDR)%>%
 dplyr::rename(`% DNAm change after training` = ES_RE,
        `SE training` = SE_RE,
        `P-value training` = pval_RE)
to_write <- left_join(to_write,
          training)%>%
    mutate(`FDR training` = p.adjust(`P-value training`))

#QQplot to show inflation of p-values
library(ggrepel)
ci = 0.95
log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

toplot <- to_write %>%
    dplyr::filter(!is.na(`P-value training`)&
                    Direction == "Hyper")%>%
    arrange(`P-value training`)
n <- nrow(toplot)

toplot <- toplot %>%
    mutate(Consistency = ifelse(sign(`% DNAm change per year of age`)==sign(`% DNAm change after training`),
                                "Aging",
                                "Rejunevating"),
           observed = -log10(`P-value training`),
           expected = -log10(ppoints(n)),
           clower  = -log10(qbeta(p = (1 - ci) / 2,
                                  shape1 = 1:n,
                                  shape2 = n:1)),
           cupper  = -log10(qbeta(p = (1 + ci) / 2,
                                  shape1 = 1:n,
                                  shape2 = n:1)))

tiff('QQplot_training_DNAm_hyperDMPs.tiff',
     width =150,
     height = 150,
     units = 'mm',
     res=300)
ggplot(toplot) +
  geom_ribbon(
    mapping = aes(x = expected,
                  ymin = clower,
                  ymax = cupper),
    alpha = 0.1
  ) +
  geom_point(aes(expected,
                 observed,
                 col = Consistency),
             #shape = 2,
             size = 2) +
  scale_color_manual(values = c("Black","Green"))+
  geom_abline(intercept = 0,
              slope = 1,
              alpha = 0.5) +
  ylim(0,5)+
  xlim(0,3.5)+
  xlab(log10Pe) +
  ylab(log10Po)+
  theme_bw(base_size = 24) +
  theme(axis.ticks = element_line(size = 0.5),
        panel.grid = element_blank(),
        legend.position = "none")
dev.off()


#Remove useless columns
to_write <- to_write %>%
    select(-c(`Number of studies`:Direction))

#Write results
write_tsv(x = to_write,
          na = "",
      file = "DNAm.txt")


#########################################
#Age-related mRNA changes
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
#Load RE meta-analysis results from metafor for age
age <- read_delim('C:/Users/e5103562/OneDrive - Victoria University/TWAS meta analysis of age/metafor_muscle_age_RE.txt')
age <- age %>%
    dplyr::rename(`log2FC per year of age` = ES_RE,
                  SE = SE_RE,
                  `P-value age` = pval_RE,
                  `FDR age` = FDR)%>%
    dplyr::mutate(`log2FC per year of age` = `log2FC per year of age`/100,
           SE = SE/100)%>%
    arrange(`P-value age`)

#Volcano plot for EWAS RE meta-analysis of age
age <- age %>%
    mutate(Significance = ifelse(`FDR age`<0.005,
                                 "Sig",
                                 "Not Sig"))
sigdir <- age$Significance
sigdir[age$Significance=="Sig"&
           age$`log2FC per year of age`<0]="Under"
sigdir[age$Significance=="Sig"&
           age$`log2FC per year of age`>0]="Over"
age$Direction <- sigdir
tiff('Volcano plot_TWAS RE meta-analysis of age.tiff',
     width =100,
     height = 80,
     units = 'mm',
     res=300)
ggplot() +
    geom_point(data = age,
               aes(`log2FC per year of age`,
                   -log10(`P-value age`),
                   col=Direction),
               size=0.5)+
    scale_color_manual(values=c("black","red","blue"))+
    labs(y="-log10(p-value)")+
    theme_classic()+
    theme(legend.position = "none")
dev.off()

#Show histogram of effect sizes for DMPs
DEGs <- subset(age,`FDR age` < 0.005)
tiff('Distribution of effect size in DEGs_TWAS RE meta-analysis of age.tiff',
     width =100,
     height = 65,
     units = 'mm',
     res=300)
ggplot(data = DEGs,
       aes(x=`log2FC per year of age`,
           fill=Direction)) +
    geom_histogram(colour = "black")+
    scale_fill_manual(values=c("red","blue"))+
    lims(x=c(min(age$`log2FC per year of age`),
             max(age$`log2FC per year of age`)))+
    theme_classic()+
    theme(legend.position = "none")
dev.off()

#Forest plot for top hypo and hyper DMP
L <- read_rds("C:/Users/e5103562/OneDrive - Victoria University/TWAS meta analysis of age/METAL_muscle_age.rds")
#Remove the NA from the list (I don't know how/why it got there)
L <- L[-which(sapply(L,is.null)==T)]

library(readxl)
dataset_summary <- read_excel('C:/Users/e5103562/OneDrive - Victoria University/mRNA datasets.xlsx',
                              na = c("","NA"),
                              n_max = 20)%>%
    dplyr::rename(prop_males = `% males`)
library(metafor)
meta_RE <- function(tib = NULL)
{
    tib <- tib %>%
        filter(Type != "Meta-analysis")
    try({
        model <- rma(yi = ES/100,
                     sei = SE/100,
                     data=tib,
                     method = "EB",
                     control=list(stepadj=0.5,
                                  maxiter=10000,
                                  threshold=1e-8))
        return(model)
    })
}

top_hypo <- DEGs %>%
    dplyr::filter(`log2FC per year of age` < 0) %>%
    dplyr::slice(1)%>%
    pull(Gene)
tib <- left_join(L[[top_hypo]],
                 dataset_summary)%>%
    arrange(-`n after filtering`)

top_hyper <- DEGs %>%
    dplyr::filter(`log2FC per year of age` > 0) %>%
    dplyr::slice(1)%>%
    pull(Gene)
tib <- left_join(L[[top_hyper]],
                 dataset_summary)%>%
    arrange(-`n after filtering`)

tiff('Forest plot_top_over.tiff',
     width =150,
     height = 180,
     units = 'mm',
     res=300)
forest(meta_RE(tib),
       header="Dataset ID",
       xlab = "log2FC per year of age",
       digits=c(2L,4L),
       slab = tib %>%
           dplyr::filter(Dataset != "Meta-analysis")%>%
           pull(Dataset)
)
dev.off()

#Load RE meta-analysis results from metafor for Vo2max
VO2max <- read_delim('C:/Users/e5103562/OneDrive - Victoria University/TWAS meta analysis of VO2max/metafor_muscle_VO2max_RE.txt')%>%
    dplyr::select(-FDR)%>%
    mutate(ES_RE = ES_RE/100,
           SE_RE = SE_RE/100)%>%
    dplyr::rename(`log2FC per unit of VO2max` = ES_RE,
                  `SE VO2max` = SE_RE,
                  `P-value VO2max` = pval_RE)
DEGs <- left_join(DEGs,
                      VO2max)%>%
    mutate(`FDR VO2max` = p.adjust(`P-value VO2max`))

#QQplot to show inflation of p-values
library(ggrepel)
ci = 0.95
log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

toplot <- DEGs %>%
    dplyr::filter(!is.na(`P-value VO2max`))%>%
    arrange(`P-value VO2max`)
n <- nrow(toplot)

toplot <- toplot %>%
    mutate(Consistency = ifelse(sign(`log2FC per year of age`)==sign(`log2FC per unit of VO2max`),
                                "Aging",
                                "Rejunevating"),
           observed = -log10(`P-value VO2max`),
           expected = -log10(ppoints(n)),
           clower  = -log10(qbeta(p = (1 - ci) / 2,
                                  shape1 = 1:n,
                                  shape2 = n:1)),
           cupper  = -log10(qbeta(p = (1 + ci) / 2,
                                  shape1 = 1:n,
                                  shape2 = n:1)))

tiff('QQplot_VO2max_mRNA.tiff',
     width =150,
     height = 150,
     units = 'mm',
     res=300)
ggplot(toplot) +
    geom_ribbon(
        mapping = aes(x = expected,
                      ymin = clower,
                      ymax = cupper),
        alpha = 0.1
    ) +
    geom_point(aes(expected,
                   observed,
                   col = Consistency),
               #shape = 2,
               size = 2) +
    scale_color_manual(values = c("Black","Green"))+
    geom_abline(intercept = 0,
                slope = 1,
                alpha = 0.5) +
    geom_label_repel(data = toplot %>% dplyr::filter(`FDR VO2max` < 0.005),
                     aes(expected,
                         observed,
                         label = Gene),
                     #size = 2,
                     #box.padding = unit(0.45, "lines"),
                     #point.padding = unit(0.45, "lines"),
                     max.overlaps = 10
    )+
    ylim(0,25)+
    xlab(log10Pe) +
    ylab(log10Po)+
    theme_bw(base_size = 24) +
    theme(axis.ticks = element_line(size = 0.5),
          panel.grid = element_blank(),
          legend.position = "none")
dev.off()


#Load ExtraMETA results for training
extrameta <- read_tsv("C:/Users/e5103562/OneDrive - Victoria University/TWAS meta analysis of training/Extrameta all genes results.txt")%>%
    dplyr::filter(Gene %in% DEGs$Gene)%>%
    dplyr::rename(`log2FC after training` = `logFC after training`)%>%
    mutate(`FDR training` = p.adjust(`P-value training`))%>%
    dplyr::select(Gene,
           `log2FC after training`,
           `SE training`,
           `P-value training`,
           `FDR training`)

#library(readxl)
#AT <- read_excel('C:/Users/e5103562/OneDrive - Victoria University/TWAS meta analysis of AT/MetaMEx genes summary.xlsx',
#                 na = "NA")%>%
#    dplyr::filter(`Gene symbol` %in% DEGs$Gene)%>%
#    dplyr::select(`Gene symbol`,
#           estimate_TA,
#           SE_TA,
#           pval_TA)%>%
#    dplyr::rename(Gene = `Gene symbol`,
#           `log2FC after AT` = estimate_TA,
#           `SE AT` = SE_TA,
#           `P-value AT` = pval_TA)%>%
#    mutate(`FDR AT` = p.adjust(`P-value AT`))

DEGs <- left_join(DEGs,
                 extrameta)

#QQplot to show inflation of p-values
library(ggrepel)
ci = 0.95
log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

toplot <- DEGs %>%
    dplyr::filter(!is.na(`P-value training`))%>%
    arrange(`P-value training`)
n <- nrow(toplot)

toplot <- toplot %>%
    mutate(Consistency = ifelse(sign(`log2FC per year of age`)==sign(`log2FC after training`),
                                "Aging",
                                "Rejunevating"),
           observed = -log10(`P-value training`),
           expected = -log10(ppoints(n)),
           clower  = -log10(qbeta(p = (1 - ci) / 2,
                                  shape1 = 1:n,
                                  shape2 = n:1)),
           cupper  = -log10(qbeta(p = (1 + ci) / 2,
                                  shape1 = 1:n,
                                  shape2 = n:1)))

tiff('QQplot_training_mRNA.tiff',
     width =150,
     height = 150,
     units = 'mm',
     res=300)
ggplot(toplot) +
    geom_ribbon(
        mapping = aes(x = expected,
                      ymin = clower,
                      ymax = cupper),
        alpha = 0.1
    ) +
    geom_point(aes(expected,
                   observed,
                   col = Consistency),
               #shape = 2,
               size = 2) +
    scale_color_manual(values = c("Black","Green"))+
    geom_abline(intercept = 0,
                slope = 1,
                alpha = 0.5) +
    geom_label_repel(data = toplot %>% dplyr::filter(`FDR training` < 0.005),
                     aes(expected,
                         observed,
                         label = Gene),
                     #size = 2,
                     #box.padding = unit(0.45, "lines"),
                     #point.padding = unit(0.45, "lines"),
                     max.overlaps = 20
    )+
    ylim(0,25)+
    xlab(log10Pe) +
    ylab(log10Po)+
    theme_bw(base_size = 24) +
    theme(axis.ticks = element_line(size = 0.5),
          panel.grid = element_blank(),
          legend.position = "none")
dev.off()

#Load MetaMEx results for inactivity
library(readxl)
disuse <- read_excel('C:/Users/e5103562/OneDrive - Victoria University/TWAS meta analysis of inactivity/MetaMEx genes summary.xlsx',
                 na = "NA")%>%
    dplyr::filter(`Gene symbol` %in% DEGs$Gene)%>%
    dplyr::select(`Gene symbol`,
           estimate_IN,
           SE_IN,
           pval_IN)%>%
    dplyr::rename(Gene = `Gene symbol`,
           `log2FC after disuse` = estimate_IN,
           `SE disuse` = SE_IN,
           `P-value disuse` = pval_IN) %>%
    mutate(`FDR disuse` = p.adjust(`P-value disuse`))

#disuse <- read_tsv("C:/Users/e5103562/OneDrive - Victoria University/TWAS meta analysis of inactivity/Full results muscle disuse.txt") %>%
#    dplyr::select(Gene.Name,
#                  LFC.mean,
#                  Stouffer.raw)%>%
#    dplyr::rename(Gene = `Gene.Name`,
#                  `log2FC after disuse`=LFC.mean,
#                  `P-value disuse`=Stouffer.raw)%>%
#    dplyr::filter(Gene %in% DEGs$Gene)%>%
#    mutate(`FDR disuse` = p.adjust(`P-value disuse`))
    
DEGs <- left_join(DEGs,
                  disuse)

#QQplot to show inflation of p-values
library(ggrepel)
ci = 0.95
log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

toplot <- DEGs %>%
    dplyr::filter(!is.na(`P-value disuse`))%>%
    arrange(`P-value disuse`)
n <- nrow(toplot)

toplot <- toplot %>%
    mutate(Consistency = ifelse(sign(`log2FC per year of age`)==sign(`log2FC after disuse`),
                                "Aging",
                                "Rejunevating"),
           observed = -log10(`P-value disuse`),
           expected = -log10(ppoints(n)),
           clower  = -log10(qbeta(p = (1 - ci) / 2,
                                  shape1 = 1:n,
                                  shape2 = n:1)),
           cupper  = -log10(qbeta(p = (1 + ci) / 2,
                                  shape1 = 1:n,
                                  shape2 = n:1)))

tiff('QQplot_disuse_mRNA.tiff',
     width =150,
     height = 150,
     units = 'mm',
     res=300)
ggplot(toplot) +
    geom_ribbon(
        mapping = aes(x = expected,
                      ymin = clower,
                      ymax = cupper),
        alpha = 0.1
    ) +
    geom_point(aes(expected,
                   observed,
                   col = Consistency),
               #shape = 2,
               size = 2) +
    scale_color_manual(values = c("Black","Green"))+
    geom_abline(intercept = 0,
                slope = 1,
                alpha = 0.5) +
    geom_label_repel(data = toplot %>% dplyr::filter(`FDR disuse` < 0.005),
                     aes(expected,
                         observed,
                         label = Gene),
                     #size = 2,
                     #box.padding = unit(0.45, "lines"),
                     #point.padding = unit(0.45, "lines"),
                     max.overlaps = 20
    )+
    xlab(log10Pe) +
    ylab(log10Po)+
    ylim(0,25)+
    theme_bw(base_size = 24) +
    theme(axis.ticks = element_line(size = 0.5),
          panel.grid = element_blank(),
          legend.position = "none")
dev.off()

DEGs <- DEGs %>%
    relocate(Gene)

#Write results
write_tsv(x = DEGs,
          na = "",
          file = "mRNA.txt")













#######################################
#Integration age
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
mRNA <- read_tsv('mRNA.txt')%>%
    dplyr::select(Gene,
           `log2FC per year of age`)
DNAm <- read_tsv('DNAm.txt') %>%
    dplyr::rename(Gene = `Annotated gene(s)`)%>%
    separate_rows(Gene,
                  sep=";")%>%
    dplyr::select(Gene,
           `% DNAm change per year of age`)

merged <- left_join(DNAm,
                    mRNA)%>%
    drop_na()

#Venn diagram
library("ggVennDiagram")
x <- list(DMGs = DNAm$Gene,
          DEGs = mRNA$Gene)

tiff('Venn diagram.tiff',
     width =100,
     height = 100,
     units = 'mm',
     res=300)
ggVennDiagram(x,
              label_alpha = 0,
              color = 1,
              lwd = 0.7,
              label = "count")+ 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
dev.off()

library(ggrepel)
tiff('Effect size DNAm vs mRNA for common genes.tiff',
     width =100,
     height = 80,
     units = 'mm',
     res=300)
ggplot(data = merged,
       aes(x = `% DNAm change per year of age`,
           y = `log2FC per year of age`,
           label = Gene)) +
    geom_point()+
    geom_label_repel(data = merged,
                     size = 2,
                     box.padding = unit(0.45, "lines"),
                     point.padding = unit(0.45, "lines"),
                     max.overlaps = 20)+
    geom_hline(yintercept = 0,
               lty="dashed")+
    geom_vline(xintercept = 0,
               lty="dashed")+
    theme_classic()
dev.off()

#Barplot
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
#Load RE meta-analysis results from metafor for age
age <- read_delim('C:/Users/e5103562/OneDrive - Victoria University/TWAS meta analysis of age/metafor_muscle_age_RE.txt')
age <- age %>%
  dplyr::rename(`log2FC per year of age` = ES_RE,
                SE = SE_RE,
                `P-value age` = pval_RE,
                `FDR age` = FDR)%>%
  arrange(`P-value age`)

DNAm <- read_tsv('C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures/DNAm.txt') %>%
  pull(`Annotated gene(s)`)
DNAm <- unique(unlist(strsplit(DNAm,
                               ";")))

#Calculate % of overlap for the different categories
n_DEG_DMG = length(intersect(age %>% dplyr::filter(`FDR age` < 0.005) %>% pull(Gene),
                      DNAm))
n_nonDEG_DMG = length(intersect(age %>% dplyr::filter(`FDR age` > 0.005) %>% pull(Gene),
                             DNAm))
n_DEG_nonDMG = length(setdiff(age %>% dplyr::filter(`FDR age` < 0.005) %>% pull(Gene),
                             DNAm))
n_nonDEG_nonDMG = length(setdiff(age %>% dplyr::filter(`FDR age` > 0.005) %>% pull(Gene),
                                DNAm))

data_summary <- tibble('Category' = c("DEGs",
                                      "non-DEGs"),
                       '% genes that are also\ndifferentially methylated' = 100*c(n_DEG_DMG/(n_DEG_DMG+n_DEG_nonDMG),
                               n_nonDEG_DMG/(n_nonDEG_DMG+n_nonDEG_nonDMG)))

tiff('Barplot DEG vs non-DEG.tiff',
     width =50,
     height = 100,
     units = 'mm',
     res=300)
ggplot(data_summary,
       aes(x = Category,
           y = `% genes that are also\ndifferentially methylated`)) + 
  geom_bar(stat = "identity")  +
  theme_bw()
dev.off()

#Same for DNAm
setwd("C:/Users/e5103562/OneDrive - Victoria University/EWAS meta analysis of age/Muscle")
library(tidyverse)
#Load METAL results that have been formatted to annotate each CpG
to_write <- read_delim('METAL_muscle_age.txt')
to_order <- tibble(CpG = to_write$CpG)

#Load RE meta-analysis results from metafor for age
age <- read_delim('metafor_muscle_age_RE.txt')
age <- left_join(to_order,
                 age)
to_write <- to_write %>%
  mutate(`Effect size` = age$ES_RE,
         SE = age$SE_RE,
         `P-value` = age$pval_RE,
         FDR = age$FDR)%>%
  dplyr::rename(`% DNAm change per year of age` = `Effect size`,
                `P-value age` = `P-value`,
                `FDR age` = FDR)%>%
  arrange(`P-value age`)

library(tidyverse)
mRNA <- read_tsv('C:/Users/e5103562/OneDrive - Victoria University/Tables & figures/mRNA.txt')$Gene

#Calculate % of overlap for the different categories
DMGs <- unique(unlist(strsplit(to_write %>% dplyr::filter(`FDR age` < 0.005) %>% pull(`Annotated gene(s)`),
                               ";")))
nonDMGs <- setdiff(unique(unlist(strsplit(to_write %>% dplyr::filter(`FDR age` > 0.005) %>% pull(`Annotated gene(s)`),
                               ";"))),
                   DMGs)
n_DMG_DEG = length(intersect(DMGs,
                             mRNA))
n_nonDMG_DEG = length(intersect(nonDMGs,
                                mRNA))
n_DMG_nonDEG = length(setdiff(DMGs,
                              mRNA))
n_nonDMG_nonDEG = length(setdiff(nonDMGs,
                                 mRNA))

data_summary <- tibble('Category' = c("DMGs",
                                      "non-DMGs"),
                       '% genes that are also\ndifferentially expressed' = 100*c(n_DMG_DEG/(n_DMG_DEG+n_DMG_nonDEG),
                                                                                 n_nonDMG_DEG/(n_nonDMG_DEG+n_nonDMG_nonDEG)))
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
tiff('Barplot DMG vs non-DMG.tiff',
     width =50,
     height = 100,
     units = 'mm',
     res=300)
ggplot(data_summary,
       aes(x = Category,
           y = `% genes that are also\ndifferentially expressed`)) + 
  geom_bar(stat = "identity")  +
  theme_bw()
dev.off()

dat <- data.frame(
  "DEG_no" = c(1784, 14543),
  "DEG_yes" = c(63, 267),
  row.names = c("Also_DMG", "Not_DMG"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("Non-DEG", "DEG")

mosaicplot(dat,
           main = "Mosaic plot",
           color = TRUE
)
#######################################
#Helper functions
library(GGally)
library(nycflights13)
library(tidyverse)
library(viridis)
#Function to do a correlation plot for the bottom half
my_fn_lower <- function(data,
                        mapping,
                        method="spearman",
                        use="pairwise"){
    
    # grab data
    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y)
    
    # calculate correlation
    corr <- cor(x, y, method=method, use=use)
    
    print(cor.test(x, y, method=method, use=use))
    
    col <- ifelse(corr < 0,
                  "blue",
                  "red")
    
    p <- ggplot(data = data,
                mapping = mapping) +
        geom_point(size=0.5)+
        geom_smooth(method="lm",
                    se = FALSE,
                    col = col,
                    lwd = 0.5)+
        theme_grey(base_size = 8)+
        geom_hline(yintercept = 0,
                   lty="dashed")+
        geom_vline(xintercept = 0,
                   lty="dashed")
    p
}

my_fn_upper <- function(data,
                        mapping,
                        method="spearman",
                        use="pairwise", ...)
{
    
    # grab data
    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y)
    
    # calculate correlation
    corr <- cor(x, y, method=method, use=use)
    
    # calculate colour based on correlation value
    # Here I have set a correlation of minus one to blue,
    # zero to white, and one to red
    # Change this to suit: possibly extend to add as an argument of `my_fn`
    colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
    fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
    
    ggally_cor(data = data,
               stars = FALSE,
               method = method,
               use = use,
               mapping = mapping,
               col="black") +
        theme_void() +
        theme(panel.background = element_rect(fill=fill))
}

#Heatmap functions
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


#################################################
#DNAm
#Comparison effect sizes for age, VO2max and AT
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
DNAm <- read_tsv('DNAm.txt')

#Graph correlation plot for AT
cleanname = function(x,
           lab="\n")
  {
  sapply(x, function(c)
    { 
    str_replace(string = x,
          pattern = "(?<=change)[:space:]",
          replacement = lab)
    }
  )
  }
DNAmtest <- DNAm
colnames(DNAmtest) = cleanname(colnames(DNAmtest))

tiff('Effect size of age vs effect size of VO2max and training_DNAm.tiff',
   width =100,
   height = 100,
   units = 'mm',
   res=300)
ggpairs(
  DNAmtest,
 columns = grep("% DNAm",colnames(DNAm)),
 switch = "both",
 lower=list(continuous=my_fn_lower),
 diag=list(continuous="blank"),
 upper=list(continuous=my_fn_upper))+
  theme(strip.text = element_text(size = 8))
dev.off()

#Heatmap to contrast effect sizes
dat_heat <- DNAm %>%
 arrange(`% DNAm change per year of age`)%>%
 dplyr::rename(`Age`=`% DNAm change per year of age`,
        `VO2max` = `% DNAm change per unit of VO2max`,
        `training`=`% DNAm change after training`)%>%
 dplyr::select(`Age`,
        `VO2max`,
        `training`) %>%
 scale(center=FALSE)%>%
 as_tibble()

myBreaks <- c(seq(min(dat_heat,na.rm=T), 0, length.out=ceiling(100/2) + 1),
       seq(max(dat_heat,na.rm=T)/100, max(dat_heat,na.rm=T), length.out=floor(100/2)))

tiff("Heatmap DNAm.tiff",
   width =50,
   height = 100,
   units = 'mm',
   res=300)
pheatmap(dat_heat,
     color = colorRampPalette(rev(brewer.pal(n = 7, name =
                          "RdBu")))(100),
     #cutree_rows = 2,
     cluster_rows = F,
     cutree_cols = 2,
     breaks = myBreaks,
     cluster_cols = TRUE,
     show_rownames = FALSE)
dev.off()

########################################################
#mRNA
#Comparison effect sizes for age, VO2max and AT
setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
library(tidyverse)
mRNA <- read_tsv('mRNA.txt')

#Graph correlation plot for AT
cleanname = function(x,
                     lab="\n")
{
    sapply(x, function(c)
    { 
        str_replace(string = x,
                    pattern = "(?<=FC)[:space:]",
                    replacement = lab)
    }
    )
}
mRNAtest <- mRNA
colnames(mRNAtest) = cleanname(colnames(mRNAtest))

tiff('Effect size of age vs effect size of VO2max, training and disuse_mRNA.tiff',
     width =120,
     height = 120,
     units = 'mm',
     res=300)
ggpairs(
    mRNAtest,
    columns = grep("log2FC",colnames(mRNA)),
    switch = "both",
    lower=list(continuous=my_fn_lower),
    diag=list(continuous="blank"),
    upper=list(continuous=my_fn_upper))+
    theme(strip.text = element_text(size = 8))
dev.off()

#Heatmap to contrast effect sizes
dat_heat <- mRNA %>%
    arrange(`log2FC per year of age`)%>%
    dplyr::rename(`Age`=`log2FC per year of age`,
                  `VO2max` = `log2FC per unit of VO2max`,
                  `training`=`log2FC after training`,
                  Disuse = `log2FC after disuse`)%>%
    dplyr::select(`Age`,
                  `VO2max`,
                  `training`,
                  Disuse) %>%
    scale(center=FALSE)%>%
    as_tibble()

myBreaks <- c(seq(min(dat_heat,na.rm=T), 0, length.out=ceiling(100/2) + 1),
              seq(max(dat_heat,na.rm=T)/100, max(dat_heat,na.rm=T), length.out=floor(100/2)))

tiff("Heatmap mRNA.tiff",
     width =80,
     height = 110,
     units = 'mm',
     res=300)
pheatmap(dat_heat,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                     "RdBu")))(100),
         #cutree_rows = 2,
         cluster_rows = F,
         cutree_cols = 2,
         breaks = myBreaks,
         cluster_cols = TRUE,
         show_rownames = FALSE)
dev.off()


######################################
#DNAm
#Volcano plot showing the number of years of rejuvenation for each unit of VO2max or after AT
vol <- DNAm %>%
    mutate(`Rejuvenation effect of VO2max (in years of age)`=`% DNAm change per unit of VO2max`/`% DNAm change per year of age`)%>%
    mutate(Significance = ifelse(`FDR VO2max`<0.005,
                                 "Sig",
                                 "Not Sig"))
sigdir <- vol$Significance
sigdir[vol$Significance=="Sig"&
           vol$`Rejuvenation effect of VO2max (in years of age)`<0]="Rejuvenating"
#sigdir[vol$Significance=="Sig"&
#     vol$`Years of age gained/lost per unit of VO2max`>0]="Aging"
vol$Direction <- sigdir

vol <- vol %>%
    filter(!is.na(`P-value VO2max`))
tiff('Volcano plot_contrast age and VO2max.tiff',
     width =105,
     height = 70,
     units = 'mm',
     res=300)
ggplot() +
    geom_point(data = vol,
               aes(`Rejuvenation effect of VO2max (in years of age)`,
                   -log10(`P-value VO2max`)),
               size=0.55)+
    #scale_color_manual(values=c("black","red"))+
    labs(y="-log10(p-value for VO2max)")+
    ylim(c(0,11))+
    theme_classic()+
    theme(legend.position = "none")+
    geom_point(data = vol %>% filter(Significance=="Sig"),
               aes(`Rejuvenation effect of VO2max (in years of age)`,
                   -log10(`P-value VO2max`)),
               col = "red",
               size=1.5)
dev.off()

#Same volcano plot for training
vol <- DNAm %>%
    mutate(`Rejuvenation effect of training (in years of age)`=`% DNAm change after training`/`% DNAm change per year of age`)

vol <- vol %>%
    dplyr::filter(!is.na(`P-value training`))

tiff('Volcano plot_contrast age and training.tiff',
     width =105,
     height = 70,
     units = 'mm',
     res=300)
ggplot() +
    geom_point(data = vol,
               aes(`Rejuvenation effect of training (in years of age)`,
                   -log10(`P-value training`)),
               size=0.55)+
    labs(y="-log10(p-value for training)")+
    ylim(c(0,11))+
    theme_classic()+
    theme(legend.position = "none")
#geom_point(data = vol %>% filter(Significance=="Sig"),
#           aes(`Years of age gained or lost after training`,
#               -log10(`P-value training`)),
#           col = "red",
#           size=1)
dev.off()



#######################
#mRNA
#Volcano plot showing the number of years of rejuvenation for each unit of VO2max or after AT
vol <- mRNA %>%
    mutate(`Rejuvenation effect of VO2max (in years of age)`=`log2FC per unit of VO2max`/`log2FC per year of age`)%>%
    mutate(Significance = ifelse(`FDR VO2max`<0.005,
                                 "Sig",
                                 "Not Sig"))
sigdir <- vol$Significance
sigdir[vol$Significance=="Sig"&
           vol$`Rejuvenation effect of VO2max (in years of age)`<0]="Rejuvenating"
#sigdir[vol$Significance=="Sig"&
#     vol$`Years of age gained/lost per unit of VO2max`>0]="Aging"
vol$Direction <- sigdir

vol <- vol %>%
    dplyr::filter(!is.na(`P-value VO2max`))

tiff('Volcano plot_contrast age and VO2max_mRNA.tiff',
     width =105,
     height = 70,
     units = 'mm',
     res=300)
ggplot() +
    geom_point(data = vol,
               aes(`Rejuvenation effect of VO2max (in years of age)`,
                   -log10(`P-value VO2max`)),
               size=0.55)+
    #scale_color_viridis()+
    #scale_color_manual(values=c("black","red"))+
    labs(y="-log10(p-value for VO2max)")+
    ylim(c(0,25))+
    theme_classic()+
    theme(legend.position = "none")+
    geom_point(data = vol %>% dplyr::filter(Significance=="Sig"),
               aes(`Rejuvenation effect of VO2max (in years of age)`,
                   -log10(`P-value VO2max`)),
               col = "red",
               size=1.5)
dev.off()

#Same volcano plot for training
vol <- mRNA %>%
    mutate(`Rejuvenation effect of training (in years of age)`=`log2FC after training`/`log2FC per year of age`)%>%
    mutate(Significance = ifelse(`FDR training`<0.005,
                                 "Sig",
                                 "Not Sig"))
sigdir <- vol$Significance
sigdir[vol$Significance=="Sig"&
           vol$`Rejuvenation effect of training (in years of age)`<0]="Rejuvenating"
#sigdir[vol$Significance=="Sig"&
#     vol$`Years of age gained/lost per unit of VO2max`>0]="Aging"
vol$Direction <- sigdir

vol <- vol %>%
    dplyr::filter(!is.na(`P-value training`))
nrow(vol %>% dplyr::filter(Significance=="Sig"))
mean(vol %>% dplyr::filter(Significance=="Sig") %>% pull(`Rejuvenation effect of training (in years of age)`))

tiff('Volcano plot_contrast age and training_mRNA.tiff',
     width =105,
     height = 70,
     units = 'mm',
     res=300)
ggplot() +
    geom_point(data = vol,
               aes(`Rejuvenation effect of training (in years of age)`,
                   -log10(`P-value training`)),
               size=0.55)+
    labs(y="-log10(p-value for training)")+
    ylim(c(0,25))+
    theme_classic()+
    theme(legend.position = "none")+
geom_point(data = vol %>% dplyr::filter(Significance=="Sig"),
           aes(`Rejuvenation effect of training (in years of age)`,
               -log10(`P-value training`)),
           col = "red",
           size=1)
dev.off()


#Same volcano plot for disuse
vol <- mRNA %>%
    mutate(`Rejuvenation effect of disuse (in years of age)`=`log2FC after disuse`/`log2FC per year of age`)%>%
    mutate(Significance = ifelse(`FDR disuse`<0.005,
                                 "Sig",
                                 "Not Sig"))
sigdir <- vol$Significance
sigdir[vol$Significance=="Sig"&
           vol$`Rejuvenation effect of disuse (in years of age)`<0]="Rejuvenating"
vol$Direction <- sigdir

vol <- vol %>%
    dplyr::filter(!is.na(`P-value disuse`))
nrow(vol %>% dplyr::filter(Significance=="Sig"))
mean(vol %>% dplyr::filter(Significance=="Sig") %>% pull(`Rejuvenation effect of disuse (in years of age)`))

tiff('Volcano plot_contrast age and disuse_mRNA.tiff',
     width =105,
     height = 70,
     units = 'mm',
     res=300)
ggplot() +
    geom_point(data = vol,
               aes(`Rejuvenation effect of disuse (in years of age)`,
                   -log10(`P-value disuse`)),
               size=0.55)+
    #scale_color_manual(values=c("black","red"))+
    labs(y="-log10(p-value for disuse)")+
    #ylim(c(0,10))+
    theme_classic()+
    theme(legend.position = "none")+
geom_point(data = vol %>%  dplyr::filter(Significance=="Sig"),
           aes(`Rejuvenation effect of disuse (in years of age)`,
               -log10(`P-value disuse`)),
           col = "red",
           size=1.5)
dev.off()








#Forest plot for top hypo and hyper DMP
L <- read_rds("C:/Users/e5103562/OneDrive - Victoria University/EWAS meta analysis of VO2max/METAL_muscle_VO2max.rds")
library(readxl)
dataset_summary <- read_excel('C:/Users/e5103562/OneDrive - Victoria University/DNAm datasets.xlsx',
               na = c("","NA"))%>%
  filter(`VO2max data`=="Yes")%>%
  dplyr::rename(prop_males = `% Male`,
         Dataset = `Dataset ID`)
library(metafor)
meta_RE <- function(tib = NULL)
{
  tib <- tib %>%
    filter(Type != "Meta-analysis")
  try({
    model <- rma(yi = ES,
           sei = SE,
           data=tib,
           method = "EB",
           control=list(stepadj=0.5,
                 maxiter=10000,
                 threshold=1e-8))
    return(model)
  })
}

top_hypo <- DNAm %>%
  arrange(`P-value VO2max`)%>%
  filter(`% DNAm change per unit of VO2max` < 0) %>%
  dplyr::slice(1)%>%
  pull(CpG)
list <- L[[top_hypo]]%>%
  mutate(Dataset = str_replace(Dataset,
                 pattern = "_VO2max",
                 replacement = ""))
tib <- left_join(list,
         dataset_summary)%>%
  arrange(-`n after preprocessing`)

top_hyper <- DNAm %>%
  arrange(`P-value VO2max`)%>%
  filter(`% DNAm change per unit of VO2max` > 0) %>%
  dplyr::slice(1)%>%
  pull(CpG)
list <- L[[top_hyper]]%>%
  mutate(Dataset = str_replace(Dataset,
                 pattern = "_VO2max",
                 replacement = ""))
tib <- left_join(list,
         dataset_summary)%>%
  arrange(-`n after preprocessing`)

tiff('Forest plot_top_hyper_VO2max.tiff',
   width =130,
   height = 140,
   units = 'mm',
   res=300)
forest(meta_RE(tib),
    header="Dataset ID",
    slab = tib %>%
      filter(Dataset != "Meta-analysis")%>%
      pull(Dataset),
    xlab = "% DNAm change per unit of VO2max",
    digits=c(2L,4L))
dev.off()




#Supplementary Table 2
memory.limit(size=100000)
library(tidyverse)
setwd("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Epigenetic aging and ")
library(readxl)
DEG_RNA <- read_excel("Differentially expressed genes with age from Su et al 2015.xlsx",
           sheet = "DEG with age")%>%
 dplyr::select(`ENSEMBL ID`,
     `Gene symbol`,
     `Beta for min(p)`,
     `min(p) (age)`)%>%
 dplyr::rename(`ENSEMBL gene ID`=`ENSEMBL ID`,
        `Gene name`=`Gene symbol`,
        `Effect size age`=`Beta for min(p)`,
        `P-value age`=`min(p) (age)`)

#Add meta-analysis of VO2max
meta <- read_tsv("Meta-analysis results mRNA expression VO2max.txt")%>%
 dplyr::select(-c(`ENSEMBL gene ID`,
      `Effect size age`,
      `P-value age`))%>%
  mutate(`Effect size VO2max` = `Effect size VO2max`/100,
      `Standard error VO2max` = `Standard error VO2max`/100)

#Add ExTraMeta and Deane et al.
MetaMEx <- read_excel("MetaMEx genes summary.xlsx")%>%
  dplyr::rename(#`log2FC after AA`=estimate_AA,
         #`Standard error AA`=SE_AA,
         #`P-value AA`=pval_AA,
         `log2FC after AT`=estimate_TA,
         `Standard error AT`=SE_TA,
         `P-value AT`=pval_TA,
         `log2FC after RT`=estimate_TR,
         `Standard error RT`= SE_TR,
         `P-value RT`= pval_TR,
         `Gene name`= `Gene symbol`)%>%
  filter(`Gene name` %in% DEG_RNA$`Gene name`)%>%
  mutate(#`FDR AA`=p.adjust(`P-value AA`),
      `FDR AT`=p.adjust(`P-value AT`),
      `FDR RT`=p.adjust(`P-value RT`))%>%
  select(`Gene name`,
      #`log2FC after AA`,
      #`Standard error AA`,
      #`P-value AA`,
      #`FDR AA`,
      `log2FC after AT`,
      `Standard error AT`,
      `P-value AT`,
      `FDR AT`,
      `log2FC after RT`,
      `Standard error RT`,
      `P-value RT`,
      `FDR RT`)

ExTraMeta <- read_tsv("Extrameta all genes results.txt")%>%
  select(-I2)%>%
  dplyr::rename(`log2FC after AT`=ES,
      `Standard error AT`=SE,
      `P-value AT`=p,
      Moderators = mod_type,
      `Gene name`=gene)%>%
  filter((`Gene name` %in% DEG_RNA$`Gene name`)&
        `P-value AT`!=0)%>%
  mutate(`FDR AT`=p.adjust(`P-value AT`))%>%
  select(`Gene name`,
      `log2FC after AT`,
      `Standard error AT`,
      `P-value AT`,
      `FDR AT`,
      Moderators)%>%
  mutate(Moderators = str_replace(Moderators,
                  "base_model",
                  ""))%>%
  mutate(Moderators = str_replace(Moderators,
                  "time",
                  "AT length"))%>%
  mutate(Moderators = str_replace(Moderators,
                  "avg_age",
                  "Age"))%>%
  mutate(Moderators = str_replace(Moderators,
                  "AT",
                  "AT type"))%>%
  mutate(Moderators = str_replace(Moderators,
                  "prop_males",
                  "Sex"))

DIS <- read_tsv("Full results muscle disuse.txt") %>%
 dplyr::select(Gene.Name,
     LFC.mean,
     Stouffer.raw)%>%
 dplyr::rename(`Gene name` = `Gene.Name`,
        `Effect size muscle disuse`=LFC.mean,
        `P-value muscle disuse`=Stouffer.raw)%>%
 mutate(`FDR muscle disuse` = p.adjust(`P-value muscle disuse`))

#Merge all
library(plyr)
merged <- join_all(list(DEG_RNA,
            meta,
            MetaMEx,
            DIS),
          by='Gene name',
          type='left')

setwd("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Epigenetic aging and ")
write.table(merged,
      file="Supplementary Table 2.txt",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      sep="\t")

#Supplementary Table 3
library(tidyverse)
setwd("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Epigenetic aging and ")
library(readxl)
DEG_prot <- read_excel("Proteins that change with age in human skeletal muscle.xlsx",
            sheet = "Significant Proteins",
            skip = 2)%>%
 select(`Gene`,`Acc No.`,`Protein ID`,`ProteinDescription`,`AgeBeta`,`Pvalue`,`BH Adjusted Pvalue`)%>%
 dplyr::rename(`Uniprot ID` = `Acc No.`,
        `Uniprot Entry Name` = `Protein ID`,
        `Gene description` = `ProteinDescription`,
        `Effect size age` = `AgeBeta`,
        `P-value age` = `Pvalue`,
        `FDR age`=`BH Adjusted Pvalue`)%>%
 filter(`FDR age`<0.05)

#Read meta-analysis of age
#DEG_prot <- read_tsv("Meta analysis proteins associated with age.txt")%>%
# dplyr::rename(`Gene` = `Gene name`)%>%
# filter(`FDR age`<0.1)

#Physical activity from Ubaida Mohien Frontiers 2021
prot_PA <- read_excel("Proteins associated with higher physical activity in Ubaida-Mohien et al 2019 Frontiers in Physiology.xlsx",
             sheet="Table S2_AllProteins",
           skip = 1)%>%
 select(`GenePrimary`,`Beta (Physical Activity)`,`P-value`)%>%
 dplyr::rename(`Gene` = `GenePrimary`,
        `Effect size` = `Beta (Physical Activity)`)
colnames(prot_PA)[-1] <- paste0(colnames(prot_PA)[-1]," PA")

#Results for VO2max
prot_VO2max <- read_tsv("GeneSMART proteins associated with VO2max.txt")%>%
 select(`Gene name`,`log2FC`,`SE`,`P.Value`)%>%
 dplyr::rename(`Gene` = `Gene name`,
        `Effect size` = `log2FC`,
        `Standard error` = `SE`,
        `P-value` = `P.Value`)
colnames(prot_VO2max)[-1] <- paste0(colnames(prot_VO2max)[-1]," VO2max")

#HIIT AT
prot_AT <- read_tsv("GeneSMART proteins Timepoint 4WP.txt")%>%
 select(`Gene name`,`log2FC`,`SE`,`P.Value`)%>%
 dplyr::rename(`Gene` = `Gene name`,
        `Effect size` = `log2FC`,
        `Standard error` = `SE`,
        `P-value` = `P.Value`)
colnames(prot_AT)[-1] <- paste0(colnames(prot_AT)[-1]," AT")

#Merge all
library(plyr)
merged <- join_all(list(DEG_prot,
            prot_PA,
            prot_VO2max,
            prot_AT),
          by='Gene',
          type='left')%>%
 mutate(`FDR PA` = p.adjust(`P-value PA`),
     `FDR VO2max` = p.adjust(`P-value VO2max`),
     `FDR AT` = p.adjust(`P-value AT`))

#Write results for paper
setwd("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Epigenetic aging and ")
write.table(merged,
      file="Supplementary Table 3.txt",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      sep="\t")

##############################Figure 1: Venn diagram
setwd("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Epigenetic aging and ")
library(tidyverse)
DNAm <- read_tsv('Supplementary Table 1.txt')
mRNA <- read_tsv('Supplementary Table 2.txt')
prot <- read_tsv('Supplementary Table 3.txt')

DMG = unlist(strsplit(DNAm%>%
            pull(`Annotated gene(s)`),split=";"))
DMG = unique(na.omit(DMG))
DEG_mRNA = mRNA %>%
 pull(`Gene name`)
DEG_prot = prot %>%
 pull(`Gene`)

L = list(
 Protein = DEG_prot,
 mRNA = DEG_mRNA,
 DMG = DMG)
genes_common_toall <- Reduce(intersect,L)

library(VennDiagram)
venn.diagram(L,
       category.names = c("",
                "",
                ""),
       filename = 'Intersection of age-related CpGs, mRNAs and proteins.tiff',
       height = 2.5 ,
       width = 2.5 ,
       units = "in",
       resolution = 600,
       lwd = 1,
       col=c("#ED7D31", '#70AD47', '#FFC000'),
       fill = c(alpha("#ED7D31",0.3), alpha('#70AD47',0.3), alpha('#FFC000',0.3)),
       cex = 0.7,
       fontfamily = "sans",
       cat.cex = 0.5,
       cat.default.pos = "outer",
       cat.pos = c(-27, 27, 180),
       cat.dist = c(0.15, 0.13, 0.05),
       cat.fontfamily = "sans",
       cat.col = c("#ED7D31", '#70AD47', '#FFC000'),
       margin = 0.2
       #rotation = 1
)


#Part B of the figure (mRNA)
setwd("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Epigenetic aging and ")
library(tidyverse)
mRNA <- read_tsv('Supplementary Table 2.txt')%>%
 dplyr::rename(`log2FC per year`=`Effect size age`,
        `log2FC per mL/min/kg`=`Effect size VO2max`,
        `log2FC after disuse`=`Effect size muscle disuse`)%>%
 mutate(VO2maxscore = -log10(`FDR VO2max`)*sign(`log2FC per mL/min/kg`),
     #AAscore = -log10(`FDR AA`)*sign(`log2FC after AA`),
    ATscore = -log10(`FDR AT`)*sign(`log2FC after AT`),
     RTscore = -log10(`FDR RT`)*sign(`log2FC after RT`),
     inactivityscore = -log10(`FDR muscle disuse`)*sign(`log2FC after disuse`))
 #mutate(VO2maxscore = -log10(`FDR VO2max`)*sign(`log2FC per mL/min/kg`),
 #   ATscore = -log10(`FDR AT`)*sign(`log2FC after AT`),
 #    inactivityscore = -log10(`FDR muscle disuse`)*sign(`log2FC after disuse`))
 
#Add column to color points differently according to FDR with VO2max or 
mRNA <- mRNA %>%
 mutate(color = ifelse(`FDR VO2max`<0.05|`FDR AT`<0.05|`FDR muscle disuse`<0.05,
            TRUE,
            FALSE))

mRNA_full <- mRNA %>%
 filter(!is.na(`log2FC after AT`))
nrow(mRNA_full %>%
    filter(sign(`log2FC per year`)!=sign(`log2FC after AT`)))/nrow(mRNA_full)

mRNA_sig <- mRNA %>%
  filter(`FDR AT`<0.05)
mean(abs(mRNA_sig$`log2FC after AT`/mRNA_sig$`log2FC per year`),
   na.rm=TRUE)

#Add topconfects
#library(topconfects)
#confects_VO2max <- normal_confects(effect = mRNA$`log2FC per mL/min/kg`,
#                  se = mRNA$`Standard error VO2max`,
#                  fdr = 0.005)$table
#confects_VO2max <- confects_VO2max %>%
# arrange(index)%>%
# pull(confect)
#confects_ <- normal_confects(effect = mRNA$`log2FC after AT`,
#                  se = mRNA$`Standard error AT`,
#                  fdr = 0.005)$table
#confects_ <- confects_ %>%
# arrange(index)%>%
# pull(confect)

#mRNA <- mRNA %>%
# mutate(confects_VO2max = confects_VO2max,
#     confects_ = confects_)

tiff('Effect size of age vs effect size of VO2max AT and disuse_mRNA.tiff',
   width =6.5,
   height = 6,
   units = 'in',
   res=600)
ggpairs(
 mRNA,
    columns = grep("log2FC",colnames(mRNA)),
    switch = "both",
    lower=list(continuous=my_fn_lower),
    diag=list(continuous="blank"),
    upper=list(continuous=my_fn_upper)
 )
dev.off()

library(plyr)
#Age vs AT
tiff('Effect size of age vs effect size of VO2max_mRNA_bubbles.tiff',
   width =3,
   height = 3,
   units = 'in',
   res=600)
ggplot(data = #subset(
     mRNA%>%
     #filter(`Gene name`!="MICALL1")%>%
     arrange(abs(VO2maxscore)),
    #,abs(VO2max_score)>7&abs(_score)>3)
    mapping = aes(x = `log2FC per year`,
           y = `log2FC per mL/min/kg`))+
 #geom_point(size=0.5,
 #      col = "grey51",
 #     data = . %>% filter(color ==FALSE))+
 #geom_point(size=0.5,
 #      col = "grey51",
 #      data = . %>% filter(color ==TRUE))+
  #geom_smooth(method=lm)+
  xlim(-0.02,0.02)+
 #ylim(-0.2,0.2)+
 geom_point(mapping = aes(#size = abs(ATscore),
              col = VO2maxscore))+
 geom_smooth(method=lm,
       se = FALSE)+ 
 #stat_density_2d(aes(fill = ..level..),
 #        geom="polygon",
 #        data = . %>% filter(color ==TRUE))+
 #scale_fill_viridis(option = "C")+
  scale_colour_gradient2(high = muted("red"),
              mid = "grey75",
              low = muted("blue"),
              breaks = c(-log10(0.05),
                   log10(0.05)),
              limits = c(-5,5)
  )+
 geom_hline(yintercept = 0,
       lty="dashed")+
 geom_vline(xintercept = 0,
       lty="dashed")+
 theme(legend.position = "none",
    panel.background = element_rect(fill = "white",
                    colour = "white",
                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                    colour = "gray75"), 
    panel.grid.minor = element_line(size = 0.125, linetype = 'solid',
                    colour = "gray75"))
dev.off()

#Heatmap
dat_heat <- mRNA %>%
 #filter(color==TRUE)%>%
 dplyr::select(`log2FC per year`,
        `log2FC per mL/min/kg`,
        `log2FC after AT`,
        `log2FC after RT`,
        `log2FC after disuse`) %>%
 dplyr::rename(`Age`=`log2FC per year`,
        `VO2max` = `log2FC per mL/min/kg`,
        `AT`=`log2FC after AT`,
        `RT`=`log2FC after RT`,
        `Disuse` = `log2FC after disuse`)%>%
 scale(center=FALSE)%>%
 as_tibble()

myBreaks <- c(seq(min(dat_heat,na.rm=T), 0, length.out=ceiling(100/2) + 1),
       seq(max(dat_heat,na.rm=T)/100, max(dat_heat,na.rm=T), length.out=floor(100/2)))

tiff("Heatmap mRNA.tiff",
   width =5,
   height = 7,
   units = 'in',
   res=600)
pheatmap(dat_heat,
     color = colorRampPalette(rev(brewer.pal(n = 7, name =
                          "RdBu")))(100),
     cutree_rows = 2,
     cutree_cols = 2,
     breaks = myBreaks,
     cluster_cols = TRUE,
     show_rownames = FALSE)
dev.off()

#Part C of the figure (protein)
setwd("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Epigenetic aging and ")
library(tidyverse)
prot <- read_tsv('Supplementary Table 3.txt')%>%
 dplyr::rename(`log2FC per year`=`Effect size age`,
        `log2FC per unit of PA` = `Effect size PA`,
        `log2FC per mL/min/kg`=`Effect size VO2max`,
        `log2FC after AT`=`Effect size AT`)
#Add empty column for inactivity
prot <- prot %>%
 mutate(`log2FC after disuse`=rep(NA,nrow(prot)))

#Add column to color points differently according to FDR with VO2max or 
prot <- prot %>%
 mutate(color = ifelse(`FDR VO2max`<0.005|`FDR AT`<0.005,TRUE,FALSE))

tiff('Effect size of age vs effect size of VO2max and AT_prot.tiff',
   width =6.5,
   height = 6,
   units = 'in',
   res=600)
ggpairs(prot,
 #subset(prot,color==TRUE),
 columns = c(5,10,13,19),
 switch = "both",
 #lower=list(continuous=my_fn_lower),
 diag=list(continuous="blank"),
 upper=list(continuous=my_fn_upper)
)
dev.off()

library(plyr)
#Age vs AT
tiff('Effect size of VO2max vs effect size of AT_prot.tiff',
   width =3,
   height = 3,
   units = 'in',
   res=600)
ggplot(data = #subset(
     prot,
    #,abs(VO2max_score)>7&abs(_score)>3)
    mapping = aes(x = `log2FC per mL/min/kg`,
           y = `log2FC after AT`))+
 #xlim(min(DNAm$`% DNAm per year`),
 #  max(DNAm$`% DNAm per year`))+
 #ylim(min(DNAm$`% DNAm after AT`,na.rm=TRUE),
 #   max(DNAm$`% DNAm after AT`,na.rm=TRUE))+
 geom_point(size=0.5,
       col = "grey",
       data = . %>% filter(color ==FALSE))+
 geom_point(size=2,
       col = "red",
       data = . %>% filter(color ==TRUE))+
 #stat_density_2d(aes(fill = ..level..),
 #        geom="polygon",
 #        data = . %>% filter(color ==TRUE))+
 #scale_fill_viridis(option = "C")+
 #scale_fill_gradientn(colours=plasma(100))+
  geom_smooth(method=lm)+ 
 geom_hline(yintercept = 0,
       lty="dashed")+
 geom_vline(xintercept = 0,
       lty="dashed")+
 theme(legend.position = "none")
dev.off()


#Heatmap
dat_heat <- prot %>%
  filter(color==TRUE)%>%
  dplyr::select(`log2FC per year`,
         `log2FC per mL/min/kg`,
         `log2FC after AT`) %>%
  dplyr::rename(`Age`=`log2FC per year`,
         `VO2max` = `log2FC per mL/min/kg`,
         `AT`=`log2FC after AT`)%>%
  scale(center=FALSE)%>%
  as_tibble()

myBreaks <- c(seq(min(dat_heat,na.rm=T), 0, length.out=ceiling(100/2) + 1),
       seq(max(dat_heat,na.rm=T)/100, max(dat_heat,na.rm=T), length.out=floor(100/2)))

tiff("Heatmap prot signif.tiff",
   width =3.5,
   height = 7,
   units = 'in',
   res=600)
pheatmap(dat_heat,
     color = colorRampPalette(rev(brewer.pal(n = 7, name =
                           "RdBu")))(100),
     cutree_rows = 2,
     cutree_cols = 2,
     breaks = myBreaks,
     cluster_cols = TRUE,
     show_rownames = FALSE)
dev.off()

############Forest plot


#########################Try interactive plot
library(plotly)
library(quantmod)

accumulate_by <- function(dat, var) {
  var <- lazyeval::f_eval(var, dat)
  lvls <- plotly:::getLevels(var)
  dats <- lapply(seq_along(lvls), function(x) {
    cbind(dat[var %in% lvls[seq(1, x)], ], frame = lvls[[x]])
  })
  dplyr::bind_rows(dats)
}

setwd("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Epigenetic aging and ")
library(tidyverse)
DNAm <- read_tsv('Supplementary Table 1.txt')%>%
 dplyr::rename(`% DNAm per year`=`Effect size age`,
        `% DNAm per mL/min/kg`=`Effect size VO2max`,
        `% DNAm after AT`=`Effect size AT`)
#Add empty column for inactivity
DNAm <- DNAm %>%
 mutate(`% DNAm after disuse`=rep(NA,nrow(DNAm)))

#Add column to color points differently according to FDR with VO2max or 
DNAm <- DNAm %>%
 mutate(color = ifelse(`FDR VO2max`<0.005|`FDR AT`<0.005
            ,TRUE,FALSE))

DNAm <- DNAm %>%
 accumulate_by(~`FDR VO2max`)

fig <- DNAm %>%
 plot_ly(
  x = ~`% DNAm per year`, 
  y = ~`% DNAm per mL/min/kg`, 
  frame = ~frame,
  type = 'scatter', 
  mode = 'markers', 
  fill = 'color', 
  #fillcolor='rgba(114, 186, 59, 0.5)',
  #line = list(color = 'rgb(114, 186, 59)'),
  text = ~`Annotated gene(s)`, 
  hoverinfo = 'text') %>%
 add_segments(x = -1, xend = 1, y = 0, yend = 0) %>%
 add_segments(x = 0, xend = 0, y = -1, yend = 1)

fig <- fig %>% layout(
  yaxis = list(
    title = "% DNAm per mL/min/kg", 
    #range = c(0,250), 
    zeroline = F
    #tickprefix = "$"
  ),
  xaxis = list(
    title = "% DNAm per year", 
    #range = c(0,30), 
    #zeroline = F, 
    showgrid = F
  )
) 
fig <- fig %>% animation_opts(
  frame = 50, 
  transition = 0, 
  redraw = FALSE
)
fig <- fig %>% animation_slider(
  currentvalue = list(
    prefix = "FDR "
  )
)

fig

#######################
#Integration all contrasts at transcriptome level with mitch

setwd("C:/Users/e5103562/OneDrive - Victoria University/Tables & Figures")
mRNA <- read_tsv('mRNA.txt')%>%
    dplyr::select(Gene,
                  `log2FC per year of age`)
DNAm <- read_tsv('DNAm.txt') %>%
    dplyr::rename(Gene = `Annotated gene(s)`)%>%
    separate_rows(Gene,
                  sep=";")%>%
    dplyr::select(Gene,
                  `% DNAm change per year of age`)

#############
#Doughnut chart
require(moonBook)
require(webr)
doughnut <- GeneSMART_results %>%
 filter(PVALUE_ADJ_CORR_basefit < 0.005) %>%
 #filter(!is.na(EFFECTSIZE_CORR_control)) %>%
 mutate(Consistency = ifelse(sign(`Effect size`)==sign(EFFECTSIZE_CORR_basefit),
               "Same direction",
               "Opposite direction"))
dir <- function(v)
{
 x <- as.numeric(v["Effect size"])
 y <- as.numeric(v["EFFECTSIZE_CORR_basefit"])
 output <- "UP with age &
 DOWN with fitness"
 if(sign(x)==sign(y)&sign(x)==-1)
 {
  output <- "DOWN with age &
 DOWN with fitness"
 }
 else if(sign(x)==sign(y)&sign(x)==1)
 {
  output <- "UP with age &
 UP with fitness"
 }
 else if(sign(x)!=sign(y)&sign(x)==1)
 {
  output <- "DOWN with age &
 UP with fitness"
 }
 return(output)
}
Directions <- apply(doughnut,1,dir)
doughnut <- doughnut %>%
 mutate(Directions = Directions)

tiff('Doughnut pie DNAm age VO2max direction of effect.tiff',
   width =3.5,
   height = 3.5,
   units = 'in',
   res=600)
PieDonut(doughnut,
     aes(Consistency,Directions),
     showPieName = FALSE,
     donutLabelSize = 4,
     labelposition = 1,
     r0=0,
     explode=1)
dev.off()

#Ratio histogram
doughnut <- doughnut %>%
 mutate(ratio = EFFECTSIZE_CORR_basefit/`Effect size`)
gg_color_hue <- function(n) {
 hues = seq(15, 375, length = n + 1)
 hcl(h = hues, l = 65, c = 100)[1:n]
}
mainCol=gg_color_hue(2)

tiff('Distribution of ratio of fitness ES and age ES DNAm.tiff',
   width =6,
   height = 3,
   units = 'in',
   res=600)
ggplot(data = doughnut,
    aes(x=ratio,
      fill=Consistency)) +
 geom_histogram(colour = "black")+
 scale_fill_manual(values=c(mainCol[1],mainCol[2]))+
 labs(x="Ratio between VO2max ES and age ES")+
 theme_bw()
dev.off()



#Doughnut
setwd("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Epigenetic aging and ")
DEG_RNA <- read_tsv("Supplementary Table 3.txt")

DEG_RNA <- DEG_RNA %>%
 filter(`FDR VO2max` < 0.005) %>%
 mutate(Consistency = ifelse(sign(`Effect size age`)==sign(`Effect size VO2max`),
               "Same direction",
               "Opposite direction"))
dir <- function(v)
{
 x <- as.numeric(v["Effect size age"])
 y <- as.numeric(v["Effect size VO2max"])
 output <- "UP with age &
 DOWN with fitness"
 if(sign(x)==sign(y)&sign(x)==-1)
 {
  output <- "DOWN with age &
 DOWN with fitness"
 }
 else if(sign(x)==sign(y)&sign(x)==1)
 {
  output <- "UP with age &
 UP with fitness"
 }
 else if(sign(x)!=sign(y)&sign(x)==1)
 {
  output <- "DOWN with age &
 UP with fitness"
 }
 return(output)
}
Directions <- apply(DEG_RNA,1,dir)
DEG_RNA <- DEG_RNA %>%
 mutate(Directions = Directions)

tiff('Doughnut pie mRNA age VO2max direction of effect.tiff',
   width =3.5,
   height = 3.5,
   units = 'in',
   res=600)
PieDonut(DEG_RNA,
     aes(Consistency,Directions),
     showPieName = FALSE,
     donutLabelSize = 4,
     labelposition = 1,
     r0=0,
     explode=1)

#########################Take intersection of DNAm and mRNA as example############
##############Intersect all significant genes##############
setwd("F:/ISEAL 2012-2013 and postdoc 2017-2019/Projects/Epigenetic aging and ")
library(tidyverse)
DNAm <- read_tsv('Supplementary Table 1.txt')
mRNA <- read_tsv('Supplementary Table 2.txt')
#prot <- read_tsv('Supplementary Table 3.txt')

DMG = unlist(strsplit(DNAm%>%
             filter(`FDR VO2max`<0.05|`FDR AT`<0.05)%>%
             pull(`Annotated gene(s)`),split=";"))
DMG = unique(na.omit(DMG))

DEG_mRNA = mRNA %>%
  filter(`FDR VO2max`<0.05|`FDR AT`<0.05|`FDR muscle disuse`<0.05)%>%
  pull(`Gene name`)

intersect(DMG,DEG_mRNA)
#DEG_prot = prot %>%
#  filter(`FDR PA`<0.05|`FDR VO2max`<0.05|`FDR AT`<0.05)%>%
#  pull(`Gene`)

L <- list(DMG = DMG,
     DEG_mRNA = DEG_mRNA,
     DEG_prot = DEG_prot)
library(ComplexHeatmap)
m1 = make_comb_mat(L)
library(UpSetR)
UpSet(m1,
   comb_order = order(comb_size(m1)))


