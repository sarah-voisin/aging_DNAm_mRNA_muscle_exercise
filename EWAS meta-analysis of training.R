#Set working directory
setwd("C:/Users/e5103562/OneDrive - Victoria University/EWAS meta analysis of training")

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
    mutate_at(vars(EFFECTSIZE_CORR:SE_CORR),
              function(x){x*100})
    
  write.table(f,
              file=paste0('./',folder,'/',dataset,".tbl"),
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE,
              sep="\t")
}

sapply(files[c(2,4)],
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
meta_res <- read_tsv("METAL_muscle_training.TBL")
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

#View results including in at least 5 studies
meta_res_robust <- meta_res_tot %>% filter(HetDf>=4)
padj <- p.adjust(meta_res_robust$`P-value`,method="BH")
meta_res_robust <- meta_res_robust %>%
  mutate(`Adjusted p-value` = padj)
meta_res_robust <- meta_res_robust %>%
  mutate(sig = ifelse(`Adjusted p-value`<0.005,
                      "Sig",
                      "Not Sig"))
nrow(meta_res_robust) #638,213

#####Arrange
sigdir <- meta_res_robust$sig
sigdir[meta_res_robust$sig=="Sig"&
         meta_res_robust$Effect<0]="Hypo"
sigdir[meta_res_robust$sig=="Sig"&
         meta_res_robust$Effect>0]="Hyper"
meta_res_robust <- meta_res_robust %>%
  mutate(sigdir = sigdir)
meta_res_robust <- meta_res_robust %>%
  arrange(desc(sigdir))

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
                HetPVal,
                sig,
                sigdir)
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
                              "Heterogeneity p-value",
                              "Significance",
                              "Direction")

nDMPs <- nrow(meta_res_robust %>% filter(Significance=="Sig"))

write.table(meta_res_robust,
            file="METAL_muscle_training.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#Run random effects meta-analysis with metafor
library(metafor)
library(readxl)
dataset_summary <- read_excel('Datasets summary.xlsx',
                          na = c("","NA"))%>%
  dplyr::rename(prop_males = `% Male`)%>%
  filter(`Training study`=="Yes")
meta_res_robust <- read_delim("METAL_muscle_training.txt",
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
        file = "METAL_muscle_training.rds")

#Set working directory
setwd("C:/Users/e5103562/OneDrive - Victoria University/EWAS meta analysis of training/")
library(tidyverse)
L <- read_rds("METAL_muscle_training.rds")
library(metafor)


meta_RE <- function(tib = tib)
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
    return(c(as.numeric(model$b),
             model$se,
             model$pval))
  })
}

#Run on age-related CpGs
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

DMPs <- subset(to_write,
               `FDR age`<0.005)

subL <- L[intersect(names(L),DMPs$CpG)]
models <- t(sapply(subL,
                 meta_RE))
CpG <- rownames(models)
models <- as_tibble(models)
colnames(models) <- c("ES_RE","SE_RE","pval_RE")
models <- models %>%
  mutate(CpG = CpG,
         FDR = p.adjust(pval_RE))
write.table(models,
            file="metafor_muscle_training_RE.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

