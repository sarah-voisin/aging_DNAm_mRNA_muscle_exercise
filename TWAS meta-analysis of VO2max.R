#Set working directory
setwd("C:/Users/e5103562/OneDrive - Victoria University/TWAS meta analysis of VO2max/")

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

sapply(files[3:4],
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
meta_res <- read_tsv("METAL_muscle_VO2max.TBL")%>%
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

#View results including at least 2 studies
meta_res_robust <- meta_res %>% dplyr::filter(HetDf>=1)
padj <- p.adjust(meta_res_robust$`P-value`,method="BH")
meta_res_robust <- meta_res_robust %>%
    mutate(`Adjusted p-value` = padj)
meta_res_robust <- meta_res_robust %>%
    mutate(sig = ifelse(`Adjusted p-value`<0.005,
                        "Sig",
                        "Not Sig"))
nrow(meta_res_robust) #17,405

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
    dplyr::rename(`log2FC per unit of VO2max` = Effect,
                  SE = StdErr,
                  FDR = `Adjusted p-value`,
                  `Number of studies` = HetDf,
                  `Heterogeneity index (I2)` = HetISq,
                  `Heterogeneity p-value` = HetPVal,
                  Significance = sig,
                  Direction = sigdir)

write.table(meta_res_robust,
            file="METAL_muscle_VO2max.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

#Run fixed effects meta-analysis with metafor to compare results with METAL
library(metafor)
library(readxl)
dataset_summary <- read_excel('C:/Users/e5103562/OneDrive - Victoria University/mRNA datasets.xlsx',
                              na = c("","NA"),
                              n_max = 20)%>%
    dplyr::rename(prop_males = `% males`)%>%
    dplyr::filter(VO2max=="Yes"&
                      !is.na(`VO2max mean ± SD`))
meta_res_robust <- read_delim("METAL_muscle_VO2max.txt",
                              lazy = FALSE)

#Create empty list
L.names <- meta_res_robust$`Gene name`
L <- vector("list", length(L.names))
names(L) <- L.names

#Initiate list
tib <- tibble(`Gene name` = meta_res_robust$`Gene name`,
              Dataset = rep("Meta-analysis",nrow(meta_res_robust)),
              ES = meta_res_robust$`log2FC per unit of VO2max`,
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
        file = "METAL_muscle_VO2max.rds")

#Set working directory
setwd("C:/Users/e5103562/OneDrive - Victoria University/TWAS meta analysis of VO2max/")
library(tidyverse)
L <- read_rds("METAL_muscle_VO2max.rds")
library(metafor)

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

#Run on all age-related genes
DEGs <- read_delim("C:/Users/e5103562/OneDrive - Victoria University/TWAS meta analysis of age/metafor_muscle_age_RE.txt")%>%
    dplyr::filter(FDR < 0.005)%>%
    pull(Gene)
subL <- L[intersect(names(L),DEGs)] #330/330 DEGs tested

models <- t(sapply(subL,
                   meta_RE))
gene <- rownames(models)
models <- as_tibble(models)
colnames(models) <- c("ES_RE","SE_RE","pval_RE")
models <- models %>%
    mutate(Gene = gene,
           FDR = p.adjust(pval_RE))
write.table(models,
            file="metafor_muscle_VO2max_RE.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")

