setwd("/Users/wrq796/OneDrive - Victoria University/TWAS meta analysis of training")

#Load extrameta results for exercise
#Extrameta
load("/Users/wrq796/OneDrive - Victoria University/TWAS meta analysis of training/meta_analysis_results.RData")

#Subset the results from chronic, muscle
chronic_muscle <- all_meta_analysis_res$`longterm,muscle`
#Have a look at number of genes in muscle for chronic training
length(chronic_muscle) #14,309 genes present
#Convert to gene symbol
library(org.Hs.eg.db)
library(annotate)
names(chronic_muscle) <- getSYMBOL(names(chronic_muscle),
                                   data='org.Hs.eg')
#Remove NA gene symbols from results
chronic_muscle <- chronic_muscle[-which(is.na(names(chronic_muscle)))]

#Extract model type (i.e. simple model or moderators), effect size, standard error and p-value for each gene.
#For each gene, the first entry is the best fit
my_fn <- function(l)
{
    tib <- tibble(mod_type = sub(".*:", "", names(l)[1]),
                  ES = l[[1]]$coeffs["intrcpt","beta"],
                  SE = l[[1]]$coeffs["intrcpt","se"],
                  p = l[[1]]$coeffs["intrcpt","pval"],
                  I2 = l[[1]]$I2)
    return(tib)
}

results <- lapply(chronic_muscle,
                  my_fn)

#Unlist and bind into a single table
all_results <- results %>%
    purrr::reduce(bind_rows)

#Add gene name
all_results <- all_results %>%
    mutate(Gene = names(chronic_muscle),
           ES = log(2^ES))%>%
    dplyr::rename(`logFC after training` = ES,
                  `SE training` = SE,
                  `P-value training` = p)

#Save results
write.table(all_results,
            file="Extrameta all genes results.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep="\t")
