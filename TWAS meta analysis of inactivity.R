setwd("C:/Users/e5103562/OneDrive - Victoria University/TWAS meta analysis of inactivity")

DIS <- read_tsv("Full results muscle disuse.txt") %>%
    dplyr::select(Gene.Name,
                  LFC.mean,
                  Stouffer.raw)%>%
    dplyr::rename(`Gene name` = `Gene.Name`,
                  `Effect size muscle disuse`=LFC.mean,
                  `P-value muscle disuse`=Stouffer.raw)%>%
    mutate(`FDR muscle disuse` = p.adjust(`P-value muscle disuse`))

