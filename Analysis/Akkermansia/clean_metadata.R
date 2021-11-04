library(here)
library(openxlsx)
library(dplyr)
library(readr)
df <- read.xlsx(here("Raw", "Akkermansia", "metadata.xlsx"), sheet = 2) %>% 
    filter(!SGB %in% c(9223, 9224, 9227, 92248),
           !is.na(Gender), !is.na(BMI), !is.na(Disease)) %>%
    select(GenomeID, Completeness, Contamination, BMI, Disease, Gender) %>%
    filter(Completeness > 95, Contamination < .5) %>%
    mutate(Disease = Disease == "healthy",
           Disease = as.numeric(Disease),
           Obese = (BMI >= 30) %>% as.numeric(),
           Gender = (Gender == "male") %>% as.numeric())

write.table(df, here("Raw", "Akkermansia", "metadata.tsv"), sep = "\t")
write.table(df %>% select(GenomeID) %>% mutate(GenomeID = paste0(GenomeID, ".fna")),
  here("Raw", "Akkermansia", "strains.tsv"), sep = "\t", 
  col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(df %>% select(Disease), here("Raw", "Akkermansia", "Disease.tsv"), sep = "\t", 
  col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(df %>% select(Gender), here("Raw", "Akkermansia", "Gender.tsv"), sep = "\t", 
  col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(df %>% select(Obese), here("Raw", "Akkermansia", "Pheno.tsv"), sep = "\t", 
  col.names = FALSE, row.names = FALSE, quote = FALSE)


write.table(df %>% select(GenomeID, Obese) %>% 
              mutate(Path = paste0("/Raw/Akkermansia/Genomes/", GenomeID, ".fna")) %>%
              rename("Pheno" = Obese) %>%
              select(GenomeID, Pheno, Path), 
  here("Raw", "Akkermansia", "strains"), sep = "\t", 
  row.names = FALSE, quote = FALSE)
