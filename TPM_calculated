##Import the count table of all virus sequences in all samples
countsMatrix <- read.table("all_samples_count_table.txt", sep = "\t", header = T)

##Import the length information of all virus sequences in all samples
len <- readxl::read_xlsx("len_43130.xlsx")

##merge files countsMatrix and len
countsMatrix_len <- dplyr::left_join(len, countsMatrix, by = "contig")

##Calculate TPM of samples in different groups 
library(tidyr) 
tpmMatrix <- countsMatrix_len %>%
  pivot_longer(c(-contig, -len),
               names_to = "Group",
               values_to = "SampleCounts") %>%
  group_by(Group) %>%
  mutate(SampleTPM = (((SampleCounts/len)*1e6)/sum(SampleCounts/len))) %>%
  pivot_wider(id_cols = "contig",
              names_from = "Group",
              values_from = "SampleTPM")

#Export file tpmMatrix, which contains TPM value of all virus sequences in all samples
write.table(tpmMatrix ,"TPM_43130.txt", sep = "\t", quote = F, row.names = F)
