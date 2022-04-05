library(tidyverse)
library(reshape2)

proteinMetrics <- read_tsv("./roc_protein.tsv", col_names = TRUE) %>% mutate (Parameter = paste0("(", MAPQ, ",", Filter_Abs, ",", Filter_Rel,")"))
rnaMetrics <- read_tsv("./roc_rna.tsv", col_names = TRUE) %>% mutate (Parameter = paste0("(", MAPQ, ",", Filter_Abs, ",", Filter_Rel,")"))


#ROC curve: X = FPR, Y=TPR
p<-ggplot(proteinMetrics, aes(x=FPR, y=TPR, label=Parameter)) + 
  geom_point()+
#  geom_text(hjust=0, vjust=0, size = 2) + 
  ylim(0,1)
ggsave("protein.png", p)

p<-ggplot(rnaMetrics, aes(x=FPR, y=TPR, label=Parameter)) + 
  geom_point()+
 # geom_text(hjust=0, vjust=0, size = 2) + 
  ylim(0,1)
ggsave("rna.png", p)

datasetMatrix <- read_tsv("MockAssessmentMatrix.tsv", col_names = TRUE) %>% select (Class, Dataset, PPV, TPR)

p <- ggplot (datasetMatrix, aes(x=TPR, y=PPV, color = Class, shape = Dataset)) +
    geom_point(size = 5) +
    scale_shape_manual(values = c(15,16,17,18))+
    xlab("Recall") + 
    ylab("Precision")

ggsave("dataset.svg", p, width = 10, height = 5)
