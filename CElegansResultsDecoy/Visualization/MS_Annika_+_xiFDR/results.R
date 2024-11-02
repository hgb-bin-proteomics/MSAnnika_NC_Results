library(tidyverse)
library(readxl)

df <- read_excel("results.xlsx")

df <- df %>% 
  mutate(Database = factor(Database, levels = c("Full Proteome", "Filtered Proteome")))
  
# export to 1200 x 800

ggplot(df, aes(x=Replicate, y=Crosslinks, fill=Database)) +
  geom_bar(stat="identity", color="black", width = 0.7, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#212529", "#ced4da"),
                    labels = c("Full Proteome", "Filtered Proteome")) +
  #ggtitle("Dataset of C. elegans nuclei by Müller et al., 2024:\nNumber of identified crosslinks per replicate at 1% estimated FDR\n(3 replicates, crosslinker: DSG)") +
  xlab("Replicate") +
  ylab("Number of identified crosslinks at 1% estimated FDR") +
  ylim(c(0, 330)) +
  labs(fill="Protein Database") +
  geom_text(aes(y=Crosslinks/2, label=round(Crosslinks, 0)), color = rep(c("black", "white"), 3), size=5.0, position = position_dodge(width = 0.7)) +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  guides(fill = guide_legend(ncol = 1, label.position = "right"))
