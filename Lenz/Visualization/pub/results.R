library(tidyverse)
library(readxl)

df <- read_excel("results.xlsx")

df <- df %>% 
  mutate(Type = factor(Type, levels = c("Intra", "Inter true", "Inter false")))
  
# export to 900 x 800

ggplot(df, aes(x=Tool, y=Crosslinks, fill=Type)) +
  geom_bar(stat="identity", color="black", width = 0.7, position = position_stack(reverse = T)) +
  scale_fill_manual(values = c("#212529", "#ced4da", "#ef233c"),
                    labels = c("Intra", "Inter true", "Inter false")) +
  #ggtitle("Dataset of Lenz et al., 2021:\nNumber of identified crosslinks at 1% estimated FDR\n(crosslinker: BS3)") +
  xlab("Pipeline") +
  ylab("Number of identified crosslinks at 1% estimated FDR") +
  ylim(c(0, 7100)) +
  labs(fill="Crosslinks") +
  geom_text(aes(y = txt_pos, label=round(Crosslinks, 0)), color = rep(c("white", "black", "#ef233c"), 2), size = 5.0, 
            vjust = 0.5, hjust = 0.5) +
  geom_text(aes(y = 7100, label = paste0(round(FDR, 2), "%")), size = 5.0, angle = 0, hjust = 0.5) +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  guides(fill = guide_legend(ncol = 1, label.position = "right"))
