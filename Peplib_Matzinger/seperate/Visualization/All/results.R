library(tidyverse)
library(readxl)

df <- read_excel("results.xlsx")

df <- df %>% 
  mutate(errorbar_position = ifelse(Type == "True XL", meanTrue, meanTrue + meanFalse)) %>% 
  mutate(txt_position = ifelse(Type == "True XL", round(meanTrue / 2), round(meanTrue + meanFalse) + 25)) %>% 
  mutate(Type = factor(Type, levels = c("True XL", "False XL")))
  
ggplot(df, aes(x=Tool, y=Mean, fill=Type)) +
  geom_bar(stat="identity", color="black", width = 0.7, position = position_stack(reverse = T)) +
  scale_fill_manual(values = c("#212529", "#ced4da"),
                    labels = c("True", "False")) +
  geom_errorbar(aes(x=Tool, ymin=errorbar_position-SD, ymax=errorbar_position+SD), 
                colour = rep(c("#40916c", "#ef233c"), 6), width = 0.3, stat = "identity", linewidth = 1.0,
                position = position_dodge(0.5)) +
  ggtitle("Dataset of synthetic peptides by Matzinger et al., 2022:\nNumber of identified crosslinks per tool at 1% estimated FDR\n(3 replicates, crosslinker: ADH)") +
  xlab("Tool") +
  ylab("Number of identified Crosslinks (mean) at 1% estimated FDR") +
  ylim(c(0, 200)) +
  labs(fill="Crosslinks") +
  geom_text(aes(y=txt_position, label=round(Mean, 0)), color = rep(c("white", "black"), 6), size=5.0) +
  geom_text(aes(y = 185, label = paste0(round(FDR, 2), "%")), size=5.0, angle = 0, hjust = 0.5) +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  guides(fill = guide_legend(ncol = 1, label.position = "right"))
