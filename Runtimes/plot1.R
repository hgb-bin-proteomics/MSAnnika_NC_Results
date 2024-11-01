library(tidyverse)
library(readxl)
library(ggsci)
# https://github.com/nanxstats/ggsci

df <- read_excel("plot1.xlsx")
df$Database <- as.factor(df$Database)

# export 1200 x 800

lw = 0.5
ggplot(df, aes(Database, Mean, color = Tool)) +
  scale_color_npg() +
  geom_point(size = 3) + 
  geom_line(aes(group = Tool), linewidth = lw) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.1, linewidth = lw) +
  xlab("Database size (number of proteins)") +
  ylab("Average runtime in seconds") +
  theme_minimal(base_size = 18)
