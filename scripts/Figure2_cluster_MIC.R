library(tidyverse)
library(cowplot)
theme_set(theme_cowplot(font_size=10))

metadata <- read.table("data/Mortimer_NYC_SupplementaryTable1.tsv", 
                       sep = "\t",
                       header = T,
                       stringsAsFactors = F,
                       comment.char = "")

metadata <- metadata %>% mutate(cefixime = ifelse(cefixime == "<=0.016" | cefixime == "<= 0.016" | cefixime == "<0.016", "0.016", cefixime),
                              ceftriaxone = ifelse(ceftriaxone == "<0.002" | ceftriaxone == "<= 0.002"| ceftriaxone == "<=0.002" | ceftriaxone == "<0.002", "0.002",
                                                   ifelse(ceftriaxone == "<0.016" | ceftriaxone ==  "<=0.016" | ceftriaxone == "<= 0.016", "0.016", 
                                                          ifelse(ceftriaxone == "0.5 NS", "0.5", ceftriaxone))),
                              azithromycin = ifelse(azithromycin == "<=0.016", "0.016",
                                                    ifelse(azithromycin == "0 .047              ", "0.047", azithromycin)),
                              ciprofloxacin = ifelse(ciprofloxacin == "<0.002" | ciprofloxacin == "<=0.002" | ciprofloxacin == ">=0.002" | ciprofloxacin == "<= 0.002" | ciprofloxacin == "0.0002", "0.002",
                                                     ifelse(ciprofloxacin == ">=32" | ciprofloxacin == ">= 32" | ciprofloxacin == ">32", "32", ciprofloxacin)))
metadata <- metadata %>% mutate(cefixime = as.numeric(cefixime),
                              ceftriaxone = as.numeric(ceftriaxone),
                              azithromycin = as.numeric(azithromycin),
                              ciprofloxacin = as.numeric(ciprofloxacin))

metadata$final_cluster <- as.factor(as.character(metadata$final_cluster))
metadata <- metadata %>% 
  add_count(final_cluster) 

clusters.common <- metadata %>% filter(n >= 10)

uncommon <- metadata %>% filter(n < 10)

clusters.common.cro <- clusters.common %>% 
  ggplot(aes(x = fct_infreq(factor(final_cluster)), y = ceftriaxone)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-8), 2^(-7), 2^(-6), 2^(-5), 2^(-4), 2^(-3), 2^(-2)),
                     labels = c("0.002", "0.004", "0.008", "0.016", "0.032", "0.064", "0.125", "0.25")) +
  coord_flip() +
  labs(x = "Cluster", y = "CRO MIC", title = "Ceftriaxone") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0.012, linetype= "dashed", color = "black", size = 1) + 
  geom_hline(yintercept = 0.25, linetype= "dashed", color = "red", size = 1)

clusters.common.azi <- clusters.common %>% 
  ggplot(aes(x = fct_infreq(factor(final_cluster)), y = azithromycin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-6), 2^(-5), 2^(-4), 2^(-3), 2^(-2), 2^(-1), 2^(0)),
                     labels = c("0.016", "0.032", "0.064", "0.125", "0.25", "0.5", "1")) +
  coord_flip() +
  labs(x = "Cluster", y = "AZM MIC", title = "Azithromycin") +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0.19, linetype= "dashed", color = "black", size = 1) +
  geom_hline(yintercept = 1, linetype= "dashed", color = "red", size = 1)

clusters.common.cip <- clusters.common %>% 
  ggplot(aes(x = fct_infreq(factor(final_cluster)), y = ciprofloxacin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-8), 2^(-7), 2^(-6), 2^(-5), 2^(-4)),
                     labels = c("0.002", "0.004", "0.008", "0.016", "0.032", "0.064")) +
  coord_flip() +
  labs(x = "Cluster", y = "CIP MIC", title = "Ciprofloxacin") +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0.004, linetype= "dashed", color = "black", size = 1) +
  geom_hline(yintercept = 0.06, linetype= "dashed", color = "red", size = 1)

plot_grid(clusters.common.cro,
          clusters.common.azi, 
          clusters.common.cip, 
          nrow = 1, 
          align = "h",
          axis = "tb")


