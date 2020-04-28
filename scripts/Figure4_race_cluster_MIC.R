library(tidyverse)
library(cowplot)
library(ggpubr)
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

valid_mics <- c(0.002, 0.003, 0.004, 0.006, 0.008, 0.012, 0.016, 0.023, 0.032, 0.047, 0.064, 0.094, 0.125, 0.19, 0.25, 0.38, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 32)

metadata_filtered <- metadata %>% mutate(cefixime = ifelse(cefixime %in% valid_mics, cefixime, NA),
                                         ceftriaxone = ifelse(ceftriaxone %in% valid_mics, cefixime, NA),
                                         ciprofloxacin = ifelse(ciprofloxacin %in% valid_mics, cefixime, NA),
                                         azithromycin = ifelse(azithromycin %in% valid_mics, cefixime, NA))

metadata <- metadata %>% 
  add_count(final_cluster) %>% 
  mutate(clustered = ifelse(n > 1, "clustered", "unclustered"))

metadata$final_cluster <- as.factor(as.character(SNPclusters.10.updated$final_cluster))

msm <- metadata %>% filter(partner == "MSM")
het <- metadata %>% filter(partner == "MSW" | partner == "WSM")

# Figure 4A

het.clustering.plot <- het %>%
  ggplot(aes(x = Race)) +
  geom_histogram(stat = "count") +
  facet_wrap(~clustered) +
  xlab("Race/Ethnicity") +
  ylab("Count")

het.chisq.clustered <- het %>%
  count(Race, clustered) %>%
  spread(Race, n, fill = 0) %>%
  column_to_rownames("clustered")

# Figure 4B

het.race.chisq <- chisq.test(het.chisq.clustered)

residuals.het.df <- data.frame(het.race.chisq$residuals)
residuals.het.df.long <- residuals.het.df %>%
  rownames_to_column(var = "clustering") %>%
  gather(Race, residual, -clustering)

residuals.het.heatmap <- ggplot(residuals.het.df.long, aes(x = Race, y = clustering)) +
  geom_tile(aes(fill = residual)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  ylab("") +
  xlab("Race/Ethnicity")

clustering.plot <- plot_grid(het.clustering.plot, residuals.het.heatmap, labels = c("A", "B"), ncol = 1, rel_heights = c(2, 1))


msm.chisq.clustered <- msm %>%
  count(Race, clustered) %>%
  spread(Race, n, fill = 0) %>%
  column_to_rownames("clustered")

msm.race.chisq <- chisq.test(msm.chisq.clustered)

# Supplementary Figure 4

msm.cro <- msm %>%
  ggplot(aes(x = Race, y = ceftriaxone)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1)),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5")) + 
  labs(x = "", y = "CRO MIC") + 
  stat_compare_means(label.y = 2)


msm.azi <- msm %>%
  ggplot(aes(x = Race, y = azithromycin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5, 2^7),
                     labels = c("0.008", "0.032", "0.125", "0.5", "2", "8", "32", "128")) +
  labs(x = "", y = "AZM MIC") +
  stat_compare_means(label.y =15)

msm.cip <- msm %>% 
  ggplot(aes(x = Race, y = ciprofloxacin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5", "2", "8", "32")) +
  labs(x = "Race/Ethnicity", y = "CIP MIC") +
  stat_compare_means(label.y = 15) +
  stat_compare_means(comparisons = list(c("NH-Black", "NH-White"), 
                                        c("NH-Other", "NH-Asian"),
                                        c("NH-White", "NH-Other")), label = "p.signif")

# Figure 4C

comparisons.het.cro <- list(c("Hispanic", "NH-White"), c("NH-Black", "NH-White"))
comparisons.het.azi <- list(c("Hispanic", "NH-White"), c("NH-Black", "NH-White"))
comparisons.het.cip <- list(c("Hispanic", "NH-Black"), c("NH-Black", "NH-White"))

het.cro <- het %>%
  ggplot(aes(x = Race, y = ceftriaxone)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1)),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5")) + 
  labs(x = "", y = "CRO MIC") + 
  stat_compare_means(label.y = 2) +
  stat_compare_means(comparisons = comparisons.het.cro, label = "p.signif")



het.azi <- het %>%
  ggplot(aes(x = Race, y = azithromycin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5, 2^7),
                     labels = c("0.008", "0.032", "0.125", "0.5", "2", "8", "32", "128")) +
  labs(x = "", y = "AZI MIC") +
  stat_compare_means(label.y = 7) +
  stat_compare_means(comparisons = comparisons.het.azi, label = "p.signif")

het.cip <- het %>% 
  ggplot(aes(x = Race, y = ciprofloxacin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5", "2", "8", "32")) +
  labs(x = "Race/Ethnicity", y = "CIP MIC") +
  stat_compare_means(label.y = 10) +
  stat_compare_means(comparisons = comparisons.het.cip, label = "p.signif")

het.mic <- plot_grid(het.cro,
                     het.azi, 
                     het.cip, 
                     ncol = 1, 
                     align = "v",
                     axis = "tb")

het.cluster.mic <- plot_grid(clustering.plot, het.mic, labels = c("", "C"), nrow = 1)


# lineage? 

clades <- read.table("data/major_clades.txt", header = F, sep = "\t")
colnames(clades) <- c("run_accession", "clade")
het <- het %>% left_join(clades)

comparisons.het.cro <- list(c("Hispanic", "NH-White"), c("NH-Black", "NH-White"))
comparisons.het.azi <- list(c("Hispanic", "NH-White"), c("NH-Black", "NH-White"))
comparisons.het.cip <- list(c("Hispanic", "NH-Black"), c("NH-Black", "NH-White"))

het.cro <- het %>% filter(clade == "MSM_clade") %>%
  ggplot(aes(x = Race, y = ceftriaxone)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1)),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5")) + 
  labs(x = "", y = "CRO MIC") + 
  stat_compare_means(label.y = 2) +
  stat_compare_means(comparisons = comparisons.het.cro, label = "p.signif")



het.azi <- het %>% filter(clade == "MSM_clade") %>%
  ggplot(aes(x = Race, y = azithromycin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5, 2^7),
                     labels = c("0.008", "0.032", "0.125", "0.5", "2", "8", "32", "128")) +
  labs(x = "", y = "AZI MIC") +
  stat_compare_means(label.y = 7) +
  stat_compare_means(comparisons = comparisons.het.azi, label = "p.signif")

het.cip <- het %>% filter(clade == "MSM_clade") %>%
  ggplot(aes(x = Race, y = ciprofloxacin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5", "2", "8", "32")) +
  labs(x = "Race/Ethnicity", y = "CIP MIC") +
  stat_compare_means(label.y = 10) +
  stat_compare_means(comparisons = comparisons.het.cip, label = "p.signif")

het.mic <- plot_grid(het.cro,
                     het.azi, 
                     het.cip, 
                     ncol = 1, 
                     align = "v",
                     axis = "tb")

het.cluster.mic <- plot_grid(clustering.plot, het.mic, labels = c("", "C"), nrow = 1)

# using filtered MICs

# Supplementary Figure 4

msm <- metadata_filtered %>% filter(partner == "MSM")
het <- metadata_filtered %>% filter(partner == "MSW" | partner == "WSM")

msm.cro <- msm %>%
  ggplot(aes(x = Race, y = ceftriaxone)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1)),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5")) + 
  labs(x = "", y = "CRO MIC") + 
  stat_compare_means(label.y = 2)


msm.azi <- msm %>%
  ggplot(aes(x = Race, y = azithromycin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5, 2^7),
                     labels = c("0.008", "0.032", "0.125", "0.5", "2", "8", "32", "128")) +
  labs(x = "", y = "AZI MIC") +
  stat_compare_means(label.y =15)

msm.cip <- msm %>% 
  ggplot(aes(x = Race, y = ciprofloxacin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5", "2", "8", "32")) +
  labs(x = "Race/Ethnicity", y = "CIP MIC") +
  stat_compare_means(label.y = 15) +
  stat_compare_means(comparisons = list(c("NH-Black", "NH-White"), 
                                        c("NH-Other", "NH-Asian"),
                                        c("NH-White", "NH-Other")), label = "p.signif")

# Figure 4C

comparisons.het.cro <- list(c("Hispanic", "NH-White"), c("NH-Black", "NH-White"))
comparisons.het.azi <- list(c("Hispanic", "NH-White"), c("NH-Black", "NH-White"))
comparisons.het.cip <- list(c("Hispanic", "NH-Black"), c("NH-Black", "NH-White"))

het.cro <- het %>%
  ggplot(aes(x = Race, y = ceftriaxone)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1)),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5")) + 
  labs(x = "", y = "CRO MIC") + 
  stat_compare_means(label.y = 2) +
  stat_compare_means(comparisons = comparisons.het.cro, label = "p.signif")



het.azi <- het %>%
  ggplot(aes(x = Race, y = azithromycin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5, 2^7),
                     labels = c("0.008", "0.032", "0.125", "0.5", "2", "8", "32", "128")) +
  labs(x = "", y = "AZI MIC") +
  stat_compare_means(label.y = 7) +
  stat_compare_means(comparisons = comparisons.het.azi, label = "p.signif")

het.cip <- het %>% 
  ggplot(aes(x = Race, y = ciprofloxacin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5", "2", "8", "32")) +
  labs(x = "Race/Ethnicity", y = "CIP MIC") +
  stat_compare_means(label.y = 10) +
  stat_compare_means(comparisons = comparisons.het.cip, label = "p.signif")

het.mic <- plot_grid(het.cro,
                     het.azi, 
                     het.cip, 
                     ncol = 1, 
                     align = "v",
                     axis = "tb")



