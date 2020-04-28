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

metadata$final_cluster <- as.factor(as.character(metadata$final_cluster))

# Figure 3A

partner.clustering.plot <- ggplot(metadata %>%
                                    filter(!is.na(partner)) %>%
                                    filter(partner != "TSM" & partner != "WSMW"), aes(x=partner))  +
  geom_histogram(stat="count") + facet_wrap(~clustered) +
  xlab("Sex of Sex Partner") +
  ylab("Count")

# Figure 3B

partner.clustered <- metadata %>%
  filter(!is.na(partner)) %>%
  filter(partner != "TSM" & partner != "WSMW") %>%
  count(partner, clustered) %>%
  spread(partner, n, fill = 0) %>%
  column_to_rownames("clustered")

partner.chisq <- chisq.test(partner.clustered)

residuals.df <- data.frame(partner.chisq$residuals)
residuals.df.long <- residuals.df %>%
  rownames_to_column(var = "clustering") %>%
  gather(partner, residual, -clustering)

residuals.heatmap <- ggplot(residuals.df.long, aes(x = partner, y = clustering)) +
  geom_tile(aes(fill = residual)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  ylab("") +
  xlab("Sex of Sex Partner")


clustering.plot <- plot_grid(partner.clustering.plot, residuals.heatmap, labels = c("A", "B"), ncol = 1, rel_heights = c(2, 1))


# Figure 3C

comparisons.cro <- list(c("MSM", "MSW"), c("MSM", "WSM"), c("MSMW", "WSM"))
comparisons.azi <- list(c("MSMW", "MSW"), c("MSMW", "WSM"), c("MSM", "MSW"), c("MSM", "WSM"))
comparisons.cip <- list(c("MSM", "MSMW"), c("MSM", "MSW"), c("MSM", "WSM"),  c("MSMW", "WSM"))
partner.cro <- metadata %>% filter(!is.na(partner)) %>%
  filter(partner != "TSM" & partner != "WSMW") %>%
  ggplot(aes(x = partner, y = ceftriaxone)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1)),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5")) + 
  labs(x = "", y = "CRO MIC") + 
  stat_compare_means(comparisons = comparisons.cro, label = "p.signif") +
  stat_compare_means(label.y = 2)


partner.azi <- metadata %>% filter(!is.na(partner)) %>%
  filter(partner != "TSM" & partner != "WSMW") %>%
  ggplot(aes(x = partner, y = azithromycin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5, 2^7),
                     labels = c("0.008", "0.032", "0.125", "0.5", "2", "8", "32", "128")) +
  labs(x = "", y = "AZM MIC") +
  stat_compare_means(comparisons = comparisons.azi, label = "p.signif") +
  stat_compare_means(label.y =15)

partner.cip <- metadata %>% filter(!is.na(partner)) %>%
  filter(partner != "TSM" & partner != "WSMW") %>%
  ggplot(aes(x = partner, y = ciprofloxacin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5", "2", "8", "32")) +
  labs(x = "Sex of Sex Partner", y = "CIP MIC") +
  stat_compare_means(comparisons = comparisons.cip, label = "p.signif") +
  stat_compare_means(label.y = 15)

partner.mic <- plot_grid(partner.cro,
                         partner.azi, 
                         partner.cip, 
                         ncol = 1, 
                         align = "v",
                         axis = "tb")

msm.cluster.mic <- plot_grid(clustering.plot, partner.mic, nrow = 1, labels = c("", "C"))


# how much of this association is due to lineage A and lineage B?

clades <- read.table("data/major_clades.txt", header = F, sep = "\t")
colnames(clades) <- c("run_accession", "clade")
metadata <- metadata %>% left_join(clades)


partner.cro.B <- metadata %>% filter(!is.na(partner), clade == "Het_clade") %>%
  filter(partner != "TSM" & partner != "WSMW") %>%
  ggplot(aes(x = partner, y = ceftriaxone)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1)),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5")) + 
  labs(x = "", y = "CRO MIC") + 

  stat_compare_means(label.y = 2)


partner.azi.B <- metadata %>% filter(!is.na(partner), clade == "Het_clade") %>%
  filter(partner != "TSM" & partner != "WSMW") %>%
  ggplot(aes(x = partner, y = azithromycin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5, 2^7),
                     labels = c("0.008", "0.032", "0.125", "0.5", "2", "8", "32", "128")) +
  labs(x = "", y = "AZI MIC") +

  stat_compare_means(label.y =15)

partner.cip.B <- metadata %>% filter(!is.na(partner), clade == "Het_clade") %>%
  filter(partner != "TSM" & partner != "WSMW") %>%
  ggplot(aes(x = partner, y = ciprofloxacin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5", "2", "8", "32")) +
  labs(x = "Sex of Sex Partner", y = "CIP MIC") +

  stat_compare_means(label.y = 15)

partner.mic.B <- plot_grid(partner.cro.B,
                         partner.azi.B, 
                         partner.cip.B, 
                         ncol = 1, 
                         align = "v",
                         axis = "tb")

partner.cro.A <- metadata %>% filter(!is.na(partner), clade == "MSM_clade") %>%
  filter(partner != "TSM" & partner != "WSMW") %>%
  ggplot(aes(x = partner, y = ceftriaxone)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1)),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5")) + 
  labs(x = "", y = "CRO MIC") + 
  
  stat_compare_means(label.y = 2)


partner.azi.A <- metadata %>% filter(!is.na(partner), clade == "MSM_clade") %>%
  filter(partner != "TSM" & partner != "WSMW") %>%
  ggplot(aes(x = partner, y = azithromycin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5, 2^7),
                     labels = c("0.008", "0.032", "0.125", "0.5", "2", "8", "32", "128")) +
  labs(x = "", y = "AZI MIC") +
  
  stat_compare_means(label.y =15)

comparisons.cip.A <- list(c("MSM", "MSMW"), c("MSM", "MSW"), c("MSM", "WSM"))

partner.cip.A <- metadata %>% filter(!is.na(partner), clade == "MSM_clade") %>%
  filter(partner != "TSM" & partner != "WSMW") %>%
  ggplot(aes(x = partner, y = ciprofloxacin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5", "2", "8", "32")) +
  labs(x = "Sex of Sex Partner", y = "CIP MIC") +
  stat_compare_means(comparisons = comparisons.cip.A, label = "p.signif") +
  stat_compare_means(label.y = 15)

partner.mic.A <- plot_grid(partner.cro.A,
                           partner.azi.A, 
                           partner.cip.A, 
                           ncol = 1, 
                           align = "v",
                           axis = "tb")



lineage.partner.mic <- plot_grid(partner.mic.A, partner.mic.B, align = "h", labels = c("A", "B"))


# analysis with filtered MICs

# Figure 3C

comparisons.cro <- list(c("MSM", "MSW"), c("MSM", "WSM"))
comparisons.azi <- list(c("MSM", "MSW"), c("MSM", "WSM"))
comparisons.cip <- list(c("MSM", "MSW"), c("MSM", "WSM"))
partner.cro <- metadata_filtered %>% filter(!is.na(partner)) %>%
  filter(partner != "TSM" & partner != "WSMW") %>%
  ggplot(aes(x = partner, y = ceftriaxone)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1)),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5")) + 
  labs(x = "", y = "CRO MIC") + 
  stat_compare_means(comparisons = comparisons.cro, label = "p.signif") +
  stat_compare_means(label.y = 2)


partner.azi <- metadata_filtered %>% filter(!is.na(partner)) %>%
  filter(partner != "TSM" & partner != "WSMW") %>%
  ggplot(aes(x = partner, y = azithromycin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5, 2^7),
                     labels = c("0.008", "0.032", "0.125", "0.5", "2", "8", "32", "128")) +
  labs(x = "", y = "AZI MIC") +
  stat_compare_means(comparisons = comparisons.azi, label = "p.signif") +
  stat_compare_means(label.y =15)

partner.cip <- metadata_filtered %>% filter(!is.na(partner)) %>%
  filter(partner != "TSM" & partner != "WSMW") %>%
  ggplot(aes(x = partner, y = ciprofloxacin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-7), 2^(-5), 2^(-3), 2^(-1), 2, 2^3, 2^5),
                     labels = c("0.002", "0.008", "0.032", "0.125", "0.5", "2", "8", "32")) +
  labs(x = "Sex of Sex Partner", y = "CIP MIC") +
  stat_compare_means(comparisons = comparisons.cip, label = "p.signif") +
  stat_compare_means(label.y = 15)

partner.mic <- plot_grid(partner.cro,
                         partner.azi, 
                         partner.cip, 
                         ncol = 1, 
                         align = "v",
                         axis = "tb")
