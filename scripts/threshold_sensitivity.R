# repeat all cluster analyses with 12 SNP cutoff

library(tidyverse)
library(phytools)

theme_set(theme_cowplot(font_size=10))

SNPclusters.12 <- read.table("data/clustering/2019-07-24_nyc_12SNP_clusters.txt", sep = "\t", header = F, comment.char = "!")
colnames(SNPclusters.12) <- c("run_accession", "cluster")
SNPclusters.12 <- SNPclusters.12 %>%
  add_row(run_accession = "ERR2631833", cluster = 482) %>%
  add_row(run_accession = "ERR3201242", cluster = 483)

# get missing isolates from clusters using MRCA and descendants

gcTree <- read.tree("data/gubbins/2019-01-23_nyc_metadataComplete.final_tree.tre")

get.missing.isolates <- function(cluster.name, clusters.df, tree) {
  original.isolates <- filter(clusters.df, cluster == cluster.name)$run_accession
  cluster.mrca <- findMRCA(tree, tips=as.character(original.isolates), type = "node")
  cluster.descendants <- getDescendants(tree, cluster.mrca)
  descendant.names <- na.omit(gcTree$tip.label[cluster.descendants])
  new.cluster.df <- data.frame("run_accession" = as.character(descendant.names), "new_cluster" = cluster.name)
  return(new.cluster.df)
}

new.clusters <- get.missing.isolates(1, SNPclusters.12, gcTree)

for (i in 2:140) {
  nc <- get.missing.isolates(i, SNPclusters.12, gcTree)
  new.clusters <- bind_rows(new.clusters, nc)
}

new.clusters <- new.clusters %>% 
  add_count(new_cluster) %>% 
  arrange(desc(n)) %>% 
  distinct(run_accession, .keep_all = TRUE) %>% 
  dplyr::select(run_accession, new_cluster)

SNPclusters.12.updated <- left_join(SNPclusters.12, new.clusters)

SNPclusters.12.updated <- SNPclusters.12.updated %>% 
  mutate(final_cluster = ifelse(!is.na(new_cluster), new_cluster, cluster)) %>% 
  add_count(final_cluster) %>% 
  mutate(clustered = ifelse(n > 1, "clustered", "unclustered"))

SNPclusters.12.updated$final_cluster <- as.factor(as.character(SNPclusters.12.updated$final_cluster))

SNPclusters.12.summary <- SNPclusters.12.updated %>% select(run_accession, final_cluster)
colnames(SNPclusters.12.summary) <- c("run_accession", "final_cluster_12")

# top clusters (>=10 samples)

library(cowplot)
library(viridis)
theme_set(theme_cowplot(font_size=10))

metadata <- read.table("data/Mortimer_NYC_SupplementaryTable1.tsv", 
                       sep = "\t",
                       header = T,
                       stringsAsFactors = F,
                       comment.char = "")
metadata <- metadata %>% left_join(SNPclusters.12.summary)

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

common_clusters <- metadata %>% count(final_cluster) %>% arrange(desc(n)) %>% filter(n >= 10)

metadata.common <- metadata %>% filter(final_cluster %in% common_clusters$final_cluster)

## cluster size



metadata$final_cluster_12 <- as.factor(as.character(metadata$final_cluster_12))
metadata <- metadata %>% 
  add_count(final_cluster_12) 

clusters.common <- metadata %>% filter(n >= 10)

clusters.common.size <- clusters.common %>%
  ggplot(aes(x = fct_infreq(factor(final_cluster_12)))) +
  geom_bar() +
  coord_flip() + 
  labs(x = "Cluster", y = "Number of Isolates", title = "Cluster Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#age

clusters.common.age <- clusters.common %>% 
  ggplot(aes(x = fct_infreq(factor(final_cluster_12)), y = Age)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) + 
  coord_flip() +
  labs(y = "Age (Years)", title = "Patient Age") +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) 

# network

clusters.common.network <- clusters.common %>%
  ggplot(aes(x = fct_infreq(factor(final_cluster_12)), fill = fct_infreq(factor(partner)))) +
  geom_bar(position = "fill") +
  scale_fill_viridis(option = "magma", discrete = T) +
  coord_flip() +
  labs(title = "Sex of Sex Partner") +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

# race/ethnicity

clusters.common.race <- clusters.common %>%
  ggplot(aes(x = fct_infreq(factor(final_cluster_12)), fill = fct_infreq(factor(Race)))) +
  geom_bar(position = "fill") +
  scale_fill_viridis(option = "viridis", discrete = T) +
  coord_flip() +
  labs(title = "Race/Ethnicity") +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

# hiv

clusters.common.hiv <- clusters.common %>%
  ggplot(aes(x = fct_infreq(factor(final_cluster_12)), fill = fct_infreq(factor(hivstatus)))) +
  geom_bar(position = "fill") +
  scale_fill_viridis(option = "plasma", discrete = T) +
  coord_flip() +
  labs(title = "HIV Status") +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

# chlamydia

clusters.common.ct <- clusters.common %>%
  ggplot(aes(x = fct_infreq(factor(final_cluster_12)), fill = fct_infreq(factor(CT_any)))) +
  geom_bar(position = "fill") +
  scale_fill_viridis(option = "inferno", discrete = T, labels = c("Negative", "Positive")) +
  coord_flip() +
  labs(title = "Chlamydia") +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

plot_grid(clusters.common.size,
          clusters.common.age, 
          clusters.common.network, 
          clusters.common.race, 
          clusters.common.hiv, 
          clusters.common.ct, 
          nrow = 1, 
          align = "h",
          axis = "tb")

# cluster MIC




clusters.common.cro <- clusters.common %>% 
  ggplot(aes(x = fct_infreq(factor(final_cluster_12)), y = ceftriaxone)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-8), 2^(-7), 2^(-6), 2^(-5), 2^(-4)),
                     labels = c("0.002", "0.004", "0.008", "0.016", "0.032", "0.064")) +
  coord_flip() +
  labs(x = "Cluster", y = "CRO MIC", title = "Ceftriaxone") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0.012, linetype= "dashed", color = "black", size = 1)

clusters.common.azi <- clusters.common %>% 
  ggplot(aes(x = fct_infreq(factor(final_cluster_12)), y = azithromycin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-6), 2^(-5), 2^(-4), 2^(-3), 2^(-2), 2^(-1), 2^(0)),
                     labels = c("0.016", "0.032", "0.064", "0.125", "0.25", "0.5", "1")) +
  coord_flip() +
  labs(x = "Cluster", y = "AZI MIC", title = "Azithromycin") +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0.19, linetype= "dashed", color = "black", size = 1)

clusters.common.cip <- clusters.common %>% 
  ggplot(aes(x = fct_infreq(factor(final_cluster_12)), y = ciprofloxacin)) + 
  geom_boxplot() + 
  scale_y_continuous(trans='log2',
                     breaks = c(2^(-9), 2^(-8), 2^(-7), 2^(-6), 2^(-5)),
                     labels = c("0.002", "0.004", "0.008", "0.016", "0.032")) +
  coord_flip() +
  labs(x = "Cluster", y = "CIP MIC", title = "Ciprofloxacin") +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 0.004, linetype= "dashed", color = "black", size = 1)

plot_grid(clusters.common.cro,
          clusters.common.azi, 
          clusters.common.cip, 
          nrow = 1, 
          align = "h",
          axis = "tb")

# clustered/unclustered demographics

metadata <- metadata %>% 
  add_count(final_cluster_12) %>% 
  mutate(clustered = ifelse(n > 1, "clustered", "unclustered"))
partner.clustering.plot <- ggplot(metadata %>%
                                    filter(!is.na(partner)) %>%
                                    filter(partner != "TSM" & partner != "WSMW"), aes(x=partner))  +
  geom_histogram(stat="count") + facet_wrap(~clustered) +
  xlab("Sex of Sex Partner") +
  ylab("Count")


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

