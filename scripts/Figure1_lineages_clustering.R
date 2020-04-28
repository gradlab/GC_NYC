library(tidyverse)
library(cowplot)
library(viridis)
theme_set(theme_cowplot(font_size=10))

metadata <- read.table("data/Mortimer_NYC_SupplementaryTable1.tsv", 
                       sep = "\t",
                       header = T,
                       stringsAsFactors = F,
                       comment.char = "")

# Figure 1A made in ITOL

# Figure 1B

nyc_clades <- read.table("data/major_clades.txt", sep = "\t", comment.char = "!", stringsAsFactors = F)
colnames(nyc_clades) <- c("run_accession", "clade")

clades <- metadata %>% left_join(nyc_clades)

clades <- clades %>% 
  dplyr::select(run_accession, Race, partner, clade) %>% 
  mutate(orientation = ifelse(partner == "MSM" | partner == "MSMW" | partner == "TSM", "MSM",
                              ifelse(partner == "WSM" | partner == "WSMW" | partner == "MSW", "heterosexual", "Unknown"))) %>%
  mutate(Lineage = ifelse(clade == "MSM_clade", "Lineage A", "Lineage B"))

clades <- clades %>%
  mutate(orientation = ifelse(is.na(orientation), "Unknown", orientation)) %>%
  mutate(partner = ifelse(is.na(partner), "Unknown", partner))


clade.partner.plot <- clades %>%
  ggplot(aes(x=fct_infreq(partner))) +
  geom_histogram(stat="count") +
  facet_wrap(~Lineage, strip.position = "top") +
  ylab("Number of Isolates") +
  xlab("Sex of Sex Partner") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

clades.chi2.table <- clades %>% 
  filter(!is.na(partner)) %>% 
  count(Lineage, partner) %>% 
  spread(partner, n, fill = 0) %>%
  column_to_rownames("Lineage")

clade.partner.chi2 <- chisq.test(clades.chi2.table)

# Figure 1C

## describe clusters by proportion of each demo

common_clusters <- metadata %>% count(final_cluster) %>% arrange(desc(n)) %>% filter(n >= 10)

metadata.common <- metadata %>% filter(final_cluster %in% common_clusters$final_cluster)

## cluster size



metadata$final_cluster <- as.factor(as.character(metadata$final_cluster))
metadata <- metadata %>% 
  add_count(final_cluster) 

clusters.common <- metadata %>% filter(n >= 10)

clusters.common.size <- clusters.common %>%
  ggplot(aes(x = fct_infreq(factor(final_cluster)))) +
  geom_bar() +
  coord_flip() + 
  labs(x = "Cluster", y = "Number of Isolates", title = "Cluster Size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#age

clusters.common.age <- clusters.common %>% 
  ggplot(aes(x = fct_infreq(factor(final_cluster)), y = Age)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) + 
  coord_flip() +
  labs(y = "Age (Years)", title = "Patient Age") +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) 

# network

clusters.common.network <- clusters.common %>%
  ggplot(aes(x = fct_infreq(factor(final_cluster)), fill = fct_infreq(factor(partner)))) +
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
  ggplot(aes(x = fct_infreq(factor(final_cluster)), fill = fct_infreq(factor(Race)))) +
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
  ggplot(aes(x = fct_infreq(factor(final_cluster)), fill = fct_infreq(factor(hivstatus)))) +
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
  ggplot(aes(x = fct_infreq(factor(final_cluster)), fill = fct_infreq(factor(CT_any)))) +
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