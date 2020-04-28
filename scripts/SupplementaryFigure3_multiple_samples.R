library(tidyverse)
library(cowplot)
library(lubridate)
theme_set(theme_cowplot(font_size=10))

pairwise_differences.recombination_removed <- read.table("data/distances/2019-02-20_nyc_gubbins_distances.txt", 
                                                         header=T, 
                                                         sep = "\t", 
                                                         comment.char = "")

pairwise.recombination_removed.long <- gather(pairwise_differences.recombination_removed, "isolate2", "snp_difference", -snp.dists.0.6)
colnames(pairwise.recombination_removed.long) <- c("isolate1", "isolate2", "snp_differences")
pairwise.recombination_removed.long <- pairwise.recombination_removed.long %>%
  mutate(isolate2 = str_replace_all(isolate2, c("X"="", "\\." = "#")))

metadata <- read.table("data/Mortimer_NYC_SupplementaryTable1.tsv", 
                       sep = "\t",
                       header = T,
                       stringsAsFactors = F,
                       comment.char = "")

metadata$Dateofcollection <- ymd(metadata$Dateofcollection)

multiple_samples <- metadata %>% select(run_accession, PatientNumber, Dateofcollection) %>%
  filter(!is.na(PatientNumber)) %>%
  group_by(PatientNumber) %>%
  mutate(sample_count = n()) %>%
  filter(sample_count > 1) %>%
  arrange(PatientNumber, Dateofcollection) %>%
  mutate(visit_count = n_distinct(Dateofcollection))

multiple_visits <- multiple_samples %>% filter(visit_count > 1)

multiple_sites <- multiple_samples %>% filter(sample_count != visit_count)



multiple_sites_pairs <- multiple_sites %>%
  group_by(PatientNumber, Dateofcollection) %>%
  summarize(samples = paste(run_accession, collapse=",")) %>%
  separate(samples, c('isolate1', 'isolate2'), sep = ",") %>%
  filter(!is.na(isolate2))


multiple_sites_pairs.distance.norecombination <- multiple_sites_pairs %>% left_join(pairwise.recombination_removed.long)

# all patients with multiple visits and multiple samples at 1 visit had the same strain at both sites
# had to add pairwise comparison for one patient who had 3 isolates

multiple_visits_pairs <- multiple_visits %>% 
  group_by(PatientNumber, Dateofcollection) %>%
  summarize(samples = paste(run_accession, collapse=",")) %>%
  separate(samples, c('isolate1_sameday', 'isolate2_sameday'), sep = ",") %>%
  ungroup() %>%
  select(-isolate2_sameday) %>%
  group_by(PatientNumber) %>%
  summarize(samples = paste(isolate1_sameday, collapse=",")) %>%
  separate(samples, c('isolate1', 'isolate2', 'isolate3'), sep = ",") %>%
  select(-isolate3) %>%
  add_row(PatientNumber = "196", isolate1 = "ERR3200897", isolate2 = "ERR1204773")


multiple_visits_pairs.distance.norecombination <- multiple_visits_pairs %>% left_join(pairwise.recombination_removed.long)


### SNP distance vs time distance

collection_dates <- metadata %>% select(run_accession, Dateofcollection)

colnames(collection_dates) <- c("isolate1", "date1")
multiple_visits_pairs.distance.norecombination.dates <- multiple_visits_pairs.distance.norecombination %>% left_join(collection_dates)

colnames(collection_dates) <- c("isolate2", "date2")
multiple_visits_pairs.distance.norecombination.dates <- multiple_visits_pairs.distance.norecombination.dates %>% left_join(collection_dates)

multiple_visits_pairs.distance.norecombination.dates <- multiple_visits_pairs.distance.norecombination.dates %>%
  mutate(date_difference = abs(date2-date1))

multi_visit_plot <- ggplot(multiple_visits_pairs.distance.norecombination.dates, aes(x=date_difference, y=snp_differences)) +
  geom_point() + 
  xlab("Temporal Distance (Days)") +
  ylab("SNP Distance") 

multi_site_plot <- ggplot(multiple_sites_pairs.distance.norecombination,
                          aes(x = snp_differences)) +
  geom_histogram(bins = 100) +
  xlab("SNP Distance") +
  ylab("Patient Count")

multi_site_zoom <- ggplot(multiple_sites_pairs.distance.norecombination %>% filter(snp_differences < 100),
                          aes(x = snp_differences)) +
  geom_histogram(bins = 100) +
  xlab("SNP Distance") +
  ylab("Patient Count")

multi_site_combine <- ggdraw() + 
  draw_plot(multi_site_plot, 0, 0, 1, 1) +
  draw_plot(multi_site_zoom, 0.5, 0.5, 0.5, 0.5)

multi_site_visit <- plot_grid(multi_site_combine, multi_visit_plot, labels = c("A", "B"), nrow = 1)
