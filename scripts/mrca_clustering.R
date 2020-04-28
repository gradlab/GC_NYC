library(tidyverse)
library(phytools)
library(lubridate)

theme_set(theme_cowplot(font_size=10))

SNPclusters.10 <- read.table("data/clustering/2019-01-28_nyc_10SNP_clusters.txt", sep = "\t", header = F, comment.char = "!")
colnames(SNPclusters.10) <- c("run_accession", "cluster")

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

new.clusters <- get.missing.isolates(1, SNPclusters.10, gcTree)

for (i in 2:147) {
  nc <- get.missing.isolates(i, SNPclusters.10, gcTree)
  new.clusters <- bind_rows(new.clusters, nc)
}

metadata <- read.table("data/Mortimer_NYC_SupplementaryTable1.tsv", 
                       sep = "\t",
                       header = T,
                       stringsAsFactors = F,
                       comment.char = "")
SNPclusters.10.metadata <- SNPclusters.10 %>% left_join(metadata) %>% mutate(CollectionDate = parse_date_time(CollectionDate, "ym"))

missing.isolates.dates <- function(cluster.name, clusters.df, tree) {
  cluster <- filter(clusters.df, cluster == cluster.name)
  original.isolates <- filter(clusters.df, cluster == cluster.name)$run_accession
  original.isolates.start <- min(cluster$CollectionDate)
  original.isolates.stop <- max(cluster$CollectionDate)
  cluster.mrca <- findMRCA(tree, tips=as.character(original.isolates), type = "node")
  cluster.descendants <- getDescendants(tree, cluster.mrca)
  descendant.names <- na.omit(gcTree$tip.label[cluster.descendants])
  new.isolates <- setdiff(descendant.names, original.isolates)
  new_isolates <- filter(clusters.df, run_accession %in% new.isolates)
  new_isolates <- new_isolates %>% 
    mutate(in_range = ifelse(CollectionDate >= original.isolates.start & CollectionDate <= original.isolates.stop, "within", "outside"),
           start = original.isolates.start,
           stop = original.isolates.stop)
  return(new_isolates)
}

new_isolates <- missing.isolates.dates(1, SNPclusters.10.metadata, gcTree)

for (i in 2:147) {
  ni <- missing.isolates.dates(i, SNPclusters.10.metadata, gcTree)
  new_isolates <- bind_rows(new_isolates, ni)
}

new.clusters <- new.clusters %>% 
  add_count(new_cluster) %>% 
  arrange(desc(n)) %>% 
  distinct(run_accession, .keep_all = TRUE) %>% 
  dplyr::select(run_accession, new_cluster)

SNPclusters.10.updated <- left_join(SNPclusters.10, new.clusters)

SNPclusters.10.updated <- SNPclusters.10.updated %>% 
  mutate(final_cluster = ifelse(!is.na(new_cluster), new_cluster, cluster)) %>% 
  add_count(final_cluster) %>% 
  mutate(clustered = ifelse(n > 1, "clustered", "unclustered"))

SNPclusters.10.updated$final_cluster <- as.factor(as.character(SNPclusters.10.updated$final_cluster))