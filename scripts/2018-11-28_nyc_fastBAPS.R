library(fastbaps)
library(phytools)

sparse.data <- import_fasta_sparse_nt("../2018-11-19_nyc_pseudogenomes_noDups_snps.fasta", prior = "baps")
gubbins.tree <- phytools::read.newick("../gubbins_subset/2018-11-19_nyc_subset.final_tree.tre")
gubbins.rooted <- phytools::midpoint.root(gubbins.tree)
best.partition <- best_baps_partition(sparse.data, gubbins.rooted)
best.partition.df <- data.frame(id = gubbins.rooted$tip.label, fastbaps = best.partition, stringsAsFactors= FALSE)
write.table(best.partition.df, "2018-11-28_nyc_fastBAPS.txt", sep = "\t", row.names = F, quote = F)
