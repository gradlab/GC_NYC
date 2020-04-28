library(caper)
library(dplyr)
library(phytools)

gcTree <- read.newick("data/gubbins/2019-01-23_nyc_metadataComplete.final_tree.tre")
gcTree <- multi2di(gcTree)
gcTree$edge.length[gcTree$edge.length==0]<-max(nodeHeights(gcTree))*1e-6
samples <- gcTree$tip.label
metadata <- read.table("data/Mortimer_NYC_SupplementaryTable1.tsv", header = T, sep = "\t", comment.char = "")
metadata <- metadata %>% mutate(orientation = ifelse(partner == "MSM" | partner == "MSMW" | partner == "TSM", "MSM",
                                                                     ifelse(partner == "MSW" | partner == "WSM" | partner == "WSMW", "heterosexual", NA)))


metadata.orientation <- metadata %>% dplyr::select(Lane, orientation)
metadata.ct <- metadata %>% dplyr::select(Lane, CT_any)
metadata.hiv <- metadata %>% dplyr::select(Lane, hivstatus) %>% 
  mutate(hivstatus = as.character(hivstatus)) %>%
  mutate(hivstatus = ifelse(hivstatus == "" | hivstatus == "Unknown", NA, hivstatus))
metadata.race <- metadata %>% dplyr::select(Lane, Race) %>%
  mutate(black = ifelse(Race == "NH-Black", "yes", "no"))
metadata.syphilis <- metadata %>% dplyr::select(Lane, Hx_Syphilis) %>%
  mutate(Hx_Syphilis = ifelse(Hx_Syphilis == 9, NA, Hx_Syphilis))
  

gcTree.orientation <- comparative.data(phy = gcTree, data = metadata.orientation, names.col = Lane, na.omit = TRUE)
gcTree.ct <- comparative.data(phy = gcTree, data = metadata.ct, names.col = Lane, na.omit = TRUE)
gcTree.hiv <- comparative.data(phy = gcTree, data = metadata.hiv, names.col = Lane, na.omit = TRUE)
gcTree.race <- comparative.data(phy = gcTree, data = metadata.race, names.col = Lane, na.omit = TRUE)
gcTree.syphilis <- comparative.data(phy = gcTree, data = metadata.syphilis, names.col = Lane, na.omit = TRUE)
#Fritz and Purvis D

orientation.d <- phylo.d(data = gcTree.orientation, binvar = orientation, permut = 1000)
ct_any.d <- phylo.d(data = gcTree.ct, binvar = CT_any, permut = 1000)
hiv.d <- phylo.d(data = gcTree.hiv, binvar = hivstatus, permut = 1000)
race.d <- phylo.d(data = gcTree.race, binvar = black, permut = 1000)
syphilis.d <- phylo.d(data = gcTree.syphilis, binvar = Hx_Syphilis, permut = 1000)
