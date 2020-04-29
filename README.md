# GC_NYC
Scripts and data for Mortimer et al. 2020 "The distribution and spread of susceptible and resistant Neisseria gonorrhoeae across demographic groups in a major metropolitan center"

## data

### clustering

#### 2019-01-24_nyc_fastBAPS.txt

FastBAPS results for NYC samples

#### 2019-01-28_nyc_10SNP_clusters.txt

Clusters defined using a 10 non-recombinant SNP distance threshold

#### 2019-01-29_nyc_MRCA_clusters.txt

Clusters defined using all descendents of the MRCA of 10 SNP clusters

#### 2019-07-24_nyc_12SNP_clusters.txt

Clusters deined using a 12 non-recombinant SNP distance threshold

### distances

#### 2019-02-20_nyc_gubbins_distances.txt

Pairwise SNP distances from Gubbins masked alignment

### gubbins

Output files from Gubbins. branch_base_reconstructions.embl and .fasta..seq.joint were not included due to GitHub file size limits.

#### Mortimer_NYC_SupplementaryTable1.tsv

Metadata associated with sample including clinical information from electronic medical record and MICs for azithromycin, cefixime, ceftriaxone, and ciprofloxacin.

#### Mortimer_NYC_SupplementaryTable2.tsv

BAPS groups from a global collection of isolates originally collected for Ma and Mortimer et al. 2020.

#### major_clades.txt

Assignments of each NYC isolate to the two main gonococcal lineages.

## scripts

#### Figure1_lineages_clustering.R

Script to create figures 1B and 1C, including performing statistical test of association between patient sexual behavior and gonocococcal lineages.

#### Figure2_cluster_MIC.R

Script to plot MICs of isolates in the ten largest transmission clusters (Figure 4 in current version of paper).

#### Figure3_partner_cluster_MIC.R

Script to plot and perform statistical tests on differences in clustering and MICs across sexual behavior groups (Figure 2 in current version of paper).

#### Figure4_race_cluster_MIC.R

Script to plot and perform statistical tests on differences in clustering and MICs across races/ethnicities among heterosexual patients (Figure 3 in current version of paper).

#### SupplementaryFigure3_multiple_samples.R

Script to plot SNP distances between isolates from the same patient collected both at the same visit or at different visits.

#### mrca_clustering.R

Script to create MRCA based clusters.

#### phyloSignal_discreteTraits.R

Script to perform statistical tests related to phylogenetic signal of discrete traits (sexual behavior, race/ethnicity, HIV status).
