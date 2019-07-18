# Evenness
The data used in the code and in the paper by Chao and Ricotta (2019, Ecology) as a worked example were collected by Caccianiga et al. (2006) and analyzed in Ricotta et al. (2016, 2018). The full data set, available from Ricotta et al. (2016, their Appendix S1), contains abundances for a total of 45 Alpine species sampled in 59 vegetation plots each of 25 m^2 along a primary succession on the Rutor glacier (northern Italy). Based on the age of the glacial deposits, plots were assigned to three successional stages: early-succession (17 plots), mid-succession (32 plots) and late-succession (10 plots). We first compute species relative abundances within each plot, and then average these relative abundances across the plots within each stage. 

"Evenness" includes four files:

(1) Alpine data: the relative abundances of the 45 Alpine species for three successional stages; see Figure 1.

(2) Alpine phylo_tree: the phylogenetic tree of the 45 Alpine species; see Figure 1. This tree was taken from Ricotta et al. (2015, Appendix A).

(3) Evenness.R: Main code for computing the profiles of six classess of evenness measures (listed in Table 1 with plots in Figure 2) and for computing the contribution of each species/node to Jaccard-type taxonomic and phylogenetic dissimilarity measures (see Figures 3 and 4). The corresponding contribution to Sorensen-type taxonomic and phylogenetic dissimilarity measures are shown in Appendix S3.

(4) Lorenz.R: R code for computing and plotting the Lorenz curves of assemblages; see Appendix S4 for examples. 
