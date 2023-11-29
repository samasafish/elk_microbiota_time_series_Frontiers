## data from dada2 in R to phyloseq
## author: Sam Pannoni
## date: 2021 Dec...
########## Project Description:   ####
#  Time series elk samples from 4 elk (biological replicates), 
#  3 replicates/elk (technical replicates) each at 5 time points
#  Goal is to determine if time has a significant effect on the microbiome
#  
# Code description:  
# This phyloseq code "elk_16s....R" is intended to be run post "...dada2.R" script, 
# the latter performs raw data filtering and DADA2 variant calling and multi-run merging to format 
# and produce an ASV table, and Taxonomy table. 
# 
# This script performs visualization and statistical summaries of the 16s gene specifically
#
#Goals: test the COMMUNITY COMPOSITION difference between samples... (tricky definition)
#   we must assess community size (total abundances in each sample {sample: sum of species or total}) and 
#   shape separately (relative abundances of species between samples {sample: species/species total})
#
#    size effects = species richness (OVERALL ABUNDANCES of samples being compared, 
#   																sum of species row compared btwn samples) != usu. in microbiome studies
#   																<<<univariate permutation test of richness>>> DA test edgeR DESeq2
#    shape effects = composition differences between samples (RELATIVE ABUNDANCES of spp. within samples) 
#   																<<<multivariate permutaion test>>>  chi-square CA or CCA and adonis
#   	shape methods : (euclidean PCA on centered log-ratio (Aitchisons), or... chi-square dist. CCA {pure shape})
#     avoid the confounding nature of Bray-Curtis dissimilarity b/c its  a poor approximation of 
#     dissimilarity (not a distance since it violates symmetry etc.) and it confounds shape and size
#     
#     Example: a difference in total abundance between samples may cause a sig. diff. w. BC because it confounds size/shape
#     					using chi-square or Aitchisons will control for this by testing only shape (composition)
#     			
#     			if you think size is important, euclidean dist on CLR data  is size and shape, 
#     			as is chi-square on raw transformed abundances
#     			
#     use vegan::ANOSIM to test differences between beta diversity groups using appropriate distances.
#     
##### Statistical comparisons Notes ####
# quick notes:
#                         2 categories         >2 categories
# Parametric..............T-test...............ANOVA
# Non-parametric (rank).. Mann-Whitney U... ...Kruskal-Wallace
# Permutation ............Permutation t-test ..perMANOVA (Adonis)
# 
# Wilcoxon rank-sum test (also called Mann–Whitney test) is 
# the non-parametric two-sample T-test
# 
# ANOSIM used to compare within and between group similarity

# Specific taxa (apriori) can be tested with a chi-square test
# 
# Welche's T 

#plan:
# We will do an Adonis to test between group differences (Population, age, sex, bf)
# significance in this indicates that the communities are different
# 
# ANOSIM for groups of variable like sex and age

#Distances:
# Unifrac and weighted unifrac under value and overvalue taxa based on abundance
# Unweighted unifrac is unreliable and will not be used (dominated by abundant taxa and 
# sensitive to sampling depth/rarefaction)
# 
# weighted UniFrac can overweight differences between large proportional abundances 
# and underweight differences between small proportional abundances.
# 
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0161196
# GUniFrac helps to aleviate this imbalance but still is a combination of the two
# devtools::install_github("ChongWu-Biostat/MiSPU")
# library("MiSPU") #new pacakge in C
# GUniFrac(otu.tab, tree, alpha = 0.5)
# 
# information UniFrac and ratio UniFrac 
# that are not as sensitive to rarefaction and allow greater separation of outliers 
# than classic unweighted and weighted UniFrac. The goal is to address the 
# limitations of unweighted UniFrac's highly sensitive to rarefaction instance and 
# to sequencing depth in uniform data sets with no clear structure or separation 
# between groups.
# install.packages("UniFrac")
# library("UniFrac")Notes
########## load libs Short   ####
# remove(list = ls())
# install.packages("pacman")
package_list_all <- c("ANCOMBC", "graphics",
									# "microbiomeSeq", 
                  "utils", "vegan", "qwraps2", "breakaway",
                  "DECIPHER", "phyloseq", "DT", "microbiome",
                  "DESeq2", "metagenomeSeq", "HMP", "selbal",
                  "tidyverse", "reshape2", "reshape", "ggplot2",
                  "ggthemes", "cowplot", "Cairo", "viridis", "ALDEx2",
                  "magrittr", "dendextend", "rms", "picante", "DivNet",
                  "wesanderson", "metagMisc", "vsn", "hexbin")

package_list <- c("Biostrings", "ShortRead",
									"ANCOMBC", "vegan", "phyloseq", 
									"microbiome","metagMisc", "DESeq2",
									"breakaway", "DivNet", 
									"tidyverse", "dplyr", "reshape2", "reshape", "magrittr",
									"ggthemes", "cowplot", "Cairo", "viridis", "wesanderson")

pacman::p_load(package_list, character.only = T, update = F, install = F) 
pacman::p_load("Biostrings", "ShortRead", update = F)
# options(DT.options = list(
#   initComplete = JS("function(settings, json) {",
#                     "$(this.api().table().header()).css({'background-color': 
#   '#000', 'color': '#fff'});","}")))
########## load libs Long  ####
## if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# install.packages("remotes")
# BiocManager::install(c("RBGL","impute", "graph","pcaExplorer", "vsn", "phyloseqGraphTest", "igraph", "ggnetwork", "DECIPHER", "biomformat"))
# devtools::install_github(repo = "malucalle/selbal")
# devtools::install_github("adw96/DivNet", build_vignettes = TRUE, force = T)
# devtools::install_github("adw96/breakaway")
# BiocManager::install("WGCNA")
# install.packages("adespatial")
# devtools::install_github("umerijaz/microbiomeSeq", force = T) #error
# install.packages("eulerr")
# BiocManager::install("microbiome")
# devtools::install_github('microsud/microbiomeutilities')
# BiocManager::install("ALDEx2")
# install.packages(c("picante","HMP", "dendextend", "rms"))
# install.packages("qwraps2")
# require(graphics); require(utils)
# install.packages("devtools")
# install.packages("qwraps2")
# install.packages("devtools")
# install_github("tobiasgf/lulu")
# install.packages("phangorn")
# install.packages("V8")
# install.packages("mvabund") 
# install_github("joey711/phyloseq") 
# install_github("benjjneb/dada2")
# install_github("mikelove/DESeq2")
# install_github("vmikk/metagMisc")
# install_github("talgalili/d3heatmap")
# install_github("FrederickHuangLin/ANCOMBC")
# install.packages("qwraps2")
# install_github("tobiasgf/lulu")
# install.packages("phangorn")
# install.packages("V8")
# install.packages("mvabund")
# install.packages("plotly")
# install.packages("wesanderson")
# install_github("js229/Vennerable"); library(Vennerable);

########## source scripts ####
source("code/phyloseq_functions2.R"); #fast plotting and data manipulations


########## load data  -- skip if RDS files exist ####
cat("\n\nRead in the files for phyloseq....\n")
# tax_gg_genus <- read.table("output/gg_taxonomy.txt", header=T, sep = "\t")
# tax_rdp <-  read.table("output/rdp_taxonomy.txt", header=T, sep = "\t")
tax_rdp_spp <-  read.table("r_output/rdp16_spp_taxonomy.txt", header=T, sep = "\t")
counts_table <- read.table("r_output/ASVs_counts.txt", header=TRUE)
metadata <- read.csv("time_series_mappingFile_R.txt", header=TRUE,  sep="\t")
ASV_fasta_st <- readDNAStringSet("r_output/ASVs.fa")
gtr_tree <- read_tree("r_output/bootstrap_GTR_tree.tre")

rownames(metadata) <- metadata$Description
psmeta <- sample_data(metadata)

#transpose
tran.OTU <- as.matrix(t(counts_table))
rownames(tran.OTU)

# tax
class(tax_rdp_spp$Genus)
unique(tax_rdp_spp$Family)
tax_rdp_spp$Family <- as.character(tax_rdp_spp$Family)
tax_rdp_spp$Family[is.na(tax_rdp_spp$Family)] <- "f__Unknown"

tax_rdp_spp$Genus <- as.character(tax_rdp_spp$Genus)
tax_rdp_spp$Genus[is.na(tax_rdp_spp$Genus)] <- "g__Unknown"
unique(tax_rdp_spp$Genus)

tax_rdp_spp$Species <- as.character(tax_rdp_spp$Species)
tax_rdp_spp$Species[is.na(tax_rdp_spp$Species)] <- "s__Unknown"
unique(tax_rdp_spp$Species)

tax_rdp_m <- as.matrix(tax_rdp_spp)
# tax_rdp_m <- gsub("__$", "__Unknown", tax_rdp_m, perl = T)
# tax_rdp_m <- gsub("_$", "__Unknown", tax_rdp_m, perl = T)

tax_rdp.clean = as.matrix(tax_rdp_m[row.names(tax_rdp_m) %in% colnames(tran.OTU),])
dim(tax_rdp.clean)==dim(tax_rdp_m)
#remove extra species levels for now
# tax.clean.7 = tax.clean[,c(-8, -9)]
# tax_m <- gsub("_$", "", tax_m, perl = T)
all(rownames(tran.OTU) %in% metadata$Description)
# tax.clean$taxa_names = row.names(tax.clean)

tran.OTU = tran.OTU[order(row.names(tran.OTU)),]
metadata = metadata[order(row.names(metadata)),]
length(dimnames(tran.OTU))
length(dimnames(tax_rdp.clean))
cat("\n\nCombine into phyloseq object....\n")

OTU = otu_table(tran.OTU, taxa_are_rows = FALSE)
TAX.rdp = tax_table(tax_rdp.clean)
meta = sample_data(metadata)
gtr_tree.phy <- phy_tree(gtr_tree)
elk_rdp.ps <- phyloseq(OTU, meta, TAX.rdp, gtr_tree.phy)

# TAX.rdp.nas = tax_table(as.matrix(tax_rdp_spp))
# elk_rdp_FS.ps <- phyloseq(OTU, meta, TAX.rdp.nas)

set.seed(3313)
# elk_rdp.ps <- merge_phyloseq(elk_rdp.ps, phy_tree(gtr_tree)
phy_tree(elk_rdp.ps) <- ape::root(phy_tree(elk_rdp.ps), 
                                  sample(taxa_names(elk_rdp.ps), 1), 
                                  resolve.root = TRUE)
ape::is.rooted(phy_tree(elk_rdp.ps))

sample_data(elk_rdp.ps)
# plot_tree(elk_rdp.ps, color="Day", shape="Phylum", label.tips="Genus", plot.margin=0.6) 
# too big
# 
# Remove stuff that is redundant at this point
rm(ASV_fasta_st, OTU, TAX.rdp, counts_table, tax_rdp_m, tax_rdp_spp, 
	 tax_rdp.clean, tran.OTU, psmeta, meta, metadata, gtr_tree, gtr_tree.phy)

########## filter raw reads and ASV summaries    ####
   # install.packages("iNEXT")
   # require(metagMisc);require(iNEXT)
   # https://rdrr.io/github/vmikk/metagMisc/man/phyloseq_coverage.html
    
    length(get_taxa_unique( elk_rdp.ps,     "Genus"))#116
    length(get_taxa_unique( elk_rdp.ps,     "Species"))#8
    get_taxa_unique( elk_rdp.ps,     "Phylum") #NA and chloroplasts
    # make data-frame of nreads for both ASVs and Samples then plot
    plot.ASV_sum_log10 <- plot_nreads_summary(elk_rdp.ps, log10y = T) #log10y
    plot.ASV_sum <- plot_nreads_summary(elk_rdp.ps, log10y = F) #log10y
    #uncomment below to save as pdf
    # system2(command = "mkdir", args = "./r_output/visualizations")
    # pdf("r_output/visualizations/ASV_read_summary_log10.pdf", 
        # width=10,height=4);plot.ASV_sum_log10;dev.off()
    # pdf("plots/rbcL_plots/rbcL_ASV_summary.pdf", width=10,height=4);plot.ASV_sum;dev.off()
    
    # remove some garbage?
    sum(sample_sums(elk_rdp.ps) < 8000) #1 sample
    min((sample_sums(elk_rdp.ps))) #7812 minumum sample count
    # sum(taxa_sums(elk_rdp.ps) == 0) # 0 taxa with 0
    table(tax_table(elk_rdp.ps)[, "Kingdom"], exclude = NULL) #11 Archea
    table(tax_table(elk_rdp.ps)[, "Phylum"], exclude = NULL) #93 NA Phyla
    table(tax_table(elk_rdp.ps)[, "Class"], exclude = NULL) #404 NA Class, 2 chloroplast
    table(tax_table(elk_rdp.ps)[, "Order"], exclude = NULL) #493 NA Phyla, 2 chloroplast
    sum(taxa_sums(elk_rdp.ps) < 1000) #4839 taxa with less than x total observations
    sum(taxa_sums(elk_rdp.ps) < 10)   #855 taxa with less than x total observations
    sum(taxa_sums(elk_rdp.ps) < 2)    #95 doubletons
    ps0 <- subset_taxa(elk_rdp.ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
    ps0 <- subset_taxa(ps0, !Class %in% c("Chloroplast"))
    ps0 <- subset_taxa(ps0, !Kingdom %in% c("Archaea"))
    sum(taxa_sums(ps0) < 10) # 835
    sum(taxa_sums(ps0) < 2) #94
    sum(taxa_sums(ps0) < 1000) #4736
    table(tax_table(ps0)[, "Kingdom"], exclude = NULL) #Bacteria only
    table(tax_table(ps0)[, "Phylum"], exclude = NULL) #0 NA Phyla
    table(tax_table(ps0)[, "Class"], exclude = NULL) #311 NA Class, 0 chloroplast
    sum(sample_sums(ps0) < 10000) #2 samples
    min(sample_sums(ps0)) #7799
    sum(taxa_sums(elk_rdp.ps) == 1) #95 taxa with 1
    
    ### revision from Frontiers Reviewer #1 ####
    
    sum(get_taxa_unique( elk_rdp.ps, "Family") == "Mitochondria")
    sum(get_taxa_unique( elk_rdp.ps, "Family") == "mitochondria")
    sum(get_taxa_unique( elk_rdp.ps, "Genus") == "Mitochondria")
    sum(get_taxa_unique( elk_rdp.ps, "Genus") == "mitochondria")
    get_taxa_unique( elk_rdp.ps, "Class") == "Alphaproteobacteria"
    ps0 <- subset_taxa(elk_rdp.ps, Class %in% "Alphaproteobacteria")
    table(tax_table(ps0)[, "Genus"], exclude = NULL)
    rm(ps0)
    #####
    
    
    # filter done for Archaea, chloroplasts and Phylum NAs
    # 
    #reset the ps object
    elk_rdp.ps <- ps0
    rm(ps0)
    
    plot.ASV_sum_log10_postFilt <- plot_nreads_summary(elk_rdp.ps, log10y = T) #log10y
    # plot.ASV_sum <- plot_nreads_summary(elk_rdp.ps, log10y = F)
    
      #### summarise measurement and counts per sample ####
	  sum(sample_data(elk_rdp.ps)$Elk == "1") #15 entries or 3 reps/day/elk
    sum(sample_data(elk_rdp.ps)$Elk == "2") #15
    sum(sample_data(elk_rdp.ps)$Elk == "3") #15
    sum(sample_data(elk_rdp.ps)$Elk == "4") #15
    length(unique(sample_data(elk_rdp.ps)$Elk)) # total 4 ELK
    length(unique(sample_data(elk_rdp.ps)$Day)) # 5 different Days
    length(unique(sample_data(elk_rdp.ps)$Rep)) # 3 experimental replicates
    
    sample_otu_dframe <- as.data.frame(t(otu_table(elk_rdp.ps)))
    asv_sums_by_sample <- colSums(sample_otu_dframe)
    summary(asv_sums_by_sample) #reads/ASVs by sample
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 7799   19849   27936   35772   42160  126043 
    summary(apply(sample_otu_dframe > 1, 2, sum)) #ASV presence by sample
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 301.0   606.2   776.5   849.9  1033.5  1965.0 
    # 
    #top 10 ASVs and their taxonomy
    top10_taxa = names(sort(taxa_sums(elk_rdp.ps), TRUE)[1:10])
    top10.ps = prune_taxa(top10_taxa, elk_rdp.ps)
    get_taxa_unique(top10.ps, taxonomic.rank = "Phylum")
    # "Firmicutes"    "Bacteroidetes"
    get_taxa_unique(top10.ps, taxonomic.rank = "Genus")
    # "Sporobacter" "g__Unknown"  "Prevotella"  "Bacteroides"
    
    # saveRDS(object = elk_rdp.ps, file = "r_output/intermediate_rds_files/elk_rdp.ps.rds")

########## make data CoDa compostitional ####
# BiocManager::install("philr") #comp phylo dist
# browseVignettes("philr")
    BiocManager::install("ALDEx2")
    devtools::install_github("tpq/propr")
    install.packages("cli", force=T)
    devtools::install_github('ggloor/CoDaSeq/CoDaSeq', dependencies = F)
    # The elements of a composition are called components or parts. In a composition the 
    # value of each component is not informative by itself and the relevant information is 
    # contained in the ratios between the components or parts.
    # Microbiome data is compositional because the information that abundance tables contain 
    # is relative. In a microbiome abundance table, the total number of counts per sample is 
    # highly variable and constrained by the maximum number of DNA reads that the sequencer 
    # can provide. This total count constraint induces strong dependencies among the 
    # abundances of the different taxa; an increase in the abundance of one taxon implies the 
    # decrease of the observed number of counts for some of the other taxa so that the total 
    # number of counts does not exceed the specified sequencing depth. Moreover, observed raw 
    # abundances and the total number of reads per sample are non-informative since they 
    # represent only a fraction or random sample of the original DNA content in the 
    # environment. These characteristics of microbiome abundance data clearly fall into the 
    # notion of compositional data.
    # 
    # the Aitchison distance, which is simply the Euclidean distance between samples after 
    # clr transformation, and the distances between samples are the same as the phylogenetic ilr
    # 
    # The replacement for β-diversity exploration of microbiome data is the variance-based 
    # compositional principal component (PCA) biplot 
    
    # a couple options. CLR cannot be subset since the denominator uses geometric mean which will change 
    # if features are removed. iqlr inter-quartile range log ratio (denom is invariant features) is a 
    # better metric for this case and is superior to alr when a spike-in abundance is unavailable.
    # https://cran.r-project.org/web/packages/propr/vignettes/g_questions.html
    
    #idea: modeling: take iqlr of time zero elk, use the denom to scale all CLR in the remaining days
    #								this didn't work no features were the same between elk!!!
    #							take the clr of time zero elk...
    #							worked.
    #			this would be a useful way to estimate change across samples according to time 0 
    #			while preserving comparability
    #			WARNING:: this is not CoDa anymore but a normalization. But this also means that you can subset
    #			
    #			Clustering: use CLR to cluster since its representative across all samples
    
    allElk_day_0 <- subset_samples(elk_rdp.ps, Day == "0")
    elk_rdp_genus.ps <- tax_glom(elk_rdp.ps, "Genus", NArm = F)
    elk_rdp_species.ps <- tax_glom(elk_rdp.ps, "Species", NArm = F)
    
    x_iqlr_spp <- ALDEx2::aldex.clr(t(elk_rdp_species.ps@otu_table), elk_rdp_species.ps@sam_data$Elk, denom = "iqlr", mc.samples = 128, useMC = T)
    x_iqlr_genus <- ALDEx2::aldex.clr(t(elk_rdp_genus.ps@otu_table), elk_rdp_genus.ps@sam_data$Elk, denom = "iqlr", mc.samples = 128, useMC = T)
    
    x_clr_spp <- ALDEx2::aldex.clr(t(elk_rdp_species.ps@otu_table), elk_rdp_species.ps@sam_data$Elk, denom = "all", mc.samples = 128, useMC = T)
    x_clr_genus <- ALDEx2::aldex.clr(t(elk_rdp_genus.ps@otu_table), elk_rdp_genus.ps@sam_data$Elk, denom = "all", mc.samples = 128, useMC = T)
   
    x_iqlr <- ALDEx2::aldex.clr(t(elk_rdp.ps@otu_table), elk_rdp.ps@sam_data$Elk, denom = "iqlr", mc.samples = 128, useMC = T)
    x_clr <- ALDEx2::aldex.clr(t(elk_rdp.ps@otu_table), elk_rdp.ps@sam_data$Elk, denom = "all", mc.samples = 128, useMC = T)
    
    #calculate Day0 CLR centered (day0_CLR)
    # x_iqlr_day0 <- ALDEx2::aldex.clr(t(allElk_day_0@otu_table), allElk_day_0@sam_data$Elk, denom = "iqlr", mc.samples = 128, useMC = T)
    #no features intersect! yikes
    x_clr_day0 <- ALDEx2::aldex.clr(t(allElk_day_0@otu_table), allElk_day_0@sam_data$Elk, denom = "all", mc.samples = 128, useMC = T, verbose = T)
    # the rows to calculate denom
    length(x_clr_day0@denom) #3523 taxa to compute day0 geom-mean
    x_clr_day0@reads$ETS_011[x_clr_day0@denom,]
    
    x_day0_CLR <- ALDEx2::aldex.clr(t(elk_rdp.ps@otu_table), elk_rdp.ps@sam_data$Elk, denom = x_clr_day0@denom, mc.samples = 128, useMC = T)
    
    #convert to table for rest of analysis
    iqlr_taxa       <- propr::aldex2propr(x_iqlr, how = "rho")
    genus_iqlr_taxa <- propr::aldex2propr(x_iqlr_genus, how = "rho")
    spp_iqlr_taxa   <- propr::aldex2propr(x_iqlr_spp, how = "rho")
    
    clr_taxa        <- propr::aldex2propr(x_clr, how = "rho")
    day0clr_taxa    <- propr::aldex2propr(x_day0_CLR, how = "rho")
    genus_clr_taxa  <- propr::aldex2propr(x_clr_genus, how = "rho")
    spp_clr_taxa    <- propr::aldex2propr(x_clr_spp, how = "rho")
    
    #convert to phylsoeq objects
    elk_rdp_ASV_iQlr.ps   <- elk_rdp.ps
    elk_rdp_genus_iQlr.ps <- elk_rdp.ps
    elk_rdp_spp_iQlr.ps   <- elk_rdp.ps
    
    elk_rdp_ASV_CLR.ps   <- elk_rdp.ps
    elk_rdp_day0_CLR.ps  <- elk_rdp.ps
    elk_rdp_genus_CLR.ps <- elk_rdp.ps
    elk_rdp_spp_CLR.ps   <- elk_rdp.ps
   
    otu_table(elk_rdp_ASV_iQlr.ps)   <- otu_table(iqlr_taxa@logratio, taxa_are_rows = F)
    otu_table(elk_rdp_genus_iQlr.ps) <- otu_table(genus_iqlr_taxa@logratio, taxa_are_rows = F)
    otu_table(elk_rdp_spp_iQlr.ps)   <- otu_table(spp_iqlr_taxa@logratio, taxa_are_rows = F)
    
    otu_table(elk_rdp_ASV_CLR.ps)    <- otu_table(clr_taxa@logratio, taxa_are_rows = F)
    otu_table(elk_rdp_day0_CLR.ps)   <- otu_table(day0clr_taxa@logratio, taxa_are_rows = F)
    otu_table(elk_rdp_genus_CLR.ps)  <- otu_table(genus_clr_taxa@logratio, taxa_are_rows = F)
    otu_table(elk_rdp_spp_CLR.ps)    <- otu_table(spp_clr_taxa@logratio, taxa_are_rows = F)
    
    # Save for midpoint reloading
    saveRDS(object = elk_rdp_ASV_iQlr.ps,   file = "r_output/intermediate_rds_files/elk_rdp_ASV_iQlr.ps.rds")
    saveRDS(object = elk_rdp_genus_iQlr.ps, file = "r_output/intermediate_rds_files/elk_rdp_genus_iQlr.ps.rds")
    saveRDS(object = elk_rdp_spp_iQlr.ps,   file = "r_output/intermediate_rds_files/elk_rdp_spp_iQlr.ps.rds")
    
    saveRDS(object = elk_rdp_ASV_CLR.ps,    file = "r_output/intermediate_rds_files/elk_rdp_ASV_CLR.ps.rds")
    saveRDS(object = elk_rdp_day0_CLR.ps,   file = "r_output/intermediate_rds_files/elk_rdp_day0_CLR.ps.rds")
    saveRDS(object = elk_rdp_genus_CLR.ps,  file = "r_output/intermediate_rds_files/elk_rdp_genus_CLR.ps.rds")
    saveRDS(object = elk_rdp_spp_CLR.ps,    file = "r_output/intermediate_rds_files/elk_rdp_spp_CLR.ps.rds")
    
    
    plot_nreads_summary(elk_rdp_day0_CLR.ps,log10y = F)
    plot_nreads_summary(elk_rdp_ASV_CLR.ps, log10y = F)
    
########## normalize [X] ANCOM-BC Bias-adjusted abundances ####
    ##### [X]ANCOMBC on Genus glom ######
    genus_data.ps = tax_glom(elk_rdp.ps, "Genus")
    otu_table(genus_data.ps)
    out_elk_genus_anc <- ancombc(genus_data.ps, 
                                 formula = "Elk", 
                                 p_adj_method = "holm", zero_cut = 0.99, lib_cut = 0, 
                                 group = "Elk", struc_zero = TRUE, 
                                 neg_lb = TRUE, tol = 1e-5, 
                                 max_iter = 100, conserve = TRUE, alpha = 0.05, 
                                 global = TRUE)
    #adjust abundances #
    samp_frac = out_elk_genus_anc$samp_frac
    # Replace NA with 0
    sum(is.na(samp_frac)) # 0 
    samp_frac[is.na(samp_frac)] = 0 
    
    round(abundances(genus_data.ps)[, 1:6], 5) %>% 
    	datatable(caption = "observed abundances")
    # Add pesudo-count (1) to avoid taking the log of 0
    log_obs_abn = log(abundances(genus_data.ps) + 1) 
    # Adjust the log observed abundances
    log_obs_abn_adj = t(t(log_obs_abn) - samp_frac)
    # Show the first 6 samples
    round(log_obs_abn_adj[, 1:6], 2) %>% 
    	datatable(caption = "Bias-adjusted log observed abundances")
    
    elk_rdp_adujsted.ps <- genus_data.ps
    otu_table(elk_rdp_adujsted.ps) <- otu_table(log_obs_abn_adj, taxa_are_rows = T)
    
    ##### [X]repeat with ASVs #####
    out_elk_asv_anc <- ancombc(elk_rdp.ps, 
    														formula = "Elk", 
    														p_adj_method = "holm", zero_cut = 0.999, 
    													  lib_cut = 0, 
    														group = "Elk", struc_zero = F, 
    														neg_lb = F, tol = 1e-5, 
    														max_iter = 100, conserve = F, alpha = 0.05, 
    														global = TRUE)
    #adjust abundances #
    samp_frac = out_elk_asv_anc$samp_frac
    # Replace NA with 0
    sum(is.na(samp_frac)) # 0 
    # samp_frac[is.na(samp_frac)] = 0 
    
    round(abundances(elk_rdp.ps)[, 1:6], 5) %>% 
      DT::datatable(caption = "observed abundances")
    # Add pesudo-count (1) to avoid taking the log of 0
    log_obs_abn = log(abundances(elk_rdp.ps) + 1) 
    # Adjust the log observed abundances
    log_obs_abn_adj = t(t(log_obs_abn) - samp_frac)
    # Show the first 6 samples
    round(log_obs_abn_adj[, 1:6], 2) %>% 
      DT::datatable(caption = "Bias-adjusted log observed abundances")
    
    #make a new phyloseq obj
    elk_rdp_ASV_adujsted.ps <- elk_rdp.ps
    otu_table(elk_rdp_ASV_adujsted.ps) <- otu_table(log_obs_abn_adj, taxa_are_rows = T)
    
    plot_nreads_summary(elk_rdp_ASV_adujsted.ps, log10y = F)
    plot_nreads_summary(elk_rdp.ps, log10y = T)
    plot_taxa_prevalence(elk_rdp_ASV_adujsted.ps, level = "Phylum", detection = 1)
    # plot_taxa_prevalence(elk_rdp.ps, level = "Phylum", detection = 1)

    #top 10 ASVs and their taxonomy
    top10_taxa_adj = names(sort(taxa_sums(elk_rdp_ASV_adujsted.ps), TRUE)[1:10])
    top10_adj.ps = prune_taxa(top10_taxa_adj, elk_rdp_ASV_adujsted.ps)

    get_taxa_unique(top10_adj.ps, taxonomic.rank = "Phylum")
    # [1] "Firmicutes"    "Bacteroidetes"
    
    get_taxa_unique(top10_adj.ps, taxonomic.rank = "Genus")
    # [1] "Sporobacter"    "g__Unknown"     "Bacteroides"    "Clostridium_IV"

    saveRDS(object = elk_rdp_ASV_adujsted.ps, file = "r_output/intermediate_rds_files/elk_rdp_ASV_ANCOM_adjusted.ps.rds")

    ##### VST DESeq2 ####
    countData <- as.matrix(t(elk_rdp.ps@otu_table))
    class(countData) <- NULL
    sam_data <- as.data.frame(elk_rdp.ps@sam_data)
    sam_data$Elk <- factor(as.character(sam_data$Elk))
    sam_data$Day <- factor(as.character(sam_data$Day))
    sam_data$Day <- scale(as.numeric(as.character(sam_data$Day)))
    
    # sam_data$Day <- as.factor(sam_data$Day)
    deseq_counts <- DESeq2::DESeqDataSetFromMatrix(countData = countData, 
    																			 colData = sam_data, 
    																			 design = ~Day*Elk) 
    deseq_counts_vst <- DESeq2::varianceStabilizingTransformation(deseq_counts)
    vst_trans_count_tab <- SummarizedExperiment::assay(deseq_counts_vst)
    
    # making our phyloseq object with transformed table
    elk_rdp_VST.ps <- elk_rdp.ps
    otu_table(elk_rdp_VST.ps) <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
    
    saveRDS(object = elk_rdp_VST.ps, file = "r_output/intermediate_rds_files/elk_rdp_VST.ps.rds")
    
######### load RDS [optional]####
    #normalizations (ok to subset)
    elk_rdp_adujsted.ps 		<-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_ANCOM_adjusted.ps.rds")
    elk_rdp_ASV_adujsted.ps <-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_ASV_ANCOM_adjusted.ps.rds")
    elk_rdp.ps							<-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp.ps.rds")
    elk_rdp_VST.ps					<-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_VST.ps.rds")
    
    elk_rdp_ASV_iQlr.ps   <-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_ASV_iQlr.ps.rds")
    elk_rdp_genus_iQlr.ps <-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_genus_iQlr.ps.rds")
    elk_rdp_spp_iQlr.ps   <-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_spp_iQlr.ps.rds")
    
    #CoDa conforming CLR (dont subset further)
    elk_rdp_ASV_CLR.ps    <-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_ASV_CLR.ps.rds")
    elk_rdp_day0_CLR.ps   <-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_day0_CLR.ps.rds")
    elk_rdp_genus_CLR.ps  <-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_genus_CLR.ps.rds")
    elk_rdp_spp_CLR.ps    <-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_spp_CLR.ps.rds")
    																	
    
#### clustering differences ####
    library(dendextend)
    library(pvclust)
    library(corrplot)
    # install.packages('circlize')
    # get data to compare
    vst_otu_table <- t(elk_rdp_VST.ps@otu_table)
    row.names(vst_otu_table) <- paste0("E",elk_rdp_VST.ps@sam_data$Elk, 
    																	 "_Day.", 
    																	 elk_rdp_VST.ps@sam_data$Day,
    																	 ".", elk_rdp_VST.ps@sam_data$Rep)
		
  	ancom_otu_table <- t(elk_rdp_ASV_adujsted.ps@otu_table)
  	row.names(ancom_otu_table) <- paste0("E",elk_rdp_ASV_adujsted.ps@sam_data$Elk, "_Day.", 
  																			 elk_rdp_ASV_adujsted.ps@sam_data$Day, ".", elk_rdp_ASV_adujsted.ps@sam_data$Rep)
  	
  	comp_otu_table <- otu_table(microbiome::transform(elk_rdp.ps, transform = "compositional"))
  	row.names(comp_otu_table) <- paste0("E",elk_rdp.ps@sam_data$Elk, "_Day.", 
  																			elk_rdp.ps@sam_data$Day, ".", elk_rdp.ps@sam_data$Rep)
  	
  	CLR_otu_table <- otu_table(elk_rdp_spp_CLR.ps)
  	row.names(CLR_otu_table) <- paste0("E",elk_rdp_spp_CLR.ps@sam_data$Elk, "_Day.", 
  																		 elk_rdp_spp_CLR.ps@sam_data$Day, ".", elk_rdp_spp_CLR.ps@sam_data$Rep)
  	
  	
  	
  	##### which cluster algorithm is best ####
    d_vst      <- dist(vst_otu_table)
    d_ancom    <- dist(ancom_otu_table)
    d_spp_CLR  <- dist(CLR_otu_table)
    d_spp_iQlr <- dist(elk_rdp_spp_iQlr.ps@otu_table)
    
    hclust_methods <- c("ward.D", "single", "complete", "average", "mcquitty", 
    										"median", "centroid", "ward.D2")
    elk_vst_dendlist <- dendlist()
    elk_ancom_dendlist <- dendlist()
    elk_spp_CLR_dendlist <- dendlist()
    
	    for(i in seq_along(hclust_methods)) {
    	hc_vst <- hclust(d_vst, method = hclust_methods[i])   
    	elk_vst_dendlist <- dendlist(elk_vst_dendlist, as.dendrogram(hc_vst))
    }
    for(i in seq_along(hclust_methods)) {
    	hc_vst <- hclust(d_ancom, method = hclust_methods[i])   
    	elk_ancom_dendlist <- dendlist(elk_ancom_dendlist, as.dendrogram(hc_vst))
    }
    for(i in seq_along(hclust_methods)) {
    	hc_vst <- hclust(d_spp_CLR, method = hclust_methods[i])   
    	elk_spp_CLR_dendlist <- dendlist(elk_spp_CLR_dendlist, as.dendrogram(hc_vst))
    }
    names(elk_vst_dendlist)   <- hclust_methods
    names(elk_ancom_dendlist) <- hclust_methods
    names(elk_spp_CLR_dendlist) <- hclust_methods
    # elk_vst_dendlist
    
    vst_dendlist_cor <- cor.dendlist(elk_vst_dendlist)
    corrplot::corrplot(vst_dendlist_cor, "pie", "lower")
    vst_dendlist_cor_spearman <- cor.dendlist(elk_vst_dendlist, method_coef = "spearman")
    corrplot::corrplot(vst_dendlist_cor_spearman, "pie", "lower")
    
    ancom_dendlist_cor <- cor.dendlist(elk_ancom_dendlist)
    corrplot::corrplot(ancom_dendlist_cor, "pie", "lower")
    ancom_dendlist_cor_spearman <- cor.dendlist(elk_ancom_dendlist, method_coef = "spearman")
    corrplot::corrplot(ancom_dendlist_cor_spearman, "pie", "lower")
    
    dev.off()
    spp_CLR_dendlist_cor <- cor.dendlist(elk_spp_CLR_dendlist)
    corrplot::corrplot(spp_CLR_dendlist_cor, "pie", "lower")
    spp_CLR_dendlist_cor_spearman <- cor.dendlist(elk_spp_CLR_dendlist, method_coef = "spearman")
    corrplot::corrplot(spp_CLR_dendlist_cor_spearman, "pie", "lower")
    
    
    ##### ward.D2 is best ####
   
    # for some reason pvclust requires a transposed matrix compared to dendogram, if not it will create an ASV cluster (bad)
    
    # use euclidean for Aitchisons dist on CLR
    vst_result_10k   <- pvclust(t(vst_otu_table), method.dist = "euclidean", 
    									    method.hclust="ward.D2", nboot=10000, parallel = T)
    
    ancom_result_10k <- pvclust(t(ancom_otu_table), method.dist = "euclidean", 
    											method.hclust="ward.D2", nboot=10000, parallel = T)
    
    comp_result_10k  <- pvclust(t(comp_otu_table), method.dist = "euclidean", 
    											method.hclust="ward.D2", nboot=10000, parallel = T)
    
    sppCLR_result_10k  <- pvclust(t(CLR_otu_table), method.dist = "euclidean", 
    														method.hclust="ward.D2", nboot=10000, parallel = T)
    
    spp_iQlr_result_10k  <- pvclust(t(elk_rdp_spp_iQlr.ps@otu_table), method.dist = "euclidean", 
    															method.hclust="ward.D2", nboot=10000, parallel = T)
    
    genusCLR_result_10k  <- pvclust(t(elk_rdp_genus_CLR.ps@otu_table), method.dist = "euclidean", 
    															method.hclust="ward.D2", nboot=10000, parallel = T)
    
    genus_iQlr_result_10k  <- pvclust(t(elk_rdp_genus_iQlr.ps@otu_table), method.dist = "euclidean", 
    																method.hclust="ward.D2", nboot=10000, parallel = T)
    
    day0_CLR_result_10k  <- pvclust(t(elk_rdp_day0_CLR.ps@otu_table), method.dist = "euclidean", 
    																	method.hclust="ward.D2", nboot=10000, parallel = T)
    
    ASVCLR_result_10k <- pvclust(t(elk_rdp_ASV_CLR.ps@otu_table), method.dist = "euclidean", 
    														 method.hclust="ward.D2", nboot=10000, parallel = T)
    
    saveRDS(ASVCLR_result_10k, file = "r_output/intermediate_rds_files/asvCLR_hier_10K.rds")
    saveRDS(vst_result_10k, file = "r_output/intermediate_rds_files/vst_hier_10K.rds")
    saveRDS(ancom_result_10k, file = "r_output/intermediate_rds_files/ancom_hier_10K.rds")
    saveRDS(comp_result_10k, file = "r_output/intermediate_rds_files/comp_hier_10K.rds")
    saveRDS(sppCLR_result_10k, file = "r_output/intermediate_rds_files/sppCLR_hier_10K.rds")
    saveRDS(spp_iQlr_result_10k, file = "r_output/intermediate_rds_files/spp_iQlr_hier_10K.rds")
    saveRDS(genusCLR_result_10k, file = "r_output/intermediate_rds_files/genusCLR_hier_10K.rds")
    saveRDS(genus_iQlr_result_10k, file = "r_output/intermediate_rds_files/genus_iQlr_hier_10K.rds")
    saveRDS(day0_CLR_result_10k, file = "r_output/intermediate_rds_files/day0_CLR_hier_10K.rds")
    
    ASVCLR_result_10k     <- readRDS(file = "r_output/intermediate_rds_files/asvCLR_hier_10K.rds")
    vst_result_10k        <- readRDS(file = "r_output/intermediate_rds_files/vst_hier_10K.rds")
    ancom_result_10k      <- readRDS(file = "r_output/intermediate_rds_files/ancom_hier_10K.rds")
    comp_result_10k       <- readRDS(file = "r_output/intermediate_rds_files/comp_hier_10K.rds")
    sppCLR_result_10k     <- readRDS(file = "r_output/intermediate_rds_files/sppCLR_hier_10K.rds")
    spp_iQlr_result_10k   <- readRDS(file = "r_output/intermediate_rds_files/spp_iQlr_hier_10K.rds")
    genusCLR_result_10k   <- readRDS(file = "r_output/intermediate_rds_files/genusCLR_hier_10K.rds")
    genus_iQlr_result_10k <- readRDS(file = "r_output/intermediate_rds_files/genus_iQlr_hier_10K.rds")
    day0_CLR_result_10k   <- readRDS(file = "r_output/intermediate_rds_files/day0_CLR_hier_10K.rds")
    # pdf("r_output/visualizations/beta_div_plots/hier_tree_boot10k_VST.pdf",         # File name
    # 		width = 14, height = 7, useDingbats = T, compress = F)
    # par(mar = c(2,3,3,2) + 0.1)
    # plot(vst_result_10k)
    # # pvrect2(vst_result_10k, alpha=0.95, pv = "au", max.only = F, lower_rect = 50)
    # dev.off()
    # 
    # pdf("r_output/visualizations/beta_div_plots/hier_tree_boot10k_ANCOM.pdf",         # File name
    # 		width = 14, height = 7, useDingbats = T, compress = F)
    # plot(ancom_result_10k)
    # # pvrect2(ancom_result_10k, alpha=0.95, pv = "au", max.only = F, lower_rect = 50)
    # dev.off()
    # 
    # plot(comp_result_10k)
    # pvrect2(comp_result_10k, alpha=0.95, pv = "au", max.only = F, xpd = t, lower_rect = -0.035)
    
    plot(genusCLR_result_10k)
    plot(sppCLR_result_10k)
    plot(genus_iQlr_result_10k)
    plot(spp_iQlr_result_10k)
    plot(day0_CLR_result_10k)
    plot(ASVCLR_result_10k)
    pvrect2(ASVCLR_result_10k, alpha=0.95, pv = "au", max.only = F, xpd = t, lower_rect = -0.035)
    
    vst_result_10k.dend      <- vst_result_10k %>% as.dendrogram 
    ancom_result_10k.dend    <- ancom_result_10k %>% as.dendrogram 
    comp_result_10k.dend     <- comp_result_10k %>% as.dendrogram 
    
    ASVCLR_result_10k.dend     <- ASVCLR_result_10k %>% as.dendrogram
    genusCLR_result_10k.dend   <- genusCLR_result_10k %>% as.dendrogram 
    sppCLR_result_10k.dend     <- sppCLR_result_10k %>% as.dendrogram 
    genus_iQlr_result_10k.dend <- genus_iQlr_result_10k %>% as.dendrogram 
    spp_iQlr_result_10k.dend   <- spp_iQlr_result_10k %>% as.dendrogram 
    day0_CLR_result_10k.dend   <- day0_CLR_result_10k %>% as.dendrogram 
    
    # make plots nicer
    dendro_pretty <- function(dend) {
  						require(viridis)
					    day_labels <- as.character(paste0("Day.", elk_rdp_genus_CLR.ps@sam_data$Day[order.dendrogram(dend)],
					    																	".", elk_rdp_genus_CLR.ps@sam_data$Rep[order.dendrogram(dend)],
					    										 "--Elk.", elk_rdp_genus_CLR.ps@sam_data$Elk[order.dendrogram(dend)]))
					    day_labels <- gsub("14--", "14-", day_labels)
					    # day_labels <- factor(day_labels)
					    
					    day_pch_labels <- paste0("Day.", elk_rdp_genus_CLR.ps@sam_data$Day[order.dendrogram(dend)])
					    day_pch_labels <- factor(day_pch_labels)
					    levels(day_pch_labels)[levels(day_pch_labels) ==  "Day.0"] <- 15#23 #21, 15, 22, 25
					    levels(day_pch_labels)[levels(day_pch_labels) ==  "Day.1"] <- 15#21
					    levels(day_pch_labels)[levels(day_pch_labels) ==  "Day.3"] <- 15#24
					    levels(day_pch_labels)[levels(day_pch_labels) ==  "Day.7"] <- 15#22
					    levels(day_pch_labels)[levels(day_pch_labels) == "Day.14"] <- 15#25
					    day_pch_labels <- as.numeric(as.character(day_pch_labels))
					    
					    day_color_labels <- paste0("Day.", elk_rdp_genus_CLR.ps@sam_data$Day[order.dendrogram(dend)])
					    day_color_labels <- factor(day_color_labels)
					    levels(day_color_labels)[levels(day_color_labels) == "Day.0"]  <- viridis(6)[2]
					    levels(day_color_labels)[levels(day_color_labels) == "Day.1"]  <- viridis(6)[3]
					    levels(day_color_labels)[levels(day_color_labels) == "Day.3"]  <- viridis(6)[4]
					    levels(day_color_labels)[levels(day_color_labels) == "Day.7"]  <- viridis(6)[5]
					    levels(day_color_labels)[levels(day_color_labels) == "Day.14"] <- viridis(6)[6]
					    
					    label_obj <- data.frame(day_labels       = day_labels,
					    												day_pch_labels   = day_pch_labels,
					    												day_color_labels = day_color_labels,
					    												stringsAsFactors = F)
					    label_obj$day_labels <- day_labels
					    label_obj$day_pch_labels <- day_pch_labels
					    label_obj$day_color_labels <- day_color_labels
					    return(label_obj)
    }
    
    asvCLR_labs <-dendro_pretty(ASVCLR_result_10k.dend)
    day0_labs <- dendro_pretty(day0_CLR_result_10k.dend)
    genusCLR_labs <- dendro_pretty(genusCLR_result_10k.dend)
    sppCLR_labs <- dendro_pretty(sppCLR_result_10k.dend)
    genusiQlr_labs <- dendro_pretty(genus_iQlr_result_10k.dend)
    sspiQlr_labs <- dendro_pretty(spp_iQlr_result_10k.dend)
    # genusCLR_result_10k.dend  
    # sppCLR_result_10k.dend    
    # genus_iQlr_result_10k.dend
    # spp_iQlr_result_10k.dend  
    
    par(mar = c(4,2,2,8) + 0.1)
    
    ASVCLR_result_10k.dend <- ASVCLR_result_10k.dend %>%
										    	dendextend::set("branches_k_color", k=4) %>%
										    	# set("hang_leaves", hang_height=80) %>% 
										    	dendextend::set("leaves_pch", asvCLR_labs$day_pch_labels) %>%  # node point type
										    	dendextend::set("leaves_cex", 1) %>%  # node point size
										    	dendextend::set("leaves_col", asvCLR_labs$day_color_labels) %>% # node point color
										    	dendextend::set("labels_cex", 0.8) %>%
										    	dendextend::set("labels",asvCLR_labs$day_labels) #%>%
										     # plot(horiz = TRUE)
   
    
    genusCLR_result_10k.dend <- genusCLR_result_10k.dend %>%
    	dendextend::set("branches_k_color", k=4) %>%
    	# set("hang_leaves", hang_height=80) %>% 
    	dendextend::set("leaves_pch", genusCLR_labs$day_pch_labels) %>%  # node point type
    	dendextend::set("leaves_cex", 1) %>%  # node point size
    	dendextend::set("leaves_col", genusCLR_labs$day_color_labels) %>% # node point color
    	dendextend::set("labels_cex", 0.8) %>%
    	dendextend::set("labels",   genusCLR_labs$day_labels)
    
    sppCLR_result_10k.dend <- sppCLR_result_10k.dend %>%
    	dendextend::set("branches_k_color", k=4) %>%
    	# set("hang_leaves", hang_height=80) %>% 
    	dendextend::set("leaves_pch", sppCLR_labs$day_pch_labels) %>%  # node point type
    	dendextend::set("leaves_cex", 1) %>%  # node point size
    	dendextend::set("leaves_col", sppCLR_labs$day_color_labels) %>% # node point color
    	dendextend::set("labels_cex", 0.8) %>%
    	dendextend::set("labels",sppCLR_labs$day_labels)
     
     genus_iQlr_result_10k.dend <- genus_iQlr_result_10k.dend %>%
    	dendextend::set("branches_k_color", k=4) %>%
    	# set("hang_leaves", hang_height=80) %>% 
    	dendextend::set("leaves_pch", genusiQlr_labs$day_pch_labels) %>%  # node point type
    	dendextend::set("leaves_cex", 1) %>%  # node point size
    	dendextend::set("leaves_col", genusiQlr_labs$day_color_labels) %>% # node point color
    	dendextend::set("labels_cex", 0.8) %>%
    	dendextend::set("labels",genusiQlr_labs$day_labels)
 
    spp_iQlr_result_10k.dend <- spp_iQlr_result_10k.dend %>%
    	dendextend::set("branches_k_color", k=4) %>%
    	# set("hang_leaves", hang_height=80) %>% 
    	dendextend::set("leaves_pch", sspiQlr_labs$day_pch_labels) %>%  # node point type
    	dendextend::set("leaves_cex", 1) %>%  # node point size
    	dendextend::set("leaves_col", sspiQlr_labs$day_color_labels) %>% # node point color
    	dendextend::set("labels_cex", 0.8) %>%
    	dendextend::set("labels",sspiQlr_labs$day_labels)
    
    day0_CLR_result_10k.dend <- day0_CLR_result_10k.dend %>%
    	dendextend::set("branches_k_color", k=4) %>%
    	# set("hang_leaves", hang_height=80) %>% 
    	dendextend::set("leaves_pch", day0_labs$day_pch_labels) %>%  # node point type
    	dendextend::set("leaves_cex", 1) %>%  # node point size
    	dendextend::set("leaves_col", day0_labs$day_color_labels) %>% # node point color
    	dendextend::set("labels_cex", 0.8) %>%
    	dendextend::set("labels",day0_labs$day_labels) #%>%
    # plot(horiz =  TRUE)
    # 
 
    ##### tanglegrams ####
    	dendlist_compare <- dendlist(vst_result_10k.dend,ancom_result_10k.dend,
    															 comp_result_10k.dend,  #3
    															 genusCLR_result_10k.dend, #4
    															 sppCLR_result_10k.dend,   #5
    															 genus_iQlr_result_10k.dend,#6
    															 spp_iQlr_result_10k.dend) #7
    
    dendlist_CLR <- dendlist(sppCLR_result_10k.dend,
    												 ASVCLR_result_10k.dend)
    
    
    dendlist_CLR %>% 	#genus to spp CLR
    	dendlist(which = c(1,2)) %>% #entanglement #0.7336125
    	dendextend::untangle("random")			  %>%  #entanglement #0.05858792
    	dendextend::untangle("step1") 				%>%  #entanglement #0.1671132
    	# dendextend::untangle("step2side") 				%>%  #entanglement #0.1543208
    	# dendextend::untangle("random")			  %>%  #entanglement #0.05858792
    	dendextend::set("rank_branches")			%>%  #entanglement #
    	dendextend::pvclust_show_signif(dend = sppCLR_result_10k.dend, pvclust_obj = sppCLR_result_10k, 
    																	signif_type = "au", alpha = 0.05, show_type = "lwd") %>%
    	dendextend::tanglegram(center = TRUE, sort = F, k_labels = 4, k_branches = 4,
    												 margin_inner = 6,
    												 common_subtrees_color_lines = T, 
    												 common_subtrees_color_branches = T,
    												 color_lines = c(rep("azure4",    60))) %>%
    	dendextend::entanglement()
    
    
    
    pdf("r_output/visualizations/beta_div_plots/hier_tree_compare_VST_ANCOM.pdf",         # File name
    		width = 7, height = 8, useDingbats = T, compress = F)
    par(mar = c(4,4,3,4) + 0.1)
    # bottom, left, top and right margins 
    dendlist_compare %>% 
    	dendextend::dendlist(which = c(1:2)) %>%
    	dendextend::untangle("step1") %>%
    	dendextend::untangle("step2") %>%
    	dendextend::set("rank_branches") %>%
    	dendextend::tanglegram(main_left = "VST Euclidean dist.",
										 main_right = "ANCOM Euclidean dist.",
										 sub = "Ward.D2 clustering",
										 cex_sub = 0.95,
										 columns_width = c(5, 2, 5),
										 margin_inner = 6, sort = F,
										 k_labels = 4, k_branches = 4, 
										 center = TRUE, lab.cex = 1.05,
										 common_subtrees_color_lines = F,
										 common_subtrees_color_branches = T,
										 lwd = 3,margin_outer = 2,
										 color_lines = c(rep("azure4",    33),    #grey bars
										 								rep(viridis(6)[5], 2), #2 green
										 								rep(viridis(10)[6], 5), #5
										 								rep("azure4",      7), #7 grey
										 								rep(viridis(10)[1], 2), 
																		rep(viridis(10)[4], 4),
    																rep("azure4",      2),
																		rep(viridis(10)[2], 3),
																		rep("azure4",      2)))
  	 dev.off()
  	 
  	 length(unique(common_subtrees_clusters(dendlist_compare[[1]], dendlist_compare[[2]]))[-1]) #11 identical subtree
  	 

  	 dendlist_compare %>% 	#genus to spp CLR
  	 	dendlist(which = c(4:5)) %>%
  	 	dendextend::untangle("step1side") 				%>%  #entanglement #0.03
  	 	dendextend::untangle("step2side") 				%>%  #entanglement #0.002
  	 	# dendextend::untangle("random")			  %>%  #entanglement #
  	 	dendextend::set("rank_branches")			%>%  #entanglement #
  	 	dendextend::tanglegram(center = TRUE, sort = F, k_labels = 4, k_branches = 4,
						  	 						 margin_inner = 6,
						  	 						 common_subtrees_color_lines = T, 
						  	 						 common_subtrees_color_branches = T) %>%
  	 	dendextend::entanglement()
  	 
  	 dendlist_compare %>% #genus & spp iQlr
  	 	dendlist(which = c(6:7)) %>% #entanglement #0.334915
  	 	# dendextend::untangle("random")			  %>%  #entanglement #0.04472843
  	 	dendextend::untangle("step1") 				%>%  #entanglement #0.01752131
  	 	# dendextend::untangle("step2") 				%>%  #entanglement #0.002571784
  	  dendextend::set("rank_branches")			%>%  #entanglement #0.01603115
  	 	dendextend::tanglegram(center = T, sort = F,  
  	 												 k_labels = 4, k_branches = 4,
						  	 						 margin_inner = 6,
						  	 						 common_subtrees_color_lines = T,
						  	 						 common_subtrees_color_branches = T) %>%
  	 	dendextend::entanglement()                 #0.01215627
  	 
  	 
  	 
  	 ##### single plots #####
  	 mypal = ggsci::pal_npg(alpha = 0.9)(9)
  	 scales::show_col(mypal)
  	 scales::show_col(sppCLR_labs$day_color_labels)
  	 
  	 par(mar = c(1,1,1,1) + 0.1)
  	 
  	 Rotate.sppCLR_result_10k.dend <- dendextend::pvclust_show_signif(sppCLR_result_10k.dend, sppCLR_result_10k, 
  	 																	signif_type = "au", alpha = 0.05, show_type = "lwd") %>%
  	 	dendextend::click_rotate()
  	 
  	  Rotate.sppCLR_result_10k.dend <- dendextend::pvclust_show_signif(Rotate.sppCLR_result_10k.dend, sppCLR_result_10k, 
  	 																																 signif_type = "au", alpha = 0.05, show_type = "lwd") %>%
  	  dendextend::click_rotate()
  	 
  	 #realign colors and leaves 	
  	 rotate_labs <- dendro_pretty(Rotate.sppCLR_result_10k.dend)
  	 
  	 par(mar = c(1,7,7,1))
  	
  	 pdf("r_output/visualizations/beta_div_plots/hier_tree_boot10k_sppCLR.pdf",         # File name
  	 				width = 9, height = 8, useDingbats = F, compress = F)
  	 
  	 dendextend::pvclust_show_signif(Rotate.sppCLR_result_10k.dend, sppCLR_result_10k, 
  	  																signif_type = "au", alpha = 0.05, show_type = "lwd") %>%
  	 	dendextend::set("branches_k_color", k=4, value = mypal[1:4]) %>%
  	 	dendextend::set("leaves_col", value = rotate_labs$day_color_labels) %>%
  	 	dendextend::set("leaves_pch", 19) %>%
  	 	dendextend::set("leaves_cex", 2) %>%
  	 	circlize_dendrogram(labels = F, dend_track_height = .7, labels_track_height = 0.2)
  	 # title("Unsupervised clustering")
  	 # plot(horiz=TRUE, axes=F,  main = "Ward2 clustering of CLR species") #%>%
  	 	# dendextend::as.ggdend()	%>% ggplot()
  	 
  	 
  	 legend("topleft", inset= 0, legend=c("Elk 1", "Elk 2","Elk 3", "Elk 4"),
  	 			 col=c(mypal[1],mypal[3],mypal[2],mypal[4]), 
  	 			 lty=1, cex=0.8, lwd = 5, box.lty=0, bty = "n", xpd = T)

  	 legend("bottomleft", inset=c(0,0.05), title = "Bootstrap support",
  	 			 legend=c("p < 0.05", "p = ns"),
  	 			 lty=1, cex=0.8, lwd = c(5, 1.5), box.lty=0, bty = "n",
  	 			 col= c("darkgrey", "darkgrey"), xpd = T)
  
  	 legend("topleft", inset=c(0.027, 0.12),
  	 			 legend=c(" Day 0", " Day 1", " Day 3", " Day 7", " Day 14"),
  	 			 box.lty=0, border = 0,bty = "n",pt.cex = 1.5,
  	 			 col= c(1:5), pch = 19, horiz=F, cex=0.8, xpd = T)
  	 dev.off()
  	 # 'arg' should be one of “labels”, “labels_colors”, “labels_cex”, 
    # “labels_to_character”, “leaves_pch”, “leaves_cex”, “leaves_col”, 
    # “nodes_pch”, “nodes_cex”, “nodes_col”, “hang_leaves”, “rank_branches”, 
    # “branches_k_color”, “branches_k_lty”, “branches_col”, “branches_lwd”, 
    # “branches_lty”, “by_labels_branches_col”, “by_labels_branches_lwd”, 
    # “by_labels_branches_lty”, “by_lists_branches_col”, “by_lists_branches_lwd”, 
    # “by_lists_branches_lty”, “highlight_branches_col”, “highlight_branches_lwd”, “clear_branches”, “clear_leaves”
    
#### Rarefaction curve ####
    # otu_table(subset_samples(elk_rdp.ps, Elk == "1"))
    library(MicrobiotaProcess)
    library(patchwork)
    # for reproducibly random number
    set.seed(1024)
    rareres <- get_rarecurve(obj=elk_rdp.ps, chunks=400)
    rareres$data$Elk <- as.factor(paste0("Elk.",rareres$data$Elk))
    rareres$data$Day <- as.factor(paste0("Elk.",rareres$data$Elk,"_day.",rareres$data$Day ))
    options(scipen = 9999)
    p_rare <- ggrarecurve(obj=rareres,factorNames="Day",linesize = 0.75,
    											indexNames=c("Observe","Chao1","ACE"),
												  ) +
												    	theme(legend.spacing.y=unit(0.01,"cm"),
												    				legend.text=element_text(size=4))
    
    prare1 <- ggrarecurve(obj=rareres, factorNames="Elk",linesize = 0.75,se = T,
    											indexNames=c("Observe", "Chao1", "ACE")
											    ) +
											    	scale_fill_manual(values=
											    											c(colorRampPalette(RColorBrewer::brewer.pal(4,"Set2"))(4))
											    										 )+
											    	scale_color_manual(values=
											    										 	c(colorRampPalette(RColorBrewer::brewer.pal(4,"Set2"))(4))
											    										 )+
											    	theme_bw()+
											    	theme(axis.text=element_text(size=8), panel.grid=element_blank(),
											    				strip.background = element_rect(colour=NA,fill="grey"),
											    				strip.text.x = element_text(face="bold"))          
    
    prare2 <- ggrarecurve(obj=rareres, factorNames="Elk",linesize = 0.75,
    											shadow=FALSE,
    											indexNames=c("Observe", "Chao1", "ACE")
											    ) +
											    	scale_color_manual(values=
											    										 	c(colorRampPalette(RColorBrewer::brewer.pal(4,"Set2"))(4))
											    										 )+
											    	theme_bw()+
											    	theme(axis.text=element_text(size=8), panel.grid=element_blank(),
											    				strip.background = element_rect(colour=NA,fill="grey"),
											    				strip.text.x = element_text(face="bold"))
    # p_rare / prare1 / prare2
    
    ggsave2("r_output/visualizations/rare_curve_day1.pdf", 
    				plot = p_rare, width = 9, height = 4, units = "in")
    ggsave2("r_output/visualizations/rare_curve_elk1.pdf", 
    				plot = prare1, width = 9, height = 4, units = "in")
    ggsave2("r_output/visualizations/rare_curve_elk2.pdf", 
    				plot = prare2, width = 9, height = 4, units = "in")
    
#### Richness #####
     #Obtain breakaway estimates
     #
    sample_data(elk_rdp.ps)$Elk <- as.factor(sample_data(elk_rdp.ps)$Elk)
    sample_data(elk_rdp.ps)$Rep <- as.factor(sample_data(elk_rdp.ps)$Rep)
    sample_data(elk_rdp.ps)$Day <- as.numeric(as.character(sample_data(elk_rdp.ps)$Day))
    set.seed(3313)
    ba_adiv <- breakaway::breakaway(elk_rdp.ps)
     ba_adiv[31]
	     ba_adiv$ETS_1431 %$% plot
	  #this one is sensitive to singletons so I'll swap it out
	   ba_adiv_nof1 <-breakaway_nof1(elk_rdp.ps)
	  ba_adiv_nof1[31]
	  ba_adiv_nof1$ETS_1431 %$% plot
	  
	  ba_adiv[31] <- ba_adiv_nof1[31]
	  #Plot estimates ugly
     # plot(ba_adiv, elk_rdp.ps, label = "Day",  color = "Elk", shape = "Day")
     # plot(ba_adiv_nof1, elk_rdp.ps, label = "Day",  color = "Elk", shape = "Day")
     #clean it up for publication
     ba_tibble <- summary(ba_adiv) %>% as_tibble
    
     # WARNING these tests have pseudo replication
     # Test if Richness changes by Day as fixed with all technical reps
     bt_rich_test_DayFixed_all_samps <- breakaway::betta(formula = estimate ~  Day,	
     																					data = ba_tibble, ses = error)
     bt_rich_test_DayFixed_all_samps$table #ns but also maybe incorporating psuedo reps
     
     bt_rich_test_DayElkFix_all_samps <- breakaway::betta(formula = estimate ~  Day + Elk,	
     																										data = ba_tibble, ses = error)
     bt_rich_test_DayElkFix_all_samps$table #ns but also maybe incorporating psuedo reps
     
     bt_rich_test_DayFixElkRandom_all_samps <- breakaway::betta_random(formula = estimate ~  Day | Elk,	
     																										 data = ba_tibble, ses = error)
     bt_rich_test_DayFixElkRandom_all_samps$table #ns but also maybe incorporating psuedo reps
     
     ##plot all estimates including technical replicates
     #reorder sample_data(elk_rdp.ps)
     ba_tibble <- ba_tibble[sample_names = (sample_data(elk_rdp.ps)$Description),]
    
     #add some factors
     ba_tibble$Elk <- sample_data(elk_rdp.ps)$Elk
     ba_tibble$Rep <- sample_data(elk_rdp.ps)$Rep
     ba_tibble$Day <- sample_data(elk_rdp.ps)$Day
     ba_tibble$SampleID <- sample_data(elk_rdp.ps)$SampleID
     ba_tibble$Elk <- as.factor(ba_tibble$Elk)
     ba_tibble$Rep <- as.factor(ba_tibble$Rep)
     # ba_tibble$Day <- as.factor(ba_tibble$Day)#wrong
     
     #reorder the factors
     ba_tibble <- ba_tibble %>%
     	arrange(Elk, Day, Rep) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
     	mutate(sample_names=factor(sample_names, levels=sample_names))
     print(ba_tibble, n = 60)
     
     options(scipen=999)
     labels_elk <- c("Elk 1", "Elk 2", "Elk 3", "Elk 4")
     names(labels_elk) <- c("1", "2", "3", "4")
     # ba_tibble[which(ba_tibble$Day==14),]
     #plot
     ba_richness_pop <-  ggplot(ba_tibble,
     													 aes(x = Day, y = estimate,
     															 ymin  = lower, ymax = upper,
     															 color = as.factor(Day), fill = Rep)) +
     														lims(y = c(0, 4000)) + 
     														ylab("breakaway richness estimate") +
     														xlab("Day / Rep") + labs(color = "Day") +
	     													theme_clean() +
     														panel_border() +
    														guides(fill = "none")+
    														theme(legend.background = element_blank(),
     																	plot.background = element_blank()) +
     														geom_pointrange(position = position_jitterdodge(dodge.width = .7, jitter.width = 0)) + 
     														scale_y_continuous(limits = c(0,max(ba_tibble$upper)), 
     																							 breaks = c(0,500,1000,2000,3000,4000), expand = c(0, 0)) +
     														coord_cartesian(ylim=c(0, 4000))+
     														scale_x_continuous(limits = c(-1,15), breaks = c(0,1,3,7,14)) +
     														theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
     														facet_wrap(~Elk, shrink = T, scales = "free_x", nrow = 1,
     																			  labeller = labeller(Elk = labels_elk))
     ba_richness_pop
     ggsave2(filename = "r_output/visualizations/brkawy_richness_elk_long.pdf", 
     				 plot = ba_richness_pop, 
     				 width = 8, height = 3, units = "in")
     
     # Because the data contains technical replicates (that are not independant) 
     # we need to average these values for each biological replicate or it will inflate our significance
     # library(dplyr)
    ba_tibble_aveTechReps <- ba_tibble %>% 	ungroup() %>%
    		mutate(Day = as.factor(Day)) %>%
     	 	group_by(Elk, Day) %>% 
     	 	dplyr::summarize(Ave_estimate = mean(estimate), 
     	 						ave_error = mean(error),
     	 						ave_lower = mean(lower),
     	 						ave_upper = mean(upper)) %>% 
     	 	ungroup() %>% 
     	 	mutate(Day = as.numeric(as.character(Day)), Elk = factor(Elk), sample_name = paste0("Elk", `Elk`, ".Day.", `Day`)) %>%
    		arrange(Elk, Day) %>%
    		mutate(sample_name=factor(sample_name, levels=sample_name))
    
    #plot
    ba_richness_AveTechReps<-  ggplot(ba_tibble_aveTechReps,
    													 aes(x = Day, y = Ave_estimate, 
    													 		ymin  = ave_lower, ymax = ave_upper,
    													 		color = as.factor(Day))) +
														    	lims(y = c(0, 4000)) + 
														    	ylab("Ave. richness estimate") +
														    	xlab("Day") + labs(color = "Day") +
    															geom_pointrange(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0)) + 
    															scale_x_continuous(limits = c(-1,15), breaks = c(0,1,3,7,14)) +														
														    	scale_y_continuous(limits = c(0,max(ba_tibble_aveTechReps$ave_upper)), breaks = c(0,500,1000,2000,3000,4000), expand = c(0, 0)) +
    															coord_cartesian(ylim=c(0, 4000))+
														    	theme_clean()+
														    	panel_border() +
														    	theme(legend.background = element_blank(),
														    				plot.background = element_blank(),
														    				axis.text.x = element_blank(),
														    				axis.ticks.x = element_blank()) +
														    	facet_wrap(~Elk, shrink = T, scales = "free_x", nrow = 1,
														    						 labeller = labeller(Elk = labels_elk))
    ba_richness_AveTechReps
    ggsave2(filename = "r_output/visualizations/brkawy_richness_elk_long_aveReps.pdf", 
    				plot = ba_richness_AveTechReps, 
    				width = 8, height = 3, units = "in")
    
   richness_plotPair <- ggpubr::ggarrange(ba_richness_pop, ba_richness_AveTechReps, nrow = 2, labels = c("A.", "B."),
    									common.legend = TRUE, legend = "right")
   richness_plotPair
   
   ggsave2(filename = "r_output/visualizations/brkawy_richness_elk_pair.pdf", 
   				plot = richness_plotPair, 
   				width = 8, height = 6, units = "in")
   
   
   #Examine models of Richness
   # Test if Richness changes by Day as fixed
    bt_rich_test_DayFixed <- breakaway::betta(formula = Ave_estimate ~  Day,	
    																					data = ba_tibble_aveTechReps, ses = ave_error)
    bt_rich_test_DayFixed$table #ns
    
    bt_rich_test_DayplusElkFixed <- breakaway::betta(formula = Ave_estimate ~  Day + Elk,	
    																					data = ba_tibble_aveTechReps, ses = ave_error)
    bt_rich_test_DayplusElkFixed$table #ns

    bt_rich_test_DayplusElkRandom <- breakaway::betta_random(formula = Ave_estimate ~  Day | Elk,	
    																												data = ba_tibble_aveTechReps, ses = ave_error)
    bt_rich_test_DayplusElkRandom$table #ns
    
    #not as informative since all comparisons are based on Elk1 intercept & slope
    
    # experiment by filtering down to 1 elk and using Day as fixed to parse out each slope
     ba_tibble_elk1 <- ba_tibble_aveTechReps %>% filter(Elk==1)
     ba_tibble_elk2 <- ba_tibble_aveTechReps %>% filter(Elk==2)
     ba_tibble_elk3 <- ba_tibble_aveTechReps %>% filter(Elk==3)
     ba_tibble_elk4 <- ba_tibble_aveTechReps %>% filter(Elk==4)
     
     bt_elk1_rich_day <-     breakaway::betta(formula = Ave_estimate ~ Day, data = ba_tibble_elk1, ses = ave_error)
     bt_elk1_rich_day$table #sig
    
    ##
     bt_elk2_rich_day <-     breakaway::betta(formula = Ave_estimate ~ Day, data = ba_tibble_elk2, ses = ave_error)
     bt_elk2_rich_day$table #ns
       
     ##
     bt_elk3_rich_day <-     breakaway::betta(formula = Ave_estimate ~ Day, data = ba_tibble_elk3, ses = ave_error)
     bt_elk3_rich_day$table #ns

     ##
     bt_elk4_rich_day <-     breakaway::betta(formula = Ave_estimate ~ Day, data = ba_tibble_elk4, ses = ave_error)
     bt_elk4_rich_day$table #ns
     
	   bt_elk4_rich_day$table[1] #intercept
	   bt_elk4_rich_day$table[2] #slope
	   bt_elk4_rich_day$table[6] #P-value
	  
	   ba_tibble_aveTechReps_wplot <- ba_tibble_aveTechReps %>% 
	   													 mutate(slope = c(rep(bt_elk1_rich_day$table[2], 5),
	   																					  rep(bt_elk2_rich_day$table[2], 5), 
	   																					  rep(bt_elk3_rich_day$table[2], 5),
	   																					  rep(bt_elk4_rich_day$table[2], 5)),
	   																 intercept = c(rep(bt_elk1_rich_day$table[1], 5),
	   																 							 rep(bt_elk2_rich_day$table[1], 5),
	   																 							 rep(bt_elk3_rich_day$table[1], 5),
	   																 							 rep(bt_elk4_rich_day$table[1], 5)),
	   																     pval = c(rep(bt_elk1_rich_day$table[6], 5),
	   																 							rep(bt_elk2_rich_day$table[6], 5),
	   																 							rep(bt_elk3_rich_day$table[6], 5),
	   																 							rep(bt_elk4_rich_day$table[6], 5)),
	   																		locX_label = c(rep(1, 20)),
	   																		locY_label = c(rep(100, 20)))

     ba_richness_AveTechReps_slopes<-  ggplot(ba_tibble_aveTechReps_wplot,
     																	        aes(x = Day, y = Ave_estimate, 
     																		        	ymin  = ave_lower, ymax = ave_upper,
     																			        color = as.factor(Day))) +
    																	geom_pointrange(inherit.aes = T) + 
																     	ylab("Ave. richness estimate") +
																     	xlab("Day") + labs(color = "Day") +
																     	scale_y_continuous(limits = c(0,max(ba_tibble_aveTechReps_wplot$ave_upper)), 
																     										 breaks = c(0,500,1000,2000, 3000, 4000), expand = c(0, 0)) +
     																	coord_cartesian(ylim=c(0, 2500))+
																     	scale_x_continuous(limits = c(0,14), breaks = c(0,1,3,7,14)) +
																     	facet_wrap(~Elk, shrink = T, scales = "free_x", nrow = 1,
																     						 labeller = labeller(Elk = labels_elk)) +
     																	geom_abline(aes(intercept = intercept, slope = slope), 
     																							color = "darkgrey") +
     																	geom_text(data = ba_tibble_aveTechReps_wplot,aes(x = as.numeric(locX_label),
    																								y = as.numeric(locY_label), label = paste("p-value =",pval)), color = "darkgrey",
     																						hjust = 0, vjust = 0) +
																     	theme_clean()+ panel_border() +
																	    theme(legend.background = element_blank(),
																     				plot.background = element_blank(),
																     				axis.text.x = element_blank(), 
																     				axis.ticks.x = element_blank()) 
     ba_richness_AveTechReps_slopes
     ggsave2(filename = "r_output/visualizations/brkawy_richness_elk_long_aveReps_slopes.pdf", 
     				plot = ba_richness_AveTechReps_slopes, 
     				width = 8, height = 3, units = "in")
     
     richness_plotPair3 <- ggpubr::ggarrange(ba_richness_pop, ba_richness_AveTechReps_slopes, 
     																			 nrow = 2, labels = c("A.", "B."),
     																			 common.legend = TRUE, legend = "right")
     richness_plotPair3
     ggsave2(filename = "r_output/visualizations/brkawy_richness_elk_pair3.pdf", 
     				plot = richness_plotPair3, 
     				width = 8, height = 6, units = "in")
     
     ### regression on all elk by Day, func estimate ~ Day |Elk
     bt_elk4_rich_day$table[1] #intercept
     bt_elk4_rich_day$table[2] #slope
     bt_rich_test_DayplusElkRandom$table[3] #error
     bt_elk4_rich_day$table[6] #P-value
     
     bt_rich_test_DayplusElkRandom
     bt_rich_test_DayplusElkRandom$cov
     bt_rich_test_DayplusElkRandom$function.args
     bt_rich_test_DayplusElkFixed$table
     ba_richness_AveTechReps_GlobalSlope <-  ggplot(ba_tibble_aveTechReps,
     																								aes(x = Day, y = Ave_estimate, 
     																				 						ymin  = ave_lower, ymax = ave_upper,
     																				 						color = as.factor(Day), shape = Elk)) +
     																					geom_pointrange(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0)) +
     																					scale_shape_manual(values=c(1, 2, 0, 5))+																	     	
     																					ylab("Ave. richness estimate") + xlab("Day") + 
     																					labs(color = "Day", shape = "Elk", subtitle = "All Elk: betta model = richness~ Day | Elk)") +
																				     	scale_y_continuous(limits = c(0,max(ba_tibble_aveTechReps$ave_upper)), 
																				     										 breaks = c(0,500,1000,2000, 3000, 4000), expand = c(0, 0)) +
																				     	coord_cartesian(ylim=c(0, 2500))+
																				     	scale_x_continuous(limits = c(-1,15), breaks = c(0,1,3,7,14)) +
																				     	geom_abline(intercept = bt_rich_test_DayplusElkRandom$table[1], 
																				     							slope = bt_rich_test_DayplusElkRandom$table[2], 
																				     							color = "darkgrey") +
																				     	geom_text(x = -0.5, y = 100, 
																				     								label = paste("p-value =",bt_rich_test_DayplusElkRandom$table[6]), 
																				     								color = "darkgrey", hjust = 0, vjust = 0) +
																				     	theme_clean()+ panel_border() +
																				     	theme(plot.subtitle = element_text(size = 11.5),
																				     				legend.background = element_blank(),
																				     				plot.background = element_blank(),
																				     				# axis.text.x = element_blank(), 
																				     				axis.ticks.x = element_blank()) 
     ba_richness_AveTechReps_GlobalSlope
		
     ggsave2(filename = "r_output/visualizations/brkawy_richness_elk_aveReps_GlobalSlope.pdf", 
					plot = ba_richness_AveTechReps_GlobalSlope, 
					width = 4, height = 3, units = "in")
     
     richness_plotPair3 <- ggpubr::ggarrange(ba_richness_pop,
     																				ba_richness_AveTechReps_slopes, 
     																				ba_richness_AveTechReps_GlobalSlope,
     																			  nrow = 3, labels = c("A.", "B.", "C."),
     																				legend.grob = get_legend(ba_richness_AveTechReps_GlobalSlope),
     																			  legend = "right")
     
     richness_plotPair3
     ggsave2(filename = "r_output/visualizations/brkawy_richness_pair3.pdf", 
     				plot = richness_plotPair3, 
     				width = 8, height = 8, units = "in")
#### Alpha Diversity #####
     # Shannon–Wiener index downweights rare species and gives more weight to evenness and common species.
     # The effect of the sample size is generally negligible for both of them. 
     # The Simpson diversity index is a similarity index where the higher the value, the lower the diversity. 
     # It measures the probability that two individuals randomly selected from a sample
     # will belong to the same species. With this index, 0 represents infinite
     # diversity and 1, no diversity.
     # 
     #only look at genus level diversity
     genus_data.ps = tax_glom(elk_rdp.ps, "Genus")
     # merged_genus_data.ps = merge_samples(genus_data.ps, "SampleType")
     meta <- genus_data.ps %>%
       sample_data %>% 
       as_tibble %>% 
       select(c("Elk", "Day", "Rep")) %>% 
       mutate("sample_names" = genus_data.ps %>% sample_names )
     # meta$sample_names
     
     #plug-in type estimates for example
     genus_data.ps %>% sample_richness %>% plot
     genus_data.ps %>% sample_inverse_simpson %>% plot
     genus_data.ps %>% sample_shannon %>% plot
     
     #run divnet for distances and errors
     #Let’s treat all samples as independent observations (X = NULL) and fit the DivNet model:
     dv <- DivNet::divnet(genus_data.ps, X = NULL, ncores = 30)
     dv %>% names
     dv$shannon %>% head
		##### Shannon #######
     #make summary table by sample
     combined_shannon <- meta %>%
     	left_join( dv$shannon %>% summary,
     						 by = "sample_names")
     
     combined_shannon <- combined_shannon %>% arrange(Elk, Day, Rep) %>%
     	mutate(sample_names=factor(sample_names, levels=sample_names))
     
     # plot with ggplot 
     # 
     Div_shann_gg <-  ggplot(combined_shannon,
    													aes(x = Day, y = estimate, 
     																ymin  = lower, ymax = upper,
     																color = as.factor(Day), fill = Rep)) +
															ylab("DivNet Shannon div. estimate") +
															xlab("Day / Rep") + labs(color = "Day") +
															theme_clean()+
															panel_border() +	guides(fill = "none")+
															theme(legend.background = element_blank(),
																		plot.background = element_blank()) +
     													geom_pointrange(position = position_jitterdodge(dodge.width = .7, jitter.width = 0)) +
															scale_y_continuous(limits = c(2,4), 
																								 expand = c(0, 0)) +
															coord_cartesian(ylim=c(2.9,3.25))+
															scale_x_discrete(labels = combined_shannon$Day, expand = c(0.1, 0)) +
															theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
															facet_wrap(~Elk, shrink = T, scales = "free_x", nrow = 1,
																				 labeller = labeller(Elk = labels_elk))
     Div_shann_gg
ggsave2(filename = "r_output/visualizations/DivNet_Shann_elk_wTechnicalReps.pdf", 
				plot = Div_shann_gg, 
				width = 8, height = 3, units = "in")

#test for group differences

# Because the data contains technical replicates (that are not independent) 
# we need to average these values for each biological replicate or it will inflate our significance

combined_shannon_aveTechReps <- combined_shannon %>% 
																		group_by(Elk, Day) %>% 
																		dplyr::summarize(Ave_estimate = mean(estimate), 
																							ave_error = mean(error),
																							ave_lower = mean(lower),
																							ave_upper = mean(upper)) %>% 
																		ungroup() %>% 
																		mutate(Elk = factor(Elk), sample_name = paste0("Elk", `Elk`, ".Day.", `Day`)) %>%
																		arrange(Elk, Day) %>%
																		mutate(sample_name=factor(sample_name, levels=sample_name))


#Examine models of Shannon
# Test if Shannon changes by Day as fixed
# Dont model with random effects because not enough data
DivN_test_DayFixed <- breakaway::betta(formula = Ave_estimate ~  Day,	
																					data = combined_shannon_aveTechReps, ses = ave_error)
DivN_test_DayFixed$table

DivN_test_DayFixedplusElkFixed <- breakaway::betta(formula = Ave_estimate ~  Day + Elk,	
																								 data = combined_shannon_aveTechReps, ses = ave_error)
DivN_test_DayFixedplusElkFixed$table

DivN_test_DayFixedplusElkRandom <- breakaway::betta_random(formula = Ave_estimate ~  Day | Elk,	
																									 data = combined_shannon_aveTechReps, ses = ave_error)
DivN_test_DayFixedplusElkRandom$table #ns


#not as informative since all comparisons are based on Elk1 intercept & slope

# experiment by filtering down to 1 elk and using Day as fixed to parse out each slope
DivN_test_elk1 <- combined_shannon_aveTechReps %>% filter(Elk==1)
DivN_test_elk2 <- combined_shannon_aveTechReps %>% filter(Elk==2)
DivN_test_elk3 <- combined_shannon_aveTechReps %>% filter(Elk==3)
DivN_test_elk4 <- combined_shannon_aveTechReps %>% filter(Elk==4)

DivN_test_elk1_Shann_day <-     breakaway::betta(formula = Ave_estimate ~ Day, data = DivN_test_elk1, ses = ave_error)
DivN_test_elk1_Shann_day$table #sig
##
DivN_test_elk2_Shann_day <-     breakaway::betta(formula = Ave_estimate ~ Day, data = DivN_test_elk2, ses = ave_error)
DivN_test_elk2_Shann_day$table #ns
##
DivN_test_elk3_Shann_day <-     breakaway::betta(formula = Ave_estimate ~ Day, data = DivN_test_elk3, ses = ave_error)
DivN_test_elk3_Shann_day$table #ns
##
DivN_test_elk4_Shann_day <-     breakaway::betta(formula = Ave_estimate ~ Day, data = DivN_test_elk4, ses = ave_error)
DivN_test_elk4_Shann_day$table #ns
DivN_test_elk4_Shann_day$table[1] #intercept
DivN_test_elk4_Shann_day$table[2] #slope
DivN_test_elk4_Shann_day$table[6] #P-value

combined_shannon_aveTechReps_wplot <- combined_shannon_aveTechReps %>% 
	mutate(slope = c(rep(DivN_test_elk1_Shann_day$table[2], 5),
									 rep(DivN_test_elk2_Shann_day$table[2], 5), 
									 rep(DivN_test_elk3_Shann_day$table[2], 5),
									 rep(DivN_test_elk4_Shann_day$table[2], 5)),
				intercept = c(rep(DivN_test_elk1_Shann_day$table[1], 5),
				 							rep(DivN_test_elk2_Shann_day$table[1], 5),
				 							rep(DivN_test_elk3_Shann_day$table[1], 5),
				 							rep(DivN_test_elk4_Shann_day$table[1], 5)),
				pval = c(rep(0.0001, 5),
				 				 rep(DivN_test_elk2_Shann_day$table[6], 5),
				 				 rep(DivN_test_elk3_Shann_day$table[6], 5),
				 				 rep(DivN_test_elk4_Shann_day$table[6], 5)),
				 locX_label = c(rep(0.5, 20)),
				 locY_label = c(rep(2.95, 20)))


DivN_Shann_AveTechReps_slopes <-  ggplot(combined_shannon_aveTechReps_wplot,
																				 aes(x = Day, y = Ave_estimate, 
																				 		ymin  = ave_lower, ymax = ave_upper,
																				 		color = as.factor(Day))) +
																				geom_pointrange() + 
																				ylab("Ave. Shannon div. estimate") +
																				xlab("Day") + labs(color = "Day") +
																				theme_clean()+ panel_border() +
																				theme(legend.background = element_blank(),
																							plot.background = element_blank()) +
																				scale_y_continuous(limits = c(2,4), 
																													 expand = c(0, 0)) +
																				coord_cartesian(ylim=c(2.9,3.25))+
																				scale_x_discrete(labels = combined_shannon_aveTechReps_wplot$Day, expand = c(0.1, 0)) +
																				theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
																				facet_wrap(~Elk, shrink = T, scales = "free_x", nrow = 1,
																									 labeller = labeller(Elk = labels_elk)) +
																				geom_abline(aes(intercept = intercept, slope = slope), 
																										color = "darkgrey") +
																				geom_text(aes(x = as.numeric(locX_label),
																											y = as.numeric(locY_label), 
																											label = paste("p-value =",pval)), color = "darkgrey",
																											hjust = 0, vjust = -.5)
																			
DivN_Shann_AveTechReps_slopes


DivN_shann_plotPair2 <- ggpubr::ggarrange(Div_shann_gg, DivN_Shann_AveTechReps_slopes, 
																				nrow = 2, labels = c("A.", "B."),
																				common.legend = TRUE, legend = "right")
DivN_shann_plotPair2
ggsave2(filename = "r_output/visualizations/DivN_shann_elk_pair2.pdf", 
				plot = DivN_shann_plotPair2, 
				width = 8, height = 6, units = "in")


#plot the global slope test

DivN_Shann_AveTechReps_GlobalSlope <-  ggplot(combined_shannon_aveTechReps,
																							 aes(x = Day, y = Ave_estimate, 
																							 		ymin  = ave_lower, ymax = ave_upper,
																							 		color = as.factor(Day), shape = Elk)) +
																						geom_pointrange(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0)) +
																						scale_shape_manual(values=c(1, 2, 0, 5))+																	     	
																						ylab("Ave. Shannon div. estimate") + xlab("Day") + 
																						labs(color = "Day", shape = "Elk", subtitle = "All Elk: betta model = Shannon ~ Day | Elk)") +
																						scale_y_continuous(limits = c(2,4), 
																															 expand = c(0, 0)) +
																						coord_cartesian(ylim=c(2.9,3.25))+
																						scale_x_continuous(limits = c(-1,15), breaks = c(0,1,3,7,14)) +
																						geom_abline(intercept = DivN_test_DayFixedplusElkRandom$table[1], 
																												slope = DivN_test_DayFixedplusElkRandom$table[2], 
																												color = "darkgrey") +
																						geom_text(x = -0.5, y = 2.91, 
																											label = paste("p-value =",DivN_test_DayFixedplusElkRandom$table[6]), 
																											color = "darkgrey", hjust = 0, vjust = 0) +
																						theme_clean()+ panel_border() +
																						theme(plot.subtitle = element_text(size = 11.5),
																									legend.background = element_blank(),
																									plot.background = element_blank(),
																									# axis.text.x = element_blank(), 
																									axis.ticks.x = element_blank()) 
DivN_Shann_AveTechReps_GlobalSlope

ggsave2(filename = "r_output/visualizations/divN_Shann_elk_aveReps_GlobalSlope.pdf", 
				plot = DivN_Shann_AveTechReps_GlobalSlope, 
				width = 4, height = 3, units = "in")

Shann_plotPair3 <- ggpubr::ggarrange(Div_shann_gg,
																				DivN_Shann_AveTechReps_slopes, 
																				DivN_Shann_AveTechReps_GlobalSlope,
																				nrow = 3, labels = c("A.", "B.", "C."),
																				legend.grob = get_legend(DivN_Shann_AveTechReps_GlobalSlope),
																				legend = "right")

Shann_plotPair3
ggsave2(filename = "r_output/visualizations/divN_Shann_pair3.pdf", 
				plot = Shann_plotPair3, 
				width = 8, height = 8, units = "in")

		##### Simpson ####

#plug-in type estimates for example
genus_data.ps %>% sample_simpson %>% plot
dv$simpson %>% head

#make summary table by sample
combined_simps <- meta %>%
	left_join( dv$simpson %>% summary,
						 by = "sample_names")

combined_simps <- combined_simps %>% arrange(Elk, Day, Rep) %>%
	mutate(sample_names=factor(sample_names, levels=sample_names))

# plot with ggplot 
Div_simps_gg <-  ggplot(combined_simps,
												aes(x = Day, y = estimate, 
														ymin  = lower, ymax = upper,
														color = as.factor(Day), fill = Rep)) +
												ylab("DivNet Simpson div. estimate") +
												xlab("Day / Rep") + labs(color = "Day") +	guides(fill = "none")+
												theme_clean()+ panel_border() +
												theme(legend.background = element_blank(),
															plot.background = element_blank()) +
												geom_pointrange(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0)) +
												scale_y_continuous(limits = c(0.060, 0.09),
																					 expand = c(0, 0)) +
												coord_cartesian(ylim=c(0.06, 0.085))+
												scale_x_continuous(limits = c(-1,15), breaks = c(0,1,3,7,14)) +
												theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
												facet_wrap(~Elk, shrink = T, scales = "free_x", nrow = 1,
																	 labeller = labeller(Elk = labels_elk))
Div_simps_gg
ggsave2(filename = "r_output/visualizations/DivNet_Simpson_elk_wTechnicalReps.pdf", 
				plot = Div_simps_gg, 
				width = 8, height = 3, units = "in")

#test for group differences

# Because the data contains technical replicates (that are not independent) 
# we need to average these values for each biological replicate or it will inflate our significance

combined_simps_aveTechReps <- combined_simps %>% 
	group_by(Elk, Day) %>% 
	dplyr::summarize(Ave_estimate = mean(estimate), 
						ave_error = mean(error),
						ave_lower = mean(lower),
						ave_upper = mean(upper)) %>% 
	ungroup() %>% 
	mutate(Elk = factor(Elk), sample_name = paste0("Elk", `Elk`, ".Day.", `Day`)) %>%
	arrange(Elk, Day) %>%
	mutate(sample_name=factor(sample_name, levels=sample_name))


#Examine models of Simpson
# Test if Simpson changes by Day as fixed
# Dont model with random effects because not enough data
DivN_test_Simps_DayFixed <- breakaway::betta(formula = Ave_estimate ~  Day,	
																			 data = combined_simps_aveTechReps, ses = ave_error)
DivN_test_Simps_DayFixed$table

DivN_test_Simps_DayFixedplusElkFixed <- breakaway::betta(formula = Ave_estimate ~  Day + Elk,	
																									 data = combined_simps_aveTechReps, ses = ave_error)
DivN_test_Simps_DayFixedplusElkFixed$table

DivN_test_Simps_DayFixedplusElkRandom <- breakaway::betta_random(formula = Ave_estimate ~  Day | Elk,	
																												 data = combined_simps_aveTechReps, ses = ave_error)
DivN_test_Simps_DayFixedplusElkRandom$table
#not informative since all comparisons are based on Elk1 intercept & slope

# experiment by filtering down to 1 elk and using Day as fixed to parse out each slope
DivN_Simpstest_elk1 <- combined_simps_aveTechReps %>% filter(Elk==1)
DivN_Simpstest_elk2 <- combined_simps_aveTechReps %>% filter(Elk==2)
DivN_Simpstest_elk3 <- combined_simps_aveTechReps %>% filter(Elk==3)
DivN_Simpstest_elk4 <- combined_simps_aveTechReps %>% filter(Elk==4)

DivN_test_elk1_Simps_day <-     breakaway::betta(formula = Ave_estimate ~ Day, data = DivN_Simpstest_elk1, ses = ave_error)
DivN_test_elk1_Simps_day$table #sig
##
DivN_test_elk2_Simps_day <-     breakaway::betta(formula = Ave_estimate ~ Day, data = DivN_Simpstest_elk2, ses = ave_error)
DivN_test_elk2_Simps_day$table #ns
##
DivN_test_elk3_Simps_day <-     breakaway::betta(formula = Ave_estimate ~ as.numeric(Day), data = DivN_Simpstest_elk3, ses = ave_error)
DivN_test_elk3_Simps_day$table #ns
##
DivN_test_elk4_Simps_day <-     breakaway::betta(formula = Ave_estimate ~ as.numeric(Day), data = DivN_Simpstest_elk4, ses = ave_error)
DivN_test_elk4_Simps_day$table #ns

DivN_test_elk4_Simps_day$table[1] #intercept
DivN_test_elk4_Simps_day$table[2] #slope
DivN_test_elk4_Simps_day$table[6] #P-value

combined_simpson_aveTechReps_wplot <- combined_simps_aveTechReps %>% 
	mutate(slope = c(rep(DivN_test_elk1_Simps_day$table[2], 5),
									 rep(DivN_test_elk2_Simps_day$table[2], 5), 
									 rep(DivN_test_elk3_Simps_day$table[2], 5),
									 rep(DivN_test_elk4_Simps_day$table[2], 5)),
				intercept = c(rep(DivN_test_elk1_Simps_day$table[1], 5),
				 							rep(DivN_test_elk2_Simps_day$table[1], 5),
				 							rep(DivN_test_elk3_Simps_day$table[1], 5),
				 							rep(DivN_test_elk4_Simps_day$table[1], 5)),
				pval = c(rep(0.0001, 5),
				 				 rep(DivN_test_elk2_Simps_day$table[6], 5),
				 				 rep(DivN_test_elk3_Simps_day$table[6], 5),
				 				 rep(0.0001, 5)),
				 locX_label = c(rep(1, 20)),
				 locY_label = c(rep(0.061, 20)))

DivN_Simps_AveTechReps_slopes <-  ggplot(combined_simpson_aveTechReps_wplot,
																				 aes(x = Day, y = Ave_estimate, 
																				 		ymin  = ave_lower, ymax = ave_upper,
																				 		color = as.factor(Day))) +
																	geom_pointrange() + 
																	ylab("Ave. Simpson div. estimate") +
																	xlab("Day") + labs(color = "Day") + guides(fill = "none")+
																	theme_clean()+ panel_border() +
																	theme(legend.background = element_blank(),
																				plot.background = element_blank()) +
																	scale_y_continuous(limits = c(0.060, 0.09),
																										 expand = c(0, 0)) +
																	coord_cartesian(ylim=c(0.06, 0.085))+
																	scale_x_continuous(limits = c(-1,15), breaks = c(0,1,3,7,14)) +
																	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
																	facet_wrap(~Elk, shrink = T, scales = "free_x", nrow = 1,
																						 labeller = labeller(Elk = labels_elk)) +
																	geom_abline(aes(intercept = intercept, slope = slope), 
																							color = "darkgrey") +
																	geom_text(aes(x = as.numeric(locX_label),
																								y = as.numeric(locY_label), 
																								label = paste("p-value =",pval)), color = "darkgrey",
																								hjust = .08, vjust = -.5)

DivN_Simps_AveTechReps_slopes


DivN_simps_plotPair2 <- ggpubr::ggarrange(Div_simps_gg, DivN_Simps_AveTechReps_slopes, 
																					nrow = 2, labels = c("A.", "B."),
																					common.legend = TRUE, legend = "right")
DivN_simps_plotPair2
ggsave2(filename = "r_output/visualizations/DivN_simps_elk_pair2.pdf", 
				plot = DivN_simps_plotPair2, 
				width = 8, height = 6, units = "in")

#plot the global slope test

DivN_Simps_AveTechReps_GlobalSlope <-  ggplot(combined_simps_aveTechReps,
																							aes(x = Day, y = Ave_estimate, 
																									ymin  = ave_lower, ymax = ave_upper,
																									color = as.factor(Day), shape = Elk)) +
																				geom_pointrange(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0)) +
																				scale_shape_manual(values=c(1, 2, 0, 5))+																	     	
																				ylab("Ave. Simpson div. estimate") + xlab("Day") + 
																				labs(color = "Day", shape = "Elk", subtitle = "All Elk: betta model = Simpson ~ Day | Elk)") +
																				scale_y_continuous(limits = c(0.060, 0.09),
																													 expand = c(0, 0)) +
																				coord_cartesian(ylim=c(0.06, 0.085))+
																				scale_x_continuous(limits = c(-1,15), breaks = c(0,1,3,7,14)) +
																				geom_abline(intercept = DivN_test_Simps_DayFixedplusElkRandom$table[1], 
																										slope = DivN_test_Simps_DayFixedplusElkRandom$table[2], 
																										color = "darkgrey") +
																				geom_text(x = -0.5, y = 0.061, 
																									label = paste("p-value =",0.0001), 
																									color = "darkgrey", hjust = 0, vjust = 0) +
																				theme_clean()+ panel_border() +
																				theme(plot.subtitle = element_text(size = 11.5),
																							legend.background = element_blank(),
																							plot.background = element_blank(),
																							# axis.text.x = element_blank(), 
																							axis.ticks.x = element_blank()) 
DivN_Simps_AveTechReps_GlobalSlope

ggsave2(filename = "r_output/visualizations/DivN_Simpson_elk_aveReps_GlobalSlope.pdf", 
				plot = DivN_Simps_AveTechReps_GlobalSlope, 
				width = 4, height = 3, units = "in")

Simps_plotPair3 <- ggpubr::ggarrange(Div_simps_gg,
																		 DivN_Simps_AveTechReps_slopes, 
																		 DivN_Simps_AveTechReps_GlobalSlope,
																		 nrow = 3, labels = c("A.", "B.", "C."),
																		 legend.grob = get_legend(DivN_Simps_AveTechReps_GlobalSlope),
																		 legend = "right")

Simps_plotPair3
ggsave2(filename = "r_output/visualizations/divN_Shann_pair3.pdf", 
				plot = Shann_plotPair3, 
				width = 8, height = 8, units = "in")



#### Combine richness, Shann, Simps plots ####

All_slopes_gg <- ggpubr::ggarrange(ba_richness_AveTechReps_slopes + theme(legend.position = "none"),
																		 DivN_Simps_AveTechReps_slopes + theme(legend.position = "none"),
																		 DivN_Shann_AveTechReps_slopes + theme(legend.position = "none"), 
																	   align = "hv",
																		 nrow = 3, labels = c("A.", "B.", "C."),
																		 common.legend = F)

All_slopes_gg

All_Global_slopes_gg <- ggpubr::ggarrange(ba_richness_AveTechReps_GlobalSlope+ 
																						labs(subtitle = "Richness") +
																						theme(plot.subtitle = element_text(size = 11), 
																									legend.position = "none"),
																					 DivN_Shann_AveTechReps_GlobalSlope+ 
																						labs(subtitle = "Shannon") +
																						theme(plot.subtitle = element_text(size = 11), 
																									legend.position = "none"),
																					 DivN_Simps_AveTechReps_GlobalSlope+ 
																					 	labs(subtitle = "Simpson") +
																					 	theme(plot.subtitle = element_text(size = 11), 
																					 				legend.position = "none"), 
																					 align = "hv",
																					 nrow = 1, labels = c("D.", "", ""))

All_Global_slopes_gg

combined_slopes_gg <- ggpubr::ggarrange(All_slopes_gg + theme(legend.position = "none"), align = "hv",
																				All_Global_slopes_gg, nrow = 2, heights = c(1,1/3),
																				legend.grob = get_legend(DivN_Simps_AveTechReps_GlobalSlope), 
																				legend = "right")

combined_slopes_gg

ggsave2(filename = "r_output/visualizations/combined_div_slopes.pdf", 
				plot = combined_slopes_gg, 
				width = 9, height = 9, units = "in")


#### Beta-Div test with DivNet ####
			#Tax glom to genus 
			genus_data.ps <- tax_glom(elk_rdp.ps, taxrank = "Genus")
			genus_data.ps@sam_data$Elk <- as.factor(as.numeric(as.character(genus_data.ps@sam_data$Elk)))
			# genus_data.ps@sam_data$Day # factor or num will change
			#when satisfied with models increase the n_boot to 10k
			
		##### Model ~Day+Elk factors #####
			##all variables will be compared using the same slope and diff intercepts
			# day is a factor here so
				genus_data.ps@sam_data$Day <- as.factor(as.numeric(genus_data.ps@sam_data$Day))
							class(genus_data.ps@sam_data$Day) #want factors
							class(genus_data.ps@sam_data$Elk)
							class(genus_data.ps@sam_data$Rep)
				
				divnet_genus <-  divnet(W = genus_data.ps,
																ncores = 20,
																formula = ~Day+Elk,
																variance ="none")
				
				bc_test <- testBetaDiversity(divnet_genus, 
																		 h0 = "bray-curtis",
																		 n_boot = 10000, 
																		 sample_specimen_matrix = divnet_genus$X,
																		 groups = genus_data.ps@sam_data$Day)
				bc_test$p_value #p=0.049 with 10k
				
				bc_test_Elk <- testBetaDiversity(divnet_genus, 
																		 h0 = "bray-curtis",
																		 n_boot = 1000,
																		 sample_specimen_matrix = divnet_genus$X,
																		 groups = genus_data.ps@sam_data$Elk)
				bc_test_Elk$p_value #p=0.034 with 1K
	
			# alternative distances
				set.seed(3532)
				euc_test <- testBetaDiversity(dv = divnet_genus, 
																			h0 = "euclidean",
																			n_boot = 1000,
																			sample_specimen_matrix = divnet_genus$X,
																			groups = genus_data.ps@sam_data$Day)
				euc_test$p_value #0.002 with 1k boots
				
				set.seed(3423)
				ait_test <- testBetaDiversity(dv = divnet_genus, 
																			sample_specimen_matrix = divnet_genus$X,
																			h0 = "aitchison",
																			n_boot = 1000,
																			groups = genus_data.ps@sam_data$Day)
				ait_test$p_value #0.002
		##### Model ~Day+Elk+Rep cont #### 
			#all variables will be compared using the same slope and diff intercepts
				
				genus_data.ps@sam_data$Day <- as.numeric(as.character(genus_data.ps@sam_data$Day))
				class(genus_data.ps@sam_data$Day)
				divnet_genusDER <-  divnet(W = genus_data.ps,
																	 ncores = 20,
																	 formula = ~Day+Elk+Rep,
																	 variance ="none")
				
				bc_testDER <- testBetaDiversity(divnet_genusDbyE, 
																				h0 = "bray-curtis",
																				n_boot = 1000, 
																				sample_specimen_matrix = divnet_genusDER$X,
																				groups = genus_data.ps@sam_data$Day)
				bc_testDER$p_value #p=0.017 with 1k
				
				bc_test_DER_elk <- testBetaDiversity(divnet_genusDER, 
																						 h0 = "bray-curtis",
																						 n_boot = 1000,
																						 sample_specimen_matrix = divnet_genusDER$X,
																						 groups = genus_data.ps@sam_data$Elk)
				bc_test_DER_elk$p_value #p=0.034 with 1K
				
				bc_test_DER_rep <- testBetaDiversity(divnet_genusDER, 
																						 h0 = "bray-curtis",
																						 n_boot = 1000,
																						 sample_specimen_matrix = divnet_genusDER$X,
																						 groups = genus_data.ps@sam_data$Rep)
				bc_test_DER_rep$p_value #p=0.034 with 1K
				
				# alternative distances
				euc_test_DER <- testBetaDiversity(divnet_genusDER, 
																					h0 = "euclidean",
																					n_boot = 1000,
																					sample_specimen_matrix = divnet_genusDER$X,
																					groups = genus_data.ps@sam_data$Day)
				euc_test_DER$p_value #0.002 with 1k boots
				
				ait_test_DER <- testBetaDiversity(dv = divnet_genusDER, 
																					sample_specimen_matrix = divnet_genusDER$X,
																					h0 = "aitchison",
																					n_boot = 1000,
																					groups = genus_data.ps@sam_data$Day)
				ait_test_DER$p_value #0.002
		##### Model ~Day*Elk     cont #### 
			# all Elk will be compared to 1 Day slope/interc (cont)
			# using different slopes and diff intercepts each
				vignette("beta_diversity", package = "DivNet")
				genus_data.ps@sam_data$Day <- as.numeric(as.character(genus_data.ps@sam_data$Day))
				genus_data.ps@sam_data$Elk <-as.factor(as.character(genus_data.ps@sam_data$Elk))
				class(genus_data.ps@sam_data$Day)
				
				
				divnet_Day_cont_xE <-  divnet(W = genus_data.ps,
																ncores = 30,
																formula = ~Day*Elk,
																variance = "none")
				#test

				bc_test_Day_con_xE_DAY <- testBetaDiversity(divnet_Day_cont_xE, 
																		 h0 = "bray-curtis",
																		 n_boot = 10000, 
																		 sample_specimen_matrix = divnet_Day_cont_xE$X,
																		 groups =as.factor(genus_data.ps@sam_data$Day))
																		 #as.factor required
				bc_test_Day_con_xE_DAY$p_value #p=0.1805 with 10k
				
				#try aitchison
				aitch_genusDcont_xE <- DivNet::fit_aitchison(W = genus_data.ps@otu_table,
																											ncores = 30, X = divnet_Day_cont_xE$X)
				testBetaDiversity()
				aitch_test_genusDcont_xE <- testBetaDiversity(divnet_Day_cont_xE, 
																										h0 = "aitchison",
																										n_boot = 10, 
																										sample_specimen_matrix = divnet_Day_cont_xE$X,
																										groups =as.factor(genus_data.ps@sam_data$Day))
				#as.factor required
				aitch_test_genusDcont_xE$p_value #p=0.2 with 10k
				
	  ##### try with a model that accounts for technical replicates ######
				# https://rdrr.io/github/adw96/DivNet/f/vignettes/beta_diversity.Rmd
				genus_data.ps@sam_data$SampleID
				unique_specimens <- unique(B2$specimen)
				
				# create "specimen matrix" to be used as design matrix X in DivNet
				specimen_matrix <- do.call(cbind,lapply(unique_specimens,function(x) 
					as.numeric(B2$specimen ==x)))
				specimen_df <- as.data.frame(specimen_matrix) %>%
					(function(x){ rownames(x) <- B2$Bioinformatics.ID; return(x)})
				# The counts data is available in DivNet as a phyloseq object called B2_phyloseq. For this analysis, we will aggregate to the species level using the taxonomic information MBQC has provided along with OTU counts.
				
				data("B2_phyloseq")
				mbqc_species <- tax_glom(B2_phyloseq, taxrank="species")
				# Now we are ready to fit a DivNet model! Note again that we are not estimating any variances in the initial fit.
				
				set.seed(43)
				
				mbqc_div <- divnet(W = mbqc_species,
													 X = as.matrix(specimen_df),
													 variance = "none")
				# Now using testBetaDiversity(), we can conduct a hypothesis test via a nonparametric cluster bootstrap. Here, we have decided to use Aitchison distance, so we are testing a hypothesis of equality of mean centered-log-ratio-transformed relative abundance across groups (colorectal cancer vs. control).
				
				set.seed(14321)
				ait_test_np <- testBetaDiversity(dv = mbqc_div,
																				 h0 = "aitchison",
																				 groups = B2$health_status,
																				 sample_specimen_matrix = as.matrix(specimen_df),
																				 n_boot = 1000)
				
				ait_test_np$p_value
				
		##### Model ~Day*Elk+Rep cont #### 
				# all Elk will be compared to 1 Day slope/interc (cont)
				# using different slopes and diff intercepts for each
				# additional Rep variable with same slope and diff intercept
				genus_data.ps@sam_data$Day <- as.numeric(as.character(genus_data.ps@sam_data$Day))
				class(genus_data.ps@sam_data$Day)
				
				divnet_Day_cont_xER <-  divnet(W = genus_data.ps,
																			ncores = 20,
																			formula = ~Day*Elk+Rep,
																			variance ="none")
				#test
				bc_test_Day_con_xER_DAY <- testBetaDiversity(divnet_Day_cont_xER, 
																										h0 = "bray-curtis",
																										n_boot = 10000, 
																										sample_specimen_matrix = divnet_Day_cont_xER$X,
																										groups =as.factor(genus_data.ps@sam_data$Day)) 
																										#as.fact required
				bc_test_Day_con_xER_DAY$p_value #p=0.2982 with 10k
				
				bc_test_Day_con_xER_ELK <- testBetaDiversity(divnet_Day_cont_xER, 
																										 h0 = "bray-curtis",
																										 n_boot = 1000, 
																										 sample_specimen_matrix = divnet_Day_cont_xER$X,
																										 groups =genus_data.ps@sam_data$Elk)
				bc_test_Day_con_xER_ELK$p_value #p=0.128 with 10k
				
				bc_test_Day_con_xER_REP <- testBetaDiversity(divnet_Day_cont_xER, 
																										 h0 = "bray-curtis",
																										 n_boot = 1000, 
																										 sample_specimen_matrix = divnet_Day_cont_xER$X,
																										 groups =genus_data.ps@sam_data$Rep)
				bc_test_Day_con_xER_REP$p_value #p=0.099 with 10k
		
		##### {skip} plot DivNet variance estimates ####
			length(DivNet::simplifyBeta(divnet_Day_cont_xE, genus_data.ps, "bray-curtis", "Day")$Covar1)
			length(DivNet::simplifyBeta(divnet_Day_cont_xER, genus_data.ps, "bray-curtis", "Elk")$Covar1)
			
			length(DivNet::simplifyBeta(divnet_Day_cont_xE, genus_data.ps, "bray-curtis", "Day")$Covar1)
			length(DivNet::simplifyBeta(divnet_Day_cont_xE, genus_data.ps, "bray-curtis", "Elk")$Covar1)
			class(divnet_Day_cont_xE)
			plot(divnet_Day_cont_xE$fitted_z[1:200], divnet_Day_cont_xE$`bray-curtis`[1:200])
			head(divnet_Day_cont_xE$`bray-curtis`)
			# DivNet::simplifyBeta(divnet_Day_cont_xER, genus_data.ps, "aitchison",   "Elk")
			# DivNet::simplifyBeta(divnet_Day_cont_xER, genus_data.ps, "aitchison",   "Day")
			# vsd <- DESeq2::vst()
			library("vegan")
			bray_dist_py <- phyloseq::ordinate(genus_data.ps, "PCoA", distance = "morisita")
			plot_ordination(genus_data.ps, bray_dist_py, color = "Day") 

			barplot(bray_dist_py$values$Relative_eig[1:10])
			Axis1.percent <- bray_dist_py$values$Relative_eig[[1]] * 100
			Axis2.percent <- bray_dist_py$values$Relative_eig[[2]] * 100
			
			bray_dist <- divnet_Day_cont_xE$`bray-curtis`
			euc_dist  <- divnet_Day_cont_xE$euclidean
			pcoa <- ape::pcoa(euc_dist)
			barplot(pcoa$values$Relative_eig[1:10])
			Axis1.percent <- pcoa$values$Relative_eig[[1]] * 100
			Axis2.percent <- pcoa$values$Relative_eig[[2]] * 100
			Axis3.percent <- pcoa$values$Relative_eig[[3]] * 100
			Axis4.percent <- pcoa$values$Relative_eig[[4]] * 100

			pcoa.axis.data <- data.frame(sample = rownames(pcoa$vectors),
							Ax1 = pcoa$vectors[, 1],
							Ax2 = pcoa$vectors[, 2],
							Ax3 = pcoa$vectors[, 3],
							Ax4 = pcoa$vectors[, 4]
						)
			
			pcoa.axis.data_py <- data.frame(sample = rownames(bray_dist_py$vectors),
							Ax1 = bray_dist_py$vectors[, 1],
							Ax2 = bray_dist_py$vectors[, 2],
							Ax3 = bray_dist_py$vectors[, 3],
							Ax4 = bray_dist_py$vectors[, 4]
						)
			
			pcoa.axis.data <- as_tibble(pcoa.axis.data) %>% 
																left_join(genus_data.ps@sam_data, 
																	by = c("sample" = "Description"))
			
			pcoa.axis.data_py <- as_tibble(pcoa.axis.data_py) %>% 
				left_join(genus_data.ps@sam_data, 
									by = c("sample" = "Description"))
			
			
			ggplot(data = pcoa.axis.data_py, aes(x = Ax1, y = Ax2)) +
				geom_point(aes(color = as.factor(Day), group = Elk)) +
				xlab(paste("PCOA1 - ", round(Axis1.percent, 2), "%", sep = "")) +
				ylab(paste("PCOA2 - ", round(Axis2.percent, 2), "%", sep = "")) +
				theme_bw() +
				theme(
					axis.title.x = element_text(face = "bold", size = 15),
					axis.title.y = element_text(face = "bold", size = 15),
					axis.text = element_text(face = "bold", size = 10),
					legend.title = element_text(face = "bold", size = 10),
					legend.text = element_text(face = "bold", size = 10),
					legend.key.size = unit(1, 'lines')
				) +
				stat_ellipse(aes(group = Elk), level = 0.95) + 
				ggtitle("DivNet Bray-Curtis Distances PCOA")
			
			
			ggplot(data = pcoa.axis.data_py, aes(x = Ax2, y = Ax1)) +
				geom_point(aes(color = as.factor(Elk), group = Elk)) +
				xlab(paste("PCOA1 - ", round(Axis1.percent, 2), "%", sep = "")) +
				ylab(paste("PCOA2 - ", round(Axis2.percent, 2), "%", sep = "")) +
				theme_bw() +
				theme(
					axis.title.x = element_text(face = "bold", size = 15),
					axis.title.y = element_text(face = "bold", size = 15),
					axis.text = element_text(face = "bold", size = 10),
					legend.title = element_text(face = "bold", size = 10),
					legend.text = element_text(face = "bold", size = 10),
					legend.key.size = unit(1, 'lines')
				) +
				stat_ellipse(aes(group = Elk), level = 0.95) + 
				ggtitle("DivNet Bray-Curtis Distances PCOA")
			
			####
			
			beta_div_DAY_bray_plot <- simplifyBeta(divnet_Day_cont_xE, genus_data.ps, "bray-curtis", "Day") %>%
				ggplot(aes(x = interaction(Covar1, Covar2), 
									 y = beta_est,
									 col = interaction(Covar1, Covar2))) +
								geom_point() +
								# geom_linerange(aes(ymin = lower, ymax = upper)) + 
								theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
								xlab("") + ylab("Estimates of Bray-Curtis distance")
			
			beta_div_ELK_bray_plot <- simplifyBeta(divnet_Day_cont_xER, genus_data.ps, "bray-curtis", "Elk") %>%
				ggplot(aes(x = interaction(Covar1, Covar2), 
									 y = beta_est,
									 col = interaction(Covar1, Covar2))) +
				geom_point() +
				# geom_linerange(aes(ymin = lower, ymax = upper)) + 
				theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
				xlab("") + ylab("Estimates of Bray-Curtis distance")
			
			
			beta_div_euc_plot <- simplifyBeta(divnet_genus_pred, genus_data.ps, "euclidean", "predator") %>%
				ggplot(aes(x = interaction(Covar1, Covar2), 
									 y = beta_est,
									 col = interaction(Covar1, Covar2))) +
								geom_point() +
								geom_linerange(aes(ymin = lower, ymax = upper)) + 
								theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
								xlab("") + ylab("Estimates of Euclidean distance")

			beta_div_ait_plot <- simplifyBeta(dv = divnet_genus_pred, genus_data.ps, "aitchison", "predator") %>%
								ggplot(aes(x = interaction(Covar1, Covar2), 
													y = beta_est,
													col = interaction(Covar1, Covar2))) +
								geom_point() +
								geom_linerange(aes(ymin = lower, ymax = upper)) + 
								theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
								xlab("") + ylab("Estimates of Aitchison distance")
			
			beta_div_bray_plot <- beta_div_bray_plot + theme_bw() + 
				theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")
			ggsave2(filename = "DivNet_model_beta_div_bray_plot.pdf", 
							plot = beta_div_bray_plot,
							path = "output/plots", height = 4, width = 3)


#### visualise Beta-Div - PCA & Adonis #####
			# When the distance metric is Euclidean, PCoA is equivalent to Principal Components Analysis
			# generating and visualizing the PCA with phyloseq
			# RDA without constraining variables is PCA
			# Unfortunately, a linear assumption causes PCA to suffer from a serious problem, 
			# the horseshoe effect, which makes it unsuitable for most ecological data sets (Gauch 1982). 
			# The PCA solution is often distorted into a horseshoe shape (with the toe either up or down)
			#  if beta diversity is moderate to high. The horseshoe can appear even if there is an 
			#  important secondary gradient.
			elk_rdp_VST.ps@sam_data$Day <- as.factor(elk_rdp_VST.ps@sam_data$Day)
			vst_pca <- ordinate(elk_rdp_VST.ps, method= "RDA", distance="euclidean")
			
			elk_rdp_ASV_adujsted.ps@sam_data$Day <- as.factor(elk_rdp_VST.ps@sam_data$Day)
			ancom_pca <- ordinate(elk_rdp_ASV_adujsted.ps, method= "RDA", distance="euclidean")
			
			# ggsci::scale_color_aaas(alpha = 0.8)
			# ggsci::pal_aaas(palette = c("default"), alpha = 0.8)
			# Day colors
			# "0"="#9E0142",	"1"="#D53E4F" ,	"3"="#F46D43",	
			# "7"="#66C2A5",	"14"="#3288BD",	"Taxa"="#5E4FA2"

			#plot scree
			plot_ordination(elk_rdp_VST.ps,            vst_pca, axes = c(1,2), type = "scree") #1-4
			plot_ordination(elk_rdp_ASV_adujsted.ps, ancom_pca, axes = c(1,2), type = "scree") #1-4
			
			# plot biplot
			plot_ordination(elk_rdp_VST.ps, vst_pca, color= "Day", axes = c(1,2), type = "split") + 
				geom_point(size=1) + labs(col="Day") + 
				# geom_text(aes(label=elk_rdp_VST.ps@sam_data$Elk, hjust=0.3, vjust=-0.4)) +
				ggtitle("PCA - Aitchison distance") +
				# coord_fixed(sqrt(vst_eigen_vals[2]/vst_eigen_vals[1])) + ggtitle("PCA") + 
				scale_color_manual(values=c("0"="#9E0142",	"1"="#D53E4F" ,	"3"="#F46D43",	
																		"7"="#66C2A5",	"14"="#3288BD",	"Taxa"="#5E4FA2")) +
				theme_bw() #+ theme(legend.position="none")
			
			plot_ordination(elk_rdp_ASV_adujsted.ps, ancom_pca, color= "Day", axes = c(2,3), type = "samples") + 
				geom_point(size=1) + labs(col="Day") +
				geom_text(aes(label=elk_rdp_VST.ps@sam_data$Elk, hjust=0.3, vjust=-0.4), size = 3) +
				ggtitle("PCA - Aitchison distance") +
				scale_color_manual(values=c("0"="#9E0142",	"1"="#D53E4F" ,	"3"="#F46D43",	
																		"7"="#66C2A5",	"14"="#3288BD",	"Taxa"="#5E4FA2")) +
				theme_bw() #+ theme(legend.position="none")
			
			######with another package#####
			# elk_rdp_ASV_adujsted.ps@sam_data$Elk <- as.factor(paste0("Elk.",elk_rdp_ASV_adujsted.ps@sam_data$Elk))
			# # elk_rdp_ASV_adujsted.ps@sam_data$Day <- as.factor(paste0("Day.",elk_rdp_ASV_adujsted.ps@sam_data$Day))
			# 
			# elk_rdp_ASV_adujsted.ps@sam_data$Day <- factor(elk_rdp_ASV_adujsted.ps@sam_data$Day, 
			# 																							 levels = c( "Day.0","Day.1","Day.3","Day.7","Day.14"))
			# 	
			pcares <- get_pca(obj=elk_rdp_ASV_adujsted.ps, method=NULL)
			# Visulizing the result
			pcaplot1 <- ggordpoint(obj=pcares, biplot=TRUE, speciesannot=TRUE,
														 factorNames=c("Elk"), ellipse=TRUE) +
				ggsci::scale_color_aaas(alpha = 0.8) +
				ggsci::scale_fill_aaas(alpha = 0.8)
			
			# pc = c(1, 3) to show the first and third principal components.
			pcaplot2 <- ggordpoint(obj=pcares, biplot=TRUE, speciesannot=TRUE,
														 factorNames=c("Elk"), ellipse=TRUE, pc = c(1, 3)) +
				ggsci::scale_color_aaas(alpha = 0.8) +
				ggsci::scale_fill_aaas(alpha = 0.8)
		pcaplot1 | pcaplot2
	
		# try with Day labels
			pcaplot3 <- ggordpoint(obj=pcares, biplot=TRUE, speciesannot=TRUE,
														 factorNames=c("Day"), ellipse=F) +
				ggsci::scale_color_aaas(alpha = 0.8) +
				ggsci::scale_fill_aaas(alpha = 0.8)
			
			# pc = c(1, 3) to show the first and third principal components.
			pcaplot4 <- ggordpoint(obj=pcares, biplot=TRUE, speciesannot=TRUE,
														 factorNames=c("Day"), ellipse=F, pc = c(1, 3)) +
				ggsci::scale_color_aaas(alpha = 0.8) +
				ggsci::scale_fill_aaas(alpha = 0.8)
			
			pcaplot3 | pcaplot4
##### RDA model####
			# install.packages("Rfast")
			devtools::install_github("andjar/ALASCA",ref="main")
			library("ALASCA");library(vegan)
			# http://ordination.okstate.edu/overview.htm
			
			#RDA or CCA
			# The choice depends on (1) the nature of your data (hetero or homogeneous) and 
			# (2) the objective (linear or unimodal link). if your data are homogeneous a method with 
			# a linear link between (variables / species), IS appropriate (rda), 
			# whereas if they are heterogeneous you should apply an ordination with a unimodal link (cca).
			# the application of a Detrended Correspondence Analysis (DCA) a priori will allow to choose 
			# the appropriate ordination where if the length of the first axis is <5 
			# unimodal ordination methods are preferable and vice versa
			
			# options
			# elk_rdp_ASV_iQlr.ps
			# elk_rdp_genus_iQlr.ps
			# elk_rdp_spp_iQlr.ps
			# elk_rdp_ASV_CLR.ps
			# elk_rdp_day0_CLR.ps
			# elk_rdp_genus_CLR.ps
			# elk_rdp_spp_CLR.ps
			
			library(vegan)
			day0CLR_aitch_dist  <- dist(elk_rdp_day0_CLR.ps@otu_table)
			DCA_ord             <- vegan::decorana(day0CLR_aitch_dist)
			plot(DCA_ord, type = "points", display = "both") #sites = circles, spp = +
			# The numeric scale on the axis is not very useful for the interpretation 
			# (an exception for this is DCA, in which the scales are in units of beta diversity).
			# axis <5 use CCA
			# scatterplot3d::
				# PCoA suffers from a number of flaws, in particular the arch effect (discussed later in the context of PCA and CA). 
				# These flaws stem, in part, from the fact that PCoA maximizes a linear correlation. 
				# Nonmetric Multidimensional Scaling (NMDS) rectifies this by maximizing the rank order correlation. 
				# The algorithm in brief outline, proceeds as follows:
				# ‘stress’ (mismatch between the rank order of distances in the data, and the rank order of distances in the ordination)
			# MDS_ord <- metaMDS(comm = elk_rdp_ASV_adujsted.ps@otu_table, distance = "euclidean")
			# MDS_ord2 <- metaMDS(k = 4, comm = elk_rdp_ASV_adujsted.ps@otu_table, distance = "euclidean", 
			# 										previous.best=MDS_ord)
			# plot(MDS_ord)
			
			
			### prc
			prc_mod <- vegan::prc(response = elk_rdp_day0_CLR.ps@otu_table, 
										 treatment = as.factor(elk_rdp_day0_CLR.ps@sam_data$Elk), 
										 time = as.factor(elk_rdp_day0_CLR.ps@sam_data$Day))
			summary(prc_mod)
			logabu <- colSums(abs(t(elk_rdp_ASV_adujsted.ps@otu_table)))
			summary(logabu)
			par()
			plot(prc_mod)
			ctrl <- how(blocks = elk_rdp_day0_CLR.ps@sam_data$Elk,
									plots = Plots(strata = elk_rdp_day0_CLR.ps@sam_data$Rep,
																type = "free"),
									within = Within(type = "series",mirror = F), nperm = 999)
			anova(prc_mod, permutations = ctrl, first= T)
			# Permutation test for rda under reduced model
			# Blocks:  elk_rdp_day0_CLR.ps@sam_data$Elk 
			# Plots: elk_rdp_day0_CLR.ps@sam_data$Rep, plot permutation: free
			# Permutation: series
			# Number of permutations: 999
			# 
			# Df Variance      F Pr(>F)  
			# RDA1      1   2922.2 13.083  0.043 *
			# 	Residual 40   8934.6                
			# ---
			# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1            

			
#### ALASCA LMM of genera #####
			# remove.packages("ALASCA")
			# devtools::install_github("andjar/ALASCA",ref="main")
			library(magrittr)
			library(tidyverse)
			library(ALASCA)
			#prebuilt
			# elk_rdp_ASV_iQlr.ps   
			# elk_rdp_genus_iQlr.ps 
			# elk_rdp_spp_iQlr.ps   
			# elk_rdp_ASV_CLR.ps    
			# elk_rdp_day0_CLR.ps   
			# elk_rdp_genus_CLR.ps  
			# elk_rdp_spp_CLR.ps    
			
			elk_rdp_ASV_CLR.ps    <-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_ASV_CLR.ps.rds")
			elk_rdp_spp_CLR.ps    <-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_spp_CLR.ps.rds")
			elk_rdp_ASV_iQlr.ps   <-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_ASV_iQlr.ps.rds")
			elk_rdp_spp_iQlr.ps   <-  readRDS(file = "r_output/intermediate_rds_files/elk_rdp_spp_iQlr.ps.rds")
			
			
			# species as columns
			# df1.CLRday0 <- as.data.frame(elk_rdp_day0_CLR.ps@otu_table)
			df2.CLRasv <- as.data.frame(elk_rdp_ASV_CLR.ps@otu_table)
			# df3.CLRgen <- as.data.frame(elk_rdp_genus_CLR.ps@otu_table)
			df4.CLRspp <- as.data.frame(elk_rdp_spp_CLR.ps@otu_table)
			
			df1.iQlrASV <- as.data.frame(elk_rdp_ASV_iQlr.ps@otu_table)
			# df2.iQlrGen  <- as.data.frame(elk_rdp_genus_iQlr.ps@otu_table)
			df3.iQlrSpp  <- as.data.frame(elk_rdp_spp_iQlr.ps@otu_table)
			
			# length(colnames(df1.CLRday0)) #5177
			# length(colnames(df2.CLRasv))  #5177
			# length(colnames(df3.CLRgen))  #158
			# length(colnames(df4.CLRspp))  #162
			# 
			# length(colnames(df1.iQlrASV)) #5177
			# length(colnames(df2.iQlrGen))  #158
			# length(colnames(df3.iQlrSpp))  #162
			
			# df1.CLRday0 <- cbind(elk_rdp_day0_CLR.ps@sam_data,  df1.CLRday0)
			df2.CLRasv  <- cbind(elk_rdp_spp_CLR.ps@sam_data,   df2.CLRasv)
			# df3.CLRgen  <- cbind(elk_rdp_genus_CLR.ps@sam_data, df3.CLRgen)
			df4.CLRspp  <- cbind(elk_rdp_spp_CLR.ps@sam_data,   df4.CLRspp)
			
			df1.iQlrASV <- cbind(elk_rdp_spp_CLR.ps@sam_data,  df1.iQlrASV)
			# df2.iQlrGen <- cbind(elk_rdp_ASV_CLR.ps@sam_data,   df2.iQlrGen)
			df3.iQlrSpp <- cbind(elk_rdp_spp_CLR.ps@sam_data, df3.iQlrSpp)
			
			rename_ALASCA_to_long <- function(ALASCA_df, metadata_col_length = 9) {
									require(dplyr)
									require(reshape2
													)
									# length <- length(colnames(ALASCA_df))
									df_rename <- ALASCA_df %>%  mutate(Elk = as.factor(Elk)) %>% 
																		mutate(Elk = recode(Elk, "1" = "Elk_b", 
																														 "2" = "Elk_c",
																														 "3" = "Elk_d", 
																														 "4" = "Elk_a")) %>%
																		mutate(Elk = factor(Elk, 
																												levels = c("Elk_a", "Elk_b", 
																																	 "Elk_c", "Elk_d"))) %>%
																		mutate(ID = as.factor(Description), 
																								 time = as.factor(Day), 
																									Rep = as.factor(Rep),
																								 time_cont = as.numeric(as.character(Day)),
																								 group = as.factor(Elk), 
																								 group_time = as.factor(paste0(Elk, ".", Day)), 
																								 sub_group_elk = as.factor(paste0(Elk,".",Rep)),			 
																								 sub_group_day = as.factor(paste0(Day,".",Rep))) 
									new.length <- length(colnames(df_rename))
									df_rename <- df_rename %>%
																		select(ID, Rep, Elk, group, group_time, sub_group_elk, sub_group_day, time, time_cont, (metadata_col_length+1):all_of(new.length)) %>%
																					 arrange(Rep, time, group) %>%  
																					 mutate(ID=factor(ID, levels=ID))
									long <- reshape2::melt(df_rename, id.vars = colnames(df_rename)[1:metadata_col_length])
									print(paste("col length", new.length, "-------- new dims", dim(long)[1]))
									return(long)
			}
			
			
			# Test nesting
				
			# we have change the factors to be implicitly nested by making hybrid variables sub_group, group_time
			# this means each sub factor is unique and cannot be found in other levels (elk.a.day0:rep1.a.0)

			
			
			# df1.CLRday0_long <- rename_ALASCA_to_long(df1.CLRday0)
			df2.CLRasv_long <- rename_ALASCA_to_long(df2.CLRasv)
			# df3.CLRgen_long <- rename_ALASCA_to_long(df3.CLRgen)
			df4.CLRspp_long <- rename_ALASCA_to_long(df4.CLRspp)
			
			df1.iQlrASV_long <- rename_ALASCA_to_long(df1.iQlrASV)
			# df2.iQlrGen_long <- rename_ALASCA_to_long(df2.iQlrGen)
			df3.iQlrSpp_long <- rename_ALASCA_to_long(df3.iQlrSpp)
			
			library(ggplot2)
			ggplot(data = df3.iQlrSpp_long, aes(x = factor(time), y = value, color = group)) +
				geom_boxplot() +
				facet_wrap(~variable, scales = "free_y") +
				theme(legend.position = "bottom") +
				labs(x = "Time", y = "Value", color = "group")
				
			do.call(
				ggpubr::ggarrange,
				c(ALASCA::plotParts(df3.iQlrSpp_long, participantColumn = "group", valueColumn = "value", addSmooth = NA), 
					common.legend = TRUE, legend = "bottom")
			)
			
			
			
		##### run and compare LMM models ####
# 			library(lme4)      # "golden standard" for mixed-effects modelling in R (no p-values)
# 			library(lmerTest)  # p-values for MEMs based on the Satterthwaite approximation
# 			library(psycho)    # mainly for an "analyze()" function
# 			library(broom)     # for tidy results
# 			library(knitr)     # beautifying tables
# 			library(sjPlot)    # for visualising MEMs
# 			library(effects)   # for visualising MEMs
# 			library(report)    # for describing models
# 			library(emmeans)   # for post-hoc analysis
# 			library(nlme)
# 			devtools::install_github("andjar/ALASCA",ref="main")
#       library(ALASCA)
			
			colnames(df4.CLRspp_long)
			df4.CLRspp_long$sub_group_day
			#linear models
			model.formula.lm      <- value ~ time_cont+ time_cont*group
      model.formula.lm2     <- value ~ time_cont
     
      #random effects
			model.formula2.0   <- value ~ time_cont          + (1|group)            #random effect for elk, intercept only
      model.formula2.1   <- value ~ time_cont          + (1|group) + (1|Rep) #with fixed effect interaction
      model.formula2.2   <- value ~ time_cont * group  + (1|group)           #with fixed effect interaction
      model.formula2.2   <- value ~ time_cont * group  + (1|group) + (1|Rep) #with fixed effect interaction
			model.formula3.0   <- value ~ time_cont          + (1|group/Rep) #nested random
			model.formula3.1   <- value ~ time               + (1|group)
			model.formula3.2   <- value ~ time * group       + (1|Rep)    #nested random
			model.formula3.3   <- value ~ time * group       + (1|group) +(1|Rep) #nested random
			model.formula3.4   <- value ~ time * group       + (1|group/Rep)
			# model.formula4.0   <- value ~ time_cont                      + (time_cont | group) # nested random for time
			# model.formula4.1   <- value ~ time_cont                      + (time_cont | group/sub_group_day) # remove correlation
			# model.formula4.2   <- value ~ time_cont + group              + (time_cont | group/sub_group_day) # remove correlation
			# model.formula4.3   <- value ~ time_cont * group              + (time_cont | group/sub_group_day) # remove correlation
			# model.formula5.0   <- value ~ time_cont + time_cont*group    + (time_cont | group/sub_group_day) # add fixed interaction
			# model.formula5.1   <- value ~ time_cont + time_cont*group    + (time_cont | group/sub_group_day)# remove corr
			
			fm2.0   <- lme4::lmer(model.formula2.0,   data = df4.CLRspp_long)
			fm2.1   <- lme4::lmer(model.formula2.1,   data = df4.CLRspp_long)
			fm2.2   <- lme4::lmer(model.formula2.2,   data = df4.CLRspp_long)
			fm3.0   <- lme4::lmer(model.formula3.0,   data = df4.CLRspp_long)
			fm3.1   <- lme4::lmer(model.formula3.1,   data = df4.CLRspp_long)
			fm3.2   <- lme4::lmer(model.formula3.2,   data = df4.CLRspp_long)
			fm3.3   <- lme4::lmer(model.formula3.3,   data = df4.CLRspp_long)

			# fm4.0   <- lme4::lmer(model.formula4.0,   data = df3.iQlrSpp_long)
			# fm4.1   <- lme4::lmer(model.formula4.1, data = df3.iQlrSpp_long)
			# fm4.2   <- lme4::lmer(model.formula4.2, data = df3.iQlrSpp_long)
			# fm4.3   <- lme4::lmer(model.formula4.3, data = df3.iQlrSpp_long)
			# fm5.0   <- lme4::lmer(model.formula5.0,   data = df3.iQlrSpp_long)
			# fm5.1   <- lme4::lmer(model.formula5.1, data = df3.iQlrSpp_long)
		
		summary(fm1)
		summary(fm2)	
		summary(fm3)	
		summary(fm3.2)	
		summary(fm3.3)	
		summary(fm4)	
		summary(fm4.1)	
		summary(fm4.2)	
		summary(fm5)	
		summary(fm5.2)	
		# random effects should be compared first in anova
		anova(fm3.3, fm3.2, fm3.1, fm3.0, fm2.1, fm2.2, fm2.0, REML= T)
		
					# refitting model(s) with ML (instead of REML)
					# Data: df3.iQlrSpp_long
					# Models:
					# 	fm2.0: value ~ time_cont + (1 | group)
					# fm3.0: value ~ time_cont + (1 | group/sub_group_day)
					# fm4.0: value ~ time_cont + (time_cont | group)
					# fm2.1: value ~ time_cont + group + (1 | group)
					# fm3.1: value ~ time_cont + group + (1 | group/sub_group_day)
					# fm4.1: value ~ time_cont + (time_cont | group/sub_group_day)
					# fm3.2: value ~ time_cont * group + (1 | group/sub_group_day)
					# fm4.2: value ~ time_cont + group + (time_cont | group/sub_group_day)
					# fm5.1: value ~ time_cont + time_cont * group + (time_cont | group/sub_group_day)
					# fm5.0: value ~ time_cont + time_cont * group + (time_cont | group/sub_group_day)
					# fm4.3: value ~ time_cont * group + (time_cont | group/sub_group_day)
					# 
					#       npar   AIC   BIC logLik deviance   Chisq Df Pr(>Chisq)    
					# fm2.0    4 56778 56806 -28385    56770                          
					# fm3.0    5 56767 56803 -28378    56757 12.5751  1  0.0003909 ***
					# fm4.0    6 56775 56818 -28382    56763  0.0000  1  1.0000000    
					# fm2.1    7 56770 56820 -28378    56756  7.2673  1  0.0070221 ** 
					# fm3.1    8 56761 56819 -28373    56745 10.3117  1  0.0013219 ** 
					# fm4.1    9 56772 56836 -28377    56754  0.0000  1  1.0000000    
					# fm3.2   11 56759 56838 -28368    56737 16.6637  2  0.0002407 ***
					# fm4.2   12 56768 56854 -28372    56744  0.0000  1  1.0000000    
					# fm5.1   15 56767 56875 -28368    56737  6.8947  3  0.0753305 .  
					# fm5.0   15 56767 56875 -28368    56737  0.0000  0               
					# fm4.3   15 56767 56875 -28368    56737  0.0000  0               
					# ---
					# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
		# 
		
		
		
		plot(fm2.2, type = c("p", "smooth"))
		plot(fm2.2, sqrt(abs(resid(.))) ~ fitted(.), type = c("p", "smooth"))
		qqmath(fm2.2, id = 0.05)
		dotplot(ranef(fm2.2, condVar=T))
		plot(allEffects(fm5)) 
		performance::check_model(fm2.2)

		
		
		##### run ALASCA #####
		# library(ALASCA)
		#scale fun
		# scaleFun <- function(df){
		# 	for(i in unique(df$variable)){
		# 		df$value[df$variable == i] <- df$value[df$variable == i]/sd(df$value[df$variable == i & df$time == "TP1"])
		# 	}
		# 	return(df)
		# }
		colnames(df4.CLRspp_long)
		class(df4.CLRspp_long$Rep)
		# model <- value ~ time * group  + (1|ID) #kinda weird singular results
		# lots of errors
		model1 <- value ~ time + (1|Rep) + (1|group) 
	
		#rfast F for only one random intercept, turn off for nested effects
			ASCA.mod.CLRspp <- ALASCA::ALASCA(df4.CLRspp_long, model1, #scaleFun = scaleFun,
																 separateTimeAndGroup = F, forceEqualBaseline = F, # method ="LMM",
																 participantColumn = "group", doDebug = T,useSumCoding = T,
																 stratificationColumn = "ID", 
																 validateRegression = F, nValRuns = 100,
																 # pAdjustMethod =  'hochberg', 
																 # validate = F, validationMethod = "bootstrap", #nValRuns = 10, #1000
																 # validate = T, validationMethod = "loo", nValFold = 59, #error
																 # validate = T, validationMethod = "permutation",  #error
																 useRfast = F, save = F,
																 plot.filetype = "pdf", 
																 plot.xlabel = "Day",plot.grouplabel = "Elk",
																 filename = "ALASCA.CLRspp", 
																 filepath ="r_output/visualizations/ALASCA/ALA_LMM_CLRspp_m2_groupValidation/" )
			# str(PE.mod.tipg)
			summary(ASCA.mod.CLRspp)
		
			# plot(ASCA.mod.CLRspp, tooDense = 5)
			# PE.mod <- flipIt(ASCA.mod.CLRspp)
			# summary(PE.mod)
			# ASCA.mod.CLRspp$RegressionCoefficients
			plotVal(ASCA.mod.CLRspp)
			screeplot(ASCA.mod.CLRspp)
			plot(ASCA.mod.CLRspp, tooDense = 5) #effect = "time" or effect = "group"
			plot(ASCA.mod.CLRspp, component = 2, tooDense = 5)
			plot(ASCA.mod.CLRspp, component = 3, tooDense = 5)
			plot(ASCA.mod.CLRspp, component = 4, tooDense = 5)
			
			
			do.call(
				ggpubr::ggarrange,
				c(plotPred(PE.mod, variable = c("Elk", "Rep")), 
					common.legend = TRUE, legend = "bottom")
			)
			# summary(res.simple$regr.model[[1]])
			summary(PE.mod$RegressionCoefficients)
			head(PE.mod$RegressionCoefficients)
			ALASCA::plotResiduals(PE.mod)#[[1]]
			plot(density(residuals(PE.mod, variable = "Elk")[[1]]), main = "Elk")
			qqnorm(residuals(res.simple, variable = "BMI")[[1]], main = "BMI") 
			qqline(residuals(res.simple, variable = "BMI")[[1]])
			
			
			#####
			model2 <- value ~ time * group + (1|Rep) #shows a strong time signal
			# this model would test each elks effect
			#no errors with either group or ID as validation, but no CI w/ group
			ASCA.mod.CLRspp_m2 <- ALASCA::ALASCA(df4.CLRspp_long, model2, #scaleFun = scaleFun,
																				separateTimeAndGroup = T, forceEqualBaseline = T, method ="LMM",
																				participantColumn = "ID", doDebug = T,
																				validate = F, validateRegression = F, #valreg not available lmm
																				nValRuns = 100,
																				# validationMethod = "bootstrap",
																				# validationMethod = "permutation",
																				# validationMethod = "loo", #nValFold = 7,
																				useRfast = T, save = F, 
																				plot.filetype = "pdf", 
																				plot.xlabel = "Day",plot.grouplabel = "Elk",
																				filename = "ALASCA.CLRspp_m2_IDval", 
																				filepath ="r_output/visualizations/ALASCA/ALA_LMM_CLRspp_m2_IDval/" )
			# str(PE.mod.tipg)
			summary(ASCA.mod.CLRspp_m2)
			# plot(ASCA.mod.CLRspp, tooDense = 5)
			PE.mod2 <- flipIt(ASCA.mod.CLRspp_m2)
			summary(PE.mod2)
			screeplot(PE.mod2)
			plot(ASCA.mod.CLRspp_m2, tooDense = 5) #effect = "time" or effect = "group"
			plot(ASCA.mod.CLRspp_m2, component = 2, tooDense = 5)
			plot(ASCA.mod.CLRspp_m2, component = 3, tooDense = 5)
			plot(ASCA.mod.CLRspp_m2, component = 4, tooDense = 5)
			
			
			#####
			model3 <- value ~ time +	time:group  + (1|ID) #this worked great! no time sig
			#singular errors
			ASCA.mod.CLRspp_m3 <- ALASCA::ALASCA(df4.CLRspp_long, model3, #scaleFun = scaleFun,
																				separateTimeAndGroup = T, forceEqualBaseline = T, method ="LMM",
																				participantColumn = "ID", doDebug = T,
																				validate = T, validateRegression = F, #valreg not available lmm
																				nValRuns = 100,
																				# validationMethod = "bootstrap",
																				# validationMethod = "permutation",
																				# validationMethod = "loo", #nValFold = 7,
																				useRfast = T, save = T, 
																				plot.filetype = "pdf", 
																				plot.xlabel = "Day",plot.grouplabel = "Elk",
																				filename = "ALASCA.CLRspp_m3", 
																				filepath ="r_output/visualizations/ALASCA/ALA_LMM_CLRspp_m3_IDval/" )
			# str(PE.mod.tipg)
			summary(ASCA.mod.CLRspp_m3)
			# plot(ASCA.mod.CLRspp, tooDense = 5)
			PE.mod3 <- flipIt(ASCA.mod.CLRspp_m3)
			summary(PE.mod3)
			screeplot(PE.mod3)
			plot(PE.mod3, tooDense = 5) #effect = "time" or effect = "group"
			plot(PE.mod3, component = 2, tooDense = 5)
			plot(PE.mod3, component = 3, tooDense = 5)
			plot(PE.mod3, component = 4, tooDense = 5)
			
			
			#######
			model4 <- value ~ time +	time*group  + (1|group) #
			#no errors but also no CIs
			ASCA.mod.CLRspp_m4 <- ALASCA::ALASCA(df4.CLRspp_long, model4, #scaleFun = scaleFun,
																					 separateTimeAndGroup = T, forceEqualBaseline = T, method ="LMM",
																					 participantColumn = "ID", doDebug = T,
																					 validate = T, validateRegression = F, #valreg not available lmm
																					 nValRuns = 100,
																					 # validationMethod = "bootstrap",
																					 # validationMethod = "permutation",
																					 # validationMethod = "loo", #nValFold = 7,
																					 useRfast = T, save = T, 
																					 plot.filetype = "pdf", 
																					 plot.xlabel = "Day",plot.grouplabel = "Elk",
																					 filename = "ALASCA.CLRspp_m4", 
																					 filepath ="r_output/visualizations/ALASCA/ALA_LMM_CLRspp_m4_IDval/" )
			# str(PE.mod.tipg)
			summary(ASCA.mod.CLRspp_m4)
			# plot(ASCA.mod.CLRspp, tooDense = 5)
			PE.mod4 <- flipIt(ASCA.mod.CLRspp_m4)
			summary(PE.mod4)
			screeplot(PE.mod4)
			plot(PE.mod4, tooDense = 5) #effect = "time" or effect = "group"
			plot(PE.mod4, component = 2, tooDense = 5)
			plot(PE.mod4, component = 3, tooDense = 5)
			plot(PE.mod4, component = 4, tooDense = 5)
			
			######
			
			model5 <- value ~ time + (1|Rep) + (1|group) 
			model5 <- value ~ Day +  (1|Elk) + (1|Elk:Day) + (1|Rep) + (1|Elk:Rep) + (1|Day:Rep)
			
			#singular errors
			rm(ASCA.mod.CLRspp_m5)
			ASCA.mod.CLRspp_m5 <- ALASCA::ALASCA(df4.CLRspp_long, model5, #scaleFun = scaleFun,
																					 separateTimeAndGroup = F, forceEqualBaseline = T, method ="Limm",
																					 participantColumn = "Day", doDebug = T,
																					 validate = T, validateRegression = F, #valreg not available lmm
																					 nValRuns = 10,
																					 # validationMethod = "bootstrap",
																					 # validationMethod = "permutation",
																					 # validationMethod = "loo", #nValFold = 7,
																					 useRfast = F, save = T, 
																					 plot.filetype = "pdf", 
																					 plot.xlabel = "Day",plot.grouplabel = "Elk",
																					 filename = "ALASCA.CLRspp_m5", 
																					 filepath ="r_output/visualizations/ALASCA/ALA_LMM_CLRspp_m5_IDval/" )
			# str(PE.mod.tipg)
			summary(ASCA.mod.CLRspp_m5)
			# plot(ASCA.mod.CLRspp, tooDense = 5)
			PE.mod5 <- flipIt(ASCA.mod.CLRspp_m5)
			summary(PE.mod5)
			screeplot(PE.mod5)
			plot(PE.mod5, tooDense = 5) #effect = "time" or effect = "group"
			plot(PE.mod5, component = 2, tooDense = 5)
			plot(PE.mod5, component = 3, tooDense = 5)
			plot(PE.mod5, component = 4, tooDense = 5)
			
			######
			donovan <- value ~ time + (1|group) + (1|Rep) #mean of time with change in general interc.
			
			
			top_loadings <- subset(getLoadings(PE.mod)$time, PC == 1)
			top_loadings <- top_loadings[order(top_loadings$loading, decreasing = TRUE),]
			top_loadings[1:5,]
			#validate the model
			ALASCA::validate(PE.mod, participantColumn = "group", validateRegression = FALSE)
			plot(PE.mod, tooDense = 5)
			
			
			#### RDA #####
			mod1 <- vegan::rda(distme ~ treatment + Condition(block), data = env, )
			
		
			# Look at eclidean distance after CLR transformation
			# since Bray-curtis cant be calculated from CLR
			
			# Change counts to relative abundances, center log-ration, and ANCOMBC 
			elk_rdp_ASV_adujsted.ps
			elk_rdp_VST.ps
			
			sample_data(elk_rdp_adujsted.ps)$Elk <- as.factor(sample_data(elk_rdp_adujsted.ps)$Elk)
			sample_data(elk_rdp_adujsted.ps)$Day <- as.factor(sample_data(elk_rdp_adujsted.ps)$Day)
			sample_data(ps4.clr)$Day <- as.factor(sample_data(ps4.clr)$Day)
			
		
			### plot dissimilarity
			groupDisim_bray_Day <- phyloseq_group_dissimilarity(elk_rdp_adujsted.ps, notch = F,
																													group = "Day", method = "bray", 
																													method_title = F, between_groups = F) + 
				theme_bw() + theme(legend.title = element_blank(),legend.text = element_text(size= 12),
													 axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1))
			
			ggsave2(filename = "group_disim_bray_Day_ANCOM_adjusted.pdf" ,
							plot = groupDisim_bray_Day, 
							path = "r_output/visualizations/beta_div_plots/",width = 4, height = 5)
			
			groupDisim_bray_elk <- phyloseq_group_dissimilarity(elk_rdp_adujsted.ps, notch = F,
																													group = "Elk", method = "bray", 
																													method_title = F, between_groups = F) + 
				theme_bw() + theme(legend.title = element_blank(),legend.text = element_text(size= 12),
													 axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1))
			
			ggsave2(filename = "group_disim_bray_elk_ANCOM_adjusted.pdf" ,
							plot = groupDisim_bray_elk, 
							path = "r_output/visualizations/beta_div_plots/",width = 4, height = 5)
			
			
		
			# heatmap of ecological distance
			table(tax_table(elk_rdp_adujsted.ps)[, "Phylum"], exclude = NULL)
			
			arrange(as_tibble(sample_data(elk_rdp_adujsted_top300.ps)), Elk, Day, Rep)$Description
			elk_rdp_top300.ps <- prune_taxa(names(sort(taxa_sums(elk_rdp.ps),TRUE)[1:100]), elk_rdp.ps)
			elk_rdp_Bactero.ps <- subset_taxa(elk_rdp.ps, Phylum=="Bacteroidetes")
			elk_rdp_Bactero.ps <- prune_taxa(names(sort(taxa_sums(elk_rdp_Bactero.ps),TRUE)[1:100]), elk_rdp_Bactero.ps)
			plot_heatmap(elk_rdp_Bactero.ps, method = "NMDS", distance = "bray", 
									 sample.label = "Elk", na.value="white",
									 sample.order = arrange(as_tibble(sample_data(elk_rdp_Bactero.ps)), Elk, Day, Rep)$Description)
			
			
			
			pheatmap::pheatmap(mat = as.matrix(phyloseq::otu_table(elk_rdp_Bactero.ps)), 
												 cluster_rows = T, cluster_cols = T, 
												 show_rownames = T, show_colnames = F, 
												 labels_row  = as.character(arrange(as_tibble(sample_data(elk_rdp_Bactero.ps)), Elk, Day, Rep)$Description))
			
#### Adonis on Aitchisons ####
			devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
			remotes::install_github("gavinsimpson/ggvegan", ref = "master", force = T)
			library(pairwiseAdonis)
			library("ggvegan");library(vegan);library(ggplot2)
			# autoplot()
			# get rid of psuedo reps?
			# Elk.rep1_ancom.ps	<- 	subset_samples(elk_rdp_ASV_adujsted.ps, Rep == "1")
			# Elk.1_ancom.ps	<- 	subset_samples(elk_rdp_ASV_adujsted.ps, Elk == "Elk.1")
			
			#aithchisons dist
			distme <- get_dist(elk_rdp_ASV_adujsted.ps, distmethod ="euclidean", method=NULL)
			sampleda <- data.frame(sample_data(elk_rdp_ASV_adujsted.ps), check.names=FALSE)
			sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
			sampleda$Rep <- factor(sampleda$Rep)
			sampleda$Elk <- factor(sampleda$Elk)
			# sampleda$Day <- factor(sampleda$Day)
			sampleda$Day <- as.numeric(as.character(sampleda$Day))
			sampleda$DayElk <- as.factor(paste0("day.",as.numeric(as.character(sampleda$Day)),"_elk.",sampleda$Elk ))
			length(sampleda$DayElk)#60
			length(levels(sampleda$DayElk))#20
			# sampleda$Rep
		  	
			xtabs(~ Elk+Day, sampleda)
		  set.seed(1024)
		  # look at trend in Day, but you have replication of Elk within Day and psuedo replication of Elk
		  # ideal formula = Day + (Elk|Day) but cant use random effects in vegan
		  ##this is pseudo rep but it should lower power so its extra conservative

		  perm <- with(sampleda, how(blocks = Elk,
		  													 plots = Plots(strata = Day, type = "series", mirror = F), #Days in Elk
		  													 within = Within(type = "free", mirror = F), #Reps in Day
		  													 nperm = 9999, minperm = 999
		  													 ))
			adores <- adonis2(distme ~ Elk * Day, data=sampleda, permutations = perm, 
												parallel = 10,
											  by = "terms" ) #terms sequential similar to random effect for the first term
			adores
			
			disper <- betadisper(distme, sampleda$Day)
			# when using sets with different beta-diversity you may
			# get false significant adonis results even when- the location/composition
			# is actually the same (!) this result is due to different
			# multivariate spread in dat1 and dat2:
			# Permdisp gives you a measure of the  dispersion of the samples 
			# (i.e. the variability in the community composition)within your groups. 
			# So you might interpret a significant permdisp result as one group of 
			# communities exhibits more or less group dispersion than the rest
			# effectively you are looking for a description of the multi-dimensional 
			# shape of the data cloud in response to your treatments. It is not unreasonable 
			# to state that there are both location (PERMANOVA) and dispersion (betadisper 
			# or more correctly PERMDISP) effects.
			# anova(disper)
			permutest(disper, permutations = perm) #day ***sig
			permutest(disper, permutations = perm, pairwise = T) # *sig
			# Pairwise comparisons:
			# 	(Observed p-value below diagonal, permuted p-value above diagonal)
			# 0           1           3           7     14
			# 0              0.164600000 0.529000000 0.007300000 0.5329
			# 1  0.042325030             0.713500000 0.005000000 0.7175
			# 3  0.382279742 0.589449073             0.058700000 0.9930
			# 7  0.000082016 0.000042733 0.012824029             0.0659
			# 14 0.369566294 0.603706263 0.985881423 0.013025408    
			
			disper <- betadisper(distme, sampleda$Day)
			permutest(disper, permutations = perm, pairwise = T) #DayElk *sig
			
			
			plot(disper)
			# autoplot.cca(disper)

			boxplot(disper)
			disper.HSD <- TukeyHSD(disper)
			disper.HSD
			plot(disper.HSD)
			
			pairwiseAdonis::pairwise.adonis2(distme ~ Day, data=sampleda, strata = as.factor(sampleda$Elk) )
			

#### plots for abundance Elk ####
     # library(ggplot2)
     library(MicrobiotaProcess)
     
     genis_taxa <- get_taxadf(obj=elk_rdp_ASV_adujsted.ps, taxlevel=6)
     # classtaxa@sam_data$Elk
     # The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`). 
     bars_genis<- ggbartax(obj=genis_taxa, count=F, facetNames="Elk", topn=25 ) +
     	xlab(NULL) +
     	ylab("relative abundance (%)") +
     	# scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
     	guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
     bars_genis
     ggsave2("r_output/visualizations/barPlot_Genus_ancom.pdf", 
     				plot = bars_genis, width = 6, height = 6, units = "in")
     
 
#### Global test HMp ####
      HMP_test.ps <- tax_glom(elk_rdp.ps, taxrank = "Genus", NArm = T)
      HMP_test.ps <- subset_taxa(HMP_test.ps, Genus!="g__Unknown" )
      HMP_test.ps <- subset_samples(HMP_test.ps, Elk == 4)
      # sample_data(HMP_test.ps)
      Day1   <- phyloseq::subset_samples(HMP_test.ps, Day == 0)
      Day2   <- phyloseq::subset_samples(HMP_test.ps, Day == 1)
      Day3   <- phyloseq::subset_samples(HMP_test.ps, Day == 3)
      Day4   <- phyloseq::subset_samples(HMP_test.ps, Day == 7)
      Day5   <- phyloseq::subset_samples(HMP_test.ps, Day == 14)
      #Output OTU tables
      Day1_otu <- data.frame(phyloseq::otu_table(Day1))
      Day2_otu <- data.frame(phyloseq::otu_table(Day2))
      Day3_otu <- data.frame(phyloseq::otu_table(Day3))
      Day4_otu <- data.frame(phyloseq::otu_table(Day4))
      Day5_otu <- data.frame(phyloseq::otu_table(Day5))
      #Group rare phyla
      #HMP test
      group_data <- list(Day1_otu, Day2_otu, Day3_otu, Day4_otu, Day5_otu)
      # The xdc test follows a Chi-square distribution with degrees of freedom equal to 
      # (J-1)*K, where J is the number of groups and K is the number of taxa.
      xdc <- HMP::Xdc.sevsample(group_data)   
      xdc
      # $`Xdc statistics`
      # [1] 400.8842
      # 
      # $`p value`
      # [1] 2.331468e-14
#### SelBal balances #####
      devtools::install_github(repo = "malucalle/selbal")
			source("code/selbal_functions.r")
      
      df1.CLRday0_long
      
      
#### Core Venn diagrams ####
			install.packages("ggVennDiagram")
      library(ggVennDiagram)
      #goal: create a venn for each elk across sample days
      ps1.comp <- microbiome::transform(elk_rdp.ps, "compositional")
      Elk1.ps.comp<- subset_samples(ps1.comp,	Elk=="1")
      Elk2.ps.comp<- subset_samples(ps1.comp,	Elk=="2")
      Elk3.ps.comp<- subset_samples(ps1.comp,	Elk=="3")
      Elk4.ps.comp<- subset_samples(ps1.comp,	Elk=="4")
      
      Day_variable <- unique(as.character(meta(Elk4.ps1.comp)$Day))
			
   ##### elk 1 ####
      #asvs
      list_core_Elk1 <- c() # an empty object to store information
      for (n in Day_variable){ # for each variable n in day
      	#print(paste0("Identifying Core Taxa for ", n))
      	ps.sub <- subset_samples(Elk1.ps.comp, Day_variable == n) # Choose sample from Day by 
      	core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
      												 detection = 0.001, # 1/1000 = 0.001 proportion for presence
      												 prevalence = 0.20) #in at least 20% of samples (1/5 samples)
      	print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each Day sample
      	list_core_Elk1[[paste0("Day ",as.character(n))]] <- core_m # add to a list core taxa for each group.
      	list_core_Elk1 <- list_core_Elk1[c("Day 0","Day 1","Day 3","Day 7","Day 14")]
      }
      #genus
      Elk1.ps.comp.gen <- tax_glom(Elk1.ps.comp, taxrank = "Genus")
      list_core_Elk1_gen <- c() # an empty object to store information
      for (n in Day_variable){ # for each variable n in day
      	#print(paste0("Identifying Core Taxa for ", n))
      	ps.sub <- subset_samples(Elk1.ps.comp.gen, Day_variable == n) # Choose sample from Day by n
      	core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
      												 detection = 0.00001, # 1/10000 = 0.00001 proportion for presence
      												 prevalence = 0.20, include.lowest = T) #in at least 20% of samples (1/5 total)
      	print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each Day sample
      	list_core_Elk1_gen[[paste0("Day ",as.character(n))]] <- core_m # add to a list core taxa for each group.
      	list_core_Elk1_gen <- list_core_Elk1_gen[c("Day 0","Day 1","Day 3","Day 7","Day 14")]
      }
      #plots
      elk1_venn_asv <-  ggVennDiagram::ggVennDiagram(list_core_Elk1, edge_lty = "solid", 
      																							 label = "count",label_size = 4, 
      																							 label_alpha = 0, label_color = "Black", edge_size = 0.9) +
																					      	scale_fill_distiller(palette = "Spectral") + 
																					      	labs(title = "Elk 1 (ASVs)",
																					      			 subtitle = "threshold = 0.1%") +
																					      	scale_x_continuous(expand = expansion(mult = .2)) + theme_void()
      elk1_venn_genus <-  ggVennDiagram::ggVennDiagram(list_core_Elk1_gen, edge_lty = "solid", 
      																							 label = "count",label_size = 4, 
      																							 label_alpha = 0, label_color = "Black", edge_size = 0.9) +
																					      	scale_fill_distiller(palette = "Spectral") + 
																					      	labs(title = "Elk 1 (Genus)",
																					      			 subtitle = "threshold = 0.001%") +
																					      	scale_x_continuous(expand = expansion(mult = .2)) + theme_void()
   ##### elk2 #####
      #asvs
      list_core_Elk2 <- c() # an empty object to store information
      for (n in Day_variable){ # for each variable n in day
      	#print(paste0("Identifying Core Taxa for ", n))
      	ps.sub <- subset_samples(Elk2.ps.comp, Day_variable == n) # Choose sample from Day by 
      	core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
      												 detection = 0.001, # 1/1000 = 0.001 proportion for presence
      												 prevalence = 0.20) #in at least 20% of samples (1/5 samples)
      	print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each Day sample
      	list_core_Elk2[[paste0("Day ",as.character(n))]] <- core_m # add to a list core taxa for each group.
      	list_core_Elk2 <- list_core_Elk2[c("Day 0","Day 1","Day 3","Day 7","Day 14")]
      }
      #genus
      Elk2.ps.comp.gen <- tax_glom(Elk2.ps.comp, taxrank = "Genus")
      list_core_Elk2_gen <- c() # an empty object to store information
      for (n in Day_variable){ # for each variable n in day
      	#print(paste0("Identifying Core Taxa for ", n))
      	ps.sub <- subset_samples(Elk2.ps.comp.gen, Day_variable == n) # Choose sample from Day by n
      	core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
      												 detection = 0.00001, # 1/10000 = 0.00001 proportion for presence
      												 prevalence = 0.20, include.lowest = T) #in at least 20% of samples (1/5 total)
      	print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each Day sample
      	list_core_Elk2_gen[[paste0("Day ",as.character(n))]] <- core_m # add to a list core taxa for each group.
      	list_core_Elk2_gen <- list_core_Elk2_gen[c("Day 0","Day 1","Day 3","Day 7","Day 14")]
      }
      #plots
      elk2_venn_asv <-  ggVennDiagram::ggVennDiagram(list_core_Elk2, edge_lty = "solid", 
      																							 label = "count",label_size = 4, 
      																							 label_alpha = 0, label_color = "Black", edge_size = 0.9) +
      	scale_fill_distiller(palette = "Spectral") + 
      	labs(title = "Elk 2 (ASVs)",
      			 subtitle = "threshold = 0.1%") +
      	scale_x_continuous(expand = expansion(mult = .2)) + theme_void()
      elk2_venn_genus <-  ggVennDiagram::ggVennDiagram(list_core_Elk2_gen, edge_lty = "solid", 
      																								 label = "count",label_size = 4, 
      																								 label_alpha = 0, label_color = "Black", edge_size = 0.9) +
      	scale_fill_distiller(palette = "Spectral") + 
      	labs(title = "Elk 2 (Genus)",
      			 subtitle = "threshold = 0.001%") +
      	scale_x_continuous(expand = expansion(mult = .2)) + theme_void()
   ##### elk3 ####
      #asvs
      list_core_Elk3 <- c() # an empty object to store information
      for (n in Day_variable){ # for each variable n in day
      	#print(paste0("Identifying Core Taxa for ", n))
      	ps.sub <- subset_samples(Elk3.ps.comp, Day_variable == n) # Choose sample from Day by 
      	core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
      												 detection = 0.001, # 1/1000 = 0.001 proportion for presence
      												 prevalence = 0.20) #in at least 20% of samples (1/5 samples)
      	print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each Day sample
      	list_core_Elk3[[paste0("Day ",as.character(n))]] <- core_m # add to a list core taxa for each group.
      	list_core_Elk3 <- list_core_Elk3[c("Day 0","Day 1","Day 3","Day 7","Day 14")]
      }
      #genus
      Elk3.ps.comp.gen <- tax_glom(Elk3.ps.comp, taxrank = "Genus")
      list_core_Elk3_gen <- c() # an empty object to store information
      for (n in Day_variable){ # for each variable n in day
      	#print(paste0("Identifying Core Taxa for ", n))
      	ps.sub <- subset_samples(Elk3.ps.comp.gen, Day_variable == n) # Choose sample from Day by n
      	core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
      												 detection = 0.00001, # 1/10000 = 0.00001 proportion for presence
      												 prevalence = 0.20, include.lowest = T) #in at least 20% of samples (1/5 total)
      	print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each Day sample
      	list_core_Elk3_gen[[paste0("Day ",as.character(n))]] <- core_m # add to a list core taxa for each group.
      	list_core_Elk3_gen <- list_core_Elk3_gen[c("Day 0","Day 1","Day 3","Day 7","Day 14")]
      }
      #plots
      elk3_venn_asv <-  ggVennDiagram::ggVennDiagram(list_core_Elk3, edge_lty = "solid", 
      																							 label = "count",label_size = 4, 
      																							 label_alpha = 0, label_color = "Black", edge_size = 0.9) +
      	scale_fill_distiller(palette = "Spectral") + 
      	labs(title = "Elk 3 (ASVs)",
      			 subtitle = "threshold = 0.1%") +
      	scale_x_continuous(expand = expansion(mult = .2)) + theme_void()
      elk3_venn_genus <-  ggVennDiagram::ggVennDiagram(list_core_Elk3_gen, edge_lty = "solid", 
      																								 label = "count",label_size = 4, 
      																								 label_alpha = 0, label_color = "Black", edge_size = 0.9) +
      	scale_fill_distiller(palette = "Spectral") + 
      	labs(title = "Elk 3 (Genus)",
      			 subtitle = "threshold = 0.001%") +
      	scale_x_continuous(expand = expansion(mult = .2)) + theme_void()
   ##### elk4 ####
      #asvs
      list_core_Elk4 <- c() # an empty object to store information
      for (n in Day_variable){ # for each variable n in day
      	#print(paste0("Identifying Core Taxa for ", n))
      	ps.sub <- subset_samples(Elk4.ps.comp, Day_variable == n) # Choose sample from Day by 
      	core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
      												 detection = 0.001, # 1/1000 = 0.001 proportion for presence
      												 prevalence = 0.20) #in at least 20% of samples (1/5 samples)
      	print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each Day sample
      	list_core_Elk4[[paste0("Day ",as.character(n))]] <- core_m # add to a list core taxa for each group.
      	list_core_Elk4 <- list_core_Elk4[c("Day 0","Day 1","Day 3","Day 7","Day 14")]
      }
      #genus
      Elk4.ps.comp.gen <- tax_glom(Elk4.ps.comp, taxrank = "Genus")
      list_core_Elk4_gen <- c() # an empty object to store information
      for (n in Day_variable){ # for each variable n in day
      	#print(paste0("Identifying Core Taxa for ", n))
      	ps.sub <- subset_samples(Elk4.ps.comp.gen, Day_variable == n) # Choose sample from Day by n
      	core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
      												 detection = 0.00001, # 1/10000 = 0.00001 proportion for presence
      												 prevalence = 0.20, include.lowest = T) #in at least 20% of samples (1/5 total)
      	print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each Day sample
      	list_core_Elk4_gen[[paste0("Day ",as.character(n))]] <- core_m # add to a list core taxa for each group.
      	list_core_Elk4_gen <- list_core_Elk4_gen[c("Day 0","Day 1","Day 3","Day 7","Day 14")]
      }
      #plots
      elk4_venn_asv <-  ggVennDiagram::ggVennDiagram(list_core_Elk4, edge_lty = "solid", 
      																							 label = "count",label_size = 4, 
      																							 label_alpha = 0, label_color = "Black", edge_size = 0.9) +
      	scale_fill_distiller(palette = "Spectral") + 
      	labs(title = "Elk 4 (ASVs)",
      			 subtitle = "threshold = 0.1%") +
      	scale_x_continuous(expand = expansion(mult = .2)) + theme_void()
      elk4_venn_genus <-  ggVennDiagram::ggVennDiagram(list_core_Elk4_gen, edge_lty = "solid", 
      																								 label = "count",label_size = 4, 
      																								 label_alpha = 0, label_color = "Black", edge_size = 0.9) +
      	scale_fill_distiller(palette = "Spectral") + 
      	labs(title = "Elk 4 (Genus)",
      			 subtitle = "threshold = 0.001%") +
      	scale_x_continuous(expand = expansion(mult = .2)) + theme_void()
      
   ##### all venns combined ####
      # plot all
      combined_venns_asv <- ggpubr::ggarrange(elk1_venn_asv, elk2_venn_asv, elk3_venn_asv, elk4_venn_asv,
      																				align = "hv", ncol = 1,
      																				common.legend = T, legend = "bottom")
      
      combined_venns_genus <- ggpubr::ggarrange(elk1_venn_genus, elk2_venn_genus, elk3_venn_genus, elk4_venn_genus,
      																					align = "hv", ncol = 1,
      																					common.legend = T, legend = "bottom")
      
      combined_venns <- ggpubr::ggarrange(combined_venns_asv,combined_venns_genus,
      																		align = "hv", ncol = 2, legend = "bottom")
      
      ggsave2(filename = "r_output/visualizations/venns/all_venns.pdf", 
      				plot = combined_venns, 
      				width = 6, height = 12, units = "in")
			
      
### Horizon plots ####
      # https://blekhmanlab.github.io/biomehorizon/
      devtools::install_github("blekhmanlab/biomehorizon")
      library(biomehorizon)
      ps1.comp <- subset_taxa(elk_rdp.ps, Genus != "g__Unknown")
      ps1.comp <- microbiome::transform(ps1.comp, "compositional")
      AllElk.ps.comp.gen <- tax_glom(ps1.comp, taxrank = "Genus")
      Elk1.ps.comp.gen <- tax_glom(subset_samples(ps1.comp,	Elk=="1"), taxrank = "Genus")
      Elk2.ps.comp.gen <- tax_glom(subset_samples(ps1.comp,	Elk=="2"), taxrank = "Genus")
      Elk3.ps.comp.gen <- tax_glom(subset_samples(ps1.comp,	Elk=="3"), taxrank = "Genus")
      Elk4.ps.comp.gen <- tax_glom(subset_samples(ps1.comp,	Elk=="4"), taxrank = "Genus")
      
      phy_to_horiz <- function(physeq) {
											      	require(dplyr)
											      	horiz_list = list();
											      	horiz_list$otudata = as.data.frame(t(otu_table(physeq))) %>%
																										      	mutate(OTU_ID = row.names(t(otu_table(physeq)))) %>%
																										      	select(OTU_ID, everything())
															horiz_list$metadata =  meta(physeq) %>% mutate(
												      																subject = paste0("Elk",as.character(Elk)), 
												      																sample = Description, 
												      																collection_date = Day, 
												      																supplement = Rep) %>%
																					      						select(subject,sample, collection_date, supplement) %>% 
																					      						arrange(subject, collection_date, supplement);
											      	horiz_list$taxonomydata  = as.data.frame(tax_table(physeq)) %>%
																										      	mutate(taxon_id = row.names(as.data.frame(tax_table(physeq))),
																				      								taxonomy = paste("F: ",Family," G: ", Genus, sep = "")) %>%
																										      	select(taxon_id, taxonomy);
											      
														return(horiz_list)
      											}
      Elk1.horiz <- phy_to_horiz(Elk1.ps.comp.gen)
      Elk2.horiz <- phy_to_horiz(Elk2.ps.comp.gen)
      Elk3.horiz <- phy_to_horiz(Elk3.ps.comp.gen)
      Elk4.horiz <- phy_to_horiz(Elk4.ps.comp.gen)
      
   ##### elk 1 plots ####
							   paramList1 <- prepanel(otudata =  Elk1.horiz$otudata,
							   											 metadata = Elk1.horiz$metadata,
							   											 taxonomydata = Elk1.horiz$taxonomydata,
							   											 subj = "Elk1",
							   											 thresh_prevalence = 100/15*14,
							   											 thresh_abundance = .1,
							   											 facetLabelsByTaxonomy = T,
							   											 regularInterval = 2,
							   											 origin = function(y) { mean(y, na.rm = TRUE) },
							   											 band.thickness = 1, #1percent
							   											 nbands = 5)
							  
							   horiz_elk1_band1_gg  <-  horizonplot(paramList1, aesthetics = horizonaes(title = "Microbiome Horizon Plot", 
							   																	 xlabel = "Samples from Elk1 (3 replicates/day)", 
							   																	 ylabel = "Taxa found in 14/15 samples at <0.1%", 
							   																	 legendTitle = "1% abund.\nQuartiles", 
							   																	 legendPosition	= "right")) +
																								   scale_x_continuous(
																								   					expand = c(0,0),
																								   					# breaks = seq(from = 1, by=1, to=16),
																								   					breaks = c(1, 2, 4, 8, 15),
																								   					labels = c("Day0", "Day1", "Day3", 
																								   										 "Day7", "Day14"))+
																								   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
							   horiz_elk1_band1_gg
							      
							   ggsave2(filename = "r_output/visualizations/horizon_AbundTime_plots/horiz_elk1_prev93_abund0.1_originMean1percBand.pdf", 
							   				plot = horiz_elk1_band1_gg, 
							   				width = 8, height = 12, units = "in")

   ##### elk 2 plots ####
   paramList2 <- prepanel(otudata =  Elk2.horiz$otudata,
   											 metadata = Elk2.horiz$metadata,
   											 taxonomydata = Elk2.horiz$taxonomydata,
   											 subj = "Elk2",
   											 thresh_prevalence = 100/15*14,
   											 thresh_abundance = .1,
   											 facetLabelsByTaxonomy = T,
   											 regularInterval = 2,
   											 # maxGap = 5, minSamplesPerFacet = 2,
   											 # maxGap = 4,
   											 # origin = 0,
   											 origin = function(y) { mean(y, na.rm = TRUE) },
   											 # band.thickness = function(y) {max((abs(y - origin(y))), na.rm=TRUE) / nbands},
   											 band.thickness = 1, #1percent
   											 nbands = 5)
   
   horiz_elk2_band1_gg  <-  horizonplot(paramList2, aesthetics = horizonaes(title = "Microbiome Horizon Plot", 
   																																					xlabel = "Samples from Elk2 (3 replicates/day)", 
   																																					ylabel = "Taxa found in 14/15 samples at <0.1%", 
   																																					legendTitle = "1% abund.\nQuartiles", 
   																																					legendPosition	= "right")) +
																													   	scale_x_continuous(
																													   		expand = c(0,0),
																													   		# breaks = seq(from = 1, by=1, to=16),
																													   		breaks = c(1, 2, 4, 8, 15),
																													   		labels = c("Day0", "Day1", "Day3", 
																													   							 "Day7", "Day14"))+
																													   	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
   horiz_elk2_band1_gg
   
   ggsave2(filename = "r_output/visualizations/horizon_AbundTime_plots/horiz_elk2_prev93_abund0.1_originMean1percBand.pdf", 
   				plot = horiz_elk2_band1_gg, 
   				width = 8, height = 12, units = "in")
   
   ##### elk 3 plots ####
   paramList3 <- prepanel(otudata =  Elk3.horiz$otudata,
   											 metadata = Elk3.horiz$metadata,
   											 taxonomydata = Elk3.horiz$taxonomydata,
   											 subj = "Elk3",
   											 thresh_prevalence = 100/15*14,
   											 thresh_abundance = .1,
   											 facetLabelsByTaxonomy = T,
   											 regularInterval = 2,
   											 # maxGap = 5, minSamplesPerFacet = 2,
   											 # maxGap = 4,
   											 # origin = 0,
   											 origin = function(y) { mean(y, na.rm = TRUE) },
   											 # band.thickness = function(y) {max((abs(y - origin(y))), na.rm=TRUE) / nbands},
   											 band.thickness = 1, #1percent
   											 nbands = 5)
   
   horiz_elk3_band1_gg  <-  horizonplot(paramList3, aesthetics = horizonaes(title = "Microbiome Horizon Plot", 
   																																					xlabel = "Samples from Elk3 (3 replicates/day)", 
   																																					ylabel = "Taxa found in 14/15 samples at <0.1%", 
   																																					legendTitle = "1% abund.\nQuartiles", 
   																																					legendPosition	= "right")) +
																																   	scale_x_continuous(
																																   		expand = c(0,0),
																																   		# breaks = seq(from = 1, by=1, to=16),
																																   		breaks = c(1, 2, 4, 8, 15),
																																   		labels = c("Day0", "Day1", "Day3", 
																																   							 "Day7", "Day14"))+
																																   	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
   horiz_elk3_band1_gg
   
   ggsave2(filename = "r_output/visualizations/horizon_AbundTime_plots/horiz_elk3_prev93_abund0.1_originMean1percBand.pdf", 
   				plot = horiz_elk1_band10_gg, 
   				width = 8, height = 12, units = "in")
   
   ##### elk 4 plots ####
   paramList4 <- prepanel(otudata =  Elk4.horiz$otudata,
   											 metadata = Elk4.horiz$metadata,
   											 taxonomydata = Elk4.horiz$taxonomydata,
   											 subj = "Elk4",
   											 thresh_prevalence = 100/15*14,
   											 thresh_abundance = .1,
   											 facetLabelsByTaxonomy = T,
   											 regularInterval = 2,
   											 # maxGap = 5, minSamplesPerFacet = 2,
   											 # maxGap = 4,
   											 # origin = 0,
   											 origin = function(y) { mean(y, na.rm = TRUE) },
   											 # band.thickness = function(y) {max((abs(y - origin(y))), na.rm=TRUE) / nbands},
   											 band.thickness = 1, #1percent
   											 nbands = 5)
   
   horiz_elk4_band1_gg  <-  horizonplot(paramList4, aesthetics = horizonaes(title = "Microbiome Horizon Plot", 
   																																					xlabel = "Samples from Elk (3 replicates/day)", 
   																																					ylabel = "Taxa found in 14/15 samples at <0.1%", 
   																																					legendTitle = "1% abund.\nQuartiles", 
   																																					legendPosition	= "right")) +
																																   	scale_x_continuous(
																																   		expand = c(0,0),
																																   		# breaks = seq(from = 1, by=1, to=16),
																																   		breaks = c(1, 2, 4, 8, 15),
																																   		labels = c("Day0", "Day1", "Day3", 
																																   							 "Day7", "Day14"))+
																																   	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
   horiz_elk4_band1_gg
   
   ggsave2(filename = "r_output/visualizations/horizon_AbundTime_plots/horiz_elk4_prev93_abund0.1_originMean1percBand.pdf", 
   				plot = horiz_elk4_band1_gg, 
   				width = 8, height = 12, units = "in")
   
   
   ##### combine plots ####
   horiz_elk4_band1_gg
   
   combined_horiz_1 <- ggpubr::ggarrange(horiz_elk1_band1_gg,horiz_elk2_band1_gg,
   																		align = "hv", ncol = 2, common.legend = T,
   																		legend = "right")
  
   combined_horiz_2 <- ggpubr::ggarrange(horiz_elk3_band1_gg,horiz_elk4_band1_gg,
   																			align = "hv", ncol = 2, common.legend = T,
   																			legend = "right")
   ggsave2(filename = "r_output/visualizations/horizon_AbundTime_plots/combined_horiz1.pdf", 
   				plot = combined_horiz_1, 
   				width = 20, height = 12, units = "in")
   ggsave2(filename = "r_output/visualizations/horizon_AbundTime_plots/combined_horiz2.pdf", 
   				plot = combined_horiz_2, 
   				width = 20, height = 12, units = "in")
   
   
   ##### plot 1 Taxa across elk #####
   bacteroides <- colnames(otu_table(subset_taxa(Elk3.ps.comp.gen, Genus=="Bacteroides")))[1]
   Sporobacter <- colnames(otu_table(subset_taxa(Elk3.ps.comp.gen, Genus=="Sporobacter")))[1]
   Intestinimonas <- colnames(otu_table(subset_taxa(Elk3.ps.comp.gen, Genus=="Intestinimonas")))[1]
   
   AllElk.ps.comp.gen
   AllElk.horiz <- phy_to_horiz(AllElk.ps.comp.gen)
   
   ###### bact plot ####

   paramList_bact <- prepanel(otudata = AllElk.horiz$otudata, 
   													 metadata = AllElk.horiz$metadata, 
   													 taxonomydata = AllElk.horiz$taxonomydata,
   													 singleVarOTU = bacteroides, 
   													 subj = c("Elk1","Elk2","Elk3","Elk4"),
													   thresh_prevalence = (100/15*14),
													   thresh_abundance = .1,
													   # regularInterval = 3,
													   # maxGap = 5, minSamplesPerFacet = 2,
													   # maxGap = 4,
													   # origin = 20,
													   origin = function(y) { mean(y, na.rm = TRUE) },
													   # band.thickness = function(y) {max((abs(y - origin(y))), na.rm=TRUE) / nbands},
													   band.thickness = 2, #1percent
													   nbands = 5)

horiz_all_elk_bact_gg  <-  horizonplot(paramList_bact, 
																		 aesthetics = horizonaes(title = "Genus:Bacteroides Horizon Plot", 
																					 xlabel = "Day of sampling", 
																					 ylabel = "Elk Sample (triplicates)", 
																					 legendTitle = "2% abund.\nQuartiles", 
																					 legendPosition	= "right")) +
																			scale_x_continuous(expand = c(0,0),
																				breaks = c(1, 2, 4, 8, 15),
																				labels = c("Day0", "Day1", "Day3", "Day7", "Day14"))+
																			theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
horiz_all_elk_bact_gg

ggsave2(filename = "r_output/visualizations/horizon_AbundTime_plots/horiz_all_elk_bacteroides.pdf", 
				plot = horiz_all_elk_bact_gg, 
				width = 6, height = 5, units = "in")

   ###### Sporo plot ####
paramList_Sporo <- prepanel(otudata = AllElk.horiz$otudata, 
													 metadata = AllElk.horiz$metadata, 
													 taxonomydata = AllElk.horiz$taxonomydata,
													 singleVarOTU = Sporobacter, 
													 subj = c("Elk1","Elk2","Elk3","Elk4"),
													 thresh_prevalence = (100/15*14),
													 thresh_abundance = .1,
													 # regularInterval = 3,
													 # maxGap = 5, minSamplesPerFacet = 2,
													 # maxGap = 4,
													 # origin = 20,
													 origin = function(y) { mean(y, na.rm = TRUE) },
													 # band.thickness = function(y) {max((abs(y - origin(y))), na.rm=TRUE) / nbands},
													 band.thickness = 2, #1percent
													 nbands = 5)

horiz_all_elk_Sporo_gg  <-  horizonplot(paramList_Sporo, 
																			 aesthetics = horizonaes(title = "Genus:Sporobacter Horizon Plot", 
																			 												xlabel = "Day of sampling", 
																			 												ylabel = "Elk Sample (triplicates)", 
																			 												legendTitle = "2% abund.\nQuartiles", 
																			 												legendPosition	= "right")) +
																						scale_x_continuous(expand = c(0,0),
																															 breaks = c(1, 2, 4, 8, 15),
																															 labels = c("Day0", "Day1", "Day3", "Day7", "Day14"))+
																						theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
horiz_all_elk_Sporo_gg

ggsave2(filename = "r_output/visualizations/horizon_AbundTime_plots/horiz_all_elk_Sporobacter.pdf", 
				plot = horiz_all_elk_Sporo_gg, 
				width = 6, height = 5, units = "in")

   ###### Intestinimonas ####
    paramList_Intest <- prepanel(otudata = AllElk.horiz$otudata, 
														metadata = AllElk.horiz$metadata, 
														taxonomydata = AllElk.horiz$taxonomydata,
														singleVarOTU = Intestinimonas, 
														subj = c("Elk1","Elk2","Elk3","Elk4"),
														thresh_prevalence = (100/15*14),
														thresh_abundance = .1,
														# regularInterval = 3,
														# maxGap = 5, minSamplesPerFacet = 2,
														# maxGap = 4,
														# origin = 20,
														origin = function(y) { mean(y, na.rm = TRUE) },
														# band.thickness = function(y) {max((abs(y - origin(y))), na.rm=TRUE) / nbands},
														band.thickness = 2, #1percent
														nbands = 5)

horiz_all_elk_Intest_gg  <-  horizonplot(paramList_Intest, 
																				aesthetics = horizonaes(title = "Genus:Intestinimonas Horizon Plot", 
																																xlabel = "Day of sampling", 
																																ylabel = "Elk Sample (triplicates)", 
																																legendTitle = "2% abund.\nQuartiles", 
																																legendPosition	= "right")) +
																							scale_x_continuous(expand = c(0,0),
																																 breaks = c(1, 2, 4, 8, 15),
																																 labels = c("Day0", "Day1", "Day3", "Day7", "Day14"))+
																							theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
horiz_all_elk_Intest_gg

ggsave2(filename = "r_output/visualizations/horizon_AbundTime_plots/horiz_all_elk_Intestinimonas.pdf", 
				plot = horiz_all_elk_Intest_gg, 
				width = 6, height = 5, units = "in")
