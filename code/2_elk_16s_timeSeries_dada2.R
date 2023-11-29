# Elk sampling time series
# Starting with Raw data (fastq's)
# read data in, trim barcodes
# FilterAndTrim with dada2, infer seq variants, seqTable
# Merge multiple runs, Remove Chimeras, Assign Taxonomy

##### Load libraries ####
  # if (!requireNamespace("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
  devtools::install_github("benjjneb/dada2") #new version 1.14.1
  BiocManager::install("ShortRead", update = F, force = T)
  BiocManager::install("Matrix", update = F, force = T)
  # BiocManager::install(c("phyloseq","metagenomeSeq","vegan"))
  # BiocManager::install(c("ggplot2","randomforests","ape","microbiome","dplyr","pheatmap","RColorBrewer","ggtree"))
  # BiocManager::install("genefilter")
  # BiocManager::install("gridExtra")
  # BiocManager::install("png")
  # BiocManager::install("DESeq2")
  # library(dada2); packageVersion("dada2")
  # library(biomformat)
  # library(mvabund);   packageVersion("mvabund")
  # library(Biostrings); packageVersion("Biostrings")
  # library(ShortRead);  packageVersion("ShortRead")

pacman::p_load("dada2", "mvabund", "Biostrings", "ShortRead", "Cairo", update = F)
packageVersion("dada2") #‘1.21.0’
packageVersion("ShortRead")#‘1.44.3’
packageVersion("Matrix") #‘1.3.4’
##### Load data ####

# Load data locations, each sequencing run should be seperate
path1 <- "./data_links" 
# CHANGE ME to the directories containing the fastq files
list.files(path1)
# root.dir = rprojroot::find_rstudio_root_file()
# setwd(root.dir)
getwd()

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R3.fastq
# save locations
fnFs.1 <- sort(list.files(path1, pattern= ".1.fq.gz", full.names = TRUE))
fnRs.1 <- sort(list.files(path1, pattern= ".2.fq.gz", full.names = TRUE))
length(fnFs.1) == length(fnRs.1)

##### REMOVE PRIMERS ######
#remove primers with cutadapt

FWD <- "CAGCMGCCGCGGTAATWC"  ## CHANGE ME to your forward primer sequence
REV <- "CCGTCAATTCMTTTRAGTTT"  ## CHANGE ME...
# FYI
# 16S_forward     CAGCMGCCGCGGTAATWC
# 16S_reverse     CCGTCAATTCMTTTRAGTTT
# rbcL_F  ATGTCACCACAAACAGAGACTAAAGC
# rbcL_R  GTAAAATCAAGTCCACCRCG
# ITS_F   ATGCGATACTTGGTGTGAAT
# ITS_R   TCCTCCGCTTATTGATATGC
# 18S_F   CAGCAGCCGCGGTAATTCC
# 18S_R   CCCGTGTTGAGTCAAATTAAGC

#function to get all possible variants of the primers
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
# FWD.orients

#make some new directory strings for filtering N's
fnFs.filtN.1 <- file.path(path1, "filtN.1", basename(fnFs.1)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN.1 <- file.path(path1, "filtN.1", basename(fnRs.1))


# Filter all the N's, save to a subdirectory named above w/in the raw data folder
OutFnt.1raw <- dada2::filterAndTrim(fnFs.1, fnFs.filtN.1, fnRs.1, fnRs.filtN.1, maxN = 0, multithread = T) #run 1

# count the primers in a single pair of reads for each read set
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN.1[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN.1[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN.1[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN.1[[1]]))
# Forward Complement Reverse RevComp
# FWD.ForwardReads      83          0       0       0
# FWD.ReverseReads       0          0       0       0
# REV.ForwardReads       0          0       0       2
# REV.ReverseReads     127          0       0       0

# Note: Orientation mixups are a common trip-up. If, for example, 
# the REV primer is matching the Reverse reads in its RevComp orientation, 
# then replace REV with its reverse-complement orientation 
# (REV <- REV.orient[["RevComp"]]) before proceeding.

# cutadapt <- Sys.which(names = "~/.local/bin/cutadapt")
cutadapt <- "/home/sam.pannoni/miniconda3/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

# cutadapt time!
path.cut.1 <- file.path(path1, "cutadapt")
if(!dir.exists(path.cut.1)) dir.create(path.cut.1)
fnFs.cut.1 <- file.path(path.cut.1, basename(fnFs.1))
fnRs.cut.1 <- file.path(path.cut.1, basename(fnRs.1))



FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# R1.flags <- paste("-g", FWD) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
# R2.flags <- paste("-G", REV) 
# Run Cutadapt
# take between 1s and 15s per read pair
#batch 1
for(i in seq_along(fnFs.1)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "--minimum-length", 100,
                             "--pair-filter=any --discard-untrimmed",
                             "-o", fnFs.cut.1[i], "-p", fnRs.cut.1[i], # output files
                             fnFs.filtN.1[i], fnRs.filtN.1[i])) # input files
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut.1[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut.1[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut.1[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut.1[[1]]))
# Forward Complement Reverse RevComp
# FWD.ForwardReads       0          0       0       0
# FWD.ReverseReads       0          0       0       0
# REV.ForwardReads       0          0       0       0
# REV.ReverseReads       0          0       0       0


##### Filter and Trim ##### 

# Forward and reverse fastq filenames have the format:
cutFs.1 <- sort(list.files(path.cut.1, pattern = ".1.fq.gz", full.names = TRUE))
cutRs.1 <- sort(list.files(path.cut.1, pattern = ".2.fq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: "ETS_012.1.fq.gz"
# strsplit(basename(cutFs.1), "\\.")[[1]][1] this cuts on "." and then selects the first instance of it
get.sample.name <- function(fname) strsplit(basename(fname), "\\.")[[1]][1]
sample.names.1 <- unname(sapply(cutFs.1, get.sample.name))

head(sample.names.1)


#plot forward read quality
plotQualityProfile(cutFs.1[2])

#plot reverse read quality
plotQualityProfile(cutRs.1[1:2])


lens.fn <- lapply(cutFs.1, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
lens.fn <- lapply(cutRs.1, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
# sample.names.1 <- sapply(strsplit(basename(fnFs.1), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs.1 <- file.path(path.cut.1, "filtered", basename(cutFs.1))
filtRs.1 <- file.path(path.cut.1, "filtered", basename(cutRs.1))

#same number of forward and Rev reads check
if(length(filtFs.1) != length(filtRs.1)) stop("Forward and reverse files do not match.")

outFnT.1 <- filterAndTrim(cutFs.1, filtFs.1, cutRs.1, filtRs.1, matchIDs = F,
                        truncLen=c(250,210), maxN=0, maxEE=c(2,2), minLen = 100,
                        truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
#
##### Learn error rates ####
set.seed(100);                                                                                               
Sys.getenv("DISPLAY") #:0 for default
Sys.setenv("DISPLAY"=":50") #for taranis x2go window
# CairoX11();

errR.1 <- learnErrors(filtRs.1, nbases = 1e8, multithread=TRUE, randomize=TRUE)
errF.1 <- learnErrors(filtFs.1, nbases = 1e8, multithread=TRUE, randomize=TRUE)
plotErrors(errR.1, nominalQ=TRUE)
plotErrors(errF.1, nominalQ=TRUE)


##### dereplicate fastq ####
derepFs.1 <- derepFastq(filtFs.1, verbose=TRUE)
derepRs.1 <- derepFastq(filtRs.1, verbose=TRUE)
derepFs.1[[1]]
derepRs.1[[1]]

#Name the derep-class objects by the sample names
names(derepFs.1) <- sample.names.1
names(derepRs.1) <- sample.names.1

#infer variants
dadaFs.1 <- dada(derepFs.1, err=errF.1, multithread=TRUE)
dadaRs.1 <- dada(derepRs.1, err=errR.1, multithread=TRUE)

##### MERGE PAIRED READS ########

#try without concat
mergers.1 <- mergePairs(dadaFs.1, derepFs.1, dadaRs.1, derepRs.1, returnRejects=F, verbose=TRUE)
head(mergers.1[[1]])


seqtab_1 <- makeSequenceTable(mergers.1)
dim(seqtab_1)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_1)))

seqtab1.noChim <- removeBimeraDenovo(seqtab_1, method="consensus", multithread=TRUE)
sum(seqtab1.noChim)/sum(seqtab_1)
# [1] 0.9050513 % retained

# This next command is analogous to “cutting a band” in-silico 
# to get amplicons of the targeted length.

#seqtab2 <- seqtab1[,nchar(colnames(seqtab1)) %in% seq(250,256)]

#count reads
getN <- function(x) sum(getUniques(x))
OutFnt.1raw # first filter
outFnT.1 # filter after cutadapt
dadaFs.1 # dada2 variants
dadaRs.1 #dada2 variants
mergers.1 #merged F+R
seqtab_1 # sequences
seqtab1.noChim #seqs w/o bimeras
#run 1
track <- cbind(OutFnt.1raw[,1], outFnT.1, sapply(dadaFs.1, getN), sapply(dadaRs.1, getN), sapply(mergers.1, getN), rowSums(seqtab1.noChim),rowSums(seqtab1.noChim)/OutFnt.1raw[,1])

colnames(track) <- c("total_raw_pairs_in", "Primer_trim_filter", "Quality_filter_trim", "denoisedF", "denoisedR", "merged", "non_chimera", "proportion")
rownames(track) <- sample.names.1
head(track)
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)

#sum into one table
options(scipen=1)
read_summary_run1 <- summary(track)
col_totals_1 <- colSums2(track)
Sum_with_totals_1 <- rbind(read_summary_run1, col_totals_1)

system2(command = "mkdir", args = "r_output/filtering_stats")
write.csv(Sum_with_totals_1, file = "r_output/filtering_stats/sequence_filtering_stats.csv")
write.csv(track,             file = "r_output/filtering_stats/sample_filtering_seq_totals.csv")


# summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1], 
#                           filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
#                           dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN), 
#                           nonchim=rowSums(seqtab.nochim), 
#                           final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

##### file making with "ASV" notation ####
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab1.noChim)
asv_headers <- vector(dim(seqtab1.noChim)[2], mode="character")

for (i in 1:dim(seqtab1.noChim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))

cat("\n\nSave unique fasta file....\n")
write(asv_fasta, "r_output/ASVs.fa")
# write(asv_fasta_nNs, "output_r/out_rbcl/ASVs_nNs.fa")
# write(asv_fasta_trunc, "output_r/out_rbcl/ASVs_truncated.fa")

# count table:
asv_tab <- t(seqtab1.noChim)
row.names(asv_tab) <- sub(">", "", asv_headers)
cat("\n\nSave ASV table...\n")
write.table(asv_tab, "r_output/ASVs_counts.txt", 
            sep="\t", quote=F, col.names=NA)



##### metadata from exisitng files *if file names are data rich ####
# sample.names.all <- c(sample.names.1,sample.names.2) #omitted sample.names.3
# ?str_split
# sample_data <- stringr::str_split(string = sample.names.all, "_",simplify = T)
# sample_data <- cbind(sample.names.all, sample_data)
# colnames(sample_data) <- c("file_name", "Grazer_Species", "Amplicon", "Sample_number")
# write.table(sample_data, "output_r/out_rbcl/metadata_rbcl.tsv", sep = "\t", col.names = T)


##### Assign taxonomy ####
tax_gg <- assignTaxonomy(seqtab1.noChim, "../taxonomy_databases/gg_13_8_train_set_97.fa.gz", 
                      multithread=TRUE, tryRC = TRUE)
tax_rdp16 <- assignTaxonomy(seqtab1.noChim, "../taxonomy_databases/rdp_train_set_16.fa.gz", 
                      multithread=TRUE, tryRC = TRUE)
# tax_rdp14 <- assignTaxonomy(seqtab.noChim, "../../R_taxonomy_dbs/for_dada2/rdp_train_set_14.fa.gz", 
#                           multithread=TRUE, tryRC = TRUE)
tax_refseq <- assignTaxonomy(seqtab1.noChim, "../taxonomy_databases/RefSeq-RDP16S_v3_May2018.fa.gz", 
                            multithread=TRUE, tryRC = TRUE)


# rm(list = c(""))
taxa_rdp16_wSpp <-addSpecies(tax_rdp16, "../taxonomy_databases/rdp_species_assignment_16.fa.gz",verbose = T )
taxa_refseq_wSpp <-addSpecies(tax_refseq, "../taxonomy_databases/RefSeq-RDP_dada2_assignment_species.fa.gz",verbose = T )


# CairoX11()

# save tax table
asv_tax_gg <- tax_gg
row.names(asv_tax_gg) <- sub(">", "", asv_headers)
cat("\n\nSave taxonomy file....\n")
write.table(asv_tax_gg, "r_output/gg_taxonomy.txt", sep="\t", quote=F, col.names=T, row.names = T)

asv_tax_refseq <- tax_refseq
row.names(asv_tax_refseq) <- sub(">", "", asv_headers)
cat("\n\nSave taxonomy file....\n")
write.table(asv_tax_refseq, "r_output/rdp_taxonomy.txt", sep="\t", quote=F, col.names=T, row.names = T)

asv_tax_rdp16 <- taxa_rdp16_wSpp
row.names(asv_tax_rdp16) <- sub(">", "", asv_headers)
cat("\n\nSave taxonomy file....\n")
write.table(asv_tax_rdp16, "r_output/rdp16_spp_taxonomy.txt", sep="\t", quote=F, col.names=T, row.names = T)


# asv_tax_open_plus <- tax_open_plus
# row.names(asv_tax_open_plus) <- sub(">", "", asv_headers)
# cat("\n\nSave taxonomy file....\n")
# write.table(asv_tax_open_plus, "output_r/out_rbcl/ASVs_rbcL_taxonomy_open_plus_d2.txt", sep="\t", quote=F, col.names=T, row.names = T)

# #resume some stuff in qiime if wanted
# uniquesToFasta(getUniques(seqtab.nochim), fout="output/18s_unique_seqs.fasta", ids=paste0("Seq", seq(length(getUniques(seqtab.nochim)))))
# ids_sam <- paste0("Seq", seq(length(getUniques(seqtab.nochim))))
# write.table(ids_matrix, "output/ids_18s_nochim.txt", sep="\t", col.names = F, row.names = F)


########### reload the files #####
library(phyloseq)
# seqtab = readRDS("data/seqtab_final.rds")
cat("\n\nRead in the files for phyloseq....\n")
tax_rdp_spp <-  read.table("r_output/rdp16_spp_taxonomy.txt", header=T, sep = "\t")
counts_table <- read.table("r_output/ASVs_counts.txt", header=TRUE)
metadata <- read.csv("./time_series_mappingFile_R.txt", header=TRUE,  sep="\t")

# metadata$SampleID <- gsub("RBCL", "rbcL",x = metadata$SampleID)
# metadata$SampleID <- gsub("\\.", "_",x = metadata$SampleID, perl = T)
rownames(metadata) <- metadata$Description
psmeta <- sample_data(metadata)
# colnames(psmeta)[7] <- "Grazer_Species"
ASV_fasta_st <- readDNAStringSet("r_output/ASVs.fa")
str(ASV_fasta_st)
names(ASV_fasta_st)
# ASV_fasta_st <- gsub("N", "", fasta_final_for_blast_q2$rbcL_sequence)
#transpose
tran.OTU <- as.matrix(t(counts_table))
rownames(tran.OTU)

#clean taxonomy
class(tax_rdp_spp)
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
elk_ts_rdp.ps <- phyloseq(OTU, meta, TAX.rdp)

#### Construct phylogenetic tree ####

# The gene for the large subunit of the ribulose-bisphosphate carboxylase (rbcL), 
# located on the chloroplast genome, is an appropriate choice for inference of phylogenetic 
# relationships at higher taxo- nomic levels (13-15). Because of its slow synonymous nucleotide 
# substitution rate in comparison with nuclear genes and its functional constraint that reduces 
# the evolutionary rate of nonsynonymous substitutions, rbcL is considered to be more useful than 
# the isozymes and the restriction fragment length polymorphisms at these taxonomic levels. 

# Phylogenetic relatedness is commonly used to inform downstream 
# analyses, especially the calculation of phylogeny-aware distances 
# between microbial communities. The DADA2 sequence inference method 
# is reference-free, so we must construct the phylogenetic tree 
# relating the inferred sequence variants de novo. We begin by 
# performing a multiple-alignment using the 
# DECIPHER R package (Wright 2015).

names(ASV_fasta_st) <- ASV_fasta_st@ranges@NAMES # This propagates to the tip labels of the tree
str(ASV_fasta_st)
head(ASV_fasta_st)

library("DECIPHER")
library(phangorn)
citation("phangorn")
numcores <- detectCores()
cat("\n\nMake alignment from fasta....\n")
alignment <- AlignSeqs(DNAStringSet(ASV_fasta_st), processors = numcores-1, anchor=NA, verbose=FALSE)

# The phangorn R package is then used to construct a phylogenetic tree. 
# https://f1000research.com/articles/5-1492/v2
# Here we first construct a neighbor-joining tree, and then fit a GTR+G+I
# (Generalized time-reversible with Gamma rate variation) maximum likelihood 
# tree using the neighbor-joining tree as a starting point.
cat("\n\nMake distance matrix with phangorn....\n")
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
# str(dm)
cat("\n\nMake NJ tree....\n")
treeNJ <- NJ(dm) # Note, tip order != sequence order
#maximum likelihood tree 
cat("\n\nFit ML model....\n")
fit = pml(treeNJ, data=phangAlign)
print(fit)
cat("\n\n optimize model....\n")
optim.pml(fit)
# plot(vert.tree,no.margin=TRUE,edge.width=2)
fitGTR <- update(fit, k=4, inv=0.2)
cat("\n\nFit GTR Tree....\n")
fitGTR.2 <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, multicore=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
#start 14:05 12/2/2021
#fin	 23:00
cat("\n\nMerge tree into phyloseq object....\n")

Tree<-fitGTR.2$tree
cat("\n\nSave tree and phyloseq object....\n")
write.tree(Tree, file="r_output/bootstrap_GTR_tree.tre")

cat("\n\nEnd of script \n")

#worked on taranis

######







