###############################################################
# Script to analyze 18S data with DADA2 - ITS workflow        #
# Based on DADA2 workflow for Big Data by Benjamin Callahan   #
# Data: Miseq-18S-demultiplexed xoxo                          #
# Mona Parizadeh - May 2025                                   #
###############################################################

# Load libraries ####
library(dada2); packageVersion("dada2") 
require(parallel)
library(phyloseq); packageVersion("phyloseq") 
library(ShortRead); packageVersion("ShortRead") 

# Define path ####
path <- "~/Documents/xoxo_article/Sequences18S_Microbiome_Insights"
list.files(path)

# Filter and trim ####
#Lists forward and reverse fastq files and sorts the files to make sure they are all in the same direction (R1 and R2)
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq.gz", full.names = TRUE)) # F 
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq.gz", full.names = TRUE)) # R 
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX
sample.names <- sapply(strsplit(basename(fnFs), "_L001"), `[`, 1)
head(sample.names)

# Identify primers ####
FWD <- "CYGCGGTAATTCCAGCTC"  
REV <- "AYGGTATCTRATCRTCTTYG"  
#Verify the presence and orientation of the primers 
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
#pre-filter the sequences to remove those w/ Ns
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#Identify and count primers
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients,
                                primerHits, fn = fnRs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#no primers found

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have the format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Examine quality profiles of forward and reverse reads
plotQualityProfile(cutFs[1:2]) # forward reads 
plotQualityProfile(cutRs[1:2]) # reverse reads 

# Perform filtering and trimming 
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

detectCores() #in order to choose multithread

# Filter the forward and reverse reads
#takes time
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(250,200), 
                     maxN=0, maxEE=c(2,2), minLen = 50, truncQ=2, rm.phix=TRUE, compress=TRUE, 
                     verbose = TRUE, matchIDs = TRUE , multithread=TRUE) 
head(out)

# Calculate error rates ####
errF <- learnErrors(filtFs, multithread=TRUE)  
errR <- learnErrors(filtRs, multithread=TRUE) 

# Dereplicate the filtered fastq files ####
#This combines all identical sequencing reads into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence, to reduce computation time by eliminating redundant comparisons
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer the sequence variants in each sample ####
#This produces an object that describes dada2 denoising results (removes all sequencing errors)
#Using the developed error model, it calculates abundance p-values for each unique sequence and checks if a given error rate is too abundant in a sample to be explained by sequencing errors
dadaFs <- dada(derepFs, err=errF, pool=TRUE, multithread=TRUE)   
dadaRs <- dada(derepRs, err=errR, pool=TRUE, multithread=TRUE)
dadaFs[[1]] 
dadaRs[[1]] 

# Merge the denoised forward and reverse reads ####
#This merges error-free forward and reverse reads (after being denoised) and removes paired reads that did not exactly overlapped
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 10, maxMismatch = 0, returnRejects = TRUE, verbose=TRUE) 

# Construct sequence table ####
seqtab <- makeSequenceTable(mergers) 
#The sequences being tabled vary in length; analogous to an OTU table
dim(seqtab) 
head(mergers[[1]])

# Inspect distribution of sequence lengths (length of ASVs):
table(nchar(getSequences(seqtab)))

#save
saveRDS(seqtab, file = "~/Documents/xoxo_article/files/18s/seqtab_xoxo18s.rds")

#### Remove chimeric sequences ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#takes time
# Chimera is a PCR artifact (human-made), while sequencing a 16S region gene, primer starts amplifying the region, but at some point it'll fall off and not go all the way through, making an oligo that can later be used to prime another completely unrelated 16S region. So we'll get the amplicons that are actually fusions of two or more parent sequences
# The algorithm will start from the least abundant reads and will test all the combinations with all the more abundant sequences to see if there's any overlap
# A bimera is a two-parent chimera, in which the left side is made up of one parent sequence, and the right-side made up of a second parent sequence.
# if half of the reads are being removed as chimeras, we probably still have primers on our samples.
#dim(seqtab.nochim) 
sum(seqtab.nochim)/sum(seqtab) 
#Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though). If most of your reads were removed as chimeric, upstream processing may need to be revisited. In almost all cases this is caused by primer sequences with ambiguous nucleotides that were not removed prior to beginning the DADA2 pipeline.

#collapse mismatches (just testing) ####
seqtab.col = collapseNoMismatch(seqtab.nochim, minOverlap = 280, verbose = TRUE)
dim(seqtab.col) 

saveRDS(seqtab.nochim, file = "~/Documents/xoxo_article/files/18s/seqtabNoChim_xoxo18s.rds")

# Track reads through the pipeline ####
getN <- function(x) sum(getUniques(x))
#track duplicate sequences and merge them
track <- cbind(out,sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.col))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "collapsed")
rownames(track) <- sample.names
head(track) # check to see the changes
tail(track)
write.table(cbind(rownames(track),track), "~/Documents/xoxo_article/files/18s/dada2_results_xoxo18s.txt", sep = "\t", row.names = F, quote = F) 
write.csv(as.matrix(seqtab.nochim),"~/Documents/xoxo_article/files/18s/asv_xoxo18s.csv")

rownames(seqtab.nochim) = sapply(strsplit(rownames(seqtab.nochim), "_S"), `[`, 1)

# Assign taxonomy using Silva v132 dada2 formatted 18s 'train sets' ####
#https://zenodo.org/records/1447330
taxaSilva <- assignTaxonomy(seqtab.nochim, "~/Documents/xoxo_article/files/18s/silva_132_18s/silva_132.18s.99_rep_set.dada2.fa.gz", multithread=TRUE)
colnames(taxaSilva) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Genbank")
unname(head(taxaSilva))
dim(taxaSilva) 

write.csv(as.matrix(taxaSilva),"~/Documents/xoxo_article/files/18s/taxa_xoxo18_genbank.csv")

##add species
#taxaSilva.spe <- addSpecies(taxaSilva, "~/Documents/xoxo_article/files/18s/silva132_db/silva_species_assignment_v132.fa.gz")
#colnames(taxaSilva.spe) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Remove Genbank
taxaSilva.spe <- taxaSilva
taxaSilva.spe <- taxaSilva.spe[,colnames(taxaSilva.spe) != "Genbank"]
dim(taxaSilva.spe)

# Inspect the taxonomic assignments 
taxaSilva.print <- taxaSilva.spe # Removing sequence rownames for display only
rownames(taxaSilva.print) <- NULL
head(taxaSilva.print)
dim(taxaSilva.print)

write.csv(as.matrix(taxaSilva.print),"~/Documents/xoxo_article/files/18s/taxa_xoxo18_march24.csv")

# Import the biosample (metadata) file ####
# It contains information about samples and treatments
map = read.csv("~/Documents/xoxo_article/files/metadata_xoxo_mi_18months_qpcr_euk_bac.csv")
rownames(map) = map$MI_ID
dim(map)

# Make phyloseq object ####
ps = phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), tax_table(taxaSilva.spe), sample_data(map))

#Corrections ####
#change month 19 to 18 (correction)
levels(as.factor(sample_data(ps)$MONTH))
which(sample_data(ps)$MONTH == "19")
sample_data(ps)[52,4]
sample_data(ps)[52,4] = stringr::str_replace(as.vector(sample_data(ps)[52,4]), '19', '18')
sample_data(ps)[52,4]

# save
saveRDS(ps, file="~/Documents/xoxo_article/files/18s/ps_xoxo_18s.rds")
save.image("~/Documents/xoxo_article/files/18s/a1_xoxo_18s_dada2.RData")

