###############################################################
# Cleaning and denoising 18S data                             #
# Data: Miseq-18S - xoxo                                      #
# Mona Parizadeh - May 2025                                   #
###############################################################

# Load libraries ####
library(phyloseq); packageVersion("phyloseq") 
library(vegan); packageVersion("vegan") 
library(ggplot2); packageVersion("ggplot2") 
library(tidyverse); packageVersion("tidyverse") 
library(BiocGenerics)
library(SummarizedExperiment)

# Import data #### 
setwd("~/Documents/xoxo_article/files/18s/")
ps = readRDS("ps_xoxo_18s.rds"); ps

# Explore data ####
nsamples(ps)
ntaxa(ps)
rank_names(ps)
sample_variables(ps)

#number of seq per sample 
summary(sample_sums(ps))  #or: summary(apply(comm,1,sum))
sd(sample_sums(ps), na.rm=TRUE)/sqrt(length(sample_sums(ps)[!is.na(sample_sums(ps))])) #SE
head(sort(sample_sums(ps),TRUE))
hist(sample_sums(ps))
hist(sample_sums(ps), main="Histogram: Read Counts", xlab="Total Reads", 
     border="blue", col="grey", las=1, breaks=10,xlim=c(0,90000),labels = TRUE)


taxa_names(ps) = paste0("ASV", seq(ntaxa(ps))) #replace sequence w/ ASV
tax_table(ps)
as.factor(tax_table(ps)[,4]) #order (fungi etc.)

#species richness per sample
summary(specnumber(otu_table(ps))) #species richness or number of ASVs per sample (or: summary(apply(comm[]>0,1,sum)))
sd(specnumber(otu_table(ps)), na.rm=TRUE)/
  sqrt(length(specnumber(otu_table(ps)[!is.na(specnumber(otu_table(ps)))]))) #SE

#number of seqs per ASV
summary(taxa_sums(ps)) 
sd(taxa_sums(ps), na.rm=TRUE)/
  sqrt(length(taxa_sums(ps)[!is.na(taxa_sums(ps))])) #SE

#Order Months
sample_data(ps)$MONTH <- factor(as.numeric(sample_data(ps)$MONTH), labels = c("1","2","3","4","5","6","7","8","9",
                                                                           "10","11","12","13","14","15","16","17","18"))
sample_data(ps)$MONTH
# Divide months into three groups (1-6, 7-12, 13-18) and create a new variable
sample_data(ps)$MONTH_GROUP <- cut(as.numeric(sample_data(ps)$MONTH), breaks = c(0, 6, 12, 18), labels = c("1-6 Months", "7-12 Months", "13-18 Months"))

#Mock samples ####
mock.ID = sample_names(ps)[grep("Mock.",sample_data(ps)$SAMPLE_ID)]
mock = subset_samples(ps, sample_names(ps) %in% mock.ID)
mock = prune_taxa(taxa_sums(mock)>0, mock); mock
mock %>% 
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% #remove NAs
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  plot_bar(fill="Genus") +
  labs(title = "Eukaryotic genera present in the mock samples") +
  xlab("Mock samples") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")
mock %>% 
  tax_glom(taxrank = "Species", NArm=TRUE) %>% #remove NAs
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  plot_bar(fill="Species") +
  labs(title = "Eukaryotic species present in the mock samples") +
  xlab("Mock samples") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")
mock %>%
  tax_glom(taxrank = "Species") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() 

#Investigating if species in mock samples are present in other samples (may check for the ASVs as well, if necessary)
ps.melt = psmelt(ps)

levels(as.factor(ps.melt[which(ps.melt$Species=="Blastocystis_sp._subtype_3"),]$OTU)) 
levels(as.factor(ps.melt[which(ps.melt$Species=="Blastocystis_sp._subtype_3"),]$INFANT_ID))

levels(as.factor(ps.melt[which(ps.melt$Species=="Candida_orthopsilosis"),]$OTU))
levels(as.factor(ps.melt[which(ps.melt$Species=="Candida_orthopsilosis"),]$INFANT_ID))

levels(as.factor(ps.melt[which(ps.melt$Species=="Candida_albicans"),]$OTU))
levels(as.factor(ps.melt[which(ps.melt$Species=="Candida_albicans"),]$INFANT_ID))

#Remove mock samples ####
ps.noMock = subset_samples(ps, !sample_names(ps) %in% mock.ID)
ps.noMock = prune_taxa(taxa_sums(ps.noMock)>0, ps.noMock); ps.noMock

#Explore samples per infant ####
#freq shows duplicates, so the real number of samples per infant is freq + 1
inf_samp = count(sample_data(ps.noMock)$INFANT_ID[duplicated(sample_data(ps.noMock)$INFANT_ID)]) 
summary(inf_samp$freq)
inf_samp[which(inf_samp$freq == summary(inf_samp$freq)[1]),]
#XC10
subset_samples(ps.noMock, sample_data(ps.noMock)$INFANT_ID == "XC10")
sample_data(subset_samples(ps.noMock, sample_data(ps.noMock)$INFANT_ID == "XC10"))$MONTH_GROUP
#Later remove infant XC10 for statistical analyses, which has only one sample during the 2nd 6-months and no samples during the 3rd 6-months of sampling and mode of birth and sex are NA

#Explore samples & taxa ####
#Explore taxanomy 
rank_names(ps.noMock)
get_taxa_unique(ps.noMock, "Kingdom") #or levels(as.factor(tax_table(ps.noMock)[,1]))
get_taxa_unique(ps.noMock, "Phylum")
get_taxa_unique(ps.noMock, "Class") #should remove Chloroplastida
#Holozoa: fungi; Nucletmycea: fungi; Chloroplastida: plants; Alveolata: protists; Rhizaria: unicellular eukaryotes; 
#Tubulinea: Amoeba; Protosporangiida: Amoeba; Stramenopiles: eukaryotes (SAR); Gracilipodida: Amoeba; Discosea: Amoeba;
#Aphelidea: fungi; Schizoplasmodiida: Amoeba; uncultured: - 
get_taxa_unique(ps.noMock, "Order") #should remove Metazoa_(Animalia)
get_taxa_unique(ps.noMock, "Family") 
get_taxa_unique(ps.noMock, "Genus")
get_taxa_unique(ps.noMock, "Species")

#Check for contaminants 
#(remove later)
unique(tax_table(ps.noMock)[,3]) #Class
subset_taxa(ps.noMock, Class=="Chloroplastida")  
subset_taxa(ps.noMock, (is.na(Class) | Class =="Chloroplastida" | Class=="uncultured"))

#explore the presence of "Zea" in samples per month
zea = subset_taxa(ps.noMock, Genus=="Zea")  
zea = prune_samples(sample_sums(zea)>0, zea); zea

zeamelt = zea %>% 
 # transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  psmelt()
ggplot(zeamelt, aes(x = MONTH, y = Abundance)) +
  geom_bar(stat = "identity") +
  labs(title = "Zea present in the samples/month",
       x = "Months",
       y = "Relative Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")
ggplot(zeamelt, aes(x = INFANT_ID, y = Abundance)) +
  geom_bar(stat = "identity") +
  labs(title = "Zea present in the samples/infant",
       x = "Infants",
       y = "Relative Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

unique(tax_table(ps.noMock)[,4]) #Order
subset_taxa(ps.noMock, Order=="Metazoa_(Animalia)") 

#number of seq per sample 
summary(sample_sums(ps.noMock))  
sd(sample_sums(ps.noMock), na.rm=TRUE)/sqrt(length(sample_sums(ps.noMock)[!is.na(sample_sums(ps.noMock))])) #SE
head(sort(sample_sums(ps.noMock),TRUE))
hist(sample_sums(ps.noMock),labels = TRUE, xlim=c(0,90000))

prune_taxa(taxa_sums(ps.noMock) > 10, ps.noMock)
prune_taxa(taxa_sums(ps.noMock) > 100, ps.noMock) 
filter_taxa(ps.noMock, function(x) sum(x > 3) > (0.1*length(x)), prune=TRUE) 
filter_taxa(ps.noMock, function (x) {sum(x > 0) > 1}, prune=TRUE)
filter_taxa(ps.noMock, function (x) {sum(x > 0) > 2}, prune=TRUE)
prune_samples(sample_sums(ps.noMock)>=1000, ps.noMock)
prune_samples(sample_sums(ps.noMock)>=500, ps.noMock)

#Remove undefined phyla ####
get_taxa_unique(ps.noMock, "Phylum") #check if there is NA phyla
subset_taxa(ps.noMock, is.na(Phylum) | Phylum %in% c("", "uncharacterized")) 
ps.noNaPhyla = subset_taxa(ps.noMock, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ps.noNaPhyla = prune_samples(sample_sums(ps.noNaPhyla)>0, ps.noNaPhyla)
ps.noNaPhyla
unique(tax_table(ps.noNaPhyla)[,2])
100-(ntaxa(ps.noNaPhyla)/ntaxa(ps.noMock))*100 
ntaxa(ps.noMock) - ntaxa(ps.noNaPhyla) #282 ASVs
nsamples(ps.noMock) - nsamples(ps.noNaPhyla)

#number of seq per sample 
summary(sample_sums(ps.noNaPhyla))
sd(sample_sums(ps.noNaPhyla), na.rm=TRUE)/sqrt(length(sample_sums(ps.noNaPhyla)[!is.na(sample_sums(ps.noNaPhyla))])) #SE
head(sort(sample_sums(ps.noNaPhyla),TRUE))
hist(sample_sums(ps.noNaPhyla))

#look for outliers ####
ps2 = ps.noNaPhyla

#nmds
sample_data(ps2)$MONTH = factor(sample_data(ps2)$MONTH, levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18))
ps.ra = transform_sample_counts(ps2, function(otu) otu/sum(otu)) #relative abudance
otu_table(ps.ra)[1:5, 1:5] 
nmds = ordinate(ps2, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps.ra, nmds, color = "MONTH", shape = "BIRTH") + 
  theme_bw() + geom_point() + ggtitle("nMDS") +
  geom_text(aes(label = SAMPLE_ID), check_overlap = TRUE, size = 4, nudge_y = -0.05) #nudge_x = -0.5
#outlier: "XC05-03"

#Alpha diversity ####
alpha.shn = estimate_richness(ps2, split=TRUE, measures="Shannon") 
summary(alpha.shn)
plot_richness(ps2, "INFANT_ID","MONTH", measures = "Shannon")
plot_richness(ps2, "INFANT_ID","MONTH", measures = "Shannon") +
  geom_text(aes(label = SAMPLE_ID), check_overlap = FALSE, size = 5)
#outlier: "XC09-09" #keep it

out = c("XC05-03")

#Remove outliers 
ps2.noOut = prune_samples(!sample_data(ps2)$SAMPLE_ID %in% out, ps2)
ps2.noOut = prune_taxa(taxa_sums(ps2.noOut)>0,ps2.noOut); ps2.noOut
ps2.noOut.ra = transform_sample_counts(ps2.noOut, function(otu) otu/sum(otu)) 
nmds1 = ordinate(ps2.noOut.ra, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps2.noOut.ra, nmds1, color = "MONTH", shape = "BIRTH") + 
  theme_bw() + geom_point() + ggtitle("nMDS") +
  geom_text(aes(label = SAMPLE_ID), check_overlap = TRUE, size = 5)

#Keep samples with at least 1000 reads (poor quality sequences) ####
ps3 = ps2.noOut
summary(sample_sums(ps3)) #1st Qu. 3960 reads/sample
ps3.1krds = prune_samples(sample_sums(ps3) >= 1000, ps3)
ps3.1krds = prune_taxa(taxa_sums(ps3.1krds) > 0, ps3.1krds); ps3.1krds
#lost samples
nsamples(ps3) - nsamples(ps3.1krds) 
ntaxa(ps3) - ntaxa(ps3.1krds) 

#Remove contaminants ####
# remove Chloroplastida 
subset_taxa(ps3.1krds, (Order=="Metazoa_(Animalia)"))
ps_clean = prune_samples(sample_sums(ps_clean) > 0, ps_clean); ps_clean
ps_clean = subset_taxa(ps3.1krds, (!is.na(Class) & Class!="Chloroplastida" & Class!="uncultured"))
ps_clean = prune_samples(sample_sums(ps_clean) > 0, ps_clean); ps_clean
ntaxa(ps3.1krds) - ntaxa(ps_clean)
nsamples(ps3.1krds) - nsamples(ps_clean)
#check
taxo_clean = as.data.frame(tax_table(ps_clean)); taxo_clean
unique(taxo_clean$Class)

# remove Metazoa_(Animalia)
ps_clean1 = subset_taxa(ps_clean, (Order!="Metazoa_(Animalia)"))
ps_clean1 = prune_samples(sample_sums(ps_clean1) > 0, ps_clean1); ps_clean1
tax_table(ps_clean1)
unique(tax_table(ps_clean1)[,4])

100-(ntaxa(ps_clean1)/ntaxa(ps3.1krds))*100 
100-(nsamples(ps_clean1)/nsamples(ps3.1krds))*100 

#Filtering ####
ps4 = ps_clean1

sort(taxa_sums(ps4));summary(taxa_sums(ps4))
sort(sample_sums(ps4));summary(sample_sums(ps4))
#Filter ASVs w/ less than 10 reads ####
ps4_seq10 = prune_taxa(taxa_sums(ps4) > 10, ps4)
ps4_seq10 = prune_samples(sample_sums(ps4_seq10)>0,ps4_seq10); ps4_seq10 
sort(taxa_sums(ps4_seq10));summary(taxa_sums(ps4_seq10))
sort(sample_sums(ps4_seq10)); summary(sample_sums(ps4))
100-(ntaxa(ps4_seq10)/ntaxa(ps4))*100
ntaxa(ps4) - ntaxa(ps4_seq10)

# Variance Stabilizing Transformation ####
# Create gm_mean function for variance stabilizing transformation (stabilize variants based on sample size)
gm_mean = function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

vst <- phyloseq_to_deseq2(ps5, ~BIRTH)
vst <- estimateSizeFactors(vst, geoMeans = apply(counts(vst), 1, gm_mean))
vst_blind <- DESeq2::varianceStabilizingTransformation(vst, blind = TRUE)
vst_blind_mat <- SummarizedExperiment::assay(vst_blind) # Extract transformed OTU table
vst_blind_mat <- t(vst_blind_mat)
vst_blind_mat[which(vst_blind_mat < 0)] <- 0
dists <- dist(t(assay(vst)))

# Create a new phyloseq object with VS-transformed counts
ps_vst <- phyloseq(otu_table(vst_blind_mat, taxa_are_rows = FALSE), sample_data(ps5), tax_table(ps5)); ps_vst

#Explore samples per infant 
#freq shows duplicates, so the real number of samples per infant is freq + 1
vst_inf_samp = count(sample_data(ps_vst)$INFANT_ID[duplicated(sample_data(ps_vst)$INFANT_ID)]) 
summary(vst_inf_samp$freq)

sample_data(subset_samples(ps_vst, sample_data(ps_vst)$BIRTH == "VAGINAL"))$MONTH #no remaining samples in the month 18 
sample_data(subset_samples(ps_vst, sample_data(ps_vst)$BIRTH == "Csec"))$MONTH

#FUNGI ####
ps.fun = subset_taxa(ps_vst, Order =="Fungi")
ps.fun = prune_samples(sample_sums(ps.fun)>0, ps.fun)
ps.fun
setdiff(sample_names(ps5),sample_names(ps.fun)) 
setdiff(as.vector(sample_data(ps5)$SAMPLE_ID),as.vector(sample_data(ps.fun)$SAMPLE_ID)) #XC05-09 has no fungi

#save
saveRDS(ps5, "~/Documents/xoxo_article/files/18s/ps_xoxo_18s_Cleaned.rds")
saveRDS(ps_vst, "~/Documents/xoxo_article/files/18s/ps_xoxo_18s_VST.rds")
save.image("~/Documents/xoxo_article/files/18s/a2_xoxo_18s_cleaning.RData")



