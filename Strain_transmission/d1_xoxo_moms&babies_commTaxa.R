###############################################################
# Find common species between moms and their infants          #
# Data: Metaphlan4 - xoxo                                     # 
# Mona Parizadeh - May 2025                                   #
###############################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") 
library(vegan); packageVersion("vegan") 
library(ggplot2); packageVersion("ggplot2") 
library(tidyverse); packageVersion("tidyverse")
library(dplyr);packageVersion("dplyr") 
library(patchwork)
library(ggh4x) #for facet_nested

#Import data ####
setwd("~/Documents/xoxo_article/files/metaphlan4/")
#import otu table
#result of runnung the metaphlanToPhyloseq_abs.ipynb script
comm = read.delim2("species_abundance_m4.csv", header = T, sep=",", check.names=FALSE) #check.names=FALSE to avoid changing dash to dot in colnames
colnames(comm); dim(comm)
comm[1:5,145:149] #check colnames
#remove clade_name column
comm1 = comm %>% select (-clade_name)
#comm1 = comm1 %>% select (-clade_taxid)
dim(comm1)
rownames(comm1) = comm1$Otu #rename rownames with otus
comm1[1:5,1:5]
comm2 = comm1 %>% select (-Otu) #remove otu colnames
is.numeric(comm2) #make sure if the dataframe is numeric
comm2 = cbind(sapply(comm2, as.numeric)) #convert to numeric
is.numeric(comm2)
rownames(comm2) = rownames(comm1)
colnames(comm2) = sapply(strsplit(colnames(comm2), "_S"), `[`, 1)
comm2[1:5,1:5]
dim(comm2)

#import taxa table
taxa = read.delim2("species_taxa_m4.csv", header = T, sep=",")
taxa[1:5,]; dim(taxa); class(taxa)
rownames(taxa) = taxa$Otu #rename ronames with otus
taxa = taxa %>% select (-Otu) %>% as("matrix")  #remove otu colnames

#import metadata
meta = read.delim2("../metadata_xoxo_mi_18months_mom.csv", header = T, sep = ",")
head(meta);dim(meta)
rownames(meta) = meta$MI_ID
diff = setdiff(rownames(meta),colnames(comm2)); diff #find differences 
meta1 = meta[!(rownames(meta) %in% diff),]; dim(meta1) #remove the extra samples

#add a column with the same IDs for infants and their correspinding mom
meta1$type = meta1$INFANT_ID; meta1$type
which(meta1$type == "Mom2")
meta1$type[139] = stringr::str_replace(meta1$type[139], 'Mom2', 'XC02')
which(meta1$type == "Mom4")
meta1$type[140] = stringr::str_replace(meta1$type[140], 'Mom4', 'XC04')
which(meta1$type == "Mom5")
meta1$type[141] = stringr::str_replace(meta1$type[141], 'Mom5', 'XC05')
which(meta1$type == "Mom6")
meta1$type[142] = stringr::str_replace(meta1$type[142], 'Mom6', 'XC06')
which(meta1$type == "Mom7")
meta1$type[143] = stringr::str_replace(meta1$type[143], 'Mom7', 'XC07')
which(meta1$type == "Mom9")
meta1$type[144] = stringr::str_replace(meta1$type[144], 'Mom9', 'XC09')
which(meta1$type == "Mom11")
meta1$type[145] = stringr::str_replace(meta1$type[145], 'Mom11', 'XC11')
which(meta1$type == "Mom15")
meta1$type[146] = stringr::str_replace(meta1$type[146], 'Mom15', 'XC15')
which(meta1$type == "Mom19")
meta1$type[147] = stringr::str_replace(meta1$type[147], 'Mom19', 'XC19')
meta1$type

#bulid phyloseq object
ps = phyloseq(otu_table(comm2, taxa_are_rows = TRUE), tax_table(taxa), sample_data(meta1))
ps = prune_taxa(taxa_sums(ps)>0, ps); ps 

#Corrections ####
#change month 19 to 18 (correction)
levels(as.factor(sample_data(ps)$MONTH))
which(sample_data(ps)$MONTH == "19")
sample_data(ps)[50,]
sample_data(ps)[50,5]
sample_data(ps)[50,5] = stringr::str_replace(as.vector(sample_data(ps)[50,5]), '19', '18')
sample_data(ps)[50,]

#number of seq per sample ####
summary(sample_sums(ps))  #or: summary(apply(comm,1,sum))
sd(sample_sums(ps), na.rm=TRUE)/sqrt(length(sample_sums(ps)[!is.na(sample_sums(ps))])) #SE
head(sort(sample_sums(ps),TRUE))
hist(sample_sums(ps))

#species richness per sample ####
summary(specnumber(otu_table(ps))) #species richness or number of ASVs per sample (or: summary(apply(comm[]>0,1,sum)))
sd(specnumber(otu_table(ps)), na.rm=TRUE)/
  sqrt(length(specnumber(otu_table(ps)[!is.na(specnumber(otu_table(ps)))]))) #SE

#number of seqs per ASV ####
summary(taxa_sums(ps)) 
sd(taxa_sums(ps), na.rm=TRUE)/
  sqrt(length(taxa_sums(ps)[!is.na(taxa_sums(ps))])) #SE

#Explore data
colnames(sample_data(ps))
sample_data(ps)[,2:3] #DATASET

#Explore taxanomy ####
rank_names(ps)
get_taxa_unique(ps, "Kingdom") #or levels(as.factor(tax_table(ps)[,1]))
subset_taxa(ps, Kingdom == "Bacteria")
tax_table(subset_taxa(ps, Kingdom == "Archaea"))
tax_table(subset_taxa(ps, Kingdom == "Eukaryota"))

ps %>% 
  tax_glom(taxrank = "Family", NArm=TRUE) %>% #remove NAs
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  plot_bar(fill="Kingdom") +
  xlab("Samples") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

#Mock samples ####
mock.ID = sample_names(ps)[grep("Mock.",sample_data(ps)$SAMPLE_ID)]
mock = subset_samples(ps, sample_names(ps) %in% mock.ID)
mock = prune_taxa(taxa_sums(mock)>0, mock); mock
mock %>% 
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% #remove NAs
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  plot_bar(fill="Genus") +
  facet_wrap(~ Kingdom) +
  labs(title = "Genera present in mock samples") +
  xlab("Mock samples") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")
mock %>% 
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% #remove NAs
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  plot_bar(fill="Genus") +
  labs(title = "Bacteria genera present in the mock samples") +
  xlab("Mock samples") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")
mock %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() 

tax_table(mock)[,7]

#Remove mock samples ####
ps.noMock = subset_samples(ps,!sample_names(ps) %in% mock.ID)
ps.noMock = prune_taxa(taxa_sums(ps.noMock)>0, ps.noMock); ps.noMock

length(setdiff(rownames(tax_table(ps)),rownames(tax_table(ps.noMock))))

#Order Months
sample_data(ps.noMock)$MONTH <- factor(as.numeric(sample_data(ps.noMock)$MONTH), 
                                       labels = c("0","1","2","3","4","5","6","7","8","9",
                                                  "10","11","12","13","14","15","16","17","18"))
sample_data(ps.noMock)$MONTH
# Divide months into four groups (0 for moms and 1-6, 7-12, 13-18 for infants) and create a new variable
sample_data(ps.noMock)$MONTH_GROUP <- cut(as.numeric(sample_data(ps.noMock)$MONTH), breaks = c(0, 6, 12, 18), labels = c("1-6 Months", "7-12 Months", "13-18 Months"))


#check for Eukaryota (only 4)
tax_table(subset_taxa(ps.noMock, tax_table(ps.noMock)[,1]=="Eukaryota"))

ps.noMock %>%
  tax_glom(taxrank = "Family", NArm=TRUE) %>% #remove NAs
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  plot_bar(fill="Phylum") +
  labs(title = "Phyla") +
  xlab("Samples") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

ps.noMock %>%
  tax_glom(taxrank = "Family", NArm=TRUE) %>% #remove NAs
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  plot_bar(fill="Phylum") +
  facet_grid(. ~ type+DATASET, scales = "free") +
  labs(title = "Phyla") +
  xlab("MONTH_GROUP") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

#keep only bacteria
ps_bac = subset_taxa(ps.noMock, Kingdom == "Bacteria")
ps_bac = prune_taxa(taxa_sums(ps_bac)>0,ps_bac); ps_bac
(1-(ntaxa(ps_bac)/ntaxa(ps.noMock)))*100

ps_bac %>%
  tax_glom(taxrank = "Family", NArm=TRUE) %>% #remove NAs
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  plot_bar(fill="Class") +
  facet_wrap(~ DATASET, scales = "free_x") +
  labs(title = "Bacteria present in samples") +
  xlab("Samples") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")

ps_bac = subset_taxa(ps_bac, Phylum != "Bacteria_unclassified"); ps_bac

#explore taxa
sort(get_taxa_unique(ps_bac, "Phylum"))  
sort(get_taxa_unique(ps_bac, "Class"))   #no NA
sort(get_taxa_unique(ps_bac, "Order"))   #no NA
sort(get_taxa_unique(ps_bac, "Family"))  #no NA
sort(get_taxa_unique(ps_bac, "Genus"))   #no NA
sort(get_taxa_unique(ps_bac, "Species")) #no NA


#number of seq per sample 
summary(taxa_sums(ps_bac))
summary(sample_sums(ps_bac))
sd(sample_sums(ps_bac), na.rm=TRUE)/sqrt(length(sample_sums(ps_bac)[!is.na(sample_sums(ps_bac))])) #SE
head(sort(sample_sums(ps_bac),TRUE))
hist(sample_sums(ps_bac)) #not much different from ps histogram

prune_taxa(taxa_sums(ps_bac) > 10, ps_bac)
prune_taxa(taxa_sums(ps_bac) > 100, ps_bac) 
filter_taxa(ps_bac, function(x) sum(x > 3) > (0.1*length(x)), prune=TRUE) 
filter_taxa(ps_bac, function (x) {sum(x > 0) > 1}, prune=TRUE)
filter_taxa(ps_bac, function (x) {sum(x > 0) > 2}, prune=TRUE)
prune_samples(sample_sums(ps_bac)>=100, ps_bac)
prune_samples(sample_sums(ps_bac)>=80, ps_bac)
prune_samples(sample_sums(ps_bac)>=50, ps_bac)

#look for outliers ####
sample_data(ps_bac)$MONTH = factor(sample_data(ps_bac)$MONTH, levels = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18))
sample_data(ps_bac)$MONTH # 0 is for moms
# Divide months into three groups (1-6, 7-12, 13-18) and create a new variable
sample_data(ps_bac)$MONTH_GROUP <- cut(as.numeric(sample_data(ps_bac)$MONTH), breaks = c(0, 1, 7, 13, 19), 
                                       labels = c("Mothers", "1-6 Month Infants", "7-12 Month Infants", "13-18 Month Infants"))
#1. nmds ####
nmds = ordinate(ps_bac, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps_bac, nmds, color = "type", shape = "DATASET") + 
  theme_bw() + geom_point(size = 3) + ggtitle("nMDS") +
  geom_text(aes(label = SAMPLE_ID), check_overlap = TRUE, size = 2, nudge_y = -0.05)

#2. alpha diversity ####
alpha.shn = estimate_richness(ps_bac, split=TRUE, measures="Shannon") 
summary(alpha.shn)
plot_richness(ps_bac, "INFANT_ID","type", measures = c("Shannon","Simpson")) +
  geom_text(aes(label = SAMPLE_ID), check_overlap = FALSE, size = 3)

adiv <- data.frame(
  "Simpson" = phyloseq::estimate_richness(ps_bac, measures = "Simpson"),
  "Shannon" = phyloseq::estimate_richness(ps_bac, measures = "Shannon"))
head(adiv)

# warning: The data you have provided does not have any singletons.
# This means you should not use the output of Metaphlan to estimate richness (eg. Chao S1). 
# However, you shouldn't have been doing that with the output of other methods either, as the high levels of FP singletons made richness estimates wrong anyway. Right now, I don't think a method exists that can make valid richness estimates from high-throughput amplicon data due to the difficulty of calling singletons accurately, and the sensitivity of richness estimation to the number of singletons.
# Other measures of diversity that aren't totally reliant on singletons, eg. Shannon/Simpson, are valid to use, and you can ignore the warning in phyloseq when calculating those measures.

adiv %>%
  gather(key = metric, value = value, c("Simpson", "Shannon")) %>%
  mutate(metric = factor(metric, levels = c("Simpson", "Shannon"))) %>%
  ggplot(aes(x = metric, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = metric), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none")

# While Simpson's index cares more about relative abundances, the Shannon index cares more about species richness
# So the importance of rare species decreases in order species richness > Shannon index > Simpson index

#Outlier removal based on alpha diversity####
plot_richness(ps_bac, x = "DATASET", color = "type", measures = c("Simpson", "Shannon")) +
  geom_boxplot(outlier.color = NA) +
  geom_text(aes(label = SAMPLE_ID)) +
  #geom_jitter(aes(color = BIRTH), height = 0, width = .2) +
  theme_bw() 
hist(adiv$Shannon)
hist(adiv$Simpson)

shn_bbs_bac = estimate_richness(ps_bac, split=TRUE, measures="Shannon") 

#samples w/ very high alpha diversity
div.out = "S0063-0016" #from infant analysis 
#Remove outliers 
ps_bac.div01 = prune_samples(!rownames(sample_data(ps_bac)) %in% div.out, ps_bac)
ps_bac.div01 = prune_taxa(taxa_sums(ps_bac.div01)>0,ps_bac.div01)
ps_bac.div01

ps2=ps_bac.div01
sample_data(ps2)[,c(2,3,25)]

#number of seq per sample ####
summary(sample_sums(ps2))  #or: summary(apply(comm,1,sum))
sd(sample_sums(ps), na.rm=TRUE)/sqrt(length(sample_sums(ps)[!is.na(sample_sums(ps))])) #SE
mean(table(sample_data(ps2)[,3]))

#samples per infant ####
#mom vs baby
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$DATASET == "MOMS"))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$DATASET == "BABIES"))$INFANT_ID)))

#Remove infants XC10 #### mom
sample_data(subset_samples(ps2, sample_data(ps2)$INFANT_ID == "XC10"))$MONTH_GROUP
#XC10 has only one sample during the 2nd 6-months and no samples during the 3rd 6-months of sampling and mode of birth and sex are NA
#remove it for any statistical analyses (any analyses other than characterization)
psnoXC10 = subset_samples(ps2, sample_data(ps2)$INFANT_ID != "XC10")
psnoXC10 = prune_taxa(taxa_sums(psnoXC10)>0, psnoXC10); psnoXC10
ntaxa(ps2) - ntaxa(psnoXC10)
nsamples(ps2) - nsamples(psnoXC10)


#Alpha diversity measures ####
ps3=psnoXC10

sample_data(ps3)$DATASET = factor(sample_data(ps3)$DATASET, levels = c("MOMS","BABIES"))
sample_data(ps3)$MONTH=factor(sample_data(ps3)$MONTH, levels=c(0,1,2,3,4,5,6,7,8,9,
                                                               10,11,12,13,14,15,16,17,18))

#shannon
shn.rich = cbind(estimate_richness(ps3, measures = 'shannon'),
                 sample_data(ps3))
ggplot(shn.rich, aes(x = MONTH_GROUP, y = Shannon, color=INFANT_ID)) +  
  geom_boxplot()

summary(shn.rich)
sd(shn.rich$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich$Shannon[!is.na(shn.rich$Shannon)])) #SE
compare_means(Shannon ~ MONTH_GROUP, data = shn.rich, method = "wilcox.test", p.adjust.method = "holm") 

#simpson
smp.rich = cbind(estimate_richness(ps3,measures = 'simpson'),
                 sample_data(ps3))
ggplot(smp.rich, aes(x = DATASET, y = Simpson, color=INFANT_ID)) +  
  geom_boxplot()

library(ggpubr); packageVersion("ggpubr") #‘0.6.0’
#shannon
#check for normality
shap.shn = shapiro.test(shn.rich$Shannon); shap.shn
shap.shn$p.value < 0.0001 #normally distributed
compare_means(Shannon ~ DATASET, data = shn.rich, method = "wilcox.test", p.adjust.method = "holm") 

#dataset
compare_means(Shannon ~ DATASET, data = shn.rich, method = "wilcox.test", p.adjust.method = "holm") 
shn.rich.mom = shn.rich[which(shn.rich$DATASET == "MOMS"),]; dim(shn.rich.mom)
shn.rich.inf = shn.rich[which(shn.rich$DATASET == "BABIES"),]; dim(shn.rich.inf)
summary(shn.rich.mom$Shannon)
sd(shn.rich.mom$Shannon, na.rm=TRUE)/sqrt(length(shn.rich.mom$Shannon[!is.na(shn.rich.mom$Shannon)])) #SE
summary(shn.rich.inf$Shannon)
sd(shn.rich.inf$Shannon, na.rm=TRUE)/sqrt(length(shn.rich.inf$Shannon[!is.na(shn.rich.inf$Shannon)])) #SE

#group
compare_means(Shannon ~ MONTH_GROUP, data = shn.rich, method = "wilcox.test", p.adjust.method = "holm") 
shn.rich.mom0 = shn.rich[which(shn.rich$MONTH_GROUP == "Mothers"),]; dim(shn.rich.mom0)
summary(shn.rich.mom0$Shannon)
sd(shn.rich.mom0$Shannon, na.rm=TRUE)/sqrt(length(shn.rich.mom0$Shannon[!is.na(shn.rich.mom0$Shannon)])) #SE
shn.rich.m6 = shn.rich[which(shn.rich$MONTH_GROUP == "1-6 Month Infants"),]; dim(shn.rich.m6)
summary(shn.rich.m6$Shannon)
sd(shn.rich.m6$Shannon, na.rm=TRUE)/sqrt(length(shn.rich.m12$Shannon[!is.na(shn.rich.m12$Shannon)])) #SE
shn.rich.m12 = shn.rich[which(shn.rich$MONTH_GROUP == "7-12 Month Infants"),]; dim(shn.rich.m12)
summary(shn.rich.m12$Shannon)
sd(shn.rich.m12$Shannon, na.rm=TRUE)/sqrt(length(shn.rich.m12$Shannon[!is.na(shn.rich.m12$Shannon)])) #SE
shn.rich.m18 = shn.rich[which(shn.rich$MONTH_GROUP == "13-18 Month Infants"),]; dim(shn.rich.m18)
summary(shn.rich.m18$Shannon)
sd(shn.rich.m18$Shannon, na.rm=TRUE)/sqrt(length(shn.rich.m18$Shannon[!is.na(shn.rich.m18$Shannon)])) #SE

#birth
compare_means(Shannon ~ BIRTH, data = shn.rich, method = "wilcox.test", p.adjust.method = "holm") 
shn.rich.vg = shn.rich[which(shn.rich$BIRTH == "VAGINAL"),]; dim(shn.rich.vg)
shn.rich.cs = shn.rich[which(shn.rich$BIRTH == "Csec"),]; dim(shn.rich.cs)
summary(shn.rich.vg$Shannon)
sd(shn.rich.vg$Shannon, na.rm=TRUE)/sqrt(length(shn.rich.vg$Shannon[!is.na(shn.rich.vg$Shannon)])) #SE
summary(shn.rich.cs$Shannon)
sd(shn.rich.cs$Shannon, na.rm=TRUE)/sqrt(length(shn.rich.cs$Shannon[!is.na(shn.rich.cs$Shannon)])) #SE

#excluded xc10 (ps3)
#shannon
img_shn_dataset = ggplot(shn.rich, aes(x = DATASET, y = Shannon, color=DATASET, fill=DATASET)) + #alpha: color intensity)) 
  theme_bw() +
  geom_violin(trim = FALSE, alpha = 0.7) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black") +  # Smaller boxplot inside the violin plot
  scale_color_manual(name = "DATASET", values = c("#acb872","#4c851b"), 
                     labels = c(MOMS = "Mothers",BABIES = "Infants")) +
  scale_fill_manual(name = "DATASET", values = c("#acb872","#4c851b"), 
                    labels = c(MOMS = "Mothers",BABIES = "Infants")) +
  scale_x_discrete(labels = c(MOMS = "Mothers",BABIES = "Infants")) +
  ggpubr::stat_compare_means(comparisons = list(c("MOMS","BABIES")), method = "wilcox.test",
                             size = 6, label.y = 4.5,aes(label = paste0("p = ", after_stat(p.signif))))+
  geom_jitter(col = "black", alpha=0.5, size = 0.5) +
  labs(x = "Source of Samples", y = "Bacterial Diversity (Shannon Index)", color = "Source of Samples") +
  theme(axis.title=element_text(face="bold",size=18, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=12, color="black", angle = 30),
        axis.text.y=element_text(size=14, color="black"),
        strip.background=element_rect(colour="black")) ; img_shn_dataset 

#FIGURE S3.A####
img_shn_brt_gr = ggplot(shn.rich, aes(x = MONTH_GROUP, y = Shannon, color=MONTH_GROUP, fill=MONTH_GROUP)) + #alpha: color intensity)) 
  theme_bw() +
  geom_violin(trim = FALSE, alpha = 0.7) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black") +
  facet_grid(~BIRTH, scales = "free",space = "free", 
             labeller = labeller(BIRTH=c(Csec = "C-section",VAGINAL= "Vaginal Birth")))+
  scale_color_manual(name = "MONTH_GROUP", values = c("#4c851b","#113a70","#5d8dc9","#a0bbde")) +#"#acb872"
  scale_fill_manual(name = "MONTH_GROUP", values = c("#acb872","#113a70","#5d8dc9","#a0bbde") )+
  geom_jitter(col = "black", alpha=0.5, size = 0.5) +
  # Add significance annotations
  geom_segment(aes(x = 1, xend = 3, y = 4.4, yend = 4.4), size = 0.2,col="black") +
  geom_segment(aes(x = 1, xend = 1, y = 4.38, yend = 4.4), size = 0.2, col="black") +
  geom_segment(aes(x = 3, xend = 3, y = 4.38, yend = 4.4), size = 0.2, col="black") +
  geom_text(aes(x = 2, y = 4.5, label = paste("****")), size = 11, col="black") +
  labs(x = "", y = "Bacterial Diversity (Shannon Index)") +
  theme(axis.title=element_text(face="bold",size=25, color="black"),
        legend.position="none",
        strip.text = element_text(size=18,face = "bold"),
        axis.text.x=element_text(size=18, color="black", angle = 30, hjust = 1, face = "bold"),
        axis.text.y=element_text(size=14, color="black"),
        strip.background=element_rect(colour="black",fill = "white")); img_shn_brt_gr 


ggplot(shn.rich, aes(x = MONTH_GROUP, y = Shannon, color=MONTH_GROUP, fill=MONTH_GROUP)) + #alpha: color intensity)) 
  theme_bw() +
  geom_boxplot(col="black",size=0.2) +
  geom_jitter(col = "black", alpha=0.5, size = 0.5) +
  geom_line(aes(group = type, colour = type)) +
  labs(x = "\nSource of Samples", y = "Bacterial Diversity (Shannon Index)", color = "Source of Samples") +
  theme(axis.title=element_text(face="bold",size=18, color="black"),
        legend.position="none",
        strip.text = element_text(size=14,face = "bold"),
        axis.text.x=element_text(size=12, color="black"),#, angle = 30
        axis.text.y=element_text(size=14, color="black"),
        strip.background=element_rect(colour="black"))

ggplot(shn.rich, aes(x = DATASET, y = Shannon, color=DATASET, fill=DATASET)) + #alpha: color intensity)) 
  theme_bw() +
  geom_violin(trim = FALSE, alpha = 0.7) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black") +
  facet_grid(~BIRTH, scales = "free",space = "free", 
             labeller = labeller(DATASET=c(MOMS = "Mothers",BABIES = "Infants"))) +
  scale_color_manual(name = "GROUP", values = c("#4c851b","#113a70","#5d8dc9","#a0bbde")) +#"#acb872"
  scale_fill_manual(name = "GROUP", values = c("#acb872","#113a70","#5d8dc9","#a0bbde") )+
  geom_jitter(col = "black", alpha=0.5, size = 0.5) +
  labs(x = "\nSource of Samples", y = "Bacterial Diversity (Shannon Index)", color = "Source of Samples") +
  theme(axis.title=element_text(face="bold",size=18, color="black"),
        legend.position="none",
        strip.text = element_text(size=14,face = "bold"),
        axis.text.x=element_text(size=12, color="black"),#, angle = 30
        axis.text.y=element_text(size=14, color="black"),
        strip.background=element_rect(colour="black"))


# PERMANOVA ####
#make dataframe
df = as(sample_data(ps3), "data.frame")
#bray-curtis distance
dis = phyloseq::distance(ps3,  method = "bray")
#permanova
set.seed(412)
adns = adonis2(dis ~ DATASET/MONTH_GROUP * BIRTH/type, df, permutations = 999, 
               by="terms", strata = df$INFANT_ID) #distance = bray , strata = df$DATASET, 
adns

#Ordination ####
#nmds
ps_bac_bray <- ordinate(ps3, "NMDS", "bray")
plot_ordination(ps3, ps_bac_bray, type="samples", color="type", shape="DATASET") + 
  geom_point(size = 3) 

nmds = ordinate(ps3,  method = "NMDS", k = 2, try = 100, distance = "bray")
# To flip the NMDS along the x-axis (Axis 1)
nmds$points[,1] <- -nmds$points[,1]
# To flip the NMDS along the y-axis (Axis 2)
nmds$points[,2] <- -nmds$points[,2]

plot_ordination(ps3, nmds, type = "samples", color = "type", shape = "DATASET") + 
  geom_point() + 
  #facet_wrap(~BIRTH) +
  stat_ellipse(aes(group = BIRTH))

#define variables as factors
gr = get_variable(ps3, "MONTH_GROUP")
sample_data(ps3)$GROUP = factor(gr)
dat = get_variable(ps3, "DATASET")
sample_data(ps3)$DATASET = factor(dat)

#group
gr.dat = paste(gr, dat, sep = "")
#order as factor
gr.dat.fac = factor(gr.dat,levels=c("MothersMOMS","1-6 Month InfantsBABIES",
                                    "7-12 Month InfantsBABIES", "13-18 Month InfantsBABIES"))

nmds_bac = plot_ordination(
  physeq = ps3,                                                        
  ordination = nmds)+       
  stat_ellipse(geom = "polygon", aes(group = GROUP, fill = GROUP), 
               level = 0.95, linetype = "dashed", colour = "black", alpha = 0.3,show.legend = FALSE) + 
  geom_point(aes(colour = GROUP, shape = DATASET), size = 5) +   
  #scale_color_manual(name = "MONTH_GROUP", values = c("#4c851b","#113a70","#5d8dc9","#a0bbde")) + #"#acb872"
  scale_color_manual(name = "GROUP", values = c("#4c851b","#113a70","#5d8dc9","#a0bbde")) + #"#acb872"  
  scale_fill_manual(name = "GROUP", values = c("#4c851b","#113a70","#5d8dc9","#a0bbde")) +
  scale_shape_manual(name="DATASET", values = c(19,17),
                     labels = c(MOMS = "Moms",BABIES = "Babies")) +
  #group by birth
  facet_grid(~BIRTH, scales = "free",space = "free", 
             labeller = labeller(BIRTH=c(Csec = "C-section",VAGINAL= "Vaginal Birth")))+
  theme_classic() +  
  theme(                             
    legend.text = element_text(size = 12),      
    legend.title = element_text(size = 12, face="bold"),
    legend.position = "right")+
  theme(axis.text.y.left = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face="bold"))+
  guides(fill = guide_legend(override.aes = list(shape = 21))) +   #fills legend points based on the fill command
  ggtitle(""); nmds_bac

pcoa = ordinate(ps3, method = "MDS", k = 2, try = 100, distance = "bray")

# Computing Bray-Curtis Dissimilarities and PCoA
comm_mat <- vegdist(otu_table(ps3), "bray")
PCoA_comm_mat <- capscale(comm_mat ~ 1, distance = "bray")
PCoA_comm_mat$CA$eig[1:3]/sum(PCoA_comm_mat$CA$eig)
PCoA_scores <- scores(PCoA_comm_mat)$sites

#FIGURE S3.B####  
gg.pcoa.gr = plot_ordination(
  physeq = ps3,                                                        
  ordination = pcoa) +
  #group by birth
  stat_ellipse(geom = "polygon", aes(group = GROUP, fill = GROUP), 
               level = 0.95, linetype = "dashed", colour = "black", alpha = 0.3,show.legend = FALSE) + 
  geom_point(aes(colour = gr.dat.fac, shape = gr.dat.fac), size = 5) +   
  scale_color_manual(name = "", values = c("#4c851b","#113a70","#5d8dc9","#a0bbde"),
                     labels = c(MothersMOMS = "Mothers","1-6 Month InfantsBABIES" = "1-6 Month Infants",
                                "7-12 Month InfantsBABIES" = "7-12 Month Infants",
                                "13-18 Month InfantsBABIES" = "13-18 Month Infants")) + #"#acb872"  
  scale_fill_manual(name = "Groups", values = c("#4c851b","#113a70","#5d8dc9","#a0bbde"),
                    labels = c(MothersMOMS = "Mothers","1-6 Month InfantsBABIES" = "1-6 Month Infants",
                               "7-12 Month InfantsBABIES" = "7-12 Month Infants",
                               "13-18 Month InfantsBABIES" = "13-18 Month Infants")) +     
  scale_shape_manual(name="", values = c(19,17,17,17),
                     labels = c(MothersMOMS = "Mothers","1-6 Month InfantsBABIES" = "1-6 Month Infants",
                                "7-12 Month InfantsBABIES" = "7-12 Month Infants",
                                "13-18 Month InfantsBABIES" = "13-18 Month Infants")) + 
  facet_grid(~BIRTH, scales = "free",space = "free", 
             labeller = labeller(BIRTH=c(Csec = "C-section",VAGINAL= "Vaginal")))+
  theme_classic() +  
  theme(axis.title=element_text(face="bold",size=25, color="black"),                             
        legend.text = element_text(size = 20),      
        legend.title = element_text(size = 22, face="bold"),
        legend.position = "right",
        strip.text = element_text(size=18,face = "bold"),
        axis.text.x=element_text(size=18, color="black"),
        axis.text.y=element_text(size=18, color="black"),
        strip.background=element_rect(colour="black"))+
  guides(fill = guide_legend(override.aes = list(shape = 21))); gg.pcoa.gr

plot_ordination(
  physeq = ps3,                                                        
  ordination = pcoa) +
  #group by birth
  stat_ellipse(aes(group = DATASET, colour = DATASET), level = 0.95, linetype = "dashed", alpha = 0.6) + 
  geom_point(aes(colour = GROUP, shape = BIRTH), size = 3) +  
  theme_classic() +  
  theme(                             
    legend.text = element_text(size = 12),      
    legend.title = element_text(size = 12, face="bold"),
    legend.position = "right")+
  theme(axis.text.y.left = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face="bold"))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))    #fills legend points based on the fill command

#RENAME infants' IDs ####
meta.ps3 <- sample_data(ps3)
sample_data(ps3)$type
# Modify INFANT_ID selectively
meta.ps3$type <- ifelse(meta.ps3$type == "XC02", "P01",
                        ifelse(meta.ps3$type == "XC04", "P02",
                               ifelse(meta.ps3$type == "XC05", "P03",
                                      ifelse(meta.ps3$type == "XC06", "P04",
                                             ifelse(meta.ps3$type == "XC07", "P05",
                                                    ifelse(meta.ps3$type == "XC09", "P06",
                                                           ifelse(meta.ps3$type == "XC11", "P08",
                                                                  ifelse(meta.ps3$type == "XC15", "P09",
                                                                         ifelse(meta.ps3$type == "XC19", "P10",
                                                                                meta.ps3$type)))))))))  # Leave others unchanged
meta.ps3$type
sample_data(ps3) <- meta.ps3 

sample_data(ps3)$type 

#RENAME bacterial phyla ####
taxa.ps3 = as.data.frame(tax_table(ps3))  

taxa.ps3$Phylum <- dplyr::case_when(taxa.ps3$Phylum == "Actinobacteria" ~ "Actinomycetota",
                                    taxa.ps3$Phylum == "Bacteroidetes" ~ "Bacteroidota", 
                                    taxa.ps3$Phylum == "Elusimicrobia" ~ "Elusimicrobiota",
                                    taxa.ps3$Phylum == "Firmicutes" ~ "Bacillota",
                                    taxa.ps3$Phylum == "Fusobacteria" ~ "Fusobacteriota",  
                                    taxa.ps3$Phylum == "Proteobacteria" ~ "Pseudomonadota",
                                    taxa.ps3$Phylum == "Tenericutes" ~ "Mycoplasmatota",
                                    taxa.ps3$Phylum == "Verrucomicrobia" ~ "Verrucomicrobiota",
                                    TRUE ~ taxa.ps3$Phylum)  # Keep the rest unchanged

taxa.ps3$Phylum = sub("_"," ", taxa.ps3$Phylum)
taxa.ps3$Phylum = sub("-"," ", taxa.ps3$Phylum); taxa.ps3$Phylum
unique(taxa.ps3$Phylum)
tax_table(ps3) <- tax_table(as.matrix(taxa.ps3))
tax_table(ps3) 
ps3
sample_data(ps3)$type

#Colors ####  
cols <- function(a) image(1:length(a), 1, as.matrix(1:length(a)),
                          col=a, axes=T , xlab="", ylab="")
phyla_colors <-c("Actinomycetota"="#5d8dc9","Bacteroidota"="#df65b0","Candidatus Melainabacteria"="#88419d", 
                 "Candidatus Saccharibacteria"="#e9c1f5","Elusimicrobiota"="#a1d99b","Bacillota"="#31a354",                 
                 "Fusobacteriota"="cyan4","Pseudomonadota"="yellow","Spirochaetes"="#fc4e2a",               
                 "Mycoplasmatota"="salmon","Verrucomicrobiota"="#980043")

cols(phyla_colors)

# Phylum ####
comm = otu_table(ps3)
taxa = as.data.frame(tax_table(ps3))
metadata = as(sample_data(ps3),"data.frame")
comm.taxo_phlm <- aggregate(comm, by=list(class=taxa$Phylum), sum)
rownames(comm.taxo_phlm) <- comm.taxo_phlm[,1]
comm.taxo_phlm <- comm.taxo_phlm[,-1]
t(comm.taxo_phlm) -> comm.taxo_phlm
colnames(comm.taxo_phlm)
dim(comm.taxo_phlm)

#Barplot 
ps3_aggregated <- ps3 %>%
  tax_glom(taxrank = "Species", NArm = TRUE) %>% # Remove NAs and group by species
  transform_sample_counts(function(x) {x / sum(x)}) %>% # Transform to relative abundance
  psmelt() %>% # Melt phyloseq object into a data frame
  group_by(INFANT_ID, Phylum, DATASET, type, MONTH_GROUP,BIRTH) %>% # Group by INFANT_ID, Phylum, DATASET, type, and MONTH_GROUP
  summarize(Relative_Abundance = sum(Abundance), .groups = "drop") %>% # Sum abundances
  group_by(DATASET, type, MONTH_GROUP,BIRTH) %>% # Group by DATASET, type, and MONTH_GROUP
  mutate(Relative_Abundance = Relative_Abundance / sum(Relative_Abundance)) # Scale to 0-1
ps3_aggregated

#Phyla
gg.phyla=ggplot(ps3_aggregated, aes(x = INFANT_ID, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(aes(MONTH_GROUP),stat = "identity") +
  facet_wrap(~ type, ncol = 9) +
  scale_fill_manual("Phyla", values = phyla_colors)+
  labs(
    x = "",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(axis.title=element_text(face="bold",size=22, color="black"),
        legend.position="right",
        strip.text = element_text(size=18,face = "bold"),
        axis.text.x=element_text(size=18, color="black", angle = 30, hjust=1),
        axis.text.y=element_text(size=14, color="black"),
        legend.title = element_text(face = "bold", size = 25),
        legend.text = element_text(size = 23,face="italic")); gg.phyla

#Vaginal ####
#-Phyla####
ps3_aggregated_vag <- ps3 %>%
  subset_samples(BIRTH == "VAGINAL") %>%  # Keep only C-section samples
  tax_glom(taxrank = "Species", NArm = TRUE) %>% # Remove NAs and group by species
  transform_sample_counts(function(x) {x / sum(x)}) %>% # Transform to relative abundance
  psmelt() %>% # Melt phyloseq object into a data frame
  group_by(INFANT_ID, Phylum, DATASET, type, MONTH_GROUP,BIRTH) %>% # Group by INFANT_ID, Phylum, DATASET, type, and MONTH_GROUP
  summarize(Relative_Abundance = sum(Abundance), .groups = "drop") %>% # Sum abundances
  group_by(DATASET, type, MONTH_GROUP,BIRTH) %>% # Group by DATASET, type, and MONTH_GROUP
  mutate(Relative_Abundance = Relative_Abundance / sum(Relative_Abundance)) # Scale to 0-1
ps3_aggregated_vag

#FIGURE S3.C - phyla ####
gg.phyla_vag=ggplot(ps3_aggregated_vag, aes(x = INFANT_ID, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(aes(MONTH_GROUP),stat = "identity",width = 0.9) +
  facet_wrap(~ type, nrow = 1) +
  scale_fill_manual("Phyla", values = phyla_colors)+
  labs(
    x = "",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(axis.title=element_text(face="bold",size=25, color="black"),
        legend.position="none",
        strip.text = element_text(size=22,face = "bold"),
        #axis.text.x=element_text(size=18, color="black", angle = 300, hjust=0,vjust = 0,face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=16, color="black"),
        strip.background=element_rect(colour="black"),
        legend.title = element_text(face = "bold", size = 25),
        legend.text = element_text(size = 23)); gg.phyla_vag

#Genera####
ps3_aggregated_vag_gen <- ps3 %>%
  subset_samples(BIRTH == "VAGINAL") %>%  # Keep only C-section samples
  tax_glom(taxrank = "Species", NArm = TRUE) %>% # Remove NAs and group by species
  transform_sample_counts(function(x) {x / sum(x)}) %>% # Transform to relative abundance
  psmelt() %>% # Melt phyloseq object into a data frame
  group_by(INFANT_ID, Genus, DATASET, type, MONTH_GROUP,BIRTH) %>% # Group by INFANT_ID, Genus, DATASET, type, and MONTH_GROUP
  summarize(Relative_Abundance = sum(Abundance), .groups = "drop") %>% # Sum abundances
  group_by(DATASET, type, MONTH_GROUP,BIRTH) %>% # Group by DATASET, type, and MONTH_GROUP
  mutate(Relative_Abundance = Relative_Abundance / sum(Relative_Abundance)) # Scale to 0-1
ps3_aggregated_vag_gen
unique(ps3_aggregated_vag_gen$Genus)

ps3_aggregated_vag_gen$Genus <- gsub("_unclassified", "", ps3_aggregated_csec_gen$Genus) # Remove "_unclassified"
unique(ps3_aggregated_csec_gen$Genus)

cols <- function(a) image(1:length(a), 1, as.matrix(1:length(a)),
                          col=a, axes=T , xlab="", ylab="")
genera = c("Akkermansia","Agathobaculum","Anaerobutyricum","Anaerostipes","Bacteroides", 
           "Bifidobacterium","Blautia", "Catenibacterium", "Clostridium", 
           "Collinsella","Dorea", "Enterococcus","Escherichia", 
           "Eubacterium", "Evtepia", "Faecalibacillus", "Faecalibacterium","Formosa", 
           "Fusicatenibacter", "Haemophilus", "Holdemanella","Intestinibacter", 
           "Klebsiella", "Lacrimispora","Lactococcus", "Lachnospira",
           "Mediterraneibacter", "Parabacteroides", "Phocaeicola", "Prevotella", "Romboutsia",
           "Roseburia", "Ruminococcus", "Sellimonas", "Streptococcus", "Turicibacter",
           "Veillonella", "Weissella")

genera_colors <- c("Akkermansia"="#67000d","Agathobaculum"="#bf4554","Anaerobutyricum"="#980043","Anaerostipes"="#7fcdbb","Bacteroides"="#1d91c0", 
                   "Bifidobacterium"="#dd1c77","Blautia"="#c7e9b4", "Catenibacterium"="#8aad76", "Clostridium"="cyan4", #"Clostridiaceae_unclassified"="cyan4",
                   "Collinsella"="#225ea8","Dorea"="cyan2", "Enterococcus"="#5d8dc9","Escherichia"="#feb24c", #"Eubacteriales_unclassified"="tomato",
                   "Eubacterium"="tomato", "Evtepia"="#ed3257", "Faecalibacillus"="#bf0429", "Faecalibacterium"="#fd8d3c","Formosa"="#006d2c", 
                   "Fusicatenibacter"="#bfde14", "Haemophilus"="#a84d78", "Holdemanella"="#f5b5c2","Intestinibacter"= "#b01eaa", #"GGB4456"=  "#74c476",
                   "Klebsiella"="#df65b0", "Lacrimispora"="#b0ae37","Lactococcus"="yellow", "Lachnospira"="#f5c731", #"Lachnospiraceae_unclassified"="#f5c731",
                   "Mediterraneibacter"="#31a354", "Parabacteroides"="#a1d99b", "Phocaeicola"="#fed976", "Prevotella"="#feaa76", "Romboutsia"="salmon",
                   "Roseburia"="#e31a1c", "Ruminococcus"="#e9c1f5", "Sellimonas"="#dd8cf5", "Streptococcus"="#88419d", "Turicibacter"="#76d1db",
                   "Veillonella"="#4e84c7", "Weissella"="#0969e0","Others"="#cfcccc")
cols(genera_colors)


ps3_aggregated_vag_gen$Genus.color <- ifelse(
  ps3_aggregated_vag_gen$Genus %in% genera,
  ps3_aggregated_vag_gen$Genus,
  "Others"
)

#FIGURE S3.C - genera ####
gg.genera_vag = ggplot(ps3_aggregated_vag_gen, aes(x = INFANT_ID, y = Relative_Abundance, fill = Genus.color)) +
  geom_bar(aes(MONTH_GROUP),stat = "identity",width = 0.95) +
  facet_wrap(~ type, nrow = 1) +
  scale_fill_manual("Genera", values = genera_colors)+
  labs(
    x = "",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(axis.title=element_text(face="bold",size=25, color="black"),
        legend.position="none",
        strip.text = element_blank(),
        axis.text.x=element_text(size=18, color="black", angle = 60, hjust=1,face = "bold"),
        axis.text.y = element_text(size=16, color="black"),
        legend.title = element_text(face = "bold", size = 25),
        legend.text = element_text(face="italic",size = 23));gg.genera_vag

#Csec ####
#-Phyla####
ps3_aggregated_csec <- ps3 %>%
  subset_samples(BIRTH == "Csec") %>%  # Keep only C-section samples
  tax_glom(taxrank = "Species", NArm = TRUE) %>% # Remove NAs and group by species
  transform_sample_counts(function(x) {x / sum(x)}) %>% # Transform to relative abundance
  psmelt() %>% # Melt phyloseq object into a data frame
  group_by(INFANT_ID, Phylum, DATASET, type, MONTH_GROUP,BIRTH) %>% # Group by INFANT_ID, Phylum, DATASET, type, and MONTH_GROUP
  summarize(Relative_Abundance = sum(Abundance), .groups = "drop") %>% # Sum abundances
  group_by(DATASET, type, MONTH_GROUP,BIRTH) %>% # Group by DATASET, type, and MONTH_GROUP
  mutate(Relative_Abundance = Relative_Abundance / sum(Relative_Abundance)) # Scale to 0-1
ps3_aggregated_csec

#FIGURE S3.D - phyla ####
gg.phyla_csec=ggplot(ps3_aggregated_csec, aes(x = INFANT_ID, y = Relative_Abundance, fill = Phylum)) +
  geom_bar(aes(MONTH_GROUP),stat = "identity",width = 0.9) +
  facet_wrap(~ type, nrow = 1) +
  scale_fill_manual("Phyla", values = phyla_colors)+
  labs(
    x = "",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(axis.title=element_text(face="bold",size=25, color="black"),
        legend.position="none",
        strip.text = element_text(size=22,face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=16, color="black"),
        strip.background=element_rect(colour="black"),
        legend.title = element_text(face = "bold", size = 27),
        legend.text = element_text(size = 25)); gg.phyla_csec

#Genera####
ps3_aggregated_csec_gen <- ps3 %>%
  subset_samples(BIRTH == "Csec") %>%  # Keep only C-section samples
  tax_glom(taxrank = "Species", NArm = TRUE) %>% # Remove NAs and group by species
  transform_sample_counts(function(x) {x / sum(x)}) %>% # Transform to relative abundance
  psmelt() %>% # Melt phyloseq object into a data frame
  group_by(INFANT_ID, Genus, DATASET, type, MONTH_GROUP,BIRTH) %>% # Group by INFANT_ID, Genus, DATASET, type, and MONTH_GROUP
  summarize(Relative_Abundance = sum(Abundance), .groups = "drop") %>% # Sum abundances
  group_by(DATASET, type, MONTH_GROUP,BIRTH) %>% # Group by DATASET, type, and MONTH_GROUP
  mutate(Relative_Abundance = Relative_Abundance / sum(Relative_Abundance)) # Scale to 0-1
ps3_aggregated_csec_gen
unique(ps3_aggregated_csec_gen$Genus)

ps3_aggregated_csec_gen$Genus.color <- ifelse(
  ps3_aggregated_csec_gen$Genus %in% genera,
  ps3_aggregated_csec_gen$Genus,
  "Others"
)

#FIGURE S3.D - genera ####
gg.genera_csec = ggplot(ps3_aggregated_csec_gen, aes(x = INFANT_ID, y = Relative_Abundance, fill = Genus.color)) +
  geom_bar(aes(MONTH_GROUP),stat = "identity", width = 0.9) +
  facet_wrap(~ type, nrow = 1) +
  scale_fill_manual("Genera", values = genera_colors)+
  labs(
    x = "",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(axis.title=element_text(face="bold",size=25, color="black"),
        legend.position="none",
        strip.text = element_blank(),
        axis.text.x=element_text(size=18, color="black", angle = 60, hjust=1,face = "bold"),
        axis.text.y = element_text(size=16, color="black"),
        legend.title = element_text(face = "bold", size = 27),
        legend.text = element_text(face="italic",size = 25));gg.genera_csec

#save ####
save.image("~/Documents/xoxo_article/files/metaphlan4/d1_xoxo_moms&babies_commTaxa.RData")

