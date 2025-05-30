###############################################################
# Cleaning and explanatory analysis of metagenomic data       #
# Data: Metaphlan4 - xoxo                                     # 
# Mona Parizadeh - May 2025                                   #
###############################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") 
library(vegan); packageVersion("vegan") 
library(ggplot2); packageVersion("ggplot2") 
library(tidyverse); packageVersion("tidyverse") 
library(dplyr);packageVersion("dplyr")
library(ggpubr); packageVersion("ggpubr")

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
#remove taxa from mom samples later when creating the ps object

#import metadata
meta = read.delim2("../metadata_xoxo_mi_18months_qpcr_euk_bac.csv", header = T, sep = ",")
head(meta);dim(meta)
rownames(meta) = meta$MI_ID
diff = setdiff(rownames(meta),colnames(comm2)); diff #find differences 
meta1 = meta[!(rownames(meta) %in% diff),]; dim(meta1) #remove the extra samples

#bulid phyloseq object
ps = phyloseq(otu_table(comm2, taxa_are_rows = TRUE), tax_table(taxa), sample_data(meta1))
ps = prune_taxa(taxa_sums(ps)>0, ps); ps #remove taxa from mom samples

#Corrections ####
#change month 19 to 18 (correction)
levels(as.factor(sample_data(ps)$MONTH))
which(sample_data(ps)$MONTH == "19")
sample_data(ps)[50,4]
sample_data(ps)[50,4] = stringr::str_replace(as.vector(sample_data(ps)[50,4]), '19', '18')
sample_data(ps)[50,4]

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
sample_data(ps)[,2] #DATASET

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
  labs(title = "Kingdoms present in mock and baby samples") +
  xlab("Baby samples") + ylab("Relative Abundance") +
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

#check for Eukaryota (only 4)
tax_table(subset_taxa(ps.noMock, tax_table(ps.noMock)[,1]=="Eukaryota"))

ps.noMock %>%
  tax_glom(taxrank = "Family", NArm=TRUE) %>% #remove NAs
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  plot_bar(fill="Phylum") +
  labs(title = "Phyla present in baby samples") +
  xlab("Baby samples") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

#keep only bacteria
ps_bbs_bac = subset_taxa(ps.noMock, Kingdom == "Bacteria")
ps_bbs_bac = prune_taxa(taxa_sums(ps_bbs_bac)>0,ps_bbs_bac); ps_bbs_bac
(1-(ntaxa(ps_bbs_bac)/ntaxa(ps.noMock)))*100

ps_bbs_bac %>%
  tax_glom(taxrank = "Family", NArm=TRUE) %>% #remove NAs
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  plot_bar(fill="Class") +
  facet_wrap(~ BIRTH, scales = "free_x") +
  labs(title = "Bacteria present in baby samples") +
  xlab("Baby samples") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")

#explore taxa
sort(get_taxa_unique(ps_bbs_bac, "Phylum"))  #no NA nor uncharacterized
sort(get_taxa_unique(ps_bbs_bac, "Class"))   #no NA
sort(get_taxa_unique(ps_bbs_bac, "Order"))   #no NA
sort(get_taxa_unique(ps_bbs_bac, "Family"))  #no NA
sort(get_taxa_unique(ps_bbs_bac, "Genus"))   #no NA
sort(get_taxa_unique(ps_bbs_bac, "Species")) #no NA

#number of seq per sample 
summary(taxa_sums(ps_bbs_bac))
summary(sample_sums(ps_bbs_bac))
sd(sample_sums(ps_bbs_bac), na.rm=TRUE)/sqrt(length(sample_sums(ps_bbs_bac)[!is.na(sample_sums(ps_bbs_bac))])) #SE
head(sort(sample_sums(ps_bbs_bac),TRUE))
hist(sample_sums(ps_bbs_bac)) #not much different from ps histogram

prune_taxa(taxa_sums(ps_bbs_bac) > 10, ps_bbs_bac)
prune_taxa(taxa_sums(ps_bbs_bac) > 100, ps_bbs_bac) 
filter_taxa(ps_bbs_bac, function(x) sum(x > 3) > (0.1*length(x)), prune=TRUE) 
filter_taxa(ps_bbs_bac, function (x) {sum(x > 0) > 1}, prune=TRUE)
filter_taxa(ps_bbs_bac, function (x) {sum(x > 0) > 2}, prune=TRUE)
prune_samples(sample_sums(ps_bbs_bac)>=100, ps_bbs_bac)
prune_samples(sample_sums(ps_bbs_bac)>=80, ps_bbs_bac)
prune_samples(sample_sums(ps_bbs_bac)>=50, ps_bbs_bac)

#look for outliers ####
sample_data(ps_bbs_bac)$MONTH = factor(sample_data(ps_bbs_bac)$MONTH, levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18))
sample_data(ps_bbs_bac)$MONTH
# Divide months into three groups (1-6, 7-12, 13-18) and create a new variable
sample_data(ps_bbs_bac)$MONTH_GROUP <- cut(as.numeric(sample_data(ps_bbs_bac)$MONTH), breaks = c(0, 6, 12, 18), labels = c("1-6 Months", "7-12 Months", "13-18 Months"))

#1. nmds ####
nmds = ordinate(ps_bbs_bac, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps_bbs_bac, nmds, color = "MONTH", shape = "BIRTH") + 
  theme_bw() + geom_point(size = 3) + ggtitle("nMDS") +
  geom_text(aes(label = SAMPLE_ID), check_overlap = TRUE, size = 2, nudge_y = -0.05)

#2. alpha diversity ####
alpha.shn = estimate_richness(ps_bbs_bac, split=TRUE, measures="Shannon") 
summary(alpha.shn)
plot_richness(ps_bbs_bac, "INFANT_ID","MONTH", measures = c("Shannon","Simpson")) +
  geom_text(aes(label = SAMPLE_ID), check_overlap = FALSE, size = 3)

adiv <- data.frame(
  "Simpson" = phyloseq::estimate_richness(ps_bbs_bac, measures = "Simpson"),
  "Shannon" = phyloseq::estimate_richness(ps_bbs_bac, measures = "Shannon"))
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
plot_richness(ps_bbs_bac, x = "BIRTH", color = "BIRTH", measures = c("Simpson", "Shannon")) +
  geom_boxplot(outlier.color = NA) +
  geom_text(aes(label = SAMPLE_ID)) +
  #geom_jitter(aes(color = BIRTH), height = 0, width = .2) +
  theme_bw() 
hist(adiv$Shannon)
hist(adiv$Simpson)

shn_bbs_bac = estimate_richness(ps_bbs_bac, split=TRUE, measures="Shannon") 

#samples w/ very high alpha diversity
div.out = rownames(shn_bbs_bac)[(which(shn_bbs_bac>4))]; div.out 
#Remove outliers 
ps_bbs_bac.div01 = prune_samples(!rownames(sample_data(ps_bbs_bac)) %in% div.out, ps_bbs_bac)
ps_bbs_bac.div01 = prune_taxa(taxa_sums(ps_bbs_bac.div01)>0,ps_bbs_bac.div01)
ps_bbs_bac.div01

ps2=ps_bbs_bac.div01

#number of seq per sample ####
summary(sample_sums(ps2))  #or: summary(apply(comm,1,sum))
sd(sample_sums(ps), na.rm=TRUE)/sqrt(length(sample_sums(ps)[!is.na(sample_sums(ps))])) #SE
mean(table(sample_data(ps2)[,3]))

saveRDS(ps2,"~/Documents/xoxo_article/files/metaphlan4/ps_xoxo_metagenome_cleaned.rds")

#samples per infant ####
#birth
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$BIRTH == "VAGINAL"))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$BIRTH == "Csec"))$INFANT_ID)))

#sex
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$SEX == "Female"))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$SEX == "Male"))$INFANT_ID)))

#birth
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$ANIMALS == "pets"))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$ANIMALS == "pets,barn"))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$ANIMALS == "rodents,cockroaches"))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$ANIMALS == "pets,rodents,cockroaches"))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$ANIMALS == "pets,barn,corral,rodents,cockroaches"))$INFANT_ID)))

#weaning time
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$END_BF == 3))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$END_BF == 4))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$END_BF == 5))$INFANT_ID)))

#household size
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$NUM_PEOPLE == 2))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$NUM_PEOPLE == 3))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$NUM_PEOPLE == 4))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$NUM_PEOPLE == 6))$INFANT_ID)))

#age
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$MONTH_GROUP == "1-6 Months"))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$MONTH_GROUP == "7-12 Months"))$INFANT_ID)))
length(levels(as.factor(sample_data(subset_samples(ps2,  sample_data(ps2)$MONTH_GROUP == "13-18 Months"))$INFANT_ID)))

#Rarefaction 
#NO NEED TO RAREFY METAPHLAN RESULTS, since they have already been normalized
#MetaPhlAn estimates the coverage of each marker and computes the clade's coverage as the robust average of the coverage across the markers of the same clade. Finally, the clade's coverages are normalized across all detected clades to obtain the relative abundance of each taxon

#Remove infant XC10 #### 
sample_data(subset_samples(ps2, sample_data(ps2)$INFANT_ID == "XC10"))$MONTH_GROUP
#XC10 has only one sample during the 2nd 6-months and no samples during the 3rd 6-months of sampling and mode of birth and sex are NA
#remove it for any statistical analyses (any analyses other than characterization)
psnoXC10 = subset_samples(ps2, sample_data(ps2)$INFANT_ID != "XC10")
psnoXC10 = prune_taxa(taxa_sums(psnoXC10)>0, psnoXC10); psnoXC10
ntaxa(ps2) - ntaxa(psnoXC10)
nsamples(ps2) - nsamples(psnoXC10)

saveRDS(psnoXC10,"ps_xoxo_metagenome_noXC10.rds")

#Alpha diversity measures ####
ps3=psnoXC10
#shannon
shn.rich = cbind(estimate_richness(ps3, measures = 'shannon'),
                 sample_data(ps3))
ggplot(shn.rich, aes(x = BIRTH, y = Shannon, color=BIRTH)) +  
  geom_boxplot()
ggplot(shn.rich, aes(x = ANIMALS, y = Shannon, color=ANIMALS)) +  
  geom_boxplot()
ggplot(shn.rich, aes(x = MONTH, y = Shannon, color=MONTH)) +  
  geom_boxplot()
ggplot(shn.rich, aes(x = SEX, y = Shannon, color=SEX)) +  
  geom_boxplot()
ggplot(shn.rich, aes(x = MONTH_GROUP, y = Shannon, color=BIRTH)) +  
  geom_boxplot()

#shannon
#check for normality
shap.shn = shapiro.test(shn.rich$Shannon); shap.shn
shap.shn$p.value < 0.01 #normally distributed
shn.rich %>%
  group_by(BIRTH) %>%
  shapiro_test(Shannon)

#LMM ####
library(lme4)
lmm.shn.brt <- lmer(Shannon ~ BIRTH + 
                      (1 | INFANT_ID), data = shn.rich)
summary(lmm.shn.brt)
plot(lmm.shn.brt)

# Check residuals
lmm.shn.brt.resid <- resid(lmm.shn.brt)
lmm.shn.brt.fitted <- fitted(lmm.shn.brt)

#Check the independence of the model residuals with the variable
# Plot residuals vs. fitted
plot(lmm.shn.brt.fitted, lmm.shn.brt.resid)
abline(h = 0, col = "red") #good model (dhomogenous dispersion of the residuals) 
#the pattern of residuals does not depend on the variable

# Histogram of residuals
hist(lmm.shn.brt.resid, breaks = 30, main = "Residuals Histogram")

# QQ plot for normality
qqnorm(lmm.shn.brt.resid) 
qqline(lmm.shn.brt.resid) #residuals are normally distributed

boxplot(lmm.shn.brt.resid ~ MONTH_GROUP, data = shn.rich, xlab = "age", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#To compare the average alpha diversity for each age group (adjusted for infant_ID variation)
library(emmeans)
# Run the pairwise comparisons (Tukey's)
lmm.shn.brt.emm <- emmeans(lmm.shn.brt, pairwise ~ BIRTH)
#p-value is adjusted (Tukey's by default)
plot(lmm.shn.brt.emm)
# View the estimated means and the comparisons
lmm.shn.brt.emm$emmeans
lmm.shn.brt.emm$contrasts

# Get the adjusted means
emm_shn.brt <- as.data.frame(emmeans(lmm.shn.brt, ~ BIRTH)); emm_shn.brt
#extract p-values
library(stringr)
contrasts_shn.brt <- as.data.frame(lmm.shn.brt.emm$contrasts)
contrasts_shn.brt <- contrasts_shn.brt %>%
  mutate(
    contrast_clean = str_remove_all(contrast, "[\\(\\)]"),  # remove parentheses
    group_split = str_split_fixed(contrast_clean, " - ", 2),
    group1 = str_trim(group_split[,1]),
    group2 = str_trim(group_split[,2])
  )
contrasts_shn.brt         

pval.symp.shn.brt <- symnum(contrasts_shn.brt$p.value, corr = FALSE,
                            cutpoints = c(0, .0001, .001,.01,.05, .1, 1),
                            symbols = c("****","***","**","*","ns"," "))

contrasts_shn.brt$pval.symp <- as.character(pval.symp.shn.brt)

# Map factor levels to numeric positions
emm_shn.brt$BIRTH <- factor(emm_shn.brt$BIRTH)
x_positions.shn <- setNames(1:length(levels(emm_shn.brt$BIRTH)), levels(emm_shn.brt$BIRTH))

# Assign numeric xmin and xmax
contrasts_shn.brt$xmin <- x_positions.shn[contrasts_shn.brt$group1]
contrasts_shn.brt$xmax <- x_positions.shn[contrasts_shn.brt$group2]

# Dynamically adjust y.position
max_y.shn <- max(shn.rich$Shannon, na.rm = TRUE)
contrasts_shn.brt$y.position <- max_y.shn + seq(0.1, 0.3, length.out = nrow(contrasts_shn.brt))

# Plot-error bar
ggplot(emm_shn.brt , aes(x = BIRTH, y = emmean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
  labs(
    x = "Birth Mode",
    y = "Adjusted Alpha Diversity",
    title = "Estimated Alpha Diversity by Birth Mode"
  ) +
  theme_minimal()

# Plot-error bar & boxplot
ggplot(emm_shn.brt , aes(x = BIRTH, y = emmean)) +
  geom_violin(data = shn.rich, aes(x = BIRTH, y = Shannon), fill = "gray", alpha = 0.5, width = 0.8) +
  geom_boxplot(data = shn.rich, aes(x = BIRTH, y = Shannon), width = 0.2, outlier.shape = NA, alpha = 0.6) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
  labs(
    x = "Birth Mode",
    y = "Adjusted Alpha Diversity",
    title = "Estimated Alpha Diversity by Birth Mode"
  ) +
  theme_minimal() +
  stat_compare_means()

#FIGURE 5.B ####
#excluded xc10 (ps3)
#shannon

img_shn_birth<- ggplot(emm_shn.brt, aes(x = BIRTH, y = emmean, fill=as.factor(BIRTH)))+
  theme_bw()+
  #geom_violin(alpha=0.9)+
  geom_violin(trim = FALSE, alpha = 0.7,data = shn.rich, aes(x = BIRTH, y = Shannon)) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black",
               data = shn.rich, aes(x = BIRTH, y = Shannon)) +  # Smaller boxplot inside the violin plot
  geom_bracket(data = contrasts_shn.brt,
               bracket.nudge.y = 0.2, tip.length = 0.01,
               aes(xmin = xmin, xmax = xmax, y.position = 4,vjust = 0, #move the text up or down relative to the bracket
                   label = pval.symp), label.size = 12,inherit.aes = FALSE) +
  scale_color_manual(name = "BIRTH", values = c("#225ea8","#df65b0"),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_fill_manual(name = "BIRTH", values = c("#225ea8","#df65b0"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_x_discrete(labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  #scale_y_continuous(limits = c(-0.5, max_y.shn + 2.5)) + # Adjust the limits of the y-axis here
  geom_jitter(data = shn.rich, aes(x = BIRTH, y = Shannon),col = "black", alpha=0.3, size = 0.7) +
  labs(x = "Mode of Birth", y = "Bacterial Shannon Diversity", color = "Mode of Birth") +
  theme(axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=27, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=22, color="black",vjust = 5,face="bold"),
        axis.text.y=element_text(size=22, color="black")); img_shn_birth

# PERMANOVA ####
#make dataframe
df = as(sample_data(ps3), "data.frame")
#bray-curtis distance
dis = phyloseq::distance(ps3,  method = "bray")
#permanova
set.seed(411)
adns = adonis2(dis ~ BIRTH*SEX*ANIMALS*END_BF*NUM_PEOPLE*MONTH_GROUP, df, permutations = 999, strata = df$INFANT_ID) #by = "terms" #distance = bray
adns

# Filter the data to only keep rows where P_value >= 0.05 and order based on p-value and R2
adns_sig <- subset(adns,  `Pr(>F)`<= 0.05) %>% 
  arrange(desc(R2),desc(`Pr(>F)`)); adns_sig
adns_sig$var <- rownames(adns_sig); adns_sig$var
adns_sig$var <- c("Age Category", "Age Category : Proximity to Animals","Proximity to Animals", "Mode of Birth", 
                  "Weaning Time", "Household", "Sex")
adns_sig

# Determine the var# Determine the significance level
significance <- ifelse(adns_sig$`Pr(>F)` <= 0.001, "***", 
                       ifelse(adns_sig$`Pr(>F)` <= 0.01, "**", 
                              ifelse(adns_sig$`Pr(>F)` <= 0.05, "*", 
                                     ifelse(adns_sig$`Pr(>F)` <= 0.1, ".", ""))))

adns_p = ggplot(adns_sig, aes(x = R2, y = reorder(var,R2))) +
  geom_bar(stat = "identity", fill = "#5d8dc9") +
  geom_label(aes(label = paste(paste0(round(R2*100,digits=1), "%;"),paste0(`Pr(>F)`,significance))), hjust = 0.25,  
             color = "#1b2a3c", fill = "#ffd966", size = 5.5, label.padding = unit(0.2, "lines"), label.size = 0.3) + 
  labs(x = "R-squared", y = "") +
  theme_bw() +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "#113a70",linewidth = 2) + 
  theme(                             
    legend.text = element_text(size = 12),      
    legend.title = element_text(size = 14, face="bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y =  element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16, face="bold"))+
  scale_y_discrete(position = "right"); adns_p 

adns_p = ggplot(adns_sig, aes(x = R2, y = reorder(var,R2))) +
  geom_bar(stat = "identity", fill = "#5d8dc9", width = 0.6) + #reduce bar width
  geom_text(aes(label = paste(significance)),  hjust = 0.5,  vjust = 1.5, size = 6, angle = 90, color = "#113a70") + 
  labs(x = "R-squared", y = "") +
  theme_bw() +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "#113a70",linewidth = 2) +
  theme(aspect.ratio = 1/2, #reduce axis scale
        axis.text.x = element_text(size = 14),
        axis.text.y =  element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 18, face="bold"))+
  scale_y_discrete(position = "right"); adns_p 

ggplot(adns_sig, aes(y = R2, x = reorder(var,R2))) +
  geom_bar(stat = "identity", fill = "#cb99cb") +
  geom_label(aes(label = paste(paste0(round(R2*100,digits=1), "%;"),paste0(`Pr(>F)`,significance))), hjust = 0.25, vjust = -0.4, 
             color = "#281e28", fill = "#ffd966", size = 5.5, label.padding = unit(0.2, "lines"), label.size = 0.3) + 
  labs(y = "R-squared", x = "") +
  theme_bw() +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "#7f007f",linewidth = 2) +
  theme(                             
    legend.text = element_text(size = 12),      
    legend.title = element_text(size = 14, face="bold"),
    axis.text.y = element_text(size = 12),
    axis.text.x =  element_text(size = 14, face = "bold", angle = 15),
    axis.title = element_text(size = 16, face="bold")) 

#Ordination ####
sample_data(ps3)$MONTH=factor(sample_data(ps3)$MONTH, levels=c(1,2,3,4,5,6,7,8,9,
                                                               10,11,12,13,14,15,16,17,18))
#nmds
ps_bbs_bac_bray <- ordinate(ps3, "NMDS", "bray")
plot_ordination(ps3, ps_bbs_bac_bray, type="samples", color="BIRTH", shape="BIRTH") + 
  geom_point(size = 3) 

nmds = ordinate(ps3,  method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps3, nmds, type = "samples", color = "MONTH", shape = "BIRTH") + 
  geom_point() + 
  facet_wrap(~MONTH)+
  stat_ellipse(aes(group = BIRTH))

nmds2 = ordinate(ps3, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps3, nmds2, color = "BIRTH", shape = "BIRTH") + 
  theme_bw() + geom_point(size = 4) + ggtitle("nMDS") +
  geom_text(aes(label = SAMPLE_ID), check_overlap = TRUE, size = 2, nudge_y = -0.05)

nmds_bac_birth = plot_ordination(
  physeq = ps3,                                                        
  ordination = nmds2)+       
  #group by birth
  stat_ellipse(aes(group = BIRTH), level = 0.95, linetype = "dashed", alpha = 0.6) + 
  geom_point(aes(fill = BIRTH, shape = BIRTH), size = 3) +    
  scale_fill_manual(name="Mode of Birth", values = c("#225ea8","#df65b0"),
                    labels = c(Csec = "C-Section", VAGINAL = "Vaginal Birth"))+
  scale_shape_manual(name="Mode of Birth", values = c(21, 22),
                     labels = c(Csec = "C-Section", VAGINAL = "Vaginal Birth"))+
  theme_classic() +  
  #stat_ellipse(aes(group = BIRTH))+
  theme(                             
    legend.text = element_text(size = 12),      
    legend.title = element_text(size = 12, face="bold"))+
  theme(axis.text.y.left = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face="bold"))+
  guides(fill = guide_legend(override.aes = list(shape = 21))) +   #fills legend points based on the fill command
  ggtitle(""); nmds_bac_birth

#pcoa
pcoa = ordinate(ps3, method = "MDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps3, pcoa, shape = "BIRTH", color = "MONTH") + geom_point() + ggtitle("PCoA") +
  geom_text(aes(label = MI_ID), check_overlap = TRUE, size = 3, nudge_y = -0.01) 
#scale_color_manual(name = "BIRTH", values = c("chartreuse4", "darkred"))
plot_ordination(ps3, pcoa, shape = "ANIMALS", color = "ANIMALS") + geom_point() + ggtitle("PCoA") +
  geom_text(aes(label = MI_ID), check_overlap = TRUE, size = 3, nudge_y = -0.01) 
plot_ordination(ps3, pcoa, shape = "BIRTH", color = "MONTH") + geom_point() + ggtitle("PCoA") +
  geom_text(aes(label = MI_ID), check_overlap = TRUE, size = 3, nudge_y = -0.01)

pcoa_bac = plot_ordination(
  physeq = ps3,                                                        
  ordination = pcoa)+       
  #group by birth
  #stat_ellipse(aes(group = BIRTH), level = 0.95, linetype = "dashed", alpha = 0.6) +
  geom_point(aes(fill = MONTH, shape = BIRTH), size = 3) +    
  scale_shape_manual(name="Mode of Birth", values = c(21, 22),
                     labels = c(Csec = "C-Section", VAGINAL = "Vaginal Birth"))+
  scale_fill_manual(values = c("#7CE3D8","darkgoldenrod1","mediumvioletred","indianred1","tan3",
                               "cornflowerblue","seagreen","red2","tan4","yellowgreen",
                               "#980043","#dd1c77","#df65b0","#fd8d3c",'#8c6bb1',
                               "#74c476",'#4d004b',"yellow"),
                    name="Age (Month)") +
  theme_classic() +                                                      
  theme(                             
    legend.text = element_text(size = 12),      
    legend.title = element_text(size = 12, face="bold"))+
  #legend.title = element_blank())+ #removes legend title
  #legend.background = element_rect(fill = "white", color = "black")  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face="bold"))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
pcoa_bac

# Computing Bray-Curtis Dissimilarities and PCoA
comm_mat <- vegdist(otu_table(ps3), "bray")
PCoA_comm_mat <- capscale(comm_mat ~ 1, distance = "bray")
PCoA_comm_mat$CA$eig[1:3]/sum(PCoA_comm_mat$CA$eig)
PCoA_scores <- scores(PCoA_comm_mat)$sites

pcoa_scores = scores(pcoa$values$Relative_eig)

#FIGURE 5.A ####
pcoa_bac_birth = plot_ordination(
  physeq = ps3,                                                        
  ordination = pcoa)+       
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = BIRTH, fill=BIRTH), level = 0.95, linetype = "dashed", #alpha = 0.4,
               alpha = 0.3,colour = "black",show.legend = FALSE) + 
  geom_point(aes(fill = BIRTH, shape = BIRTH), size = 3) +    
  scale_fill_manual(name="Mode of Birth", values = c("#225ea8","#df65b0"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal"))+
  scale_shape_manual(name="Mode of Birth", values = c(21, 24),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal"))+
  theme_classic() +  
  # labs(x = paste0("PCoA1 (", round(PCoA_comm_mat$CA$eig[1:1]/sum(PCoA_comm_mat$CA$eig)*100,digits=1), "%)"), 
  #      y = paste0("PCoA2 (", round(PCoA_comm_mat$CA$eig[2:1]/sum(PCoA_comm_mat$CA$eig)*100,digits=1), "%)"))+ 
  labs(
    x = paste0("PCoA1 (", round(pcoa_scores[1:1]/sum(pcoa_scores)*100,digits=1), "%)"), 
    y = paste0("PCoA2 (", round(pcoa_scores[2:2]/sum(pcoa_scores)*100,digits=1), "%)"))+ 
  theme(                             
    legend.text = element_text(size = 25),      
    legend.title = element_text(size = 27, face="bold"),
    legend.position = "right",
    axis.text = element_text(size = 21),
    axis.title = element_text(size = 27, face="bold")) +
  annotate("text", x = 0.05, y = 0.38, label = "paste(italic(R) ^ 2, \" = 0.019; \", italic(p), \" < 0.001***\")", parse = TRUE, 
           size = 11, color = "black"); pcoa_bac_birth

#birth - age
pcoa_bac_birth_age = plot_ordination(
  physeq = ps3,                                                        
  ordination = pcoa)+       
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = BIRTH, fill=BIRTH), level = 0.95, linetype = "dashed", alpha = 0.4) + 
  geom_point(aes(fill = BIRTH, shape = BIRTH), size = 2) +    
  scale_fill_manual(name="Mode of Birth", values = c("#225ea8","#df65b0"),
                    labels = c(Csec = "C-Section", VAGINAL = "Vaginal Birth"))+
  scale_shape_manual(name="Mode of Birth", values = c(21, 22),
                     labels = c(Csec = "C-Section", VAGINAL = "Vaginal Birth"))+
  facet_wrap(~MONTH_GROUP, ncol = 1, strip.position = "right") +
  theme_classic() +  
  labs(x = paste0("PCoA1 (", round(PCoA_comm_mat$CA$eig[1:1]/sum(PCoA_comm_mat$CA$eig)*100,digits=1), "%)"), 
       y = paste0("PCoA2 (", round(PCoA_comm_mat$CA$eig[2:1]/sum(PCoA_comm_mat$CA$eig)*100,digits=1), "%)"))+ 
  theme(                             
    legend.text = element_text(size = 18),      
    legend.title = element_text(size = 16, face="bold"),
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 16, face="bold"),
    strip.text = element_text(size = 14, face="bold")); pcoa_bac_birth_age

#Temporal Shannon & Simpson Diversity ####
# Transformation to achieve normality 
shapiro.test(shn.rich$Shannon)
shn.rich %>%
  group_by(BIRTH) %>%
  shapiro_test(Shannon)
#shn.rich$Shannon_norm = sqrt(shn.rich$Shannon) #OR log(metadata$Shannon)
library(readxl)
#LLM shannon - infant ID is random effect
shn.rich$MONTH <- as.numeric(as.character(shn.rich$MONTH))

lmm.shn.mnt.brt <- lmer(Shannon ~ MONTH + BIRTH + (1 | INFANT_ID), data = shn.rich)
# Show correlation matrix
print(summary(lmm.shn.mnt.brt), correlation = TRUE)
# Or extract variance-covariance matrix
vcov(lmm.shn.mnt.brt)
plot(lmm.shn.mnt.brt) #check residuals & model fit

library(performance)
lmm.shn.mnt.brt.r = r2(lmm.shn.mnt.brt)  # returns marginal and conditional R²
lmm.shn.mnt.brt.r 
#R²m = marginal R² (variance explained by fixed effects)
#R²c = conditional R² (variance explained by fixed + random effects)

lmm.shn.mnt.brt.resid = resid(lmm.shn.mnt.brt)
qqnorm(lmm.shn.mnt.brt.resid); qqline(lmm.shn.mnt.brt.resid)
summary(lmm.shn.mnt.brt.resid)
plotNormalHistogram(lmm.shn.mnt.brt.resid)

shap.lmm.shn.mnt.brt.resid = shapiro.test(lmm.shn.mnt.brt.resid); shap.lmm.shn.mnt.brt.resid #test for normality
shap.lmm.shn.mnt.brt.resid$p.value < 0.05 #normally distributed

#pairwise comparison
library(emmeans)
lmm.shn.mnt.brt.emm = emmeans(lmm.shn.mnt.brt, ~ MONTH + BIRTH)
pairs(lmm.shn.mnt.brt.emm)
lmm.shn.mnt.brt.emm.df = as.data.frame(lmm.shn.mnt.brt.emm)

#overall effect of MONTH
anova(lmm.shn.mnt.brt)

ggplot(shn.rich, aes(x = MONTH, y = Shannon, color = BIRTH)) +
  geom_point(alpha = 0.5) +  # plot raw data points
  geom_smooth(method = "lm", se = TRUE) +  # linear regression line with confidence interval
  theme_minimal() +
  labs(title = "Shannon Diversity over Time", x = "Month", y = "Shannon Diversity")

#FIGURE 5.C ####
img_shn_mnt_birth = 
  ggplot(data = shn.rich %>%
           arrange(MONTH, Shannon) %>%
           dplyr::slice(1:nrow(.)), 
         mapping = aes(x = as.factor(MONTH), y = Shannon, color = BIRTH, fill = BIRTH)) +
  #ggplot(shn.rich, aes(x = MONTH, y = Shannon)) +
  geom_smooth(aes(group = BIRTH),method = "lm", se = TRUE) +
  theme_light() +
  geom_jitter(width = 0.1,alpha=0.3,colour = "black")+
  scale_color_manual(name = "Mode of Birth", values = c("#225ea8","#df65b0"),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_fill_manual(name = "Mode of Birth", values = c("#225ea8","#df65b0"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_x_discrete(breaks = c(1,6,12,18)) +
  annotate("text", x = 1, y = 1.9,
           label = paste0("R²m = ", round(lmm.shn.mnt.brt.r$R2_marginal, 3),
                          "\nR²c = ", round(lmm.shn.mnt.brt.r$R2_conditional, 3),
                          "; CI = ", "95%"), hjust = 0, size = 9) +
  ylim(0,2)+
  geom_jitter(alpha=0.5, size = 1) +
  labs(#title = paste(strwrap("Temporal Bacterial Shannon Diversity"), collapse = "\n"),
    x = "Months", y = "Temporal Changes of Bacterial Shannon Diversity", color = "Mode of Birth") +
  theme(legend.text = element_text(size = 25),      
        legend.title = element_text(size = 27, face="bold"),
        axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=23, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=22, color="black",vjust = 2,face="bold"),
        axis.text.y=element_text(size=22, color="black")); img_shn_mnt_birth #result.shn$pval

# Homogeneity of dispersion test ####
comm = as.matrix(otu_table(ps3))
metadata = data.frame(sample_data(ps3))
#var_per_asv = apply(comm, 1, var)
#distances_to_centroid <- as.integer(apply(comm, 2, var)) 

#birth ####
disper.brt = betadisper(dis, df$BIRTH, type = "centroid"); disper.brt
plot(disper.brt, hull = FALSE, ellipse = TRUE)
# Extract the distances to the centroid for each group
distances_to_centroid_brt <- disper.brt$distance
metadata$distances_to_centroid_brt <- distances_to_centroid_brt

lmm.disp.brt <- lmer(distances_to_centroid_brt ~ BIRTH +
                       (1 | INFANT_ID), data = metadata)

summary(lmm.disp.brt)
# Check residuals
plot(residuals(lmm.disp.brt))
# Check random effects
ranef(lmm.disp.brt)
#normality test
library(rcompanion)
#p-value < 0.05 is not normally distributed
shapiro.test(distances_to_centroid_brt) #not normally distributed
plotNormalHistogram(distances_to_centroid_brt)
shapiro_test(sqrt(distances_to_centroid_brt)) 
plotNormalHistogram(sqrt(distances_to_centroid_brt))
#still not normally distributed: p-value < 0.05
#so we go with GLMM instead of LMM

#GLMM ####
glmm.disp.brt <- glmmTMB(distances_to_centroid_brt ~ BIRTH + (1 | INFANT_ID), 
                         data = metadata, family = Gamma(link = "log")); summary(glmm.disp.brt)
glmm.disp.brt.emm <- emmeans(glmm.disp.brt, pairwise ~ BIRTH); summary(glmm.disp.brt.emm)

# Get the adjusted means
emm_disp.brt <- as.data.frame(emmeans(glmm.disp.brt, ~ BIRTH)); emm_disp.brt
#extract p-values
library(stringr)
contrasts_disp.brt <- as.data.frame(glmm.disp.brt.emm$contrasts)
contrasts_disp.brt <- contrasts_disp.brt %>%
  mutate(
    contrast_clean = str_remove_all(contrast, "[\\(\\)]"),  # remove parentheses
    group_split = str_split_fixed(contrast_clean, " - ", 2),
    group1 = str_trim(group_split[,1]),
    group2 = str_trim(group_split[,2])
  )
contrasts_disp.brt         

pval.symp.disp.brt <- symnum(contrasts_disp.brt$p.value, corr = FALSE,
                             cutpoints = c(0, .0001, .001,.01,.05, .1, 1),
                             symbols = c("****","***","**","*"," ","ns"))

contrasts_disp.brt$pval.symp <- as.character(pval.symp.disp.brt)

# Map factor levels to numeric positions
emm_disp.brt$BIRTH <- factor(emm_disp.brt$BIRTH)
x_positions.disp.brt <- setNames(1:length(levels(emm_disp.brt$BIRTH)), levels(emm_disp.brt$BIRTH))

# Assign numeric xmin and xmax
contrasts_disp.brt$xmin <- x_positions.disp.brt[contrasts_disp.brt$group1]
contrasts_disp.brt$xmax <- x_positions.disp.brt[contrasts_disp.brt$group2]

# Dynamically adjust y.position
max_y.disp.brt <- max(distances_to_centroid_brt, na.rm = TRUE)
contrasts_disp.brt$y.position <- max_y.disp + seq(0.1, 0.3, length.out = nrow(contrasts_disp.brt))

#plot
ggplot(metadata2, aes(x = BIRTH, y = distances_to_centroid_brt)) +
  #geom_violin(trim = FALSE, fill = "gray80") +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  labs(
    x = "Mode of Birth",
    y = "Dispersion",
    title = "Group-wise Dispersion in Bacterial Composition"
  ) +
  theme_minimal()

#plot
# Get the adjusted means
#Compute mean and standard error by group
summary_metadata <- metadata %>%
  group_by(BIRTH) %>%
  summarise(
    mean_disp = mean(distances_to_centroid_brt),
    se_disp = sd(distances_to_centroid_brt) / sqrt(n())
  )

#Plot with error bars
ggplot(summary_metadata, aes(x = BIRTH, y = mean_disp)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_disp - se_disp, ymax = mean_disp + se_disp), width = 0.2) +
  labs(
    x = "Mode of Birth",
    y = "Mean Distance to Centroid",
    title = "Dispersion by Mode of Birth with Error Bars"
  ) +
  theme_minimal()

#FIGURE 5.D ####
p.disper.birth <- ggplot(emm_disp.brt, aes(x = BIRTH, y = emmean, fill=as.factor(BIRTH)))+
  theme_bw()+
  #geom_violin(alpha=0.9)+
  geom_violin(trim = FALSE, alpha = 0.7,data = metadata, aes(x = BIRTH, y = distances_to_centroid_brt)) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black",
               data = metadata, aes(x = BIRTH, y = distances_to_centroid_brt)) +  # Smaller boxplot inside the violin plot
  geom_bracket(data = contrasts_disp.brt,
               bracket.nudge.y = 0.2, tip.length = 0.01,
               aes(xmin = xmin, xmax = xmax, y.position = 0.9,vjust = 0, #move the text up or down relative to the bracket
                   label = pval.symp), label.size = 12,inherit.aes = FALSE,type="text") +
  scale_color_manual(name = "BIRTH", values = c("#225ea8","#df65b0"),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_fill_manual(name = "BIRTH", values = c("#225ea8","#df65b0"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_x_discrete(labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  geom_jitter(data = metadata, aes(x = BIRTH, y = distances_to_centroid_brt),col = "black", alpha=0.3, size = 0.7) +
  xlab("Mode of Birth")+
  ylab("Bacterial Dispersion")+
  theme(axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=25, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=22, color="black",vjust = 5,face="bold"),
        axis.text.y=element_text(size=22, color="black")); p.disper.birth

#save
save.image("~/Documents/xoxo_article/files/metaphlan4/b1_xoxo_metagenome_cleaning_expl.RData")


