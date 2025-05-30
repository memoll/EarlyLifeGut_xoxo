###############################################################
# Functional analysis of metagenomic data                     #
# Data: HUMANN - xoxo                                         # 
# Mona Parizadeh - May 2025                                   #
###############################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") 
library(vegan); packageVersion("vegan") 
library(ggplot2); packageVersion("ggplot2") 
library(tidyverse); packageVersion("tidyverse") 
library(dplyr); packageVersion("dplyr") #1.1.4
library(Maaslin2); packageVersion("Maaslin2") 

setwd('~/Documents/xoxo_article/')

#import data
meta <- read.delim("files/metadata_xoxo_mi_18months_mom_update.csv", header = T, sep = "," ); meta; dim(meta)
rownames(meta) = meta$MI_ID
func = read.delim("06b_humann_split/_unstratified.NORM_TABLE.tsv", header = T, sep = "\t" ); func; dim(func)
rownames(func) = func$X..Pathway
func$X..Pathway <- NULL
colnames(func) = sub("_.*", "",colnames(func))
colnames(func) = sub("\\.", "-",colnames(func))

ps = phyloseq(otu_table(func, taxa_are_rows = TRUE), sample_data(meta))

#Explore data ####
summary(taxa_sums(ps))
summary(sample_sums(ps))

#Check for outliers ####
#% NMDS####
#relative abundance
ps.ra = transform_sample_counts(ps, function(otu) otu/sum(otu)) 
#ordinate
nmds = ordinate(ps, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps.ra, nmds, color = "DATASET", shape = "BIRTH") + 
  theme_bw() + geom_point(size=3) + ggtitle("nMDS") +
  geom_text(aes(label = SAMPLE_ID), check_overlap = FALSE, size = 5) + 
  scale_shape_manual(values = c(19, 23))
# XC05-06

#% alpha diversity ####
shn = estimate_richness(ps, split=TRUE, measures="Shannon") 
plot_richness(ps, "DATASET","BIRTH", measures = "Shannon") +
  geom_text(aes(label = SAMPLE_ID), size = 3)
# XC06-01-R

# 1.Filter functions with less than 5 reads (for low number of reads) #### 
summary(taxa_sums(ps))
ps.rdrare = prune_taxa(taxa_sums(ps) > 5, ps)
ps.rdrare = prune_samples(sample_sums(ps.rdrare)>0,ps.rdrare)
ps.rdrare

#Clean ####
#Remove mock samples 
mock.ID = sample_names(ps.rdrare)[grep("Mock.",sample_data(ps.rdrare)$SAMPLE_ID)]
mock = subset_samples(ps.rdrare, sample_names(ps.rdrare) %in% mock.ID)
mock = prune_taxa(taxa_sums(mock)>0, mock); mock
otu_table(mock)

ps.noMock = subset_samples(ps.rdrare, !sample_names(ps.rdrare) %in% mock.ID)
ps.noMock = prune_taxa(taxa_sums(ps.noMock)>0, ps.noMock); ps.noMock

#remove mom samples
ps.noMom = subset_samples(ps.noMock, sample_data(ps.noMock)$DATASET != "MOMS")
ps.noMom = prune_taxa(taxa_sums(ps.noMom)>0, ps.noMom); ps.noMom

#XC10 has only one sample during the 2nd 6-months and no samples during the 3rd 6-months of sampling and mode of birth and sex are NA
#remove it for any statistical analyses (any analyses other than characterization)
ps.noXC10 = subset_samples(ps.noMom, sample_data(ps.noMom)$INFANT_ID != "XC10")
ps.noXC10 = prune_taxa(taxa_sums(ps.noXC10)>0, ps.noXC10); ps.noXC10

#clean functions
func1 = otu_table(ps.noXC10)
rownames(func1)[grep("UNMAPPED",rownames(func1))]
rownames(func1)[grep("UNINTEGRATED",rownames(func1))]
func2 <- func1[rownames(func1) != "UNMAPPED", ]
func3 <- func2[rownames(func2) != "UNINTEGRATED", ]
ps2 = phyloseq(sample_data(meta), otu_table(func3, taxa_are_rows = TRUE))
ps2

is.numeric(sample_data(ps2)$END_BF)
is.numeric(sample_data(ps2)$NUM_PEOPLE)

#no need to vst (metaphlan4 performs it)

#PERMANOVA ####
#make dataframe
df = as(sample_data(ps2), "data.frame")
df$MONTH_GROUP = factor(df$MONTH_GROUP, levels = c("1-6 Months","7-12 Months","13-18 Months"))

#bray-curtis distance
dis = phyloseq::distance(ps2,  method = "bray")
sample_variables(ps2)
set.seed(112)
adns = adonis2(dis ~ BIRTH*END_BF*SEX*NUM_PEOPLE*ANIMALS*MONTH_GROUP, df, 
               permutations = 999, strata = df$INFANT_ID, by = "terms") #distance = bray
adns

set.seed(222)
adns.age = adonis2(dis ~ MONTH_GROUP, df, permutations = 999, strata = df$INFANT_ID,
                   by = "terms") #distance = bray
adns.age

adonis2(dis ~ MONTH_GROUP, df, permutations = 999, strata = df$INFANT_ID,
        by = "terms") 
adonis2(dis ~ BIRTH, df, permutations = 999, strata = df$INFANT_ID,
        by = "terms") 
adonis2(dis ~ SEX, df, permutations = 999, strata = df$INFANT_ID,
        by = "terms") 
adonis2(dis ~ ANIMALS, df, permutations = 999, strata = df$INFANT_ID,
        by = "terms") 
adonis2(dis ~ END_BF, df, permutations = 999, strata = df$INFANT_ID,
        by = "terms") 
adonis2(dis ~ NUM_PEOPLE, df, permutations = 999, strata = df$INFANT_ID,
        by = "terms") 

# Homogeneity of dispersion test ####
metadata2 = data.frame(sample_data(ps2))
metadata2$BIRTH = factor(metadata2$BIRTH)
metadata2$SEX = factor(metadata2$SEX)
metadata2$ANIMALS = factor(metadata2$ANIMALS)
metadata2$END_BF = factor(metadata2$END_BF)
metadata2$NUM_PEOPLE = factor(metadata2$NUM_PEOPLE)
metadata2$MONTH_GROUP = factor(metadata2$MONTH_GROUP, levels = c("1-6 Months","7-12 Months","13-18 Months"))

#age
disper.age = betadisper(dis, df$MONTH_GROUP, type = "centroid"); disper.age
plot(disper.age, hull = FALSE, ellipse = TRUE)
# Extract the distances to the centroid for each group
distances_to_centroid_age <- disper.age$distance
metadata2$distances_to_centroid_age <- distances_to_centroid_age

#normality test
library(rcompanion)
#p-value < 0.05 is not normally distributed
shapiro.test(distances_to_centroid_age) #normally distributed
plotNormalHistogram(distances_to_centroid_age)
metadata2$distances_to_centroid_age_norm = sqrt(metadata2$distances_to_centroid_age) #OR log(metadata$Shannon)
shapiro.test(metadata2$distances_to_centroid_age_norm)
plotNormalHistogram(metadata2$distances_to_centroid_age_norm)
#still not normally distributed

#LMM####
lmm.disp.age <- lmer(distances_to_centroid_age ~ MONTH_GROUP + 
                       (1 | INFANT_ID), data = metadata2)

summary(lmm.disp.age)

lmm.resid.age <- residuals(lmm.disp.age)
lmm.fitted.age <- fitted(lmm.disp.age)
ggplot(data = data.frame(fitted =lmm.fitted.age, residuals = lmm.resid.age)) +
  geom_point(aes(x = fitted, y = residuals)) +
  geom_hline(yintercept = 0, color = "red") +
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals")

qqnorm(lmm.resid.age); qqline(lmm.resid.age)
summary(lmm.resid.age)
plotNormalHistogram(lmm.resid.age)
#still not normally distributed. GLMM would be a more robust method

#GLMM ####
glmm.disp.age <- glmmTMB(distances_to_centroid_age ~ MONTH_GROUP + (1 | INFANT_ID), 
                         data = metadata2, family = Gamma(link = "log")); summary(glmm.disp.age)

summary(glmm.disp.age)

glmm.resid.age <- residuals(glmm.disp.age)
glmm.fitted.age <- fitted(glmm.disp.age)
ggplot(data = data.frame(fitted =glmm.fitted.age, residuals = glmm.resid.age)) +
  geom_point(aes(x = fitted, y = residuals)) +
  geom_hline(yintercept = 0, color = "red") +
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals")

qqnorm(glmm.resid.age); qqline(glmm.resid.age)
summary(glmm.resid.age)
plotNormalHistogram(glmm.resid.age)

#check for overdispersion
library(performance)
check_overdispersion(glmm.disp.age) #No overdispersion detected.

glmm.disp.age.emm <- emmeans(glmm.disp.age, pairwise ~ MONTH_GROUP); summary(glmm.disp.age.emm)

# Get the adjusted means
emm_disp.age <- as.data.frame(emmeans(glmm.disp.age.emm, ~ MONTH_GROUP)); emm_disp.age
#extract p-values
library(stringr)
contrasts_disp.age <- as.data.frame(glmm.disp.age.emm$contrasts)
contrasts_disp.age <- contrasts_disp.age %>%
  mutate(
    contrast_clean = str_remove_all(contrast, "[\\(\\)]"),  # remove parentheses
    group_split = str_split_fixed(contrast_clean, " - ", 2),
    group1 = str_trim(group_split[,1]),
    group2 = str_trim(group_split[,2])
  )
contrasts_disp.age         

pval.symp.disp.age <- symnum(contrasts_disp.age$p.value, corr = FALSE,
                             cutpoints = c(0, .0001, .001,.01,.05, .1, 1),
                             symbols = c("****","***","**","*"," ","ns"))

contrasts_disp.age$pval.symp <- as.character(pval.symp.disp.age)

# Map factor levels to numeric positions
emm_disp.age$MONTH_GROUP <- factor(emm_disp.age$MONTH_GROUP)
x_positions.disp.age <- setNames(1:length(levels(emm_disp.age$MONTH_GROUP)), levels(emm_disp.age$MONTH_GROUP))

# Assign numeric xmin and xmax
contrasts_disp.age$xmin <- x_positions.disp.age[contrasts_disp.age$group1]
contrasts_disp.age$xmax <- x_positions.disp.age[contrasts_disp.age$group2]

# Dynamically adjust y.position
max_y.disp.age <- max(distances_to_centroid_age, na.rm = TRUE)
contrasts_disp.age$y.position <- max_y.disp.age + seq(0.1, 0.3, length.out = nrow(contrasts_disp.age))

#plot
ggplot(metadata2, aes(x = MONTH_GROUP, y = distances_to_centroid_age)) +
  #geom_violin(trim = FALSE, fill = "gray80") +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  labs(
    x = "Age",
    y = "Dispersion",
    title = "Group-wise Dispersion in Bacterial Functional Composition"
  ) +
  theme_minimal()

#plot
# Get the adjusted means
#Compute mean and standard error by group
summary_metadata2 <- metadata2 %>%
  group_by(MONTH_GROUP) %>%
  summarise(
    mean_disp = mean(distances_to_centroid_age),
    se_disp = sd(distances_to_centroid_age) / sqrt(n())
  )

#Plot with error bars
ggplot(summary_metadata2, aes(x = MONTH_GROUP, y = mean_disp)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_disp - se_disp, ymax = mean_disp + se_disp), width = 0.2) +
  labs(
    x = "Infants' Age",
    y = "Mean Distance to Centroid",
    title = "Dispersion by Infants' Age with Error Bars"
  ) +
  theme_minimal()

# #FIGURE 6.C####
p.disper.age<- ggplot(emm_disp.age, aes(x = MONTH_GROUP, y = emmean, fill=as.factor(MONTH_GROUP)))+
  theme_bw()+
  #geom_violin(alpha=0.9)+
  geom_violin(trim = FALSE, alpha = 0.7,data = metadata2, aes(x = MONTH_GROUP, y = distances_to_centroid_age)) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black",
               data = metadata2, aes(x = MONTH_GROUP, y = distances_to_centroid_age)) +  # Smaller boxplot inside the violin plot
  # annotate("text",  x = 2, y = 0.6, label = "paste(italic(p), \" = ns \")",  
  #          parse = TRUE, size =9, color = "black") +
  geom_bracket(data = contrasts_disp.age,
               bracket.nudge.y = 0.2, tip.length = 0.01,
               aes(xmin = xmin, xmax = xmax, y.position = c(0.2,0.42,0.15),vjust = 0.2, #move the text up or down relative to the bracket
                   label = pval.symp), label.size = 10,inherit.aes = FALSE,type="text") +
  scale_fill_manual(name="Infants' Age (Months)", values = c("#113a70","#5d8dc9","#a0bbde")) +
  scale_x_discrete(name="Infants' Age (Months)", 
                   labels = c("1-6 Months"="1-6", "7-12 Months"="7-12", 
                              "13-18 Months"="13-18"))+
  #scale_y_continuous(limits = c(0, max_y.disp + 2.5)) + # Adjust the limits of the y-axis here
  geom_jitter(data = metadata2, aes(x = MONTH_GROUP, y = distances_to_centroid_age),col = "black", alpha=0.3, size = 0.7) +
  xlab("Infants' Age (Months)")+
  ylab("Bacterial Functional Dispersion")+
  #ggtitle("C") +
  theme(axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=25, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=20, color="black",vjust = 2, face = "bold"),
        axis.text.y=element_text(size=20, color="black")); p.disper.age

#birth
disper.brt = betadisper(dis, df$BIRTH, type = "centroid"); disper.brt
plot(disper.brt, hull = FALSE, ellipse = TRUE)
distances_to_centroid_brt <- disper.brt$distance
metadata2$distances_to_centroid_brt <- distances_to_centroid_brt
glmm.disp.brt <- glmmTMB(distances_to_centroid_brt ~ BIRTH + (1 | INFANT_ID), data = metadata2, family = Gamma(link = "log")); summary(glmm.disp.brt)
glmm.disp.brt.emm <- emmeans(glmm.disp.brt, pairwise ~ BIRTH); summary(glmm.disp.brt.emm)

#sex
disper.sx = betadisper(dis, df$SEX, type = "centroid"); disper.sx
plot(disper.sx, hull = FALSE, ellipse = TRUE)
distances_to_centroid_sx <- disper.sx$distance
metadata2$distances_to_centroid_sx <- distances_to_centroid_sx
glmm.disp.sx <- glmmTMB(distances_to_centroid_sx ~ SEX + (1 | INFANT_ID), data = metadata2, family = Gamma(link = "log")); summary(glmm.disp.sx)
glmm.disp.sx.emm <- emmeans(glmm.disp.sx, pairwise ~ SEX); summary(glmm.disp.sx.emm)

#animals
disper.anm = betadisper(dis, df$ANIMALS, type = "centroid"); disper.anm
plot(disper.anm, hull = FALSE, ellipse = TRUE)
distances_to_centroid_anm <- disper.anm$distance
metadata2$distances_to_centroid_anm <- distances_to_centroid_anm
glmm.disp.anm <- glmmTMB(distances_to_centroid_anm ~ ANIMALS + (1 | INFANT_ID), data = metadata2, family = Gamma(link = "log")); summary(glmm.disp.anm)
glmm.disp.anm.emm <- emmeans(glmm.disp.anm, pairwise ~ ANIMALS); summary(glmm.disp.anm.emm)

#weaning
disper.bf = betadisper(dis, df$END_BF, type = "centroid"); disper.bf
plot(disper.bf, hull = FALSE, ellipse = TRUE)
distances_to_centroid_bf <- disper.bf$distance
glmm.disp.bf <- glmmTMB(distances_to_centroid_bf ~ END_BF + (1 | INFANT_ID), data = metadata2, family = Gamma(link = "log")); summary(glmm.disp.bf)
glmm.disp.bf.emm <- emmeans(glmm.disp.bf, pairwise ~ END_BF); summary(glmm.disp.bf.emm)

#household size
disper.hs = betadisper(dis, df$NUM_PEOPLE, type = "centroid"); disper.hs
plot(disper.hs, hull = FALSE, ellipse = TRUE)
distances_to_centroid_hs <- disper.hs$distance
glmm.disp.hs <- glmmTMB(distances_to_centroid_hs ~ NUM_PEOPLE + (1 | INFANT_ID), data = metadata2, family = Gamma(link = "log")); summary(glmm.disp.hs)
glmm.disp.hs.emm <- emmeans(glmm.disp.hs, pairwise ~ NUM_PEOPLE); summary(glmm.disp.hs.emm)

# #FIGURE 6.G####
#birth ####
disper.brt = betadisper(dis, df$BIRTH, type = "centroid"); disper.brt
plot(disper.brt, hull = FALSE, ellipse = TRUE)
# Extract the distances to the centroid for each group
distances_to_centroid_brt <- disper.brt$distance
metadata2$distances_to_centroid_brt <- distances_to_centroid_brt
glmm.disp.brt <- glmmTMB(distances_to_centroid_brt ~ BIRTH + (1 | INFANT_ID), 
                         data = metadata2, family = Gamma(link = "log")); summary(glmm.disp.brt)
summary(glmm.disp.brt)

glmm.resid.brt <- residuals(glmm.disp.brt)
glmm.fitted.brt <- fitted(glmm.disp.brt)
ggplot(data = data.frame(fitted =glmm.fitted.brt, residuals = glmm.resid.brt)) +
  geom_point(aes(x = fitted, y = residuals)) +
  geom_hline(yintercept = 0, color = "red") +
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals")

qqnorm(glmm.resid.brt); qqline(glmm.resid.brt)
summary(glmm.resid.brt)
plotNormalHistogram(glmm.resid.brt)
#not norma;;y distributed

#check for overdispersion
library(performance)
check_overdispersion(glmm.disp.brt) #No overdispersion detected.

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
contrasts_disp.brt$y.position <- max_y.disp.brt + seq(0.1, 0.3, length.out = nrow(contrasts_disp.fun.brt))

#plot
ggplot(metadata2, aes(x = BIRTH, y = distances_to_centroid_brt)) +
  #geom_violin(trim = FALSE, fill = "gray80") +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  labs(
    x = "Mode of Birth",
    y = "Dispersion",
    title = "Group-wise Dispersion in Bacterial Functional Composition"
  ) +
  theme_minimal()

#plot
# Get the adjusted means
#Compute mean and standard error by group
summary_metadata <- metadata2 %>%
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

p.disper.birth <- ggplot(emm_disp.brt, aes(x = BIRTH, y = emmean, fill=as.factor(BIRTH)))+
  theme_bw()+
  #geom_violin(alpha=0.9)+
  geom_violin(trim = FALSE, alpha = 0.7,data = metadata2, aes(x = BIRTH, y = distances_to_centroid_brt)) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black",
               data = metadata2, aes(x = BIRTH, y = distances_to_centroid_brt)) +  # Smaller boxplot inside the violin plot
  geom_bracket(data = contrasts_disp.brt,
               bracket.nudge.y = 0.2, tip.length = 0.01,
               aes(xmin = xmin, xmax = xmax, y.position = 0.43,vjust = 0, #move the text up or down relative to the bracket
                   label = pval.symp), label.size = 12,inherit.aes = FALSE,type="text") +
  scale_color_manual(name = "BIRTH", values = c("#225ea8","#df65b0"),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_fill_manual(name = "BIRTH", values = c("#225ea8","#df65b0"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_x_discrete(labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  geom_jitter(data = metadata2, aes(x = BIRTH, y = distances_to_centroid_brt),col = "black", alpha=0.3, size = 0.7) +
  xlab("Mode of Birth")+
  ylab("Bacterial Functional Dispersion")+
  theme(axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=25, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=22, color="black",vjust = 5,face="bold"),
        axis.text.y=element_text(size=22, color="black")); p.disper.birth

#Alpha diversity measures ####
#shannon
shn.rich = cbind(estimate_richness(ps2, measures = 'shannon'),
                 sample_data(ps2))

hist(shn.rich$Shannon)

#Age ####
library(rcompanion)
# Check if data is normally distributed (p-value < 0.05 is not normally distributed).
#we only need to check for normality of Shannon diversity if you plan to use statistical models that assume normally distributed residuals (e.g., linear regression, ANOVA, LMM).
shap.shn = shapiro.test(shn.rich$Shannon); shap.shn #test for normality
shap.shn$p.value < 0.05 #normally distributed
shn.rich %>%
  group_by(MONTH_GROUP) %>%
  shapiro_test(Shannon)
plotNormalHistogram(shn.rich$Shannon)

shn.rich %>%
  group_by(MONTH_GROUP) %>%
  summarise(
    mean_sh = mean(Shannon),
    se_shn = sd(Shannon) / sqrt(n())
  )

#LMM ####
library(lme4)
lmm.shn <- lmer(Shannon ~ MONTH_GROUP + BIRTH + SEX + 
                  ANIMALS + END_BF + NUM_PEOPLE +  
                  (1 | INFANT_ID), data = shn.rich)
summary(lmm.shn)
plot(lmm.shn)

lmm.resid <- residuals(lmm.shn)
lmm.fitted <- fitted(lmm.shn)
ggplot(data = data.frame(fitted =lmm.fitted, residuals = lmm.resid)) +
  geom_point(aes(x = fitted, y = residuals)) +
  geom_hline(yintercept = 0, color = "red") +
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals")

#age
shn.rich$MONTH_GROUP = factor(shn.rich$MONTH_GROUP, levels = c("1-6 Months","7-12 Months","13-18 Months"))
lmm.shn.age <- lmer(Shannon ~ MONTH + MONTH_GROUP + (1 | INFANT_ID), data = shn.rich); summary(lmm.shn.age)
library(emmeans)
# Run the pairwise comparisons (Tukey's)
lmm.shn.age.emm <- emmeans(lmm.shn.age, pairwise ~ MONTH_GROUP); summary(lmm.shn.age.emm)
r2(lmm.shn.age)

#birth
shn.rich %>%
  group_by(BIRTH) %>%
  shapiro_test(Shannon)
lmm.shn.brt <- lmer(Shannon ~ BIRTH + (1 | INFANT_ID), data = shn.rich); summary(lmm.shn.brt)
lmm.shn.brt.emm <- emmeans(lmm.shn.brt, pairwise ~ BIRTH); summary(lmm.shn.brt.emm)
r2(lmm.shn.brt)

#sex
lmm.shn.sx <- lmer(Shannon ~ SEX + (1 | INFANT_ID), data = shn.rich); summary(lmm.shn.sx)
lmm.shn.sx.emm <- emmeans(lmm.shn.sx, pairwise ~ SEX); summary(lmm.shn.sx.emm)
r2(lmm.shn.sx)

#animals
lmm.shn.anm <- lmer(Shannon ~ ANIMALS + (1 | INFANT_ID), data = shn.rich); summary(lmm.shn.anm)
lmm.shn.anm.emm <- emmeans(lmm.shn.anm, pairwise ~ ANIMALS); summary(lmm.shn.anm.emm)
r2(lmm.shn.anm)

#weaning
shn.rich$END_BF <- as.factor(shn.rich$END_BF)
lmm.shn.bf <- lmer(Shannon ~ END_BF + (1 | INFANT_ID), data = shn.rich); summary(lmm.shn.bf)
lmm.shn.bf.emm <- emmeans(lmm.shn.bf, pairwise ~ END_BF); summary(lmm.shn.bf.emm)
r2(lmm.shn.bf)

#household size
shn.rich$NUM_PEOPLE <- as.factor(shn.rich$NUM_PEOPLE)
lmm.shn.hs <- lmer(Shannon ~ NUM_PEOPLE + (1 | INFANT_ID), data = shn.rich); summary(lmm.shn.hs)
lmm.shn.hs.emm <- emmeans(lmm.shn.hs, pairwise ~ NUM_PEOPLE); summary(lmm.shn.hs.emm)
r2(lmm.shn.hs)

#age
lmm.shn.age <- lmer(Shannon ~ MONTH_GROUP + (1 | INFANT_ID), data = shn.rich)
summary(lmm.shn.age)
plot(lmm.shn.age)

# Check residuals
lmm.shn.age.resid <- resid(lmm.shn.age)
lmm.shn.age.fitted <- fitted(lmm.shn.age)

#Check the independence of the model residuals with the variable
# Plot residuals vs. fitted
plot(lmm.shn.age.fitted, lmm.shn.age.resid)
abline(h = 0, col = "red") #good model (dhomogenous dispersion of the residuals) 
#the pattern of residuals does not depend on the variable

# Histogram of residuals
hist(lmm.shn.age.resid, breaks = 30, main = "Residuals Histogram")

# QQ plot for normality
qqnorm(lmm.shn.age.resid) 
qqline(lmm.shn.age.resid) #residuals are normally distributed

boxplot(lmm.shn.age.resid ~ MONTH_GROUP, data = shn.rich, xlab = "age", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#To compare the average alpha diversity for each age group (adjusted for infant_ID variation)
library(emmeans)
# Run the pairwise comparisons (Tukey's)
lmm.shn.age.emm <- emmeans(lmm.shn.age, pairwise ~ MONTH_GROUP)
#p-value is adjusted (Tukey's by default)
plot(lmm.shn.age.emm)
# View the estimated means and the comparisons
lmm.shn.age.emm$emmeans
lmm.shn.age.emm$contrasts

# Get the adjusted means
emm_shn.age <- as.data.frame(emmeans(lmm.shn.age, ~ MONTH_GROUP)); emm_shn.age
#extract p-values
library(stringr)
contrasts_shn.age <- as.data.frame(lmm.shn.age.emm$contrasts)
contrasts_shn.age <- contrasts_shn.age %>%
  mutate(
    contrast_clean = str_remove_all(contrast, "[\\(\\)]"),  # remove parentheses
    group_split = str_split_fixed(contrast_clean, " - ", 2),
    group1 = str_trim(group_split[,1]),
    group2 = str_trim(group_split[,2])
  )
contrasts_shn.age         

#y.position = max(emm_shn.age$emmean) + 0.2)  # adjust height

pval.symp.shn.age <- symnum(contrasts_shn.age$p.value, corr = FALSE,
                            cutpoints = c(0, .0001, .001,.01,.05, .1, 1),
                            symbols = c("****","***","**","*"," ","ns"))

contrasts_shn.age$pval.symp <- as.character(pval.symp.shn.age)

# Map factor levels to numeric positions
emm_shn.age$MONTH_GROUP <- factor(emm_shn.age$MONTH_GROUP)
x_positions.shn.age <- setNames(1:length(levels(emm_shn.age$MONTH_GROUP)), levels(emm_shn.age$MONTH_GROUP))

# Assign numeric xmin and xmax
contrasts_shn.age$xmin <- x_positions.shn.age[contrasts_shn.age$group1]
contrasts_shn.age$xmax <- x_positions.shn.age[contrasts_shn.age$group2]

# Dynamically adjust y.position
max_y.shn.age <- max(shn.rich$Shannon, na.rm = TRUE)
contrasts_shn.age$y.position <- max_y.shn.age + seq(0.1, 0.3, length.out = nrow(contrasts_shn.age))

# Plot-error bar
ggplot(emm_shn.age , aes(x = MONTH_GROUP, y = emmean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
  geom_bracket(
    data = contrasts_shn.age,
    aes(xmin = xmin, xmax = xmax,
        y.position = c(3.3,3.8,4),
        label = pval.symp
    ),
    inherit.aes = FALSE
  ) +
  labs(
    x = "Age Group",
    y = "Adjusted Alpha Diversity",
    title = "Estimated Alpha Diversity by Age Group"
  ) +
  theme_minimal()

# Plot-error bar & boxplot
ggplot(emm_shn.age , aes(x = MONTH_GROUP, y = emmean)) +
  geom_violin(data = shn.rich, aes(x = MONTH_GROUP, y = Shannon), fill = "gray", alpha = 0.5, width = 0.8) +
  geom_boxplot(data = shn.rich, aes(x = MONTH_GROUP, y = Shannon), width = 0.2, outlier.shape = NA, alpha = 0.6) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
  labs(
    x = "Age Group",
    y = "Adjusted Alpha Diversity",
    title = "Estimated Alpha Diversity by Age Group"
  ) +
  theme_minimal() +
  stat_compare_means()

#FIGURE 6.A ####
p1_shn_age<- ggplot(emm_shn.age, aes(x = MONTH_GROUP, y = emmean, fill=as.factor(MONTH_GROUP)))+
  theme_bw()+
  #geom_violin(alpha=0.9)+
  geom_violin(trim = FALSE, alpha = 0.7,data = shn.rich, aes(x = MONTH_GROUP, y = Shannon)) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black",
               data = shn.rich, aes(x = MONTH_GROUP, y = Shannon)) +  # Smaller boxplot inside the violin plot
  #scale_x_continuous(breaks = 1:18) +
  geom_bracket(data = contrasts_shn.age,
               bracket.nudge.y = 0.2, tip.length = 0.01,
               aes(xmin = xmin, xmax = xmax, y.position = c(5.7,6,5.5),vjust = 0.1, #move the text up or down relative to the bracket
                   label = pval.symp), label.size = 12,inherit.aes = FALSE) +
  scale_fill_manual(name="Infants' Age (Months)", values = c("#113a70","#5d8dc9","#a0bbde")) +
  scale_x_discrete(name="Infants' Age (Months)", 
                   labels = c("1-6 Months"="1-6", "7-12 Months"="7-12", 
                              "13-18 Months"="13-18"))+
  #scale_y_continuous(limits = c(-0.5, max_y.shn.age + 1)) + # Adjust the limits of the y-axis here
  geom_jitter(data = shn.rich, aes(x = MONTH_GROUP, y = Shannon),col = "black", alpha=0.3, size = 0.7) +
  xlab("Infants' Age (Months)")+
  ylab("Bacterial Functional Shannon Diversity")+
  #ggtitle("C") +
  theme(axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=25, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=20, color="black",vjust = 2, face = "bold"),
        axis.text.y=element_text(size=20, color="black")); p1_shn_age

#Birth ####
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
                            symbols = c("****","***","**","*"," ","ns"))

contrasts_shn.brt$pval.symp <- as.character(pval.symp.shn.brt)

# Map factor levels to numeric positions
emm_shn.brt$BIRTH <- factor(emm_shn.brt$BIRTH)
x_positions.shn.brt <- setNames(1:length(levels(emm_shn.brt$BIRTH)), levels(emm_shn.brt$BIRTH))

# Assign numeric xmin and xmax
contrasts_shn.brt$xmin <- x_positions.shn.brt[contrasts_shn.brt$group1]
contrasts_shn.brt$xmax <- x_positions.shn.brt[contrasts_shn.brt$group2]

# Dynamically adjust y.position
max_y.shn.brt <- max(shn.rich$Shannon, na.rm = TRUE)
contrasts_shn.brt$y.position <- max_y.shn.brt + seq(0.1, 0.3, length.out = nrow(contrasts_shn.brt))

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

#FIGURE 6.E ####
img_shn_birth<- ggplot(emm_shn.brt, aes(x = BIRTH, y = emmean, fill=as.factor(BIRTH)))+
  theme_bw()+
  #geom_violin(alpha=0.9)+
  geom_violin(trim = FALSE, alpha = 0.7,data = shn.rich, aes(x = BIRTH, y = Shannon)) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black",
               data = shn.rich, aes(x = BIRTH, y = Shannon)) +  # Smaller boxplot inside the violin plot
  geom_bracket(data = contrasts_shn.brt,
               bracket.nudge.y = 0.2, tip.length = 0.01,
               aes(xmin = xmin, xmax = xmax, y.position = 5.85,vjust = 0, #move the text up or down relative to the bracket
                   label = pval.symp), label.size = 12,inherit.aes = FALSE) +
  scale_color_manual(name = "BIRTH", values = c("#225ea8","#df65b0"),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_fill_manual(name = "BIRTH", values = c("#225ea8","#df65b0"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_x_discrete(labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  #scale_y_continuous(limits = c(-0.5, max_y.shn.brt + 2.5)) + # Adjust the limits of the y-axis here
  geom_jitter(data = shn.rich, aes(x = BIRTH, y = Shannon),col = "black", alpha=0.3, size = 0.7) +
  labs(x = "Mode of Birth", y = "Bacterial Functional Shannon Diversity", color = "Mode of Birth") +
  theme(axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=27, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=22, color="black",vjust = 5,face="bold"),
        axis.text.y=element_text(size=22, color="black")); img_shn_birth

#Age ####
shn.rich$MONTH <- factor(as.numeric(as.character(shn.rich$MONTH)))
lmm.shn.mnt <- lmer(Shannon ~ MONTH + (1 | INFANT_ID), data = shn.rich)
# Show correlation matrix
print(summary(lmm.shn.mnt), correlation = TRUE)
# Or extract variance-covariance matrix
vcov(lmm.shn.mnt)
plot(lmm.shn.mnt) #check residuals & model fit

# # Get model-predicted values
# lmm.shn.mnt.pred <- predict(lmm.shn.mnt); lmm.shn.mnt.pred

library(performance)
lmm.shn.mnt.r = r2(lmm.shn.mnt)  # returns marginal and conditional R²
lmm.shn.mnt.r
#R²m = marginal R² (variance explained by fixed effects)
#R²c = conditional R² (variance explained by fixed + random effects)

#check for normality
# Check if data is normally distributed (p-value < 0.05 is not normally distributed).
lmm.shn.mnt.resid = resid(lmm.shn.mnt)
qqnorm(lmm.shn.mnt.resid); qqline(lmm.shn.mnt.resid)
summary(lmm.shn.mnt.resid)
plotNormalHistogram(lmm.shn.mnt.resid)

shap.lmm.shn.mnt.resid = shapiro.test(lmm.shn.mnt.resid); shap.lmm.shn.mnt.resid #test for normality
shap.lmm.shn.mnt.resid$p.value < 0.05 #normally distributed

#pairwise comparison
library(emmeans)
lmm.shn.mnt.emm = emmeans(lmm.shn.mnt, ~ MONTH)
pairs(lmm.shn.mnt.emm)
lmm.shn.mnt.emm.df = as.data.frame(lmm.shn.mnt.emm)

#overall effect of MONTH
anova(lmm.shn.mnt)

ggplot(shn.rich, aes(x = MONTH, y = Shannon)) +
  geom_point(alpha = 0.5) +  # plot raw data points
  geom_smooth(method = "lm", se = TRUE) +  # linear regression line with confidence interval
  theme_minimal() +
  labs(title = "Shannon Diversity over Time", x = "Month", y = "Shannon Diversity")

#FIGURE 6.B ####
p1_lm_shn_mnt = ggplot(shn.rich, aes(x = MONTH, y = Shannon)) +
  geom_smooth(aes(group = 1),method = "lm", se = TRUE, fill="#a0bbde",colour = "#5d8dc9") +
  geom_point(data = lmm.shn.mnt.emm.df, aes(x = MONTH, y = emmean), color = "#df65b0", size = 4, inherit.aes = FALSE, alpha = 0.9) +
  geom_errorbar(data = lmm.shn.mnt.emm.df, aes(x = MONTH, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "#df65b0", inherit.aes = FALSE) + ##68477a
  theme_minimal() +
  annotate("text", x = 1.5, y = 5.6,
           label = paste0("R²m = ", round(lmm.shn.mnt.r$R2_marginal, 3),
                          "\nR²c = ", round(lmm.shn.mnt.r$R2_conditional, 3),
                          "; CI = ", "95%"), hjust = 0, size = 10) +
  geom_jitter(width = 0.1,alpha=0.3,colour = "black")+
  scale_x_discrete(breaks = c(1,6,12,18)) +
  xlab("Months")+
  ylab("Temporal Changes of \nBacterial Functional Shannon Diversity")+
  #ggtitle("D") +
  theme(axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=23, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=20, color="black",face="bold"),
        axis.text.y=element_text(size=20, color="black"));p1_lm_shn_mnt

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

#FIGURE 6.F ####
img_shn_mnt_birth = 
  ggplot(data = shn.rich %>%
           arrange(MONTH) %>%
           dplyr::slice(1:nrow(.)), 
         mapping = aes(x = as.factor(MONTH), y = Shannon, color = BIRTH, fill = BIRTH)) +
  geom_smooth(aes(group = BIRTH),method = "lm", se = TRUE) +
  theme_light() +
  geom_jitter(width = 0.1,alpha=0.3,colour = "black")+
  scale_color_manual(name = "Mode of Birth", values = c("#225ea8","#df65b0"),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_fill_manual(name = "Mode of Birth", values = c("#225ea8","#df65b0"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_x_discrete(breaks = c(1,6,12,18)) +
  annotate("text", x = 1, y = 5.6,
           label = paste0("R²m = ", round(lmm.shn.mnt.brt.r$R2_marginal, 3),
                          "\nR²c = ", round(lmm.shn.mnt.brt.r$R2_conditional, 3),
                          "; CI = ", "95%"), hjust = 0, size = 10) +
  geom_jitter(alpha=0.5, size = 1) +
  labs(
    x = "Months", y = "Temporal Changes of \nBacterial Functional Shannon Diversity", color = "Mode of Birth") +
  theme(legend.text = element_text(size = 25),      
        legend.title = element_text(size = 27, face="bold"),
        axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=23, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=22, color="black",vjust = 2,face="bold"),
        axis.text.y=element_text(size=22, color="black")); img_shn_mnt_birth 

# Homogeneity of sample variance (p-value < 0.05 indicates variance is not equal).
library(car); package.version("car") #3.1-3
lev.shn = leveneTest(Shannon ~ BIRTH, data = shn.rich); lev.shn
lev.shn$`Pr(>F)`[1]<0.05 # variances are NOT equal

#Ordination ####
pcoa = ordinate(ps2, method = "MDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps2, pcoa, shape = "BIRTH", color = "MONTH") + geom_point() + ggtitle("PCoA") +
  geom_text(aes(label = MI_ID), check_overlap = TRUE, size = 3, nudge_y = -0.01) 
#scale_color_manual(name = "BIRTH", values = c("chartreuse4", "darkred"))
plot_ordination(ps2, pcoa, shape = "ANIMALS", color = "ANIMALS") + geom_point() + ggtitle("PCoA") +
  geom_text(aes(label = MI_ID), check_overlap = TRUE, size = 3, nudge_y = -0.01) 
plot_ordination(ps2, pcoa, shape = "BIRTH", color = "MONTH") + geom_point() + ggtitle("PCoA") +
  geom_text(aes(label = MI_ID), check_overlap = TRUE, size = 3, nudge_y = -0.01)
# Computing Bray-Curtis Dissimilarities and PCoA
comm_mat <- vegdist(otu_table(ps2), "bray")
PCoA_comm_mat <- capscale(comm_mat ~ 1, distance = "bray")
PCoA_comm_mat$CA$eig[1:3]/sum(PCoA_comm_mat$CA$eig)
PCoA_scores <- scores(PCoA_comm_mat)$sites

pcoa_scores = scores(pcoa$values$Relative_eig)

#FIGURE 6.D####
#age
pcoa_bac_age = plot_ordination(
  physeq = ps2,                                                        
  ordination = pcoa)+       
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = MONTH_GROUP, fill=MONTH_GROUP), level = 0.95, linetype = "dashed", alpha = 0.4) + 
  geom_point(aes(fill = MONTH_GROUP, shape = MONTH_GROUP), size = 3) +    
  scale_color_manual(name="Infants' Age", values = c("#113a70","#5d8dc9","#a0bbde")) +
  scale_fill_manual(name="Infants' Age", values = c("#113a70","#5d8dc9","#a0bbde")) +
  scale_shape_manual(name="Infants' Age", values = c(21, 22, 24))+
  theme_classic() +  
  labs(
    x = paste0("PCoA1 (", round(pcoa_scores[1:1]/sum(pcoa_scores)*100,digits=1), "%)"), 
    y = paste0("PCoA2 (", round(pcoa_scores[2:2]/sum(pcoa_scores)*100,digits=1), "%)"))+ 
  theme(                             
    legend.text = element_text(size = 25),      
    legend.title = element_text(size = 27, face="bold"),
    legend.position = "right",
    axis.text = element_text(size = 21),
    axis.title = element_text(size = 27, face="bold")) +
  annotate("text", x = 0.06, y = 0.21, label = "paste(italic(R) ^ 2, \" = 0.067; \", italic(p), \" < 0.001***\")", #permanova for age category (MONTH_GROUP)
           parse = TRUE, size = 10, color = "black"); pcoa_bac_age

#birth
#FIGURE 6.H####
pcoa_bac_birth = plot_ordination(
  physeq = ps2,                                                        
  ordination = pcoa)+       
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = BIRTH, fill=BIRTH), level = 0.95, linetype = "dashed", alpha = 0.4) + 
  geom_point(aes(fill = BIRTH, shape = BIRTH), size = 3) +    
  scale_fill_manual(name="Mode of Birth", values = c("#225ea8","#df65b0"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal"))+
  scale_shape_manual(name="Mode of Birth", values = c(21, 24),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal"))+
  theme_classic() +  
  labs(
    x = paste0("PCoA1 (", round(pcoa_scores[1:1]/sum(pcoa_scores)*100,digits=1), "%)"), 
    y = paste0("PCoA2 (", round(pcoa_scores[2:2]/sum(pcoa_scores)*100,digits=1), "%)"))+   
  theme(                             
    legend.text = element_text(size = 25),      
    legend.title = element_text(size = 27, face="bold"),
    legend.position = "right",
    axis.text = element_text(size = 21),
    axis.title = element_text(size = 27, face="bold")) +
  annotate("text", x = 0.05, y = 0.21, label = "paste(italic(R) ^ 2, \" = 0.017; \", italic(p), \" < 0.05*\")", #permanova for age category (MONTH_GROUP)
           parse = TRUE, size = 10, color = "black"); pcoa_bac_birth

plot_ordination(
  physeq = ps2,                                                        
  ordination = pcoa)+       
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = SEX, fill=SEX), level = 0.95, linetype = "dashed", alpha = 0.4) + 
  geom_point(aes(fill = SEX, shape = SEX), size = 3)

plot_ordination(
  physeq = ps2,                                                        
  ordination = pcoa)+       
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = ANIMALS, fill=ANIMALS), level = 0.95, linetype = "dashed", alpha = 0.4) + 
  geom_point(aes(fill = ANIMALS, shape = ANIMALS), size = 3)

plot_ordination(
  physeq = ps2,                                                        
  ordination = pcoa)+       
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = END_BF, fill=END_BF), level = 0.95, linetype = "dashed", alpha = 0.4) + 
  geom_point(aes(fill = END_BF, shape = END_BF), size = 3)

plot_ordination(
  physeq = ps2,                                                        
  ordination = pcoa)+       
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = NUM_PEOPLE, fill=NUM_PEOPLE), level = 0.95, linetype = "dashed", alpha = 0.4)  
geom_point(aes(fill = NUM_PEOPLE, shape = NUM_PEOPLE), size = 3)



#FIGURE ALL ####
library(patchwork)
p1_shn_age + p1_lm_shn_mnt + p.disper.age + pcoa_bac_age +
  img_shn_birth + img_shn_mnt_birth + p.disper.birth + pcoa_bac_birth + 
  #  plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() +
  plot_layout(widths = c(1.2,1.7,1.2,1.1), guides = 'collect') +
  plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D', 'E','F','G','H'), '1')) & 
  theme(legend.justification = c("top"),
        plot.tag = element_text(size = 25, face = "bold"))


#TOP 20 Functions ####
#keep xc10
#clean functions
func1.1 = otu_table(ps.noMom)
rownames(func1.1)[grep("UNMAPPED",rownames(func1.1))]
rownames(func1.1)[grep("UNINTEGRATED",rownames(func1.1))]
func2.1 <- func1.1[rownames(func1.1) != "UNMAPPED", ]
func3.1 <- func2.1[rownames(func2.1) != "UNINTEGRATED", ]
ps3 = phyloseq(sample_data(meta), otu_table(func3.1, taxa_are_rows = TRUE))
ps3

# relative abundance
ps3.ra = transform_sample_counts(ps3, function(otu) otu/sum(otu)) 

#save as table
ps3.df=as.data.frame(taxa_sums(ps3.ra))
is.numeric(ps3.df$`taxa_sums(ps3.ra)`)
ps3.df$`taxa_sums(ps3.ra)`= round(ps3.df$`taxa_sums(ps3.ra)`,digits = 4)
write.csv2(ps3.df,file="~/Documents/xoxo_article/files/bac_func_relab.csv",row.names = TRUE)

# Determine 20 most abundant bacterial genera
top20 <- names(sort(taxa_sums(ps3.ra), decreasing = TRUE))[1:20] # Select top 20 genera
ps.top20<- prune_taxa(top20, ps3.ra) # Create new ps object for top 10 genera
ps.top20
pmelt_top20 <- psmelt(ps.top20) # Melt to dataframe

top.df = as.data.frame(taxa_sums(ps.top20))
top.df$PTW = rownames(top.df)
top.df$`taxa_sums(ps.top20)`= round(top.df$`taxa_sums(ps.top20)`,digits = 4)
#sort by relative abundance
top.df.sort <- top.df[order(top.df$`taxa_sums(ps.top20)`, decreasing = TRUE), ]
write.csv2(top.df.sort,file="~/Documents/xoxo_article/files/bac_func_relab_top.csv",row.names = TRUE)

# Determine percent of community top 20 genera account for mode of births (mean and standard deviation)
ps.vg <- subset_samples(ps3, BIRTH == "VAGINAL"); ps.vg
ps.vg.ra <- transform_sample_counts(ps.vg, function(otu) otu/sum(otu)) 
top20.vg <- names(sort(taxa_sums(ps.vg.ra), decreasing = TRUE))[1:20] 
ps.top20.vg<- prune_taxa(top20.vg, ps.vg.ra) # Create new ps object for top 10 genera
ps.top20.vg
pmelt_top20.vg <- psmelt(ps.top20.vg) # Melt to dataframe

ps.cs <- subset_samples(ps3, BIRTH == "Csec"); ps.cs
ps.cs.ra <- transform_sample_counts(ps.cs, function(otu) otu/sum(otu)) 
top20.cs <- names(sort(taxa_sums(ps.cs.ra), decreasing = TRUE))[1:20] 
ps.top20.cs<- prune_taxa(top20.cs, ps.cs.ra) # Create new ps object for top 10 genera
ps.top20.cs
pmelt_top20.cs <- psmelt(ps.top20.cs) # Melt to dataframe

p_top20 <- ggplot(pmelt_top20, aes(x = BIRTH, y = Abundance, fill = OTU, color = OTU)) + 
  geom_bar(stat = "identity", position = "fill") + 
  xlab("") + 
  ylab("Bacterial Relative Abundance") +
  scale_color_viridis_d(name="Top Functions", option = "E")+
  scale_fill_viridis_d(name="Top Functions", option = "E")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c(Csec = "C-Section", VAGINAL = "Vaginal Birth")) +
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(face = "bold", size = 9),
        legend.text = element_text(size = 9, face = "italic"), 
        legend.text.align = 0,
        legend.position = "right"); p_top20 

# Function to retrieve colors for a given order of functions
get_func_col <- function(func_order) {
  return(func_col[func_order])
}
#library(RColorBrewer)
cols <- function(a) image(1:length(a), 1, as.matrix(1:length(a)),
                          col=a, axes=T , xlab="", ylab="")
a <- c("#980043",'#e31a1c',"#df65b0","#006d2c","#31a354",
       "#74c476","#a1d99b","yellow" ,'#edf8b1','#c7e9b4',
       '#7fcdbb','#41b6c4','#1d91c0',"#5d8dc9",'#253494',
       '#081d58','#ffeda0','#fed976','#feb24c','#fd8d3c',
       '#fc4e2a',"#dd1c77",'#bd0026','#88419d','#225ea8',
       '#67000d','#810f7c','#4d004b','#737373','#525252',
       '#252525','#000000')
#cols(a)
col<-a[1:20]
# new_func_order <- as.character(unique(pmelt_top20$OTU))
# col2 <- get_func_col(new_func_order)

#FIGURE #### 
#age
p_top20_age <- ggplot(pmelt_top20, aes(x = MONTH_GROUP, y = Abundance, fill = OTU, colour = OTU)) + 
  geom_bar(stat = "identity", position = "fill") + 
  xlab("\nInfants' Age") + 
  ylab("Relative Abundance") +
  scale_color_manual(name="Top 20 Bacterial Functions", values = col)+ #change colors based on the first plot
  scale_fill_manual(name="Top 20 Bacterial Functions", values = col)+
  scale_x_discrete(labels = c(Csec = "C-Section", VAGINAL = "Vaginal Birth")) +
  theme_bw()+
  theme(axis.title=element_text(face="bold",size=18, color="black"),
        legend.position="right",
        axis.text.x=element_text(size=10, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12)); p_top20_age

# MaAsLin2 ####
# MaAsLin2 using a compound Poisson linear model 
# Run the MaAsLin2 models for normalized and relative abundance
# Outcome: normalized abundance of each function
# Independent variable: age and birth mode

df_input_data <- data.frame(otu_table(ps2)); dim(df_input_data); head(df_input_data)
df_input_metadata <- data.frame(sample_data(ps2)); dim(df_input_metadata); df_input_metadata

#AGE ####
#ref 13-18 Months ####
fit_data_age_ref18m = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/humann/maaslin_xoxo_metagenome_ref18m",
  normalization = "NONE", #[choices: "TSS", "CLR", "CSS", "NONE", "TMM"]
  #data is already normalized, if not for count data: TMM or CSS; for non-count data (CPLM): use TSS or CLR
  #min_prevalence = 0.1, #a 10% prevelance filter 
  #min_abundance = 0.001,
  transform = "LOG", #default [ Choices: LOG, LOGIT, AST (arcsine square-root transformation), NONE ]
  correction = "BH", #default [choices: "BH", "holm", "hochberg", "hommel", "bonferroni", "BY"]
  analysis_method = "LM", # [choices: "LM", "CPLM", "NEGBIN", "ZINB"] 
  #default; for count data: negative binomial distribution; for non-count input: LM or CPLM (compound poisson linear model)
  plot_heatmap = TRUE, #only for variables with more than 2 levels
  fixed_effects = c("MONTH_GROUP"),
  random_effects = c("INFANT_ID"),
  reference =  c("MONTH_GROUP,13-18 Months")) #only for variables with more than 2 levels; NO space between variable and ref category

# Create a new column in our results data.frame to show the stars for the significant findings
fit_data_age_ref18m$results$sig <- ifelse(fit_data_age_ref18m$results$qval >= 0.05, "",
                                          ifelse(fit_data_age_ref18m$results$qval < 0.05 & fit_data_age_ref18m$results$qval >= 0.01, "*",
                                                 ifelse(fit_data_age_ref18m$results$qval < 0.01 & fit_data_age_ref18m$results$qval >= 0.001, "**",
                                                        ifelse(fit_data_age_ref18m$results$qval < 0.001 & fit_data_age_ref18m$results$qval >= 0.0001, "***",
                                                               ifelse(fit_data_age_ref18m$results$qval < 0.0001, "****", "error")))))


#Table 4.C-3####  
#keep only significant results
age_ref18m_res.sig = fit_data_age_ref18m$results[which(fit_data_age_ref18m$results$sig != ""),]; age_ref18m_res.sig; dim(age_ref18m_res.sig)

length(age_ref18m_res.sig$coef[which(age_ref18m_res.sig$coef<0)]) #increased
unique(age_ref18m_res.sig$feature[which(age_ref18m_res.sig$coef<0)])

length(age_ref18m_res.sig$coef[which(age_ref18m_res.sig$coef>0)]) #decreased
unique(age_ref18m_res.sig$feature[which(age_ref18m_res.sig$coef>0)])

age_ref18m_res.sig_top = fit_data_age_ref18m$results %>%
  filter(sig %in% c("***","****")); age_ref18m_res.sig_top; dim(age_ref18m_res.sig_top)

#FIGURE 6.I####
p.age.ref18m_vsig = ggplot(age_ref18m_res.sig_top, aes(x=value, y=feature, fill=coef)) +
  geom_tile(aes(fill=coef),width=1)+
  #geom_dendro(res_forDendro_age) +
  geom_text(aes(label=sig), size=7, nudge_y=-0.3)+
  scale_y_discrete(position = "right") +
  scale_fill_gradient2(mid="white",low="#253494", high="#df65b0", na.value=NA, name="Coefficient\n")+ #,limits=c(-7,7)
  theme_classic() +
  theme(axis.title.x =element_text(face="bold",size=15, color="black"),
        axis.title.y =element_text(face="bold",size=17, color="black"),
        axis.text.x=element_text(colour = "black", size = 16, vjust = 2, face = "bold"),
        #axis.text.y.left =element_text(colour = "black", size = 22),
        axis.text.y.right =element_text(colour = "black", size = 14,
                                        margin = margin(b = 1, unit = "in")),
        legend.title = element_text(colour="black", size=15, face="bold"),
        legend.text = element_text(colour="black", size = 13),
        legend.position = "left")+
  scale_x_discrete(name="Infants' Age (Months)", 
                   labels = c("1-6 Months"="1-6","7-12 Months"="7-12"))+
  xlab("Infants' Age (Months)")+
  ylab("Bacterial Functional Pathways") +
  # scale_x_discrete(labels = c("7-12 Months" = "7-12 Months", "13-18 Months" = "13-18 Months")) +
  coord_cartesian(expand=FALSE) +
  ggtitle("Reference: 1-6 Months") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 17)); p.age.ref18m_vsig

# 18 vs 6 ####
age_ref18m_res.sig_6m = age_ref18m_res.sig[which(age_ref18m_res.sig$value == "1-6 Months"),]; age_ref18m_res.sig_6m; dim(age_ref18m_res.sig_6m)
age_ref18m_res.sig_6m_neg = age_ref18m_res.sig_6m[which(age_ref18m_res.sig_6m$coef <0),]; dim(age_ref18m_res.sig_6m_neg)
length(age_ref18m_res.sig_6m_neg$feature)
age_ref18m_res.sig_6m_pos = age_ref18m_res.sig_6m[which(age_ref18m_res.sig_6m$coef >0),]; dim(age_ref18m_res.sig_6m_pos)
length(age_ref18m_res.sig_6m_pos$feature)
# 18 vs 12 ####
age_ref18m_res.sig_12m = age_ref18m_res.sig[which(age_ref18m_res.sig$value == "7-12 Months"),]; age_ref18m_res.sig_12m; dim(age_ref18m_res.sig_12m)
age_ref18m_res.sig_12m_neg = age_ref18m_res.sig_12m[which(age_ref18m_res.sig_12m$coef <0),]; dim(age_ref18m_res.sig_12m_neg)
length(age_ref18m_res.sig_12m_neg$feature)
age_ref18m_res.sig_12m_pos = age_ref18m_res.sig_12m[which(age_ref18m_res.sig_12m$coef >0),]; dim(age_ref18m_res.sig_12m_pos)
length(age_ref18m_res.sig_12m_pos$feature)

#plot
p.age.ref18m = ggplot(age_ref18m_res.sig, aes(x=value, y=feature, fill=coef))+
  geom_tile()+
  #geom_dendro(res_forDendro_age) +
  geom_text(aes(label=sig), color="black", size=4.5, nudge_y=-0.3)+
  scale_fill_gradient2(mid="white", low="#253494", high="#df65b0", na.value="grey85", name="Coefficient")+
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", angle=40, vjust=.95, hjust=0.95),
        #legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.8, "cm"))+
  xlab(NULL)+
  ylab(NULL) +
  scale_x_discrete(labels = c("7-12 Months" = "7-12 Months", "13-18 Months" = "13-18 Months")) +
  coord_cartesian(expand=FALSE); p.age.ref18m 

#ref 1-6 Months ####
fit_data_age_ref6m = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/humann/maaslin_xoxo_metagenome_ref6m",
  normalization = "NONE", #[choices: "TSS", "CLR", "CSS", "NONE", "TMM"]
  #data is already normalized, if not for count data: TMM or CSS; for non-count data (CPLM): use TSS or CLR
  #min_prevalence = 0.1, #a 10% prevelance filter 
  #min_abundance = 0.001,
  transform = "LOG", #default [ Choices: LOG, LOGIT, AST (arcsine square-root transformation), NONE ]
  correction = "BH", #default [choices: "BH", "holm", "hochberg", "hommel", "bonferroni", "BY"]
  analysis_method = "LM", # [choices: "LM", "CPLM", "NEGBIN", "ZINB"] 
  #default; for count data: negative binomial distribution; for non-count input: LM or CPLM (compound poisson linear model)
  plot_heatmap = TRUE, #only for variables with more than 2 levels
  fixed_effects = c("MONTH_GROUP"),
  random_effects = c("INFANT_ID"),
  reference =  c("MONTH_GROUP,1-6 Months")) #only for variables with more than 2 levels; NO space between variable and ref category

# Create a new column in our results data.frame to show the stars for the significant findings
fit_data_age_ref6m$results$sig <- ifelse(fit_data_age_ref6m$results$qval >= 0.05, "",
                                         ifelse(fit_data_age_ref6m$results$qval < 0.05 & fit_data_age_ref6m$results$qval >= 0.01, "*",
                                                ifelse(fit_data_age_ref6m$results$qval < 0.01 & fit_data_age_ref6m$results$qval >= 0.001, "**",
                                                       ifelse(fit_data_age_ref6m$results$qval < 0.001 & fit_data_age_ref6m$results$qval >= 0.0001, "***",
                                                              ifelse(fit_data_age_ref6m$results$qval < 0.0001, "****", "error")))))
#Table 4.C-1#### 
#keep only significant results
age_ref6m_res.sig = fit_data_age_ref6m$results[which(fit_data_age_ref6m$results$sig != ""),]; age_ref6m_res.sig; dim(age_ref6m_res.sig)

# 6 vs 12 ####
age_ref6m_res.sig_12m = age_ref6m_res.sig[which(age_ref6m_res.sig$value == "7-12 Months"),]; age_ref6m_res.sig_12m; dim(age_ref6m_res.sig_12m)
age_ref6m_res.sig_12m_neg = age_ref6m_res.sig_12m[which(age_ref6m_res.sig_12m$coef <0),]; dim(age_ref6m_res.sig_12m_neg)
age_ref6m_res.sig_12m_pos = age_ref6m_res.sig_12m[which(age_ref6m_res.sig_12m$coef >0),]; dim(age_ref6m_res.sig_12m_pos)

# 6 vs 18 ####
age_ref6m_res.sig_18m = age_ref6m_res.sig[which(age_ref6m_res.sig$value == "13-18 Months"),]; age_ref6m_res.sig_18m; dim(age_ref6m_res.sig_18m)
age_ref6m_res.sig_18m_neg = age_ref6m_res.sig_18m[which(age_ref6m_res.sig_18m$coef <0),]; dim(age_ref6m_res.sig_18m_neg)
age_ref6m_res.sig_18m_pos = age_ref6m_res.sig_18m[which(age_ref6m_res.sig_18m$coef >0),]; dim(age_ref6m_res.sig_18m_pos)

#plot
p.age.ref6m = ggplot(age_ref6m_res.sig, aes(x=value, y=feature, fill=coef))+
  geom_tile()+
  #geom_dendro(res_forDendro_age) +
  geom_text(aes(label=sig), color="black", size=4.5, nudge_y=-0.3)+
  scale_fill_gradient2(mid="white", low="#253494", high="#df65b0", na.value="grey85", name="Coefficient")+
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", angle=40, vjust=.95, hjust=0.95),
        #legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.8, "cm"))+
  xlab(NULL)+
  ylab(NULL) +
  scale_x_discrete(labels = c("7-12 Months" = "7-12 Months", "13-18 Months" = "13-18 Months")) +
  coord_cartesian(expand=FALSE); p.age.ref6m 

#ref 7-12 Months ####
fit_data_age_ref12m = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/humann/maaslin_xoxo_metagenome_ref12m",
  normalization = "NONE", #[choices: "TSS", "CLR", "CSS", "NONE", "TMM"]
  #data is already normalized, if not for count data: TMM or CSS; for non-count data (CPLM): use TSS or CLR
  #min_prevalence = 0.1, #a 10% prevelance filter 
  #min_abundance = 0.001,
  transform = "LOG", #default [ Choices: LOG, LOGIT, AST (arcsine square-root transformation), NONE ]
  correction = "BH", #default [choices: "BH", "holm", "hochberg", "hommel", "bonferroni", "BY"]
  analysis_method = "LM", # [choices: "LM", "CPLM", "NEGBIN", "ZINB"] 
  #default; for count data: negative binomial distribution; for non-count input: LM or CPLM (compound poisson linear model)
  plot_heatmap = TRUE, #only for variables with more than 2 levels
  fixed_effects = c("MONTH_GROUP"),
  random_effects = c("INFANT_ID"),
  reference =  c("MONTH_GROUP,7-12 Months")) #only for variables with more than 2 levels; NO space between variable and ref category

# Create a new column in our results data.frame to show the stars for the significant findings
fit_data_age_ref12m$results$sig <- ifelse(fit_data_age_ref12m$results$qval >= 0.05, "",
                                          ifelse(fit_data_age_ref12m$results$qval < 0.05 & fit_data_age_ref12m$results$qval >= 0.01, "*",
                                                 ifelse(fit_data_age_ref12m$results$qval < 0.01 & fit_data_age_ref12m$results$qval >= 0.001, "**",
                                                        ifelse(fit_data_age_ref12m$results$qval < 0.001 & fit_data_age_ref12m$results$qval >= 0.0001, "***",
                                                               ifelse(fit_data_age_ref12m$results$qval < 0.0001, "****", "error")))))

#Table 4.C-2#### 
#keep only significant results
age_ref12m_res.sig = fit_data_age_ref12m$results[which(fit_data_age_ref12m$results$sig != ""),]; age_ref12m_res.sig; dim(age_ref12m_res.sig)

# 12 vs 18 ####
age_ref12m_res.sig_18m = age_ref12m_res.sig[which(age_ref12m_res.sig$value == "13-18 Months"),]; age_ref12m_res.sig_18m; dim(age_ref12m_res.sig_18m)
age_ref12m_res.sig_18m_pos = age_ref12m_res.sig_18m[which(age_ref12m_res.sig_18m$coef >0),]; dim(age_ref12m_res.sig_18m_pos)
age_ref12m_res.sig_18m_neg = age_ref12m_res.sig_18m[which(age_ref12m_res.sig_18m$coef <0),]; dim(age_ref12m_res.sig_18m_neg)

#plot
age_ref12m_res.sig = fit_data_age_ref12m$results[which(fit_data_age_ref12m$results$sig == "***"),]; age_ref12m_res.sig; dim(age_ref12m_res.sig)
p.age.ref12m_vsig = ggplot(age_ref12m_res.sig, aes(x=value, y=feature, fill=coef)) +
  geom_tile(aes(fill=coef),width=0.6)+
  #geom_point(aes(size = coef),shape = 16) +
  #geom_dendro(res_forDendro_age) +
  geom_text(aes(label=sig), size=7, nudge_y=-0.3)+
  scale_fill_gradient2(mid="white", low="#253494", high="#df65b0", na.value="grey85", name="Coefficient")+
  theme_classic() +
  theme(axis.title=element_text(face="bold",size=16, color="black"),
        axis.text.x=element_text(colour = "black",  size = 12),
        axis.text.y.left =element_text(colour = "black", size = 13, face="italic"),
        axis.text.y.right =element_text(colour = "black", size = 15),
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.width = unit(0.8, "cm"),
        legend.title = element_text(colour="black", size=18, face="bold"),
        legend.text = element_text(colour="black", size = 16),
        legend.position = "right") +
  xlab("Infants' Age")+
  ylab("") +
  scale_x_discrete(labels = c("1-6 Months" = "1-6 Months", "13-18 Months" = "13-18 Months")) +
  #coord_cartesian(expand=TRUE) +
  ggtitle("Reference: 7-12 Months\n") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)); p.age.ref12m_vsig

#BIRTH, with no ref ####
fit_data_brt = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/humann/maaslin_xoxo_metagenome_refBirth",
  normalization = "NONE", #[choices: "TSS", "CLR", "CSS", "NONE", "TMM"]
  #data is already normalized, if not for count data: TMM or CSS; for non-count data (CPLM): use TSS or CLR
  #min_prevalence = 0.1, #a 10% prevelance filter 
  #min_abundance = 0.001,
  transform = "LOG", #default [ Choices: LOG, LOGIT, AST (arcsine square-root transformation), NONE ]
  correction = "BH", #default [choices: "BH", "holm", "hochberg", "hommel", "bonferroni", "BY"]
  analysis_method = "LM", # [choices: "LM", "CPLM", "NEGBIN", "ZINB"] 
  #default; for count data: negative binomial distribution; for non-count input: LM or CPLM (compound poisson linear model)
  plot_heatmap = TRUE, #only for variables with more than 2 levels
  fixed_effects = c("BIRTH"), 
  random_effects = c("INFANT_ID"))

#ref 13-18 Months ####
fit_data_age_ref18m_brt = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/humann/maaslin_xoxo_metagenome_ref18mBirth",
  normalization = "NONE", #[choices: "TSS", "CLR", "CSS", "NONE", "TMM"]
  #data is already normalized, if not for count data: TMM or CSS; for non-count data (CPLM): use TSS or CLR
  #min_prevalence = 0.1, #a 10% prevelance filter 
  #min_abundance = 0.001,
  transform = "LOG", #default [ Choices: LOG, LOGIT, AST (arcsine square-root transformation), NONE ]
  correction = "BH", #default [choices: "BH", "holm", "hochberg", "hommel", "bonferroni", "BY"]
  analysis_method = "LM", # [choices: "LM", "CPLM", "NEGBIN", "ZINB"] 
  #default; for count data: negative binomial distribution; for non-count input: LM or CPLM (compound poisson linear model)
  plot_heatmap = TRUE, #only for variables with more than 2 levels
  fixed_effects = c("MONTH_GROUP","BIRTH"), 
  random_effects = c("INFANT_ID"),
  reference =  c("MONTH_GROUP,13-18 Months")) #only for variables with more than 2 levels; NO space between variable and ref category

# Create a new column in our results data.frame to show the stars for the significant findings
fit_data_age_ref18m_brt$results$sig <- ifelse(fit_data_age_ref18m_brt$results$qval >= 0.1, "",
                                              ifelse(fit_data_age_ref18m_brt$results$qval < 0.1 & fit_data_age_ref18m_brt$results$qval >= 0.05, "*",
                                                     ifelse(fit_data_age_ref18m_brt$results$qval < 0.05 & fit_data_age_ref18m_brt$results$qval >= 0.01, "**",
                                                            ifelse(fit_data_age_ref18m_brt$results$qval < 0.01 & fit_data_age_ref18m_brt$results$qval >= 0.001, "***",
                                                                   ifelse(fit_data_age_ref18m_brt$results$qval < 0.001, "****", "error")))))

#keep only significant results
age_ref18m_brt_res.sig = fit_data_age_ref18m_brt$results[which(fit_data_age_ref18m_brt$results$sig != ""),]; age_ref18m_brt_res.sig; dim(age_ref18m_brt_res.sig)

# Csec vs vaginal ####
age_ref18m_brt_res.sig_vag = age_ref18m_brt_res.sig[which(age_ref18m_brt_res.sig$value == "VAGINAL"),]; age_ref18m_brt_res.sig_vag; dim(age_ref18m_brt_res.sig_vag)
age_ref18m_brt_res.sig_vag_neg = age_ref18m_brt_res.sig_vag[which(age_ref18m_brt_res.sig_vag$coef <0),]; dim(age_ref18m_brt_res.sig_vag_neg)
length(age_ref18m_brt_res.sig_vag_neg$feature)

age_ref18m_brt_res.sig_vag_pos = age_ref18m_brt_res.sig_vag[which(age_ref18m_brt_res.sig_vag$coef >0),]; dim(age_ref18m_brt_res.sig_vag_pos)
length(age_ref18m_brt_res.sig_vag_pos$feature)

#ref 1-6 Months ####
fit_data_age_ref6m_brt = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/humann/maaslin_xoxo_metagenome_ref6mBirth",
  normalization = "NONE", #[choices: "TSS", "CLR", "CSS", "NONE", "TMM"]
  #data is already normalized, if not for count data: TMM or CSS; for non-count data (CPLM): use TSS or CLR
  #min_prevalence = 0.1, #a 10% prevelance filter 
  #min_abundance = 0.001,
  transform = "LOG", #default [ Choices: LOG, LOGIT, AST (arcsine square-root transformation), NONE ]
  correction = "BH", #default [choices: "BH", "holm", "hochberg", "hommel", "bonferroni", "BY"]
  analysis_method = "LM", # [choices: "LM", "CPLM", "NEGBIN", "ZINB"] 
  #default; for count data: negative binomial distribution; for non-count input: LM or CPLM (compound poisson linear model)
  plot_heatmap = TRUE, #only for variables with more than 2 levels
  fixed_effects = c("MONTH_GROUP","BIRTH"),
  random_effects = c("INFANT_ID"),
  reference =  c("MONTH_GROUP,1-6 Months")) #only for variables with more than 2 levels; NO space between variable and ref category

# Create a new column in our results data.frame to show the stars for the significant findings
fit_data_age_ref6m_brt$results$sig <- ifelse(fit_data_age_ref6m_brt$results$qval >= 0.1, "",
                                             ifelse(fit_data_age_ref6m_brt$results$qval < 0.1 & fit_data_age_ref6m_brt$results$qval >= 0.05, "*",
                                                    ifelse(fit_data_age_ref6m_brt$results$qval < 0.05 & fit_data_age_ref6m_brt$results$qval >= 0.01, "**",
                                                           ifelse(fit_data_age_ref6m_brt$results$qval < 0.01 & fit_data_age_ref6m_brt$results$qval >= 0.001, "***",
                                                                  ifelse(fit_data_age_ref6m_brt$results$qval < 0.001, "****", "error")))))

#keep only significant results
age_ref6m_brt_res.sig = fit_data_age_ref6m_brt$results[which(fit_data_age_ref6m_brt$results$sig != ""),]; age_ref6m_brt_res.sig; dim(age_ref6m_brt_res.sig)

# Csec vs vaginal ####
age_ref6m_brt_res.sig_vag = age_ref6m_brt_res.sig[which(age_ref6m_brt_res.sig$value == "VAGINAL"),]; age_ref6m_brt_res.sig_vag; dim(age_ref6m_brt_res.sig_vag)

age_ref6m_brt_res.sig_vag_neg = age_ref6m_brt_res.sig_vag[which(age_ref6m_brt_res.sig_vag$coef <0),]; dim(age_ref6m_brt_res.sig_vag_neg)

age_ref6m_brt_res.sig_vag_pos = age_ref6m_brt_res.sig_vag[which(age_ref6m_brt_res.sig_vag$coef >0),]; dim(age_ref6m_brt_res.sig_vag_pos)
length(age_ref6m_brt_res.sig_vag_pos$feature)

#ref 7-12 Months ####
fit_data_age_ref12m_brt = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/humann/maaslin_xoxo_metagenome_ref12mBirth",
  normalization = "NONE", #[choices: "TSS", "CLR", "CSS", "NONE", "TMM"]
  #data is already normalized, if not for count data: TMM or CSS; for non-count data (CPLM): use TSS or CLR
  #min_prevalence = 0.1, #a 10% prevelance filter 
  #min_abundance = 0.001,
  transform = "LOG", #default [ Choices: LOG, LOGIT, AST (arcsine square-root transformation), NONE ]
  correction = "BH", #default [choices: "BH", "holm", "hochberg", "hommel", "bonferroni", "BY"]
  analysis_method = "LM", # [choices: "LM", "CPLM", "NEGBIN", "ZINB"] 
  #default; for count data: negative binomial distribution; for non-count input: LM or CPLM (compound poisson linear model)
  plot_heatmap = TRUE, #only for variables with more than 2 levels
  fixed_effects = c("MONTH_GROUP","BIRTH"),
  random_effects = c("INFANT_ID"),
  reference =  c("MONTH_GROUP,7-12 Months")) #only for variables with more than 2 levels; NO space between variable and ref category

# Create a new column in our results data.frame to show the stars for the significant findings
fit_data_age_ref12m_brt$results$sig <- ifelse(fit_data_age_ref12m_brt$results$qval >= 0.1, "",
                                              ifelse(fit_data_age_ref12m_brt$results$qval < 0.1 & fit_data_age_ref12m_brt$results$qval >= 0.05, "*",
                                                     ifelse(fit_data_age_ref12m_brt$results$qval < 0.05 & fit_data_age_ref12m_brt$results$qval >= 0.01, "**",
                                                            ifelse(fit_data_age_ref12m_brt$results$qval < 0.01 & fit_data_age_ref12m_brt$results$qval >= 0.001, "***",
                                                                   ifelse(fit_data_age_ref12m_brt$results$qval < 0.001, "****", "error")))))

#keep only significant results
age_ref12m_brt_res.sig = fit_data_age_ref12m_brt$results[which(fit_data_age_ref12m_brt$results$sig != ""),]; age_ref12m_brt_res.sig; dim(age_ref12m_brt_res.sig)

# Csec vs vaginal ####
age_ref12m_brt_res.sig_vag = age_ref12m_brt_res.sig[which(age_ref12m_brt_res.sig$value == "VAGINAL"),]; age_ref12m_brt_res.sig_vag; dim(age_ref12m_brt_res.sig_vag)

#save
save.image("~/Documents/xoxo_article/files/b3_xoxo_function.RData")
