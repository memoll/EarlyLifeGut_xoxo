###############################################################
# Differential analysis of metagenomic data - MaAsLin2        #
# Data: Metaphlan4 - xoxo                                     # 
# Mona Parizadeh - May 2025                                   #
###############################################################

# MaAsLin2 using a compound Poisson linear model 
# Run the MaAsLin2 models for normalized and relative abundance
# Outcome: normalized abundance of each ASV/species
# Independent variable: age

#load libraries ####
library(phyloseq); packageVersion("phyloseq") 
library(vegan); packageVersion("vegan") 
library(ggplot2); packageVersion("ggplot2") 
library(dplyr); packageVersion("dplyr") 
library(Maaslin2); packageVersion("Maaslin2") 

#import data ####
ps = readRDS("~/Documents/xoxo_article/files/metaphlan4/ps_xoxo_metagenome_noXC10.rds");ps 
length(unique(tax_table(ps)[,6])) 

df_input_data <- data.frame(otu_table(ps)); dim(df_input_data); head(df_input_data)
df_input_metadata <- data.frame(sample_data(ps)); dim(df_input_metadata); df_input_metadata
df_taxa <- data.frame(tax_table(ps)); dim(df_taxa); df_taxa
#add otus from rownames as a new column called feature
df_taxa$feature <- rownames(df_taxa)

#AGE ####
#ref 13-18 Months ####
fit_data_age_ref18m = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/metaphlan4/maaslin_xoxo_metagenome_ref18m",
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

# add taxa to fit data
fit_data_age_ref18m$results <- fit_data_age_ref18m$results %>%
  dplyr::left_join(df_taxa, by = "feature"); dim(fit_data_age_ref18m$results)
head(fit_data_age_ref18m$results)

#rename phyla names
fit_data_age_ref18m$results$Phylum = ifelse(fit_data_age_ref18m$results$Phylum == "Actinobacteria" , "Actinomycetota",
                                            ifelse(fit_data_age_ref18m$results$Phylum == "Bacteroidetes" , "Bacteroidota", 
                                                   ifelse(fit_data_age_ref18m$results$Phylum == "Firmicutes" , "Bacillota",
                                                          ifelse(fit_data_age_ref18m$results$Phylum == "Proteobacteria" , "Pseudomonadota",
                                                                 ifelse(fit_data_age_ref18m$results$Phylum == "Verrucomicrobia" , "Verrucomicrobiota",
                                                                        fit_data_age_ref18m$results$Phylum)))))
#keep only significant results
age_ref18m_res.sig = fit_data_age_ref18m$results[which(fit_data_age_ref18m$results$sig != ""),]; age_ref18m_res.sig; dim(age_ref18m_res.sig)
unique(age_ref18m_res.sig$Genus); length(unique(age_ref18m_res.sig$Genus))
unique(age_ref18m_res.sig$Species)
age_ref18m_res.sig %>%
  group_by(Phylum) %>%
  summarize(n = n()) %>%
  arrange(desc(n))
#count(Phylum) %>%
#arrange(desc(n))

age_ref18m_res.sig = fit_data_age_ref18m$results[which(fit_data_age_ref18m$results$sig != ""),]; age_ref18m_res.sig; dim(age_ref18m_res.sig)

age_ref18m_res.sig_top = fit_data_age_ref18m$results %>%
  filter(sig %in% c("***","****")); age_ref18m_res.sig_top; dim(age_ref18m_res.sig_top)

unique(age_ref18m_res.sig$Species)
unique(age_ref18m_res.sig_top$Species)
summary(age_ref18m_res.sig_top$coef[which(age_ref18m_res.sig_top$coef<0)])
summary(age_ref18m_res.sig_top$coef[which(age_ref18m_res.sig_top$coef>0)])

# 18 vs 6 ####
age_ref18m_res.sig_6m = age_ref18m_res.sig[which(age_ref18m_res.sig$value == "1-6 Months"),]; age_ref18m_res.sig_6m; dim(age_ref18m_res.sig_6m)
unique(age_ref18m_res.sig_6m$Genus); length(unique(age_ref18m_res.sig_6m$Genus))

age_ref18m_res.sig_6m_neg = age_ref18m_res.sig_6m[which(age_ref18m_res.sig_6m$coef <0),]; dim(age_ref18m_res.sig_6m_neg)
length(age_ref18m_res.sig_6m_neg$feature)
unique(age_ref18m_res.sig_6m_neg$Genus); length(unique(age_ref18m_res.sig_6m_neg$Genus))
length(unique(age_ref18m_res.sig_6m_neg$Species))

age_ref18m_res.sig_6m_pos = age_ref18m_res.sig_6m[which(age_ref18m_res.sig_6m$coef >0),]; dim(age_ref18m_res.sig_6m_pos)
length(age_ref18m_res.sig_6m_pos$feature)
unique(age_ref18m_res.sig_6m_pos$Genus); length(unique(age_ref18m_res.sig_6m_pos$Genus))

# 18 vs 12 ####
age_ref18m_res.sig_12m = age_ref18m_res.sig[which(age_ref18m_res.sig$value == "7-12 Months"),]; age_ref18m_res.sig_12m; dim(age_ref18m_res.sig_12m)
unique(age_ref18m_res.sig_12m$Genus); length(unique(age_ref18m_res.sig_12m$Genus))

age_ref18m_res.sig_12m_neg = age_ref18m_res.sig_12m[which(age_ref18m_res.sig_12m$coef <0),]; dim(age_ref18m_res.sig_12m_neg)
length(age_ref18m_res.sig_12m_neg$feature)
unique(age_ref18m_res.sig_12m_neg$Genus); length(unique(age_ref18m_res.sig_12m_neg$Genus))

age_ref18m_res.sig_12m_pos = age_ref18m_res.sig_12m[which(age_ref18m_res.sig_12m$coef >0),]; dim(age_ref18m_res.sig_12m_pos)
length(age_ref18m_res.sig_12m_pos$feature)
unique(age_ref18m_res.sig_12m_pos$Genus); length(unique(age_ref18m_res.sig_12m_pos$Genus))

#plot-OTU
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

#plot-Species
p.age.ref18m_spe = ggplot(age_ref18m_res.sig, aes(x=value, y=Species, fill=coef)) +
  geom_tile()+
  #geom_dendro(res_forDendro_age) +
  facet_grid(Phylum ~ ., scales = "free", space = "free") + #switch = "y"
  geom_text(aes(label=sig), size=7, nudge_y=-0.3)+
  scale_fill_gradient2(mid="white", low="#5d8dc9", high="#df65b0", na.value="grey85", name="Coefficient",limits=c(-7,7))+
  theme_classic() +
  theme(axis.title=element_text(face="bold",size=18, color="black"),
        axis.text.x=element_text(colour = "black", face = "bold", size = 10),
        axis.text.y=element_text(colour = "black", size = 15, face = "italic"),
        legend.title = element_text(colour="black", size=18, face="bold"),
        legend.text = element_text(colour="black", size = 16),
        legend.position = "bottom",
        strip.text.y = element_text(size = 16, face = "bold.italic", margin = margin(0.2,0.2,0.2,0.2, "cm"), angle = 0),
        strip.background=element_rect(colour = "#5d8dc9",fill="#5d8dc9")) +
  xlab("Infants' Age")+
  ylab("Bacterial Species\n") +
  coord_cartesian(expand=FALSE) +
  ggtitle("Reference: 13-18 Months") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 25)); p.age.ref18m_spe

#ref 1-6 Months ####
fit_data_age_ref6m = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/metaphlan4/maaslin_xoxo_metagenome_ref6m",
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

# add taxa to fit data
fit_data_age_ref6m$results <- fit_data_age_ref6m$results %>%
  dplyr::left_join(df_taxa, by = "feature"); dim(fit_data_age_ref6m$results)
head(fit_data_age_ref6m$results)

#rename phyla names
fit_data_age_ref6m$results$Phylum = ifelse(fit_data_age_ref6m$results$Phylum == "Actinobacteria" , "Actinomycetota",
                                           ifelse(fit_data_age_ref6m$results$Phylum == "Bacteroidetes" , "Bacteroidota", 
                                                  ifelse(fit_data_age_ref6m$results$Phylum == "Firmicutes" , "Bacillota",
                                                         ifelse(fit_data_age_ref6m$results$Phylum == "Proteobacteria" , "Pseudomonadota",
                                                                ifelse(fit_data_age_ref6m$results$Phylum == "Verrucomicrobia" , "Verrucomicrobiota",
                                                                       fit_data_age_ref6m$results$Phylum)))))

#keep only significant results
age_ref6m_res.sig = fit_data_age_ref6m$results[which(fit_data_age_ref6m$results$sig != ""),]; age_ref6m_res.sig; dim(age_ref6m_res.sig)
unique(age_ref6m_res.sig$Genus); length(unique(age_ref6m_res.sig$Genus))

# 6 vs 12 ####
age_ref6m_res.sig_12m = age_ref6m_res.sig[which(age_ref6m_res.sig$value == "7-12 Months"),]; age_ref6m_res.sig_12m; dim(age_ref6m_res.sig_12m)
unique(age_ref6m_res.sig_12m$Genus); length(unique(age_ref6m_res.sig_12m$Genus))
# age_ref6m_res.sig_12m %>%
#   count(Phylum) %>%
#   arrange(desc(n))

age_ref6m_res.sig_12m_neg = age_ref6m_res.sig_12m[which(age_ref6m_res.sig_12m$coef <0),]; dim(age_ref6m_res.sig_12m_neg)
length(age_ref6m_res.sig_12m_neg$feature)
unique(age_ref6m_res.sig_12m_neg$Genus); length(unique(age_ref6m_res.sig_12m_neg$Genus))

age_ref6m_res.sig_12m_pos = age_ref6m_res.sig_12m[which(age_ref6m_res.sig_12m$coef >0),]; dim(age_ref6m_res.sig_12m_pos)
length(age_ref6m_res.sig_12m_pos$feature)
unique(age_ref6m_res.sig_12m_pos$Genus); length(unique(age_ref6m_res.sig_12m_pos$Genus))
unique(age_ref6m_res.sig_12m_pos$Species); length(unique(age_ref6m_res.sig_12m_pos$Species))
View(age_ref6m_res.sig_12m_pos %>%
       group_by(Phylum, Species) %>%                    
       summarize(count = n(), .groups = 'drop') %>%     
       arrange(Phylum, desc(count)) )

# 6 vs 18 ####
age_ref6m_res.sig_18m = age_ref6m_res.sig[which(age_ref6m_res.sig$value == "13-18 Months"),]; age_ref6m_res.sig_18m; dim(age_ref6m_res.sig_18m)
unique(age_ref6m_res.sig_18m$Genus); length(unique(age_ref6m_res.sig_18m$Genus))

age_ref6m_res.sig_18m_neg = age_ref6m_res.sig_18m[which(age_ref6m_res.sig_18m$coef <0),]; dim(age_ref6m_res.sig_18m_neg)
length(age_ref6m_res.sig_18m_neg$feature)
unique(age_ref6m_res.sig_18m_neg$Genus); length(unique(age_ref6m_res.sig_18m_neg$Genus))

age_ref6m_res.sig_18m_pos = age_ref6m_res.sig_18m[which(age_ref6m_res.sig_18m$coef >0),]; dim(age_ref6m_res.sig_18m_pos)
length(age_ref6m_res.sig_18m_pos$feature)
unique(age_ref6m_res.sig_18m_pos$Genus); length(unique(age_ref6m_res.sig_18m_pos$Genus))

#plot-OTU
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

#FIGURE 1.A####
#plot-Species
#age_ref6m_res.sig = fit_data_age_ref6m$results[which(c(fit_data_age_ref6m$results$sig == "****"),fit_data_age_ref6m$results$sig == "***"),]; age_ref6m_res.sig; dim(age_ref6m_res.sig)
age_ref6m_res.sig = fit_data_age_ref6m$results %>%
  filter(sig %in% c("***","****")); age_ref6m_res.sig
unique(age_ref6m_res.sig$Genus); length(unique(age_ref6m_res.sig$Genus))
age_ref6m_res.sig$Species = sub("_"," ", age_ref6m_res.sig$Species)
#plot-Species
p.age.ref6m_spe = ggplot(age_ref6m_res.sig, aes(x=value, y=Species, fill=coef)) +
  geom_tile(aes(fill=coef),width=0.6)+
  #geom_dendro(res_forDendro_age) +
  facet_grid(Phylum ~ ., scales = "free", space = "free") + #switch = "y"
  geom_text(aes(label=sig), size=7, nudge_y=-0.2)+
  scale_fill_gradient2(mid="white", low="#5d8dc9", high="#df65b0", na.value="grey85", name="Coefficient")+
  theme_classic() +
  theme(axis.title=element_text(face="bold",size=16, color="black"),
        axis.text.x=element_text(colour = "black",  size = 12, vjust = 2, face = "bold"),
        axis.text.y.left =element_text(colour = "black", size = 13, face="italic"),
        axis.text.y.right =element_text(colour = "black", size = 14),
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.width = unit(0.8, "cm"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 14),
        legend.position = "bottom",
        #strip.clip = "off",
        strip.text.y = element_text(size = 14,face = "bold.italic",margin = margin(0.2,0.2,0.2,0.2, "cm"), angle = 0), 
        strip.background=element_rect(colour = "#5d8dc9",fill="#5d8dc9")) +
  scale_x_discrete(name="Infants' Age (Months)", 
                   labels = c("7-12 Months"="7-12","13-18 Months"="13-18"))+
  xlab("Infants' Age (Months)")+
  ylab("Bacterial Species") +
  # scale_x_discrete(labels = c("7-12 Months" = "7-12 Months", "13-18 Months" = "13-18 Months")) +
  coord_cartesian(expand=FALSE) +
  ggtitle("Reference: 1-6 Months") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)); p.age.ref6m_spe

#for adjustment: https://forum.posit.co/t/geom-tile-adjustments/35128

#ref 7-12 Months ####
fit_data_age_ref12m = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/metaphlan4/maaslin_xoxo_metagenome_ref12m",
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

# add taxa to fit data
fit_data_age_ref12m$results <- fit_data_age_ref12m$results %>%
  dplyr::left_join(df_taxa, by = "feature"); dim(fit_data_age_ref12m$results)
head(fit_data_age_ref12m$results)

#rename phyla names
fit_data_age_ref12m$results$Phylum = ifelse(fit_data_age_ref12m$results$Phylum == "Actinobacteria" , "Actinomycetota",
                                            ifelse(fit_data_age_ref12m$results$Phylum == "Bacteroidetes" , "Bacteroidota", 
                                                   ifelse(fit_data_age_ref12m$results$Phylum == "Firmicutes" , "Bacillota",
                                                          ifelse(fit_data_age_ref12m$results$Phylum == "Proteobacteria" , "Pseudomonadota",
                                                                 ifelse(fit_data_age_ref12m$results$Phylum == "Verrucomicrobia" , "Verrucomicrobiota",
                                                                        fit_data_age_ref12m$results$Phylum)))))


#keep only significant results
age_ref12m_res.sig = fit_data_age_ref12m$results[which(fit_data_age_ref12m$results$sig != ""),]; age_ref12m_res.sig; dim(age_ref12m_res.sig)
age_ref12m_res.sig = fit_data_age_ref12m$results %>%
  filter(sig %in% c("***","****")); age_ref12m_res.sig
unique(age_ref12m_res.sig$Genus); length(unique(age_ref12m_res.sig$Genus))
# age_ref12m_res.sig %>%
#   count(Phylum) %>%
#   arrange(desc(n))
age_ref12m_res.sig$Species = sub("_"," ", age_ref12m_res.sig$Species)

#FIGURE 1.B####
#plot-Species
p.age.ref12m_spe = ggplot(age_ref12m_res.sig, aes(x=value, y=Species, fill=coef)) +
  geom_tile(aes(fill=coef),width=0.6)+
  #geom_dendro(res_forDendro_age) +
  facet_grid(Phylum ~ ., scales = "free", space = "free") + #switch = "y"
  geom_text(aes(label=sig), size=7, nudge_y=-0.1)+
  scale_fill_gradient2(mid="white", low="#5d8dc9", high="#df65b0", na.value="grey85", name="Coefficient")+
  theme_classic() +
  theme(axis.title=element_text(face="bold",size=16, color="black"),
        axis.text.x=element_text(colour = "black",  size = 12, vjust = 2, face = "bold"),
        axis.text.y.left =element_text(colour = "black", size = 13, face="italic"),
        axis.text.y.right =element_text(colour = "black", size = 14),
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.width = unit(0.8, "cm"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 14),
        legend.position = "bottom",
        #strip.clip = "off",
        strip.text.y = element_text(size = 14,face = "bold.italic",margin = margin(0.2,0.2,0.2,0.2, "cm"), angle = 0), 
        strip.background=element_rect(colour = "#5d8dc9",fill="#5d8dc9")) +
  scale_x_discrete(name="Infants' Age (Months)", 
                   labels = c("1-6 Months"="1-6","13-18 Months"="13-18"))+
  xlab("Infants' Age (Month)")+
  ylab("") +
  coord_cartesian(expand=FALSE) +
  ggtitle("Reference: 7-12 Months") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)); p.age.ref12m_spe

# 12 vs 18 ####
age_ref12m_res.sig_18m = age_ref12m_res.sig[which(age_ref12m_res.sig$value == "13-18 Months"),]; age_ref12m_res.sig_18m; dim(age_ref12m_res.sig_18m)
unique(age_ref12m_res.sig_18m$Genus); length(unique(age_ref12m_res.sig_18m$Genus))

age_ref12m_res.sig_18m_pos = age_ref12m_res.sig_18m[which(age_ref12m_res.sig_18m$coef >0),]; dim(age_ref12m_res.sig_18m_pos)

age_ref12m_res.sig_18m_neg = age_ref12m_res.sig_18m[which(age_ref12m_res.sig_18m$coef <0),]; dim(age_ref12m_res.sig_18m_neg)

length(age_ref18m_res.sig$Species)
unique(age_ref18m_res.sig$Species)

#BIRTH ####
#only BIRTH, with no ref ####
fit_data_age_brt = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/metaphlan4/maaslin_xoxo_metagenome_refBirth",
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
  output = "~/Documents/xoxo_article/files/metaphlan4/maaslin_xoxo_metagenome_ref18m",
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

# add taxa to fit data
fit_data_age_ref18m_brt$results <- fit_data_age_ref18m_brt$results %>%
  dplyr::left_join(df_taxa, by = "feature"); dim(fit_data_age_ref18m_brt$results)
head(fit_data_age_ref18m_brt$results)

#keep only significant results
age_ref18m_brt_res.sig = fit_data_age_ref18m_brt$results[which(fit_data_age_ref18m_brt$results$sig != ""),]; age_ref18m_brt_res.sig; dim(age_ref18m_brt_res.sig)
unique(age_ref18m_brt_res.sig$Genus); length(unique(age_ref18m_brt_res.sig$Genus))
age_ref18m_brt_res.sig %>%
  group_by(Phylum) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

length(age_ref18m_brt_res.sig$coef[which(age_ref18m_brt_res.sig$coef<0)]) #increased
unique(age_ref18m_brt_res.sig$Genus[which(age_ref18m_brt_res.sig$coef<0)])

length(age_ref18m_brt_res.sig$coef[which(age_ref18m_brt_res.sig$coef>0)]) #decreased
unique(age_ref18m_brt_res.sig$Genus[which(age_ref18m_brt_res.sig$coef>0)])

# Csec vs vaginal ####
age_ref18m_brt_res.sig_vag = age_ref18m_brt_res.sig[which(age_ref18m_brt_res.sig$value == "VAGINAL"),]; age_ref18m_brt_res.sig_vag; dim(age_ref18m_brt_res.sig_vag)
unique(age_ref18m_brt_res.sig_vag$Species); length(unique(age_ref18m_brt_res.sig_vag$Species))

age_ref18m_brt_res.sig[which(age_ref18m_brt_res.sig$value == "Csec"),]

age_ref18m_brt_res.sig_vag_neg = age_ref18m_brt_res.sig_vag[which(age_ref18m_brt_res.sig_vag$coef <0),]; dim(age_ref18m_brt_res.sig_vag_neg)
length(age_ref18m_brt_res.sig_vag_neg$feature)

age_ref18m_brt_res.sig_vag_pos = age_ref18m_brt_res.sig_vag[which(age_ref18m_brt_res.sig_vag$coef >0),]; dim(age_ref18m_brt_res.sig_vag_pos)
length(age_ref18m_brt_res.sig_vag_pos$feature)
unique(age_ref18m_brt_res.sig_vag_pos$Genus); length(unique(age_ref18m_brt_res.sig_vag_pos$Genus))

#ref 1-6 Months ####
fit_data_age_ref6m_brt = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/metaphlan4/maaslin_xoxo_metagenome_ref6m",
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

# add taxa to fit data
fit_data_age_ref6m_brt$results <- fit_data_age_ref6m_brt$results %>%
  dplyr::left_join(df_taxa, by = "feature"); dim(fit_data_age_ref6m_brt$results)
head(fit_data_age_ref6m_brt$results)

#keep only significant results
age_ref6m_brt_res.sig = fit_data_age_ref6m_brt$results[which(fit_data_age_ref6m_brt$results$sig != ""),]; age_ref6m_brt_res.sig; dim(age_ref6m_brt_res.sig)
unique(age_ref6m_brt_res.sig$Genus); length(unique(age_ref6m_brt_res.sig$Genus))
age_ref6m_brt_res.sig %>%
  count(Phylum) %>%
  arrange(desc(n))

# Csec vs vaginal ####
age_ref6m_brt_res.sig_vag = age_ref6m_brt_res.sig[which(age_ref6m_brt_res.sig$value == "VAGINAL"),]; age_ref6m_brt_res.sig_vag; dim(age_ref6m_brt_res.sig_vag)
unique(age_ref6m_brt_res.sig_vag$Genus); length(unique(aage_ref6m_brt_res.sig_vag$Genus))
age_ref6m_brt_res.sig_vag %>%
  group_by(Phylum) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

age_ref6m_brt_res.sig_vag_neg = age_ref6m_brt_res.sig_vag[which(age_ref6m_brt_res.sig_vag$coef <0),]; dim(age_ref6m_brt_res.sig_vag_neg)

age_ref6m_brt_res.sig_vag_pos = age_ref6m_brt_res.sig_vag[which(age_ref6m_brt_res.sig_vag$coef >0),]; dim(age_ref6m_brt_res.sig_vag_pos)
length(age_ref6m_brt_res.sig_vag_pos$feature)
unique(age_ref6m_brt_res.sig_vag_pos$Genus); length(unique(age_ref6m_brt_res.sig_vag_pos$Genus))

View(fit_data_age_ref6m_brt$results[which(fit_data_age_ref6m_brt$results$sig != "****"),])

#ref 7-12 Months ####
fit_data_age_ref12m_brt = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/metaphlan4/maaslin_xoxo_metagenome_ref12m",
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

# add taxa to fit data
fit_data_age_ref12m_brt$results <- fit_data_age_ref12m_brt$results %>%
  dplyr::left_join(df_taxa, by = "feature"); dim(fit_data_age_ref12m_brt$results)
head(fit_data_age_ref12m_brt$results)

#keep only significant results
age_ref12m_brt_res.sig = fit_data_age_ref12m_brt$results[which(fit_data_age_ref12m_brt$results$sig != ""),]; age_ref12m_brt_res.sig; dim(age_ref12m_brt_res.sig)
unique(age_ref12m_brt_res.sig$Genus); length(unique(age_ref12m_brt_res.sig$Genus))
age_ref12m_brt_res.sig %>%
  count(Phylum) %>%
  arrange(desc(n))

# Csec vs vaginal ####
age_ref12m_brt_res.sig_vag = age_ref12m_brt_res.sig[which(age_ref12m_brt_res.sig$value == "VAGINAL"),]; age_ref12m_brt_res.sig_vag; dim(age_ref12m_brt_res.sig_vag)
unique(age_ref12m_brt_res.sig_vag$Genus); length(unique(age_ref12m_brt_res.sig_vag$Genus))
age_ref12m_brt_res.sig_vag %>%
  group_by(Phylum) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

age_ref12m_brt_res.sig_vag_neg = age_ref12m_brt_res.sig_vag[which(age_ref12m_brt_res.sig_vag$coef <0),]; dim(age_ref12m_brt_res.sig_vag_neg)

age_ref12m_brt_res.sig_vag_pos = age_ref12m_brt_res.sig_vag[which(age_ref12m_brt_res.sig_vag$coef >0),]; dim(age_ref12m_brt_res.sig_vag_pos)

#save
save.image("~/Documents/xoxo_article/files/metaphlan4/b2_xoxo_metagenome_maaslin.RData")



