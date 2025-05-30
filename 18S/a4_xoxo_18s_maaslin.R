###############################################################
# Differential analysis of 18S data - MaAsLin2                #
# Data: Miseq-18S- xoxo.                                      # 
# Mona Parizadeh - May 2025                                   #
###############################################################

# MaAsLin2 using a compound Poisson linear model 
# Run the MaAsLin2 models for normalized and relative abundance
# Outcome: normalized abundance of each ASV/genera
# Independent variable: age 

#load libraries ####
library(phyloseq); packageVersion("phyloseq") 
library(vegan); packageVersion("vegan") 
library(ggplot2); packageVersion("ggplot2") 
library(dplyr); packageVersion("dplyr") 
library(Maaslin2); packageVersion("Maaslin2") 

#import data ####
ps = readRDS("~/Documents/xoxo_article/files/18s/ps_xoxo_18s_noXC10.rds");ps 
length(unique(tax_table(ps)[,6])) 
ps.gen = tax_glom(ps, taxrank = "Genus"); ps.gen 

df_input_data <- data.frame(otu_table(ps.gen)); dim(df_input_data); head(df_input_data)
df_input_metadata <- data.frame(sample_data(ps.gen)); dim(df_input_metadata); df_input_metadata
df_taxa <- data.frame(tax_table(ps.gen)); dim(df_taxa); df_taxa
#add otus from rownames as a new column called feature
df_taxa$feature <- rownames(df_taxa)

#AGE ####
#ref 13-18 Months ####
fit_data_age_ref18m = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/18s/maaslin_xoxo_18s_ref18m",
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
fit_data_age_ref18m$results$sig <- ifelse(fit_data_age_ref18m$results$qval >= 0.1, "",
                                          ifelse(fit_data_age_ref18m$results$qval < 0.1 & fit_data_age_ref18m$results$qval >= 0.05, "*",
                                                 ifelse(fit_data_age_ref18m$results$qval < 0.05 & fit_data_age_ref18m$results$qval >= 0.01, "**",
                                                        ifelse(fit_data_age_ref18m$results$qval < 0.01 & fit_data_age_ref18m$results$qval >= 0.001, "***",
                                                               ifelse(fit_data_age_ref18m$results$qval < 0.001, "****", "error")))))

# add taxa to fit data
fit_data_age_ref18m$results <- fit_data_age_ref18m$results %>%
  left_join(df_taxa, by = "feature"); dim(fit_data_age_ref18m$results)
head(fit_data_age_ref18m$results)

#FINAL RESULTS ####
#keep only significant results
age_ref18m_res.sig = fit_data_age_ref18m$results[which(fit_data_age_ref18m$results$sig != ""),]; age_ref18m_res.sig; dim(age_ref18m_res.sig)
unique(age_ref18m_res.sig$Genus); length(unique(age_ref18m_res.sig$Genus))
age_ref18m_res.sig %>%
  group_by(Phylum) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

#RESULTS#### 
# 18 vs 6 ####
age_ref18m_res.sig_6m = age_ref18m_res.sig[which(age_ref18m_res.sig$value == "1-6 Months"),]; age_ref18m_res.sig_6m; dim(age_ref18m_res.sig_6m)
unique(age_ref18m_res.sig_6m$Genus); length(unique(age_ref18m_res.sig_6m$Genus))

age_ref18m_res.sig_6m_neg = age_ref18m_res.sig_6m[which(age_ref18m_res.sig_6m$coef <0),]; dim(age_ref18m_res.sig_6m_neg)
length(age_ref18m_res.sig_6m_neg$feature)
unique(age_ref18m_res.sig_6m_neg$Genus); length(unique(age_ref18m_res.sig_6m_neg$Genus))
length(which(age_ref18m_res.sig_6m_neg$Genus ==  "Clavispora-Candida_clade"))
length(which(age_ref18m_res.sig_6m_neg$Genus == "Cryptosporidium"))
length(which(age_ref18m_res.sig_6m_neg$Genus == "Geotrichum"))
length(which(age_ref18m_res.sig_6m_neg$Genus ==  "Kurtzmaniella-Candida_clade"))
length(which(age_ref18m_res.sig_6m_neg$Genus == "Saccharomyces"))

age_ref18m_res.sig_6m_pos = age_ref18m_res.sig_6m[which(age_ref18m_res.sig_6m$coef >0),]; dim(age_ref18m_res.sig_6m_pos)
length(age_ref18m_res.sig_6m_pos$feature)
unique(age_ref18m_res.sig_6m_pos$Genus); length(unique(age_ref18m_res.sig_6m_pos$Genus))
length(which(age_ref18m_res.sig_6m_pos$Genus == "Malassezia"))

# 18 vs 12 ####
age_ref18m_res.sig_12m = age_ref18m_res.sig[which(age_ref18m_res.sig$value == "7-12 Months"),]; age_ref18m_res.sig_12m; dim(age_ref18m_res.sig_12m)
unique(age_ref18m_res.sig_12m$Genus); length(unique(age_ref18m_res.sig_12m$Genus))

age_ref18m_res.sig_12m_neg = age_ref18m_res.sig_12m[which(age_ref18m_res.sig_12m$coef <0),]; dim(age_ref18m_res.sig_12m_neg)
length(age_ref18m_res.sig_12m_neg$feature)
unique(age_ref18m_res.sig_12m_neg$Genus); length(unique(age_ref18m_res.sig_12m_neg$Genus))
length(which(age_ref18m_res.sig_12m_neg$Genus ==  "Geotrichum"))
length(which(age_ref18m_res.sig_12m_neg$Genus == "Kurtzmaniella-Candida_clade"))

age_ref18m_res.sig_12m_pos = age_ref18m_res.sig_12m[which(age_ref18m_res.sig_12m$coef >0),]; dim(age_ref18m_res.sig_12m_pos)
length(age_ref18m_res.sig_12m_pos$feature)

#plot-OTU
p.age.ref18m = ggplot(age_ref18m_res.sig, aes(x=value, y=feature, fill=coef))+
  geom_tile()+
  #geom_dendro(res_forDendro_age) +
  geom_text(aes(label=sig), color="black", size=4.5, nudge_y=-0.3)+
  scale_fill_gradient2(mid="white", low="#59A1A0", high="#FFB90D", na.value="grey85", name="Coefficient")+
  theme_bw() +
  theme(axis.text.x=element_text(color = "black", angle=40, vjust=.95, hjust=0.95),
        #legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.8, "cm"))+
  xlab(NULL)+
  ylab(NULL) +
  scale_x_discrete(labels = c("7-12 Months" = "7-12 Months", "13-18 Months" = "13-18 Months")) +
  coord_cartesian(expand=FALSE); p.age.ref18m 

#FIGURE 
#plot-Genus
p.age.ref18m_gen = ggplot(age_ref18m_res.sig, aes(x=value, y=Genus, fill=coef)) +
  geom_tile()+
  #geom_dendro(res_forDendro_age) +
  facet_grid(Phylum ~ ., scales = "free", space = "free") + #switch = "y"
  geom_text(aes(label=sig), size=7, nudge_y=-0.3)+
  scale_fill_gradient2(mid="white", low="#59A1A0", high="#FFB90D", na.value="grey85", name="Coefficient")+
  theme_classic() +
  theme(axis.title=element_text(face="bold",size=18, color="black"),
        axis.text.x=element_text(colour = "black", face = "bold", size = 12),
        axis.text.y=element_text(colour = "black", size = 15, face = "italic"),
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.width = unit(0.8, "cm"),
        legend.title = element_text(colour="black", size=18, face="bold"),
        legend.text = element_text(colour="black", size = 16),
        legend.position = "bottom",
        #strip.clip = "off",
        strip.text.y = element_text(size = 16, face = "bold.italic", margin = margin(0,0.2,0,0.2, "cm"), angle = 0),
        strip.background=element_rect(colour = "#59A1A0",fill="#59A1A0")) +
  xlab("Infants' Age")+
  ylab("Microbial Eukaryotic Genera\n") +
  # scale_x_discrete(labels = c("7-12 Months" = "7-12 Months", "13-18 Months" = "13-18 Months")) +
  coord_cartesian(expand=FALSE) +
  ggtitle("Reference: 13-18 Months") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 25)); p.age.ref18m_gen

#ref 1-6 Months ####
fit_data_age_ref6m = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "~/Documents/xoxo_article/files/18s/maaslin_xoxo_18s_ref6m",
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
  left_join(df_taxa, by = "feature"); dim(fit_data_age_ref6m$results)
head(fit_data_age_ref6m$results)

#keep only significant results
age_ref6m_res.sig = fit_data_age_ref6m$results[which(fit_data_age_ref6m$results$sig != ""),]; age_ref6m_res.sig; dim(age_ref6m_res.sig)
unique(age_ref6m_res.sig$Genus); length(unique(age_ref6m_res.sig$Genus))
age_ref6m_res.sig %>%
  count(Phylum) %>%
  arrange(desc(n))

#RESULTS#### 
# 6 vs 12 ####
age_ref6m_res.sig_12m = age_ref6m_res.sig[which(age_ref6m_res.sig$value == "7-12 Months"),]; age_ref6m_res.sig_12m; dim(age_ref6m_res.sig_12m)
unique(age_ref6m_res.sig_12m$Genus); length(unique(age_ref6m_res.sig_12m$Genus))
age_ref6m_res.sig_12m %>%
  count(Phylum) %>%
  arrange(desc(n))

age_ref6m_res.sig_12m_neg = age_ref6m_res.sig_12m[which(age_ref6m_res.sig_12m$coef <0),]; dim(age_ref6m_res.sig_12m_neg)
length(age_ref6m_res.sig_12m_neg$feature)
unique(age_ref6m_res.sig_12m_neg$Genus); length(unique(age_ref6m_res.sig_12m_neg$Genus))
length(which(age_ref6m_res.sig_12m_neg$Genus == "Malassezia"))

age_ref6m_res.sig_12m_pos = age_ref6m_res.sig_12m[which(age_ref6m_res.sig_12m$coef >0),]; dim(age_ref6m_res.sig_12m_pos)
length(age_ref6m_res.sig_12m_pos$feature)
unique(age_ref6m_res.sig_12m_pos$Genus); length(unique(age_ref6m_res.sig_12m_pos$Genus))
length(which(age_ref6m_res.sig_12m_pos$Genus == "Pichia"))
length(which(age_ref6m_res.sig_12m_pos$Genus == "Saccharomyces"))

# 6 vs 18 ####
age_ref6m_res.sig_18m = age_ref6m_res.sig[which(age_ref6m_res.sig$value == "13-18 Months"),]; age_ref6m_res.sig_18m; dim(age_ref6m_res.sig_18m)
unique(age_ref6m_res.sig_18m$Genus); length(unique(age_ref6m_res.sig_18m$Genus))

age_ref6m_res.sig_18m_neg = age_ref6m_res.sig_18m[which(age_ref6m_res.sig_18m$coef <0),]; dim(age_ref6m_res.sig_18m_neg)
length(age_ref6m_res.sig_18m_neg$feature)
unique(age_ref6m_res.sig_18m_neg$Genus); length(unique(age_ref6m_res.sig_18m_neg$Genus))
length(which(age_ref6m_res.sig_18m_neg$Genus == "Malassezia"))

age_ref6m_res.sig_18m_pos = age_ref6m_res.sig_18m[which(age_ref6m_res.sig_18m$coef >0),]; dim(age_ref6m_res.sig_18m_pos)
length(age_ref6m_res.sig_18m_pos$feature)
unique(age_ref6m_res.sig_18m_pos$Genus); length(unique(age_ref6m_res.sig_18m_pos$Genus))
length(which(age_ref6m_res.sig_18m_pos$Genus == "Clavispora-Candida_clade")) 
length(which(age_ref6m_res.sig_18m_pos$Genus == "Cryptosporidium")) 
length(which(age_ref6m_res.sig_18m_pos$Genus == "Geotrichum")) 
length(which(age_ref6m_res.sig_18m_pos$Genus == "Kurtzmaniella-Candida_clade")) 
length(which(age_ref6m_res.sig_18m_pos$Genus == "Saccharomyces")) 

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

#FIGURE 2.C ####
age_ref6m_res.sig = fit_data_age_ref6m$results %>%
  filter(sig %in% c("***","****")); age_ref6m_res.sig
unique(age_ref6m_res.sig$Genus); length(unique(age_ref6m_res.sig$Genus))
age_ref6m_res.sig$Genus = sub("_"," ", age_ref6m_res.sig$Genus)
#plot-Genus
p.age.ref6m_gen = ggplot(age_ref6m_res.sig, aes(x=value, y=Genus, fill=coef)) +
  geom_tile()+
  #geom_dendro(res_forDendro_age) +
  facet_grid(Phylum ~ ., scales = "free", space = "free") + #switch = "y"
  geom_text(aes(label=sig), size=7, nudge_y=0)+
  scale_fill_gradient2(mid="white", low="#59A1A0", high="#FFB90D", na.value="grey85", name="Coefficient")+
  theme_classic() +
  theme(axis.title=element_text(face="bold",size=16, color="black"),
        axis.text.x=element_text(colour = "black",  size = 12, vjust = 2, face = "bold",margin = margin(0.2,0.2,0.2,0.2, "cm")),
        axis.text.y.left =element_text(colour = "black", size = 13, face="italic"),
        axis.text.y.right =element_text(colour = "black", size = 14),
        # legend.key.height = unit(0.3, "cm"),
        # legend.key.width = unit(0.8, "cm"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 14),
        legend.position = "bottom",
        #strip.clip = "off",
        strip.text.y.left = element_text(size = 14,face = "bold.italic", angle = 0), 
        strip.text.y.right = element_text(size = 14,face = "bold",angle = 0), 
        strip.background=element_rect(colour = "#59A1A0",fill="#59A1A0")) +  
  scale_x_discrete(name="Infants' Age (Months)", 
                   labels = c("7-12 Months"="7-12","13-18 Months"="13-18"))+
  xlab("Infants' Age (Months)")+
  ylab("Eukaryotic Genera") +
  coord_cartesian(expand=FALSE) +
  ggtitle("Reference: 1-6 Months") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)); p.age.ref6m_gen

#ref 12-18 Months ####
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

#keep only significant results
age_ref12m_res.sig = fit_data_age_ref12m$results[which(fit_data_age_ref12m$results$sig != ""),]; age_ref12m_res.sig; dim(age_ref12m_res.sig)
unique(age_ref12m_res.sig$Genus); length(unique(age_ref12m_res.sig$Genus))
age_ref12m_res.sig %>%
  count(Phylum) %>%
  arrange(desc(n))

# 12 vs 18 ####
age_ref12m_res.sig_18m = age_ref12m_res.sig[which(age_ref12m_res.sig$value == "13-18 Months"),]; age_ref12m_res.sig_18m; dim(age_ref12m_res.sig_18m)
unique(age_ref12m_res.sig_18m$Genus); length(unique(age_ref12m_res.sig_18m$Genus))
age_ref12m_res.sig_18m %>%
  count(Phylum) %>%
  arrange(desc(n))

#FIGURE 2.D ####
age_ref12m_res.sig = fit_data_age_ref12m$results %>%
  filter(sig %in% c("***","****")); age_ref12m_res.sig
unique(age_ref12m_res.sig$Genus); length(unique(age_ref12m_res.sig$Genus))
age_ref12m_res.sig$Genus = sub("_"," ", age_ref12m_res.sig$Genus)

p.age.ref12m_gen = ggplot(age_ref12m_res.sig, aes(x=value, y=Genus, fill=coef)) +
  geom_tile(aes(fill=coef),width=0.6)+
  #geom_dendro(res_forDendro_age) +
  facet_grid(Phylum ~ ., scales = "free", space = "free") + #switch = "y"
  geom_text(aes(label=sig), size=7, nudge_y=0)+
  scale_fill_gradient2(mid="white", low="#59A1A0", high="#FFB90D", na.value="grey85", name="Coefficient")+
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
        strip.background=element_rect(colour = "#59A1A0",fill="#59A1A0")) +  
  scale_x_discrete(name="Infants' Age (Months)", 
                   labels = c("1-6 Months"="1-6","13-18 Months"="13-18"))+
  xlab("Infants' Age (Months)")+        
  ylab("") +
  coord_cartesian(expand=FALSE) +
  ggtitle("Reference: 7-12 Months") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)); p.age.ref12m_gen

#Age & BIRTH ####
#only BIRTH, with no ref ####
fit_data_age_brt = Maaslin2(
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

# Csec vs vaginal ####
age_ref18m_brt_res.sig_vag = age_ref18m_brt_res.sig[which(age_ref18m_brt_res.sig$value == "VAGINAL"),]; age_ref18m_brt_res.sig_vag; dim(age_ref18m_brt_res.sig_vag)
age_ref18m_brt_res.sig[which(age_ref18m_brt_res.sig$value == "Csec"),]

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
  group_by(Phylum) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

# Csec vs vaginal ####
age_ref6m_brt_res.sig_cs = age_ref6m_brt_res.sig[which(age_ref6m_brt_res.sig$value == "Csec"),]; age_ref6m_brt_res.sig_cs; dim(age_ref6m_brt_res.sig_cs)
age_ref6m_brt_res.sig[which(age_ref6m_brt_res.sig$value == "VAGINAL"),]

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
  group_by(Phylum) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

# Csec vs vaginal ####
age_ref12m_brt_res.sig_vag = age_ref12m_brt_res.sig[which(age_ref12m_brt_res.sig$value == "VAGINAL"),]; age_ref12m_brt_res.sig_vag; dim(age_ref12m_brt_res.sig_vag)
age_ref12m_brt_res.sig[which(age_ref12m_brt_res.sig$value == "Csec"),]

#save
save.image("~/Documents/xoxo_article/files/18s/a4_xoxo_18s_maaslin.RData")

