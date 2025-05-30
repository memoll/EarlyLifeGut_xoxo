###############################################################
# Explanatory analysis of 18S data                            #
# Data: Miseq-18S- xoxo                                       #
# Mona Parizadeh - May 2025                                   #
###############################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan") 
library(ggplot2); packageVersion("ggplot2") 
library(tidyverse); packageVersion("tidyverse") 
library(dplyr);packageVersion("dplyr") 
library(ggpubr); packageVersion("ggpubr") 

# Import data (rarefied) #### 
setwd("~/Documents/xoxo_article/files/18s/")
ps.vst = readRDS("ps_xoxo_18s_VST.rds") 
ps.vst 

#Remove infant XC10 #### 
sample_data(subset_samples(ps.vst, sample_data(ps.vst)$INFANT_ID == "XC10"))$MONTH_GROUP
#XC10 has only one sample during the 2nd 6-months and no samples during the 3rd 6-months of sampling and mode of birth and sex are NA
#remove it for any statistical analyses (any analyses other than characterization)
ps = subset_samples(ps.vst, sample_data(ps.vst)$INFANT_ID != "XC10")
ps = prune_taxa(taxa_sums(ps)>0, ps); ps
ntaxa(ps.vst) - ntaxa(ps)
nsamples(ps.vst) - nsamples(ps)

saveRDS(ps,"ps_xoxo_18s_noXC10.rds")

#Richness ####
#Alpha diversity measures ####
adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"))
head(adiv)

#shannon
shn.rich = cbind(estimate_richness(ps,measures = 'shannon'),
                 sample_data(ps))
ggplot(shn.rich, aes(x = BIRTH, y = Shannon, color=BIRTH)) +  
  geom_boxplot()
ggplot(shn.rich, aes(x = ANIMALS, y = Shannon, color=ANIMALS)) +  
  geom_boxplot()
ggplot(shn.rich, aes(x = MONTH, y = Shannon, color=MONTH)) +  
  geom_boxplot()
ggplot(shn.rich, aes(x = MONTH_GROUP, y = Shannon, color=MONTH_GROUP)) +  
  geom_boxplot()
ggplot(shn.rich, aes(x = SEX, y = Shannon, color=SEX)) +  
  geom_boxplot()
ggplot(shn.rich, aes(x = INFANT_ID, y = Shannon, color=SEX)) +  
  geom_boxplot()

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
x_positions.shn <- setNames(1:length(levels(emm_shn.brt$BIRTH)), levels(emm_shn.brt$BIRTH))

# Assign numeric xmin and xmax
contrasts_shn.brt$xmin <- x_positions.shn[contrasts_shn.brt$group1]
contrasts_shn.brt$xmax <- x_positions.shn[contrasts_shn.brt$group2]

# Dynamically adjust y.position
max_y.shn <- max(shn.rich$Shannon_norm, na.rm = TRUE)
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

#FIGURE 5.F ####
img_shn_birth<- ggplot(emm_shn.brt, aes(x = BIRTH, y = emmean, fill=as.factor(BIRTH)))+
  theme_bw()+
  #geom_violin(alpha=0.9)+
  geom_violin(trim = FALSE, alpha = 0.7,data = shn.rich, aes(x = BIRTH, y = Shannon)) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black",
               data = shn.rich, aes(x = BIRTH, y = Shannon)) +  # Smaller boxplot inside the violin plot
  geom_bracket(data = contrasts_shn.brt,
               bracket.nudge.y = 0.2, tip.length = 0.01,
               aes(xmin = xmin, xmax = xmax, y.position = 4.2,vjust = 0, #move the text up or down relative to the bracket
                   label = pval.symp), label.size = 12,inherit.aes = FALSE) +
  scale_color_manual(name = "BIRTH", values = c("#59A1A0","#8c6bb1"),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_fill_manual(name = "BIRTH", values = c("#59A1A0","#8c6bb1"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_x_discrete(labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_y_continuous(limits = c(-0.5, max_y.shn + 2.5)) + # Adjust the limits of the y-axis here
  geom_jitter(data = shn.rich, aes(x = BIRTH, y = Shannon),col = "black", alpha=0.3, size = 0.7) +
  labs(x = "Mode of Birth", y = "Eukaryotic Shannon Diversity", color = "Mode of Birth") +
  theme(axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=27, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=22, color="black",vjust = 5,face="bold"),
        axis.text.y=element_text(size=22, color="black")); img_shn_birth

#PERMANOVA ####
#make dataframe
df = as(sample_data(ps), "data.frame")
#bray-curtis distance
dis = phyloseq::distance(ps,  method = "bray")
sample_variables(ps)
set.seed(111)
adns = adonis2(dis ~ BIRTH*SEX*ANIMALS*END_BF*NUM_PEOPLE*MONTH_GROUP, df, permutations = 999, by = "terms", strata = df$INFANT_ID) #distance = bray
adns

# Filter the data to only keep rows where P_value >= 0.05 and order based on p-value and R2
adns_sig <- subset(adns,  `Pr(>F)`<= 0.05) %>% 
  arrange(desc(R2),desc(`Pr(>F)`)); adns_sig
adns_sig$var <- rownames(adns_sig); adns_sig$var
adns_sig$var <- c("Proximity to Animals", "Proximity to Animals: \nAge Category", "Age Category",  "Weaning Time",
                  "Sex", "Mode of Birth", "Household Size")
adns_sig

# Determine the var# Determine the significance level
significance <- ifelse(adns_sig$`Pr(>F)` <= 0.001, "***", 
                       ifelse(adns_sig$`Pr(>F)` <= 0.01, "**", 
                              ifelse(adns_sig$`Pr(>F)` <= 0.05, "*", 
                                     ifelse(adns_sig$`Pr(>F)` <= 0.1, ".", ""))))

adns_p = ggplot(adns_sig, aes(x = R2*100, y = reorder(var,R2))) +
  #geom_bar(stat = "identity", fill = "#5d8dc9") +
  geom_bar(stat = "identity", fill = "#59A1A0") +
  geom_text(aes(label = paste(significance)),  hjust = 0.5,  vjust = 1.5, size = 5, angle = 90) + 
  labs(x = "R-squared", y = "") +
  theme_bw() +
  geom_vline(xintercept = c(5,10), linetype = "dashed", color = "#094242",linewidth = 0.5) +
  theme(                             
    axis.text.x = element_text(size = 10),
    axis.text.y =  element_text(size = 14, face = "bold",colour="black"),
    axis.title = element_text(size = 14, face="bold")); adns_p 
adns_p +
  #plot_layout(ncol = 2, nrow = 1) +
  plot_annotation(tag_levels = list(c('G'), '1')) & 
  theme(plot.tag = element_text(size = 40, face = "bold"))

adns_p = ggplot(adns_sig, aes(x = R2, y = reorder(var,R2))) +
  geom_bar(stat = "identity", fill = "#59A1A0") + #width = 0.6
  geom_text(aes(label = paste(significance)),  hjust = 0.5,  vjust = 1.5, size = 7, angle = 90, color = "#094242") + 
  labs(x = "R-squared", y = "") +
  scale_y_discrete(position = "right") +
  theme_bw() +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "#094242",linewidth = 1) +
  theme(#aspect.ratio = 1/2, #reduce axis scale
    axis.text.x = element_text(size = 12),
    axis.text.y =  element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18, face="bold"))+
  scale_y_discrete(position = "right"); adns_p 

adns_p = ggplot(adns_sig, aes(x = R2, y = (var))) +
  geom_bar(stat = "identity", fill = "#59a1a0") +
  geom_label(aes(label = paste(paste0(round(R2*100,digits=1), "%;"),paste0(`Pr(>F)`,significance))), hjust = 0.25,  
             color = "#1a3030", fill = "#ffd966", size = 5.5, label.padding = unit(0.2, "lines"), label.size = 0.3) + 
  labs(x = "R-squared", y = "") +
  theme_bw() +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "#2c5050",linewidth = 2) +
  theme(                             
    legend.text = element_text(size = 12),      
    legend.title = element_text(size = 14, face="bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y =  element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16, face="bold"))+
  scale_y_discrete(position = "right"); adns_p 

#PCoA ordination ####
pcoa = ordinate(ps, method = "MDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps, pcoa, shape = "BIRTH", color = "MONTH") + geom_point() + ggtitle("PCoA") +
  geom_text(aes(label = MI_ID), check_overlap = TRUE, size = 3, nudge_y = -0.01) 
#scale_color_manual(name = "BIRTH", values = c("chartreuse4", "darkred"))
plot_ordination(ps, pcoa, shape = "ANIMALS", color = "ANIMALS") + geom_point() + ggtitle("PCoA") +
  geom_text(aes(label = MI_ID), check_overlap = TRUE, size = 3, nudge_y = -0.01) 
plot_ordination(ps, pcoa, shape = "BIRTH", color = "MONTH") + geom_point() + ggtitle("PCoA") +
  geom_text(aes(label = MI_ID), check_overlap = TRUE, size = 3, nudge_y = -0.01)

pcoa_euk = plot_ordination(
  physeq = ps,                                                        
  ordination = pcoa)+       
  #group by birth
  #stat_ellipse(aes(group = BIRTH), level = 0.95, linetype = "dashed", alpha = 0.6) +
  geom_point(aes(fill = INFANT_ID, shape = SEX), size = 3) +    
  scale_shape_manual(name="Mode of Birth", values = c(21, 22),
                     labels = c(Csec = "C-Section", VAGINAL = "Vaginal Birth"))+
  scale_fill_manual(values = c("#7CE3D8","darkgoldenrod1","mediumvioletred","indianred1","tan3",
                               "cornflowerblue","seagreen","red2","tan4","yellowgreen",
                               "#980043","#dd1c77","#df65b0","#fd8d3c",'#8c6bb1',
                               "#74c476",'#4d004b',"yellow"),
                    name="Age (Month)") +
  theme_classic() +  
  #stat_ellipse(aes(group = SEX))+
  theme(                             
    legend.text = element_text(size = 12),      
    legend.title = element_text(size = 12, face="bold"))+
  #facet_wrap(~BIRTH)+
  #legend.title = element_blank())+ #removes legend title
  #legend.background = element_rect(fill = "white", color = "black")  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12, face="bold"))+
  guides(fill = guide_legend(override.aes = list(shape = 21))) +   #fills legend points based on the fill command
  ggtitle("PCoA - Eukaryotes")
pcoa_euk

# Computing Bray-Curtis Dissimilarities and PCoA
comm_mat <- vegdist(otu_table(ps), "bray")
PCoA_comm_mat <- capscale(comm_mat ~ 1, distance = "bray")
PCoA_comm_mat$CA$eig[1:3]/sum(PCoA_comm_mat$CA$eig)
PCoA_scores <- scores(PCoA_comm_mat)$sites

pcoa_scores = scores(pcoa$values$Relative_eig)
#FIGURE 5.E ####
pcoa_euk_birth = plot_ordination(
  physeq = ps,                                                        
  ordination = pcoa)+       
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = BIRTH, fill=BIRTH), level = 0.95, linetype = "dashed", alpha = 0.4) + 
  geom_point(aes(fill = BIRTH, shape = BIRTH), size = 3) +    
  scale_fill_manual(name="Mode of Birth", values = c("#59A1A0","#8c6bb1"),
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
  annotate("text", x = 0, y = 0.52, label = "paste(italic(R) ^ 2, \" = 0.008; \", italic(p), \" < 0.001***\")", parse = TRUE, 
           size = 11, color = "black"); pcoa_euk_birth

#birth - age
pcoa_euk_birth_age = plot_ordination(
  physeq = ps,                                                        
  ordination = pcoa)+       
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = BIRTH, fill=BIRTH), level = 0.95, linetype = "dashed", alpha = 0.4) + 
  geom_point(aes(fill = BIRTH, shape = BIRTH), size = 2) +    
  scale_fill_manual(name="Mode of Birth", values = c("#59A1A0","#8c6bb1"),
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
    strip.text = element_text(size = 13, face="bold")); pcoa_euk_birth_age
#save_plot("~/Documents/xoxo_article/images/xoxo_18s_pcoa_birth_age.pdf",pcoa_euk_birth_age)

# Homogeneity of dispersion test ####
comm = as.matrix(otu_table(ps))
metadata = data.frame(sample_data(ps))
#birth ####
disper.brt = betadisper(dis, df$BIRTH, type = "centroid"); disper.brt
plot(disper.brt, hull = FALSE, ellipse = TRUE)
# Extract the distances to the centroid for each group
distances_to_centroid_brt <- disper.brt$distance
metadata$distances_to_centroid_brt <- distances_to_centroid_brt
glmm.disp.brt <- glmmTMB(distances_to_centroid_brt ~ BIRTH + (1 | INFANT_ID), 
                         data = metadata, family = Gamma(link = "log")); summary(glmm.disp.brt)
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
ggplot(metadata, aes(x = BIRTH, y = distances_to_centroid_brt)) +
  #geom_violin(trim = FALSE, fill = "gray80") +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  labs(
    x = "Mode of Birth",
    y = "Dispersion",
    title = "Group-wise Dispersion in Eukaryotic Composition"
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

#FIGURE 5.H ####
p.disper.birth <- ggplot(emm_disp.brt, aes(x = BIRTH, y = emmean, fill=as.factor(BIRTH)))+
  theme_bw()+
  #geom_violin(alpha=0.9)+
  geom_violin(trim = FALSE, alpha = 0.7,data = metadata, aes(x = BIRTH, y = distances_to_centroid_brt)) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black",
               data = metadata, aes(x = BIRTH, y = distances_to_centroid_brt)) +  # Smaller boxplot inside the violin plot
  geom_bracket(data = contrasts_disp.fun.brt,
               bracket.nudge.y = 0.2, tip.length = 0.01,
               aes(xmin = xmin, xmax = xmax, y.position = 0.595,vjust = 0, #move the text up or down relative to the bracket
                   label = pval.symp), label.size = 12,inherit.aes = FALSE,type="text") +
  scale_color_manual(name = "BIRTH", values = c("#59A1A0","#8c6bb1"),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_fill_manual(name = "BIRTH", values = c("#59A1A0","#8c6bb1"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_x_discrete(labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  geom_jitter(data = metadata, aes(x = BIRTH, y = distances_to_centroid_brt),col = "black", alpha=0.3, size = 0.7) +
  xlab("Mode of Birth")+
  ylab("Eukaryotic Dispersion")+
  theme(axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=25, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=22, color="black",vjust = 5,face="bold"),
        axis.text.y=element_text(size=22, color="black")); p.disper.birth

#Temporal Shannon & Diversity ####
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
Anova(lmm.shn.mnt.brt)

ggplot(shn.rich, aes(x = MONTH, y = Shannon, color = BIRTH)) +
  geom_point(alpha = 0.5) +  # plot raw data points
  geom_smooth(method = "lm", se = TRUE) +  # linear regression line with confidence interval
  theme_minimal() +
  labs(title = "Shannon Diversity over Time", x = "Month", y = "Shannon Diversity")

#FIGURE 5.G ####
img_shn_mnt_birth = 
  ggplot(data = shn.rich %>%
           arrange(MONTH, Shannon_norm) %>%
           dplyr::slice(1:nrow(.)), 
         mapping = aes(x = as.factor(MONTH), y = Shannon, color = BIRTH, fill = BIRTH)) +
  #ggplot(shn.rich, aes(x = MONTH, y = Shannon)) +
  geom_smooth(aes(group = BIRTH),method = "lm", se = TRUE) +
  #geom_point(data = lmm.shn.mnt.brt.emm.df, aes(x = MONTH, y = emmean), color = "#df65b0", size = 4, inherit.aes = FALSE, alpha = 0.9) +
  #geom_errorbar(data = lmm.shn.mnt.brt.emm.df, aes(x = MONTH, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "#df65b0", inherit.aes = FALSE) + ##68477a
  #geom_smooth(aes(group=BIRTH), alpha=0.3) +
  theme_light() +
  geom_jitter(width = 0.1,alpha=0.3,colour = "black")+
  scale_color_manual(name = "Mode of Birth", values = c("#59A1A0","#8c6bb1"),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_fill_manual(name = "Mode of Birth", values = c("#59A1A0","#8c6bb1"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_x_discrete(breaks = c(1,6,12,18)) +
  annotate("text", x = 1, y = 3.6,
           label = paste0("R²m = ", round(lmm.shn.mnt.brt.r$R2_marginal, 3),
                          "\nR²c = ", round(lmm.shn.mnt.brt.r$R2_conditional, 3),
                          "; CI = ", "95%"), hjust = 0, size = 9) +
  #ylim(0,2)+
  geom_jitter(alpha=0.5, size = 1) +
  labs(
    x = "Months", y = "Temporal Changes of Eukaryotic Shannon Diversity", color = "Mode of Birth") +
  theme(legend.text = element_text(size = 25),      
        legend.title = element_text(size = 27, face="bold"),
        axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=23, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=22, color="black",vjust = 2,face="bold"),
        axis.text.y=element_text(size=22, color="black")); img_shn_mnt_birth #result.shn$pval

#shannon - sex
set.seed(112)
result.shn.sex <- permuspliner(data = shn.rich, x = 'MONTH', y = 'Shannon',
                               cases = 'MI_ID', category = 'SEX', groups = c('Female',"Male"))
result.shn.sex$pval < 0.05 #ns
permuspliner.plot.permsplines(result.shn.sex, xvar="MONTH", yvar="Shannon")

#%FUNGI ####
ps.fun = subset_taxa(ps, Order == "Fungi")
ps.fun = prune_samples(sample_sums(ps.fun)>0,ps.fun);ps.fun

#Alpha diversity measures ####
adiv.fun <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps.fun, measures = "Shannon"))
head(adiv.fun)

#shannon
shn.fun.rich = cbind(estimate_richness(ps.fun,measures = 'shannon'),
                     sample_data(ps.fun))

#LMM ####
library(lme4)
lmm.shn.fun.brt <- lmer(Shannon ~ BIRTH + 
                          (1 | INFANT_ID), data = shn.fun.rich)
summary(lmm.shn.fun.brt)
plot(lmm.shn.fun.brt)

# Check residuals
lmm.shn.fun.brt.resid <- resid(lmm.shn.fun.brt)
lmm.shn.fun.brt.fitted <- fitted(lmm.shn.fun.brt)

#Check the independence of the model residuals with the variable
# Plot residuals vs. fitted
plot(lmm.shn.fun.brt.fitted, lmm.shn.fun.brt.resid)
abline(h = 0, col = "red") #good model (dhomogenous dispersion of the residuals) 
#the pattern of residuals does not depend on the variable

# Histogram of residuals
hist(lmm.shn.fun.brt.resid, breaks = 30, main = "Residuals Histogram")

# QQ plot for normality
qqnorm(lmm.shn.fun.brt.resid) 
qqline(lmm.shn.fun.brt.resid) #residuals are normally distributed

boxplot(lmm.shn.fun.brt.resid ~ MONTH_GROUP, data = shn.fun.rich, xlab = "age", ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#To compare the average alpha diversity for each age group (adjusted for infant_ID variation)
library(emmeans)
# Run the pairwise comparisons (Tukey's)
lmm.shn.fun.brt.emm <- emmeans(lmm.shn.fun.brt, pairwise ~ BIRTH)
#p-value is adjusted (Tukey's by default)
plot(lmm.shn.fun.brt.emm)
# View the estimated means and the comparisons
lmm.shn.fun.brt.emm$emmeans
lmm.shn.fun.brt.emm$contrasts

# Get the adjusted means
emm_shn.fun.brt <- as.data.frame(emmeans(lmm.shn.fun.brt, ~ BIRTH)); emm_shn.fun.brt
#extract p-values
library(stringr)
contrasts_shn.fun.brt <- as.data.frame(lmm.shn.fun.brt.emm$contrasts)
contrasts_shn.fun.brt <- contrasts_shn.fun.brt %>%
  mutate(
    contrast_clean = str_remove_all(contrast, "[\\(\\)]"),  # remove parentheses
    group_split = str_split_fixed(contrast_clean, " - ", 2),
    group1 = str_trim(group_split[,1]),
    group2 = str_trim(group_split[,2])
  )
contrasts_shn.fun.brt         

pval.symp.shn.fun.brt <- symnum(contrasts_shn.fun.brt$p.value, corr = FALSE,
                                cutpoints = c(0, .0001, .001,.01,.05, .1, 1),
                                symbols = c("****","***","**","*"," ","ns"))

contrasts_shn.fun.brt$pval.symp <- as.character(pval.symp.shn.fun.brt)

# Map factor levels to numeric positions
emm_shn.fun.brt$BIRTH <- factor(emm_shn.fun.brt$BIRTH)
x_positions.shn.fun <- setNames(1:length(levels(emm_shn.fun.brt$BIRTH)), levels(emm_shn.fun.brt$BIRTH))

# Assign numeric xmin and xmax
contrasts_shn.fun.brt$xmin <- x_positions.shn.fun[contrasts_shn.fun.brt$group1]
contrasts_shn.fun.brt$xmax <- x_positions.shn.fun[contrasts_shn.fun.brt$group2]

# Dynamically adjust y.position
max_y.shn.fun <- max(shn.fun.rich$Shannon_norm, na.rm = TRUE)
contrasts_shn.fun.brt$y.position <- max_y.shn + seq(0.1, 0.3, length.out = nrow(contrasts_shn.fun.brt))

# Plot-error bar
ggplot(emm_shn.fun.brt , aes(x = BIRTH, y = emmean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
  labs(
    x = "Birth Mode",
    y = "Adjusted Alpha Diversity",
    title = "Estimated Alpha Diversity by Birth Mode"
  ) +
  theme_minimal()

# Plot-error bar & boxplot
ggplot(emm_shn.fun.brt , aes(x = BIRTH, y = emmean)) +
  geom_violin(data = shn.fun.rich, aes(x = BIRTH, y = Shannon), fill = "gray", alpha = 0.5, width = 0.8) +
  geom_boxplot(data = shn.fun.rich, aes(x = BIRTH, y = Shannon), width = 0.2, outlier.shape = NA, alpha = 0.6) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
  labs(
    x = "Birth Mode",
    y = "Adjusted Alpha Diversity",
    title = "Estimated Alpha Diversity by Birth Mode"
  ) +
  theme_minimal() +
  stat_compare_means()

#FIGURE S5.F ####
img_shn.fun_birth<- ggplot(emm_shn.fun.brt, aes(x = BIRTH, y = emmean, fill=as.factor(BIRTH)))+
  theme_bw()+
  #geom_violin(alpha=0.9)+
  geom_violin(trim = FALSE, alpha = 0.7,data = shn.fun.rich, aes(x = BIRTH, y = Shannon)) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black",
               data = shn.fun.rich, aes(x = BIRTH, y = Shannon)) +  # Smaller boxplot inside the violin plot
  geom_bracket(data = contrasts_shn.fun.brt,
               bracket.nudge.y = 0.2, tip.length = 0.01,
               aes(xmin = xmin, xmax = xmax, y.position = 4.2,vjust = 0, #move the text up or down relative to the bracket
                   label = pval.symp), label.size = 12,inherit.aes = FALSE) +
  scale_color_manual(name = "BIRTH", values = c("#59A1A0","#8c6bb1"),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_fill_manual(name = "BIRTH", values = c("#59A1A0","#8c6bb1"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_x_discrete(labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_y_continuous(limits = c(-1, max_y.shn + 2.5)) + # Adjust the limits of the y-axis here
  geom_jitter(data = shn.fun.rich, aes(x = BIRTH, y = Shannon),col = "black", alpha=0.3, size = 0.7) +
  labs(x = "Mode of Birth", y = "Fungal Shannon Diversity", color = "Mode of Birth") +
  theme(axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=27, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=22, color="black",vjust = 5,face="bold"),
        axis.text.y=element_text(size=22, color="black")); img_shn.fun_birth

#Temporal Shannon Diversity ####
shapiro.test(shn.fun.rich$Shannon)
shn.fun.rich %>%
  group_by(BIRTH) %>%
  shapiro_test(Shannon)
#shn.fun.rich$Shannon_norm = sqrt(shn.fun.rich$Shannon) #OR log(metadata$Shannon)
library(readxl)
#LLM shannon - infant ID is random effect
shn.fun.rich$MONTH <- as.numeric(as.character(shn.fun.rich$MONTH))

lmm.shn.fun.mnt.brt <- lmer(Shannon ~ MONTH + BIRTH + (1 | INFANT_ID), data = shn.fun.rich)
# Show correlation matrix
print(summary(lmm.shn.fun.mnt.brt), correlation = TRUE)
# Or extract variance-covariance matrix
vcov(lmm.shn.fun.mnt.brt)
plot(lmm.shn.fun.mnt.brt) #check residuals & model fit

library(performance)
lmm.shn.fun.mnt.brt.r = r2(lmm.shn.fun.mnt.brt)  # returns marginal and conditional R²
#R²m = marginal R² (variance explained by fixed effects)
#R²c = conditional R² (variance explained by fixed + random effects)

lmm.shn.fun.mnt.brt.resid = resid(lmm.shn.fun.mnt.brt)
qqnorm(lmm.shn.fun.mnt.brt.resid); qqline(lmm.shn.fun.mnt.brt.resid)
summary(lmm.shn.fun.mnt.brt.resid)
plotNormalHistogram(lmm.shn.fun.mnt.brt.resid)

shap.lmm.shn.fun.mnt.brt.resid = shapiro.test(lmm.shn.fun.mnt.brt.resid); shap.lmm.shn.fun.mnt.brt.resid #test for normality
shap.lmm.shn.fun.mnt.brt.resid$p.value < 0.05 #normally distributed

#pairwise comparison
library(emmeans)
lmm.shn.fun.mnt.brt.emm = emmeans(lmm.shn.fun.mnt.brt, ~ MONTH + BIRTH)
pairs(lmm.shn.fun.mnt.brt.emm)
lmm.shn.fun.mnt.brt.emm.df = as.data.frame(lmm.shn.fun.mnt.brt.emm)

#overall effect of MONTH
anova(lmm.shn.fun.mnt.brt)
Anova(lmm.shn.fun.mnt.brt)

ggplot(shn.fun.rich, aes(x = MONTH, y = Shannon, color = BIRTH)) +
  geom_point(alpha = 0.5) +  # plot raw data points
  geom_smooth(method = "lm", se = TRUE) +  # linear regression line with confidence interval
  theme_minimal() +
  labs(title = "Shannon Diversity over Time", x = "Month", y = "Shannon Diversity")

#FIGURE S5.G ####
img_shn.fun_mnt_birth = 
  ggplot(data = shn.fun.rich %>%
           arrange(MONTH, Shannon_norm) %>%
           dplyr::slice(1:nrow(.)), 
         mapping = aes(x = as.factor(MONTH), y = Shannon, color = BIRTH, fill = BIRTH)) +
  #ggplot(shn.fun.rich, aes(x = MONTH, y = Shannon)) +
  geom_smooth(aes(group = BIRTH),method = "lm", se = TRUE) +
  theme_light() +
  geom_jitter(width = 0.1,alpha=0.3,colour = "black")+
  scale_color_manual(name = "Mode of Birth", values = c("#59A1A0","#8c6bb1"),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_fill_manual(name = "Mode of Birth", values = c("#59A1A0","#8c6bb1"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_x_discrete(breaks = c(1,6,12,18)) +
  annotate("text", x = 1, y = 3.2,
           label = paste0("R²m = ", round(lmm.shn.fun.mnt.brt.r$R2_marginal, 3),
                          "\nR²c = ", round(lmm.shn.fun.mnt.brt.r$R2_conditional, 3),
                          "; CI = ", "95%"), hjust = 0, size = 10) +
  geom_jitter(alpha=0.5, size = 1) +
  labs(
    x = "Months", y = "Temporal Changes of Fungal Shannon Diversity", color = "Mode of Birth") +
  theme(legend.text = element_text(size = 25),      
        legend.title = element_text(size = 27, face="bold"),
        axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=23, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=22, color="black",vjust = 2,face="bold"),
        axis.text.y=element_text(size=22, color="black")); img_shn.fun_mnt_birth #result.shn$pval

#PERMANOVA ####
#make dataframe
df.fun = as(sample_data(ps.fun), "data.frame")
#bray-curtis distance
dis.fun = phyloseq::distance(ps.fun,  method = "bray")
sample_variables(ps.fun)
set.seed(1113)
adns.fun = adonis2(dis.fun ~ MONTH_GROUP*BIRTH*SEX*NUM_PEOPLE*ANIMALS, df.fun, by="terms",strata = df.fun$INFANT_ID) #distance = bray
adns.fun

#PCoA ordination ####
pcoa.fun = ordinate(ps.fun, method = "MDS", k = 2, try = 100, distance = "bray")
# Computing Bray-Curtis Dissimilarities and PCoA
comm.fun_mat <- vegdist(otu_table(ps.fun), "bray")
PCoA_comm.fun_mat <- capscale(comm.fun_mat ~ 1, distance = "bray")
PCoA_comm.fun_mat$CA$eig[1:3]/sum(PCoA_comm.fun_mat$CA$eig)
PCoA_scores.fun <- scores(PCoA_comm.fun_mat)$sites

pcoa.fun_scores = scores(pcoa.fun$values$Relative_eig)
#FIGURE S2.I ####
pcoa_fun_birth = plot_ordination(
  physeq = ps.fun,                                                        
  ordination = pcoa.fun)+       
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = BIRTH, fill=BIRTH), level = 0.95, linetype = "dashed", alpha = 0.4) + 
  geom_point(aes(fill = BIRTH, shape = BIRTH), size = 3) +    
  scale_fill_manual(name="Mode of Birth", values = c("#59A1A0","#8c6bb1"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal"))+
  scale_shape_manual(name="Mode of Birth", values = c(21, 24),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal"))+
  theme_classic()+  
  labs(
    x = paste0("PCoA1 (", round(pcoa.fun_scores[1:1]/sum(pcoa.fun_scores)*100,digits=1), "%)"), 
    y = paste0("PCoA2 (", round(pcoa.fun_scores[2:2]/sum(pcoa.fun_scores)*100,digits=1), "%)"))+ 
  theme(                             
    legend.text = element_text(size = 25),      
    legend.title = element_text(size = 27, face="bold"),
    legend.position = "right",
    axis.text = element_text(size = 21),
    axis.title = element_text(size = 27, face="bold")) +
  annotate("text", x = 0, y = 0.55, label = "paste(italic(R) ^ 2, \" = 0.009; \", italic(p), \" < 0.001***\")", parse = TRUE, 
           size = 10, color = "black"); pcoa_fun_birth

#birth - age
pcoa_fun_birth_age = plot_ordination(
  physeq = ps.fun,                                                        
  ordination = pcoa.fun)+       
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = BIRTH, fill=BIRTH), level = 0.95, linetype = "dashed", alpha = 0.4) + 
  geom_point(aes(fill = BIRTH, shape = BIRTH), size = 2) +    
  scale_fill_manual(name="Mode of Birth", values = c("#59A1A0","#8c6bb1"),
                    labels = c(Csec = "C-Section", VAGINAL = "Vaginal Birth"))+
  scale_shape_manual(name="Mode of Birth", values = c(21, 22),
                     labels = c(Csec = "C-Section", VAGINAL = "Vaginal Birth"))+
  facet_wrap(~MONTH_GROUP, ncol = 1, strip.position = "right") +
  theme_classic() +  
  labs(x = paste0("PCoA1 (", round(PCoA_comm.fun_mat$CA$eig[1:1]/sum(PCoA_comm.fun_mat$CA$eig)*100,digits=1), "%)"), 
       y = paste0("PCoA2 (", round(PCoA_comm.fun_mat$CA$eig[2:1]/sum(PCoA_comm.fun_mat$CA$eig)*100,digits=1), "%)"))+ 
  theme(                             
    legend.text = element_text(size = 18),      
    legend.title = element_text(size = 16, face="bold"),
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 16, face="bold"),
    strip.text = element_text(size = 13, face="bold")); pcoa_fun_birth_age

# Homogeneity of dispersion test ####
comm.fun = as.matrix(otu_table(ps.fun))
metadata.fun = data.frame(sample_data(ps.fun))
#birth ####
disper.fun.brt = betadisper(dis.fun, df.fun$BIRTH, type = "centroid"); disper.fun.brt
plot(disper.fun.brt, hull = FALSE, ellipse = TRUE)
# Extract the distances to the centroid for each group
distances_to_centroid_fun.brt <- disper.fun.brt$distance
metadata.fun$distances_to_centroid_fun.brt <- distances_to_centroid_fun.brt
glmm.disp.fun.brt <- glmmTMB(distances_to_centroid_fun.brt ~ BIRTH + (1 | INFANT_ID), 
                             data = metadata, family = Gamma(link = "log")); summary(glmm.disp.fun.brt)
summary(glmm.disp.fun.brt)

glmm.resid.fun.brt <- residuals(glmm.disp.fun.brt)
glmm.fitted.fun.brt <- fitted(glmm.disp.fun.brt)
ggplot(data = data.frame(fitted =glmm.fitted.fun.brt, residuals = glmm.resid.fun.brt)) +
  geom_point(aes(x = fitted, y = residuals)) +
  geom_hline(yintercept = 0, color = "red") +
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals")

qqnorm(glmm.resid.fun.brt); qqline(glmm.resid.fun.brt)
summary(glmm.resid.fun.brt)
plotNormalHistogram(glmm.resid.fun.brt)

#check for overdispersion
library(performance)
check_overdispersion(glmm.disp.fun.brt) #No overdispersion detected.

glmm.disp.fun.brt.emm <- emmeans(glmm.disp.fun.brt, pairwise ~ BIRTH); summary(glmm.disp.fun.brt.emm)

# Get the adjusted means
emm_disp.fun.brt <- as.data.frame(emmeans(glmm.disp.fun.brt, ~ BIRTH)); emm_disp.fun.brt
#extract p-values
library(stringr)
contrasts_disp.fun.brt <- as.data.frame(glmm.disp.fun.brt.emm$contrasts)
contrasts_disp.fun.brt <- contrasts_disp.fun.brt %>%
  mutate(
    contrast_clean = str_remove_all(contrast, "[\\(\\)]"),  # remove parentheses
    group_split = str_split_fixed(contrast_clean, " - ", 2),
    group1 = str_trim(group_split[,1]),
    group2 = str_trim(group_split[,2])
  )
contrasts_disp.fun.brt         

pval.symp.disp.fun.brt <- symnum(contrasts_disp.fun.brt$p.value, corr = FALSE,
                                 cutpoints = c(0, .0001, .001,.01,.05, .1, 1),
                                 symbols = c("****","***","**","*"," ","ns"))

contrasts_disp.fun.brt$pval.symp <- as.character(pval.symp.disp.fun.brt)

# Map factor levels to numeric positions
emm_disp.fun.brt$BIRTH <- factor(emm_disp.fun.brt$BIRTH)
x_positions.disp.fun.brt <- setNames(1:length(levels(emm_disp.fun.brt$BIRTH)), levels(emm_disp.fun.brt$BIRTH))

# Assign numeric xmin and xmax
contrasts_disp.fun.brt$xmin <- x_positions.disp.fun.brt[contrasts_disp.fun.brt$group1]
contrasts_disp.fun.brt$xmax <- x_positions.disp.fun.brt[contrasts_disp.fun.brt$group2]

# Dynamically adjust y.position
max_y.disp.fun.brt <- max(distances_to_centroid_fun.brt, na.rm = TRUE)
contrasts_disp.fun.brt$y.position <- max_y.disp.fun.brt + seq(0.1, 0.3, length.out = nrow(contrasts_disp.fun.brt))

#plot
ggplot(metadata.fun, aes(x = BIRTH, y = distances_to_centroid_fun.brt)) +
  #geom_violin(trim = FALSE, fill = "gray80") +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  labs(
    x = "Mode of Birth",
    y = "Dispersion",
    title = "Group-wise Dispersion in Fungal Composition"
  ) +
  theme_minimal()

#plot
# Get the adjusted means
#Compute mean and standard error by group
summary_metadata.fun <- metadata.fun %>%
  group_by(BIRTH) %>%
  summarise(
    mean_disp = mean(distances_to_centroid_fun.brt),
    se_disp = sd(distances_to_centroid_fun.brt) / sqrt(n())
  )

#Plot with error bars
ggplot(summary_metadata.fun, aes(x = BIRTH, y = mean_disp)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_disp - se_disp, ymax = mean_disp + se_disp), width = 0.2) +
  labs(
    x = "Mode of Birth",
    y = "Mean Distance to Centroid",
    title = "Dispersion by Mode of Birth with Error Bars"
  ) +
  theme_minimal()

#FIGURE S2.H ####
p.disper.fun.birth <- ggplot(emm_disp.fun.brt, aes(x = BIRTH, y = emmean, fill=as.factor(BIRTH)))+
  theme_bw()+
  #geom_violin(alpha=0.9)+
  geom_violin(trim = FALSE, alpha = 0.7,data = metadata.fun, aes(x = BIRTH, y = distances_to_centroid_fun.brt)) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "black",
               data = metadata.fun, aes(x = BIRTH, y = distances_to_centroid_fun.brt)) +  # Smaller boxplot inside the violin plot
  geom_bracket(data = contrasts_disp.fun.brt,
               bracket.nudge.y = 0.2, tip.length = 0.01,
               aes(xmin = xmin, xmax = xmax, y.position = 0.62,vjust = 0, #move the text up or down relative to the bracket
                   label = pval.symp), label.size = 12,inherit.aes = FALSE,type="text") +
  scale_color_manual(name = "BIRTH", values = c("#59A1A0","#8c6bb1"),
                     labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_fill_manual(name = "BIRTH", values = c("#59A1A0","#8c6bb1"),
                    labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  scale_x_discrete(labels = c(Csec = "C-section", VAGINAL = "Vaginal")) +
  geom_jitter(data = metadata.fun, aes(x = BIRTH, y = distances_to_centroid_fun.brt),col = "black", alpha=0.3, size = 0.7) +
  xlab("Mode of Birth")+
  ylab("Fungal Dispersion")+
  theme(axis.title.x=element_text(face="bold",size=25, color="black"),
        axis.title.y=element_text(face="bold",size=25, color="black"),
        legend.position="none",
        axis.text.x=element_text(size=22, color="black",vjust = 5,face="bold"),
        axis.text.y=element_text(size=22, color="black")); p.disper.fun.birth

#save
save.image("~/Documents/xoxo_article/files/18s/a3_xoxo_18s_expl.RData")
