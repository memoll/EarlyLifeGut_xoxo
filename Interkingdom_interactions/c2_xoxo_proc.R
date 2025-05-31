###############################################################
# Procrustes analysis                                         #
# Data: Metagenomics and 18S - xoxo                           # 
# Mona Parizadeh - May 2025                                   #
###############################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") 
library(vegan); packageVersion("vegan") 
library(ade4);packageVersion("ade4") 
library(ggplot2); packageVersion("ggplot2") 

# Import non-rarefied data #### 
setwd("~/Documents/xoxo_article/files/")
#Load data ####
ps = readRDS("ps_bac_euk.rds")
ps_euk = subset_taxa(ps, tax_table(ps)[,1] == "Eukaryota"); ps_euk
ps_bac = subset_taxa(ps, tax_table(ps)[,1] == "Bacteria"); ps_bac


#PCoA ####
pcoa_euk = ordinate(ps_euk, method = "MDS", k = 2, try = 100, distance = "bray")
pcoa_bac = ordinate(ps_bac, method = "MDS", k = 2, try = 100, distance = "bray")

pcoa_euk_scores = scores(pcoa_euk$values$Relative_eig)
pcoa_bac_scores = scores(pcoa_bac$values$Relative_eig)

#Procrustes ####
proc_euk_bac = procrustes(pcoa_euk$vectors, pcoa_bac$vectors, scale = TRUE, symmetric = TRUE, scores = "sites"); proc_euk_bac 
#symmetric=TRUE: the configurations are scaled to equal dispersions and a symmetric version of the Procrustes statistic is computed.
summary(proc_euk_bac)
proc_euk_bac$scale
sort(residuals(proc_euk_bac),TRUE)
#plot
plot(proc_euk_bac, pch = 19, xlab = "Eukaryotes", ylab = "Bacteria", kind = 1) 
plot(proc_euk_bac, kind = 1, type = "text")
plot(proc_euk_bac, kind = 2) #This allows identification of samples with the worst fit. 
#The horizontal lines, from bottom to top, are the 25% (dashed), 50% (solid), and 75% (dashed) quantiles of the residuals.

#Test of significance
set.seed(443)
prot = protest(X = pcoa_euk$vectors, Y = pcoa_bac$vectors, scores = "sites", permutations = 999); prot #more ASVs in bacteria than eukaryotes
#performs symmetric Procrustes analysis repeatedly to estimate the significance of the Procrustes statistic

ctest <- data.frame(rda1=prot$Yrot[,1],
                    rda2=prot$Yrot[,2],
                    xrda1=prot$X[,1],
                    xrda2=prot$X[,2])


# Extract age values from ps based on row names of df
age_values <- sample_data(ps)[rownames(ctest), "MONTH_GROUP"]  
# Add age column to the dataframe 
ctest$MONTH_GROUP <- age_values$MONTH_GROUP

# Computing Bray-Curtis Dissimilarities and PCoA
comm_mat <- vegdist(otu_table(ps), "bray")
PCoA_comm_mat <- capscale(comm_mat ~ 1, distance = "bray")
PCoA_comm_mat$CA$eig[1:3]/sum(PCoA_comm_mat$CA$eig)
PCoA_scores <- scores(PCoA_comm_mat)$sites

#FIGURE 4.C ####
#plot
library(patchwork)
gg.age=ggplot(ctest) +
  theme_classic() +
  stat_ellipse(geom = "polygon",aes(x=xrda1, y=xrda2, group = MONTH_GROUP,  fill = MONTH_GROUP), #colour = MONTH_GROUP,
               level = 0.95, linetype = "dashed", alpha = 0.5, color = "black",show.legend = TRUE) + #color = "black", 
  scale_fill_manual(name="Infants' Age", values = c("#59A1A0","#c0a9cc","#5d8dc9")) +
  geom_point(aes(x=rda1, y=rda2), colour = "#797272", size = 3) + 
  geom_point(aes(x=xrda1, y=xrda2), colour = "#d95f02", size = 3) + #"#d35687"
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2), size=0.1,
               arrow = arrow(length = unit(0.01, "npc")),colour="black",inherit.aes=FALSE) +
  labs(x="RDA1",y="RDA2")+
  annotate("text", x = -0.01, y = 0.08, label = "paste(italic(r), \" = 0.76; \", italic(P), \" < 0.001***\")", 
           parse = TRUE, size = 7, color = "black") +
  theme(                             
    legend.text = element_text(size = 18),      
    legend.title = element_text(size = 20, face="bold"),
    legend.position = "right",
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20, face="bold")); gg.age

#Canonical Analysis ####
#analyze relationships between two rectangular data tables
pca.euk <- dudi.pca(pcoa_euk$vectors, scale=T, scan=F, nf=5)
pca.bac <- dudi.pca(pcoa_bac$vectors, scale=T, scan=F, nf=5)
co.in.data <- coinertia(pca.euk, pca.bac, scan=F, nf=5)
summary(co.in.data)
set.seed(111)
rtest = randtest(co.in.data, nrepet = 999);rtest #Monte-Carlo test 
#p-values for Procrustes Analysis are generated using a Monte Carlo simulation. 
#Sample identifiers are shuffled in one of the PC matrices, and the M2 value is re-computed nrepet times. 
#The proportion of M2 values that are equal to or lower than the actual M2 value is the Monte Carlo p-value.
plot(co.in.data)
#OR
set.seed(222)
rv1 <- RV.rtest(pca.euk$tab, pca.bac$tab, 99); rv1
plot(rv1)

#subset age categories ####
#% 1-6 months ####
ps_6m = subset_samples(ps, sample_data(ps)$MONTH_GROUP == "1-6 Months") 
ps_6m = prune_taxa(taxa_sums(ps_6m)>0, ps_6m); ps_6m

ps_euk_6m = subset_taxa(ps_6m, tax_table(ps_6m)[,1] == "Eukaryota") 
ps_euk_6m = prune_samples(sample_sums(ps_euk_6m)>0, ps_euk_6m); ps_euk_6m
ps_bac_6m = subset_taxa(ps_6m, tax_table(ps_6m)[,1] == "Bacteria") 
ps_bac_6m = prune_taxa(taxa_sums(ps_bac_6m)>0, ps_bac_6m); ps_bac_6m

#% 7-12 months ####
ps_12m = subset_samples(ps, sample_data(ps)$MONTH_GROUP == "7-12 Months") 
ps_12m = prune_taxa(taxa_sums(ps_12m)>0, ps_12m); ps_12m

ps_euk_12m = subset_taxa(ps_12m, tax_table(ps_12m)[,1] == "Eukaryota") 
ps_euk_12m = prune_taxa(taxa_sums(ps_euk_12m)>0, ps_euk_12m); ps_euk_12m
ps_bac_12m = subset_taxa(ps_12m, tax_table(ps_12m)[,1] == "Bacteria") 
ps_bac_12m = prune_taxa(taxa_sums(ps_bac_12m)>0, ps_bac_12m); ps_bac_12m

#% 13-18 months ####
ps_18m = subset_samples(ps, sample_data(ps)$MONTH_GROUP == "13-18 Months") 
ps_18m = prune_taxa(taxa_sums(ps_18m)>0, ps_18m); ps_18m

ps_euk_18m = subset_taxa(ps_18m, tax_table(ps_18m)[,1] == "Eukaryota") 
ps_euk_18m = prune_taxa(taxa_sums(ps_euk_18m)>0, ps_euk_18m); ps_euk_18m
ps_bac_18m = subset_taxa(ps_18m, tax_table(ps_18m)[,1] == "Bacteria") 
ps_bac_18m = prune_taxa(taxa_sums(ps_bac_18m)>0, ps_bac_18m); ps_bac_18m

#PCoA ####
#%1-6 months ####
pcoa_euk_6m = ordinate(ps_euk_6m, method = "MDS", k = 2, try = 100, distance = "bray")
pcoa_bac_6m = ordinate(ps_bac_6m, method = "MDS", k = 2, try = 100, distance = "bray")
#%7-12 months ####
pcoa_euk_12m = ordinate(ps_euk_12m, method = "MDS", k = 2, try = 100, distance = "bray")
pcoa_bac_12m = ordinate(ps_bac_12m, method = "MDS", k = 2, try = 100, distance = "bray")
#%13-18 months ####
pcoa_euk_18m = ordinate(ps_euk_18m, method = "MDS", k = 2, try = 100, distance = "bray")
pcoa_bac_18m = ordinate(ps_bac_18m, method = "MDS", k = 2, try = 100, distance = "bray")

#Procrustes ####
#%1-6 months ####
proc_euk_bac_6m = procrustes(pcoa_euk_6m$vectors, pcoa_bac_6m$vectors, scale = TRUE, symmetric = TRUE, scores = "sites"); proc_euk_bac_6m 
summary(proc_euk_bac_6m)
proc_euk_bac_6m$scale
sort(residuals(proc_euk_bac_6m),TRUE)
#plot
plot(proc_euk_bac_6m, pch = 19, xlab = "Bacteria", ylab = "Eukaryotes", kind = 1) 
plot(proc_euk_bac_6m, kind = 2)

#Test of significance
set.seed(4431)
prot_6m = protest(X = pcoa_euk_6m$vectors, Y = pcoa_bac_6m$vectors, scores = "species", permutations = 999); prot_6m
#str(prot_6m)

ctest_6m <- data.frame(rda1=prot_6m$Yrot[,1],
                       rda2=prot_6m$Yrot[,2],
                       xrda1=prot_6m$X[,1],
                       xrda2=prot_6m$X[,2])

# Create a data frame for legend
legend_df <- data.frame(
  group = c("Eukaryotes", "Bacteria"),
  color = c("#797272", "#d95f02")
  #color = c("#69dd6a", "#2e5894")
)

#FIGURE 4.D ####
# Computing Bray-Curtis Dissimilarities 
comm_mat_6m <- vegdist(otu_table(ps_6m), "bray")
PCoA_comm_mat_6m <- capscale(comm_mat_6m ~ 1, distance = "bray")
PCoA_comm_mat_6m$CA$eig[1:3]/sum(PCoA_comm_mat_6m$CA$eig)
PCoA_scores_6m <- scores(PCoA_comm_mat_6m)$sites

#plot
g6m=ggplot(ctest_6m) +
  theme_bw() +
  geom_point(aes(x = rda1, y = rda2, color = "Eukaryotes"), size = 3) + 
  geom_point(aes(x = xrda1, y = xrda2, color = "Bacteria"), size = 3) +
  scale_color_manual(name = "Microbial Communities",
                     values = legend_df$color,
                     labels = legend_df$group) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2), size=0.15,
               arrow = arrow(length = unit(0.01, "npc")),colour="#555555",inherit.aes=FALSE) +
  labs(x="RDA1",y="RDA2")+
  annotate("text", x = -0.02, y = 0.1, label = "paste(italic(r), \" = 0.708; \", italic(P), \" < 0.01**\")", #prot_6m
           parse = TRUE, size = 4.5, color = "black") +
  theme(                             
    legend.text = element_text(size = 18),      
    legend.title = element_text(size = 20, face="bold"),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16, face="bold")); g6m 

p.g6m = ggplot(ctest_6m) +
  theme_bw() +
  geom_point(aes(x = rda1, y = rda2, color = "Eukaryotes"), size = 3) + 
  geom_point(aes(x = xrda1, y = xrda2, color = "Bacteria"), size = 3) +
  scale_color_manual(name = "Microbial Communities",
                     values = legend_df$color,
                     labels = legend_df$group) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2), size=0.15,
               arrow = arrow(length = unit(0.01, "npc")),colour="#555555",inherit.aes=FALSE) +
  # labs(
  #   x = paste0("PCoA1 (", round(PCoA_comm_mat_6m$CA$eig[1:1]/sum(PCoA_comm_mat_6m$CA$eig)*100,digits=1), "%)"), 
  #   y = paste0("PCoA2 (", round(PCoA_comm_mat_6m$CA$eig[2:1]/sum(PCoA_comm_mat_6m$CA$eig)*100,digits=1), "%)"))+ 
  labs(x="RDA1",y="RDA2")+
  annotate("text", x = -0.02, y = 0.1, label = "paste(italic(r), \" = 0.708; \", italic(P), \" < 0.01**\")", #prot_6m
           parse = TRUE, size = 7, color = "black") +
  theme(                             
    legend.text = element_text(size = 18),      
    legend.title = element_text(size = 20, face="bold"),
    legend.position = "right",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20, face="bold")); p.g6m 

#%7-12 months ####
proc_euk_bac_12m = procrustes(pcoa_euk_12m$vectors, pcoa_bac_12m$vectors, scale = TRUE, symmetric = TRUE, scores = "sites"); proc_euk_bac_12m 
summary(proc_euk_bac_12m)
proc_euk_bac_12m$scale
sort(residuals(proc_euk_bac_12m),TRUE)
#plot
plot(proc_euk_bac_12m, pch = 19, xlab = "Eukaryotes", ylab = "Bacteria", kind = 1) 
plot(proc_euk_bac_12m, kind = 2)

#Test of significance
set.seed(4432)
prot_12m = protest(X = pcoa_euk_12m$vectors, Y = pcoa_bac_12m$vectors, scores = "species", permutations = 999); prot_12m

ctest_12m <- data.frame(rda1=prot_12m$Yrot[,1],
                        rda2=prot_12m$Yrot[,2],
                        xrda1=prot_12m$X[,1],
                        xrda2=prot_12m$X[,2])

#FIGURE 4.E ####
# Computing Bray-Curtis Dissimilarities 
comm_mat_12m <- vegdist(otu_table(ps_12m), "bray")
PCoA_comm_mat_12m <- capscale(comm_mat_12m ~ 1, distance = "bray")
PCoA_comm_mat_12m$CA$eig[1:3]/sum(PCoA_comm_mat_12m$CA$eig)
PCoA_scores_12m <- scores(PCoA_comm_mat_12m)$sites

#plot
g12m=ggplot(ctest_12m) +
  theme_bw() +
  geom_point(aes(x = rda1, y = rda2, color = "Eukaryotes"), size = 3) + 
  geom_point(aes(x = xrda1, y = xrda2, color = "Bacteria"), size = 3) +
  scale_color_manual(name = "Microbial Communities",
                     values = legend_df$color,
                     labels = legend_df$group) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2), size=0.15,
               arrow = arrow(length = unit(0.01, "npc")),colour="#555555",inherit.aes=FALSE) +
  labs(x="RDA1",y="RDA2")+
  annotate("text", x = -0.02, y = 0.1, label = "paste(italic(r), \" = 0.708; \", italic(P), \" = ns \")",
           parse = TRUE, size = 4.5, color = "black") +
  theme(                             
    legend.text = element_text(size = 18),      
    legend.title = element_text(size = 20, face="bold"),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16, face="bold")); g12m 

p.g12m = ggplot(ctest_12m) +
  theme_bw() +
  geom_point(aes(x = rda1, y = rda2, color = "Eukaryotes"), size = 3) + 
  geom_point(aes(x = xrda1, y = xrda2, color = "Bacteria"), size = 3) +
  scale_color_manual(name = "Microbial Communities",
                     values = legend_df$color,
                     labels = legend_df$group) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2), size=0.15,
               arrow = arrow(length = unit(0.01, "npc")),colour="#555555",inherit.aes=FALSE) +
  labs(
    x = paste0("\nPCoA1 (", round(PCoA_comm_mat_12m$CA$eig[1:1]/sum(PCoA_comm_mat_12m$CA$eig)*100,digits=1), "%)"), 
    y = paste0("PCoA2 (", round(PCoA_comm_mat_12m$CA$eig[2:1]/sum(PCoA_comm_mat_12m$CA$eig)*100,digits=1), "%)"))+ 
  annotate("text", x = -0.02, y = 0.1, label = "paste(italic(r), \" = 0.708; \", italic(P), \" : ns \")",
           parse = TRUE, size = 5.5, color = "black") +
  theme(                             
    legend.text = element_text(size = 18),      
    legend.title = element_text(size = 20, face="bold"),
    legend.position = "right",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20, face="bold")); p.g12m 

#%13-18 months ####
proc_euk_bac_18m = procrustes(pcoa_euk_18m$vectors, pcoa_bac_18m$vectors, scale = TRUE, symmetric = TRUE, scores = "sites"); proc_euk_bac_18m 
summary(proc_euk_bac_18m)
proc_euk_bac_18m$scale
sort(residuals(proc_euk_bac_18m),TRUE)
#plot
plot(proc_euk_bac_18m, pch = 19, xlab = "Eukaryotes", ylab = "Bacteria", kind = 1) 
plot(proc_euk_bac_18m, kind = 2)

#Test of significance
set.seed(4433)
prot_18m = protest(X = pcoa_euk_18m$vectors, Y = pcoa_bac_18m$vectors, scores = "species", permutations = 999); prot_18m

ctest_18m <- data.frame(rda1=prot_18m$Yrot[,1],
                        rda2=prot_18m$Yrot[,2],
                        xrda1=prot_18m$X[,1],
                        xrda2=prot_18m$X[,2])

#FIGURE 4.F ####
# Computing Bray-Curtis Dissimilarities 
comm_mat_18m <- vegdist(otu_table(ps_18m), "bray")
PCoA_comm_mat_18m <- capscale(comm_mat_18m ~ 1, distance = "bray")
PCoA_comm_mat_18m$CA$eig[1:3]/sum(PCoA_comm_mat_18m$CA$eig)
PCoA_scores_18m <- scores(PCoA_comm_mat_18m)$sites

#plot
g18m = ggplot(ctest_18m) +
  theme_bw() +
  geom_point(aes(x = rda1, y = rda2, color = "Eukaryotes"), size = 3) + 
  geom_point(aes(x = xrda1, y = xrda2, color = "Bacteria"), size = 3) +
  scale_color_manual(name = "Microbial Communities",
                     values = legend_df$color,
                     labels = legend_df$group) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2), size=0.15,
               arrow = arrow(length = unit(0.01, "npc")),colour="#555555",inherit.aes=FALSE) +
  labs(x="RDA1",y="RDA2")+
  annotate("text", x = 0.01, y = 0.13, label = "paste(italic(r), \" = 0.556; \", italic(P), \" < 0.05*\")",
           parse = TRUE, size = 4.5, color = "black") +
  theme(                             
    legend.text = element_text(size = 18),      
    legend.title = element_text(size = 20, face="bold"),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16, face="bold")); g18m 

p.g18m = ggplot(ctest_18m) +
  theme_bw() +
  geom_point(aes(x = rda1, y = rda2, color = "Eukaryotes"), size = 3) + 
  geom_point(aes(x = xrda1, y = xrda2, color = "Bacteria"), size = 3) +
  scale_color_manual(name = "Microbial Communities",
                     values = legend_df$color,
                     labels = legend_df$group) +
  geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2), size=0.15,
               arrow = arrow(length = unit(0.01, "npc")),colour="#555555",inherit.aes=FALSE) +
  labs(
    x = paste0("\nPCoA1 (", round(PCoA_comm_mat_18m$CA$eig[1:1]/sum(PCoA_comm_mat_18m$CA$eig)*100,digits=1), "%)"), 
    y = paste0("PCoA2 (", round(PCoA_comm_mat_18m$CA$eig[2:1]/sum(PCoA_comm_mat_18m$CA$eig)*100,digits=1), "%)"))+ 
  annotate("text", x = 0.01, y = 0.11, label = "paste(italic(r), \" = 0.556; \", italic(P), \" < 0.05*\")",
           parse = TRUE, size = 5.5, color = "black") +
  theme(                             
    legend.text = element_text(size = 18),      
    legend.title = element_text(size = 20, face="bold"),
    legend.position = "right",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20, face="bold"));p.g18m

save.image("~/Documents/xoxo_article/files/c2_xoxo_proc.RData")

