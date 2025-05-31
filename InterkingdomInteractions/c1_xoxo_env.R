###############################################################
# Fitness analysis                                            #
# Data: Metagenomics and 18S - xoxo                           # 
# Mona Parizadeh - May 2025                                   #
###############################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") 
library(vegan); packageVersion("vegan") 
library(ggplot2); packageVersion("ggplot2") 
library(tidyverse); packageVersion("tidyverse") 
library(dplyr);packageVersion("dplyr") 

# Import non-rarefied data #### 
setwd("~/Documents/xoxo_article/files/")
#Load data ####
ps = readRDS("ps_bac_euk_GEN.rds"); ps

#PERMANOVA ####
#make dataframe
df = as(sample_data(ps), "data.frame")
#bray-curtis distance
dis = phyloseq::distance(ps,  method = "bray")
set.seed(111)
adns = adonis2(dis ~ MONTH_GROUP, df, permutations = 999, strata = df$INFANT_ID) #distance = bray
adns

#Ordination ####
# Computing Bray-Curtis Dissimilarities and PCoA
comm_mat <- vegdist(otu_table(ps), "bray")
PCoA_comm_mat <- capscale(comm_mat ~ 1, distance = "bray")
PCoA_comm_mat$CA$eig[1:3]/sum(PCoA_comm_mat$CA$eig)
PCoA_scores <- scores(PCoA_comm_mat)$sites

pcoa_scores = scores(pcoa$values$Relative_eig)

#plot
pcoa = ordinate(ps, method = "MDS", k = 2, try = 100, distance = "bray")
pcoa1 = plot_ordination(ps,pcoa)
#Empty points of the PCoA replace them with the desired shapes 
pcoa1$layers
pcoa1$layers = pcoa1$layers[-1] #remove the original points to add the desired colors and shapes
library(ggrepel); packageVersion("ggrepel") #‘0.9.5’
p = pcoa1 +
  theme_bw() +      
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = MONTH_GROUP, fill = MONTH_GROUP, colour = MONTH_GROUP), 
               level = 0.95, alpha = 0.01, show.legend = FALSE,linetype = "dashed") + #linetype = "dashed" #alpha=0.7
  geom_point(aes(fill = MONTH_GROUP, shape = MONTH_GROUP), size = 3, alpha=0.8) +   
  scale_color_manual(name="Infants' Age", values = c("black","black","black")) +
  scale_fill_manual(name="Infants' Age", values = c("#59A1A0","#d95f02","#225ea8")) +
  scale_shape_manual(name="Age Category", values = c(21,24,22))+
  theme_classic() + 
  labs(x="PCoA1",y="PCoA2")+
  theme(                             
    legend.text = element_text(size = 12),      
    legend.title = element_text(size = 12, face="bold"),
    legend.position = "right",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18, face="bold"),
    title = element_text(size = 10, face="bold")) +
  annotate("text", x = 0, y = 0.5, label = "paste(italic(R) ^ 2, \" = 0.063; \", italic(p), \" < 0.001***\")", #permanova for age category (MONTH_GROUP)
           parse = TRUE, size = 5, color = "black"); p

p1 = pcoa1 +
  theme_bw() +      
  #group by birth
  stat_ellipse(geom = "polygon",aes(group = MONTH_GROUP, fill = MONTH_GROUP, colour = MONTH_GROUP), 
               level = 0.95, alpha = 0.01, show.legend = TRUE,linetype = "dashed") + #linetype = "dashed" #alpha=0.7
  geom_point(aes(fill = MONTH_GROUP, shape = MONTH_GROUP), size = 3, alpha=0.8) +   
  scale_color_manual(name="Age Category", values = c("#59A1A0","#c0a9cc","#5d8dc9")) +
  scale_fill_manual(name="Age Category", values = c("#59A1A0","#c0a9cc","#5d8dc9")) +
  scale_shape_manual(name="Age Category", values = c(21,24,22))+
  theme_classic() + 
  labs(
    x = paste0("PCoA1 (", round(PCoA_comm_mat$CA$eig[1:1]/sum(PCoA_comm_mat$CA$eig)*100,digits=1), "%)"), 
    y = paste0("PCoA2 (", round(PCoA_comm_mat$CA$eig[2:1]/sum(PCoA_comm_mat$CA$eig)*100,digits=1), "%)"))+ 
  theme(                             
    legend.position = "right",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18, face="bold"),
    title = element_text(size = 10, face="bold")); p1

# Fitness ####
comm = otu_table(ps)
taxo = tax_table(ps)

#Genus
library(seqtools); packageVersion("seqtools") #‘0.2.0’
comm.gen = taxocomm(comm, taxo, "Genus")
comm.gen.ra = decostand(comm.gen, method="total") #relative abundance
gen01 = sort(apply(comm.gen.ra[,apply(comm.gen.ra,2,mean)>0.01], 2, mean),TRUE)
gen01.ra = comm.gen.ra[,which(colnames(comm.gen.ra) %in% names(gen01))]

#TABLE ####
#fit vectors on pcoa
env.gen = envfit(as.data.frame(pcoa$vectors), gen01.ra)
env.gen
sort(env.gen$vectors$r, decreasing = TRUE)
#dataframe
env.gen.df = as.data.frame(env.gen$vectors$arrows*sqrt(env.gen$vectors$r))

#find kingdom
kingdoms_for_genera <- as.data.frame(tax_table(ps))[tax_table(ps)[,6] %in% 
                                                      rownames(env.gen.df), c("Genus", "Kingdom")]; kingdoms_for_genera
matched_kingdoms <- kingdoms_for_genera[match(rownames(env.gen.df), kingdoms_for_genera$Genus), ]; matched_kingdoms

# filter significant results (p = 0.001)
#env.gen.df[env.gen$vectors$pvals < 0.05, ]
env.gen.df.sig = env.gen.df[env.gen$vectors$pvals <= 0.001, ]; env.gen.df.sig

env.gen.df.sig$genus = rownames(env.gen.df.sig)

# Check if the value from the 'genus' column of df exists in tax_table(ps)$genus
matching_rows <- match(env.gen.df.sig$genus, tax_table(ps)[,6])

# Extract the corresponding Kingdom values
kingdom_values <- tax_table(ps)[,1][matching_rows]

# Add the 'kingdom' column to df
env.gen.df.sig$Kingdom <- kingdom_values

kingdom_colors <- c("Bacteria" = "#CCCCCC", "Eukaryota" = "#f1b197") 
env.gen.df.sig$kingdom_color <- kingdom_colors[env.gen.df.sig$Kingdom]

#plot (bacterial genilies on PCoA)
p.env.gen = plot_ordination(ps, pcoa) 
#Empty points of the PCoA replace them with the desired shapes 
p.env.gen$layers
p.env.gen$layers = p.env.gen$layers[-1] #remove the original points to add the desired colors and shapes


#FIGURE 4.A ####
#rename genera (correction):
env.gen.df.sig$genus = sub("_"," ", env.gen.df.sig$genus)
env.gen.df.sig$genus = sub("-"," ", env.gen.df.sig$genus)
gg = p1 +
  theme_classic() +      
  geom_segment(data=env.gen.df.sig,aes(x=0,xend=0.4*Axis.1,y=0,yend=0.4*Axis.2), size=0.4,
               arrow = arrow(length = unit(0.01, "npc")),colour="black",inherit.aes=FALSE) +
  ggrepel::geom_label_repel(data=env.gen.df.sig,aes(x=0.4*Axis.1,y=0.4*Axis.2,label=genus),size=6, fontface="bold.italic",
                            color = "black",fill = env.gen.df.sig$kingdom_color, label.padding = unit(0.4, "lines"), 
                            label.size = 0.1, show.legend = FALSE) +
  #scale_fill_discrete(guide = guide_legend(order = env.gen.df.sig$Kingdom, title = "Microbial Community", title.position = "top")) +
  labs(
    x = paste0("PCoA1 (", round(pcoa_scores[1:1]/sum(pcoa_scores)*100,digits=1), "%)"), 
    y = paste0("PCoA2 (", round(pcoa_scores[2:2]/sum(pcoa_scores)*100,digits=1), "%)"))+ 
  theme(                             
    legend.text = element_text(size = 18),      
    legend.title = element_text(size = 20, face="bold"),
    legend.position = "right",
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    title = element_text(size=18, face="bold")); gg

gg = p1 +
  #theme_classic() +     
  geom_segment(data=env.gen.df.sig,aes(x=0,xend=0.4*Axis.1,y=0,yend=0.4*Axis.2), size=0.4,
               arrow = arrow(length = unit(0.01, "npc")),colour="black",inherit.aes=FALSE) +
  ggrepel::geom_label_repel(data=env.gen.df.sig,aes(x=0.4*Axis.1,y=0.4*Axis.2,label=genus),size=6, fontface="bold.italic",
                            color = "black",fill = env.gen.df.sig$kingdom_color,alpha = 0.9, label.padding = unit(0.4, "lines"), 
                            label.size = 0.1, show.legend = FALSE) +
  theme(                             
    legend.position = "right",
    legend.text = element_text(size = 12),      
    legend.title = element_text(size = 14, face="bold")); gg

#subset samples - age ####
ps.6m = subset_samples(ps, sample_data(ps)$AGE_GROUP == "1_6Months")
ps.6m = prune_taxa(taxa_sums(ps.6m)>0, ps.6m); ps.6m
comm.6m = otu_table(ps.6m)
taxo.6m = tax_table(ps.6m)
ps.12m = subset_samples(ps, sample_data(ps)$AGE_GROUP == "7_12Months")
ps.12m = prune_taxa(taxa_sums(ps.12m)>0, ps.12m); ps.12m
comm.12m = otu_table(ps.12m)
taxo.12m = tax_table(ps.12m)
ps.18m = subset_samples(ps, sample_data(ps)$AGE_GROUP == "13_18Months")
ps.18m = prune_taxa(taxa_sums(ps.18m)>0, ps.18m); ps.18m
comm.18m = otu_table(ps.18m)
taxo.18m = tax_table(ps.18m)

# Fitness - age ####
#1-6 months
comm.gen.6m = taxocomm(comm.6m, taxo.6m, "Genus")
comm.gen.6m.ra = decostand(comm.gen.6m, method="total") #relative abundance
gen.6m.01 = sort(apply(comm.gen.6m.ra[,apply(comm.gen.6m.ra,2,mean)>0.01], 2, mean),TRUE)
gen.6m.01.ra = comm.gen.6m.ra[,which(colnames(comm.gen.6m) %in% names(gen.6m.01))]
#pcoa
pcoa.6m = ordinate(ps.6m, method = "MDS", k = 2, try = 100, distance = "bray")
#fit vectors on pcoa
env.gen.6m = envfit(as.data.frame(pcoa.6m$vectors), gen.6m.01.ra)
env.gen.6m
# Extract the environmental variable fit results
env.gen.6m.df <- as.data.frame(scores(env.gen.6m, "vectors"))
env.gen.6m.df$pval <- env.gen.6m$vectors$pvals  # Add p-values
env.gen.6m.df$r2 <- env.gen.6m$vectors$r # Add R^2 values
# Filter significant variables (p-value < 0.05)
env.gen.6m.df[env.gen.6m.df$pval <= 0.05, ]
env.gen.6m.sig <- env.gen.6m.df[env.gen.6m.df$pval <= 0.001, ]; env.gen.6m.sig
# # Order by R² in descending order (or by p-value in ascending order)
# env.gen.6m.sig.ord <- env.gen.6m.sig[order(-env.gen.6m.sig$r2, env.gen.6m.sig$pval), ]

#7-12 months
comm.gen.12m = taxocomm(comm.12m, taxo.12m, "Genus")
comm.gen.12m.ra = decostand(comm.gen.12m, method="total") #relative abundance
gen.12m.01 = sort(apply(comm.gen.12m.ra[,apply(comm.gen.12m.ra,2,mean)>0.01], 2, mean),TRUE)
gen.12m.01.ra = comm.gen.12m.ra[,which(colnames(comm.gen.12m) %in% names(gen.12m.01))]
#pcoa
pcoa.12m = ordinate(ps.12m, method = "MDS", k = 2, try = 100, distance = "bray")
#fit vectors on pcoa
env.gen.12m = envfit(as.data.frame(pcoa.12m$vectors), gen.12m.01.ra)
env.gen.12m
# Extract the environmental variable fit results
env.gen.12m.df <- as.data.frame(scores(env.gen.12m, "vectors"))
env.gen.12m.df$pval <- env.gen.12m$vectors$pvals  # Add p-values
env.gen.12m.df$r2 <- env.gen.12m$vectors$r # Add R^2 values
# Filter significant variables (p-value < 0.05)
env.gen.12m.df[env.gen.12m.df$pval <= 0.05, ]
env.gen.12m.sig <- env.gen.12m.df[env.gen.12m.df$pval <= 0.001, ]; env.gen.12m.sig
# # Order by R² in descending order (or by p-value in ascending order)
# env.gen.12m.sig.ord <- env.gen.12m.sig[order(-env.gen.12m.sig$r2, env.gen.12m.sig$pval), ]

#13-18 months
comm.gen.18m = taxocomm(comm.18m, taxo.18m, "Genus")
comm.gen.18m.ra = decostand(comm.gen.18m, method="total") #relative abundance
gen.18m.01 = sort(apply(comm.gen.18m.ra[,apply(comm.gen.18m.ra,2,mean)>0.01], 2, mean),TRUE)
gen.18m.01.ra = comm.gen.18m.ra[,which(colnames(comm.gen.18m) %in% names(gen.18m.01))]
#pcoa
pcoa.18m = ordinate(ps.18m, method = "MDS", k = 2, try = 100, distance = "bray")
#fit vectors on pcoa
env.gen.18m = envfit(as.data.frame(pcoa.18m$vectors), gen.18m.01.ra)
env.gen.18m
# Extract the environmental variable fit results
env.gen.18m.df <- as.data.frame(scores(env.gen.18m, "vectors"))
env.gen.18m.df$pval <- env.gen.18m$vectors$pvals  # Add p-values
env.gen.18m.df$r2 <- env.gen.18m$vectors$r # Add R^2 values
# Filter significant variables (p-value < 0.05)
env.gen.18m.df[env.gen.18m.df$pval <= 0.05, ]
env.gen.18m.sig <- env.gen.18m.df[env.gen.18m.df$pval <= 0.001, ]; env.gen.18m.sig
# # Order by R² in descending order (or by p-value in ascending order)
# env.gen.18m.sig.ord <- env.gen.18m.sig[order(-env.gen.18m.sig$r2, env.gen.18m.sig$pval), ]

# Combine all unique rownames
length(unique(c(rownames(env.gen.6m.sig), rownames(env.gen.12m.sig), rownames(env.gen.18m.sig))))
length(duplicated(c(rownames(env.gen.6m.sig), rownames(env.gen.12m.sig), rownames(env.gen.18m.sig))))
# List of all dataframes
env.gen.all.df <- list(env.gen.6m.sig,env.gen.12m.sig,env.gen.18m.sig)
# Inspect the structure of each dataframe
for (i in seq_along(env.gen.all.df)) {
  cat("Dataframe", i, "structure:\n")
  str(env.gen.all.df[[i]])
  cat("\n")
}
# Create a full list of all unique rownames
env.gen.all <- unique(unlist(lapply(env.gen.all.df, rownames)))
# Initialize the merged dataframe with default values (pval = 1, r2 = 0)
env.gen.all.mrg <- data.frame(rowname = env.gen.all); env.gen.all.mrg
# Add columns for pval and r2 from each dataframe
for (i in seq_along(env.gen.all.df)) {
  # Extract the dataframe for this iteration
  env.df <- env.gen.all.df[[i]]
  
  # Convert to dataframe if it's not already
  env.df <- as.data.frame(env.df)
  
  # Check if the dataframe contains the expected columns
  expected_cols <- c('pval', 'r2')
  
  # If the expected columns are not present, display an informative error message
  if (!all(expected_cols %in% colnames(env.df))) {
    cat("Dataframe", i, "column names are:", colnames(env.df), "\n")
    stop(paste("Dataframe", i, "does not contain the expected columns:", 
               paste(setdiff(expected_cols, colnames(env.df)), collapse = ", ")))
  }
  
  # Subset the dataframe to only the 'pval' and 'r2' columns
  env.df <- env.df[, c('pval', 'r2')]
  
  # Add missing rows with default values (pval = 1, r2 = 0)
  missing_rows <- setdiff(env.gen.all, rownames(env.df))
  if (length(missing_rows) > 0) {
    # Create a dataframe for missing rows
    missing_df <- data.frame(pval = rep(1, length(missing_rows)), 
                             r2 = rep(0, length(missing_rows)), 
                             row.names = missing_rows)
    
    # Combine with the existing dataframe
    env.df <- rbind(env.df, missing_df)
  }
  
  # Ensure the rownames are in the same order as 'env.gen.all'
  env.df <- env.df[env.gen.all, ]
  
  # Rename columns for each dataframe (e.g., pval_df1, r2_df1, etc.)
  colnames(env.df) <- c(paste0("pval_df", i), paste0("r2_df", i))
  
  # Add the dataframe columns to the merged dataframe
  env.gen.all.mrg <- cbind(env.gen.all.mrg, env.df)
}
env.gen.all.mrg; dim(env.gen.all.mrg)

# plot
#Reshape the data to long format
env_long <- env.gen.all.mrg %>%
  pivot_longer(cols = starts_with("r2_df"), 
               names_to = "source", 
               values_to = "r2") %>%
  mutate(source = sub("r2_df", "", source)) %>%
  pivot_longer(cols = starts_with("pval_df"), 
               names_to = "pval_source", 
               values_to = "pval") %>%
  mutate(pval_source = sub("pval_df", "", pval_source))

# Merge the two sources if they should be the same
env_long <- env_long %>%
  filter(source == pval_source) %>%
  select(-pval_source)

#Add kingdom and phylum
#Extract the taxonomic information from the phyloseq object
taxo_df <- as.data.frame(tax_table(ps))  # Convert tax_table to a DataFrame
taxo_df$rowname <- taxo_df$Genus  # Add rownames as a column
#Merge env_long with the taxonomic dataFrame 
env_long_taxo <- env_long %>%
  left_join(taxo_df[, c("rowname", "Kingdom", "Phylum")], by = "rowname")  
env_long_taxo

env_long_taxo_ord <- env_long_taxo %>%
  mutate(
    Kingdom = factor(Kingdom),  # Ensure Kingdom is a factor
    Phylum = factor(Phylum, levels = unique(Phylum[order(Kingdom, Phylum)]))  # Order Phylum by Kingdom, then alphabetically
  )

# # Assign colors to Phyla based on their Kingdom
# env_long_taxo_ord <- env_long_taxo_ord %>%
#   mutate(strip_color = phyla_colors[Phylum]) # Map colors to the strip_color column

# Determine the significance level based on p-values
env_long_taxo_ord$significance <- ifelse(env_long_taxo_ord$pval <= 0.001, "***", 
                                         ifelse(env_long_taxo_ord$pval <= 0.01, "**", 
                                                ifelse(env_long_taxo_ord$pval <= 0.05, "*", 
                                                       ifelse(env_long_taxo_ord$pval <= 0.1, ".", ""))))
env_long_taxo_ord


# FIGURE 4.B####
library(cowplot)
library(ggtext)
library(purrr)
#rename genera (correction):
env_long_taxo_ord$rowname = sub("_"," ", env_long_taxo_ord$rowname)
env_long_taxo_ord$rowname = sub("-"," ", env_long_taxo_ord$rowname)
#rename phyla
env_long_taxo_ord = env_long_taxo_ord %>%
  mutate(Phylum = recode(Phylum, 
                         "Actinobacteria" = "Actinomycetota", 
                         "Bacteroidetes" = "Bacteroidota", 
                         "Firmicutes" = "Bacillota",
                         "Proteobacteria" = "Pseudomonadota",
                         "Verrucomicrobia" = "Verrucomicrobiota")) 

g_bar = ggplot(env_long_taxo_ord, aes(x = r2, y = rowname, fill = source)) +  # Use 'source' for filling
  geom_text(aes(label=significance), hjust=0, vjust=1,color="black",angle = 0,
            position = position_dodge2(width = 0.7), size=7)+ #nudge_x= 0.02, nudge_y = c(0.3,0,-0.2)
  geom_bar(stat = "identity", width=0.7, position=position_dodge(0.7)) +  # Side-by-side bars
  labs(x = "R-squared", y = "Microbial Genera", fill = "Infants' Age") +
  #scale_y_discrete(position = "right") + #brings y-axis to the right
  theme_bw() +
  scale_fill_manual(values = c("#59A1A0","#c0a9cc","#5d8dc9"), labels = c("1-6 Months", "7-12 Months", "13-18 Months")) + 
  facet_grid(Phylum ~., scales = "free_y", space = "free_y") +
  background_grid() +
  theme(
    legend.title = element_text(size = 22, face = "bold"), 
    legend.text = element_text(size = 20), 
    axis.text.x = element_text(size = 14, colour = "black"),
    #axis.ticks.y.left = element_line(colour = "red",linewidth = 0.1),  # Add ticks to the right y-axis
    axis.text.y =  element_text(size = 18, face = "italic",vjust = 1,colour = "black"),
    axis.title = element_text(size = 20, face="bold"),
    strip.placement = "outside",  # Place strips outside the plot area
    strip.background = element_blank(),
    strip.text = element_textbox_highlight(
      size = 16,face = "bold", #bold.italic
      # unnamed set (all facet windows except named sets below)
      color = "black", fill = '#CCCCCC', box.color = '#CCCCCC',
      halign = 0.5, linetype = 1, r = unit(0, "pt"), width = unit(0.5, "npc"), 
      padding =margin(2, 0.5, 2, 0.5, "cm"), margin = unit(1, "pt"),
      # first named set
      hi.labels = c("Opisthokonta", "SAR"),
      hi.fill = "#eb8c44", hi.box.col = "#eb8c44", hi.col = "black"), #"#eb8c44"
    legend.position="top",
    panel.spacing.y = unit(0.5, "lines")); g_bar   # Adjust spacing between facets

library(patchwork)
g_bar +
  plot_annotation(tag_levels = list(c('B'), '1')) & 
  theme(plot.tag = element_text(size = 35, face = "bold"))

#save
save.image("~/Documents/xoxo_article/files/c1_xoxo_env.RData")
