#################################################################################
# Network analysis - NetCoMi                                                    #
# Data: Metagenomics and 18S (absolute abundance) - xoxo - age                  # 
# Mona Parizadeh - May 2025                                                     #
#################################################################################

#load libraries
library(phyloseq); packageVersion("phyloseq") 
library(DESeq2); packageVersion("DESeq2") 
library(ggplot2); packageVersion("ggplot2") 
library(SpiecEasi); packageVersion("SpiecEasi") 
library(metagMisc); packageVersion("metagMisc") 
library(dplyr); packageVersion("dplyr") 
library(SPRING); packageVersion("SPRING") 
library(NetCoMi); packageVersion("NetCoMi") 
library(limma); packageVersion("limma") 
library(LaplacesDemon)

# Import non-rarefied data #### 
setwd("~/Documents/xoxo_article/files/")
#Load data ####
ps = readRDS("ps_bac_euk_abs.rds"); ps

#Replace NAs ####
comm = otu_table(ps)
meta = sample_data(ps)
taxo = tax_table(ps)

#remove "_" and "-" from the Species names
taxo <- gsub("_", " ", taxo)
taxo <- gsub("-", " ", taxo)

#replace unassigned by the lowest  assigned taxonomic level
for(i in c(1:length(colnames(tax_table(taxo))))){
  for(j in rownames(tax_table(taxo))){
    if(is.na(tax_table(taxo)[j,i]) == TRUE){
      if(substr(tax_table(taxo)[j,i-1], start=1,stop=2)=="UA"){
        tax_table(taxo)[j,i] <- tax_table(taxo)[j,i-1]
      }
      else{
        tax_table(taxo)[j,i] <- paste0("UA_",tax_table(taxo)[j,i-1])
      }
    }
  }
}

#put them back together
ps.all = phyloseq(comm, taxo, meta)

#subset samples - age ####
ps.6m = subset_samples(ps.all, sample_data(ps.all)$AGE_GROUP == "1_6Months")
ps.6m = prune_taxa(taxa_sums(ps.6m)>0, ps.6m); ps.6m
ps.12m = subset_samples(ps.all, sample_data(ps.all)$AGE_GROUP == "7_12Months")
ps.12m = prune_taxa(taxa_sums(ps.12m)>0, ps.12m); ps.12m
ps.18m = subset_samples(ps.all, sample_data(ps.all)$AGE_GROUP == "13_18Months")
ps.18m = prune_taxa(taxa_sums(ps.18m)>0, ps.18m); ps.18m

#% 1-6 months ####
# Network construction ####
#signed transformation (transforms the estimated associations into dissimilarities)
#This leads to a network where strongly positive correlated taxa have a high edge weight (1 if the correlation equals 1) 
#and strongly negative correlated taxa have a low edge weight (0 if the correlation equals -1).

#unsigned transformation (the edge weight between strongly correlated taxa is high, no matter of the sign)
#so a correlation of -1 and 1 would lead to an edge weight of 1

net_age_6m = netConstruct(data = ps.6m, 
                          filtTax = "highestVar",
                          filtTaxPar = list(highestVar = 50),
                          filtSamp = "totalReads",
                          filtSampPar = list(totalReads = 0),
                          measure = "sparcc",
                          normMethod = "none", 
                          zeroMethod = "none",
                          sparsMethod = "none", 
                          dissFunc = "unsigned", #low distance between strongly associated taxa (positively as well as negatively)
                          thresh = 0.3,
                          #p.adjust.method = "fdr",  # Adjusted p-value method
                          alpha = 0.01,   # Change significance threshold to 0.01
                          verbose = 2,
                          cores = 20,
                          seed = 123456); net_age_6m

# Network analysis ####
props_age_6m <- netAnalyze(net_age_6m, 
                           centrLCC = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           #hubPar = c("degree", "betweenness", "closeness","eigenvector"),
                           hubQuant = 0.95,
                           weightDeg = TRUE, normDeg = FALSE)

summary(props_age_6m)
props_age_6m$hubs$hubs1
sort(props_age_6m$centralities$degree1, decreasing = TRUE)
sort(props_age_6m$centralities$close1, decreasing = TRUE)
sort(props_age_6m$centralities$eigenv1, decreasing = TRUE)

#plot
taxo.6m_net = tax_table(ps.6m)[props_age_6m$lccNames1, ]
ps.6m_net = phyloseq(otu_table(ps.6m),tax_table(taxo.6m_net),sample_data(ps.6m))
# Create label vector
labels_6m <- as.vector(tax_table(ps.6m_net)[, "Species"])
names(labels_6m) <- rownames(tax_table(ps.6m_net))
labels_6m <- gsub(" ", "\n", labels_6m) #change line after space
shapeVec_6m <- as.vector(tax_table(ps.6m_net)[, "Kingdom"])
names(shapeVec_6m) <- rownames(tax_table(ps.6m_net)[, "Kingdom"])
featVec_6m <- as.vector(tax_table(ps.6m_net)[, "Phylum"])
names(featVec_6m) <- rownames(tax_table(ps.6m_net)[,"Phylum"])

#asv
plot(props_age_6m, 
     repulsion = 0.9, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#genus
plot(props_age_6m, 
     repulsion = 0.9, labels = labels_6m, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#FIGURE 7.A####
unique(as.vector(tax_table(ps.6m_net)[, "Phylum"]))
length(unique(as.vector(tax_table(ps.6m_net)[, "Phylum"])))
p.net_age_6m = plot(props_age_6m, 
                    edgeInvisPar = 1, 
                    negDiffCol = TRUE, #if TRUE, show both pos and neg
                    posCol = "#6fc7c7", 
                    negCol = "#c4b3b9",
                    edgeTranspLow = 0,
                    edgeTranspHigh = 20,
                    nodeColor = "feature", 
                    featVecCol = featVec_6m, 
                    colorVec = c("#225ea8","#dd1c77","#31a354","#59A1A0","#eded02","#FFB90D","#980043"),
                    nodeTransp = 20,
                    labels = labels_6m,
                    nodeSize = "normCounts",#"eigenvector" 
                    repulsion = 1.2, #Place the nodes further apart
                    rmSingles = TRUE, #removes singletons
                    labelScale = FALSE,#labelLength = 14, #labelScale = TRUE, 
                    #cexLabels = 1.2, cexHubLabels = 1.5,
                    cexLabels = 1, cexHubLabels = 1.2, nodeSizeSpread = 7, cexNodes = 5,
                    usePCH = TRUE, #avoid rescaling (prevent circles of becoming ovals)
                    hubBorderCol = "mediumvioletred", borderWidth = 3,  
                    labelFont = 3, hubLabelFont = 4, #1:plain,2:bold,3:italic,4:bold&italic
                    #highlightHubs = TRUE, hubTransp = 0.5*5, 
                    nodeShape = c("circle", "triangle"), featVecShape = shapeVec_6m, 
                    title1 = "1-6 Months\n",
                    mar = c(2,7,2,7), #margin
                    showTitle = TRUE, cexTitle = 2); p.net_age_6m 

legend(x = -1.55, y = 1.2, cex = 3,bty = "n", y.intersp=0.5,x.intersp=0.5,
       legend = (bquote(bold("A"))), horiz = FALSE)

legend(0.6, 1.1, cex = 1.2, title = "Estimated Correlation:", text.font = 2, y.intersp=0.7,x.intersp=0.4,
       legend = c("+","-"), lty = c(1,1), lwd = 1.5, col = c("#6fc7c7","#c4b3b9"),
       bty = "n", horiz = FALSE)
legend(x = 0.9, y = 0.6, pch=c(1,2), pt.cex = 2, cex = 1.2, bty = "n", y.intersp=0.7,x.intersp=0.5,#text.width=0.2,border = "black",
       legend = c("Eukaryotic Species","Bacterial Species"), horiz = FALSE)

p.net_age_6m$q1$Arguments$cut

#% 7-12 months ####
net_age_12m = netConstruct(data = ps.12m, 
                           filtTax = "highestVar",
                           filtTaxPar = list(highestVar = 50),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 0),
                           measure = "sparcc",
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "unsigned", 
                           thresh = 0.3,
                           alpha = 0.01,
                           verbose = 2,
                           cores = 20,
                           seed = 123456); net_age_12m

# Network analysis ####
props_age_12m <- netAnalyze(net_age_12m, 
                            centrLCC = TRUE,
                            clustMethod = "cluster_fast_greedy",
                            #hubPar = c("degree", "betweenness", "closeness","eigenvector"),
                            hubQuant = 0.95,
                            weightDeg = TRUE, normDeg = FALSE)

summary(props_age_12m)
props_age_12m$hubs$hubs1
sort(props_age_12m$centralities$degree1, decreasing = TRUE)
sort(props_age_12m$centralities$close1, decreasing = TRUE)
sort(props_age_12m$centralities$eigenv1, decreasing = TRUE)

#plot
# Create label vector
taxo.12m_net = tax_table(ps.12m)[props_age_12m$lccNames1, ]
ps.12m_net = phyloseq(otu_table(ps.12m),tax_table(taxo.12m_net),sample_data(ps.12m))
labels_12m <- as.vector(tax_table(ps.12m_net)[, "Species"])
names(labels_12m) <- rownames(tax_table(ps.12m_net))
labels_12m <- gsub(" ", "\n", labels_12m) #change line after space
shapeVec_12m <- as.vector(tax_table(ps.12m_net)[, "Kingdom"])
names(shapeVec_12m) <- rownames(tax_table(ps.12m_net)[, "Kingdom"])
featVec_12m <- as.vector(tax_table(ps.12m_net)[, "Phylum"])
names(featVec_12m) <- rownames(tax_table(ps.12m_net)[, "Phylum"])

#asv
plot(props_age_12m, 
     repulsion = 0.9, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#genus
plot(props_age_12m, 
     repulsion = 0.9, labels = labels_12m, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#FIGURE 7.B####
unique(as.vector(tax_table(ps.12m_net)[, "Phylum"]))
length(unique(as.vector(tax_table(ps.12m_net)[, "Phylum"])))
p.net_age_12m = plot(props_age_12m, 
                     edgeInvisPar = 1, 
                     negDiffCol = TRUE, #if TRUE, show both pos and neg
                     posCol = "#6fc7c7",
                     negCol = "#c4b3b9",
                     edgeTranspLow = 0,
                     edgeTranspHigh = 20,
                     nodeColor = "feature", 
                     featVecCol = featVec_12m, 
                     colorVec = c("#225ea8","#dd1c77","#31a354","#88419d","#59A1A0","#eded02","#980043"),
                     nodeTransp = 20,
                     labels = labels_12m,
                     nodeSize = "normCounts",#"eigenvector" 
                     repulsion = 1.2, #Place the nodes further apart
                     rmSingles = TRUE, #removes singletons
                     labelScale = FALSE,#labelLength = 14, #labelScale = TRUE,  
                     cexLabels = 1, cexHubLabels = 1.2, nodeSizeSpread = 7, cexNodes = 5, 
                     usePCH = TRUE, #avoid rescaling (prevent circles of becoming ovals)
                     hubBorderCol = "mediumvioletred", borderWidth = 3,
                     labelFont = 3, hubLabelFont = 4, #1:plain,2:bold,3:italic,4:bold&italic
                     nodeShape = c("circle", "triangle"), featVecShape = shapeVec_12m, 
                     title1 = "7-12 Months\n",
                     mar = c(2,7,2,7), #margin
                     showTitle = TRUE, cexTitle = 2); p.net_age_12m 

legend(x = -1.4, y = 1.23, cex = 3,bty = "n", y.intersp=0.5,x.intersp=0.5,
       legend = (bquote(bold("B"))), horiz = FALSE)

p.net_age_12m$q1$Arguments$cut

#% 13-18 months ####
net_age_18m = netConstruct(data = ps.18m, 
                           filtTax = "highestVar",
                           filtTaxPar = list(highestVar = 50),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 0),
                           measure = "sparcc",
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "unsigned", 
                           thresh = 0.3,
                           alpha = 0.01,
                           verbose = 2,
                           cores = 20,
                           seed = 123456); net_age_18m

# Network analysis ####
props_age_18m <- netAnalyze(net_age_18m, 
                            centrLCC = TRUE,
                            clustMethod = "cluster_fast_greedy",
                            #hubPar = c("degree", "betweenness", "closeness","eigenvector"),
                            hubQuant = 0.95,
                            weightDeg = TRUE, normDeg = FALSE)

summary(props_age_18m)
props_age_18m$hubs$hubs1
sort(props_age_18m$centralities$degree1, decreasing = TRUE)
sort(props_age_18m$centralities$close1, decreasing = TRUE)
sort(props_age_18m$centralities$eigenv1, decreasing = TRUE)

#plot
# Create label vector
taxo.18m_net = tax_table(ps.18m)[props_age_18m$lccNames1, ]
ps.18m_net = phyloseq(otu_table(ps.18m),tax_table(taxo.18m_net),sample_data(ps.18m))
labels_18m <- as.vector(tax_table(ps.18m_net)[, "Species"])
names(labels_18m) <- rownames(tax_table(ps.18m_net))
labels_18m <- gsub(" ", "\n", labels_18m) #change line after space
shapeVec_18m <- as.vector(tax_table(ps.18m_net)[, "Kingdom"])
names(shapeVec_18m) <- rownames(tax_table(ps.18m_net)[, "Kingdom"])
featVec_18m <- as.vector(tax_table(ps.18m_net)[, "Phylum"])
names(featVec_18m) <- rownames(tax_table(ps.18m_net)[, "Phylum"])

#asv
plot(props_age_18m, 
     repulsion = 0.9, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#genus
plot(props_age_18m, 
     repulsion = 0.9, labels = labels_18m, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))


#FIGURE 7.C####
unique(as.vector(tax_table(ps.18m_net)[, "Phylum"]))
length(unique(as.vector(tax_table(ps.18m_net)[, "Phylum"])))
p.net_age_18m = plot(props_age_18m, 
                     edgeInvisPar = 1, 
                     negDiffCol = TRUE, #if TRUE, show both pos and neg
                     posCol = "#6fc7c7", 
                     negCol = "#c4b3b9",
                     edgeTranspLow = 0,
                     edgeTranspHigh = 20,
                     nodeColor = "feature", 
                     featVecCol = featVec_18m, 
                     colorVec = c("#225ea8","#dd1c77","#31a354","#eded02"),
                     nodeTransp = 20,
                     labels = labels_18m,
                     nodeSize = "normCounts",#"eigenvector" 
                     repulsion = 1.2, #Place the nodes further apart
                     rmSingles = TRUE, #removes singletons
                     labelScale = FALSE,#labelLength = 14, #labelScale = TRUE,  
                     cexLabels = 1, cexHubLabels = 1.2, nodeSizeSpread = 7, cexNodes = 5, 
                     usePCH = TRUE, #avoid rescaling (prevent circles of becoming ovals)
                     hubBorderCol = "mediumvioletred", ,borderWidth = 3,
                     labelFont = 3, hubLabelFont = 4, #1:plain,2:bold,3:italic,4:bold&italic
                     nodeShape = c("circle", "triangle"), featVecShape = shapeVec_18m, 
                     title1 = "13-18 Months\n",
                     mar = c(2,7,2,7), #margin
                     showTitle = TRUE, cexTitle = 2); p.net_age_18m 

legend(x = -1.55, y = 1.2, cex = 3,bty = "n", y.intersp=0.5,x.intersp=0.5,
       legend = (bquote(bold("C"))), horiz = FALSE)

legend(x = 0.9, y = 0.15, cex = 1.4,bty = "n", y.intersp=0.5,x.intersp=0.5,
       legend = (bquote(bold("Phyla"))), horiz = FALSE)
legend(x = 1, y = 0, cex = 1.2,bty = "n", y.intersp=0.35,x.intersp=0.5,
       #new names
       legend = c("Actinomycetota","Bacillota","Bacteroidota","Fusobacteriota",
                  "Opisthokonta","Pseudomonadota","SAR","Verrucomicrobiota"), horiz = FALSE,
       pch=c(16,16,16,16,17,16,17,16), pt.cex = 2.5,
       col=c("#225ea8","#31a354","#dd1c77","#88419d","#59A1A0","#eded02","#FFB90D","#980043"))

p.net_age_18m$q1$Arguments$cut

#Compare Networks ####
#% 6 vs 12 months ####
# Network construction ####
net_age_6vs12m = netConstruct(data = ps.6m, 
                              data2 = ps.12m,  
                              filtTax = "highestVar",
                              filtTaxPar = list(highestVar = 50),
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 0),
                              measure = "sparcc",
                              normMethod = "none", 
                              zeroMethod = "none",
                              sparsMethod = "none", 
                              dissFunc = "unsigned", 
                              thresh = 0.3,
                              verbose = 2,
                              cores = 20,
                              seed = 123456); net_age_6vs12m

# Network analysis ####
props_age_6vs12m <- netAnalyze(net_age_6vs12m, 
                               centrLCC = TRUE,
                               clustMethod = "cluster_fast_greedy",
                               hubQuant = 0.95,
                               weightDeg = TRUE, normDeg = FALSE)

summary(props_age_6vs12m)

#plot
# Create label vector
taxo.6vs12m_net = tax_table(ps.all)[props_age_6vs12m$lccNames1, ]
ps.6vs12m_net = phyloseq(otu_table(ps.all),tax_table(taxo.6vs12m_net),sample_data(ps.all))
labels_6vs12m <- as.vector(tax_table(ps.6vs12m_net)[, "Species"])
names(labels_6vs12m) <- rownames(tax_table(ps.6vs12m_net))
labels_6vs12m <- gsub(" ", "\n", labels_6vs12m) #change line after space
shapeVec_6vs12m <- as.vector(tax_table(ps.6vs12m_net)[, "Kingdom"])
names(shapeVec_6vs12m) <- rownames(tax_table(ps.6vs12m_net)[, "Kingdom"])
featVec_6vs12m <- as.vector(tax_table(ps.6vs12m_net)[, "Phylum"])
names(featVec_6vs12m) <- rownames(tax_table(ps.6vs12m_net)[, "Phylum"])

#asv
plot(props_age_6vs12m, 
     repulsion = 0.9, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     groupNames = c("1-6 Months", "7-12 Months"),
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#genus
plot(props_age_6vs12m, 
     repulsion = 0.9, labels = labels_6vs12m, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     groupNames = c("1-6 Months", "7-12 Months"),nodeColor = "cluster", sameColThresh = 2, #sameColThresh = 3, #to find out number of cluster
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#genus - eigenvector
#FIGURE S4.A####
# #number of clusters (max number of clusters)
# unique(props_age_6vs12m$clusteringLCC$clust1)
# unique(props_age_6vs12m$clusteringLCC$clust2)
p.net_age_6vs12m = plot(props_age_6vs12m, 
                        edgeInvisPar = 1, 
                        negDiffCol = TRUE, #if TRUE, show both pos and neg
                        posCol = "#6fc7c7", 
                        negCol = "#c4b3b9",
                        edgeTranspLow = 0,
                        edgeTranspHigh = 20,
                        nodeColor = "cluster", 
                        colorVec = c("#e0a00b","#59A1A0","#fc085d","#68477a","#225ea8","#980043"),#mutual nodes between clusters:2
                        sameClustCol = TRUE,#clusters having at least sameColThresh nodes in common have the same color
                        sameColThresh = 2,#number of nodes a cluster must have in common in the two groups to have the same color
                        nodeTransp = 20,
                        labels = labels_6vs12m,
                        nodeSize = "normCounts",#"eigenvector" 
                        repulsion = 1.2, 
                        rmSingles = TRUE, #removes singletons
                        labelScale = FALSE,#labelLength = 14, #labelScale = TRUE,  
                        cexLabels = 1, cexHubLabels = 1.2, nodeSizeSpread = 7, cexNodes = 5,
                        usePCH = TRUE, #avoid rescaling (prevent circles of becoming ovals)
                        hubBorderCol = "#484d4d",borderWidth = 1,
                        labelFont = 3, hubLabelFont = 4, 
                        #highlightHubs = TRUE, hubTransp = 0.5*5, 
                        nodeShape = c("circle", "triangle"), featVecShape = shapeVec_6vs12m, 
                        title1 = "1-6 Months",
                        title2 = "7-12 Months",
                        mar = c(2,3,2,3), #margin
                        showTitle = TRUE, cexTitle = 2); p.net_age_6vs12m 

legend(-0.6, 1.1, cex = 1.2, title = "Estimated Correlation:", text.font = 2, y.intersp=0.7,x.intersp=0.4,
       legend = c("+","-"), lty = c(1,1), lwd = 1.5, col = c("#6fc7c7","#c4b3b9"),
       bty = "n", horiz = FALSE)

legend(x = -1.5, y = 1.3, cex = 3,bty = "n", y.intersp=0.5,x.intersp=0.5,
       legend = (bquote(bold("A"))), horiz = FALSE)

p.net_age_6vs12m$q1$Arguments$cut

#Compares network properties for microbial networks #### 
comp_age_6vs12m <- netCompare(props_age_6vs12m, 
                              adjust = "BH",        
                              permTest = TRUE, 
                              nPerm = 5000L, 
                              verbose = TRUE,
                              cores = 20,
                              seed = 123456)

summary(comp_age_6vs12m, 
        groupNames = c("1-6 Months", "7-12 Months"))

# Differential network construction (identifying differentially associated taxa) ####
diff_age_6vs12m <- diffnet(net_age_6vs12m,
                           diffMethod = "fisherTest", #for Fisher's z-test
                           #diffMethod = "permute", nPerm = 1000L,
                           discordThresh = 0.8, #a threshold for the posterior probability a pair of taxa is differentially correlated between the groups. 
                           #Taxa pairs with a posterior above this threshold are connected in the network.
                           adjust = "lfdr", #for local false discovery rate correction
                           cores = 8L,
                           seed = 123456)

#% 12 vs 18 months ####
# Network construction ####
net_age_12vs18m = netConstruct(data = ps.12m, 
                               data2 = ps.18m,  
                               filtTax = "highestVar",
                               filtTaxPar = list(highestVar = 50),
                               filtSamp = "totalReads",
                               filtSampPar = list(totalReads = 0),
                               measure = "sparcc",
                               normMethod = "none", 
                               zeroMethod = "none",
                               sparsMethod = "none", 
                               dissFunc = "unsigned", 
                               thresh = 0.4,
                               verbose = 2,
                               cores = 20,
                               seed = 123456); net_age_12vs18m

# Network analysis ####
props_age_12vs18m <- netAnalyze(net_age_12vs18m, 
                                centrLCC = TRUE,
                                clustMethod = "cluster_fast_greedy",
                                hubQuant = 0.95,
                                weightDeg = TRUE, normDeg = FALSE)

summary(props_age_12vs18m)

#plot
# Create label vector
taxo.12vs18m_net = tax_table(ps.all)[props_age_12vs18m$lccNames1, ]
ps.12vs18m_net = phyloseq(otu_table(ps.all),tax_table(taxo.12vs18m_net),sample_data(ps.all))
labels_12vs18m <- as.vector(tax_table(ps.12vs18m_net)[, "Species"])
names(labels_12vs18m) <- rownames(tax_table(ps.12vs18m_net))
labels_12vs18m <- gsub(" ", "\n", labels_12vs18m) #change line after space
shapeVec_12vs18m <- as.vector(tax_table(ps.12vs18m_net)[, "Kingdom"])
names(shapeVec_12vs18m) <- rownames(tax_table(ps.12vs18m_net)[, "Kingdom"])
featVec_12vs18m <- as.vector(tax_table(ps.12vs18m_net)[, "Phylum"])
names(featVec_12vs18m) <- rownames(tax_table(ps.12vs18m_net)[, "Phylum"])

#asv
plot(props_age_12vs18m, 
     repulsion = 0.9, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     groupNames = c("7-12 Months", "13-18 Months"),
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#genus
plot(props_age_12vs18m, 
     repulsion = 0.9, labels = labels_12vs18m, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     groupNames = c("7-12 Months", "13-18 Months"),nodeColor = "cluster", sameColThresh = 2,
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#genus - eigenvector
#FIGURE S4.B####
p.net_age_12vs18m = plot(props_age_12vs18m, 
                         edgeInvisPar = 1, 
                         negDiffCol = TRUE, #if TRUE, show both pos and neg
                         posCol = "#6fc7c7", 
                         negCol = "#c4b3b9",
                         edgeTranspLow = 0,
                         edgeTranspHigh = 20,
                         nodeColor = "cluster", 
                         colorVec = c("#e0a00b","#59A1A0","#68477a","#225ea8","#fc085d"),#mutual nodes between clusters:2
                         sameClustCol = TRUE,#clusters having at least sameColThresh nodes in common have the same color
                         sameColThresh = 2,#number of nodes a cluster must have in common in the two groups to have the same color
                         nodeTransp = 20,
                         labels = labels_12vs18m,
                         nodeSize = "normCounts",#"eigenvector" 
                         repulsion = 1.2, 
                         rmSingles = TRUE, #removes singletons
                         labelScale = FALSE,#labelLength = 14, #labelScale = TRUE,  
                         cexLabels = 1, cexHubLabels = 1.2, nodeSizeSpread = 7, cexNodes = 5,
                         usePCH = TRUE, #avoid rescaling (prevent circles of becoming ovals)
                         hubBorderCol = "#484d4d",borderWidth = 1,
                         labelFont = 3, hubLabelFont = 4, 
                         #highlightHubs = TRUE, hubTransp = 0.5*5, 
                         nodeShape = c("circle", "triangle"), featVecShape = shapeVec_12vs18m, 
                         title1 = "7-12 Months",
                         title2 = "13-18 Months",
                         mar = c(2,3,2,3), #margin
                         showTitle = TRUE, cexTitle = 2); p.net_age_12vs18m 

legend(x = -1.5, y = 1.2, cex = 3,bty = "n", y.intersp=0.5,x.intersp=0.5,
       legend = (bquote(bold("B"))), horiz = FALSE)

p.net_age_12vs18m$q1$Arguments$cut

#Compares network properties for microbial networks #### 
comp_age_12vs18m <- netCompare(props_age_12vs18m, 
                               adjust = "BH",        
                               permTest = TRUE, 
                               nPerm = 5000L, 
                               verbose = TRUE,
                               seed = 123456)

summary(comp_age_12vs18m, 
        groupNames = c("7-12 Months", "13-18 Months"))

# Differential network construction (identifying differentially associated taxa) ####
diff_age_12vs18m <- diffnet(net_age_12vs18m,
                            diffMethod = "fisherTest", #for Fisher's z-test
                            discordThresh = 0.8, #a threshold for the posterior probability a pair of taxa is differentially correlated between the groups. 
                            #Taxa pairs with a posterior above this threshold are connected in the network.
                            adjust = "lfdr", #for local false discovery rate correction
                            cores = 8L,
                            seed = 123456)

#% 6 vs 18 months ####
# Network construction ####
net_age_6vs18m = netConstruct(data = ps.6m, 
                              data2 = ps.18m,  
                              filtTax = "highestVar",
                              filtTaxPar = list(highestVar = 50),
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 0),
                              measure = "sparcc",
                              normMethod = "none", 
                              zeroMethod = "none",
                              sparsMethod = "none", 
                              dissFunc = "unsigned", 
                              thresh = 0.4,
                              verbose = 2,
                              cores = 20,
                              seed = 123456); net_age_6vs18m

# Network analysis ####
props_age_6vs18m <- netAnalyze(net_age_6vs18m, 
                               centrLCC = TRUE,
                               clustMethod = "cluster_fast_greedy",
                               hubQuant = 0.95,
                               weightDeg = TRUE, normDeg = FALSE)

summary(props_age_6vs18m)

#plot
# Create label vector
taxo.6vs18m_net = tax_table(ps.all)[props_age_6vs18m$lccNames1, ]
ps.6vs18m_net = phyloseq(otu_table(ps.all),tax_table(taxo.6vs18m_net),sample_data(ps.all))
labels_6vs18m <- as.vector(tax_table(ps.6vs18m_net)[, "Species"])
names(labels_6vs18m) <- rownames(tax_table(ps.6vs18m_net))
labels_6vs18m <- gsub(" ", "\n", labels_6vs18m) #change line after space
shapeVec_6vs18m <- as.vector(tax_table(ps.6vs18m_net)[, "Kingdom"])
names(shapeVec_6vs18m) <- rownames(tax_table(ps.6vs18m_net)[, "Kingdom"])
featVec_6vs18m <- as.vector(tax_table(ps.6vs18m_net)[, "Phylum"])
names(featVec_6vs18m) <- rownames(tax_table(ps.6vs18m_net)[, "Phylum"])

#asv
plot(props_age_6vs18m, 
     repulsion = 0.9, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     groupNames = c("1-6 Months", "13-18 Months"),
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#genus
plot(props_age_6vs18m, 
     repulsion = 0.9, labels = labels_6vs18m, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     groupNames = c("1-6 Months", "13-18 Months"), nodeColor = "cluster", sameColThresh = 2,
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#genus - eigenvector
#FIGURE S4.C####
length(levels(as.factor(tax_table(ps.all)[, "Phylum"]))) #no. of colors needed
p.net_age_6vs18m = plot(props_age_6vs18m, 
                        edgeInvisPar = 1, 
                        negDiffCol = TRUE, #if TRUE, show both pos and neg
                        posCol = "#6fc7c7", 
                        negCol = "#c4b3b9",
                        edgeTranspLow = 0,
                        edgeTranspHigh = 20,
                        nodeColor = "cluster", 
                        colorVec = c("#e0a00b","#59A1A0","#68477a","#225ea8","#fc085d","#980043"),#mutual nodes between clusters:2
                        sameClustCol = TRUE,#clusters having at least sameColThresh nodes in common have the same color
                        sameColThresh = 2,#number of nodes a cluster must have in common in the two groups to have the same color
                        nodeTransp = 20,
                        labels = labels_6vs18m,
                        nodeSize = "normCounts",#"eigenvector" 
                        repulsion = 1.2, 
                        rmSingles = TRUE, #removes singletons
                        labelScale = FALSE,#labelLength = 14, #labelScale = TRUE,  
                        cexLabels = 1, cexHubLabels = 1.2, nodeSizeSpread = 7, cexNodes = 5,
                        usePCH = TRUE, #avoid rescaling (prevent circles of becoming ovals)
                        hubBorderCol = "#484d4d",borderWidth = 1,
                        labelFont = 3, hubLabelFont = 4, 
                        #highlightHubs = TRUE, hubTransp = 0.5*5, 
                        nodeShape = c("circle", "triangle"), featVecShape = shapeVec_6vs18m, 
                        title1 = "1-6 Months",
                        title2 = "13-18 Months",
                        mar = c(2,3,2,3), #margin
                        showTitle = TRUE, cexTitle = 2); p.net_age_6vs18m 

legend(x = -1.4, y = 1.4, cex = 3,bty = "n", y.intersp=0.5,x.intersp=0.5,
       legend = (bquote(bold("C"))), horiz = FALSE)

p.net_age_6vs18m$q1$Arguments$cut

#Compares network properties for microbial networks #### 
comp_age_6vs18m <- netCompare(props_age_6vs18m, 
                              adjust = "BH",  #adjust for multiple testing      
                              permTest = TRUE, 
                              nPerm = 5000L, #try with 1000 perm (defaul)
                              verbose = TRUE,
                              seed = 123456)

summary(comp_age_6vs18m, 
        groupNames = c("1-6 Months", "13-18 Months"))


# Differential network construction (identifying differentially associated taxa) ####
diff_age_6vs18m <- diffnet(net_age_6vs18m,
                           diffMethod = "fisherTest", #for Fisher's z-test
                           discordThresh = 0.8, #a threshold for the posterior probability a pair of taxa is differentially correlated between the groups. 
                           #Taxa pairs with a posterior above this threshold are connected in the network.
                           adjust = "lfdr", #for local false discovery rate correction
                           cores = 8L,
                           seed = 123456)


#save ####
save.image("~/Documents/xoxo_article/files/c3_xoxo_NetCoMi_age.RData")

