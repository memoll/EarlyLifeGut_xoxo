#################################################################################
# Network analysis - NetCoMi                                                    #
# Data: Metagenomics and 18S (absolute abundance) - xoxo - birth                # 
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

#subset samples - birth ####
ps.vag = subset_samples(ps.all, sample_data(ps.all)$BIRTH == "VAGINAL")
ps.vag = prune_taxa(taxa_sums(ps.vag)>0, ps.vag); ps.vag
ps.cs = subset_samples(ps.all, sample_data(ps.all)$BIRTH == "Csec")
ps.cs = prune_taxa(taxa_sums(ps.cs)>0, ps.cs); ps.cs

#% vaginal ####
net_brt_vag = netConstruct(data = ps.vag, 
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
                           seed = 123456); net_brt_vag

# Network analysis ####
props_brt_vag <- netAnalyze(net_brt_vag, 
                            centrLCC = TRUE,
                            clustMethod = "cluster_fast_greedy",
                            #hubPar = c("degree", "betweenness", "closeness","eigenvector"),
                            hubQuant = 0.95,
                            weightDeg = TRUE, normDeg = FALSE)

summary(props_brt_vag)

props_brt_vag$hubs$hubs1
sort(props_brt_vag$centralities$degree1, decreasing = TRUE)
sort(props_brt_vag$centralities$close1, decreasing = TRUE)
sort(props_brt_vag$centralities$eigenv1, decreasing = TRUE)

#plot
# Create label vector
taxo.vag_net = tax_table(ps.vag)[props_brt_vag$lccNames1, ]
ps.vag_net = phyloseq(otu_table(ps.vag),tax_table(taxo.vag_net),sample_data(ps.vag))
labels_vag <- as.vector(tax_table(ps.vag_net)[, "Species"])
names(labels_vag) <- rownames(tax_table(ps.vag_net))
labels_vag <- gsub(" ", "\n", labels_vag) #change line after space
shapeVec_vag <- as.vector(tax_table(ps.vag_net)[, "Kingdom"])
names(shapeVec_vag) <- rownames(tax_table(ps.vag_net)[, "Kingdom"])
featVec_vag <- as.vector(tax_table(ps.vag_net)[, "Phylum"])
names(featVec_vag) <- rownames(tax_table(ps.vag_net)[, "Phylum"])

#asv
plot(props_brt_vag, 
     repulsion = 0.9, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#genus
plot(props_brt_vag, 
     repulsion = 0.9, labels = labels_vag, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#FIGURE 7.D####
unique(as.vector(tax_table(ps.vag_net)[, "Phylum"]))
length(unique(as.vector(tax_table(ps.vag_net)[, "Phylum"])))
p.net_brt_vag = plot(props_brt_vag, 
                     edgeInvisPar = 1, 
                     negDiffCol = TRUE, #if TRUE, show both pos and neg
                     posCol = "#6fc7c7", 
                     negCol = "#c4b3b9",
                     edgeTranspLow = 0,
                     edgeTranspHigh = 20,
                     nodeColor = "feature", 
                     featVecCol = featVec_vag, 
                     colorVec = c("#225ea8","#31a354","#88419d","#59A1A0","#eded02","#980043"),
                     nodeTransp = 20,
                     labels = labels_vag,
                     nodeSize = "normCounts",#"eigenvector" 
                     repulsion = 1.2, #Place the nodes further apart
                     rmSingles = TRUE, #removes singletons
                     labelScale = FALSE,#labelLength = 14, #labelScale = TRUE,  
                     cexLabels = 1, cexHubLabels = 1.2, nodeSizeSpread = 7, cexNodes = 5, 
                     usePCH = TRUE, #avoid rescaling (prevent circles of becoming ovals)
                     hubBorderCol = "mediumvioletred",borderWidth = 3,#borderWidth = 1,
                     labelFont = 3, hubLabelFont = 4, #1:plain,2:bold,3:italic,4:bold&italic
                     #highlightHubs = TRUE, hubTransp = 0.5*5, 
                     nodeShape = c("circle", "triangle"), featVecShape = shapeVec_vag, 
                     title1 = "Vaginal Birth\n",
                     mar = c(2,7,2,7), #margin
                     showTitle = TRUE, cexTitle = 2); p.net_brt_vag 

legend(x = -1.5, y = 1.3, cex = 3,bty = "n", y.intersp=0.5,x.intersp=0.5,
       legend = (bquote(bold("D"))), horiz = FALSE)

p.net_brt_vag$q1$Arguments$cut

#% C-section ####
net_brt_cs = netConstruct(data = ps.cs, 
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
                          seed = 123456); net_brt_cs

# Network analysis ####
props_brt_cs <- netAnalyze(net_brt_cs, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           #hubPar = c("degree", "betweenness", "closeness","eigenvector"),
                           hubQuant = 0.95,
                           weightDeg = TRUE, normDeg = FALSE)

summary(props_brt_cs)
props_brt_cs$hubs$hubs1
sort(props_brt_cs$centralities$degree1, decreasing = TRUE)
sort(props_brt_cs$centralities$close1, decreasing = TRUE)
sort(props_brt_cs$centralities$eigenv1, decreasing = TRUE)

#plot
# Create label vector
taxo.cs_net = tax_table(ps.cs)[props_brt_cs$lccNames1, ]
ps.cs_net = phyloseq(otu_table(ps.cs),tax_table(taxo.cs_net),sample_data(ps.cs))
labels_cs <- as.vector(tax_table(ps.cs_net)[, "Species"])
names(labels_cs) <- rownames(tax_table(ps.cs_net))
labels_cs <- gsub(" ", "\n", labels_cs) #change line after space
shapeVec_cs <- as.vector(tax_table(ps.cs_net)[, "Kingdom"])
names(shapeVec_cs) <- rownames(tax_table(ps.cs_net)[, "Kingdom"])
featVec_cs <- as.vector(tax_table(ps.cs_net)[, "Phylum"])
names(featVec_cs) <- rownames(tax_table(ps.cs_net)[, "Phylum"])

#asv
plot(props_brt_cs, 
     repulsion = 0.9, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#genus
plot(props_brt_cs, 
     repulsion = 0.9, labels = labels_cs, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#FIGURE 7.E#####
unique(as.vector(tax_table(ps.cs_net)[, "Phylum"]))
length(unique(as.vector(tax_table(ps.cs_net)[, "Phylum"])))
p.net_brt_cs = plot(props_brt_cs, 
                    edgeInvisPar = 1, 
                    negDiffCol = TRUE, #if TRUE, show both pos and neg
                    posCol = "#6fc7c7", 
                    negCol = "#c4b3b9",
                    edgeTranspLow = 0,
                    edgeTranspHigh = 20,
                    nodeColor = "feature", 
                    featVecCol = featVec_cs, 
                    colorVec = c("#225ea8","#dd1c77","#31a354","#59A1A0","#eded02"),
                    nodeTransp = 20,
                    labels = labels_cs,
                    nodeSize = "normCounts",#"eigenvector" 
                    repulsion = 1.2, #Place the nodes further apart
                    rmSingles = TRUE, #removes singletons
                    labelScale = FALSE,#labelLength = 14, #labelScale = TRUE,  
                    cexLabels = 1, cexHubLabels = 1.2, nodeSizeSpread = 7, cexNodes = 5, 
                    usePCH = TRUE, #avoid rescaling (prevent circles of becoming ovals)
                    hubBorderCol = "mediumvioletred",borderWidth = 3,#616666",borderWidth = 1,
                    labelFont = 3, hubLabelFont = 4, #1:plain,2:bold,3:italic,4:bold&italic
                    #highlightHubs = TRUE, hubTransp = 0.5*5, 
                    nodeShape = c("circle", "triangle"), featVecShape = shapeVec_cs, 
                    title1 = "C-section Birth\n",
                    mar = c(2,7,2,7), #margin
                    showTitle = TRUE, cexTitle = 2); p.net_brt_cs 

legend(x = -1.55, y = 1.2, cex = 3,bty = "n", y.intersp=0.5,x.intersp=0.5,
       legend = (bquote(bold("E"))), horiz = FALSE)

p.net_brt_cs$q1$Arguments$cut

#Compare Networks ####
#% vaginal vs Csec ####
net_brt = netConstruct(data = ps.vag, 
                       data2 = ps.cs,  
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
                       seed = 123456); net_brt 


# Network analysis ####
props_brt <- netAnalyze(net_brt, 
                        centrLCC = TRUE,
                        clustMethod = "cluster_fast_greedy",
                        hubQuant = 0.95,
                        weightDeg = TRUE, normDeg = FALSE)

summary(props_brt)

#plot
# Create label vector
taxo_net = tax_table(ps.cs)[props_brt$lccNames1, ]
ps_net = phyloseq(otu_table(ps.all),tax_table(taxo_net),sample_data(ps.all))
labels <- as.vector(tax_table(ps_net)[, "Species"])
names(labels) <- rownames(tax_table(ps_net))
labels <- gsub(" ", "\n", labels) #change line after space
shapeVec <- as.vector(tax_table(ps_net)[, "Kingdom"])
names(shapeVec) <- rownames(tax_table(ps_net)[, "Kingdom"])
featVec <- as.vector(tax_table(ps_net)[, "Phylum"])
names(featVec) <- rownames(tax_table(ps_net)[, "Phylum"])

#asv
plot(props_brt, 
     repulsion = 0.9, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     groupNames = c("VAGINAL", "Csec"),
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#genus
plot(props_brt, 
     repulsion = 0.9, labels = labels, sameLayout = TRUE, layoutGroup = "union", rmSingles = "inboth", shortenLabels = "simple",
     labelLength = 9, nodeSize = "mclr", labelScale = FALSE, cexNodes = 1.5, cexLabels = 0, cexHubLabels = 1, cexTitle = 1.2,
     groupNames = c("VAGINAL", "Csec"),nodeColor = "cluster", sameColThresh = 2,
     hubBorderCol  = "gray40", mar = c(2, 6, 4, 6))

#genus - eigenvector
#FIGURE S4.D####
#number of clusters (max number of clusters)
unique(props_brt$clusteringLCC$clust1) 
unique(props_brt$clusteringLCC$clust2)

## Number of clusters detected by Fast Greedy method
#num_clusters_fg <- length(unique(cl_fast_greedy$membership))
p.net_brt = plot(props_brt, 
                 edgeInvisPar = 1, 
                 negDiffCol = TRUE, #if TRUE, show both pos and neg
                 posCol = "#6fc7c7", 
                 negCol = "#c4b3b9",
                 edgeTranspLow = 0,
                 edgeTranspHigh = 20,
                 nodeColor = "cluster", 
                 colorVec = c("#e0a00b","#acb872","#59A1A0","#e04f0b","#bd7abf","#4e74a3"),,#mutual nodes between clusters:2
                 sameClustCol = TRUE,#clusters having at least sameColThresh nodes in common have the same color
                 sameColThresh = 2,#number of nodes a cluster must have in common in the two groups to have the same color
                 nodeTransp = 20,
                 labels = labels,
                 nodeSize = "normCounts",#"eigenvector" 
                 repulsion = 1.2, 
                 rmSingles = TRUE, #removes singletons
                 labelScale = FALSE,#labelLength = 14, #labelScale = TRUE,  
                 cexLabels = 1, cexHubLabels = 1.2, nodeSizeSpread = 7, cexNodes = 5, 
                 usePCH = TRUE, #avoid rescaling (prevent circles of becoming ovals)
                 hubBorderCol = "#484d4d",borderWidth = 1,
                 labelFont = 3, hubLabelFont = 4, 
                 #highlightHubs = TRUE, hubTransp = 0.5*5, 
                 nodeShape = c("circle", "triangle"), featVecShape = shapeVec, 
                 title1 = "Vaginal",
                 title2 = "C-section",
                 mar = c(2,7,2,7), #margin
                 showTitle = TRUE, cexTitle = 2); p.net_brt 

legend(x = -1.9, y = 1.3, cex = 3,bty = "n", y.intersp=0.5,x.intersp=0.5,
       legend = (bquote(bold("D"))), horiz = FALSE)

legend(-0.8, 1.1, cex = 1.2, title = "Estimated Correlation:", text.font = 2, y.intersp=0.7,x.intersp=0.4,
       legend = c("+","-"), lty = c(1,1), lwd = 1.5, col = c("#6fc7c7","#c4b3b9"),
       bty = "n", horiz = FALSE)
legend(x = -0.4, y = -0.4, pch=c(1,2), pt.cex = 2, cex = 1.2, bty = "n", y.intersp=0.7,x.intersp=0.5,#text.width=0.2,border = "black",
       legend = c("Eukaryotic Species","Bacterial Species"), horiz = FALSE)

p.net_brt$q1$Arguments$cut

#Compares network properties for microbial networks #### 
comp_brt <- netCompare(props_brt, 
                       adjust = "BH",        
                       permTest = TRUE, 
                       nPerm = 5000L, 
                       verbose = TRUE,
                       seed = 123456)
summary(comp_brt, 
        groupNames = c("VAGINAL", "Csec"))

# Differential network construction
diff_brt <- diffnet(net_brt,
                    diffMethod = "fisherTest", #for Fisher's z-test
                    discordThresh = 0.8, #a threshold for the posterior probability a pair of taxa is differentially correlated between the groups. 
                    #Taxa pairs with a posterior above this threshold are connected in the network.
                    adjust = "lfdr", #for local false discovery rate correction
                    cores = 8L,
                    seed = 123456)

#save ####
save.image("~/Documents/xoxo_article/files/c4_xoxo_NetCoMi_birth.RData")

