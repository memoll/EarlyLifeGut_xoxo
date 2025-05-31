###############################################################
# Calculating Absolute Abundance                              #
# Data: Bacteria and Eukaryotes - xoxo                        # 
# Mona Parizadeh - May 2025                                   #
###############################################################

library(phyloseq); packageVersion("phyloseq") 

# Import data #### 
setwd("~/") 

#Load data ####
ps.bac = readRDS("~/Documents/xoxo_article/files/metaphlan4/ps_xoxo_metagenome_cleaned.rds"); ps.bac 
ps.euk = readRDS("~/Documents/xoxo_article/files/18s/ps_xoxo_18s_VST.rds"); ps.euk

#Keep mutual samples ####
# keep only the samples in common between both becatria and eukaryotes datasets
ps.bac2 = subset_samples(ps.bac, sample_names(ps.bac) %in% sample_names(ps.euk))
ps.bac2 = prune_taxa(taxa_sums(ps.bac2)>0, ps.bac2); ps.bac2
taxa_names(ps.bac2) = paste0("bac_", taxa_names(ps.bac2));taxa_names(ps.bac2)
ps.euk2 = subset_samples(ps.euk, sample_names(ps.euk) %in% sample_names(ps.bac))
ps.euk2 = prune_taxa(taxa_sums(ps.euk2)>0, ps.euk2); ps.euk2
taxa_names(ps.euk2) = paste0("euk_", taxa_names(ps.euk2));taxa_names(ps.euk2)

#Absolute abundance ####
#bacteria
# Normalize ASV data by combining the ASV relative abundances with the ng / g of stool data
ps.bac2.ra = ps.bac2 %>% transform_sample_counts(function(x) {x/sum(x)} )
otu.bac2.ra = otu_table(ps.bac2.ra)
tax.bac2.ra = tax_table(ps.bac2.ra)
meta.bac2.ra = sample_data(ps.bac2.ra)
# match the order of bac_QUANTITY_MEAN to the column order of otu.bac2.ra
bac2.ra.match = match(colnames(otu.bac2.ra), rownames(meta.bac2.ra))
# extract the 'bac_QUANTITY_MEAN' in the correct order
bac2_quantity_mean = as.numeric(meta.bac2.ra$bac_QUANTITY_MEAN[bac2.ra.match])
# multiply each column of otu.bac2.ra by the corresponding value in bac2_quantity_mean
otu.bac2.abs <- sweep(otu.bac2.ra, 2, bac2_quantity_mean, `*`); dim(otu.bac2.abs)
#re-built the phyloseq object with absolute abundance
ps.bac2.abs = phyloseq(otu_table(otu.bac2.abs, taxa_are_rows = TRUE), tax_table(ps.bac2.ra), sample_data(meta.bac2.ra)); ps.bac2.abs

#eukaryotes
# Normalize ASV data by combining the ASV relative abundances with the ng / g of stool data
ps.euk2.ra = ps.euk2 %>% transform_sample_counts(function(x) {x/sum(x)} )
otu.euk2.ra = otu_table(ps.euk2.ra)
tax.euk2.ra = tax_table(ps.euk2.ra)
meta.euk2.ra = sample_data(ps.euk2.ra)
# match the order of euk_QUANTITY_MEAN to the column order of otu.bac2.ra
euk2.ra.match = match(rownames(otu.euk2.ra), rownames(meta.euk2.ra))
# extract the 'euk_QUANTITY_MEAN' in the correct order
euk2_quantity_mean = as.numeric(meta.euk2.ra$euk_QUANTITY_MEAN[euk2.ra.match])
# multiply each column of otu.euk2.ra by the corresponding value in euk2_quantity_mean
otu.euk2.abs <- sweep(otu.euk2.ra, 1, euk2_quantity_mean, `*`); dim(otu.euk2.abs)
#re-built the phyloseq object with absolute abundance
ps.euk2.abs = phyloseq(otu_table(otu.euk2.abs, taxa_are_rows = FALSE), tax_table(ps.euk2.ra), sample_data(meta.euk2.ra)); ps.euk2.abs

#merge
otu.abs = rbind(otu.bac2.abs,t(otu.euk2.abs)); dim(otu.abs)
# Remove columns that contain any NA values
otu.abs_noNA <- otu.abs[, colSums(is.na(otu.abs)) == 0]
tax.abs = rbind(tax.bac2.ra,tax.euk2.ra); dim(tax.abs) #the taxa table the same as the one from relative abundance table
identical(colnames(otu.abs),rownames(meta.bac2.ra)) #TRUE
identical(colnames(otu.abs),rownames(meta.euk2.ra)) #TRUE

ps.all.abs = phyloseq(otu_table(otu.abs_noNA, taxa_are_rows = TRUE), sample_data(meta.euk2.ra), 
                      tax_table(tax.abs))
ps.all.abs

#correct age categories ####
sample_data(ps.all.abs)$AGE_GROUP <- cut(as.numeric(sample_data(ps.all.abs)$MONTH), breaks = c(0, 6, 12, 18), labels = c("1_6Months", "7_12Months", "13_18Months"))

saveRDS(ps.all.abs, "~/Documents/xoxo_article/files/ps_bac_euk_abs.rds")

#save####
save.image("~/Documents/xoxo_article/files/absAbun_xoxo_bac_euk.RData")
