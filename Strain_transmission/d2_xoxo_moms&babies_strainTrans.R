###############################################################
# Strain Transmission - clade tree                            #
# Data: StrainPhlan4 - xoxo                                   # 
# Mona Parizadeh - May 2025                                   #
###############################################################

# load libraries 
library(ggplot2)
library(ggtree)
library(RColorBrewer)

setwd('~/Documents/xoxo_article/')

# Read in metadata file
meta <- read.delim("files/metadata_xoxo_mi_18months_mom.csv", header = T, sep = "," )
rownames(meta)

# Corrections - change month 19 to 18 (correction)
levels(as.factor(md$MONTH))
which(md$MONTH == "19")
md[52,5]
md[52,5] = as.numeric(sub("19", "18", as.character(md[52, 5])))
md[52,]

# Remove mock samples
md = md %>% filter(md$DATASET != "MOCK")

# Add age category
md$MONTH_GROUP <- cut(as.numeric(md$MONTH), breaks = c(-1,0, 6, 12, 18), labels = c("birth","1-6 Months", "7-12 Months", "13-18 Months"))

# Add a column with the same IDs for infants and their correspinding mom
md$type = md$INFANT_ID; md$type
which(md$type == "Mom2")
md$type[146] = stringr::str_replace(md$type[146], 'Mom2', 'XC02')
which(md$type == "Mom4")
md$type[147] = stringr::str_replace(md$type[147], 'Mom4', 'XC04')
which(md$type == "Mom5")
md$type[148] = stringr::str_replace(md$type[148], 'Mom5', 'XC05')
which(md$type == "Mom6")
md$type[149] = stringr::str_replace(md$type[149], 'Mom6', 'XC06')
which(md$type == "Mom7")
md$type[150] = stringr::str_replace(md$type[150], 'Mom7', 'XC07')
which(md$type == "Mom9")
md$type[151] = stringr::str_replace(md$type[151], 'Mom9', 'XC09')
which(md$type == "Mom11")
md$type[152] = stringr::str_replace(md$type[152], 'Mom11', 'XC11')
which(md$type == "Mom15")
md$type[153] = stringr::str_replace(md$type[153], 'Mom15', 'XC15')
which(md$type == "Mom19")
md$type[154] = stringr::str_replace(md$type[154], 'Mom19', 'XC19')
md$type

# For each clade marker from StrainPhlan analysis
# Example: t__SGB17248 Bifidobacterium longum (detected in 122 samples) ####

# Read newick tree
tre <- read.tree("05e_strainphlan4_tre/RAxML_bestTree.t__SGB17248.StrainPhlAn4.tre")
tre$tip.label

# Match tree tip labels and meta rownames
tre$tip.label = sub("_.*", "",tre$tip.label)
rownames(meta) = meta$MI_ID

library(tibble)
# Add a column in your metadata for the tip labels
meta <- meta %>% 
  rownames_to_column(var = "tip.label")

# Create ggtree object
gg <- ggtree( tre )
get_taxa_name(gg)

# Add metadata to dendrogram plot
gg <- gg %<+% meta

# Draw StrainPhlan tree
ggtree(tre, layout="circular") %<+% meta +
  geom_tippoint( size = 3, aes( color = INFANT_ID) ) +
  aes( branch.length = 'length' ) +
  theme_tree2() + 
  theme(legend.position="right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))


# Setting strain identity thresholds ####
# Check for strain transmission
# Read the tsv table of pairwise distances
nGD_17248 <- read_tsv(file = "05g_strainphlan4_MST/t__SGB17248_nGD.tsv", col_names = F, show_col_types = F)

# Make sample names identical to MI_ID in metadata
nGD_17248$X1 <- sub("_S.*", "", nGD_17248$X1)
nGD_17248$X2 <- sub("_S.*", "", nGD_17248$X2)

intersect(md$MI_ID,nGD_17248$X1)
intersect(md$MI_ID,nGD_17248$X2)

# Add metadata to the table
nGD_17248 <- left_join(nGD_17248 %>% select(MI_ID_1 = X1, everything()),
                       md %>% select(MI_ID_1 = MI_ID,
                                     DATASET_1 = DATASET,
                                     SAMPLE_ID_1 = SAMPLE_ID,
                                     INFANT_ID_1 = INFANT_ID,
                                     MONTH_1 = MONTH,
                                     SEX_1 = SEX,
                                     BIRTH_1 = BIRTH,
                                     BIRTH_DATE_1 = BIRTH_DATE,
                                     AGE_AT_SAMPLINGdays_1 = AGE_AT_SAMPLINGdays,
                                     SAMPLE_DATE_1 = SAMPLE_DATE,
                                     WEIGHTg_1 = WEIGHTg,
                                     HEIGHTcm_1 = HEIGHTcm,
                                     PARASITES_1 = PARASITES,
                                     VACCINES_1 = VACCINES,
                                     END_BF_1 = END_BF,
                                     MOM_EDUC_1 = MOM_EDUC,
                                     HOUSE_MATERIAL_1 = HOUSE_MATERIAL,
                                     NUM_ROOMS_1 = NUM_ROOMS,
                                     WATER_ACCESS_1 = WATER_ACCESS,
                                     WATER_TYPE_1 = WATER_TYPE,
                                     SERVICES_1 = SERVICES,
                                     ENERGY_1 = ENERGY,
                                     ANIMALS_1 = ANIMALS,
                                     NUM_PEOPLE_1 = NUM_PEOPLE,
                                     MONTH_GROUP_1 = MONTH_GROUP,
                                     type_1 = type))

nGD_17248 <- left_join(nGD_17248 %>% select(MI_ID_2 = X2, everything()),
                       md %>% select(MI_ID_2 = MI_ID,
                                     DATASET_2 = DATASET,
                                     SAMPLE_ID_2 = SAMPLE_ID,
                                     INFANT_ID_2 = INFANT_ID,
                                     MONTH_2 = MONTH,
                                     SEX_2 = SEX,
                                     BIRTH_2 = BIRTH,
                                     BIRTH_DATE_2 = BIRTH_DATE,
                                     AGE_AT_SAMPLINGdays_2 = AGE_AT_SAMPLINGdays,
                                     SAMPLE_DATE_2 = SAMPLE_DATE,
                                     WEIGHTg_2 = WEIGHTg,
                                     HEIGHTcm_2 = HEIGHTcm,
                                     PARASITES_2 = PARASITES,
                                     VACCINES_2 = VACCINES,
                                     END_BF_2 = END_BF,
                                     MOM_EDUC_2 = MOM_EDUC,
                                     HOUSE_MATERIAL_2 = HOUSE_MATERIAL,
                                     NUM_ROOMS_2 = NUM_ROOMS,
                                     WATER_ACCESS_2 = WATER_ACCESS,
                                     WATER_TYPE_2 = WATER_TYPE,
                                     SERVICES_2 = SERVICES,
                                     ENERGY_2 = ENERGY,
                                     ANIMALS_2 = ANIMALS,
                                     NUM_PEOPLE_2 = NUM_PEOPLE,
                                     MONTH_GROUP_2 = MONTH_GROUP,
                                     type_2 = type))

nGD_17248 = nGD_17248 %>%
  filter(!is.na(nGD_17248$DATASET_1)) %>%
  filter(!is.na(nGD_17248$DATASET_2))

# Compute time difference between sample (important for longitudinal samples)
nGD_17248$timeDiff <- abs(nGD_17248$MONTH_1 - nGD_17248$MONTH_2)

# Annotate pairs of samples. Are they related? Are they from the same individual?
nGD_17248$same_individual <- ifelse(nGD_17248$INFANT_ID_1 == nGD_17248$INFANT_ID_2, "same_individual", "different_individual")
nGD_17248$related <- ifelse(nGD_17248$type_1 == nGD_17248$type_2, "related", "unrelated")
nGD_17248$same_dataset <- ifelse(nGD_17248$DATASET_1 == nGD_17248$DATASET_2, "same_dtset", "different_dtset")

# Detect strain sharing events part in the Methods of Valles-Colomer et al. (2023)
# Keeping only the training data
nGD_17248_training <- nGD_17248 %>% filter(type_1 %in% c("XC02","XC04","XC05","XC06","XC07","XC09","XC11","XC15","XC19")) %>%
  filter(type_2 %in% c("XC02","XC04","XC05","XC06","XC07","XC09","XC11","XC15","XC19"))

nGD_17248_training <- rbind(nGD_17248_training %>% filter(same_individual == "same_individual") %>%
                              filter(timeDiff <= 15) %>%
                              group_by(type_1) %>% arrange(type_1, timeDiff) %>% slice_head(n = 1) %>% ungroup(),
                            nGD_17248_training %>% filter(same_individual != "same_individual",  related == "unrelated", same_dataset == "different_dtset" ) %>%
                              group_by(type_1, type_2) %>% slice_head(n = 1) %>% ungroup())

table(nGD_17248_training$same_individual)

# Repeat this for each marker clade
