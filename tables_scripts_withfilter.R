#Load phyloseq
library("phyloseq")

# Créer donnés avec 1e-5 data Filtrer comme c’est 1e-5 et que j’ai 1484 OTU il faut calculer la moyenne de total counts pour savoir la probabilité que la motu se trouve seulement une fois

periobaseline_physeq_norm_filter<- filter_taxa(periobaseline_physeq_norm, function(x) mean(x) > 1.225, TRUE)

#Export relative abundance 1.PHYLUM tables of patients with 1.shallow, 2.intermediate and 3.deep sites;
#Export relative abundance 2.GENERA tables for all 1.cases, 2.all controls, 3.all cases shallow, 4.intermediate and 5.deep sites.


#Already normalized object to use: periobaseline_physeq_normf

#Create objects to separate patients from controls
#Create object with patients only
periobaseline_physeq_norm_filter_onlypatients <- subset_samples(periobaseline_physeq_norm_filter, abs_status!="C")

#Create objects with controls only
periobaseline_physeq_norm_filter_onlycontrols <- subset_samples(periobaseline_physeq_norm_filter, abs_status=="C")

##Create tables for only different depth sites for PATIENTS

periobaseline_physeq_norm_filter_onlypatients_shallow <- subset_samples(periobaseline_physeq_norm_filter_onlypatients, Depth=="S")

periobaseline_physeq_norm_filter_onlypatients_intermediate <- subset_samples(periobaseline_physeq_norm_filter_onlypatients, Depth=="I")

periobaseline_physeq_norm_filter_onlypatients_deep <- subset_samples(periobaseline_physeq_norm_filter_onlypatients, Depth=="D")


#1.1 Create PHYLUM table for controls

#Add all OTUs belonging to the same Phylum
periobaseline_physeq_norm_filter_onlycontrols_phylum= tax_glom(periobaseline_physeq_norm_filter_onlycontrols, "Phylum")

#Extract data frame containing OTU counts
periobaseline_physeq_norm_filter_onlycontrols_phylum_otutable <- as.data.frame(otu_table(periobaseline_physeq_norm_filter_onlycontrols_phylum))

#Rename rows
rownames(periobaseline_physeq_norm_filter_onlycontrols_phylum_otutable) <- as.data.frame(tax_table(periobaseline_physeq_norm_filter_onlycontrols_phylum))$Phylum

#Round table to integer values
periobaseline_physeq_norm_filter_onlycontrols_phylum_otutable  <- round(periobaseline_physeq_norm_filter_onlycontrols_phylum_otutable, 0)

#Add rows to get the total number of counts across samples
periobaseline_physeq_norm_filter_onlycontrols_totals <- as.data.frame(rowSums(periobaseline_physeq_norm_filter_onlycontrols_phylum_otutable))

#Transform to percentages
periobaseline_physeq_norm_filter_onlycontrols_totals <- (periobaseline_physeq_norm_filter_onlycontrols_totals /sum(periobaseline_physeq_norm_filter_onlycontrols_totals))*100

#Round off to 2 decimal places
periobaseline_physeq_norm_filter_onlycontrols_totals <- round(periobaseline_physeq_norm_filter_onlycontrols_totals, 2)

#Rename Columns_
colnames(periobaseline_physeq_norm_filter_onlycontrols_totals) <- "Percentage"

#Export
write.table(periobaseline_physeq_norm_filter_onlycontrols_totals , file="periobaseline_physeq_norm_filter_onlycontrols_totals.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")


#1.1 Create PHYLUM table for SHALLOW patients

#Add all OTUs belonging to the same Phylum
periobaseline_physeq_norm_filter_onlypatients_shallow_phylum = tax_glom(periobaseline_physeq_norm_filter_onlypatients_shallow, "Phylum")

#Extract data frame containing OTU counts
periobaseline_physeq_norm_filter_onlypatients_shallow_phylum_otutable <- as.data.frame(otu_table(periobaseline_physeq_norm_filter_onlypatients_shallow_phylum))

#Rename rows
rownames(periobaseline_physeq_norm_filter_onlypatients_shallow_phylum_otutable) <- as.data.frame(tax_table(periobaseline_physeq_norm_filter_onlypatients_shallow_phylum))$Phylum

#Round table to integer values
periobaseline_physeq_norm_filter_onlypatients_shallow_phylum_otutable <- round(periobaseline_physeq_norm_filter_onlypatients_shallow_phylum_otutable, 0)

#Add rows to get the total number of counts across samples
periobaseline_physeq_norm_filter_onlypatients_shallow_totals <- as.data.frame(rowSums(periobaseline_physeq_norm_filter_onlypatients_shallow_phylum_otutable))

#Transform to percentages
periobaseline_physeq_norm_filter_onlypatients_shallow_totals <- (periobaseline_physeq_norm_filter_onlypatients_shallow_totals/sum(periobaseline_physeq_norm_filter_onlypatients_shallow_totals))*100

#Round off to 2 decimal places
periobaseline_physeq_norm_filter_onlypatients_shallow_totals <- round(periobaseline_physeq_norm_filter_onlypatients_shallow_totals, 2)

#Rename Columns_
colnames(periobaseline_physeq_norm_filter_onlypatients_shallow_totals) <- "Percentage"

#Export
write.table(periobaseline_physeq_norm_filter_onlypatients_shallow_totals, file="periobaseline_physeq_norm_filter_onlypatients_shallow_totals.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

#1.1 Create PHYLUM table for INTERMEDIATE patients

#Add all OTUs belonging to the same Phylum
periobaseline_physeq_norm_filter_onlypatients_intermediate_phylum = tax_glom(periobaseline_physeq_norm_filter_onlypatients_intermediate, "Phylum")

#Extract data frame containing OTU counts
periobaseline_physeq_norm_filter_onlypatients_intermediate_phylum_otutable<- as.data.frame(otu_table(periobaseline_physeq_norm_filter_onlypatients_intermediate_phylum))

#Rename rows
rownames(periobaseline_physeq_norm_filter_onlypatients_intermediate_phylum_otutable)<- as.data.frame(tax_table(periobaseline_physeq_norm_filter_onlypatients_intermediate_phylum))$Phylum

#Round table to integer values
periobaseline_physeq_norm_filter_onlypatients_intermediate_phylum_otutable<- round(periobaseline_physeq_norm_filter_onlypatients_intermediate_phylum_otutable, 0)

#Add rows to get the total number of counts across samples
periobaseline_physeq_norm_filter_onlypatients_intermediate_totals<- as.data.frame(rowSums(periobaseline_physeq_norm_filter_onlypatients_intermediate_phylum_otutable))

#Transform to percentages
periobaseline_physeq_norm_filter_onlypatients_intermediate_totals<- (periobaseline_physeq_norm_filter_onlypatients_intermediate_totals/sum(periobaseline_physeq_norm_filter_onlypatients_intermediate_totals))*100

#Round off to 2 decimal places
periobaseline_physeq_norm_filter_onlypatients_intermediate_totals <- round(periobaseline_physeq_norm_filter_onlypatients_intermediate_totals, 2)

#Rename Columns
colnames(periobaseline_physeq_norm_filter_onlypatients_intermediate_totals) <- "Percentage"

#Export
write.table(periobaseline_physeq_norm_filter_onlypatients_intermediate_totals, file="periobaseline_physeq_norm_filter_onlypatients_intermediate_totals.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")


#1.1 Create PHYLUM table for Deep patients

#Add all OTUs belonging to the same Phylum
periobaseline_physeq_norm_filter_onlypatients_deep_phylum = tax_glom(periobaseline_physeq_norm_filter_onlypatients_deep, "Phylum")

#Extract data frame containing OTU counts
periobaseline_physeq_norm_filter_patients_deep_phylum_otutable <- as.data.frame(otu_table(periobaseline_physeq_norm_filter_onlypatients_deep_phylum))

#Rename rows
rownames(periobaseline_physeq_norm_filter_patients_deep_phylum_otutable) <- as.data.frame(tax_table(periobaseline_physeq_norm_filter_onlypatients_deep_phylum))$Phylum

#Round table to integer values
periobaseline_physeq_norm_filter_onlypatients_deep_phylum_otutable <- round(periobaseline_physeq_norm_filter_patients_deep_phylum_otutable, 0)

#Add rows to get the total number of counts across samples
periobaseline_physeq_norm_filter_onlypatients_deep_totals<- as.data.frame(rowSums(periobaseline_physeq_norm_filter_onlypatients_deep_phylum_otutable))

#Transform to percentages
periobaseline_physeq_norm_filter_onlypatients_deep_totals<- (periobaseline_physeq_norm_filter_onlypatients_deep_totals/sum(periobaseline_physeq_norm_filter_onlypatients_deep_totals))*100

#Round off to 2 decimal places
periobaseline_physeq_norm_filter_onlypatients_deep_totals <- round(periobaseline_physeq_norm_filter_onlypatients_deep_totals, 2)

#Rename Columns
colnames(periobaseline_physeq_norm_filter_onlypatients_deep_totals) <- "Percentage"

#Export
write.table(periobaseline_physeq_norm_filter_onlypatients_deep_totals, file="periobaseline_physeq_norm_filter_onlypatients_deep_totals.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")



#1.1 Create GENUS table for SHALLOW patients#Add all OTUs belonging to the same Genus

periobaseline_physeq_norm_filter_onlypatients_shallow_genus = tax_glom(periobaseline_physeq_norm_filter_onlypatients_shallow, "Genus")

#Extract data frame containing OTU counts
periobaseline_physeq_norm_filter_onlypatients_shallow_genus_otutable <- as.data.frame(otu_table(periobaseline_physeq_norm_filter_onlypatients_shallow_genus))

#Rename rows

genusnames <- as.data.frame(tax_table(periobaseline_physeq_norm_filter_onlypatients_shallow_genus))$Genus

make.unique(as.character(genusnames))

rownames(periobaseline_physeq_norm_filter_onlypatients_shallow_genus_otutable) <- make.unique(as.character(genusnames))


#Round table to integer values
periobaseline_physeq_norm_filter_onlypatients_shallow_genus_otutable <- round(periobaseline_physeq_norm_filter_onlypatients_shallow_genus_otutable, 0)

#Add rows to get the total number of counts across samples
periobaseline_physeq_norm_filter_onlypatients_shallow_totals <- as.data.frame(rowSums(periobaseline_physeq_norm_filter_onlypatients_shallow_genus_otutable))

#Transform to percentages
periobaseline_physeq_norm_filter_onlypatients_shallow_totals <- (periobaseline_physeq_norm_filter_onlypatients_shallow_totals/sum(periobaseline_physeq_norm_filter_onlypatients_shallow_totals))*100

#Round off to 2 decimal places
periobaseline_physeq_norm_filter_onlypatients_shallow_totals <- round(periobaseline_physeq_norm_filter_onlypatients_shallow_totals, 2)

#Rename Columns
colnames(periobaseline_physeq_norm_filter_onlypatients_shallow_totals) <- "Percentage"

#Export
write.table(periobaseline_physeq_norm_filter_onlypatients_shallow_totals, file="periobaseline_physeq_norm_filter_onlypatients_shallow_totals.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

#1.1 Create GENUS table for INTERMEDIATE patients#Add all OTUs belonging to the same Genus
periobaseline_physeq_norm_filter_onlypatients_intermediate_genus = tax_glom(periobaseline_physeq_norm_filter_onlypatients_intermediate, "Genus")

#Extract data frame containing OTU counts
periobaseline_physeq_norm_filter_onlypatients_intermediate_genus_otutable <- as.data.frame(otu_table(periobaseline_physeq_norm_filter_onlypatients_intermediate_genus))

#Rename rows

genusnames <- as.data.frame(tax_table(periobaseline_physeq_norm_filter_onlypatients_intermediate_genus))$Genus

make.unique(as.character(genusnames))

rownames(periobaseline_physeq_norm_filter_onlypatients_intermediate_genus_otutable) <- make.unique(as.character(genusnames))


#Round table to integer values
periobaseline_physeq_norm_filter_onlypatients_intermediate_genus_otutable <- round(periobaseline_physeq_norm_filter_onlypatients_intermediate_genus_otutable, 0)

#Add rows to get the total number of counts across samples
periobaseline_physeq_norm_filter_onlypatients_intermediate_totals <- as.data.frame(rowSums(periobaseline_physeq_norm_filter_onlypatients_intermediate_genus_otutable))

#Transform to percentages
periobaseline_physeq_norm_filter_onlypatients_intermediate_totals <- (periobaseline_physeq_norm_filter_onlypatients_intermediate_totals/sum(periobaseline_physeq_norm_filter_onlypatients_intermediate_totals))*100

#Round off to 2 decimal places
periobaseline_physeq_norm_filter_onlypatients_intermediate_totals <- round(periobaseline_physeq_norm_filter_onlypatients_intermediate_totals, 2)

#Rename Columns
colnames(periobaseline_physeq_norm_filter_onlypatients_intermediate_totals) <- "Percentage"

#Export
write.table(periobaseline_physeq_norm_filter_onlypatients_intermediate_totals, file="periobaseline_physeq_norm_filter_onlypatients_intermediate_totals.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

#1.1 Create GENUS table for DEEP patients#Add all OTUs belonging to the same Genus
periobaseline_physeq_norm_filter_onlypatients_deep_genus = tax_glom(periobaseline_physeq_norm_filter_onlypatients_deep, "Genus")

#Extract data frame containing OTU counts
periobaseline_physeq_norm_filter_onlypatients_deep_genus_otutable <- as.data.frame(otu_table(periobaseline_physeq_norm_filter_onlypatients_deep_genus))

#Rename rows

genusnames <- as.data.frame(tax_table(periobaseline_physeq_norm_filter_onlypatients_deep_genus))$Genus

make.unique(as.character(genusnames))

rownames(periobaseline_physeq_norm_filter_onlypatients_deep_genus_otutable) <- make.unique(as.character(genusnames))


#Round table to integer values
periobaseline_physeq_norm_filter_onlypatients_deep_genus_otutable <- round(periobaseline_physeq_norm_filter_onlypatients_deep_genus_otutable, 0)

#Add rows to get the total number of counts across samples
periobaseline_physeq_norm_filter_onlypatients_deep_totals <- as.data.frame(rowSums(periobaseline_physeq_norm_filter_onlypatients_deep_genus_otutable))

#Transform to percentages
periobaseline_physeq_norm_filter_onlypatients_deep_totals <- (periobaseline_physeq_norm_filter_onlypatients_deep_totals/sum(periobaseline_physeq_norm_filter_onlypatients_deep_totals))*100

#Round off to 2 decimal places
periobaseline_physeq_norm_filter_onlypatients_deep_totals <- round(periobaseline_physeq_norm_filter_onlypatients_deep_totals, 2)

#Rename Columns
colnames(periobaseline_physeq_norm_filter_onlypatients_deep_totals) <- "Percentage"

#Export
write.table(periobaseline_physeq_norm_filter_onlypatients_deep_totals, file="periobaseline_physeq_norm_filter_onlypatients_deep_totals.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

#END
