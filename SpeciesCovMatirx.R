#####################################################
### Code for phylogeny matrix in random effects   ###
### Lynn et al. 2022 Ecography                    ###
### DOI: 10.1111/ecog.06114                       ###
#####################################################

# load packages
library(rotl); library(ape); library(devtools); library(stringi); 
library(tidyverse)

# data
herbdat <- read.csv("LynnETAL_2022_Ecography_sodCaseStudy.csv")

# pull species names
sppname <- herbdat %>% dplyr::select (species_name) %>% unique() %>%
  as.matrix() 

head(sppname)

# need to add underscores between genus_species
# check scientific names
names <- base::as.vector(sppname)
resolved_name <- rotl::tnrs_match_names(names = names)
resolved_name <- resolved_name %>% 
  mutate(unique_name = str_replace(unique_name, " ", "_"),
         unique_name = dplyr::recode(unique_name, 
                                     "Andropogon_gerardi" = "Andropogon_gerardii")) 


# load the big tree from Smith and Brown 2018
# can be found at https://github.com/FePhyFoFum/big_seed_plant_trees
bigtree <- read.tree("v0.1/ALLMB.tre")

# prune to needed tips
prunedtree <- keep.tip(bigtree, tip=resolved_name$unique_name)
plot(prunedtree) # visualize

#calculate the variance covariance matrix
vcv_tree <- vcv.phylo(prunedtree, corr=TRUE)

# fix speciees names to get them to match
vcvdat <- as.data.frame(vcv_tree) # 
vcvdat <- vcvdat %>% rownames_to_column(var="species_name") %>%
  mutate(species_name = str_replace(species_name, "_", " "),
         species_name = dplyr::recode(species_name, 
                               "Erythranthe guttata" = "Mimulus guttatus",
                               "Solidago canadensis" = "Solidago altissima")) %>%
  rename("Mimulus_guttatus" = "Erythranthe_guttata")

# ensure matrix is symmetric and convert to matrix
vcvdat <- vcvdat[order(vcvdat$species_name), order(names(vcvdat))]
vcvdat <-   vcvdat %>% rownames_to_column(var="X") %>% 
  dplyr::select(-X) %>% column_to_rownames(var="species_name") %>% as.matrix()
  