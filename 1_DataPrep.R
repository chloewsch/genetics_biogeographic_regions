#### NOTE go through and remove all the intermediate write.csvs
# replace file names with final archived data

## Code accompanying: Population demography maintains biogeographic boundaries (2022)

# 1 - Data preparation #

library(tidyverse)
library(sf) # v1.0.1
library(sp) #v1.4.5

# Distance calculation
library(geosphere) # v1.5.10

# Taxonomy
library(taxize) # v0.9.99


# Cleaning data ####

## 1) Bioregion data ####
# Distribution-based regions from Holt et al: https://www.science.org/doi/10.1126/science.1228282
mamm_regions <- read_sf("Wallace Mammal Ecoregions/MammalEcoregions_WGS84.shp", crs = 4326)

# Crop to the Americas
mb <- st_buffer(mamm_regions, dist = 0)
reg <- st_crop(mb, xmin = -180, ymin = -56, xmax = -26.5, ymax = 83)

# Reproject into equal area projection (South America azimuthal equal area) 
az <- "+proj=laea +lat_0=-10 +lon_0=-70 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"
reg_AEA <- st_transform(reg, az)

# create an ID column in biogeo regions:
reg_AEA$ID <- as.numeric(c(1:nrow(reg_AEA)))
reg_AEA$biogeo_region <- as.factor(paste("BGR", reg_AEA$ID, sep=""))

### Remove islands:
## Galapagos and Hawaii (part of region 2, which are attached to a bigger region in SA)
BGR2 <- reg_AEA %>% 
  filter(biogeo_region=="BGR2") %>% 
  st_cast("POLYGON") %>% 
  mutate(sub_id = as.factor(paste("piece", 1:nrow(.), sep = "")))

BGR2_2_13 <- BGR2 %>% 
  filter(sub_id %in% c("piece2", "piece13")) %>% 
  st_combine() %>% 
  st_as_sf()

# Replace BGR2 in region data with new one:
st_geometry(reg_AEA[2,]) <- st_geometry(BGR2_2_13)

## Falkland Islands (region 1; 154 polygons; crop by longitude)
BGR1 <- reg_AEA %>% 
  filter(biogeo_region=="BGR1") %>% 
  st_cast("POLYGON")

BGR1k <- st_crop(BGR1, xmin = -458858.0, ymin = -4972366.6, xmax = 630000, ymax = -2649209.5) %>% 
  st_combine() %>% 
  st_as_sf()

# Replace BGR1 in region data with new one:
st_geometry(reg_AEA[1,]) <- st_geometry(BGR1k)

## Carribbean islands that are part of BGR4:
BGR4 <- reg_AEA %>% 
  filter(biogeo_region=="BGR4") %>% 
  st_cast("POLYGON") %>% 
  mutate(sub_id = as.factor(paste("piece", 1:nrow(.), sep = "")))

BGR4_1 <- BGR4 %>% 
  filter(sub_id %in% c("piece1")) %>% 
  st_combine() %>% 
  st_as_sf()

# Replace BGR4 in region data with new one:
st_geometry(reg_AEA[4,]) <- st_geometry(BGR4_1)


## Remove region 6 (more Carribbean), edge of Alaska, some island in the North Atlantic
reg_AEA <- reg_AEA %>%
  filter(!biogeo_region %in% c("BGR6", "BGR9", "BGR10", "BGR11"))
reg_AEA$biogeo_region <- droplevels(reg_AEA$biogeo_region)

## Plot & check regions:
ggplot() +
  geom_sf(data=reg_AEA, aes(fill = biogeo_region))

## Remove un-needed data:
rm(BGR1, BGR1k, BGR2, BGR2_2_13, BGR4, BGR4_1, mb, reg, az)

## convert to spatial points:
reg_WGS <- st_transform(reg_AEA, crs = 4326)
reg_WGS_sp <- as(reg_WGS, 'Spatial')

## Write to shapefile:
#rgdal::writeOGR(reg_WGS_sp, dsn = 'Wallace Mammal Ecoregions', layer = 'MammalEcoregions_WGS84_clean', driver = "ESRI Shapefile")

## 2) Genetic diversity (expected heterozygosity) #####
# Macropopgen database (Lawrence et al): https://figshare.com/articles/dataset/MacroPopGen_Database_Geo-referenced_population-specific_microsatellite_data_across_the_American_continents/7207514
macropopgen <- read.csv("Macropopgen.csv")

# Mammal subset:
mammal_gd <- macropopgen %>% 
  filter(TaxaClass == "Mammalia", Long < 0, Continent != "Carr") %>% # 1 point not in western hemisphere; remove islands
  drop_na(Long, He)

# Create a site ID column
mammal_gd$siteID <- as.factor(paste("GD_", c(1:nrow(mammal_gd)), sep = ""))

gd_sf <- st_as_sf(mammal_gd, coords = c("Long", "Lat"), crs = 4326) # convert to sf

rm(macropopgen)

## 3) Population density data ####
# TetraDensity v1 database (Santini et al): https://figshare.com/articles/dataset/TetraDENSITY_Population_Density_dataset/5371633

tetra <- read.csv("TetraDENSITY_v.1.csv", header = TRUE) # all density measured in ind/km2

# Mammal subset
# Note all density measured in ind/km2
tetra <- tetra %>% 
  filter(Class == "Mammalia") %>% 
  drop_na(Longitude, Density)
tetra$siteID <- as.factor(paste("TD_", c(1:nrow(tetra)), sep = ""))

tetra_sf0 <- st_as_sf(tetra, coords = c("Longitude", "Latitude"), crs = 4326) # convert to sf

# Crop to Americas:
tetra_am <- st_crop(tetra_sf0, xmin = -180, ymin = -56, xmax = -26.5, ymax = 83)

# Find and remove islands!
galapagos <- grep("Gal", tetra_am$Country)
car <- which(tetra_am$siteID %in% c("TD_340", "TD_341", "TD_342", "TD_343", "TD_344", "TD_345", "TD_346")) 

nomoreislands <- c(galapagos, car) ## 11 sites

tetra_am <- tetra_am[-c(nomoreislands),]

# remove hawaii
hawaii_pts <- data.frame(lon = -156, lat = 20)
hawaii_pts <- st_as_sf(hawaii_pts, coords = c("lon", "lat"), crs = 4326) 
hawaii_box <- st_as_sfc(st_bbox(st_buffer(hawaii_pts, dist = 500*1000)))
ggplot() + geom_sf(data = tetra_am) + geom_sf(data = hawaii_box, fill = NA)

tetra_noisl <- st_difference(tetra_am, hawaii_box)

## Average within years
# Convert to dataframe:
tetra_amdf <- as.data.frame(tetra_noisl)
tetra_amdf$lon <- st_coordinates(tetra_noisl)[,1]
tetra_amdf$lat <- st_coordinates(tetra_noisl)[,2]

td_dist_yearavg <- tetra_amdf %>% 
  group_by(Genus, Species, lat, lon) %>% 
  mutate(mean_density_yr = mean(Density)) %>% 
  distinct(Genus, Species, lat, lon, .keep_all = TRUE)

tetra_sf0 <- st_as_sf(td_dist_yearavg, coords = c("lon", "lat"), crs = 4326) # convert to sf

## 4) Population differentiation (FST) and effective population size (Ne) data ####
# from Schmidt et al. https://datadryad.org/stash/dataset/doi:10.5061/dryad.cz8w9gj0c

sgdata <- read.csv("sgdata.csv", head = T)

## Population differentiation
sgdataf <- sgdata %>% 
  drop_na(global_fst) # remove rows with no FST

# Convert to sf
fst_sf <- st_as_sf(sgdataf, coords = c("lon", "lat"), crs = 4326)

## Effective population size
sgdatan <- sgdata %>% 
  drop_na(Ne) # remove rows with no Ne

# Convert to sf
ne_sf <- st_as_sf(sgdatan, coords = c("lon", "lat"), crs = 4326)


# Distance to nearest edge ####
## Genetic diversity ####
gd_sf_sp <- as(gd_sf, 'Spatial') # convert to sp
dist_gd <- dist2Line(gd_sf_sp, reg_WGS_sp, distfun = distGeo) # geodesic distance

mammal_gd$dist_nearest_BGR <- (dist_gd[,1])
mammal_gd$dist_nearest_BGR_km <- (dist_gd[,1])/1000 #convert to km
mammal_gd <- mammal_gd %>%
  select(-Location) # Remove location column for writing csv (has commas)

gd_dist <- mammal_gd

## Population density ####
tetra_sf_sp <- as(tetra_sf0, 'Spatial')
dist_td <- dist2Line(tetra_sf_sp, reg_WGS_sp, distfun = distGeo)

tetra_am$dist_nearest_BGR <- (dist_td[,1])
tetra_am$dist_nearest_BGR_km <- (dist_td[,1])/1000

## convert to dataframe for export:
tetdat <- as.data.frame(tetra_am)

# reattach coordinates:
tetdatxy <- st_coordinates(tetra_am)
tetdat$lon <- tetdatxy[,1]
tetdat$lat <- tetdatxy[,2]

tetdat <- tetdat %>% 
  select(-c(geometry, Locality)) # Remove Locality column for writing csv (has commas)

td_dist <- tetdat

## FST ####
fst_sf_sp <- as(fst_sf, 'Spatial')

# Crop regions to avoid errors:
reg2 <- st_crop(reg_WGS, xmin = -168, ymin = 15, xmax = -26.5, ymax = 83)
reg2_sp <- as(reg2, 'Spatial')
# or alternatively (to avoid geometrycollection):
reg2 <- mamm_regions %>% 
  filter(ID %in% c(7, 8, 12))
reg2_sp <- as(reg2, 'Spatial')

dist_fst <- dist2Line(fst_sf_sp, reg2_sp, distfun = distGeo)

sgdataf$dist_nearest_BGR <- (dist_fst[,1])
sgdataf$dist_nearest_BGR_km <- (dist_fst[,1])/1000
fst_dist <- sgdataf

## Ne ####
ne_sf_sp <- as(ne_sf, 'Spatial')

dist_ne <- dist2Line(ne_sf_sp, reg2_sp, distfun = distGeo) # use North America regions

sgdatan$dist_nearest_BGR <- (dist_ne[,1])
sgdatan$dist_nearest_BGR_km <- (dist_ne[,1])/1000
ne_dist <- sgdatan

# Edge classification ####
# Identify coasts:
coastlines <- mamm_regions %>%
  st_union() # dissolve all inner boundaries

# Convert to sp
coast_sp <- as(coastlines, 'Spatial')

## Genetic diversity ####
gd_coast <- dist2Line(gd_sf_sp, coast_sp, distfun = distGeo)

gd_dist$edge_type <- ifelse(round(gd_dist$dist_nearest_BGR, 2) == round(gd_coast[,1], 2), "coast", "interior")

## Population density ####
td_coast <- dist2Line(tetra_sf_sp, coast_sp, distfun = distGeo) # geodesic distance
td_dist$edge_type <- ifelse(round(td_dist$dist_nearest_BGR, 2) == round(td_coast[,1], 2), "coast", "interior")

## FST ####
coastN <-  reg2 %>%
  st_union()
coastN_sp <- as(coastN, 'Spatial')

fst_coast <- dist2Line(fst_sf_sp, coastN_sp, distfun = distGeo)
fst_dist$edge_type <- ifelse(round(fst_dist$dist_nearest_BGR, 2) == round(fst_coast[,1], 2), "coast", "interior")

## Ne ####
ne_coast <- dist2Line(ne_sf_sp, coastN_sp, distfun = distGeo)
ne_dist$edge_type <- ifelse(round(ne_dist$dist_nearest_BGR, 2) == round(ne_coast[,1], 2), "coast", "interior")

# Taxonomy ####

## Taxonomic data for all species ##
species <- unique(c(gd_dist$G_s, td_dist$G_s, fst_dist$species, ne_dist$species))
species_list <- as.list(species)

# Get NCBI unique identifier (UID) for each species:
uids <- lapply(species_list, function(x) get_uid(x, messages = FALSE)[1])

uid_strip <- unlist(unique(uids))

tax_df <- data.frame(uid = uid_strip,
                     species = species)
missingtax <- (tax_df[is.na(tax_df$uid),])

missingtax$family <- c("Callitrichidae","Sciuridae", "Sciuridae", "Soricidae", "Callitrichidae", "Callitrichidae", 
                       "Pitheciidae", "Pitheciidae", "Pitheciidae", "Pitheciidae", "Pitheciidae", "Mustelidae")
missingtax$order <- c("Primates", "Rodentia", "Rodentia", "Eulipotyphla", "Primates", "Primates", 
                      "Primates", "Primates", "Primates", "Primates", "Primates", "Carnivora") 

missingtax <- missingtax %>% 
  select(c(uid, family, order, species))

# Assign classifications:
taxize_class <- classification(uids, db = "ncbi")

pulltax <- lapply(taxize_class, function(x) as.data.frame(x))

pulltax1 <- tibble(uid = names(pulltax), pulltax) %>% 
  unnest(cols = c(pulltax)) %>% 
  filter(rank %in% c("order","family", "species")) %>% 
  select(-c(id, x)) %>% 
  distinct() %>%
  spread(rank, name) %>% 
  mutate(species = gsub(" ", "_", species))

taxoinfo <- rbind(pulltax1, missingtax)

## Genetic diversity ####
gd_dist <- gd_dist %>% select(-c(Family, Genus, Species)) %>% 
  rename(class = TaxaClass)

gd_dist_order <- merge(taxoinfo %>% select(species, family, order), gd_dist, by.x = "species", by.y = "G_s",
                       all.x = FALSE, all.y = TRUE)
# Update species names
gd_dist_order$species <- gsub("Xerospermosphilus_polionotus", "Xerospermophilus_polionotus", gd_dist_order$species)
gd_dist_order$species <- gsub("Xerospermosphilus_perotensis", "Xerospermophilus_perotensis", gd_dist_order$species)
gd_dist_order$species <- gsub("Saguinus_nigricollis", "Leontocebus_nigricollis", gd_dist_order$species)
gd_dist_order$species <- gsub("Callicebus_lucifer", "Cheracebus_lucifer", gd_dist_order$species)
gd_dist_order$species <- gsub("Callicebus_pallescens", "Plecturocebus_pallescens", gd_dist_order$species)
gd_dist_order$species <- gsub("Callicebus_regulus", "Cheracebus_regulus", gd_dist_order$species)


## Population density ####
# TD has taxonomomic info, but change G_s column to species to match other data
td_dist <- td_dist %>% 
  unite(species, c("Genus", "Species")) %>% 
  rename(class = Class, order = Order, family = Family)

## FST ####
fst_dist <- merge(taxoinfo %>% select(species, family, order), fst_dist, by = "species",
                  all.x = FALSE, all.y = TRUE)

## Ne ####
ne_dist <- merge(taxoinfo %>% select(species, family, order), ne_dist, by = "species",
                 all.x = FALSE, all.y = TRUE)

# Write data ####
#write.csv(gd_dist, "gd_dist.csv", row.names = F, quote = F)
#write.csv(td_dist, "td_dist.csv", row.names = F, quote = F)
#write.csv(fst_dist, "fst_dist.csv", row.names = F, quote = F)
#write.csv(ne_dist, "ne_dist.csv", row.names = F, quote = F)