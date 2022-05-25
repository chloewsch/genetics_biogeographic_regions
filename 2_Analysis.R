## Distance from sites to nearest bioregion edge ##

## 2. Analysis
library(tidyverse)
library(brms)

# Spatial regs
library(sf)
library(adespatial)
library(vegan)
library(spdep)

# Plots
library(tidybayes)
library(extrafont)
library(viridis)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)

# Load & prep data ####
gd_dist <- read.csv("gd_dist.csv", header = TRUE)
td_dist <- read.csv("td_dist.csv", header = TRUE)
fst_dist <- read.csv("fst_dist.csv", header = TRUE)
ne_dist <- read.csv("ne_dist.csv", header = TRUE)

taxoinfo <- read.csv("taxoinfo.csv", header = TRUE)

## scale & log variables
gd_dist$scale_He     <- scale(gd_dist$He)
td_dist$scale_logDen <- scale(log(td_dist$mean_density_yr))
fst_dist$scale_fst   <- scale(fst_dist$global_fst)
ne_dist$scale_logNe   <- scale(log(ne_dist$Ne))
gd_dist$scale_logdistk   <- scale(log(gd_dist$dist_nearest_BGR_km + 1))
td_dist$scale_logdistk   <- scale(log(td_dist$dist_nearest_BGR_km + 1))
fst_dist$scale_logdistk  <- scale(log(fst_dist$dist_nearest_BGR_km + 1))
ne_dist$scale_logdistk  <- scale(log(ne_dist$dist_nearest_BGR_km + 1))

# Preliminary models (not spatial regressions) ####
## genetic diversity #####
m_gdpr <- brm(bf(scale_He ~ scale_logdistk + (scale_logdistk|species)), 
            cores = 4, chains = 4, iter = 4000,
            prior = prior(normal(0,1), class = b, coef = scale_logdistk),
            control = list(adapt_delta = 0.95, max_treedepth = 10),
            data = gd_dist)

# Check spatial autocorrelation:
gd_sites <- st_as_sf(gd_dist, coords = c("Long", "Lat"), crs = 4326)
GD_dupid <- which(duplicated(st_geometry(gd_sites)))
GD_dupes <- gd_sites[GD_dupid,]
GD_dupesj <- st_jitter(GD_dupes, 0.01)

GD_jit <- gd_sites
st_geometry(GD_jit)[GD_dupid] <- st_geometry(GD_dupesj)

xygd <- st_coordinates(GD_jit)
nb1 <- chooseCN(xygd, ask = FALSE, type = 2, result.type = "listw") # minimum distance (km) to not get an empty neighbor set

gdres <- residuals(m_gdpr)[,1]
moran.test(gdres, nb1) # yes, I = 0.13

## population density #####
m_tdpr <- brm(bf(scale_logDen ~ scale_logdistk + (scale_logdistk|species)), 
              cores = 4, chains = 4, iter = 2000,
              prior = prior(normal(0,1), class = b, coef = scale_logdistk),
              control = list(adapt_delta = 0.95, max_treedepth = 10),
              data = td_dist)

# Check spatial autocorrelation:
td_sites <- st_as_sf(td_dist, coords = c("lon", "lat"), crs = 4326)
TD_dupid <- which(duplicated(st_geometry(td_sites)))
TD_dupes <- td_sites[TD_dupid,]
TD_dupesj <- st_jitter(TD_dupes, 0.1)

TD_jit <- td_sites
st_geometry(TD_jit)[TD_dupid] <- st_geometry(TD_dupesj)

xytd <- st_coordinates(TD_jit)
nb2 <- chooseCN(xytd, type = 5, d1 = 9, d2 = 15, result.type = "listw")

tdres <- residuals(m_tdpr)[,1]
moran.test(tdres, nb2) # none


## FST #####
m_fstpr <- brm(bf(scale_fst ~ scale_logdistk + (scale_logdistk|species)), 
             cores = 4, chains = 4, iter = 2000, 
             control = list(adapt_delta = 0.95, max_treedepth = 10),
             data = fst_dist)

## Check spatial autocorrelation:
fst_sites <- st_as_sf(fst_dist, coords = c("lon", "lat"), crs = 4326)
fst_dupid <- which(duplicated(st_geometry(fst_sites)))
fst_dupes <- fst_sites[fst_dupid,]
fst_dupesj <- st_jitter(fst_dupes, 0.01)

fst_jit <- fst_sites
st_geometry(fst_jit)[fst_dupid] <- st_geometry(fst_dupesj)

xyfst <- st_coordinates(fst_jit)
nb3 <- chooseCN(xyfst, ask = FALSE, type = 2, result.type = "listw") # minimum distance (km) to not get an empty neighbor set

fstres <- residuals(m_fstpr)[,1]
moran.test(fstres, nb3) # yes, I = 0.08

## effective population size #####
m_nepr <- brm(bf(scale_logNe ~ scale_logdistk + (scale_logdistk|species)), 
            cores = 4, chains = 4, iter = 2000,
            prior = prior(normal(0,1), class = b, coef = scale_logdistk),
            control = list(adapt_delta = 0.95, max_treedepth = 10),
            data = ne_dist)

# Check spatial autocorrelation:
ne_sites <- st_as_sf(ne_dist, coords = c("lon", "lat"), crs = 4326)
ne_dupid <- which(duplicated(st_geometry(ne_sites)))
ne_dupes <- ne_sites[ne_dupid,]
ne_dupesj <- st_jitter(ne_dupes, 0.01)
ne_jit <- ne_sites
st_geometry(ne_jit)[ne_dupid] <- st_geometry(ne_dupesj)

xyne <- st_coordinates(ne_jit)
nb4 <- chooseCN(xyne, ask = FALSE, type = 2, result.type = "listw") # minimum distance (km) to not get an empty neighbor set

neres <- residuals(m_nepr)[,1]
moran.test(neres, nb4) # yes, I = 0.07


# SAR models ####

## genetic diversity #####
# nb network:
m_gdSARpr <- brm(scale_He ~ scale_logdistk + (scale_logdistk|species) + sar(nb1, type = "lag"),
                 cores = 4, chains = 4, iter = 3000,
                 prior = prior(normal(0,1), class = b, coef = scale_logdistk),
                 control = list(adapt_delta = 0.95, max_treedepth = 10),
                 data = gd_dist, data2 = list(nb1 = nb1))

## FST #####
m_fstSARpr <- brm(scale_fst ~ scale_logdistk + (scale_logdistk|species) + sar(nb3, type = "lag"),
                  cores = 4, chains = 4, iter = 2000,
                  prior = prior(normal(0,1), class = b, coef = scale_logdistk),
                  control = list(adapt_delta = 0.95, max_treedepth = 10),
                  data = fst_dist, data2 = list(nb3 = nb3))

## effective population size  #####
m_neSARpr <- brm(scale_logNe ~ scale_logdistk + (scale_logdistk|species) + sar(nb4, type = "lag"),
                 cores = 4, chains = 4, iter = 2000,
                 prior = prior(normal(0,1), class = b, coef = scale_logdistk),
                 control = list(adapt_delta = 0.95, max_treedepth = 10),
                 data = ne_dist, data2 = list(nb4 = nb4))

# Edge interaction ####
## gene diversity ####
m_gdSARpredge1 <- brm(scale_He ~ scale_logdistk + edge_type + scale_logdistk*edge_type + 
                        (scale_logdistk + edge_type + scale_logdistk*edge_type|species) + 
                        sar(nb1, type = "lag"),
                      cores = 4, chains = 4, iter = 3000,
                      prior = prior(normal(0,1), class = b),
                      control = list(adapt_delta = 0.95, max_treedepth = 10),
                      data = gd_dist, data2 = list(nb1 = nb1))

## population density ####
m_tdpredge1 <- brm(bf(scale_logDen ~ scale_logdistk + edge_type + scale_logdistk*edge_type + 
                        (scale_logdistk + edge_type + scale_logdistk*edge_type|species)), 
                   cores = 4, chains = 4, iter = 4000,
                   prior = prior(normal(0,1), class = b),
                   control = list(adapt_delta = 0.95, max_treedepth = 10),
                   data = td_dist)


## FST ####
m_fstSARpredge1 <- brm(scale_fst ~ scale_logdistk + edge_type + scale_logdistk*edge_type + 
                         (scale_logdistk + edge_type + scale_logdistk*edge_type|species) + 
                         sar(nb3, type = "lag"),
                       cores = 4, chains = 4, iter = 2000,
                       prior = prior(normal(0,1), class = b),
                       control = list(adapt_delta = 0.95, max_treedepth = 10),
                       data = fst_dist, data2 = list(nb3 = nb3))

## effective population size ####
m_neSARpredge1 <- brm(scale_logNe ~ scale_logdistk + edge_type + scale_logdistk*edge_type + 
                        (scale_logdistk + edge_type + scale_logdistk*edge_type|species) + 
                        sar(nb4, type = "lag"),
                      cores = 4, chains = 4, iter = 2000,
                      prior = prior(normal(0,1), class = b),
                      control = list(adapt_delta = 0.95, max_treedepth = 10),
                      data = ne_dist, data2 = list(nb4 = nb4))