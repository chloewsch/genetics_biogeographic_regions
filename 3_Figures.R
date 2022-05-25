library(brms)
library(tidyverse)
library(tidybayes)
library(extrafont)
library(viridis)
library(patchwork)

# Coefficient plot (Fig. 2) ####
gdmod.df <- data.frame(Variable = rownames(summary(m_gdSARpr)$fixed)[2],
                       Coefficient = summary(m_gdSARpr)$fixed[2, 1],
                       CIlo95 = summary(m_gdSARpr)$fixed[2,3],
                       CIup95 = summary(m_gdSARpr)$fixed[2,4],
                       CIlo90 = posterior_interval(m_gdSARpr, pars = "b_scale_logdistk", prob = 0.90)[,1],
                       CIup90 = posterior_interval(m_gdSARpr, pars = "b_scale_logdistk", prob = 0.90)[,2],
                       Response_var = "genetic diversity")

tdmod.df <- data.frame(Variable = rownames(summary(m_tdpr)$fixed)[2],
                       Coefficient = summary(m_tdpr)$fixed[2, 1],
                       CIlo95 = summary(m_tdpr)$fixed[2,3],
                       CIup95 = summary(m_tdpr)$fixed[2,4],
                       CIlo90 = posterior_interval(m_tdpr, pars = "b_scale_logdistk", prob = 0.90)[,1],
                       CIup90 = posterior_interval(m_tdpr, pars = "b_scale_logdistk", prob = 0.90)[,2],
                       Response_var = "population density")

fstmod.df <- data.frame(Variable = rownames(summary(m_fstSARpr)$fixed)[2],
                        Coefficient = summary(m_fstSARpr)$fixed[2, 1],
                        CIlo95 = summary(m_fstSARpr)$fixed[2,3],
                        CIup95 = summary(m_fstSARpr)$fixed[2,4],
                        CIlo90 = posterior_interval(m_fstSARpr, pars = "b_scale_logdistk", prob = 0.90)[,1],
                        CIup90 = posterior_interval(m_fstSARpr, pars = "b_scale_logdistk", prob = 0.90)[,2],
                        Response_var = "genetic differentiation")

nemod.df <- data.frame(Variable = rownames(summary(m_neSARpr)$fixed)[2],
                       Coefficient = summary(m_neSARpr)$fixed[2, 1],
                       CIlo95 = summary(m_neSARpr)$fixed[2,3],
                       CIup95 = summary(m_neSARpr)$fixed[2,4],
                       CIlo90 = posterior_interval(m_neSARpr, pars = "b_scale_logdistk", prob = 0.90)[,1],
                       CIup90 = posterior_interval(m_neSARpr, pars = "b_scale_logdistk", prob = 0.90)[,2],
                       Response_var = "effective population size")

allModelFrame <- data.frame(rbind(gdmod.df, tdmod.df, fstmod.df, nemod.df))

allModelFrame$Response_var <- factor(allModelFrame$Response_var, levels = c("population density",
                                                              "genetic differentiation",
                                                              "genetic diversity",
                                                              "effective population size"))

# Species-specific effect data
taxoinfo <- read.csv("taxoinfo.csv", header = TRUE)

spp_plotdat <- function(model){
  model %>% 
    spread_draws(b_scale_logdistk, r_species[species, term]) %>%
    mutate(species_mean = b_scale_logdistk + r_species) %>%
    filter(term == "scale_logdistk") %>% 
    full_join(., taxoinfo, by = "species") %>%
    drop_na() %>% 
    mutate(species = gsub("_", " ", species)) %>% 
    group_by(order) %>% 
    arrange(species, .by_group = TRUE) %>% 
    mutate(species = factor(species, levels = unique(.$species)))
}

## gene diversity
gdplotdat <- spp_plotdat(m_gdSARpr) %>% 
  group_by(uid, family, order, species) %>% 
  summarise(sp_effect = mean(species_mean)) %>% 
  ungroup() %>% 
  mutate(Response_var = "genetic diversity")

gd_dist$species <- gsub("_", " ", gd_dist$species)
sitesg <- gd_dist %>% 
  group_by(species) %>% 
  summarise(sites = n()) %>% 
  ungroup()
gdplotdat <- merge(gdplotdat, sitesg, by = "species")

## population density
tdplotdat <- spp_plotdat(m_tdpr) %>% 
  group_by(uid, family, order, species) %>% 
  summarise(sp_effect = mean(species_mean)) %>% 
  ungroup() %>% 
  mutate(Response_var = "population density")

td_dist$species <- gsub("_", " ", td_dist$species)
sitest <- td_dist %>% 
  group_by(species) %>% 
  summarise(sites = n()) %>% 
  ungroup()
tdplotdat <- merge(tdplotdat, sitest, by = "species")

## FST
fstplotdat <- spp_plotdat(m_fstSARpr) %>% 
  group_by(uid, family, order, species) %>% 
  summarise(sp_effect = mean(species_mean)) %>% 
  ungroup() %>% 
  mutate(Response_var = "genetic differentiation")

fst_dist$species <- gsub("_", " ", fst_dist$species)
sitesf <- fst_dist %>% 
  group_by(species) %>% 
  summarise(sites = n()) %>% 
  ungroup()
fstplotdat <- merge(fstplotdat, sitesf, by = "species")

## effective population size
neplotdat <- spp_plotdat(m_neSARpr) %>% 
  group_by(uid, family, order, species) %>% 
  summarise(sp_effect = mean(species_mean)) %>% 
  ungroup() %>% 
  mutate(Response_var = "effective population size")

ne_dist$species <- gsub("_", " ", ne_dist$species)
sitesn <- ne_dist %>% 
  group_by(species) %>% 
  summarise(sites = n()) %>% 
  ungroup()
neplotdat <- merge(neplotdat, sitesn, by = "species")

## Plot ####
speciesfx <- rbind(gdplotdat, tdplotdat, fstplotdat, neplotdat)
speciesfx$Response_var <- factor(speciesfx$Response_var, levels = c("population density",
                                                                    "genetic differentiation",
                                                                    "genetic diversity",
                                                                    "effective population size"))

ggplot(allModelFrame) + 
  geom_hline(yintercept=seq(-1.5, 0.5, 0.25),  # these are the x axis lines (horizontal but will be flipped later)
             lwd=0.25, colour="#DEDEDE") +
  geom_hline(yintercept = 0, colour = "black", lty = 2) +
  geom_jitter(data = speciesfx, aes(x = Response_var, y = sp_effect, size = sites), color = "#CBAECB", width = 0.1, alpha = 0.6) +
  #geom_text(data = speciesfx, aes(x = Response_var, y = sp_effect, label = species), size = 2, angle = 45) +
  geom_linerange(aes(x = Response_var, ymin = CIlo90,
                     ymax = CIup90),
                 lwd = 2.5, position = position_dodge(width = 1), color = "#8E0190") +
  geom_pointrange(aes(x = Response_var, y = Coefficient, ymin = CIlo95,
                      ymax = CIup95),
                  lwd = 1, position = position_dodge(width = 1),
                  shape = 21, fill = "white", stroke = 3, color = "#8E0190") +
  geom_vline(xintercept=seq(1.5, length(unique(allModelFrame$Response_var))-0.5, 1), # these are the y axis lines that appear btwn groups
             lwd=0.25, colour="#DEDEDE") +
  scale_x_discrete(limits = levels(allModelFrame$Response_var)) +
  coord_flip() + 
  theme_classic(base_size = 14) +
  theme(axis.ticks.y = element_blank()) +
  labs(x= "", y = "model coefficients", title = "") +
  theme(text=element_text(family="Roboto Medium"))
dev.off()

# species coefficients ####
spp_ploth <- function(data, title){
  ggplot(data=data, aes(y = species, x = species_mean, color=order)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    stat_pointinterval(.width = c(0.9, 0.95)) +
    theme_minimal(base_size = 9) +
    labs(title = title, 
         x = "\u03b2 species",
         y = "") +
    scale_color_viridis(option="inferno", discrete=TRUE, begin = 0.3, end=0.8) +
    #scale_y_discrete(limits = rev(levels(data$species))) +
    coord_flip() +
    theme(text=element_text(family="Roboto"),
          axis.text.x = element_text(face = "italic")) +
    theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0, hjust=1))
}

gdph <- spp_ploth(spp_plotdat(m_gdSARpr), "genetic diversity")
tdph <- spp_ploth(spp_plotdat(m_tdpr), "population density")
fstph <- spp_ploth(spp_plotdat(m_fstSARpr), "population differentiation")
neph <- spp_ploth(spp_plotdat(m_neSARpr), "effective population size")

# interaction plots ######
p <- conditional_effects(m_tdpredge1, effects = 'scale_logdistk:edge_type')

p1 <- plot(p, plot = FALSE)[[1]] +
  xlab('log distance to region boundary') +
  ylab('log population density') +
  scale_colour_manual(values = c('#C689B2', '#000000'),
                      name = 'edge type') +
  scale_fill_manual(values = c('#C689B2', '#000000'),
                    name = 'edge type') +
  theme_minimal() +
  theme(text=element_text(family="Roboto"))

# GD
pgd <- conditional_effects(m_gdSARpredge1, effects = 'scale_logdistk:edge_type')

p2 <- plot(pgd, plot = FALSE)[[1]] +
  xlab('log distance to region boundary') +
  ylab('genetic diversity') +
  scale_colour_manual(values = c('#C689B2', '#000000'),
                      name = 'edge type') +
  scale_fill_manual(values = c('#C689B2', '#000000'),
                    name = 'edge type') +
  theme_minimal() +
  theme(text=element_text(family="Roboto"))

# Ne
pne <- conditional_effects(m_neSARpredge1, effects = 'scale_logdistk:edge_type')

p3 <- plot(pne, plot = FALSE)[[1]] +
  xlab('log distance to region boundary') +
  ylab('effective population size') +
  scale_colour_manual(values = c('#C689B2', '#000000'),
                      name = 'edge type') +
  scale_fill_manual(values = c('#C689B2', '#000000'),
                    name = 'edge type') +
  theme_minimal() +
  theme(text=element_text(family="Roboto"))

# FST
pfst <- conditional_effects(m_fstSARpredge1, effects = 'scale_logdistk:edge_type')

p4 <- plot(pfst, plot = FALSE)[[1]] +
  xlab('log distance to region boundary') +
  ylab('population differentiation') +
  scale_colour_manual(values = c('#C689B2', '#000000'),
                      name = 'edge type') +
  scale_fill_manual(values = c('#C689B2', '#000000'),
                    name = 'edge type') +
  theme_minimal() +
  theme(text=element_text(family="Roboto"))

# Combine plots:
(p2+p3)/(p4+p1) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')


# Scatterplot raw data ####
# GD
gdscat <- ggplot(data = gd_dist, aes(y = He, x = log(dist_nearest_BGR_km), 
                           color = species)) +
  geom_point() +
  scale_color_viridis(option = 'inferno', 
                      discrete = T,
                      guide = "none") +
  labs(y = "genetic diversity", x = "log distance to nearest edge")+
  theme_minimal() +
  theme(text=element_text(family="Roboto"))

# Ne
nescat <- ggplot(data = ne_dist, aes(y = log(Ne), x = log(dist_nearest_BGR_km),
                           color = species)) +
  geom_point() +
  scale_color_viridis(option = 'inferno', 
                      discrete = T,
                      guide = "none") +
  labs(y = "log effective population size", x = "log distance to nearest edge")+
  theme_minimal() +
  theme(text=element_text(family="Roboto"))


# FST
fstscat <- ggplot(data = fst_dist, aes(y = global_fst, x = log(dist_nearest_BGR_km),
                           color = species)) +
  geom_point() +
  scale_color_viridis(option = 'inferno', 
                      discrete = T,
                      guide = "none") +
  labs(y = 'population differentiation', x = "log distance to nearest edge")+
  theme_minimal() +
  theme(text=element_text(family="Roboto"))

# TD
tdscat <- ggplot(data = td_dist, aes(y = log(mean_density_yr), x = log(dist_nearest_BGR_km),
                           color = species)) +
  geom_point() +
  scale_color_viridis(option = 'inferno', 
                      discrete = T,
                      guide = "none") +
  labs(y = 'log population density', x = "log distance to nearest edge")+
  theme_minimal() +
  theme(text=element_text(family="Roboto"))

# Combine plots:
(gdscat + nescat)/(fstscat + tdscat)

# Boxplot distance to edge for edge types ####
gd_dist_box <- ggplot(gd_dist, aes(y = log(dist_nearest_BGR_km), x = edge_type)) + 
  geom_boxplot() + 
  labs(y = "log distance to nearest edge", x = "", 
       title = "Genetic diversity")+
  theme_minimal() +
  theme(text=element_text(family="Roboto"))

ne_dist_box <- ggplot(ne_dist, aes(y = log(dist_nearest_BGR_km), x = edge_type)) + 
  geom_boxplot() + 
  labs(y = "log distance to nearest edge", x = "", 
       title = "Effective population size")+
  theme_minimal() +
  theme(text=element_text(family="Roboto"))

fst_dist_box <- ggplot(fst_dist, aes(y = log(dist_nearest_BGR_km), x = edge_type)) + 
  geom_boxplot() + 
  labs(y = "log distance to nearest edge", x = "", 
       title = "Population differentiation")+
  theme_minimal() +
  theme(text=element_text(family="Roboto"))

td_dist_box <- ggplot(td_dist, aes(y = log(dist_nearest_BGR_km), x = edge_type)) + 
  geom_boxplot() + 
  labs(y = "log distance to nearest edge", x = "", 
       title = "Population density")+
  theme_minimal() +
  theme(text=element_text(family="Roboto"))

# Combine plots:
(gd_dist_box + ne_dist_box)/(fst_dist_box + td_dist_box)

# histograms of number of biogeographic regions species are sampled in ####
mamm_regions <- read_sf("Wallace Mammal Ecoregions/MammalEcoregions_WGS84_clean.shp", crs = 4326)
ne_sf <- st_as_sf(ne_dist, coords = c("lon", "lat"), crs = 4326)
td_sf <- st_as_sf(td_dist, coords = c("lon", "lat"), crs = 4326)
fst_sf <- st_as_sf(fst_dist, coords = c("lon", "lat"), crs = 4326)
gd_sf <- st_as_sf(gd_dist, coords = c("lon", "lat"), crs = 4326)

gd_bound <- st_join(gd_sf, mamm_regions) 
ne_bound <- st_join(ne_sf, mamm_regions) 
fst_bound <- st_join(fst_sf, mamm_regions) 
td_bound <- st_join(td_sf, mamm_regions) 

gd_reg <- gd_bound %>% 
  group_by(species, big_rgn) %>% 
  summarise() %>% 
  drop_na(big_rgn) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  summarise(num_regions = n())

ne_reg <- ne_bound %>% 
  group_by(species, big_rgn) %>% 
  summarise() %>%
  drop_na(big_rgn) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  summarise(num_regions = n())

fst_reg <- fst_bound %>% 
  group_by(species, big_rgn) %>% 
  summarise() %>%
  drop_na(big_rgn) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  summarise(num_regions = n())

td_reg <- td_bound %>% 
  group_by(species, big_rgn) %>% 
  summarise() %>% 
  drop_na(big_rgn) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  summarise(num_regions = n())

a <- ggplot(data = gd_reg, aes(x = num_regions)) + 
  geom_histogram(bins = 5) + 
  labs(title = "Genetic diversity", x = "number of regions") +
  theme_classic() +
  theme(text=element_text(family="Roboto"))

b <- ggplot(data = ne_reg, aes(x = num_regions)) + 
  geom_histogram(bins = 4) + 
  labs(title = "Effective population size", x = "number of regions") +
  scale_x_continuous(breaks=c(1,2,3)) +
  theme_classic() +
  theme(text=element_text(family="Roboto"))

c <- ggplot(data = fst_reg, aes(x = num_regions)) + 
  geom_histogram(bins = 4) + 
  labs(title = "Population differentiation", x = "number of regions") +
  scale_x_continuous(breaks=c(1,2,3)) +
  theme_classic() +
  theme(text=element_text(family="Roboto"))

d <- ggplot(data = td_reg, aes(x = num_regions)) + 
  geom_histogram(bins = 3) + 
  labs(title = "Population density", x = "number of regions") +
  theme_classic() +
  theme(text=element_text(family="Roboto"))

# Combine plots:
(a+b)/(c+d)