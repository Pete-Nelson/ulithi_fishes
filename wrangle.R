# Peter Nelson
# Institute of Marine Sciences
# University of California, Santa Cruz
# pnelson1@ucsc.edu

# Ulithi fishes
# read in and wrangle fish, site and environmental data
# created 26 May 2023
# last modified 6 June 2023

# maybe...
setwd("G:/My Drive/GitHub Gdrive Projects/ulithi_fishes")

# get started ----
library(tidyverse)
library(readxl)
library(janitor)

# SO many issues with variables in the original...changed in Excel!
# See now Tableau_Ulithi_Fish_2012_2023_Final_short_PN.xlsx

temp <- read_excel("data/fishes.xlsx",
                   sheet = "transects")
warnings()
glimpse(temp) # 'IP/TP' is a bad variable name and shouldn't be lgl anyway

# sites+ -----
# sites, location, habitat
temp2 <- read_excel("data/fishes.xlsx",
                    sheet = "Sites")
glimpse(temp2)
sites <- temp2 %>% 
  mutate(across(Site:Habitat, factor)) %>% 
  clean_names()
glimpse(sites)

# org parameters ----
# parameters for biomass estimation, trophic categories
temp3 <- read_excel("data/fishes.xlsx",
                    sheet = "Data")
glimpse(temp3)
trophic <- temp3 %>% 
  mutate(across(c(family, trophic_coarse:trophic_coarse_MP), factor)) %>% 
  select(-c('Trophic coarse MP', 'Trophic_finescale GB')) %>% 
  clean_names()


# counts -----
count1 <- rename(temp, phase = 'IP_TP') %>% 
  clean_names() %>% 
  mutate(across(c(site:depth, family, species, phase, morph, feeding), factor)) %>% 
  group_by(year, site, transect, depth, species) %>% # check with Giacomo!
  summarise(total = sum(counts))
glimpse(count1)

counts <- 
  count1 %>% 
  pivot_wider(names_from = species, values_from = total, values_fill = 0) %>% 
  clean_names(case = "mixed") %>% 
  mutate(count_total = sum(c_across(Acanthurus_lineatus:Plectorhynchus_albovittatus))) %>% 
  relocate(count_total, .after = depth)

# biomass -----
## to do -----
# Develop biomass df with higher level location data (eg Ulithi, Ifalik--island?) and with two trophic level designations only--confer w Giacomo and Michelle. 'biom2' is getting closer.
# Work on figures, probably in a separate code file.

# calculated here with all spp (including sharks, large schools)
biom1 <- rename(temp, phase = 'IP_TP') %>% 
  clean_names() %>% 
  mutate(across(c(site:depth, family, species, phase, morph, feeding), factor)) %>% 
  group_by(year, site, transect, depth, species) %>% # check with Giacomo!
  summarise(biom_total = sum(biomass))

## biom_site -----
biom_site <- biom1 %>% 
  group_by(year, site) %>% 
  summarise(total_site_biom = sum(biom_total)) %>% 
  left_join(., sites, by = "site")
biom2 <- biom1 %>% 
  left_join(., trophic, by = "species") %>% 
  relocate(family, .after = depth) %>% 
  select(c(year:biom_total, trophic_coarse:trophic_coarse_mp))

  

## biom_trans -----
biom_trans <- biom1 %>% 
  group_by(year, site, transect) %>% 
  summarise(trans_biom = sum(biom_total))


# clean-up ----
rm(temp, temp2, temp3)

# figures -----
## biomass -----
# raincloud plots
library(ggplot2)
library(ggdist)
library(tidyquant)

biom_site %>% 
  filter(island %in% c("Asor", "Falalop", "Federai", "Mogmog")) %>% 
  ggplot(aes(x = factor(island), y = total_site_biom, fill = factor(island))) +
  # add half-violin from {ggdist} pkg
  stat_halfeye(
    # customize smoothing ("bandwidth")
    adjust = 1,
    # move geom to right of (above) boxplot
    justification = -0.18,
    # remove slab interval (line through box)
    .width = 0,
    point_color = NA
  ) +
  geom_boxplot(
    # width of the boxes
    width = 0.24,
    # distance btwn boxes
    position = "dodge",
    # remove outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  stat_dots(
    shape = 21, size = 3, color = "darkgrey", fill = "darkgrey",
    # position to the left
    side = "right", # "left"
    # move geom to the left (below)
    justification = -0.2, # 1.15
    # binning for y values
    binwidth = 5
  ) +
  scale_fill_tq() +
  theme_tq() +
  labs(
    title = "Reef Fish Biomass",
    subtitle = "Main Ulithi Islands",
    x = "Islands",
    y = expression(Biomass ~ (g/m^2)),
    fill = "Islands"
  ) + 
  coord_flip()

base <- biom_site %>% 
  filter(island %in% c("Asor", "Falalop", "Federai", "Mogmog")) %>% 
  ggplot(
  aes(x = factor(island), 
      y = total_site_biom, 
      fill = factor(island))
)
base + geom_boxplot() +
  scale_fill_tq() +
  theme_tq()

# SCRATCH ###########################################
## example data -----
(temp <- counts %>% 
  filter(year == 2012) %>% 
  ungroup() %>% 
  select(-c(1, 3:5)) %>% 
  column_to_rownames(var = "site") %>% 
  select(c(1, 15:35)))

write_csv(temp, file.choose())

## biomass ----
# alternative code
# biomass calculated with all spp recorded
test <-
  biom1 %>% 
  pivot_wider(names_from = species, values_from = total, values_fill = 0) %>% 
  clean_names(case = "mixed") %>% 
  mutate(biomass_total = sum(c_across(Acanthurus_lineatus:Plectorhynchus_albovittatus))) %>% 
  relocate(biomass_total, .after = depth)

## spp abbrev ----
library(vegan)
colnames(df) <- make.cepnames(df) # ha ha...doesn't remotely work!
abbr <- make.cepnames(names(df[,-c(1:5)]))

colnames(df[-c(1:4)]) <- abbr

# these aren't working for some vile reason
spp <- c("Carcharhinus melanopterus", "Aulostomus chinensis", "Cephalopholis urodeta", "Lutjanus fulvus", "Lutjanus kasmira", "Macolor macularis", "Gnathodentex aureolineatus", "Monotaxis grandoculis", "Monotaxis grandoculis", "Pempheris oualensis")
(df <- tibble(obs = c(1:10), species = spp, count = sample(1:100, 10)))
dfw <- df %>% pivot_wider(names_from = species, values_from = count, values_fill = 0)
make.cepnames(names(dfw[-1]))
colnames(dfw[-1]) <- make.cepnames(names(dfw[-1]))

# raincloud plots
library(ggplot2)
library(ggdist)
library(tidyquant)
mpg %>% 
  filter(cyl %in% c(4,6,8)) %>% 
  ggplot(aes(x=factor(cyl), y=hwy, fill=factor(cyl))) +
  stat_halfeye(
    adjust = 0.5,
    justification = -0.2,
    .width = 0,
    point_color = NA
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.color = NA,
    alpha = 0.5
  ) +
  stat_dots(
    side = "left",
    justification = 1.1,
    binwidth = 0.25
  )
