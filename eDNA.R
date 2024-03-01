# Peter Nelson
# Institute of Marine Sciences
# University of California, Santa Cruz
# pnelson1@ucsc.edu

# Ulithi fishes
# read in and wrangle 2024 eDNA data
# created 22 February 2024

# get started ----
library(tidyverse)
library(readxl)
library(janitor)

temp <- read_excel("data/Summary_Fishes_eDNA_Ulithi_2024_PN.xlsx", sheet = "eDNA") %>% 
  mutate(found = factor(found),
         rank = factor(rank),
         tax_ID = as.integer(tax_ID),
         group = factor(group) # summary data includes only fish
         )
glimpse(temp)
summary(temp)

fish <- temp %>% 
  filter(found == "y", # limits those to spp previously observed by Giacomo
         rank == "species") %>% # observations limited to those at the species level
  select(-c(found, tax_ID, common, group))

fish1 <- temp %>% 
  filter(found != "x",
         rank == "species",
         name != "Coryphaena hippurus") %>% 
  select(-c(found, tax_ID, common, group))

# community ecology -----
library(vegan)
df <- t(fish[,3:16])
colnames(df) <- fish$name
sites <- clean_names(df, case = "none") %>% print()

df1 <- t(fish1[,3:16])
colnames(df1) <- fish1$name
sites1 <- clean_names(df1, case = "none") %>% print()

## diversity -----
# analyses with more conservative data set
# limited to fishes identified to species AND previously seen by Giacomo
fdiv_simpson <- diversity(sites, index = "simpson")
fdiv_shannon <- diversity(sites, index = "shannon")
barplot(fdiv_simpson, las = 2)
barplot(fdiv_shannon, las = 2)
par(mfrow = c(1, 2)) 
hist(fdiv_simpson)
hist(fdiv_shannon)

## hierarchical clustering ----
# note that "bray" is the site dissimilarity using Bray-Curtis; "gower" uses Gower...
# everything (at present) uses Bray-Curtis, but it'd be worth looking at the Gower results

# compute pair-wise dissimilarity measurements
bray <- vegdist(sites, "bray") 
gower <- vegdist(sites, "gower")
cao <- vegdist(sites, "cao")
hist(bray)
hist(gower)
hist(cao)

csin <- hclust(bray, method = "single") # single linkage clustering
ccom <- hclust(bray, method = "complete") # complete linkage
caver <- hclust(bray, method = "aver") # average linkage

par(mfrow = c(3, 1))
plot(csin, hang = -1) # hang forces branches down to base line
rect.hclust(csin, 3) # identify classes or groups, here n=3
plot(ccom, hang = -1)
rect.hclust(ccom, 3)
plot(caver, hang = -1)
rect.hclust(caver, 3)

# repeat using Gower
csin_g <- hclust(gower, method = "single") # single linkage clustering
ccom_g <- hclust(gower, method = "complete") # complete linkage
caver_g <- hclust(gower, method = "aver") # average linkage

# what's the best number of groups to identify? not necessarily 3...
plot(csin_g, hang = -1) # hang forces branches down to base line
rect.hclust(csin_g, 3) # identify classes or groups, here n=3
plot(ccom_g, hang = -1)
rect.hclust(ccom_g, 3)
plot(caver_g, hang = -1)
rect.hclust(caver_g, 3)

par(mfrow = c(1, 1))

# force clustering into larger number of classes
plot(csin, hang = -1)
rect.hclust(csin, 5)

plot(csin_g, hang = -1)
rect.hclust(csin_g, 5)

# how do spp differ in classification?
vegemite(sites, caver, scale = "Hill")

## ordination -----
ord <- cmdscale(bray)
ordiplot(ord) # suggests to me we want...4 or 5 classes (ignore warning)
plot(caver, hang = -1)
rect.hclust(caver, 5) # ugh! maybe not!
c1 <- cutree(caver, 3) # extract classification at n=3
table(c1) # shows number of obs in each class (cluster)

ord_g <- cmdscale(gower)
ordiplot(ord_g) # suggests to me we want...4 or 5 classes (ignore warning)
plot(caver_g, hang = -1)
rect.hclust(caver_g, 5) # ugh! maybe not!
c1g <- cutree(caver_g, 3) # extract classification at n=3
table(c1g) # shows number of obs in each class (cluster)

# confusion matrix: rows give first classification, columns the second
# if the classifications match there is no "confusion"...each row & each col has only one non-zero entry
table(c1, cutree(csin, 3))
table(c1, cutree(ccom, 3))

ordiplot(ordihull(ord, c1)) # seems to remove some sites...which? why?
ordiplot(ord, dis = "sites")
ordihull(ord,cutree(caver,5))
ordicluster(ord,csin)

ordiplot(ord, dis = "sites")
ordicluster(ord, caver, prune = 2) # prune top level fusions to highlight clustering

# with Gower...
table(c1g, cutree(csin_g, 3))
table(c1g, cutree(ccom_g, 3))

ordiplot(ordihull(ord_g, c1g)) # seems to remove some sites...which? why?
ordiplot(ord_g, dis = "sites")
ordihull(ord_g,cutree(caver_g,5))
ordicluster(ord_g,csin_g)

ordiplot(ord_g, dis = "sites")
ordicluster(ord_g, caver_g, prune = 2) # prune top level fusions to highlight clustering

## optimized clustering -----
# K-means clustering looks for an optimal solution given a specified number of classes
# ideally, that might in this case be based on site characteristics (eg in/out of lagoon, at or away from human settlement, part of/not part of atoll, etc)
# works in euclidean metric (not meaningful for community data)--requires standardization;
# Jari recs Hellinger transformation using decostand...

ckm <- kmeans(decostand(sites, "hell"), 3)
ckm$cluster # possibly makes much better ecological sense than prior clusters!

ordiplot(ord, dis = "sites", type = "text")
ordihull(ord, ckm$cluster, col = "blue")

# again, with Gower
ckm <- kmeans(decostand(sites, "hell"), 3)
ckm$cluster # possibly makes much better ecological sense than prior clusters!

ordiplot(ord_g, dis = "sites", type = "text")
ordihull(ord_g, ckm$cluster, col = "green")

# NOT impressed with Gower here!

## ordination w NMDS ----

# w complete eDNA dataset -----
aggr <- read_excel("data/Summary_Fishes_eDNA_Ulithi_2024_PN.xlsx", sheet = "aggregated") %>% 
  mutate(name = ScientificName,
         rank = factor(Rank),
         tax_ID = as.integer(TaxID),
         common = CommonName,
         group = factor(Group), # summary data includes only fish
         .keep = "unused") %>% 
  select(-c(tax_ID)) %>% 
  relocate(., group, rank, name, common)

glimpse(aggr)
summary(aggr)

# community ecology -----
library(vegan)

taggr <- t(aggr[,5:18])
# colnames(taggr) <- make.cepnames(as.matrix(aggr[, 3]))

temp <- t(aggr[1:5,5:10])
colnames(temp) <- make.cepnames(as.matrix(aggr[1:5, 3]))
temp

## hierarchical clustering ----
# note that "bray" is the site dissimilarity using Bray-Curtis; "gower" uses Gower...
# everything (at present) uses Bray-Curtis, but it'd be worth looking at the Gower results
# compute pair-wise dissimilarity measurements
bray <- vegdist(taggr, "bray") 
gower <- vegdist(taggr, "gower")
cao <- vegdist(taggr, "cao")
hist(bray)
hist(gower)
hist(cao)

csin <- hclust(bray, method = "single") # single linkage clustering
ccom <- hclust(bray, method = "complete") # complete linkage
caver <- hclust(bray, method = "aver") # average linkage
cward <- hclust(bray, method = "ward.D2") # Ward's 1963 criterion (dissim squared before updating cluster)

par(mfrow = c(4, 1))
plot(csin, hang = -1) # hang forces branches down to base line
plot(ccom, hang = -1)
plot(caver, hang = -1)
plot(cward, hang = -1)
# rect.hclust(caver, 3)

# repeat using Gower
csin_g <- hclust(gower, method = "single") # single linkage clustering
ccom_g <- hclust(gower, method = "complete") # complete linkage
caver_g <- hclust(gower, method = "aver") # average linkage
cward_g <- hclust(gower, method = "ward.D2")

# what's the best number of groups to identify? not necessarily 3...
plot(csin_g, hang = -1) # hang forces branches down to base line
plot(ccom_g, hang = -1)
plot(caver_g, hang = -1)
plot(cward_g, hang = -1)
# rect.hclust(caver_g, 3)

par(mfrow = c(1, 1))

# force clustering into larger number of classes
plot(csin, hang = -1)
rect.hclust(csin, 5)

plot(csin_g, hang = -1)
rect.hclust(csin_g, 5)
