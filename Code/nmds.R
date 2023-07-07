### NMDS ###
#Libraries
library(tidyverse)
library(vegan)

#Data manipulation
data <- read.csv("https://raw.githubusercontent.com/cmlglvz/Remember-R/main/Data/abundancia_taxonomia.csv", header = TRUE, sep = ";", skip = 0, row.names = 1)
short <- dplyr::filter(data,Total >= 20)
abund <- short[,c(3:36)] %>% t() %>% as.data.frame()
colnames(abund) <- short$ASV

#Abundance normalization
hell <- vegan::decostand(abund, method = "hellinger") #save file

#Dissimilarity indices
dis.indx <- vegan::vegdist(x = hell, method = "bray", diag = TRUE) #save file

#NMDS
set.seed(19910420) #reproducible results
nmds <- metaMDS(comm = hell, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)nmds <- metaMDS(comm = hell, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)
nmds <- vegan::metaMDS(comm = hell, #or dissimilarity indices
                       distance = "bray", 
                       k = 2, 
                       trymax = 100, 
                       autotransform = FALSE #we already normalized the data
                       )
str(nmds)
scrs <- vegan::scores(nmds)
nmds$points
plot(nmds)

#Plot as you  prefer