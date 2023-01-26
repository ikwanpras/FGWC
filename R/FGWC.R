# Load Data Sovi
load("G:/My Drive/SKRIPSI/EXPERIMENT SKRIPSI/FGWC/DATA/census2010.rdata")
# Load Data Matrix Distance
load("G:/My Drive/SKRIPSI/EXPERIMENT SKRIPSI/FGWC/DATA/census2010dist.rdata")
#Load Data Populasi
load("G:/My Drive/SKRIPSI/EXPERIMENT SKRIPSI/FGWC/DATA/census2010pop.rdata")

#Plot Disatance
library(ggplot2)
library(factoextra)
distance <- get_dist(census2010)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# Library yag dibutuhka
library(naspaclust)
library(rdist)

# Pentuan Parameter fgwc
param_fgwc <- c(kind='v',ncluster=3,m=2,distance='euclidean',order=3,
                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=0)

# FGWC Standart
Res_fgwc <- fgwc(data=census2010, pop=census2010pop, distmat=census2010dist, algorithm = "classic", param_fgwc,1)

fviz_cluster(Res_fgwc, data = census2010)
