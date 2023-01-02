# mengubah data populasi menjadi vektor
Sovi_Pop <- sovi_data_pop_skripsi$...1
print(Sovi_Pop)

aa <- as.matrix(sovi_distance)

# Library yag dibutuhka
library(naspaclust)
library(rdist)

# Pentuan Parameter fgwc
param_fgwc <- c(kind='v',ncluster=3,m=2,distance='euclidean',order=3,
                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=0)

# FGWC Standart
Res_fgwc <- fgwc(data=sovi_data_skripsi, pop=aa, distmat=sovi_distance, algorithm = "classic", param_fgwc,1)
