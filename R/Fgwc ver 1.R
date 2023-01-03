# Load Data Sovi
load("D:/My Drive/deden-RMarkdown & Github/FGWC/Data/sovi_data.rdata")
# Load Data Matrix Distance
load("D:/My Drive/deden-RMarkdown & Github/FGWC/Data/sovi_dist.rdata")
#Load Data Populasi
load("D:/My Drive/deden-RMarkdown & Github/FGWC/Data/sovi_pop.rdata")

#Mengubah data frame Sovi distance menjadi matrik distance
mat_dist= data.matrix(sovi_distance)


# Library yag dibutuhkan
library(naspaclust)
library(rdist)

# Pentuan Parameter fgwc
param_fgwc <- c(kind='v',ncluster=4,m=2,distance='euclidean',order=3,
                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=0)

# FGWC Standart
Res_fgwc <- fgwc(data=sovi_data, pop=Sovi_Pop, distmat=mat_dist, algorithm = "classic", param_fgwc,1)

# Menggabungkan hasil cluster algoritma FGWC dengan data
library(xlsx)
cluster.output <- cbind(data_kab ,Res_fgwc$cluster)
write.xlsx(cluster.output, file = "Final Cluster.xlsx", row.names = TRUE)
