# Load Data Sovi
load("G:/My Drive/SKRIPSI/EXPERIMENT SKRIPSI/FGWC/DATA/sovi_data.rdata")
# Load Data Matrix Distance
load("G:/My Drive/SKRIPSI/EXPERIMENT SKRIPSI/FGWC/DATA/sovi_dist.rdata")
#Load Data Populasi
load("G:/My Drive/SKRIPSI/EXPERIMENT SKRIPSI/FGWC/DATA/sovi_pop.rdata")

#Mengubah data frame Sovi distance menjadi matrik distance
mat_dist= data.matrix(sovi_distance)

# Library yag dibutuhka
library(naspaclust)
library(rdist)

# Pentuan Parameter fgwc
param_fgwc <- c(kind='v',ncluster=4,m=2,distance='euclidean',order=3,
                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=10)

# parameter TLBO 
tlbo_param <- c(vi.dist="uniform",nstud=10, tlbo.same=10,
                nselection=10,elitism=FALSE,n.elite=2)

#TLBO-FGWC
Res_tlbofgwc <- fgwc(data=sovi_data, pop=Sovi_Pop, distmat=mat_dist,
                     algorithm = "tlbo", param_fgwc,tlbo_param)

# Menggabungkan hasil cluster algoritma FGWC dengan data
cluster.output <- cbind(data_kab ,Res_tlbofgwc$cluster)
write.csv(cluster.output, file = "TLBO Final Cluster.csv", row.names = TRUE)

