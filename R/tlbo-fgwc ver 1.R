# Load Data Sovi
load("D:/My Drive/deden-RMarkdown & Github/FGWC/Data/sovi_data.rdata")
# Load Data Matrix Distance
load("D:/My Drive/deden-RMarkdown & Github/FGWC/Data/sovi_dist.rdata")
#Load Data Populasi
load("D:/My Drive/deden-RMarkdown & Github/FGWC/Data/sovi_pop.rdata")

#Mengubah data frame Sovi distance menjadi matrik distance
mat_dist= data.matrix(sovi_distance)

# Library yag dibutuhka
library(naspaclust)
library(rdist)

# Pentuan Parameter fgwc
param_fgwc <- c(kind='v',ncluster=4,m=2,distance='euclidean',order=3,
                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=0)

# parameter TLBO 
tlbo_param <- c(vi.dist="uniform",nstud=10, tlbo.same=10,
                nselection=10,elitism=FALSE,n.elite=2)

#TLBO-FGWC
Res_tlbofgwc <- fgwc(data=sovi_data, pop=Sovi_Pop, distmat=mat_dist, algorithm = "tlbo", param_fgwc,tlbo_param)
