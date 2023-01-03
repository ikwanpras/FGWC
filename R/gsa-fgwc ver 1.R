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

# Parameter GSA 
gsa_param <- c(vi.dist='normal',npar=5,same=15,G=1,vmax=0.7,new=FALSE) 

#GSA-FGWC
Res_gsafgwc <- fgwc(data=sovi_data, pop=Sovi_Pop, distmat=mat_dist, 
                    algorithm = "gsa", param_fgwc,gsa_param)
