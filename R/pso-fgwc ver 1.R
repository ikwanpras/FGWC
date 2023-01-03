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

# Parameter PSO 
pso_param <- c(vi.dist='uniform',npar=15,
               vmax=0.8, pso.same=10, c1=0.7, c2=0.6, type='chaotic',
               wmax=0.8,wmin=0.3,map=0.3)

#PSO-FGWC
Res_psofgwc <- fgwc(data=sovi_data, pop=Sovi_Pop, distmat=mat_dist, 
                    algorithm = "pso", param_fgwc,pso_param)
