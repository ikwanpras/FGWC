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
                alpha=0.7,a=1,b=1,max.iter=1000,error=1e-5,randomN=0)

# Parameter IFA
ifa_param <- c(vi.dist='uniform', ei.distr='logchaotic',
               fa.same=10, npar=15, par.no=3, par.dist='minkowski', 
               par.order=4, gamma=1, beta=1,
               alpha=1, chaos=4,update_type=4)

#IFA-FGWC
Res_ifafgwc <- fgwc(data=sovi_data, pop=Sovi_Pop, distmat=mat_dist, 
                    algorithm = "ifa", param_fgwc,ifa_param)

# Menggabungkan hasil cluster algoritma FGWC dengan data
library(xlsx)
cluster.output <- cbind(data_kab ,Res_ifafgwc$cluster)
write.xlsx(cluster.output, file = "Final Cluster.xlsx", row.names = TRUE)
