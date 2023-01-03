# Load Data
load("D:/My Drive/deden-RMarkdown & Github/FGWC/Data/census2010.rdata")
load("D:/My Drive/deden-RMarkdown & Github/FGWC/Data/census2010dist.rdata")
load("D:/My Drive/deden-RMarkdown & Github/FGWC/Data/census2010pop.rdata")

# Library yag dibutuhka
library(naspaclust)
library(rdist)

# Pentuan Parameter fgwc
param_fgwc <- c(kind='v',ncluster=3,m=2,distance='euclidean',order=3,
                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=0)

# Parameter IFA
ifa_param <- c(vi.dist='uniform', ei.distr='logchaotic',
               fa.same=10, npar=15, par.no=3, par.dist='minkowski', 
               par.order=4, gamma=1, beta=1.5,
               alpha=1, chaos=4,update_type=4)

#IFA-FGWC
Res_ifafgwc <- fgwc(data=census2010, pop=census2010pop, distmat=census2010dist, 
                algorithm = "ifa", param_fgwc,abc_param)
