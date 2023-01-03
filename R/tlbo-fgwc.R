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

# parameter TLBO 
tlbo_param <- c(vi.dist="uniform",nstud=10, tlbo.same=10,
                nselection=10,elitism=FALSE,n.elite=2)

#TLBO-FGWC
Res_tlbofgwc <- fgwc(data=census2010, pop=census2010pop, distmat=census2010dist, algorithm = "tlbo", param_fgwc,abc_param)
