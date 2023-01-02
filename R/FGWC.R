# Library yag dibutuhka
library(naspaclust)
library(rdist)

# Pentuan Parameter fgwc
param_fgwc <- c(kind='v',ncluster=3,m=2,distance='euclidean',order=3,
                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=0)

# parameter ABC
abc_param <- c(vi.dist='normal',npar=5,pso=FALSE,same=15,n.onlooker=5,limit=5)

# Parameter FPA
fpa_param <- c(vi.dist='normal',npar=5,same=15,p=0.7,gamma=1.2,lambda=1.5,ei.distr='logchaotic',chaos=3)

# Parameter GSA 
gsa_param <- c(vi.dist='normal',npar=5,same=15,G=1,vmax=0.7,new=FALSE) 

# Parameter HHO
hho_param <- c(vi.dist='normal',npar=5,same=15,algo='bairathi',a1=3,a2=1,a3=0.4)

# Parameter IFA
ifa_param <- c(vi.dist='uniform', ei.distr='logchaotic',
               fa.same=10, npar=15, par.no=3, par.dist='minkowski', 
               par.order=4, gamma=1, beta=1.5,
               alpha=1, chaos=4,update_type=4)

# Paraneter PSO 
pso_param <- c(vi.dist='uniform',npar=15,
               vmax=0.8, pso.same=10, c1=0.7, c2=0.6, type='chaotic',
               wmax=0.8,wmin=0.3,map=0.3)

# parameter TLBO 
tlbo_param <- c(vi.dist="uniform",nstud=10, tlbo.same=10,
                nselection=10,elitism=FALSE,n.elite=2)

# FGWC Standart
Res_fgwc <- fgwc(data=census2010, pop=census2010pop, distmat=census2010dist, algorithm = "classic", param_fgwc,1)

#ABC-FGWC
Res_abcfgwc <- fgwc(data=census2010, pop=census2010pop, distmat=census2010dist, algorithm = "abc", param_fgwc,abc_param)

#FPA-FGWC
Res_fpafgwc <- fgwc(data=census2010, pop=census2010pop, distmat=census2010dist, algorithm = "fpa", param_fgwc,abc_param)

#GSA-FGWC
Res_gsafgwc <- fgwc(data=census2010, pop=census2010pop, distmat=census2010dist, algorithm = "gsa", param_fgwc,abc_param)

#HHO-FGWC
Res_hhofgwc <- fgwc(data=census2010, pop=census2010pop, distmat=census2010dist, algorithm = "hho", param_fgwc,abc_param)

#IFA-FGWC
Res_ifafgwc <- fgwc(data=census2010, pop=census2010pop, distmat=census2010dist, algorithm = "ifa", param_fgwc,abc_param)

#PSO-FGWC
Res_psofgwc <- fgwc(data=census2010, pop=census2010pop, distmat=census2010dist, algorithm = "pso", param_fgwc,abc_param)

#TLBO-FGWC
Res_tlbofgwc <- fgwc(data=census2010, pop=census2010pop, distmat=census2010dist, algorithm = "tlbo", param_fgwc,abc_param)
