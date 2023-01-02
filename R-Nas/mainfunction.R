#' Fuzzy Geographicaly Weighted Clustering
#' @description Fuzzy clustering with addition of spatial configuration of membership matrix
#' @param data an object of data with d>1. Can be \code{matrix} or \code{data.frame}. If your data is univariate, bind it with \code{1} to get a 2 columns.
#' @param pop an n*1 vector contains population.
#' @param distmat an n*n distance matrix between regions.
#' @param algorithm algorithm used for FGWC
#' @param fgwc_param a vector that consists of FGWC parameter (see \code{\link{fgwcuv}} for parameter details)
#' @param opt_param a vector that consists of optimization algorithm parameter (see \code{\link{fgwcuv}} for parameter details)


#' @return an object of class \code{"fgwc"}.\cr
#' An \code{"fgwc"} object contains as follows:
#' \itemize{
#' \item \code{converg} - the process convergence of objective function
#' \item \code{f_obj} - objective function value
#' \item \code{membership} - membership matrix
#' \item \code{centroid} - centroid matrix
#' \item \code{validation} - validation indices (there are partition coefficient (\code{PC}), classification entropy (\code{CE}), 
#' SC index (\code{SC}), separation index (\code{SI}), Xie and Beni's index (\code{XB}), IFV index (\code{IFV}), and Kwon index (\code{Kwon}))
#' \item \code{max.iter} - Maximum iteration
#' \item \code{cluster} - the cluster of the data
#' \item \code{finaldata} - The final data (with the cluster)
#' \item \code{call} - the syntax called previously
#' \item \code{time} - computational time.
#' }

#' @details Fuzzy Geographically Weighted Clustering (FGWC) was developed by \insertCite{fgwc;textual}{naspaclust} by adding 
#' neighborhood effects and population to configure the membership matrix in Fuzzy C-Means. There are seven optimisation algorithms that currently
#' provided in this package, mainly from the \insertCite{yang2014;textual}{naspaclust}. The optimization algorithm uses the centroid as the parameter to be optimized. Here are the
#' algorithm that can be used:
#' \itemize{
#' \item \code{"classic"} - The classical algorithm of FGWC based on \insertCite{fgwc;textual}{naspaclust} for centroid optimisation 
#' and \insertCite{Runkler2006;textual}{naspaclust} for membership optimization.
#' \item \code{"abc"} - Optimization using artificial bee colony algorithm based on \insertCite{Karaboga2007;textual}{naspaclust} 
#' \insertCite{@see also @fgwcabc1 and @fgwcabc2 for FGWC implementation}{naspaclust}.
#' \item \code{"fpa"} - Optimization using flower pollination algorithm based on \insertCite{Yang2012}{naspaclust}.
#' \item \code{"gsa"} - Optimization using gravitational search algorithm based on \insertCite{rashedi2009;textual}{naspaclust} and 
#' \insertCite{Li2017gsa;textual}{naspaclust} \insertCite{@see also @fgwcgsa for FGWC implementation}{naspaclust}.
#' \item \code{"hho"} - Optimization using harris-hawk optimization with \code{"heidari"} \insertCite{Heidari2019}{naspaclust} (default).
#' and \code{"bairathi"} \insertCite{Bairathi2018}{naspaclust}.
#' \item \code{"ifa"} - Optimization using intelligent firefly algorithm based on \insertCite{Yang2009;textual}{naspaclust}, 
#' as well as the intelligent improvement by \insertCite{intfa;textual}{naspaclust} \insertCite{@see also @Nasution2020 for FGWC implementation}{naspaclust}.
#' \item \code{"pso"} - Optimization using particle swarm optimization based on \insertCite{Runkler2006;textual}{naspaclust} and 
#' \insertCite{Bansal2011;textual}{naspaclust} for inertia option \insertCite{@see also @fgwcpso; @putra2017; @Abdussamad for FGWC implementation}{naspaclust}.
#' \item \code{"tlbo"} - Optimization using teaching - learning based optimization based on \insertCite{Rao2012;textual}{naspaclust} and 
#' elitism improvement by \insertCite{Rao2012b;textual}{naspaclust}.
#' }
#' Furthermore, there are 10 distance that can be used to calculate the membership (see \code{\link{cdist}} for details).
#' the default parameter of FGWC (in case you do not want to tune anything) is \cr \code{
#' c(kind='u',ncluster=2,m=2,distance='euclidean',order=2,alpha=0.7,a=1,b=1,}\cr
#' \code{max.iter=500,error=1e-5,randomN=1)}.\cr
#' There is also a universal parameter to the optimization algorithm as well as the details. The default parameter
#' for the optimization algorithm is \cr
#' \code{c(vi.dist='uniform',npar=10,par.no=2,par.dist='euclidean',par.order=2,pso=TRUE,}\cr
#' \code{same=10,type='sim.annealing',ei.distr='normal',vmax=0.7,wmax=0.9,wmin=0.4,}\cr
#' \code{chaos=4,x0='F',map=0.7,ind=1,skew=0,sca=1)} \cr
#' If you do not define a certain parameter, the parameter will be set to its default value

#' @references
#' \insertAllCited{}

#' @seealso \code{\link{fgwcuv}}, \code{\link{abcfgwc}}, \code{\link{fpafgwc}},
#' \code{\link{gsafgwc}}, \code{\link{hhofgwc}}, \code{\link{ifafgwc}}, \code{\link{psofgwc}}, \code{\link{tlbofgwc}}
#' @examples
#' data('census2010')
#' data('census2010dist')
#' data('census2010pop')
#' # initiate parameter
#' param_fgwc <- c(kind='v',ncluster=3,m=2,distance='minkowski',order=3,
#'                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=10)
#' ## FGWC with classical algorithm
#' res1 <- fgwc(census2010,census2010pop,census2010dist,'classic',param_fgwc,1)
#' ## tune the ABC parameter
#' abc_param <- c(vi.dist='normal',npar=5,pso=FALSE,same=15,n.onlooker=5,limit=5) 
#' ## FGWC with ABC optimization algorithm
#' res2 <- fgwc(census2010,census2010pop,census2010dist,'abc',param_fgwc,abc_param) 
#' @export


#param_opt <- c(vi.dist='normal',npar=5,par.no=2,par.dist='euclidean',par.order=2,pso=FALSE,
#same=15,type='sim.annealing',ei.distr='normal',vmax=0.7,wmax=0.95,wmin=0.35,
#chaos=4,x0='F',map=0.7,ind=1,skew=0,sca=1,n.onlooker=5,limit=5) ## tune the optimization algorithm, add no. of onlooker bee as abc parameter


fgwc <- function(data,pop,distmat,algorithm='classic',fgwc_param,opt_param){
    if (!fgwc_param['kind']%in%c('u','v','NA')) kind <- NA else kind <- fgwc_param['kind']
    if (is.na(fgwc_param['ncluster'])|as.numeric(fgwc_param['ncluster'])<2) ncluster <- 2 else ncluster <- as.numeric(fgwc_param['ncluster'])
    if ((is.na(fgwc_param['m'])|as.numeric(fgwc_param['m'])<1)) m <- 2 else m <- as.numeric(fgwc_param['m'])
    if ((is.na(fgwc_param['distance']))) distance <- 'euclidean' else distance <- fgwc_param['distance']
    if ((is.na(fgwc_param['order'])|as.numeric(fgwc_param['order'])<=0)) order <- 2 else order <- as.numeric(fgwc_param['order'])
    if ((is.na(fgwc_param['alpha'])|as.numeric(fgwc_param['alpha'])>1|fgwc_param['alpha']<0)) alpha <- 0.7 else alpha <- as.numeric(fgwc_param['alpha'])
    if ((is.na(fgwc_param['a'])|as.numeric(fgwc_param['a'])<0)) a <- 1 else a <- as.numeric(fgwc_param['a'])
    if ((is.na(fgwc_param['b'])|as.numeric(fgwc_param['b'])<0)) b <- 1 else b <- as.numeric(fgwc_param['b'])
    if ((is.na(fgwc_param['max.iter'])|as.numeric(fgwc_param['max.iter'])<0)) max.iter <- 500 else max.iter <- as.numeric(fgwc_param['max.iter'])
    if ((is.na(fgwc_param['error'])|as.numeric(fgwc_param['error'])<0)) error <- 1e-5 else error <- as.numeric(fgwc_param['error'])
    if ((is.na(fgwc_param['randomN'])|as.numeric(fgwc_param['randomN'])<0)) randomN <- 1 else randomN <- as.numeric(fgwc_param['randomN'])
    
    if(algorithm!='classic'){
        if(is.na(opt_param['vi.dist'])) vi.dist <- 'uniform' else vi.dist <- opt_param['vi.dist']
        if(is.na(opt_param['npar'])|as.numeric(opt_param['npar'])<0) npar <- 10 else npar <- as.numeric(opt_param['npar'])
        if(is.na(opt_param['par.no'])|as.numeric(opt_param['par.no'])<0) par.no <- 2 else par.no <- as.numeric(opt_param['par.no'])
        if(is.na(opt_param['par.dist'])) par.dist <- 'euclidean' else par.dist <- opt_param['par.dist']
        if(is.na(opt_param['par.order'])|as.numeric(opt_param['par.order'])<0) par.order <- 2 else par.order <- as.numeric(opt_param['par.order'])
        if(is.na(opt_param['pso'])) pso <- TRUE else pso <- as.logical(opt_param['pso'])
        if(is.na(opt_param['same'])|as.numeric(opt_param['same'])<0) same <- 10 else same <- as.numeric(opt_param['same'])
        if(is.na(opt_param['type'])) type <- 'sim.annealing' else type <- opt_param['type']
        if(is.na(opt_param['ei.distr'])) ei.distr <- 'normal' else ei.distr <- opt_param['ei.distr']
        if(is.na(opt_param['vmax'])|as.numeric(opt_param['vmax'])<0) vmax <- 0.7 else vmax <- as.numeric(opt_param['vmax'])
        if(is.na(opt_param['wmax'])|as.numeric(opt_param['wmax'])<0) wmax <- 0.9 else wmax <- as.numeric(opt_param['wmax'])
        if(is.na(opt_param['wmin'])|as.numeric(opt_param['wmin'])<0) wmin <- 0.4 else wmin <- as.numeric(opt_param['wmin'])
        if(is.na(opt_param['chaos'])|as.numeric(opt_param['chaos'])<0) chaos <- 4 else chaos <- as.numeric(opt_param['chaos'])
        if(is.na(opt_param['x0'])) x0 <- 'F' else x0 <- opt_param['x0']
        if(is.na(opt_param['map'])|as.numeric(opt_param['map'])<0) map <- 0.7 else map <- as.numeric(opt_param['map'])
        if(is.na(opt_param['ind'])|as.numeric(opt_param['ind'])<0) ind <- 1 else ind <- as.numeric(opt_param['ind'])
        if(is.na(opt_param['skew'])|as.numeric(opt_param['skew'])<0) skew <- 0 else skew <- as.numeric(opt_param['skew'])
        if(is.na(opt_param['sca'])|as.numeric(opt_param['sca'])<0) sca <- 1 else sca <- as.numeric(opt_param['sca'])
    }

    if (algorithm=='classic'){
        return(fgwcuv(data, pop, distmat, kind=kind,ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, max.iter=max.iter, error=error,randomN=randomN, uij=NA, vi=NA))
    }
    else if(algorithm=='abc'){
        opt_param <- get_param_abc(opt_param)
        return(abcfgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter, randomN=randomN, vi.dist=vi.dist, 
                      nfood=npar, n.onlooker=as.numeric(opt_param['n.onlooker']), 
                      limit=as.numeric(opt_param['limit']), pso=pso, abc.same=same))
    }
    else if(algorithm=='fpa'){
        opt_param <- get_param_fpa(opt_param)
        return(fpafgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter, randomN=randomN,
                      vi.dist=vi.dist, nflow=npar, p=as.numeric(opt_param['p']), gamma=as.numeric(opt_param['gamma']), 
                      lambda=as.numeric(opt_param['lambda']), delta=as.numeric(opt_param['delta']),
                      ei.distr=ei.distr,flow.same=same,r=chaos,m.chaotic=map,skew=skew,sca=sca))
    }
    else if(algorithm=='gsa'){
        opt_param <- get_param_gsa(opt_param)
        print(opt_param)
        return(gsafgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter,randomN=randomN,
                      vi.dist=vi.dist,npar=npar,par.no=par.no,par.dist=par.dist, par.order=par.order,
                      gsa.same=same, G=as.numeric(opt_param['G']), vmax=vmax, new=as.logical(opt_param['new'])))
    }
    else if(algorithm=='hho'){
        opt_param <- get_param_hho(opt_param)
        print(as.numeric(opt_param[c('a1','a2','a3')]))
        return(hhofgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter, randomN=randomN,
                      vi.dist=vi.dist,nhh=npar,hh.alg=opt_param['algo'],
                      A=as.numeric(opt_param[c('a1','a2','a3')]),p=as.numeric(opt_param['p']),
                      hh.same=same,levy.beta=as.numeric(opt_param['beta']),update.type=as.numeric(opt_param['update_type'])))
    }
    else if(algorithm=='ifa'){
        opt_param <- get_param_ifa(opt_param)
        print(opt_param)
        return(ifafgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter, randomN=randomN,
                      vi.dist=vi.dist, ei.distr=ei.distr,fa.same=same, nfly=npar, ffly.no=par.no, 
                      ffly.dist=par.dist, ffly.order=par.order, gamma=as.numeric(opt_param['gamma']), 
                      ffly.beta=as.numeric(opt_param['beta']),ffly.alpha=as.numeric(opt_param['alpha']),
                      r.chaotic=chaos,m.chaotic=map,ind.levy=ind,skew.levy=skew,
                      scale.levy=sca,ffly.alpha.type=as.numeric(opt_param['update_type'])))
    }
    else if(algorithm=='pso'){
        opt_param <- get_param_pso(opt_param)
        return(psofgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter, randomN=randomN,
                      vi.dist=vi.dist,npar=npar,vmax=vmax, pso.same=same, c1=as.numeric(opt_param['c1']),
                      c2=as.numeric(opt_param['c2']), w.inert=type,
                      wmax=wmax,wmin=wmin,map=map))
    }
    else if(algorithm=='tlbo'){
        opt_param <- get_param_tlbo(opt_param)
        return(tlbofgwc(data, pop, distmat, ncluster=ncluster, m=m, distance=distance, order=order,
                      alpha=alpha, a=a, b=b, error=error,max.iter=max.iter, randomN=randomN,
                      vi.dist=vi.dist,nstud=npar,tlbo.same=same,
                      nselection=as.numeric(opt_param['nselection']),
                      elitism=as.logical(opt_param['elitism']),n.elite=as.numeric(opt_param['n.elite'])))
    }
    
}

get_param_abc <- function(param){
    paramx <- c()
    if(is.na(param['n.onlooker'])|param['n.onlooker']<0) paramx['n.onlooker'] <- 5 else paramx['n.onlooker'] <- param['n.onlooker']
    if(is.na(param['limit'])|param['limit']<0) paramx['limit'] <- 4 else paramx['limit'] <- param['limit']
    return(paramx)
}

get_param_fpa <- function(param){
    paramx <- c()
    if(is.na(param['p'])|param['p']<0|param['p']>1) paramx['p'] <- 0.8 else paramx['p'] <- param['p']
    if(is.na(param['gamma'])|param['gamma']<0) paramx['gamma'] <- 1 else paramx['gamma'] <- param['gamma']
    if(is.na(param['lambda'])|param['lambda']<0) paramx['lambda'] <- 1.5 else paramx['lambda'] <- param['lambda'] 
    if(is.na(param['delta'])|param['delta']<0) paramx['delta'] <- 0 else paramx['delta'] <- param['delta']
    return(paramx)
}

get_param_gsa <- function(param){
    paramx <- c()
    if(is.na(param['G'])|param['G']<0) paramx['G'] <- 1 else paramx['G'] <- param['G']
    if(is.na(param['new'])|!param['new']%in%c(0,1)) paramx['new'] <- TRUE else paramx['new'] <- as.logical(param['new'])
    return(paramx)
}

get_param_hho <- function(param){
    paramx <- c()
    if(is.na(param['algo'])|!param['algo']%in%c('heidari','bairathi')) paramx['algo'] <- 'heidari' else paramx['algo'] <- param['algo']
    if(is.na(param['a1'])|param['a1']<0) paramx['a1'] <- 2 else paramx['a1'] <- param['a1']
    if(is.na(param['a2'])|param['a2']<0) paramx['a2'] <- 1 else paramx['a2'] <- param['a2']
    if(is.na(param['a3'])|param['a3']<0) paramx['a3'] <- 0.5 else paramx['a3'] <- param['a3']
    if(is.na(param['p'])|param['p']<0|param['p']>1) paramx['p'] <- 0.5 else paramx['p'] <- param['p']
    if(is.na(param['beta'])|param['beta']<0) paramx['beta'] <- 1.5 else paramx['beta'] <- param['beta']
    if(is.na(param['update_type'])|param['update_type']<5) paramx['update_type'] <- 5 else paramx['update_type'] <- param['update_type']
    return(paramx)
}

get_param_ifa <- function(param){
    paramx <- c()
    if(is.na(param['gamma'])|param['gamma']<0) paramx['gamma'] <- 1 else paramx['gamma'] <- param['gamma']
    if(is.na(param['beta'])|param['beta']<0) paramx['beta'] <- 1 else paramx['beta'] <- param['beta']
    if(is.na(param['alpha'])|param['alpha']<0) paramx['alpha'] <- 1 else paramx['alpha'] <- param['alpha']
    if(is.na(param['update_type'])|param['update_type']<5) paramx['update_type'] <- 4 else paramx['update_type'] <- param['update_type']
    return(paramx)
}

get_param_pso <- function(param){
    paramx <- c()
    if(is.na(param['c1'])|param['c1']<0) paramx['c1'] <- 0.49 else paramx['c1'] <- param['c1']
    if(is.na(param['c2'])|param['c2']<0) paramx['c2'] <- 0.49 else paramx['c2'] <- param['c2']
    return(paramx)
}

get_param_tlbo <- function(param){
    paramx <- c()
    if(is.na(param['nselection'])|param['nselection']<0) paramx['nselection'] <- 10 else paramx['nselection'] <- param['nselection']
    if(is.na(param['elitism'])|param['elitism']<0) paramx['elitism'] <- F else paramx['elitism'] <- param['elitism']
    if(is.na(param['n.elite'])|param['n.elite']<0) paramx['n.elite'] <- 2 else paramx['n.elite'] <- param['n.elite']
    return(paramx)
}