#' Fuzzy Geographicaly Weighted Clustering with Flower Pollination Algorithm
#' @description Fuzzy clustering with addition of spatial configuration of membership matrix with centroid optimization using Flower Pollination Algorithm
#' @param data an object of data with d>1. Can be \code{matrix} or \code{data.frame}. If your data is univariate, bind it with \code{1} to get a 2 columns.
#' @param pop an n*1 vector contains population.
#' @param distmat an n*n distance matrix between regions.
#' @param ncluster an integer. The number of clusters.
#' @param m degree of fuzziness or fuzzifier. Default is 2.
#' @param distance the distance metric between data and centroid, the default is euclidean, see \code{\link{cdist}} for details.
#' @param order, minkowski order. default is 2.
#' @param alpha the old membership effect with [0,1], if \code{alpha} equals 1, it will be same as fuzzy C-Means, if 0, it equals to neighborhood effect.
#' @param a spatial magnitude of distance. Default is 1.
#' @param b spatial magnitude of population. Default is 1.
#' @param max.iter maximum iteration. Default is 500.
#' @param error error tolerance. Default is 1e-5.
#' @param randomN random seed for initialisation (if uij or vi is NA). Default is 0.
#' @param vi.dist a string of centroid population distribution between \code{"uniform"} (default) and \code{"normal"}. Can be defined as \code{vi.dist=} in \code{opt_param}.
#' @param nflow number of flowers population. Can be defined as \code{npar=} in \code{opt_param}. Default is 10.
#' @param p switch probability between global and local pollination, Can be defined as \code{p} in \code{opt_param}. default is 0.8.
#' @param gamma Step size scaling factor. Can be defined as \code{gamma} in \code{opt_param}. Default is 1.
#' @param lambda Levy flights index parameter between [0,2]. Can be defined as \code{lambda} in \code{opt_param}. Default is 1.5.
#' @param delta Levi flights shift. Can be defined as \code{delta} in \code{opt_param}. Default is 0.
#' @param ei.distr distribution of random walk parameter. Can be defined as \code{ei.distr} in \code{opt_param}.
#' @param flow.same number of consecutive unchange to stop the iteration. Can be defined as \code{same=} in \code{opt_param}.
#' @param r weight in logistic chaotic between [0,4]. Can be used when \code{ei.distr='logchaotic'}. Can be defined as \code{chaos} in \code{opt_param}.
#' @param m.chaotic mapping parameter in kent chaotic between [0,1]. Can be used when \code{ei.distr='kentchaotic'}. Can be defined as \code{map} in \code{opt_param}.
#' @param skew Levy distribution skewness for random walk. Can be used when \code{ei.distr='levy'}. Can be defined as \code{skew} in \code{opt_param}.
#' @param sca Levy distribution scale for random walk. Can be used when \code{ei.distr='levy'}. Can be defined as \code{sca} in \code{opt_param}.

#' @return an object of class \code{"fgwc"}.\cr
#' An \code{"fgwc"} object contains as follows:
#' \itemize{
#' \item \code{converg} - the process convergence of objective function
#' \item \code{f_obj} - objective function value
#' \item \code{membership} - membership matrix
#' \item \code{centroid} - centroid matrix
#' \item \code{validation} - validation indices (there are partition coefficient (\code{PC}), classification entropy (\code{CE}), 
#' SC index (\code{SC}), separation index (\code{SI}), Xie and Beni's index (\code{XB}), IFV index (\code{IFV}), and Kwon index (Kwon))
#' \item \code{max.iter} - Maximum iteration
#' \item \code{cluster} - the cluster of the data
#' \item \code{finaldata} - The final data (with the cluster)
#' \item \code{call} - the syntax called previously
#' \item \code{time} - computational time.
#' }

#' @details Fuzzy Geographically Weighted Clustering (FGWC) was developed by \insertCite{fgwc;textual}{naspaclust} by adding 
#' neighborhood effects and population to configure the membership matrix in Fuzzy C-Means. Furthermore,
#' the Flower Pollination Algorithm was developed by \insertCite{Yang2012;textual}{naspaclust} in order to get a more optimal
#' solution of a certain complex function.

#' @seealso \code{\link{fpafgwc}} \code{\link{gsafgwc}}
#' @examples
#' data('census2010')
#' data('census2010dist')
#' data('census2010pop')
#' # First way
#' res1 <- fpafgwc(census2010,census2010pop,census2010dist,3,2,'euclidean',4,nflow=10)
#' # Second way
#' # initiate parameter
#' param_fgwc <- c(kind='v',ncluster=3,m=2,distance='minkowski',order=3,
#'                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=10)
#' ## tune the FPA parameter
#' fpa_param <- c(vi.dist='normal',npar=5,same=15,p=0.7,
#'                gamma=1.2,lambda=1.5,ei.distr='logchaotic',chaos=3) 
#' ##FGWC with FPA
#' res2 <- fgwc(census2010,census2010pop,census2010dist,'fpa',param_fgwc,fpa_param)

#' @export
#' @import beepr
#' @import stabledist
#' @import rdist
#' @import stats

fpafgwc <- function(data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1,
					error=1e-5, max.iter=100, randomN=0, vi.dist="uniform", nflow=10, p=0.8, gamma=1, lambda=1.5, delta=0,
          ei.distr='normal', flow.same=10,r=4,m.chaotic=0.7,skew=0,sca=1){
  # require(stabledist)
  randomnn <- randomN
  ptm<-proc.time()
  n <- nrow(data)
  d <- ncol(data)
  iter=0
  beta <- 1-alpha
  same=0
  data <- as.matrix(data)
  if (alpha ==1) {
    pop <- rep(1,n)
    distmat <- matrix(1,n,n)
  }
  datax <- data
  pop <- matrix(pop,ncol=1)
  mi.mj <- pop%*%t(pop)

  flow <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster,
                                m, alpha, a, b, randomN, nflow)
  flow.swarm <- flow$centroid
  flow.other <- flow$membership
  flow.fit <- flow$I

  flow.finalpos <- flow$centroid[[which.min(flow.fit)]]
  flow.finalpos.other <- flow$membership[[which.min(flow.fit)]]
  flow.fit.finalbest <- flow$I[[which.min(flow.fit)]]
  conv <- c(flow.fit[which.min(flow.fit)])
  repeat{
  	minmax <- c(which.min(flow.fit)[1],which.max(flow.fit)[1])
  	best <- minmax[1]
  	worst <- minmax[2]
  	pollen <- flow.finalpos
    flow.swarm <- pollination(flow.swarm,p,pollen,gamma,lambda,delta,randomN,ei.distr,r,m.chaotic,skew,sca)
    flow.other <- lapply(1:nflow, function(x) uij(data,flow.swarm[[x]],m,distance,order))
  	flow.other <- lapply(1:nflow, function(x) renew_uij(data,flow.other[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
  	flow.swarm <- lapply(1:nflow, function(x) vi(data,flow.other[[x]],m))
  	flow.fit <- sapply(1:nflow, function(x) jfgwcv(data,flow.swarm[[x]],m,distance,order))
    best <- which(flow.fit==min(flow.fit))[1]
    flow.curbest <- flow.swarm[[best]]
    flow.curbest.other <- flow.other[[best]]
    flow.fit.curbest <- flow.fit[best]
    conv <- c(conv,flow.fit.finalbest)
    iter <- iter+1
    if (abs(conv[iter+1]-conv[iter])<error) same <- same+1
    else same <- 0
    if (flow.fit.curbest<=flow.fit.finalbest) {
      flow.finalpos <- flow.curbest
      flow.finalpos.other <- flow.curbest.other
      flow.fit.finalbest <- flow.fit.curbest
    }
    randomN <- randomN+nflow
    if (iter==max.iter || same==flow.same) break
  }
  finaldata=determine_cluster(datax,flow.finalpos.other)
  cluster=finaldata[,ncol(finaldata)]
  print(c(order, ncluster,m, randomN))
  fpa <- list("converg"=conv,"f_obj"=jfgwcv(data,flow.finalpos,m,distance,order),"membership"=flow.finalpos.other,"centroid"=flow.finalpos,
              "validation"=index_fgwc(data,cluster,flow.finalpos.other,flow.finalpos,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"iteration"=iter,"same"=same,"time"=proc.time()-ptm)
  class(fpa) <- 'fgwc'
  return(fpa)
}

pollination <- function(flow,p,pollen,gamma,lambda,delta,seed,ei.distr,r,m,skew,sca){
  set.seed(seed<-seed+10)
  rand <- runif(length(flow))
  dd <- dim(pollen)
  return(lapply(1:length(flow),function(x){
    if(rand[x]<p){ ##global pollination
      set.seed(seed<-seed+1)
      flow[[x]]+gamma*matrix(rstable(dd[1]*dd[2],lambda,skew,sca,delta),ncol=dd[2])*(pollen-flow[[x]])
    }
    else{ ##local pollination
      ei <- matrix(eiDist(ei.distr,dd[1]*dd[2],seed+2,r,m,lambda,skew,sca),ncol=dd[2])
      no <- 1:length(flow)
      sample <- sample(no[-x],2)
      a = sample[1]
      b = sample[2]
      flow[[x]]+ei*(flow[[a]]-flow[[b]])
    }
  }))
}
