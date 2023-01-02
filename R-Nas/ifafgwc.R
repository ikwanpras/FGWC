#' Fuzzy Geographicaly Weighted Clustering with (Intelligent) Firefly Algorithm
#' @description Fuzzy clustering with addition of spatial configuration of membership matrix with centroid optimization using (Intelligent) Firefly Algorithm.
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
#' @param vi.dist a string of centroid population distribution between \code{'uniform'} (default) and \code{'normal'}. Can be defined as \code{vi.dist=} in \code{opt_param}.
#' @param ei.distr distribution of random walk parameter. Can be defined as \code{ei.distr} in \code{opt_param}.
#' @param fa.same number of consecutive unchange to stop the iteration. Can be defined as \code{same=} in \code{opt_param}.
#' @param nfly number of fireflies. Can be defined as \code{npar=} in \code{opt_param}. Default is 10.
#' @param ffly.no The number of selected best fireflies for intelligent firefly algorithm. Can be defined as \code{par.no=} in \code{opt_param}. Default is 2
#' @param ffly.dist The distance between fireflies. Can be defined as \code{par.dist=} in \code{opt_param}. Default is \code{'euclidean'}, 
#' @param ffly.order The minkowski order of the \code{par.dist} if \code{par.dist='minkowski'}. Can be defined as \code{par.order=} in \code{opt_param}. Default is 2
#' @param gamma distance scaling factor. Can be defined as \code{gamma} in \code{opt_param}. Default is 1.
#' @param ffly.beta Attractiveness constant. Can be defined as \code{beta} in \code{opt_param}. Default is 1.
#' @param ffly.alpha Randomisation constant. Can be defined as \code{alpha=} in \code{opt_param}.
#' @param r.chaotic weight in logistic chaotic between [0,4]. Can be used when \code{ei.distr='logchaotic'}. Can be defined as \code{chaos} in \code{opt_param}. Default is 4.
#' @param m.chaotic mapping parameter in kent chaotic between [0,1]. Can be used when \code{ei.distr='kentchaotic'}. Can be defined as \code{map} in \code{opt_param}. Default is 0.7.
#' @param ind.levy Levy distribution index for random walk. Can be used when \code{ei.distr='levy'}. Can be defined as \code{ind} in \code{opt_param}. Default is 1.
#' @param skew.levy Levy distribution skewness for random walk. Can be used when \code{ei.distr='levy'}. Can be defined as \code{skew} in \code{opt_param}. Default is 0.
#' @param scale.levy Levy distribution scale for random walk. Can be used when \code{ei.distr='levy'}. Can be defined as \code{sca} in \code{opt_param}. Default is 1.
#' @param ffly.alpha.type An integer. The type of \code{ffly.alpha} update. Can be selected from 1 to 5. Can be defined as \code{update_type} in \code{opt_param}. Default is 4.

#' @return an object of class \code{'fgwc'}.\cr
#' An \code{'fgwc'} object contains as follows:
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
#' the Firefly Algorithm was developed by \insertCite{Yang2009;textual}{naspaclust} and the technique is also upgraded by 
#' \insertCite{intfa;textual}{naspaclust} by adding the intelligent phase (choosing the best firefly based on the intensity) in order to get a more optimal
#' solution of a certain complex function. FGWC using IFA has been implemented previously by \insertCite{Nasution2020;textual}{naspaclust}.

#' @references
#' \insertAllCited{}

#' @seealso \code{\link{fpafgwc}} \code{\link{gsafgwc}}
#' @examples
#' data('census2010')
#' data('census2010dist')
#' data('census2010pop')
#' # First way
#' res1 <- ifafgwc(census2010,census2010pop,census2010dist,3,2,'minkowski',4,nfly=10)
#' # Second way
#' # initiate parameter
#' param_fgwc <- c(kind='v',ncluster=3,m=2,distance='minkowski',order=3,
#'                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=10)
#' ## tune the IFA parameter
#' ifa_param <- c(vi.dist='uniform', ei.distr='logchaotic',
#'						fa.same=10, npar=15, par.no=3, par.dist='minkowski', 
#'            par.order=4, gamma=1, beta=1.5,
#'            alpha=1, chaos=4,update_type=4) 
#' ##FGWC with IFA
#' res2 <- fgwc(census2010,census2010pop,census2010dist,'ifa',param_fgwc,ifa_param)

#' @export


########################################################
#############INTELLIGENT FIREFLY ALGORITHM##############
########################################################
ifafgwc <- function (data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1,
						error=1e-5, max.iter=100,randomN=0,vi.dist="uniform", ei.distr="normal",
						fa.same=10, nfly=10, ffly.no=2, ffly.dist='euclidean', ffly.order=2, gamma=1, ffly.beta=1,
            ffly.alpha=1, r.chaotic=4,m.chaotic=0.7,ind.levy=1,skew.levy=0,scale.levy=1,ffly.alpha.type=4) {
  #require(beepr)
  randomnn <- randomN
  ptm<-proc.time()
  n <- nrow(data)
  d <- ncol(data)
  iter=0
  gen=1
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
  ffly.finalbest <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster,
                                m, alpha, a, b, randomN, 1)
  inten.finalbest <- ffly.finalbest$I
  conv <- c(inten.finalbest)
  ffly.finalpos <- ffly.finalbest$centroid
  ffly.finalpos.other <- ffly.finalbest$membership
  ffly.new <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster,
                        m, alpha, a, b, randomN, nfly)
  ffly.swarm <- ffly.new$centroid
  ffly.other <- ffly.new$membership
  inten <- ffly.new$I
  repeat {
    set.seed(randomN)
    ffly.alpha <- update_alpha(ffly.alpha,gen,max.iter,ffly.alpha.type)
    set.seed(randomN)
    ffly.swarm <- ffly.new$centroid <- moving(ffly.new,ffly.no,ffly.beta,gamma,ffly.alpha,ffly.dist,ffly.order,ei.distr,
      r.chaotic,m.chaotic,ind.levy,skew.levy,scale.levy,data,m,distance,order,mi.mj,distmat,alpha,beta,a,b,randomN)
    ffly.other <- ffly.new$membership <- lapply(1:nfly, function(x) uij(data,ffly.swarm[[x]],m,distance,order))
    ffly.other <- ffly.new$membership <- lapply(1:nfly, function(x) renew_uij(data,ffly.new$membership[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
    ffly.swarm <- ffly.new$centroid <- lapply(1:nfly, function(x) vi(data,ffly.new$membership[[x]],m))
    inten <- ffly.new$I <- sapply(1:nfly, function(x) jfgwcv(data,ffly.new$centroid[[x]],m,distance,order))
    best <- which(inten==min(inten))[1]
    ffly.curbest <- ffly.swarm[[best]]
    ffly.curbest.other <- ffly.other[[best]]
    inten.curbest <- inten[best]
    conv <- c(conv,inten.finalbest)
    if (abs(conv[gen+1]-conv[gen])<error) {
      same <- same+1
    }
    else {
      same <- 0
    }
    if (inten.curbest<=inten.finalbest) {
      ffly.finalpos <- ffly.curbest
      ffly.finalpos.other <- ffly.curbest.other
      inten.finalbest <- inten.curbest
    }
    gen <- gen+1
    randomN <- randomN+nfly
    if (gen==max.iter || same==fa.same) break
  }##end repeat
  # print(class(ffly.finalpos.other))
  if (any(class(ffly.finalpos.other)=="list")) {
  	new_uij <- ffly.finalpos.other[[1]]
  	vi <- ffly.finalpos[[1]]
  }
  else {
  	new_uij <- ffly.finalpos.other
  	vi <- ffly.finalpos
  }
  finaldata=determine_cluster(datax,new_uij)
  cluster=finaldata[,ncol(finaldata)]
  ifa <- list("converg"=conv,"f_obj"=jfgwcv(data,vi,m,distance,order),"membership"=new_uij,"centroid"=vi,
              "validation"=index_fgwc(data,cluster,new_uij,vi,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"iteration"=gen,"same"=same,"time"=proc.time()-ptm)
  print(c(order, ncluster,m, randomN))
  #result <- list(ifa=ifa,fgwc=fgwc)
  class(ifa) <- 'fgwc'
  return (ifa)
}

init.swarm <- function(data, pop, distmat, distance, order, vi.dist, ncluster,
                        m, alpha, a, b, randomN, nfly) {
  inten <- rep(0,nfly)
  beta <- 1-alpha
  start.uij <- lapply(1:nfly, function(x) gen_uij(data,ncluster,nrow(data),randomN+x))
	start.uij <- lapply(1:nfly, function (x) renew_uij(data,start.uij[[x]],pop,distmat,alpha,beta,a,b))
	start.vi <- lapply(1:nfly, function (x) vi(data,start.uij[[x]],m))
  for(i in 1:nfly) {
    inten[i] <- jfgwcv(data,start.vi[[i]],m,distance,order)
  }
  result <- list("membership"=start.uij,"centroid"=start.vi,"I"=inten)
  return(result)
}

intel.ffly <- function(ffly.list,no) {
  best <- order(ffly.list$I,decreasing = F)[1:no]
  intel.uij <- lapply(best, function(x) ffly.list$membership[[x]])
  intel.vi <- lapply(best, function(x) ffly.list$centroid[[x]])
  inten <- ffly.list$I[best]
  result <- list("membership"=intel.uij,"centroid"=intel.vi,"I"=inten)
  return(result)
}

swarm_dist <- function (swarm1,swarm2,distance,order) {
  # jarak<-rep(0,nrow(swarm1))
  # for (i in 1:nrow(swarm1)) {
  #     diff <- abs(swarm1[i,]-swarm2[i,])^order
  #     jarak[i] <- sum(diff)^(1/order)
  # }
  return(diag(dist(swarm1,swarm2,distance,order)))
}

moving <- function(ffly.all,no,ff.beta,gamma,ff.alpha,ffly.dist,ffly.order,ei.distr,r.chaotic,m.chaotic,ind.levy,skew.levy,sca.levy,
                  data,m,distance,order,mi.mj,dist,alpha,beta,a,b,randomN){##menggerakkan firefly
  times <- 0
  intel.ffly <- intel.ffly(ffly.all,no)
  dd <- dim(ffly.all$centroid[[1]])
  ffly <- ffly.all$centroid
  fit <- ffly.all$I
  for(i in 1:length(intel.ffly$centroid)){
    for(j in 1:length(ffly)){
      r <- diag(cdist(ffly[[j]],intel.ffly$centroid[[i]],ffly.dist,ffly.order))
      ei <- matrix(eiDist(ei.distr,dd[1]*dd[2],randomN+i+j,r.chaotic,m.chaotic,ind.levy,skew.levy,sca.levy),ncol=dd[2])
      if (fit[j] > intel.ffly$I[i]){
        ffly[[j]]+beta*exp(-gamma*r^2)*(intel.ffly$centroid[[i]]-ffly[[j]])+(ff.alpha*ei)
      }
      else{
        times <- times+1
        if(times==no){
          ffly[[j]]+(ff.alpha*ei)
        }
      }
      fit[j] <- jfgwcv2(data,ffly[[j]],m,distance,order,mi.mj,dist,alpha,beta,a,b)
    }
  }
  return(ffly)
}
