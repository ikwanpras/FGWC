#' Fuzzy Geographicaly Weighted Clustering with Gravitational Search Algorithm
#' @description Fuzzy clustering with addition of spatial configuration of membership matrix with centroid optimization using Gravitational Search Algorithm
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
#' @param npar number of particle. Can be defined as \code{npar=} in \code{opt_param}. Default is 10.
#' @param par.no The number of selected best particle. Can be defined as \code{par.no=} in \code{opt_param}. Default is 2
#' @param par.dist The distance between particles. Can be defined as \code{par.dist=} in \code{opt_param}. Default is \code{'euclidean'}, 
#' @param par.order The minkowski order of the \code{par.dist} if \code{par.dist='minkowski'}. Can be defined as \code{par.order=} in \code{opt_param}. Default is 2
#' @param gsa.same number of consecutive unchange to stop the iteration. Can be defined as \code{same=} in \code{opt_param}.
#' @param G initial gravitatioal constant, Can be defined as \code{G} in \code{opt_param}. default is 1.
#' @param vmax maximum velocity to be tolerated. Can be defined as \code{vmax} in \code{opt_param}. Default is 0.7
#' @param new Boolean that represents whether to use the new algorithm by Li and Dong (2017). Can be defined as \code{new} in \code{opt_param}. Default is \code{FALSE}

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

#' @details Fuzzy Geographically Weighted Clustering (FGWC) was developed by Mason and Jacobson (2007) by adding 
#' neighborhood effects and population to configure the membership matrix in Fuzzy C-Means. Furthermore,
#' the Gravitational Search Algorithm was developed by \insertCite{rashedi2009;textual}{naspaclust} and 
#' and the technique is also upgraded by \insertCite{Li2017gsa;textual}{naspaclust} in order to get a more optimal
#' solution of a certain complex function. FGWC using GSA has been implemented previously by \insertCite{fgwcgsa;textual}{naspaclust}.


#' @references
#' \insertAllCited{}

#' @seealso \code{\link{fpafgwc}} \code{\link{gsafgwc}}
#' @examples
#' data('census2010')
#' data('census2010dist')
#' data('census2010pop')
#' # First way
#' res1 <- gsafgwc(census2010,census2010pop,census2010dist,3,2,'euclidean',4,npar=10)
#' # Second way
#' # initiate parameter
#' param_fgwc <- c(kind='v',ncluster=3,m=2,distance='minkowski',order=3,
#'                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=10)
#' ## tune the GSA parameter
#' gsa_param <- c(vi.dist='normal',npar=5,same=15,G=1,vmax=0.7,new=FALSE) 
#' ##FGWC with GSA
#' res2 <- fgwc(census2010,census2010pop,census2010dist,'gsa',param_fgwc,gsa_param)

#' @export


gsafgwc <- function(data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1, 
					error=1e-5, max.iter=100,randomN=0,vi.dist="uniform",npar=10,par.no=2,par.dist='euclidean', par.order=2,
          gsa.same=10, G=1, vmax=0.7, new=F){
  # require(beepr)
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

  par <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster, 
                                m, alpha, a, b, randomN, npar)
  par.swarm <- par$centroid
  par.other <- par$membership
  par.fit <- par$I

  par.finalpos <- par$centroid[[which.min(par.fit)]]
  par.finalpos.other <- par$membership[[which.min(par.fit)]]
  par.fit.finalbest <- par$I[[which.min(par.fit)]]
  v <- lapply(1:npar, function(x) matrix(0, ncluster, d))
  #{set.seed(randomN+x+100); matrix(runif(ncluster*d, 0,1), ncluster, d)}
  pbest <- par$centroid
  pfit <- par$I
  conv <- c(par.fit[which.min(par.fit)])
  repeat{
  	minmax <- c(which.min(par.fit)[1],which.max(par.fit)[1])
  	best <- minmax[1]
  	worst <- minmax[2]
    G <- G*runif(1,0.95,1)
    v <- force_v(par,par.no,G,v,vmax,par.dist,par.order,randomN)
    par.swarm <- lapply(1:npar, function (x) v[[x]] + par.swarm[[x]])
    if(new==TRUE){
      par.swarm <- lapply(1:npar,function(x) new.move(par.swarm[[x]],pbest[[x]],par.finalpos,randomN+x))
    }
    par.other <- lapply(1:npar, function(x) uij(data,par.swarm[[x]],m,distance,order))
  	par.other <- par$membership <- lapply(1:npar, function(x) renew_uij(data,par.other[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
  	par.swarm <- par$centroid <- lapply(1:npar, function(x) vi(data,par.other[[x]],m))
    par.fit <- par$I <- sapply(1:npar, function(x) jfgwcv(data,par.swarm[[x]],m,distance,order))

    if(new==TRUE){ 
      pbest.ind <- which(par.fit<pfit)
      if(length(pbest.ind)>0){
        for(i in pbest.ind){
          pbest[[i]] <- par.swarm[[i]]
          pfit[i] <- par.fit[i]
        }
      }
    }

    best <- which(par.fit==min(par.fit))[1]
    par.curbest <- par.swarm[[best]]
    par.curbest.other <- par.other[[best]]
    par.fit.curbest <- par.fit[best]
    conv <- c(conv,par.fit.finalbest)
    iter <- iter+1
    if (abs(conv[iter+1]-conv[iter])<error) same <- same+1
    else same <- 0
    if (par.fit.curbest<=par.fit.finalbest) {
      par.finalpos <- par.curbest
      par.finalpos.other <- par.curbest.other
      par.fit.finalbest <- par.fit.curbest
    }
    randomN <- randomN+npar
    if (iter==max.iter || same==gsa.same) break
  }
  finaldata=determine_cluster(datax,par.finalpos.other)
  cluster=finaldata[,ncol(finaldata)]
  print(c(order, ncluster,m, randomN))
  gsa <- list("converg"=conv,"f_obj"=jfgwcv(data,par.finalpos,m,distance,order),"membership"=par.finalpos.other,"centroid"=par.finalpos,
              "validation"=index_fgwc(data,cluster,par.finalpos.other,par.finalpos,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"iteration"=iter,"same"=same,"time"=proc.time()-ptm)
  class(gsa) <- 'fgwc'
  return(gsa)
}

force_v <- function(par,no,G,v,vmax,par.dist,par.order,randomN){
  dd <- dim(par$centroid[[1]])
  intel.par <- intel.ffly(par,no)
  mass <- (par$I-max(par$I))/(min(par$I)-max(par$I))
  Mass <- mass/sum(mass)
  Mass.intel <- sort(Mass,decreasing=T)[1:no]
  v1 <- v
  for(i in 1:length(par$centroid)){
    Fij <- lapply(1:no,c)
    for(j in 1:no){
      r <- diag(cdist(par$centroid[[i]],intel.par$centroid[[j]],par.dist,par.order))
      set.seed(randomN <- randomN+1)
      eps <- runif(length(r),0,1e-6)
      set.seed(randomN <- randomN+1)
      rand <- matrix(runif(dd[1]*dd[2]),ncol=dd[2])
      Fij[[j]] <- rand*G*Mass[i]*Mass.intel[j]*(intel.par$centroid[[j]]-par$centroid[[i]])/(r+eps)
    }
    Fi <- Reduce("+",Fij)
    a <- Fi/Mass[i]
    set.seed(randomN <- randomN+1)
    # rand <- matrix(runif(dd[1]*dd[2]),ncol=dd[2])
    
    v1[[i]] <- rand*v1[[i]]+a
    # for(i in 1:nrow(v1[[i]])){
    #   v1[[i]][v1[[i]]< -vmax,] <- -vmax
    #   v1[[i]][v1[[i]]> vmax,] <- vmax
    # }
  }
  return(v)
}

new.move <- function(par,pbest,gbest,randomN){ ##Li dan Dong, 2017 GSA new technique
  dd <- dim(par)
  mu <- (par+pbest+gbest)/3
  sigma <- sqrt(((par-mu)^2+(pbest-mu)^2+(gbest-mu)^2)/3)
  set.seed(randomN+100)
  c1 <- matrix(runif(dd[1]*dd[2], 0,1), ncol=dd[2])
  set.seed(randomN+101)
  c2 <- matrix(runif(dd[1]*dd[2], 0,1), ncol=dd[2])
  z <- sqrt(-2*log(c1))*cos(2*pi*c2)
  return(mu+sigma*z)
}