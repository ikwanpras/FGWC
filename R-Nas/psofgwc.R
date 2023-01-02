#' Fuzzy Geographicaly Weighted Clustering with Particle Swarm Optimization
#' @description Fuzzy clustering with addition of spatial configuration of membership matrix with centroid optimization using Particle Swarm Algorithm.
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
#' @param npar number of particles. Can be defined as \code{npar=} in \code{opt_param}. Default is 10.
#' @param vmax maximum velocity to be tolerated. Can be defined as \code{vmax} in \code{opt_param}. Default is 0.7
#' @param pso.same number of consecutive unchange to stop the iteration. Can be defined as \code{same=} in \code{opt_param}.
#' @param c1 Cognitive scaling parameters. Can be defined as \code{c1=} in \code{opt_param}. Default is 0.49
#' @param c2 Social scaling parameters. Can be defined as \code{c2=} in \code{opt_param}. Default is 0.49, 
#' @param w.inert The inertia weight update method between \code{"constant"}, \code{"chaotic"}, \code{"sim.annealing"}, \code{"nat.exponent1"}, \code{"nat.exponent2"} based on Bansal (2011). 
#' Can be defined as \code{type=} in \code{opt_param}. Default is \code{'sim.annealing'}
#' @param wmax Maximum inertia weight. Can be defined as \code{wmax} in \code{opt_param}. Default is 0.9.
#' @param wmin Minimum inertia weight. Can be defined as \code{wmin} in \code{opt_param}. Default is 0.4.
#' @param map Chaotic mapping parameter. Userful when \code{w.inert='chaotic'}. Can be defined as \code{map} in \code{opt_param}. Default is 0.4.

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
#' the Particle Swarm Optimization was developed by \insertCite{Kennedy1995;textual}{naspaclust} in order to get a more optimal solution of a certain complex function.
#' PSO was also improved by \insertCite{Bansal2011;textual}{naspaclust} by modifying the inertia weight. 
#' FGWC using PSO has been implemented previously by some studies \insertCite{fgwcpso,putra2017,Abdussamad}{naspaclust}.

#' @references
#' \insertAllCited{}

#' @seealso \code{\link{fpafgwc}} \code{\link{gsafgwc}}
#' @examples
#' data('census2010')
#' data('census2010dist')
#' data('census2010pop')
#' # First way
#' res1 <- psofgwc(census2010,census2010pop,census2010dist,3,2,'minkowski',4,npar=10)
#' # Second way
#' # initiate parameter
#' param_fgwc <- c(kind='v',ncluster=3,m=2,distance='minkowski',order=3,
#'                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=10)
#' ## tune the PSO parameter
#' pso_param <- c(vi.dist='uniform',npar=15,
#'          vmax=0.8, pso.same=10, c1=0.7, c2=0.6, type='chaotic',
#'                      wmax=0.8,wmin=0.3,map=0.3)
#' ##FGWC with PSO
#' res2 <- fgwc(census2010,census2010pop,census2010dist,'pso',param_fgwc,pso_param)

#' @export

psofgwc <- function(data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1, 
					error=1e-5, max.iter=100,randomN=0,vi.dist="uniform",npar=10,
          vmax=0.7, pso.same=10, c1=0.49, c2=0.49, w.inert='sim.annealing',
                      wmax=0.9,wmin=0.4,map=0.4){
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
  v <- lapply(1:npar, function(x) {set.seed(randomN+x+100); matrix(rnorm(ncluster*d, 0,1), ncluster, d)})
  conv <- c(par.fit[which.min(par.fit)])
  pbest <- par$centroid
  pfit <- par$I
  repeat{
  	minmax <- c(which.min(par.fit)[1],which.max(par.fit)[1])
  	best <- minmax[1]
  	worst <- minmax[2]
    theta <- update_inertia(w.inert,wmax,wmin,map,iter,max.iter)
    v <- lapply(1:npar,function(x) update_v(theta,v[[x]],vmax,c1,c2,par.finalpos,pbest[[x]],par.swarm[[x]],randomN+x))
    par.swarm <- lapply(1:npar, function (x) v[[x]] + par.swarm[[x]])
    par.other <- lapply(1:npar, function(x) uij(data,par.swarm[[x]],m,distance,order))
  	par.other <- lapply(1:npar, function(x) renew_uij(data,par.other[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
  	par.swarm <- lapply(1:npar, function(x) vi(data,par.other[[x]],m))
  	par.fit <- sapply(1:npar, function(x) jfgwcv(data,par.swarm[[x]],m,distance,order))
    pbest.ind <- which(par.fit<pfit)
    if(length(pbest.ind)>0){
      for(i in pbest.ind){
        pbest[[i]] <- par.swarm[[i]]
        pfit[i] <- par.fit[i]
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
    if (iter==max.iter || same==pso.same) break
  }
  finaldata=determine_cluster(datax,par.finalpos.other)
  cluster=finaldata[,ncol(finaldata)]
  print(c(order, ncluster,m, randomN))
  pso <- list("converg"=conv,"f_obj"=jfgwcv(data,par.finalpos,m,distance,order),"membership"=par.finalpos.other,"centroid"=par.finalpos,
              "validation"=index_fgwc(data,cluster,par.finalpos.other,par.finalpos,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"iteration"=iter,"same"=same,"time"=proc.time()-ptm)
  class(pso) <- 'fgwc'
  return(pso)
}

update_v <- function(theta,v0,vmax,c1,c2,gbest,pbest,particle,randomN) {
  n <- nrow(particle)
  d <- ncol(particle)
  set.seed(randomN <- randomN+1)
  e1 <- matrix(runif(n*d),ncol=d,nrow=n)
  set.seed(randomN <- randomN+1)
  e2 <- matrix(runif(n*d),ncol=d,nrow=n)
  v_new <- theta*v0+c1*e1*(gbest-particle)+c2*e2*(pbest-particle)
  for(i in 1:ncol(v_new)) {
    x <- which(v_new[,i]<(-vmax))
    v_new[x,i] <- -vmax
    x <- which(v_new[,i]>(vmax))
    v_new[x,i] <- vmax
  }
  return(v_new)
}

# weight candidate
# constant, chaotic, simulated annealing, natural exp 1 2, exp decreasing
update_inertia <- function(w.inert, wmax, wmin, z, iter, maxiter) {
  if(w.inert=="constant") {
    return(wmax)
  }
  else if(w.inert=="chaotic") {
    return((wmax-wmin)*(1-iter/maxiter)+(wmin*z))
  }
  else if(w.inert=="sim.annealing") {
    return(wmin+((wmax-wmin)*0.95^(iter-1)))
  }
  else if(w.inert=="nat.exponent1") {
    return(wmin+((wmax-wmin)*exp(iter/(maxiter/10))))
  }
  else if(w.inert=="nat.exponent2") {
    return(wmin+((wmax-wmin)*exp((iter/(maxiter/10)^2))))
  }
}