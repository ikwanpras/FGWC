#' Fuzzy Geographicaly Weighted Clustering with Harris-Hawk Optimization
#' @description Fuzzy clustering with addition of spatial configuration of membership matrix with centroid optimization using Harris-Hawk Algorithm.
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
#' @param nhh number of harris-hawk eagles. Can be defined as \code{npar=} in \code{opt_param}. Default is 10.
#' @param hh.alg String between default is \code{'heidari'} and default is \code{'bairathi'}. The algorithm for HHO, Can be defined as \code{algo} in \code{opt_param}. default is \code{'heidari'}.
#' @param A a 3 vectors which represents initial energy and cut-off for exploitation and exploration. In \code{opt_param}, they can be defined as \code{'a1'} for initial energy, \code{'a2'} for exploitation cut-off and \code{'a3'} for exploration cut-off respectively. default is \code{c("a1"=2,"a2"=1,"a3"=0.5)}.
#' @param p a real number between 0 and 1. The eagle's movement probability
#' @param hh.same number of consecutive unchange to stop the iteration. Can be defined as \code{same=} in \code{opt_param}.
#' @param levy.beta The skewness of levy flight. Can be defined as \code{beta} in \code{opt_param}. Default is 1.5
#' @param update.type An integer. The type of energy \code{A[1]} update. Can be selected from 1 to 5. Can be defined as \code{update.type} in \code{opt_param}. Default is 5.

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
#' the Harris-Hawk Optimization was developed by \insertCite{Bairathi2018;textual}{naspaclust} and the technique is also upgraded by 
#' \insertCite{Heidari2019;textual}{naspaclust} by adding progressive rapid dives in order to get a more optimal solution of a certain complex function.

#' @references
#' \insertAllCited{}

#' @seealso \code{\link{fpafgwc}} \code{\link{gsafgwc}}
#' @examples
#' data('census2010')
#' data('census2010dist')
#' data('census2010pop')
#' # First way
#' res1 <- hhofgwc(census2010,census2010pop,census2010dist,3,2,'euclidean',4,nhh=10)
#' # Second way
#' # initiate parameter
#' param_fgwc <- c(kind='v',ncluster=3,m=2,distance='minkowski',order=3,
#'                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=10)
#' ## tune the HHO parameter
#' hho_param <- c(vi.dist='normal',npar=5,same=15,algo='bairathi',a1=3,a2=1,a3=0.4)
#' ##FGWC with HHO
#' res2 <- fgwc(census2010,census2010pop,census2010dist,'hho',param_fgwc,hho_param)

#' @export


hhofgwc <- function(data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1, 
					error=1e-5, max.iter=100,randomN=0,vi.dist="uniform",nhh=10,hh.alg='heidari',
          A=c(2,1,0.5),p=0.5,hh.same=10,levy.beta=1.5,update.type=5){
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

  hh <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster, 
                                m, alpha, a, b, randomN, nhh)
  hh.swarm <- hh$centroid
  hh.other <- hh$membership
  hh.fit <- hh$I
  
  hh.finalpos <- hh$centroid[[which.min(hh.fit)]]
  hh.finalpos.other <- hh$membership[[which.min(hh.fit)]]
  hh.fit.finalbest <- hh$I[[which.min(hh.fit)]]

  conv <- c(hh.fit[which.min(hh.fit)])
  repeat{
  	minmax <- c(which.min(hh.fit)[1],which.max(hh.fit)[1])
  	best <- minmax[1]
  	worst <- minmax[2]
  	rabbit <- hh.finalpos
    sample <- sample(1:nhh,nhh,T)
    set.seed(randomN+iter)
    if(tolower(hh.alg) == 'bairathi'){
      A[1] <- update_alpha(A[1],iter,max.iter,update.type)
      hh.swarm <- lapply(1:nhh, function(x) hh.attack.bairathi(hh.swarm[[x]],hh.swarm,hh.finalpos,A,p,sample[x],
        randomN+x,best))
    }
  	else {
      E0 <- A[1]*(2*runif(nhh)-1)
      E <- update_alpha(E0,iter,max.iter,update.type)
      hh.swarm <- lapply(1:nhh, function(x) hh.attack.heidari(hh.swarm[[x]],hh.swarm,rabbit,E[x],A,p,
      	sample[x],levy.beta,randomN+x,best,worst,hh.fit[x],data,m,distance,order,mi.mj,distmat,alpha,beta,a,b))
    }
    set.seed(randomN+iter)
  	hh.other <- lapply(1:nhh, function(x) uij(data,hh.swarm[[x]],m,distance,order))
  	hh.other <- lapply(1:nhh, function(x) renew_uij(data,hh.other[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
  	hh.swarm <- lapply(1:nhh, function(x) vi(data,hh.other[[x]],m))
  	hh.fit <- sapply(1:nhh, function(x) jfgwcv(data,hh.swarm[[x]],m,distance,order))
    best <- which(hh.fit==min(hh.fit))[1]
    hh.curbest <- hh.swarm[[best]]
    hh.curbest.other <- hh.other[[best]]
    hh.fit.curbest <- hh.fit[best]
    conv <- c(conv,hh.fit.finalbest)
    iter <- iter+1
    if (abs(conv[iter+1]-conv[iter])<error) same <- same+1
    else same <- 0
    if (hh.fit.curbest<=hh.fit.finalbest) {
      hh.finalpos <- hh.curbest
      hh.finalpos.other <- hh.curbest.other
      hh.fit.finalbest <- hh.fit.curbest
    }
    randomN <- randomN+nhh
    if (iter==max.iter || same==hh.same) break
  }
  finaldata=determine_cluster(datax,hh.finalpos.other)
  cluster=finaldata[,ncol(finaldata)]
  print(c(order, ncluster,m, randomN))
  hho <- list("converg"=conv,"f_obj"=jfgwcv(data,hh.finalpos,m,distance,order),"membership"=hh.finalpos.other,"centroid"=hh.finalpos,
              "validation"=index_fgwc(data,cluster,hh.finalpos.other,hh.finalpos,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"iteration"=iter,"same"=same,"time"=proc.time()-ptm)
  class(hho) <- 'fgwc'
  return(hho)
}

hh.attack.bairathi <- function(hawk,hawks,gbest,A,p,rand,seed,best){
  dd <- dim(hawk)
  set.seed(seed<-seed+10)
  rs <- runif(1)
  set.seed(seed<-seed+10)
  c1 <- A[1]*(2*runif(1)-1)
  set.seed(seed)
  c2 <- matrix(2*(1-runif(dd[1]*dd[2])),ncol=dd[2])
  if(rs<p && c1>=A[2]) return(hawks[[rand]]-c1*(c2*hawks[[rand]]-hawk)) ##exploration
  else if(c1>=A[3]) return(hawks[[best]]-c1*(c2*hawks[[best]]-hawk)) ##local exploitation
  else return((gbest-hawk)-c1*(c2*gbest-hawk)) ##global exploitation
}

hh.attack.heidari <- function(hawk,hawks,rabbit,E,A,p,rand,levy.beta,seed,best,worst,fithawk,data,m,
	distance,order,mi.mj,dist,alpha,beta,a,b){
	set.seed(seed<-seed+10)
	rs <- runif(3)
	dd <- dim(hawk)
	hawks.m <- Reduce('+',hawks)/length(hawks)
	if(E>=A[2]){ ##exploration phase
		if(rs[1]>=p){ ##q
			return(hawks[[rand]]-rs[2]*(hawks[[rand]]-2*rs[3]*hawk))
		}
		else{
			return(rabbit-hawks.m-rs[2]*(hawks[[worst]]+rs[3]*(hawks[[best]]-hawks[[worst]])))
		}
	}
	else{ ##exploitation phase
		set.seed(seed<-seed+10)
		J <- matrix(2*(1-runif(dd[1]*dd[2])),ncol=dd[2])
		if(rs[1]>=p){ 
			if(E>=A[3]) return((rabbit-hawk)-E*(J*rabbit-hawk)) ##soft besiege
			else return(rabbit-E*(rabbit-hawk)) ##hard besiege
		}
		else{ 
			S <- matrix(rnorm(dd[1]*dd[2]),ncol=dd[2])
			LF <- matrix(rlevy(dd[1]*dd[2],levy.beta,seed),ncol=dd[2])
			if(E>=A[3]){ ##soft besiege with progressive rapid dives
				y <- rabbit-E*(J*rabbit-hawk)
			}
			else{ ##hard besiege with progressive rapid dives
				y <- rabbit-E*(J*rabbit-hawks.m)
			}
			z <- y+S*LF
			fity <- jfgwcv2(data,y,m,distance,order,mi.mj,dist,alpha,beta,a,b)
			fitz <- jfgwcv2(data,z,m,distance,order,mi.mj,dist,alpha,beta,a,b)
			if(fity<fithawk) return(y)
			else if(fitz<fithawk) return(z)
			else return(hawk)
		}
	}	
}

rlevy <- function(n,beta,seed){
  set.seed(seed+100)
  u <- runif(n)
  set.seed(seed+101)
  v <- runif(n)
  sigma1 <- gamma(1+beta)*sin(pi*beta/2)
  sigma2 <- gamma((1+beta)/2)*beta*2^((beta-1)/2)
  return(u*((sigma1/sigma2)/v)^(1/beta))
}