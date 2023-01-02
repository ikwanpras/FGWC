#' Fuzzy Geographicaly Weighted Clustering with Teaching - Learning Based Optimization
#' @description Fuzzy clustering with addition of spatial configuration of membership matrix with centroid optimization using Teaching - Learning Based Algorithm.
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
#' @param nstud number of students. Can be defined as \code{npar=} in \code{opt_param}. Default is 10.
#' @param tlbo.same number of consecutive unchange to stop the iteration. Can be defined as \code{same=} in \code{opt_param}. Default is 10.
#' @param nselection number of teachers based on selected students. Can be defined as \code{nselection=} in \code{opt_param}. Default is equal to \code{nstud}.
#' @param elitism wheter to use elitism algorithm or not. Either \code{TRUE} or \code{FALSE}. Can be defined as \code{elitism=} in \code{opt_param}. Default is \code{FALSE}.
#' @param n.elite Number of elitist students. Can be defined as \code{n.elite=} in \code{opt_param}. Default is 2.

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
#' the Teaching - Learning Based Optimization was developed by \insertCite{Rao2012;textual}{naspaclust} and Developed by \insertCite{Rao2012b;textual}{naspaclust} 
#' by adding the elitism algorithm in order to get a more optimal solution of a certain complex function.

#' @references
#' \insertAllCited{}

#' @seealso \code{\link{fpafgwc}} \code{\link{gsafgwc}}
#' @examples
#' data('census2010')
#' data('census2010dist')
#' data('census2010pop')
#' # First way
#' res1 <- tlbofgwc(census2010,census2010pop,census2010dist,3,2,'minkowski',4,nstud=10)
#' # Second way
#' # initiate parameter
#' param_fgwc <- c(kind='v',ncluster=3,m=2,distance='minkowski',order=3,
#'                alpha=0.5,a=1.2,b=1.2,max.iter=1000,error=1e-6,randomN=10)
#' ## tune the TLBO parameter
#' tlbo_param <- c(vi.dist="uniform",nstud=10, tlbo.same=10,
#'          nselection=10,elitism=FALSE,n.elite=2)
#' ##FGWC with TLBO
#' res2 <- fgwc(census2010,census2010pop,census2010dist,'tlbo',param_fgwc,tlbo_param)

#' @export


tlbofgwc <- function(data, pop=NA, distmat=NA, ncluster=2, m=2, distance='euclidean', order=2, alpha=0.7, a=1, b=1, 
					error=1e-5, max.iter=100,randomN=0,vi.dist="uniform",nstud=10, tlbo.same=10,
          nselection=10,elitism=F,n.elite=2){
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

  stud <- init.swarm(data, mi.mj, distmat, distance, order, vi.dist, ncluster, 
                                m, alpha, a, b, randomN, nstud)
  stud.swarm <- stud$centroid
  stud.other <- stud$membership
  stud.fit <- stud$I

  stud.finalpos <- stud$centroid[[which.min(stud.fit)]]
  stud.finalpos.other <- stud$membership[[which.min(stud.fit)]]
  stud.fit.finalbest <- stud$I[[which.min(stud.fit)]]
  conv <- c(stud.fit[which.min(stud.fit)])
  repeat{
  	minmax <- c(which.min(stud.fit)[1],which.max(stud.fit)[1])
  	best <- minmax[1]
  	worst <- minmax[2]
  	teacher <- stud$centroid[[best]]
    class.ave <- Reduce('+',stud.swarm)/nstud
    studs.new <- teacher.phase(stud.swarm,stud.fit,teacher,class.ave,randomN,data,m,distance,order,mi.mj,distmat,alpha,beta,a,b)
    studs.new <- learner.phase(studs.new$studs,studs.new$fit,randomN+5,data,m,distance,order,mi.mj,distmat,alpha,beta,a,b)

    if(elitism==TRUE){
    	stud.swarm <- elitism(studs.new$studs,studs.new$fit,stud.finalpos,stud.fit.finalbest,n.elite,randomN+6)
    }
    else{
    	stud.swarm <- studs.new$studs
    }
    stud.other <- lapply(1:nstud, function(x) uij(data,stud.swarm[[x]],m,distance,order))
  	stud.other <- lapply(1:nstud, function(x) renew_uij(data,stud.other[[x]]$u,mi.mj,distmat,alpha,beta,a,b))
  	stud.swarm <- lapply(1:nstud, function(x) vi(data,stud.other[[x]],m))
  	stud.fit <- sapply(1:nstud, function(x) jfgwcv(data,stud.swarm[[x]],m,distance,order))
    best <- which(stud.fit==min(stud.fit))[1]
    stud.curbest <- stud.swarm[[best]]
    stud.curbest.other <- stud.other[[best]]
    stud.fit.curbest <- stud.fit[best]
    conv <- c(conv,stud.fit.finalbest)
    iter <- iter+1
    if (abs(conv[iter+1]-conv[iter])<error) same <- same+1
    else same <- 0
    if (stud.fit.curbest<=stud.fit.finalbest) {
      stud.finalpos <- stud.curbest
      stud.finalpos.other <- stud.curbest.other
      stud.fit.finalbest <- stud.fit.curbest
    }
    randomN <- randomN+nstud
    if (iter==max.iter || same==tlbo.same) break
  }
  finaldata=determine_cluster(datax,stud.finalpos.other)
  cluster=finaldata[,ncol(finaldata)]
  print(c(order, ncluster,m, randomN))
  tlbo <- list("converg"=conv,"f_obj"=jfgwcv(data,stud.finalpos,m,distance,order),"membership"=stud.finalpos.other,"centroid"=stud.finalpos,
              "validation"=index_fgwc(data,cluster,stud.finalpos.other,stud.finalpos,m,exp(1)), "cluster"=cluster,
              "finaldata"=finaldata, "call"=match.call(),"iteration"=iter,"same"=same,"time"=proc.time()-ptm)
  class(tlbo) <- 'fgwc'
  return(tlbo)
}

teacher.phase <- function(studs,studs.fit,teacher,average,seed,data,m,distance,order,mi.mj,dist,alpha,beta,a,b){
  set.seed(seed <- seed+2)
  tf <- matrix(1+runif(ncol(teacher)*nrow(teacher)),ncol=ncol(teacher))
  set.seed(seed <- seed+1)
  r <- matrix(runif(ncol(teacher)*nrow(teacher)),ncol=ncol(teacher))
  diff <- r*(teacher-tf*average)
  nstud <- length(studs)
  stud2 <- lapply(1:nstud, function (x) studs[[x]]+diff)
  stud.fit2 <- sapply(1:nstud, function (x) jfgwcv2(data,stud2[[x]],m,distance,order,mi.mj,dist,alpha,beta,a,b))
  whichone <- which(stud.fit2<studs.fit)
  for(i in whichone){
    studs[[i]] <- stud2[[i]]
    studs.fit[i] <- stud.fit2[i]
  }
  return(list(studs=studs,fit=studs.fit))
}

learner.phase <- function(studs,studs.fit,seed,data,m,distance,order,mi.mj,dist,alpha,beta,a,b){
	set.seed(seed <- seed+2)
	sample1 <- sample(1:length(studs),length(studs))
	set.seed(seed <- seed+2)
	sample2 <- sample(1:length(studs),length(studs))
	for(i in 1:length(studs)){
		set.seed(seed <- seed+1)
		r <- matrix(runif(ncol(studs[[1]])*nrow(studs[[1]])),ncol=ncol(studs[[1]]))	
		if(sample1[i]!=sample2[i]){
	  	a <- sample1[i]
	    b <- sample2[i]
	    if(studs.fit[a]<studs.fit[b]) stud2 <- studs[[a]]+r*(studs[[a]]-studs[[b]])
		  else stud2 <- studs[[a]]+r*(studs[[b]]-studs[[a]])
		  stud.fit2 <- jfgwcv2(data,stud2,m,distance,order,mi.mj,dist,alpha,beta,a,b)
			if(stud.fit2<studs.fit[a]){
				studs[[a]] <- stud2
    		studs.fit[a] <- stud.fit2
			}
	  }
	}
	return(list(studs=studs,fit=studs.fit))
}

elitism <- function(swarm,fit,gbest,gfit,n.elite,seed){
	worst <- order(fit,decreasing = T)[1:n.elite]
	for(i in worst){
		swarm[[i]] <- gbest
		fit[i] <- gfit
	}
	dup <- which(duplicated(fit))
	for(i in dup){
		set.seed(seed <- seed+i)
		r <- matrix(rnorm(ncol(swarm[[1]])*nrow(swarm[[1]])),ncol=ncol(swarm[[1]]))
		swarm[[i]] <- swarm[[i]]+r
	}
	return(swarm)
}