########################################################
#####################CLASSICAL FGWC#####################
########################################################

#' Classical Fuzzy Geographicaly Weighted Clustering
#' @description Fuzzy clustering with addition of spatial configuration of membership matrix
#' @param data an object of data with d>1. Can be \code{matrix} or \code{data.frame}. If your data is univariate, bind it with \code{1} to get a 2 columns.
#' @param pop an n*1 vector contains population.
#' @param distmat an n*n distance matrix between regions.
#' @param kind use \code{'u'} if you want to use membership approach and \code{'v'} for centroid approach.
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
#' @param uij membership matrix initialisation.
#' @param vi centroid matrix initialisation.

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
#' neighborhood effects and population to configure the membership matrix in Fuzzy C-Means. There are two kinds of options in doing classical FGWC.
#' The first is using \code{"u"} \insertCite{Runkler2006}{naspaclust} (default) for membership optimization 
#' and \code{"v"} \insertCite{fgwc}{naspaclust} for centroid optimisation. 

#' @references
#' \insertAllCited{}

#' @seealso \code{\link{abcfgwc}} \code{\link{fpafgwc}} \code{\link{gsafgwc}} \code{\link{hhofgwc}} \code{\link{ifafgwc}} \code{\link{psofgwc}} \code{\link{tlbofgwc}}
#' @examples
#' data('census2010')
#' data('census2010dist')
#' data('census2010pop')
#' res1 <- fgwcuv(census2010,census2010pop,census2010dist,'u',3,2,'euclidean',4)
#'

#' @export
fgwcuv <- function(data, pop, distmat, kind=NA,ncluster=2, m=2, distance='euclidean', order=2,
                      alpha=0.7, a=1, b=1, max.iter=500, error=1e-5,
                     randomN=0, uij=NA, vi=NA) {
  ##populasi berupa matriks n x 1
  ##jarak berupa matriks jarak antar daerah
  ##alpha + beta = 1
  ##m = fuzzifier
  ptm <- proc.time()
  stopifnot(kind=="v" || kind=="u" || any(is.na(kind))==TRUE)
  n <- nrow(data)
  d <- ncol(data)
  beta = 1-alpha
  iter = 0
  conv <- c(0)
  if (is.matrix(data)==FALSE) {
    data <- as.matrix(data)
  }
  ##jika alfa =1, akan menjadi fuzzy c-means,
  ##populasi dan matriks jarak dianggap matriks 1
  if (alpha==1) {
    pop <- rep(1,n)
    distmat <- matrix(1,n,n)
  }
  ##membaca matriks populasi dan menjadikannya matriks dengan 1 kolom
  pop <- matrix(pop,ncol=1)
  mi.mj <- pop%*%t(pop)
  ##membaca pendekatan FGWC jika dikosongkan atau u,
  ##pendekatannya akan menjadi matriks keanggotaan, inisialisasi awal matriks keanggotaan
  if (any(is.na(kind))==TRUE || kind=="u"){ ##fgwc biasa = fgwc u
    if (any(is.na(uij))==TRUE) {
      set.seed(randomN)
      uij <- matrix(runif(n*ncluster,0,1),ncol=ncluster)
      new_uij  <- uij/rowSums(uij)
    }
    else {
      new_uij <- uij
    }
    old_uij <- new_uij+1
    while (max(abs(new_uij-old_uij))>error && iter<max.iter) {
      old_uij <- new_uij
      vi <- vi (data,old_uij,m) ##mengubah matriks keanggotaan menjadi centroid
      uij <- uij(data,vi,m,distance,order) ##centroid yang didapat diubah menjadi matriks keanggotaan
      new_uij <- renew_uij(data,uij$u,mi.mj,distmat,alpha,beta,a,b) ##memodifikasi matriks keanggotaan dengan memanfaatkan matriks jarak dan populasi
      iter = iter+1
      conv <- c(conv,sum(new_uij^m*uij$d))
    }
  }
  else { ## fgwc.v
    ##stopifnot(any(is.na(vi))==F)
    if (is.na(vi)) { ##jika centroid tidak diinisialisasi, akan dilakukan generate dengan menggunakan distribusi uniform
      vi <- gen_vi(data,ncluster,"uniform",randomN)
    }
    v_new <- vi
    uij <- uij(data,v_new,m,'euclidean',2) ##menghitung matriks keanggotaan dengan menggunakan jarak euclidean untuk awal
    new_uij <- uij$u ##matriks keanggotaan yang disimpan dimasukkan ke yang baru
    v_old <- vi+1
    while (abs(sum(v_new-v_old))>error && iter<max.iter) {
      v_old <- v_new
      uij <- uij(data,v_old,m,distance,order)##memperbarui matriks keanggotaan dengan menggunakan jarak minkowski sesuai orde
      new_uij <- renew_uij(data,uij$u,mi.mj,distmat,alpha,beta,a,b) ##memodifikasi matriks keanggotaan dengan memanfaatkan matriks jarak dan populasi
      v_new <- vi(data,new_uij,m) ##menghitung centroid yang baru
      vi <- v_new ##centroid yang lama diperbarui
      conv <- c(conv,jfgwcv(data,vi,m,distance,order)) ##menghitung fungsi objektif dengan menggunakan pendekatan v
      iter = iter+1
    }
  }
  fgwc_obj <- sum(new_uij^m*uij$d) ##fungsi objektif akhir
  finaldata <- determine_cluster(data,new_uij)
  cluster <- finaldata[,ncol(finaldata)]
  result <- list("converg"=conv[-1],"f_obj"=fgwc_obj, "membership"=new_uij,"centroid"=vi, "validation"=index_fgwc(data,cluster,new_uij,vi,m,exp(1)) ,"iteration" = iter,
                   "cluster"=cluster, "finaldata"=finaldata, "call"=match.call(), "time" = proc.time()-ptm)
  print(c(order, ncluster,m, randomN))
  class(result) <- 'fgwc'
  return (result)
}

##vi dari nilai keanggotaan
vi <- function(data,uij,m) {
  return (t(uij^m)%*%data/colSums(uij^m))
}

##uij dari centroid
uij <- function(data,vi,m,distance,order=2) {
  u <- matrix(0,nrow(data),nrow(vi))
  d <- cdist(data,vi,distance,order)^2
  x <- (d)^(1/(m - 1))
  u <-  (1/x)/rowSums(1/x)
  res <- list("d"=d, "u"=u)
  return(res)
}

##menentukan cluster dari data
determine_cluster <- function(data,uij) {
	clust = apply(uij,1,which.max)
	return(cbind.data.frame(data,cluster=clust))
}

##memodifikasi matriks keanggotaan dengan memanfaatkan matriks jarak dan populasi
renew_uij <- function(data,old_uij,mi.mj,dist,alpha,beta,a,b) {
  diag(dist) <- Inf
  wij <- mi.mj^b/dist^a
  ##new uij
  wijmuj <- wij %*% old_uij
  A <- rowSums(wijmuj)
  new_uij <-  alpha*old_uij + (beta/A)*wijmuj
  return(new_uij)
}

##generate matrik keanggotaan
gen_uij <- function(data,ncluster,n,randomN) {
  set.seed(randomN)
  uij <- matrix(runif(ncluster*n,0,1),n,ncluster)
  return(uij/rowSums(uij))
}

##generate pusat cluster
gen_vi <- function(data,ncluster,gendist,randomN) {##generate centroid
  p <- ncol(data)
  piclass <- matrix(0,ncluster,p)
  for (i in 1:p) {
    set.seed(randomN)
    if (gendist=="normal"){
      piclass[,i] <- rnorm(ncluster,mean(data[,i]),sd(data[,i]))
    }
    else if (gendist=="uniform"){
      piclass[,i] <- runif(ncluster,min(data[,i]),max(data[,i]))
    }
  }
  return(piclass)
}

##fungsi objektif
jfgwcu <- function(data,uij,m,distance,order) { ##fungsi objektif fgwc-u
  vi <- (t(uij^m)%*%data/colSums(uij^m))
  d <- cdist(data,vi,distance,order)^2
  return(sum((uij^m)*d))
}

##fungsi objektif
jfgwcu2 <- function(data,uij,m,distance,order,mi.mj,dist,alpha,beta,a,b) { ##fungsi objektif fgwc-u
  u <- renew_uij(data,u$u,mi.mj,dist,alpha,beta,a,b)
  vi <- (t(uij^m)%*%data/colSums(uij^m))
  d <- cdist(data,vi,distance,order)^2
  return(sum((uij^m)*d))
}

jfgwcv  <- function(data,vi,m,distance,order) { ##fungsi objektif fgwc-v
  u <- matrix(0,nrow(data),nrow(vi))
  d <- cdist(data,vi,distance,order)^2
  x <- (d)^(1/(m - 1))
  for (i in 1:nrow(d)) {
    for (j in 1:ncol(d)) {
      temp <- (d[i,j]/d[i,])^(1/(m-1))
      u[i,j] <- 1/sum(temp)
    }
  }
  return(sum((u^m)*d))
}

jfgwcv2  <- function(data,vi,m,distance,order,mi.mj,dist,alpha,beta,a,b) { ##fungsi objektif fgwc-v
  u <- uij(data,vi,m,distance,order)
  u <- renew_uij(data,u$u,mi.mj,dist,alpha,beta,a,b)
  vi <- vi(data,u,m)
  d <- cdist(data,vi,distance,order)^2
  for (i in 1:nrow(d)) {
    for (j in 1:ncol(d)) {
      temp <- (d[i,j]/d[i,])^(1/(m-1))
      u[i,j] <- 1/sum(temp)
    }
  }
  return(sum((u^m)*d))
}