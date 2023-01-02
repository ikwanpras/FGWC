
##penghitungan index
index_fgwc <- function(data,cluster,uij,vi,m,a=exp(1)) {
  result<-list()
  result$PC <- PC1(uij)
  result$CE <- CE1(uij,a)
  result$SC <- SC1(data,cluster,uij,vi,m)
  result$SI <- SI1(data,uij,vi)
  result$XB <- XB1(data,uij,vi,m)
  result$IFV <- IFV1(data,uij,vi,m)
  result$Kwon <- Kwon1(data,uij,vi,m)
  return(result)
}

########################################################
#################VALIDATION MEASUREMENT#################
########################################################

##kelompok yang optimum dinyatakan dengan nilai PC yang maksimum.
PC1 <- function(uij) {
  return(sum(uij^2)/nrow(uij))
}

##kelompok yang optimum dinyatakan dengan nilai indeks CE yang minimum.
CE1 <- function(uij,a=exp(1)) {##
  return(sum(uij*log(uij,a))/(-nrow(uij)))
}

##Partisi yang optimum dinyatakan dengan nilai indeks SC yang minimum.
SC1 <- function(data,cluster,uij,vi,m) {
  d <- matrix(0,nrow(data),nrow(vi))
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i,j] <- sum((vi[j,]-data[i,])^2)
    }
  }
  pt1 <- colSums((uij^m)*d)
  vkvi <- matrix(0,nrow(vi),nrow(vi))
  for (i in 1:nrow(vi)) {
    for (k in 1:nrow(vi)) {
      vkvi[i,k] <- sum((vi[i,]-vi[k,])^2)
    }
    Ni <- length(which(cluster==i))
    vkvi[i,]<-Ni*vkvi[i,]
  }
  pt2 <- colSums(vkvi)
  return(sum(pt1/pt2))
}

##Jumlah kelompok yang optimum dinyatakan dengan nilai indeks S yang minimum.
SI1 <- function(data,uij,vi) {
  d <- matrix(0,nrow(data),nrow(vi))
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i,j] <- sum((vi[j,]-data[i,])^2)
    }
  }
  vkvi <- matrix(0,nrow(vi),nrow(vi))
  for (i in 1:nrow(vi)) {
    for (k in 1:nrow(vi)) {
      vkvi[i,k] <- sum((vi[k,]-vi[i,])^2)
    }
  }
  diag(vkvi) <- Inf
  pt1 <- sum((uij^2)*d)
  pt2 <- nrow(data)*min(vkvi)
  return(sum(pt1/pt2))
}

##Jumlah kelompok yang optimal dinyatakan dengan nilai XB yang minimum.
XB1 <- function(data,uij,vi,m) {
  d <- matrix(0,nrow(data),nrow(vi))
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i,j] <- sum((vi[j,]-data[i,])^2)
      if (d[i,j]==0) {
        d[i,j]==Inf
      }
    }
  }
  pt1 <- sum((uij^m)*d)
  pt2 <- nrow(data)*min(d)
  return(pt1/pt2)
}

##Ketika nilai IFV maksimum maka kualitas cluster semakin baik.
IFV1 <- function(data,uij,vi,m) {
  vkvi <- matrix(0,nrow(vi),nrow(vi))
  for (i in 1:nrow(vi)) {
    for (k in 1:nrow(vi)) {
      vkvi[i,k] <- sum((vi[k,]-vi[i,])^2)
    }
  }
  d <- matrix(0,nrow(data),nrow(vi))
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i,j] <- sum((vi[j,]-data[i,])^2)
    }
  }
  sigmaD <- sum(d)/(nrow(data)*nrow(vi))
  SDmax <- max(vkvi)
  log2u <- colSums(log(uij,2))/nrow(data)
  u2ij <- colSums(uij^2)
  inside <- sum(u2ij*(log(nrow(vi),2)-log2u))
  return(sum(u2ij*((log(nrow(vi),2)-log2u)^2)/nrow(data)*(SDmax/sigmaD)))
}

##Ketika nilai Kwon minimum maka kualitas cluster semakin baik.
Kwon1 <- function(data,uij,vi,m) {
  d <- matrix(0,nrow(data),nrow(vi))
  s <- matrix(0,nrow(vi))
  vivj <- matrix(0,nrow(vi),nrow(vi))
  for (j in 1:nrow(vi)) {
    s[j,] <- (sum(vi[j,]-colMeans(data))^2)
  }
  for (i in 1:nrow(data)) {
    for (j in 1:nrow(vi)) {
      d[i,j] <- sum((vi[j,]-data[i,])^2)
    }
  }
  for (i in 1:nrow(vi)) {
    for (j in 1:nrow(vi)) {
      vivj[i,j] <- sum((vi[i,]-vi[j,])^2)
    }
  }
  diag(vivj) <- Inf
  pt1 <- colSums((uij^m)*d)
  pt2 <- sum(s)/nrow(vi)
  pt3 <- min(vivj)
  return(sum((pt1+pt2)/pt3))
}
