eiDist <- function(distr='normal',n,randomN=40,r=4,m=0.7,ind=1,skew=0,sca=1) {
  set.seed(randomN)
  if(distr=='uniform') {
    return(runif(n,-1,1))
  }
  else if(distr=='normal') {
    return(rnorm(n,0,1))
  }
  else if (distr=="levy") {
    # require(stabledist)
    return(rstable(n,ind,skew,sca))
  }
  else if (distr=="logchaotic") {
    return(logchaotic(n,r,randomN))
  }
  else if (distr=="kentchaotic") {
    return(kentchaotic(n,m,randomN))
  }
  else if (distr=='sinechaotic'){
    return(sinechaotic(n,m,randomN))
  }
  else if (distr=='dyadchaotic'){
    return(dyadchaotic(n,randomN))
  }
  else if (distr=='chebychaotic'){
    return(chebychaotic(n,randomN))
  }
  else if (distr=='circhaotic'){
    return(circhaotic(n,randomN))
  }

}

logchaotic <- function(n,r=4,seed=1) {
  set.seed(seed)
  x0 <- runif(1)
  x <- c()
  for(i in 1:n) {
    x0 <- r*x0*(1-x0)
    x <- c(x,x0)
  }
  return(x)
}

kentchaotic <- function(n,m=0.7,seed) {
  set.seed(seed)
  x0 <- runif(1)
  x <- c(x0)
  for(i in 2:n) {
    if(x0>0 && x0<=m) {
      x0 <- x0/m
    }
    else {
      x0 <- (1-x0)/(1-m)
    }
    x <- c(x,x0)
  }
  return(x)
}

sinechaotic <- function(n,m,seed){
  set.seed(seed)
  x0 <- runif(1)
  x <- c()
  for(i in 1:n) {
    x0 <- m/4*sin(x0)
    x <- c(x,x0)
  }
  return(x)
}

dyadchaotic <- function(n,seed){
  set.seed(seed)
  x0 <- runif(1)
  x <- c()
  for(i in 1:n) {
    x0 <- (2*x0)%%1
    x <- c(x,x0)
  }
  return(x)
}

chebychaotic <- function(n,seed){
  set.seed(seed)
  x0 <- runif(1)
  x <- c()
  for(i in 1:n) {
    x0 <- cos((i)*acos(x0))
    x <- c(x,x0)
  }
  return(x)
}

circhaotic <- function(n,seed){
  set.seed(seed)
  x0 <- runif(1)
  x <- c()
  for(i in 1:n) {
    x0 <- x0+0.2-(0.5-2*pi)*sin(2*pi*x0)
    x <- c(x,x0)
  }
  return(x)
}

update_alpha <- function(alpha, iter, maxiter, type) {
  if(type==1) {
    return(1e-5+(alpha-(1e-5))*exp(-iter))
  }
  else if(type==2) {
    return(alpha*runif(1,0.95,0.99)^iter)
  }
  else if(type==3) {
    delta <- 1-(10^(-4)/9^(1/maxiter))
    return(1-delta*alpha)
  }
  else if (type==4) {
    return((1.11*10^(-4))^(5/maxiter)*alpha)
  }
  else if (type==5){
  	return(alpha*(1-(iter/maxiter)))
  }
}