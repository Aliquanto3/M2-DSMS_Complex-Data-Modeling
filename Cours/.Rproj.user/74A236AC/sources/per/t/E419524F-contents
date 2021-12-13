# QUESTION 1
x <- rnorm(1000)
hist(x,prob=TRUE,ylim=c(0,0.5))
lines(density(x),col="red")
curve(dnorm,add=T,col="blue")

EP <- function(x)
{
  return(3*(1-x^2)*(abs(x)<=1)/4)
}


verif <- function(x)
{
  return(x*3*(1-x^2)*(abs(x)<=1)/4)
}
sinc <- function(x)
{
  return(sin(x)/(pi*x))
}
integrate(sinc,-2*pi,2*pi)

#########################################
# QUESTION 2
#########################################
z=seq(-5,5,.01)
plot(z,3*(1-z^2)*(abs(z)<=1)/4,type='l',ylim=c(-0.5,1),xlim=c(-5,5),xlab="Noyau")
lines(z,(abs(z)<=1)/2,col='red')
lines(z,(1-abs(z))*(abs(z)<=1),col="blue")
lines(z,dnorm(z,0,sd=1/2.2),col='orange')
lines(z,(1/pi)*sin(z)/z,col='violet')

#########################################
# QUESTION 3
#########################################
a <- min(x)
b <- max(x)
temp <- b-a
a <- a-temp
b <- b+temp

abs <- a+(b-a)*(1:500)/500

KernelEst <- function(x,h,kern,abs)
{
  n <- length(x)
  m <- length(abs)
  # on crée les points y où on calcule la densité
  yy <- array(data=c(abs), dim=c(m,n))
  
  xx <- t(array(data=c(x), dim=c(n,m))) 
  
  z <- (xx-yy)/h
  
  if (kern=='EP')
  {
    hatf <- (1/(n*h))*as.vector(rowSums(3*(1-z^2)*(abs(z)<=1)/4))
  }
  if (kern=='Rect')
  {
    hatf <- (1/(n*h))*as.vector(rowSums((abs(z)<=1)/2))
  }
  if (kern=='Tri')
  {
    hatf <- (1/(n*h))*as.vector(rowSums((1-abs(z))*(abs(z)<=1)))
  }
  if (kern=='Gaus')
  {
    hatf <- (1/(n*h))*as.vector(rowSums(dnorm(z,0,sd=1/2.2)))
  }
  if (kern=='sinc')
  {
    hatf <- (1/(pi*n*h))*as.vector(rowSums(sin(z)/z))
  }
  plot(abs,hatf,ylab="densité",type='l',col='blue',lwd=2,ylim=c(0,0.5))	
  return(list(hatf=hatf,abs=abs))
}

#########################################
# QUESTION 4
#########################################
h <- 0.6
x <- rnorm(1000)

fhat <- KernelEst(x,h,"EP",abs)
lines(fhat$abs,dnorm(fhat$abs),lty=1,col="red",lwd=2)

h <- 0.4
fhat <- KernelEst(x,h,"sinc",abs)
lines(fhat$abs,dnorm(fhat$abs),lty=1,col="red",lwd=2)

#########################################
# QUESTION 5
#########################################
h <- 0.8

par(mfcol=c(2,2))

f1 <- KernelEst(x,h,"Rect",abs)$hatf
lines(abs,dnorm(abs),lty=1,col="red",lwd=2)

f2 <- KernelEst(x,h,"Gaus",abs)$hatf
lines(abs,dnorm(abs),lty=1,col="red",lwd=2)

f3 <- KernelEst(x,h,"EP",abs)$hatf
lines(abs,dnorm(abs),lty=1,col="red",lwd=2)

f4 <- KernelEst(x,h,"Tri",abs)$hatf
lines(abs,dnorm(abs),lty=1,col="red",lwd=2)


#########################################
# QUESTION 6
#########################################
MSE <- array(data=c(0),dim=c(5,500))
par(mfrow=c(3,2))
N <- 100

h <- 0.6

abs <- seq(-5,5-10/500,10/500)
length(abs)

T1par<-Sys.time() 
for(i in 1:N)
{
  #print(i)
  estim <- KernelEst(rnorm(1000),h,"Rect",abs)
  MSE[1,] <- MSE[1,]+(estim$hatf-dnorm(estim$abs))^2/N	
  estim <- KernelEst(rnorm(1000),h,"Tri",abs)
  MSE[2,] <- MSE[2,]+(estim$hatf-dnorm(estim$abs))^2/N
  estim <- KernelEst(rnorm(1000),h,"Gaus",abs)
  MSE[3,] <- MSE[3,]+(estim$hatf-dnorm(estim$abs))^2/N
  estim <- KernelEst(rnorm(1000),h,"EP",abs)
  MSE[4,] <- MSE[4,]+(estim$hatf-dnorm(estim$abs))^2/N
  estim <- KernelEst(rnorm(1000),0.2,"sinc",abs)
  MSE[5,] <- MSE[5,]+(estim$hatf-dnorm(estim$abs))^2/N
}
T2par<-Sys.time() 
Tdiffpar= difftime(T2par, T1par) 
summary(MSE[1,])
summary(MSE[2,])
summary(MSE[3,])
summary(MSE[4,])
summary(MSE[5,])

plot(abs, MSE[1,],col="black")
lines(abs, MSE[2,],col="blue")
lines(abs, MSE[3,],col="red")
lines(abs, MSE[4,],col="green")
lines(abs, MSE[5,],col="pink")

# Question supplémentaire - calcul du h par la méthode de l'erreur quadratique moyenne intégrée
EP <- function(x)
{
  return(3*(1-x^2)*(abs(x)<=1)/4)
}

EP2 <- function(x)
{
  return((EP(x))^2)
}

temp <- function(x)
{
  return(x^2*EP(x))
}

Quad <- function(x)
{
  return((15/16)*(1-x^2)^2*(abs(x)<=1))
}

Quad2 <- function(x)
{
  return(Quad(x)^2)
}

temp2 <- function(x)
{
  return(x^2*Quad(x))
}


Gauss <- function(x)
{
  return(dnorm(x))
}

Gauss2 <- function(x)
{
  return(Gauss(x)^2)
}

temp3 <- function(x)
{
  return(x^2*Gauss(x))
}

# Noyau Epanechnikov
Num <- 8*sqrt(pi)*integrate(EP2,-Inf,Inf)$value
Denom <- 3*(integrate(temp,-Inf,Inf)$value)^2
# Vérification calcul du cours page 17
resu <- (Num/Denom)^(1/5)
h_EQMI <- (Num/(length(x)*Denom))^(1/5)*sd(x)

# Noyau quadratiaue
Num <- 8*sqrt(pi)*integrate(Gauss2,-Inf,Inf)$value
Denom <- 3*(integrate(temp3,-Inf,Inf)$value)^2
# Vérification calcul du cours page 17
resu <- (Num/Denom)^(1/5)
h_EQMI <- (Num/(length(x)*Denom))^(1/5)*sd(x)

# Noyau Gaussien
Num <- 8*sqrt(pi)*integrate(Gauss2,-Inf,Inf)$value
Denom <- 3*(integrate(temp3,-Inf,Inf)$value)^2
# Vérification calcul du cours page 17
resu <- (Num/Denom)^(1/5)
resu
h_EQMI <- (Num/(length(x)*Denom))^(1/5)*sd(x)
h_EQMI


fhat <- KernelEst(x,h_EQMI,"EP",abs)$hatf

#########################################
# QUESTION 7
########################################
MISE <- rowSums(MSE)*(5-(-5))/500
print(MISE)