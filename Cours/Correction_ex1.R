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
#TODO : Error in integrate(sinc, -2 * pi, 2 * pi) : non-finite function value
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

#TODO : ce n'est pas l'espérance E, mais la moyenne du min et du max
temp <- b-a
a <- a-temp
b <- b+temp

#abscisses
abs <- a+(b-a)*(1:500)/500

KernelEst <- function(x,h,kern,abs)
{
  n <- length(x)
  m <- length(abs)
  # on crée les points y où on calcule la densité
  #abscisses
  yy <- array(data=c(abs), dim=c(m,n))
  
  #observations
  xx <- t(array(data=c(x), dim=c(n,m))) 
  
  #(x-Xi)/h
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
    #TODO : pourquoi divisé par pi ? Résultat théorique, pas pour les autres ?
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

#TODO : pourquoi 0.4 alors que 0.6 dans la question ?
#Afin de mieux matcher la courbe théorique ?
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
#On fait N simulations différentes, donc dans une boucle
for(i in 1:N){
  #print(i)
  #on prend rnorm(1000) comme x
  
  x=rnorm(1000) 
  #TODO:C'est légèrement plus long de générer rnorm(1000) avant de le mettre dans
  #les kernels que de générer un échantillon différent pour chaque KernelEst ?!
  
  estim <- KernelEst(x,h,"Rect",abs)
  MSE[1,] <- MSE[1,]+(estim$hatf-dnorm(estim$abs))^2/N
  
  estim <- KernelEst(x,h,"Tri",abs)
  MSE[2,] <- MSE[2,]+(estim$hatf-dnorm(estim$abs))^2/N
  
  estim <- KernelEst(x,h,"Gaus",abs)
  MSE[3,] <- MSE[3,]+(estim$hatf-dnorm(estim$abs))^2/N
  
  estim <- KernelEst(x,h,"EP",abs)
  MSE[4,] <- MSE[4,]+(estim$hatf-dnorm(estim$abs))^2/N
  
  #TODO : encore un h différent ici ?
  estim <- KernelEst(x,0.2,"sinc",abs)
  MSE[5,] <- MSE[5,]+(estim$hatf-dnorm(estim$abs))^2/N
}
T2par<-Sys.time() 
Tdiffpar= difftime(T2par, T1par); Tdiffpar

summary(MSE[1,])
summary(MSE[2,])
summary(MSE[3,])
summary(MSE[4,])
summary(MSE[5,])

par(mfcol=c(1,1))

plot(abs, MSE[1,],col="black")
lines(abs, MSE[2,],col="blue")
lines(abs, MSE[3,],col="red")
lines(abs, MSE[4,],col="green")
lines(abs, MSE[5,],col="pink")

# Question supplémentaire - calcul du h par la méthode de l'erreur quadratique 
#moyenne intégrée
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
#TODO : comprendre implémentation de temp. C'est quoi mu2 ?
Denom <- 3*(integrate(temp,-Inf,Inf)$value)^2
# Vérification calcul du cours page 17
resu <- (Num/Denom)^(1/5);resu
#TODO : où est la formule correspondante ?
h_EQMI <- (Num/(length(x)*Denom))^(1/5)*sd(x);h_EQMI

# Noyau quadratiaue
Num <- 8*sqrt(pi)*integrate(Gauss2,-Inf,Inf)$value
Denom <- 3*(integrate(temp3,-Inf,Inf)$value)^2
# Vérification calcul du cours page 17
resu <- (Num/Denom)^(1/5);resu
h_EQMI <- (Num/(length(x)*Denom))^(1/5)*sd(x);h_EQMI

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
#TODO : formule MISE par rapport à MSE ?
MISE <- rowSums(MSE)*(5-(-5))/500
print(MISE)

#########################################
# QUESTION 8
########################################
x <- rnorm(1000)
a <- min(x)
b <- max(x)
temp <- b-a
a <- a-temp
b <- b+temp
abs <- a+(b-a)*(1:500)/500

#pourquoi 91 ? => n variant de 10 à 1000 par pas de 10
resultats <- array(data=c(0), dim=c(91))
for(i in 1:91)
{
  n <- 100+(i-1)*10
  estim <- KernelEst(x[1:n],n^(1/5),"Gaus",abs)
  resultats[i] <- sum((estim$hatf-dnorm(estim$abs))^2*(b-a)/500)
}

plot(seq(100,1000,10),resultats)

###################################
# QUESTION 9
########################################

CV_kern=function(x)
{
  
  n <- length(x)
  a <- min(x)
  b <- max(x)
  temp <- b-a
  a <- a-temp
  b <- b+temp
  abs <- a+(b-a)*(1:500)/500
  
  N <- floor(n/2)*(n<100)+floor(n/4)*(n>=100)
  
  resultatsCV1 <- array(data=c(0),dim=c(N))
  resultatsCV2 <- array(data=c(0),dim=c(N))
  resultatsCV3 <- array(data=c(0),dim=c(N))
  
  for(m in 1:N)
  {
    h <- m*(b-a)/(N*10)
    
    esth <- KernelEst(x,h,"Gaus",abs)
    #matrice des K((X_i-X_j)/h)
    x2 <- array(data=c(x),dim=c(n,n))
    xx <- dnorm((x2-t(x2))/h)
    temp <- rowSums(xx)-diag(xx)
    
    resultatsCV1[m] <- sum((esth$hatf)^2)*((b-a)/500)
    resultatsCV2[m] <- 2/(n*(n-1)*h)*sum(temp)
    resultatsCV3[m] <- resultatsCV1[m]-resultatsCV2[m]
  }
  
  #	plot((1:N)*(b-a)/(10*N),resultatsCV1,type='l',lwd=2,col='darkred',xlab='fenêtre',ylab='CV')
  #	lines((1:N)*(b-a)/(10*N),resultatsCV2, col='blue')
  #	lines((1:N)*(b-a)/(10*N),resultatsCV3, lwd=2,col='red')
  plot((1:N)*(b-a)/(10*N),resultatsCV3,type='l',lwd=2,col='darkred',xlab='fenêtre',ylab='CV')
  
  plot(1:N*(b-a)/(N*10),resultatsCV3)
  return(c(min(resultatsCV3),which.min(resultatsCV3)*(b-a)/(10*N)))
}

CV_kern(rnorm(400))

####################### 

x <- rnorm(1000)

resultats <- array(data=c(0), dim=c(10))

for(i in seq(1:10))
{	
  print(i)
  resultats[i] <- CV_kern(x[1:(100+(i-1)*100)])[2]
}

plot(seq(100,1000,100), resultats, xlim=c(0,1000),type='l')

####################### 

library(MASS)
help(galaxies)
boxplot(galaxies)

gal <- galaxies/1000
c(width.SJ(gal, method = "dpi"), width.SJ(gal))

# On peut aussi utiliser la fonction R density() 
# qui détermine l'estimateur à noyau de la densité
# Help : ?density
a <- min(gal)
b <- max(gal)
abs <- a+(b-a)*(1:500)/500
res <- KernelEst(gal,3.25/2,"Gaus",abs)
res <- KernelEst(gal,0.55,"Gaus",abs)


plot(x = c(0, 40), y = c(0, 0.3), type = "n", bty = "l",
     xlab = "velocity of galaxy (1000km/s)", ylab = "density")
rug(gal)
lines(density(gal, width = 3.25, n = 200), lty = 1)
lines(density(gal, width = 2.56, n = 200), lty = 3)
lines(res$abs,res$hatf,col='red')

CV_kern(gal)