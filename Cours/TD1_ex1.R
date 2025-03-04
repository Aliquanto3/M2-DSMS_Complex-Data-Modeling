# Exercice 1
# 
# Dans un premier temps, on va générer un échantillon X1, . . . , Xn de densité connue (une loi N(0, 1))
# que l’on utilisera pour reconstruire la densité et la comparer au résultat théorique. Dans un second temps,
# on considèrera l’estimation de densité sur un jeu de vraies données.
# 
# Partie 1

# 1. Générer un échantillon x de n = 1000 v.a. i.i.d X1, . . . , Xn de loi N(0, 1).

set.seed(1234)
X=rnorm(1000)
hist(X,proba=T)
lines(density(X),col="red")

# 2. On considèrera dans la suite 5 noyaux : la densité de la loi N(0, 1) et les 
# 4 autres noyaux suivants :
# Représenter graphiquement ces noyaux.

#Densité de la loi N(0,1)
#dnorm(x)

Tri=function(x){
  (1-abs(x))*(abs(x)<=1)
}

Rect=function(x){
  (1/2*abs(x))*(abs(x)<=1)
}

EP=function(x){
  return((3/4*(1-x**2))*(abs(x)<=1))

}

sinc=function(x){
  sin(x)/x
}

#####################################
#Vérification condition (C2)
integrate(dnorm,-Inf,Inf)
integrate(Tri,-Inf,Inf)
integrate(Rect,-Inf,Inf)
integrate(EP,-Inf,Inf)
integrate(sinc,-Inf,Inf)

Trib=function(x){
  x*Tri(x)
}
integrate(Trib,-Inf,Inf)

Tric=function(x){
  Tri(x)**2
}
tau2=integrate(Tric,-Inf,Inf);tau2

Trid=function(x){
  abs(x)*Tri(x)
}
integrate(Trid,-Inf,Inf)
#####################################

#Représenter graphiquement ces noyaux
par(mfrow=c(2,3))

curve(dnorm(x), from=-5, to=5, , xlab="x", ylab="y",main="Loi N(0,1)")
curve(Tri(x), from=-5, to=5, , xlab="x", ylab="y",main="Tri")
curve(Rect(x), from=-5, to=5, , xlab="x", ylab="y",main="Rect")
curve(EP(x), from=-5, to=5, , xlab="x", ylab="y",main="EP")
curve(sinc(x), from=-50, to=50, , xlab="x", ylab="y",main="sinc")

dev.off()

z=seq(-5,5,.01)
plot(z,dnorm(z),type="l",ylim=c(-0.5,1),xlim=c(-5,5),xlab="x",ylab="y",col="red",main="Noyaux")
lines(z,Tri(z),col="yellow")
lines(z,Rect(z),col="green")
lines(z,EP(z),col="blue")
lines(z,sinc(z),col="orange")
legend(4, 0.95, legend=c("N(0,1)", "Tri","Rect","EP","sinc"),
       col=c("red", "yellow","green","blue","orange"),lty=1)

dev.off()

# 3. On définit une grille 500 pas sur laquelle sera calculée la densité. On 
# considère l’estimation sur l’intervalle [a, b] où a = min(Xi) − E, 
# b = max(Xi) + E où E = max(Xi) − min(Xi) est l’étendue des valeurs de 
# l’échantillon.
# Coder une fonction R KernelEst qui :
#   — prend en arguments : 1) x l’échantillon dont il faut reconstruire la 
# densité f, 2) h la fenêtre, 3) une chaîne de caractère, Tri, Rect, EP, Gaus, 
# sinc qui indiquera quel noyau utiliser pour l’estimation de la densité, 
# 4) le vecteur abs des abscisses des 500 points où l’estimateur de la
# densité sera calculé.
# — calcule l’estimateur à noyaux de la densité ˆf, dont on rappelle la définition
# — renvoie les valeurs ˆf(x) pour les points x de la grille définie ci-dessus, 
# ainsi que les valeurs de ces points x.

Grid500=expand.grid(X=seq(-5,5,length=500));Grid500

X=rnorm(1000)

a=min(X);a
b=max(X);b
tmp=b-a;tmp
a=a-tmp;a
b=b+tmp;b

abs=a+(b-a)*(1:500)/500;abs

KernelEst=function(x,h,Kern,abs,plotGraph){
  
  n=length(x)
  m=length(abs)
  
  xx=t(array(data=c(x),dim=c(n,m)))
  yy=array(data=c(abs),dim=c(m,n))
  
  z=(xx-yy)/h
  
  Kres=NULL
  
  if(Kern=="Tri"){
    Kres=(1-abs(z))*(abs(z)<=1)
    k2=function(u){
      return(((1-abs(u))*(abs(u)<=1))**2)
    }
    #faire pareil pour les autres cas
    
  }else if(Kern=="Rect"){
    Kres=(1/2*abs(z))*(abs(z)<=1)
    k2=function(u){
      return(((1/2*abs(u))*(abs(u)<=1))**2)
    }
    
  }else if(Kern=="EP"){
    Kres=(3/4*(1-z**2))*(abs(z)<=1)
    k2=function(u){
      return(((3/4*(1-u**2))*(abs(u)<=1))**2)
    }
    
  }else if(Kern=="Gaus"){
    Kres=dnorm(z,0,sd=1/2.2)
    k2=function(u){
      return((dnorm(u,0,sd=1/2.2))**2)
    }
    
  }else if(Kern=="sinc"){
    Kres=sin(z)/z
    k2=function(u){
      return((sin(u)/u)**2)
    }
  }
  
  f_hat=1/(n*h)*as.vector(rowSums(Kres))
  
  #Calcul de l'intervalle de confiance
  IC_inf=abs-qnorm(1-alpha/2)*sqrt(abs*integrate(k2,-Inf,Inf))/(n*h)
  IC_sup=abs+qnorm(1-alpha/2)*sqrt(abs*integrate(k2,-Inf,Inf))/(n*h)
  
  if(plotGraph==TRUE){
    plot(abs,f_hat,ylab="densité",type="l",col="blue",lwd=2,
         main=paste("Noyau :",Kern,"- h :",h,sep = " "))
  }
  
  return(list(hatf=f_hat,abs=abs))
}

# 4. Utiliser la fonction précédente pour estimer la densité de x en utilisant 
# les noyaux EP et sinc avec h = 0.6. Commenter. Comment améliorer l’estimation 
# avec le noyau sinc ?

h = 0.6
X=rnorm(1000)

dev.off()

par(mfrow=c(2,1))

head(KernelEst(X,h,"EP",abs))
lines(abs,EP(abs),col="red",lwd=2)
lines(abs,dnorm(abs),col="green",lwd=2)

head(KernelEst(X,h,"sinc",abs))
lines(abs,sinc(abs),col="red",lwd=2)
lines(abs,dnorm(abs),col="green",lwd=2)

dev.off()


# 5. Avec h = 0.8, comparer l’estimation de la densité de x en utilisant les 
# noyaux EP, Rect, Tri et le noyau gaussien.

h=0.8

dev.off()

par(mfrow=c(2,3))

head(KernelEst(X,h,"Tri",abs))
lines(abs,Tri(abs),col="red",lwd=2)
lines(abs,dnorm(abs),col="green",lwd=2)
head(KernelEst(X,h,"Rect",abs))
lines(abs,Rect(abs),col="red",lwd=2)
lines(abs,dnorm(abs),col="green",lwd=2)
head(KernelEst(X,h,"EP",abs))
lines(abs,EP(abs),col="red",lwd=2)
lines(abs,dnorm(abs),col="green",lwd=2)
head(KernelEst(X,h,"Gaus",abs))
lines(abs,dnorm(abs),col="red",lwd=2)
lines(abs,dnorm(abs),col="green",lwd=2)

dev.off()

par(mfrow=c(1,1))


#6. Calculer un estimateur Monte-Carlo du MSE en x, basé sur N simulations 
#indépendantes de x :
#
#Pour chacun des noyaux, donner des statistiques résumant la distribution des 
#MSE(x).

library(MLmetrics)

h=0.2

resK=KernelEst(X,h,"Tri",abs)$hatf
MSE(resK,dnorm(X))
# print(paste("La MSE pour Tri est :",
#             round(1/length(X)*sum((resK-dnorm(X))**2),digit=4)))
summary((1/length(X)*(resK-dnorm(X))**2))

resK=KernelEst(X,h,"Rect",abs)$hatf
MSE(resK,dnorm(X))
# print(paste("La MSE pour Rect est :",
#             round(1/length(X)*sum((resK-dnorm(X))**2),digit=4)))
summary((1/length(X)*(resK-dnorm(X))**2))

resK=KernelEst(X,h,"EP",abs)$hatf
MSE(resK,dnorm(X))
# print(paste("La MSE pour EP est :",
#             round(1/length(X)*sum((resK-dnorm(X))**2),digit=4)))
summary((1/length(X)*(resK-dnorm(X))**2))

resK=KernelEst(X,h,"Gaus",abs)$hatf
MSE(resK,dnorm(X))
# print(paste("La MSE pour Gauss est :",
#             round(1/length(X)*sum((resK-dnorm(X))**2),digit=4)))
summary((1/length(X)*(resK-dnorm(X))**2))

resK=KernelEst(X,h,"sinc",abs)$hatf
MSE(resK,dnorm(X))
# print(paste("La MSE pour sinc est :",
#             round(1/length(X)*sum((resK-dnorm(X))**2),digit=4)))
summary((1/length(X)*(resK-dnorm(X))**2))

N=100
h=0.6
MSE=array(data=c(0),dim=c(5,500))

for (i in 1:N){
  estim=KernelEst(rnorm(1000),h,"Rect",abs)
  MSE[1,]=MSE[1,]+(estim$hatf-dnorm(estim$abs))**2
  estim=KernelEst(rnorm(1000),h,"Tri",abs)
  MSE[2,]=MSE[2,]+(estim$hatf-dnorm(estim$abs))**2
  estim=KernelEst(rnorm(1000),h,"Gaus",abs)
  MSE[3,]=MSE[3,]+(estim$hatf-dnorm(estim$abs))**2
  estim=KernelEst(rnorm(1000),h,"EP",abs)
  MSE[4,]=MSE[4,]+(estim$hatf-dnorm(estim$abs))**2
  estim=KernelEst(rnorm(1000),h,"sinc",abs)
  MSE[5,]=MSE[5,]+(estim$hatf-dnorm(estim$abs))**2
}

summary(MSE[1,])
summary(MSE[2,])
summary(MSE[3,])
summary(MSE[4,])
summary(MSE[5,])

# 7. En déduire le MISE pour chacun des noyaux.

MISE=rowSums((5-(-5))/500*(MSE));MISE

#Question supplémentaire : vérification calcul du cours p17

EP2=function(x){
  return(EP(x)^2)
}

temp=function(x){
  return((x^2)*EP(x))
}

Num=8*sqrt(pi)*integrate(EP2,-Inf,Inf)$value
Denom=3*(integrate(temp,-Inf,Inf)$value)^2
h_EQMI=(Num/(length(x)*Denom))^(1/5)*sd(x);h_EQMI

Quad=function(x){
  15/16*(1-x**2)**2*(abs(x)<=1)
}
Quad2=function(x){Quad(x)**2}
temp2=function(x){
  return((x^2)*Quad(x))
}
Num=8*sqrt(pi)*integrate(Quad2,-Inf,Inf)$value
Denom=3*(integrate(temp2,-Inf,Inf)$value)^2
h_EQMI=(Num/(length(x)*Denom))^(1/5)*sd(x);h_EQMI

# 8. Tracer l’évolution du MISE pour h = n**(−1/5) en fonction de n. 
# Pour cela, on simulera un échantillon de 1000 v.a. i.i.d. N(0, 1) et on en 
# considèrera les sous échantillons X1, . . . , Xn pour n variant de 100 à 1000 
# par pas de 10.

X=rnorm(1000)

a=min(X);a
b=max(X);b
tmp=b-a;tmp
a=a-tmp;a
b=b+tmp;b

abs=a+(b-a)*(1:500)/500;abs

nn=seq(100,1000,10)
h = nn**(-1/5)

Xlist=list()

MSE=array(data=c(0),dim=c(5,500))
MISE=c()

for (i in 1:length(nn)){
  Xn=sample(x=X,size=nn[i],replace=FALSE)
  estim=KernelEst(Xn,h[i],"Rect",abs)
  MSE[1,]=MSE[1,]+(estim$hatf-dnorm(estim$abs))**2
  MISE[i]=sum((5-(-5))/500*(MSE))
}

plot(nn,MISE)

#####Corrigé 8.
MISE=array(data=c(0),dim=c(91))
for (i in 1:91){
  n=100+(i-1)*10
  estim=KernelEst(X[1:n],n^(-1/5),"Gaus",abs)
  MISE[i]=sum((estim$hatf-dnorm(estim$abs))^2*(b-a)/500)
} 
taille=seq(100,1000,10)
plot(taille,MISE)
plot(taille,MISE,ylim=c(0,0.01),type='l')


# 9. On souhaite maintenant minimiser le MISE en h, pour un choix de noyau fixé. 
# Rappeler quelle est la fonction J(h) de h qu’il suffit de minimiser pour 
# trouver arg minh MISE(h) ? Donner un estimateur sans biais CV(h) de J(h). 
# Programmer une fonction qui renvoie la valeur de CV (h) et la valeur de h qui 
# minimise cette fonction lorsque l’on fait varier h sur une grille de N pas 
# espacés de (b−a)/10N où a et b sont définis à la question 3) et où N = n/2 si 
# n ≤ 100 et n/4 si n > 100.

#Calcul du h par méthode de l'erreur quadratique intégrée 

J=function(h){
  
}

KernelIMoins1=function(x,h,Kern,abs,i){
  
  n=length(x)
  m=length(abs)
  
  xx=t(array(data=c(x),dim=c(n,m)))
  yy=array(data=c(abs),dim=c(m,n))
  
  z=(xx-yy)/h
  
  Kres=NULL
  
  if(Kern=="Tri"){
    Kres=(1-abs(z))*(abs(z)<=1)
    
  }else if(Kern=="Rect"){
    Kres=(1/2*abs(z))*(abs(z)<=1)
    
  }else if(Kern=="EP"){
    Kres=(3/4*(1-z**2))*(abs(z)<=1)
    
  }else if(Kern=="Gaus"){
    Kres=dnorm(z)
    
  }else if(Kern=="sinc"){
    Kres=sin(z)/z
  }
  
  f_hat=1/(n-1)*as.vector(rowSums(Kres[-i]))
  
  return(list(hatf=f_hat,abs=abs))
}

CV=function(x,h,Kern,abs){
  
  SumFIMoins1=0
  for (i in 1:n){
    SumFIMoins1=SumFIMoins1+sum(KernelIMoins1(x,h,Kern,i))
  }
  integrate(KernelEst,lower = -Inf, upper = Inf) - 
    2/n*SumFIMoins1
  
}

a=min(X);a
b=max(X);b
tmp=b-a;tmp
a=a-tmp;a
b=b+tmp;b

MinCV=function(n){

  N=floor(n/2)*(n<100)+floor(n/4)*(n>=100)
  GridCV=expand.grid(X=seq(a,b,by=(b-a)/(10*N)))
  
  h=optim(fn=j,par=c(GridCV))
  
  c(h,CV(h))
}

###Correction 9.

CV_Kern=function(x,plotGraph){
  a=min(X)
  b=max(X)
  tmp=b-a
  a=a-tmp
  b=b+tmp
  
  abs=a+(b-a)*(1:500)/500;abs
  
  N=floor(n/2)*(n<100)+floor(n/4)*(n>=100)
  
  resultatsCV1=array(data=c(0),dim=c(N))
  resultatsCV2=array(data=c(0),dim=c(N))
  resultatsCV3=array(data=c(0),dim=c(N))
  
  for (m in 1:N){
    h=m*(b-a)/(N*10)
    esth=KernelEst(x,h,"Gaus",abs,plotGraph)
    x2=array(data=c(x),dim=c(n,n))
    xx=dnorm(x2-t(x2)/h)
    temp=rowSums(xx)-diag(xx)
    
    resultatsCV1[m]=sum((esth$hatf)^2)*((b-a)/500)
    resultatsCV2[m]=2/(n*(n-1)*h)*sum(temp)
    resultatsCV3[m]=resultatsCV1[m]-resultatsCV2[m]
  }
  
  if(plotGraph==TRUE){
    plot((1:N)*(b-a)/(10*N),resultatsCV3,type="l",lwd=2,col="darkred",
         xlab = "fenêtre",ylab="CV")
    
    plot((1:N)*(b-a)/(10*N),resultatsCV3)
  }
  return(c(min(resultatsCV3),which.min(resultatsCV3)*(b-a)/(10*N)))
}
CV_Kern(X,T)
h_CV=CV_Kern(X,plotGraph=F)[2];h_CV

CV_Kern(rnorm(1000),F)

### Partie 2
# On cherche maintenant à étudier les données correspondant au mouvement de 82 
# galaxies.
# 1. Charger les données en tapant les commandes

library(MASS)
data(galaxies)

# 2. Donner quelques statistiques descriptives.
summary(galaxies)
str(galaxies)
hist(galaxies)
gal=galaxies/1000

# 3. Obtenir les estimations de la densité en utilisant les fonctions 
# programmées de la partie 1. Quelle est la fenêtre optimale ?
CV_Kern(gal,F)

head(KernelEst(gal,h,"EP",abs,F))

MISE=array(data=c(0),dim=c(91))
for (i in 1:length(gal)){
  n=100+(i-1)*10
  estim=KernelEst(gal,n^(-1/5),"Gaus",abs,F)
  MISE[i]=sum((estim$hatf-dnorm(estim$abs))^2*(b-a)/500)
} 
taille=seq(100,1000,10)
plot(taille,MISE)
plot(taille,MISE,ylim=c(0,0.01),type='l')


a=min(gal)
b=max(gal)
abs=a+(b-a)*(1:500)/500
res=KernelEst(gal,h=3.25/2,'Gaus',abs,T)
res=KernelEst(gal,0.55,'Gaus',abs,T)

plot(x=c(0,40),y=c(0,0.3),type="n",bty="l",
     xlab = "velocity of galaxy (1000km/s)",ylab = "density")
rug(gal)
lines(density(gal,width = 3.25, n=200),lty=1)
lines(density(gal,width = 2.56, n=200),lty=3)
lines(res$abs,res$hatf,col='red')

CV_Kern(gal,T)
