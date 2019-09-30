source("PLASSO.R")


## Case 1: Quantitative traits
n=1000 
p=1000 
k=20
rho=0.4  
mu = rep(0, p)
sigma = matrix(0,p,p)
for(i in 1:p/20){
for(j in 1:p/20){
sigma[i,j]=rho
}}
for(i in 1:nrow(sigma)){sigma[i,i] = 1}
ep=1
MAF=runif(p,min=0.05,max=0.5)

W= mvrnorm(n, mu, sigma, tol = 1e-6)
X=W
for(j in 1:p){X[,j]=SNPf(W[,j],MAF[j])}
X=scale(X)
beta=rep(0, p)
bbb = sample(1:p,k)
beta[bbb]= 0.2*c(rep(1, k))*(2*rbinom(k, 1, 0.5)-1)
y = 1+ X %*% beta + rnorm(n,0,ep^2)

select_plasso=plassoc(X,y,pB=10,SS=0.9)
select_plasso
length(which(select_plasso %in% bbb))/length(select_plasso) ## Precision
length(which(select_plasso %in% bbb))/length(bbb)		## Recall



## Case 2: Dichotomous traits
n=1000 
n0=n1=n/2
p=1000 
k=20
rho=0.4
mu = rep(0, p)
sigma = matrix(0,p,p)
for(i in 1:p/20){
for(j in 1:p/20){
sigma[i,j]=rho
}}
for(i in 1:nrow(sigma)){sigma[i,i] = 1}
R=sigma
svd.R<-svd(R)
R1<-svd.R$u %*% diag(sqrt(svd.R$d))
ep=1
MAF=runif(p,min=0.05,max=0.5)
cutoff=qnorm(MAF)
invlogit <- function(x){ exp(x)/(exp(x)+1)}
beta0=log(0.05/0.95)

beta=rep(0, p)
bbb = sample(1:p,k)
beta[bbb]= 1*c(rep(1, k))*(2*rbinom(k, 1, 0.5)-1)

X<-matrix(0, nrow=n0+n1, ncol=p)
y<-rep(0, n0+n1); y[(n0+1):(n0+n1)]<-1
i<-1
#sampling controls:
while ( i <= n0){
  X0<-rnorm(p, 0, 1) 	#: X0 ~ MVN(0, I)
  X1<-R1 %*% X0   	#: X1 ~ MVN(0, sigma)
  X2<-ifelse(X1<cutoff, 1, 0)
  X0<-rnorm(p, 0, 1) 	#: X0 ~ MVN(0, I)
  X1<-R1 %*% X0   	#: X1 ~ MVN(0, sigma)
  X3<-ifelse(X1<cutoff, 1, 0)
  X4<-X2+ X3
  pr=invlogit(beta0+t(X4) %*% beta)
  Y1<-sample(c(0, 1), 1, prob=c(1-pr, pr))
  if (Y1==0){
    X[i, ]<-X4
    i<-i+1
    }
}

#sampling cases:
while ( i <= n0+n1){
  X0<-rnorm(p, 0, 1) 	#: X0 ~ MVN(0, I)
  X1<-R1 %*% X0   	#: X1 ~ MVN(0, sigma)
  X2<-ifelse(X1<cutoff, 1, 0)
  X0<-rnorm(p, 0, 1) 	#: X0 ~ MVN(0, I)
  X1<-R1 %*% X0   	#: X1 ~ MVN(0, sigma)
  X3<-ifelse(X1<cutoff, 1, 0)
  X4<-X2+ X3
  pr=invlogit(beta0+t(X4) %*% beta)
  Y1<-sample(c(0, 1), 1, prob=c(1-pr, pr))
  if (Y1==1){
    X[i, ]<-X4
    i<-i+1
    }
}
## Note: this sampling process can be slow if the value of beta0+t(X4) %*% beta is either too close to 0 or 1

select_plasso=plassob(X,y,pB=10,SS=0.9)
select_plasso
length(which(select_plasso %in% bbb))/length(select_plasso) ## Precision
length(which(select_plasso %in% bbb))/length(bbb)		## Recall


