## Title: Permutation-assisted Lasso Tuning (plasso)
## Version: 0.1
## Author: Xiang Zhan (xyz5074@psu.edu)
## Date: 09/29/2019

library(glmnet) 
library(MASS)

## plasso for a continuous phenotype:
plassoc=function(X,y,pB,SS){
n=nrow(X)
p=ncol(X)
X=scale(X)
Spi=NULL
for(tt in 1:pB){
X_ko1=X[c(sample(nrow(X))), ]
a=glmnet(cbind(X, X_ko1),y)
lasso=as.matrix(a$beta)
ii=1
while (max(abs(lasso[(p+1):(2*p),ii]))==0 & ii<dim(lasso)[2]){
    ii=ii+1 
}
selected_lasso = which(abs(lasso[,ii-1])>0)
Spi=c(Spi,selected_lasso)
}
freq=tabulate(Spi)/pB
out=which(freq>SS | freq==SS)
return(out)
}


## plasso for a binary phenotype:
plassob=function(X,y,pB,SS){
n=nrow(X)
p=ncol(X)
X=scale(X)
Spi=NULL
for(k in 1:pB){
X_ko1=X[c(sample(nrow(X))), ]
a=glmnet(cbind(X, X_ko1),y,family="binomial")
lasso=as.matrix(a$beta)
ii=1
while (max(abs(lasso[(p+1):(2*p),ii]))==0 & ii<dim(lasso)[2]){
    ii=ii+1 
}
selected_lasso = which(abs(lasso[,ii-1])>0)
Spi=c(Spi,selected_lasso)
}
freq=tabulate(Spi)/pB
out=which(freq>SS | freq==SS)
return(out)
}


## Auxillary function to transfer a normal vector a to a 0/1/2 SNP variable with specifed MAF
SNPf <- function(a,maf){
p2=maf^2
p1=2*maf*(1-maf)
p0=(1-maf)^2
cutoff1= quantile(a,p0)
cutoff2= quantile(a,(1-p2))
b=a  ## Need a new vector to store the results
b[a<cutoff1|a==cutoff1]=0
b[a>cutoff1&a<cutoff2]=1
b[a>cutoff2|a==cutoff2]=2
return(b)
}


