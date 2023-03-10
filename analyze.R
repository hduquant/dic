


#########################
## Model 1
###############################
library("parallel")

#start_time <- Sys.time()
setwd("/Users/handu/Desktop/multiple_imputation/DIC/real data")
###################################################
# read data
###################################################
data <- read.table("employee.dat")
id<-data[,2]
x<-data[,6]
m<-data[,5]
y<-data[,7]
library(mvtnorm)

###################################################
# read posterior summaries
# the structure of the output file is identical to printed tables on output
# row numbers correspond to parameter labels on blimp output
###################################################

summaries <- read.table("summaries.model1.dat")
names(summaries) <- c("mean","mdn","sd","lcl","ucl","psr")

# posterior summaries

# m model 
model.1 <- summaries[c(3,4,1,2),]
row.names(model.1) <- c("m0","mx1","mv2i","mv1")

# y model
model.2 <- summaries[c(11:13,9:10),]
row.names(model.2) <- c("y0","y1","y2","yv2i","yv1")


###################################################

ll.summary <- read.csv("ll.model1.csv")

###################################################
# read posterior means of random effects
###################################################
random.summary <- read.table("avgimp.model1.dat")

x.im<- random.summary[,6] 
m.im<- random.summary[,5] 
y.im<- random.summary[,7] 
  
# m model 
mv2i <- random.summary[,10] 
var(mv2i)

# y model
yv <- random.summary[,11] 
names(yv) <-c("yv2i")
var(yv)

############## cluster means 
mu.m<- random.summary[,13] 
mu.x<- random.summary[,12] 

#########
###### marginal
y.m<- function(beta=model.2[c(1:3),1],
               varu=matrix(c(model.2[4,1],0,0,0,0,0,0,0,0),nrow=3),
               x=x.im,m=m.im,y=y.im,id=id,npercluster=6,numclusterL2 = 105,
               vare=model.2[5,1]) {
  sigma=diag(npercluster)* vare
  L=lapply(1:numclusterL2,
           function(j) dmvnorm(y[id==j], mean=cbind(1,(m[id==j]-mu.m[id==j]),(x[id==j]-mu.x[id==j]))%*%matrix(as.matrix(beta),nrow=3), 
  ( sigma+cbind(1,(m[id==j]-mu.m[id==j]),(x[id==j]-mu.x[id==j]))%*%varu%*%t(cbind(1,(m[id==j]-mu.m[id==j]),(x[id==j]-mu.x[id==j])))), log=TRUE))
  return(L)
}
m.m<- function(beta=model.1[c(1,2),1],
               varu=matrix(c(model.1[3,1],0,0,0),nrow=2),
               m=m.im,x=x.im,id=id,npercluster=6,numclusterL2 = 105,
               vare=model.1[4,1]) {
  sigma=diag(npercluster)* vare
  L=lapply(1:numclusterL2,
           function(j) dmvnorm(m[id==j], mean=cbind(1,(x[id==j]-mu.x[id==j]))%*%matrix(as.matrix(beta),nrow=2), 
                               ( sigma+cbind(1,(x[id==j]-mu.x[id==j]))%*%varu%*%t(cbind(1,(x[id==j]-mu.x[id==j])))), log=TRUE))
  return(L)
}
################ from LL file
#dic1

dic1.my.m <-(sum(ll.summary[ll.summary[,1]==1,4])+sum(ll.summary[ll.summary[,1]==2,4]))*(-4)+
  2*sum(unlist(y.m(beta=model.2[c(1:3),1],
                   varu=matrix(c(model.2[4,1],0,0,0,0,0,0,0,0),nrow=3),
                   x=x.im,m=m.im,y=y.im,id=id,npercluster=6,numclusterL2 = 105,
                   vare=model.2[5,1])) + 
          unlist( m.m(beta=model.1[c(1,2),1],
                      varu=matrix(c(model.1[3,1],0,0,0),nrow=2),
                      m=m.im,x=x.im,id=id,npercluster=6,numclusterL2 = 105,
                      vare=model.1[4,1]))   )

#   5612.406

#dic2
dic2.my.m <-(sum(ll.summary[ll.summary[,1]==1,4])+sum(ll.summary[ll.summary[,1]==2,4]))*(-4)+
  2*sum(ll.summary[ll.summary[,1]==0,2])
# 5518.996

#waic1
waic1.my.m<- -2*sum(ll.summary[ll.summary[,1]==0,2])+2*sum(ll.summary[ll.summary[,1]==0,6])
#   5576.719

#########
###### conditional

y.c<- function(beta=model.2[c(1:3),1],
               uj=cbind(yv,0,0),x=x.im,m=m.im,y=y.im,id=id,npercluster=6,numclusterL2 = 105,
               vare=model.2[5,1]) {
  sigma=diag(npercluster)* vare
  L=lapply(1:numclusterL2,
           function(j) dmvnorm(y[id==j], mean=cbind(1,(m[id==j]-mu.m[id==j]),(x[id==j]-mu.x[id==j]))%*%matrix(as.numeric(beta+colMeans(uj[id==j,])),nrow=3), sigma, log=TRUE))
  return(L)
}

m.c<- function(beta=model.1[c(1,2),1],
               uj=cbind(mv2i,0) ,m=m,x=x,id=id,npercluster=6,numclusterL2 = 105,
               vare=model.1[4,1]) {
  sigma=diag(npercluster)* vare
  L=lapply(1:numclusterL2,
           function(j) dmvnorm(m[id==j], mean=cbind(1,(x[id==j]-mu.x[id==j]))%*%matrix(as.numeric(beta+colMeans(uj[id==j,])),nrow=2), 
                               sigma, log=TRUE))
  return(L)
}

################ from LL file
#dic1
dic1.my.c <- (sum(ll.summary[ll.summary[,1]==1,5])+sum(ll.summary[ll.summary[,1]==2,5]))*(-4)+
  2*sum(unlist(y.c(beta=model.2[c(1:3),1],
                   uj=cbind(yv,0,0),x=x.im,m=m.im,y=y.im,id=id,npercluster=6,numclusterL2 = 105,
                   vare=model.2[5,1]))+ 
          unlist(m.c(beta=model.1[c(1,2),1],
                     uj=cbind(mv2i,0) ,m=m.im,x=x.im,id=id,npercluster=6,numclusterL2 = 105,
                     vare=model.1[4,1]))  )
#     5578.686

#dic2
dic2.my.c <-  (sum(ll.summary[ll.summary[,1]==1,5])+sum(ll.summary[ll.summary[,1]==2,5]))*(-4)+
  2*sum(ll.summary[ll.summary[,1]==0,3])
#  5448.113

#waic1
waic1.my.c <- -2*sum(ll.summary[ll.summary[,1]==0,3])+2*sum(ll.summary[ll.summary[,1]==0,7])
#   5541.298

#########################
## Model 2
###############################

library("parallel")

#start_time <- Sys.time()
setwd("/Users/handu/Desktop/multiple_imputation/DIC/real data")
###################################################
# read data
###################################################
data <- read.table("employee.dat")
id<-data[,2]
x<-data[,6]
m<-data[,5]
y<-data[,7]
library(mvtnorm)

###################################################
# read posterior summaries
# the structure of the output file is identical to printed tables on output
# row numbers correspond to parameter labels on blimp output
###################################################

summaries <- read.table("summaries.model2.dat")
names(summaries) <- c("mean","mdn","sd","lcl","ucl","psr")

# posterior summaries

# m model 
model.1 <- summaries[c(5,6,1:4),]
row.names(model.1) <- c("m0","mx1","mv2i","mv2is","mv2s","mv1")

# y model
model.2 <- summaries[c(19:21,12:18),]
row.names(model.2) <- c("y0","y1","y2","yv2i","yv2iem","yv2em","yv2ilm","yv2lmem","yv2lm","yv1")


###################################################

ll.summary <- read.csv("ll.model2.csv")

###################################################
# read posterior means of random effects
###################################################
random.summary <- read.table("avgimp.model2.dat")

x.im<- random.summary[,6] 
m.im<- random.summary[,5] 
y.im<- random.summary[,7] 

# m model 
mv2 <- random.summary[,c(10,11)] 
var(mv2)

# y model
yv <- random.summary[,c(12:14)] 
var(yv)

############## cluster means 
mu.m<- random.summary[,16] 
mu.x<- random.summary[,15] 

#########
###### marginal
y.m<- function(beta=model.2[c(1:3),1],
               varu=matrix(c(model.2[4,1],model.2[5,1],model.2[7,1],model.2[5,1],model.2[6,1],
                             model.2[8,1],model.2[7,1],model.2[8,1],model.2[9,1]),nrow=3),
               x=x.im,m=m.im,y=y.im,id=id,npercluster=6,numclusterL2 = 105,
               vare=model.2[10,1]) {
  sigma=diag(npercluster)* vare
  L=lapply(1:numclusterL2,
           function(j) dmvnorm(y[id==j], mean=cbind(1,(m[id==j]-mu.m[id==j]),(x[id==j]-mu.x[id==j]))%*%matrix(as.matrix(beta),nrow=3), 
                               ( sigma+cbind(1,(m[id==j]-mu.m[id==j]),(x[id==j]-mu.x[id==j]))%*%varu%*%t(cbind(1,(m[id==j]-mu.m[id==j]),(x[id==j]-mu.x[id==j])))), log=TRUE))
  return(L)
}
m.m<- function(beta=model.1[c(1,2),1],
               varu=matrix(c(model.1[3,1],model.1[4,1],model.1[4,1],model.1[5,1]),nrow=2),
               m=m.im,x=x.im,id=id,npercluster=6,numclusterL2 = 105,
               vare=model.1[6,1]) {
  sigma=diag(npercluster)* vare
  L=lapply(1:numclusterL2,
           function(j) dmvnorm(m[id==j], mean=cbind(1,(x[id==j]-mu.x[id==j]))%*%matrix(as.matrix(beta),nrow=2), 
                               ( sigma+cbind(1,(x[id==j]-mu.x[id==j]))%*%varu%*%t(cbind(1,(x[id==j]-mu.x[id==j])))), log=TRUE))
  return(L)
}
################ from LL file
#dic1

dic1.my.m <-(sum(ll.summary[ll.summary[,1]==1,4])+sum(ll.summary[ll.summary[,1]==2,4]))*(-4)+
  2*sum(unlist(y.m(beta=model.2[c(1:3),1],
                   varu=matrix(c(model.2[4,1],model.2[5,1],model.2[7,1],model.2[5,1],model.2[6,1],
                                 model.2[8,1],model.2[7,1],model.2[8,1],model.2[9,1]),nrow=3),
                   x=x.im,m=m.im,y=y.im,id=id,npercluster=6,numclusterL2 = 105,
                   vare=model.2[10,1])) + 
          unlist( m.m(beta=model.1[c(1,2),1],
                      varu=matrix(c(model.1[3,1],model.1[4,1],model.1[4,1],model.1[5,1]),nrow=2),
                      m=m.im,x=x.im,id=id,npercluster=6,numclusterL2 = 105,
                      vare=model.1[6,1]))   )

#  5604.815

#dic2
dic2.my.m <-(sum(ll.summary[ll.summary[,1]==1,4])+sum(ll.summary[ll.summary[,1]==2,4]))*(-4)+
  2*sum(ll.summary[ll.summary[,1]==0,2])
#  5507.32

#waic1
waic1.my.m<- -2*sum(ll.summary[ll.summary[,1]==0,2])+2*sum(ll.summary[ll.summary[,1]==0,6])
#   5566.584

#########
###### conditional

y.c<- function(beta=model.2[c(1:3),1],
               uj=yv,x=x.im,m=m.im,y=y.im,id=id,npercluster=6,numclusterL2 = 105,
               vare=model.2[10,1]) {
  sigma=diag(npercluster)* vare
  L=lapply(1:numclusterL2,
           function(j) dmvnorm(y[id==j], mean=cbind(1,(m[id==j]-mu.m[id==j]),(x[id==j]-mu.x[id==j]))%*%matrix(as.numeric(beta+colMeans(uj[id==j,])),nrow=3), sigma, log=TRUE))
  return(L)
}

m.c<- function(beta=model.1[c(1,2),1],
               uj=mv2 ,m=m.im,x=x.im,id=id,npercluster=6,numclusterL2 = 105,
               vare=model.1[6,1]) {
  sigma=diag(npercluster)* vare
  L=lapply(1:numclusterL2,
           function(j) dmvnorm(m[id==j], mean=cbind(1,(x[id==j]-mu.x[id==j]))%*%matrix(as.numeric(beta+colMeans(uj[id==j,])),nrow=2), 
                               sigma, log=TRUE))
  return(L)
}

################ from LL file
#dic1
dic1.my.c <- (sum(ll.summary[ll.summary[,1]==1,5])+sum(ll.summary[ll.summary[,1]==2,5]))*(-4)+
  2*sum(unlist(y.c(beta=model.2[c(1:3),1],
                   uj=yv,x=x.im,m=m.im,y=y.im,id=id,npercluster=6,numclusterL2 = 105,
                   vare=model.2[10,1]))+ 
          unlist(m.c(beta=model.1[c(1,2),1],
                     uj=mv2 ,m=m.im,x=x.im,id=id,npercluster=6,numclusterL2 = 105,
                     vare=model.1[6,1]))  )
#     5513.055

#dic2
dic2.my.c <-  (sum(ll.summary[ll.summary[,1]==1,5])+sum(ll.summary[ll.summary[,1]==2,5]))*(-4)+
  2*sum(ll.summary[ll.summary[,1]==0,3])
#  5358.068

#waic1
waic1.my.c <- -2*sum(ll.summary[ll.summary[,1]==0,3])+2*sum(ll.summary[ll.summary[,1]==0,7])
#   5478.1


#########################
## Model 3
###############################

library("parallel")

#start_time <- Sys.time()
setwd("/Users/handu/Desktop/multiple_imputation/DIC/real data")
###################################################
# read data
###################################################
data <- read.table("employee.dat")
id<-data[,2]
x<-data[,6]
m<-data[,5]
y<-data[,7]
library(mvtnorm)

###################################################
# read posterior summaries
# the structure of the output file is identical to printed tables on output
# row numbers correspond to parameter labels on blimp output
###################################################

summaries <- read.table("summaries.model3.dat")
names(summaries) <- c("mean","mdn","sd","lcl","ucl","psr")

# posterior summaries

# m model 
model.1 <- summaries[c(5:7,1:4),]
row.names(model.1) <- c("m0","mx","mx.mean","mv2i","mv2is","mv2s","mv1")

# y model
model.2 <- summaries[c(21:25,14:20),]
row.names(model.2) <- c("y0","y1","yem,mean","y2","ylm.mean","yv2i","yv2iem","yv2em","yv2ilm","yv2lmem","yv2lm","yv1")


###################################################

ll.summary <- read.csv("ll.model3.csv")

###################################################
# read posterior means of random effects
###################################################
random.summary <- read.table("avgimp.model3.dat")

x.im<- random.summary[,6] 
m.im<- random.summary[,5] 
y.im<- random.summary[,7] 

# m model 
mv2 <- random.summary[,c(10,11)] 
var(mv2)

# y model
yv <- random.summary[,c(12:14)] 
var(yv)

############## cluster means 
mu.m<- random.summary[,16] 
mu.x<- random.summary[,15] 

#########
###### marginal
library(magic)
y.m<- function(beta=model.2[c(1,2,4,3,5),1],
               varu=adiag(matrix(c(model.2[6,1],model.2[7,1],model.2[9,1],model.2[7,1],model.2[8,1],
                                   model.2[10,1],model.2[9,1],model.2[10,1],model.2[11,1]),nrow=3),0,0),
               x=x.im,m=m.im,y=y.im,id=id,npercluster=6,numclusterL2 = 105,
               vare=model.2[12,1]) {
  sigma=diag(npercluster)* vare
  L=lapply(1:numclusterL2,
           function(j) dmvnorm(y[id==j], mean=cbind(1,(m[id==j]-mu.m[id==j]),(x[id==j]-mu.x[id==j]),mu.m[id==j],mu.x[id==j])%*%
                                 matrix(as.matrix(beta),nrow=5), 
                               ( sigma+cbind(1,(m[id==j]-mu.m[id==j]),(x[id==j]-mu.x[id==j]),mu.m[id==j],mu.x[id==j])%*%
                                   varu%*%t(cbind(1,(m[id==j]-mu.m[id==j]),(x[id==j]-mu.x[id==j]),mu.m[id==j],mu.x[id==j]))), log=TRUE))
  return(L)
}
m.m<- function(beta=model.1[c(1:3),1],
               varu=adiag(matrix(c(model.1[4,1],model.1[5,1],model.1[5,1],model.1[6,1]),nrow=2),0),
               m=m.im,x=x.im,id=id,npercluster=6,numclusterL2 = 105,
               vare=model.1[7,1]) {
  sigma=diag(npercluster)* vare
  L=lapply(1:numclusterL2,
           function(j) dmvnorm(m[id==j], mean=cbind(1,(x[id==j]-mu.x[id==j]),mu.x[id==j])%*%matrix(as.matrix(beta),nrow=3), 
                               ( sigma+cbind(1,(x[id==j]-mu.x[id==j]),mu.x[id==j])%*%varu%*%
                                   t(cbind(1,(x[id==j]-mu.x[id==j]),mu.x[id==j]))), log=TRUE))
  return(L)
}
################ from LL file
#dic1

dic1.my.m <-(sum(ll.summary[ll.summary[,1]==1,4])+sum(ll.summary[ll.summary[,1]==2,4]))*(-4)+
  2*sum(unlist(y.m(beta=model.2[c(1,2,4,3,5),1],
                   varu=adiag(matrix(c(model.2[6,1],model.2[7,1],model.2[9,1],model.2[7,1],model.2[8,1],
                                       model.2[10,1],model.2[9,1],model.2[10,1],model.2[11,1]),nrow=3),0,0),
                   x=x.im,m=m.im,y=y.im,id=id,npercluster=6,numclusterL2 = 105,
                   vare=model.2[12,1])) + 
          unlist( m.m(beta=model.1[c(1:3),1],
                      varu=adiag(matrix(c(model.1[4,1],model.1[5,1],model.1[5,1],model.1[6,1]),nrow=2),0),
                      m=m.im,x=x.im,id=id,npercluster=6,numclusterL2 = 105,
                      vare=model.1[7,1]))   )

#  5604.185

#dic2
dic2.my.m <-(sum(ll.summary[ll.summary[,1]==1,4])+sum(ll.summary[ll.summary[,1]==2,4]))*(-4)+
  2*sum(ll.summary[ll.summary[,1]==0,2])
#  5511.442

#waic1
waic1.my.m <- -2*sum(ll.summary[ll.summary[,1]==0,2])+2*sum(ll.summary[ll.summary[,1]==0,6])
#    5570.962

#########
###### conditional

y.c<- function(beta=model.2[c(1,2,4,3,5),1],
               uj=cbind(yv,0,0),x=x.im,m=m.im,y=y.im,id=id,npercluster=6,numclusterL2 = 105,
               vare=model.2[12,1]) {
  sigma=diag(npercluster)* vare
  L=lapply(1:numclusterL2,
           function(j) dmvnorm(y[id==j], mean=cbind(1,(m[id==j]-mu.m[id==j]),(x[id==j]-mu.x[id==j]),mu.m[id==j],mu.x[id==j])%*%
                                 matrix(as.numeric(beta+colMeans(uj[id==j,])),nrow=5), sigma, log=TRUE))
  return(L)
}

m.c<- function(beta=model.1[c(1:3),1],
               uj=cbind(mv2,0) ,m=m.im,x=x.im,id=id,npercluster=6,numclusterL2 = 105,
               vare=model.1[7,1]) {
  sigma=diag(npercluster)* vare
  L=lapply(1:numclusterL2,
           function(j) dmvnorm(m[id==j], mean=cbind(1,(x[id==j]-mu.x[id==j]),mu.x[id==j])%*%matrix(as.numeric(beta+colMeans(uj[id==j,])),nrow=3), 
                               sigma, log=TRUE))
  return(L)
}

################ from LL file
#dic1
dic1.my.c <- (sum(ll.summary[ll.summary[,1]==1,5])+sum(ll.summary[ll.summary[,1]==2,5]))*(-4)+
  2*sum(unlist(y.c(beta=model.2[c(1,2,4,3,5),1],
                   uj=cbind(yv,0,0),x=x.im,m=m.im,y=y.im,id=id,npercluster=6,numclusterL2 = 105,
                   vare=model.2[12,1]))+ 
          unlist(m.c(beta=model.1[c(1:3),1],
                     uj=cbind(mv2,0) ,m=m.im,x=x.im,id=id,npercluster=6,numclusterL2 = 105,
                     vare=model.1[7,1]))  )
#     5513.22

#dic2
dic2.my.c <-  (sum(ll.summary[ll.summary[,1]==1,5])+sum(ll.summary[ll.summary[,1]==2,5]))*(-4)+
  2*sum(ll.summary[ll.summary[,1]==0,3])
#   5366.736

#waic1
waic1.my.c <- -2*sum(ll.summary[ll.summary[,1]==0,3])+2*sum(ll.summary[ll.summary[,1]==0,7])
#  5488.689
