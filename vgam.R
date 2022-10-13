###############################################################################
#############--------------- MyCopula.R ------------################
###############################################################################
#
# Des:
# Ref:
#
# Team  : AI Team
# Author: Ding Weijie 
# Email : 380666903@qq.com
# Date  : 2021-03-09
###

#引入copula包
library(copula)
library(psych)
library(VineCopula)

#随机生成数据集
X_1 <- runif(100, 0, 100)
X_2 <-X_1+ runif(100, 0, 50)

#初步研究两个数据的关系
plot(X_1,X_2)
abline(lm(X_2~X_1),col='red',lwd=1)
cor(X_1,X_2,method='spearman')

#以下4行代码是Copula函数的非参数估计，在不知道哪一个Copula函数下使用。
#u <- pobs(as.matrix(cbind(X_1,X_2)))[,1]
#v <- pobs(as.matrix(cbind(X_1,X_2)))[,2]
#selectedCopula <- BiCopSelect(u,v,familyset=NA)
#selectedCopula


#将数据映射到累计概率分布值上，因为copula函数的自变量是对应边际分布的累计概率分布值。这个量是服从[0,1]上的均匀分布的。

##########################################################################################
#高斯Copula函数
gaussian.cop <- normalCopula(dim=2) # 生成二维高斯Copula函数
set.seed(500)



#使用pobs函数将数据映射到累计概率分布值上，因为copula函数的自变量是对应边际分布的累计概率分布值。这个量是服从[0,1]上的均匀分布的。
m <- pobs(as.matrix(cbind(X_1,X_2)))
fit <- fitCopula(gaussian.cop,m,method='ml')#用数据集拟合高斯Copula函数，对其参数进行估计
# 获得参数估计
coef(fit) 
rho <- coef(fit)[1]

#画出拟合后的联合概率密度函数
persp(normalCopula(dim=2,rho),dCopula)

#画出拟合后的联合概率密度函数在二维的投影。
u <- rCopula(3965,normalCopula(dim=2,rho))
plot(u[,1],u[,2],pch='.',col='blue')


#使用以上高斯copula函数，对数据的分布进行联合（这里假设原始数据的分布是付出正太分布的）。
cor(u,method='spearman')
X_1_mu <- mean(X_1)
X_1_sd <- sd(X_1)
X_2_mu <- mean(X_2)
X_2_sd <- sd(X_2)
copula_dist <- mvdc(copula=normalCopula(rho,dim=2), margins=c("norm","norm"),
                    paramMargins=list(list(mean=X_1_mu, sd=X_1_sd),
                                      list(mean=X_2_mu, sd=X_2_sd)))
# 对联合后的分布，进行模拟，生成模拟数据，并与远数据比较。
sim <- rMvdc(3965,copula_dist)
plot(X_1,X_2,main='relation')
points(sim[,1],sim[,2],col='red',pch='.')
legend('bottomright',c('Observed','Simulated'),col=c('black','red'),pch=21)


##########################################################################################
t.cop <- tCopula(dim=2)
set.seed(500)
m <- pobs(as.matrix(cbind(X_1,X_2)))
fit <- fitCopula(t.cop,m,method='ml')
coef(fit)
rho <- coef(fit)[1]
df <- coef(fit)[2]
persp(tCopula(dim=2,rho,df=df),dCopula)


u <- rCopula(3965,tCopula(dim=2,rho,df=df))
plot(u[,1],u[,2],pch='.',col='blue')
cor(u,method='spearman')
X_1_mu <- mean(X_1)
X_1_sd <- sd(X_1)
X_2_mu <- mean(X_2)
X_2_sd <- sd(X_2)

copula_dist <- mvdc(copula=tCopula(rho,dim=2,df=df), margins=c("norm","norm"),
                    paramMargins=list(list(mean=X_1_mu, sd=X_1_sd),
                                      list(mean=X_2_mu, sd=X_2_sd)))
sim <- rMvdc(3965,copula_dist)
plot(X_1,X_2,main='relation')
points(sim[,1],sim[,2],col='red',pch='.')
legend('bottomright',c('Observed','Simulated'),col=c('black','red'),pch=21)





##########################################################################################
clayton.cop <- claytonCopula(dim=2)
set.seed(500)
m <- pobs(as.matrix(cbind(X_1,X_2)))
fit <- fitCopula(clayton.cop,m,method='ml')
coef(fit)
alpha <- coef(fit)[1]

persp(claytonCopula(dim=2,alpha),dCopula)


u <- rCopula(3965,claytonCopula(dim=2,alpha))
plot(u[,1],u[,2],pch='.',col='blue')
cor(u,method='spearman')
X_1_mu <- mean(X_1)
X_1_sd <- sd(X_1)
X_2_mu <- mean(X_2)
X_2_sd <- sd(X_2)

copula_dist <- mvdc(copula=claytonCopula(dim=2,alpha), margins=c("norm","norm"),
                    paramMargins=list(list(mean=X_1_mu, sd=X_1_sd),
                                      list(mean=X_2_mu, sd=X_2_sd)))
sim <- rMvdc(3965,copula_dist)
plot(X_1,X_2,main='relation')
points(sim[,1],sim[,2],col='red',pch='.')
legend('bottomright',c('Observed','Simulated'),col=c('black','red'),pch=21)




##########################################################################################
gumbel.cop <- gumbelCopula(dim=2)
set.seed(500)
m <- pobs(as.matrix(cbind(X_1,X_2)))
fit <- fitCopula(gumbel.cop,m,method='ml')
coef(fit)
alpha <- coef(fit)[1]

persp(gumbelCopula(dim=2,alpha),dCopula)


u <- rCopula(3965,gumbelCopula(dim=2,alpha))
plot(u[,1],u[,2],pch='.',col='blue')
cor(u,method='spearman')
X_1_mu <- mean(X_1)
X_1_sd <- sd(X_1)
X_2_mu <- mean(X_2)
X_2_sd <- sd(X_2)

copula_dist <- mvdc(copula=gumbelCopula(dim=2,alpha), margins=c("norm","norm"),
                    paramMargins=list(list(mean=X_1_mu, sd=X_1_sd),
                                      list(mean=X_2_mu, sd=X_2_sd)))
sim <- rMvdc(3965,copula_dist)
plot(X_1,X_2,main='relation')
points(sim[,1],sim[,2],col='red',pch='.')
legend('bottomright',c('Observed','Simulated'),col=c('black','red'),pch=21)

####


## qita
library(VGAM)
nn <- 1000; set.seed(1)
bdata <- data.frame(x2 = runif(nn))
bdata <- transform(bdata,
                   pobs0 = logitlink(-2 + x2, inverse = TRUE),
                   pobs1 = logitlink(-2 + x2, inverse = TRUE))
bdata <- transform(bdata,
                   y1 = rzoabeta(nn, shape1 = exp(1 + x2), shape2 = exp(2 - x2),
                                 pobs0 = pobs0, pobs1 = pobs1))
summary(bdata)
fit2 <- vglm(y1 ~ x2, zoabetaR(parallel.pobs = FALSE),
             data = bdata, trace = TRUE)
coef(fit2, matrix = TRUE)
summary(fit2)

data1 <- data0[[2]]
data2 <- data1[,c(101,301)]
fit2 <- vglm(X1 ~ X1.1, zoabetaR(parallel.pobs = TRUE),
             data = data2, trace = TRUE)


nrow(bdata[bdata$y1 == 0,])

library(gamlss)
gam <- gamlss(data2[,1] ~ data2[,2], data = data2, family = BEZI(), trace = FALSE)
a <- summary(gam)
library(VGAM)
vglm <- vglm(data1[,1] ~ data1[,201], zoabetaR(),trace = FALSE, data = data1)
summary(vglm)
a <- summaryvglm(vglm,signif.stars = TRUE)
coef(vglm, matrix = TRUE)
hdeff.vglm(vglm)
##
library(zoib)
m1 <- zoib(data1[,1] ~ data1[,201], random = 1,    
           EUID = OPdat$species,  zero.inflation = TRUE,  
           one.inflation = TRUE, data = OPdat,   n.iter=50, 
           n.thin=20, n.burn=50)

library(copula)

shape1=0.69647651
rate1=0.010729039
shape2=0.9635694881
rate2=0.0021278487
n<-500
cop<-claytonCopula(5,dim=2)
u<-rgamma(n,shape=shape1,scale=rate1)
v<-rgamma(n,shape=shape2,scale=rate2)
g<-fitCopula(cbind(pgamma(u,shape=shape1,scale=rate1),pgamma(v,shape=shape2,scale=rate2)),cop,optim.control=list(trace=TRUE,maxit=2000),method="ml")
iTau(gumbel.cop, cor(cbind(u,v), method="kendall"))
Cop.gumbel<-archmCopula(family="gumbel",dim=2,param=cg)  #gumbel
myMvdg<-mvdc(copula=Cop.gumbel,margins=c("gamma","gamma"),paramMargins=list(list(shape=shape1,rate=rate1),list(shape=shape2,rate=rate2)))

