library(timeSeries)
library(zoo)
library(tseries)
library(xts)
library(vars)
library(BEKKs)
library(MTS)
library(ggplot2)
library(MASS)
library(tidyr)
#library(aod)



your_wd = "/.../project"

Rcpp::sourceCpp(paste0(your_wd,"/src/avirf.cpp"))

data = read.csv(paste0(your_wd,"/data/data.csv"))


#VAR-BEKK-GARCH
##Estimating conditional mean with VAR
VARselect(data)
###AIC(n)  HQ(n)  SC(n) FPE(n) 
###3      1      1      3 
###LR test (PAC test Lütkepohl, H. (1985))
####df = 4**2 = 16
####LR1= 32
####LR5= 26.296
sigma1 = log(det(cov(residuals(VAR(data,p=1)))))
sigma2 = log(det(cov(residuals(VAR(data,p=2)))))
sigma3 = log(det(cov(residuals(VAR(data,p=3)))))
sigma4 = log(det(cov(residuals(VAR(data,p=4)))))
sigma5 = log(det(cov(residuals(VAR(data,p=5)))))
sigma6 = log(det(cov(residuals(VAR(data,p=6)))))
sigma7 = log(det(cov(residuals(VAR(data,p=7)))))

iT= dim(data)[1]
k = dim(data)[2]

LR_2 = iT * (sigma1 - sigma2)
LR_3 = iT * (sigma2 - sigma3)#***
LR_4 = iT * (sigma3 - sigma4)
LR_5 = iT * (sigma4 - sigma5)
LR_6 = iT * (sigma5 - sigma6)
LR_7 = iT * (sigma6 - sigma7)

###Granger causality test
colnames(data) = c("ncets","gdcets","gdele","GECs")

create_lag <- function(x, k){
  c(rep(NA, k), x[1:(length(x)-k)])
}

max_lag <- 4
for(j in 1:max_lag){
  data <- cbind(data, create_lag(data[,"ncets"], k=j))
  data <- cbind(data, create_lag(data[,"gdcets"], k=j))
  data <- cbind(data, create_lag(data[,"gdele"], k=j))
  data <- cbind(data, create_lag(data[,"GECs"], k=j))
}
colnames(data) <- c(
  "ncets","gdcets","gdele","GECs",
  paste(rep(c("ncets","gdcets","gdele","GECs"), 4), rep(1:4, each=4), sep="_"))

data = as.data.frame(data)

lm_ncets_1 = lm(ncets ~ ncets_1 + ncets_2 + ncets_3 + gdcets_1 + gdcets_2 + gdcets_3 + gdele_1 + gdele_2 + gdele_3 + GECs_1 + GECs_2 + GECs_3, data=data)
lm_ncets_2 = lm(ncets ~ ncets_1 + ncets_2 + ncets_3 + gdele_1 + gdele_2 + gdele_3 + GECs_1 + GECs_2 + GECs_3, data=data)
anova(lm_ncets_1, lm_ncets_2)#***
lm_ncets_3 = lm(ncets ~ ncets_1 + ncets_2 + ncets_3 + gdcets_1 + gdcets_2 + gdcets_3 + GECs_1 + GECs_2 + GECs_3, data=data)
anova(lm_ncets_1, lm_ncets_3)
lm_ncets_4 = lm(ncets ~ ncets_1 + ncets_2 + ncets_3 + gdcets_1 + gdcets_2 + gdcets_3 + gdele_1 + gdele_2 + gdele_3, data=data)
anova(lm_ncets_1, lm_ncets_4)

lm_gdcets_1 = lm(gdcets ~ ncets_1 + ncets_2 + ncets_3 + gdcets_1 + gdcets_2 + gdcets_3 + gdele_1 + gdele_2 + gdele_3 + GECs_1 + GECs_2 + GECs_3, data=data)
lm_gdcets_2 = lm(gdcets ~ gdcets_1 + gdcets_2 + gdcets_3 + gdele_1 + gdele_2 + gdele_3 + GECs_1 + GECs_2 + GECs_3, data=data)
anova(lm_gdcets_1, lm_gdcets_2)#*
lm_gdcets_3 = lm(gdcets ~ ncets_1 + ncets_2 + ncets_3 + gdcets_1 + gdcets_2 + gdcets_3 + GECs_1 + GECs_2 + GECs_3, data=data)
anova(lm_gdcets_1, lm_gdcets_3)
lm_gdcets_4 = lm(gdcets ~ ncets_1 + ncets_2 + ncets_3 + gdcets_1 + gdcets_2 + gdcets_3 + gdele_1 + gdele_2 + gdele_3, data=data)
anova(lm_gdcets_1, lm_gdcets_4)

lm_gdele_1 = lm(gdele ~ ncets_1 + ncets_2 + ncets_3 + gdcets_1 + gdcets_2 + gdcets_3 + gdele_1 + gdele_2 + gdele_3 + GECs_1 + GECs_2 + GECs_3, data=data)
lm_gdele_2 = lm(gdele ~ gdcets_1 + gdcets_2 + gdcets_3 + gdele_1 + gdele_2 + gdele_3 + GECs_1 + GECs_2 + GECs_3, data=data)
anova(lm_gdele_1, lm_gdele_2)
lm_gdele_3 = lm(gdele ~ ncets_1 + ncets_2 + ncets_3 + gdele_1 + gdele_2 + gdele_3 + GECs_1 + GECs_2 + GECs_3, data=data)
anova(lm_gdele_1, lm_gdele_3)
lm_gdele_4 = lm(gdele ~ ncets_1 + ncets_2 + ncets_3 + gdcets_1 + gdcets_2 + gdcets_3 + gdele_1 + gdele_2 + gdele_3, data=data)
anova(lm_gdele_1, lm_gdele_4)

lm_GECs_1 = lm(GECs ~ ncets_1 + ncets_2 + ncets_3 + gdcets_1 + gdcets_2 + gdcets_3 + gdele_1 + gdele_2 + gdele_3 + GECs_1 + GECs_2 + GECs_3, data=data)
lm_GECs_2 = lm(GECs ~ gdcets_1 + gdcets_2 + gdcets_3 + gdele_1 + gdele_2 + gdele_3 + GECs_1 + GECs_2 + GECs_3, data=data)
anova(lm_GECs_1, lm_GECs_2) #**
lm_GECs_3 = lm(GECs ~ ncets_1 + ncets_2 + ncets_3 + gdele_1 + gdele_2 + gdele_3 + GECs_1 + GECs_2 + GECs_3, data=data)
anova(lm_GECs_1, lm_GECs_3)
lm_GECs_4 = lm(GECs ~ ncets_1 + ncets_2 + ncets_3 + gdcets_1 + gdcets_2 + gdcets_3 + GECs_1 + GECs_2 + GECs_3, data=data)
anova(lm_GECs_1, lm_GECs_4)

###Construct VAR based on Granger causality test result
lm_ncets = lm(ncets ~ ncets_1  + ncets_2 + ncets_3 + gdcets_1 + gdcets_2 + gdcets_3, data=data)
lm_gdcets= lm(gdcets~ gdcets_1 + gdcets_2+ gdcets_3, data=data)
lm_gdele = lm(gdele ~ gdele_1  + gdele_2 + gdele_3, data=data)
lm_GECs  = lm(GECs  ~  ncets_1 + ncets_2 + ncets_3 + GECs_1 + GECs_2 + GECs_3, data=data)

mq(resi, adj=6*2+3*2)

###Extracting residuals
res_ncets = residuals(lm_ncets)
res_gdcets= residuals(lm_gdcets)
res_gdele = residuals(lm_gdele)
res_GECs  = residuals(lm_GECs)

resi = cbind(res_ncets,res_gdcets,res_gdele,res_GECs)
####Multivariate Ljung-Box test
mq(resi, adj=6*2+3*2)

#BEKK-GARCH
obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = F))

x1 <- bekk_fit(obj_spec, resi, QML_t_ratios = T, crit = 1e-09,max_iter = 500)

print(x1$A_t)

portmanteau.test(x1,5)

portmanteau.test(x1,15)

portmanteau.test(x1,20)

portmanteau.test(x1,50)

#Wald test
#coef_cov <- function(theta, r, signs) {
#  l1 = BEKKs:::score_asymm_bekk(theta, r, signs)
#  l1 = crossprod(l1)
#  J1 = BEKKs:::hesse_asymm_bekk(theta, r, signs)
#  cov_mat = solve(J1) %*% l1 %*% solve(J1)
#  return(cov_mat)
#}
#cov_matrix = coef_cov(x1$theta,x1$data,matrix(rep(-1, N), ncol = 1))
#Indexes of the A, G, and B matrices
#A_index_1 = 10 + 12
#A_index_2 = A_index_1 + 3

#G_index_1 = 10 + 16 * 2 + 12
#G_index_2 = G_index_1 + 3

#B_index_1 = 10 + 16 * 1 + 12
#B_index_2 = B_index_1 + 3

#H_0
#wald.test(Sigma=cov_matrix, b=x1$theta, Terms = c(A_index_2,A_index_1,G_index_2,G_index_1))
#H_1
#wald.test(Sigma=cov_matrix, b=x1$theta, Terms = c(B_index_1,B_index_2))

##co-volatility
eps_mean = apply(x1$e_t,MARGIN=2,mean)
x1$A[3,3] * x1$A[4,4] * eps_mean[3]
x1$A[3,3] * x1$A[4,4] * eps_mean[4]

#(A)VIRF
#Configuring some parameters
#Horizon
symmetric = TRUE
N=dim(resi)[2]
x=x1
n.ahead <- 100
#Trajectories
n.iterations <- 100000
#Percentile
q = 0.99
#1-4 denotes the national CETs, GECs, Hubei, and Guangdong CETs, respectively.
index_series = 4
#Starting time
time = 1
#Duration
span = dim(x1$H_t)[1]
#span = 2
#Generate shock matrices, \(z_1\) and \(z_2\)
shocks_mat <- function(x,q,index_series,n.ahead,n.iterations){
  z  = x$e_t
  N  = dim(z)[2]
  z1 = matrix(nrow=n.ahead,ncol=N*n.iterations)
  z2 = matrix(nrow=n.ahead,ncol=N*n.iterations)
  shocks = matrix(c(rep(0,index_series-1),quantile(z[,index_series],q),rep(0,N-index_series)),nrow=1)
  z1[1,] = c(rep(0,(index_series-1)*n.iterations),rep(quantile(z[,index_series],q),n.iterations),rep(0,(N-index_series)*n.iterations))
  for (i in 1:N){
    z1[2:n.ahead,((i-1)*n.iterations+1):(i*n.iterations)] = matrix(sample(z[,i],size=(n.ahead-1)*n.iterations,replace=T),nrow=n.ahead-1)
    z2[,((i-1)*n.iterations+1):(i*n.iterations)] = matrix(sample(z[,i],size=n.ahead*n.iterations,replace=T),nrow=n.ahead)
  }
  return(list(shocks=shocks,z1=z1,z2=z2))
}
#Estimating AVIRF
avirf_bekk <- function(x,time,span=1,q,index_series,n.ahead,n.iterations,symmetric){
  N = dim(x$A)[1]
  H_1 = matrix(NA,N,N*dim(x$H_t)[1])
  for(i in 1:dim(x$H_t)[1]){
    H_1[,(N*i-N+1):(N*i)] = matrix(x$H_t[i,],N,N)
  }
  if(symmetric==TRUE){
    shocks = matrix(c(rep(0,index_series-1),quantile(x$e_t[,index_series],q),rep(0,N-index_series)),nrow=1)
    return(virf_bekk(time,time+span,H_1,x$A,x$G,t(x$C0)%*%x$C0,shocks,n.ahead,N))
  } else {
    z = shocks_mat(x,q,index_series,n.ahead,n.iterations)
    return(virf_bekka(time,time+span,H_1,x$A,x$B,x$G,t(x$C0)%*%x$C0,z$shocks,z$z1,z$z2,n.ahead,n.iterations,N,matrix(rep(1, N),ncol = 1),matrix(rep(0, N), ncol = 1)))
  }
}

start_time <- Sys.time()
V          <- avirf_bekk(x1,time,span-1,q,index_series,n.ahead,n.iterations,symmetric)
end_time   <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

idx = c(1,rep(0,N-1))
var_names = paste0('Response in the variance \n of M',1:N,' to ',ifelse(q>0.5,'positive','negative'),' news in M',index_series)
for (i in 2:N){
  idx[i] = idx[i-1]+(N-i+2)
}

df <- as.data.frame(V[,idx])
colnames(df) = var_names
df['periods'] = 1:n.ahead
df_long_tmp <- tidyr::pivot_longer(df, cols = var_names,
                                 names_to = "variable", values_to = "value")

p <- ggplot(df_long_tmp, aes(x = periods, y = value)) +
  geom_line(color = "black") +
  geom_hline(yintercept=0,size=0.75) +
  facet_wrap(~ variable,ncol = 1,scales="free_y") +
  scale_x_continuous(limits=c(1,n.ahead),expand=c(0,0))+
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
    aspect.ratio = 0.075,
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(x = NULL, y = NULL)
print(p)

#write.csv(df,file=paste0(your_wd,"/data/",ifelse(symmetric==TRUE,'symmetric','asymmetric'),"/Response_to_",ifelse(q>0.5,'positive','negative'),"_news_in_M",index_series,".csv"),row.names = F)
#write.csv(data,file=paste0(your_wd,"/data/data.csv"),row.names = F)

