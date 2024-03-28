library(timeSeries)
library(zoo)
library(tseries)
library(xts)
library(openxlsx)
library(vars)
library(BEKKs)
library(MTS)
library(ggplot2)
library(MASS)
library(tidyr)
library(aod)



your_wd = "/.../project"

Rcpp::sourceCpp(paste0(your_wd,"/src/avirf.cpp"))

#data = read.csv(paste0(your_wd,"/data/data.csv"))

#VAR-BEKK-GARCH
#Estimating conditional mean with VAR
data_m1 = VAR(data,p=1)
data_m2 = VAR(data,p=2)
data_m3 = VAR(data,p=3)
data_m4 = VAR(data,p=4)
data_m5 = VAR(data,p=5)
data_m6 = VAR(data,p=6)
data_m7 = VAR(data,p=7)
data_m8 = VAR(data,p=8)
data_m9 = VAR(data,p=9)


data_m10 = VAR(data,p=10)



print(c(data_m1$aic,data_m1$bic,data_m1$hq))
print(c(data_m2$aic,data_m2$bic,data_m2$hq))
print(c(data_m3$aic,data_m3$bic,data_m3$hq))
print(c(data_m4$aic,data_m4$bic,data_m4$hq))
print(c(data_m5$aic,data_m5$bic,data_m5$hq))
print(c(data_m6$aic,data_m6$bic,data_m6$hq))
print(c(data_m7$aic,data_m7$bic,data_m7$hq))
print(c(data_m8$aic,data_m8$bic,data_m8$hq))
print(c(data_m9$aic,data_m9$bic,data_m9$hq))
print(c(data_m10$aic,data_m10$bic,data_m10$hq))

summary(data_m2)
#Extracting residuals
resi = residuals(data_m2)
mq(resi, adj=4^2 *2)

obj_spec <- bekk_spec(model = list(type = "dbekk", asymmetric = T))

x1 <- bekk_fit(obj_spec, resi, QML_t_ratios = T, crit = 1e-09,max_iter = 500)

#x2 = virf(x1,index_series = 1,q = 0.5)

#plot(x2)

summary(x1)

print(x1$A_t)

print(x1$B_t)

print(x1$G_t)


portmanteau.test(x1,5)

portmanteau.test(x1,15)

portmanteau.test(x1,20)


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


#(A)VIRF
#Configuring some parameters
#Horizon
N=dim(data)[2]
x=x1
n.ahead <- 50
#Trajectories
n.iterations <- 200000
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
avirf_bekk <- function(x,time,span=1,q,index_series,n.ahead,n.iterations){
  N = dim(x$A)[1]
  H_1 = matrix(NA,N,N*dim(x$H_t)[1])
  for(i in 1:dim(x$H_t)[1]){
    H_1[,(N*i-N+1):(N*i)] = matrix(x$H_t[i,],N,N)
  }
  z = shocks_mat(x,q,index_series,n.ahead,n.iterations)
  return(virf_bekka(time,time+span,H_1,x$A,x$B,x$G,t(x$C0)%*%x$C0,z$shocks,z$z1,z$z2,n.ahead,n.iterations,N,matrix(rep(1, N),ncol = 1),matrix(rep(0, N), ncol = 1)))
}

start_time <- Sys.time()
V <- avirf_bekk(x1,time,span-1,q,index_series,n.ahead,n.iterations)
end_time <- Sys.time()
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
  geom_line(color = "blue") +
  facet_wrap(~ variable,ncol = 1,scales="free_y") +
  scale_x_continuous(limits=c(0,n.ahead),expand=c(0,0))+
  scale_y_continuous(expand=c(0.0001,0.0001))+
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
    aspect.ratio = 0.35,
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(x = NULL, y = NULL)
print(p)

write.csv(df,file=paste0(your_wd,"/data/Response_to_",ifelse(q>0.5,'positive','negative'),"_news_in_M",index_series,".csv"),row.names = F)
