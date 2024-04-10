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
library(bvhar)



your_wd = "/.../project"

Rcpp::sourceCpp(paste0(your_wd,"/src/avirf.cpp"))

data = read.csv(paste0(your_wd,"/data/data.csv"))


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
data_m11 = VAR(data,p=11)
data_m12 = VAR(data,p=12)
data_m13 = VAR(data,p=13)
data_m14 = VAR(data,p=14)
data_m15 = VAR(data,p=15)
data_m16 = VAR(data,p=16)
data_m17 = VAR(data,p=17)
data_m18 = VAR(data,p=18)
data_m19 = VAR(data,p=19)
data_m20 = VAR(data,p=20)
data_m21 = VAR(data,p=21)
data_m22 = VAR(data,p=22)
data_m23 = VAR(data,p=23)




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
print(c(data_m11$aic,data_m11$bic,data_m11$hq))
print(c(data_m12$aic,data_m12$bic,data_m12$hq))
print(c(data_m13$aic,data_m13$bic,data_m13$hq))
print(c(data_m14$aic,data_m14$bic,data_m14$hq))
print(c(data_m15$aic,data_m15$bic,data_m15$hq))
print(c(data_m16$aic,data_m16$bic,data_m16$hq))
print(c(data_m17$aic,data_m17$bic,data_m17$hq))
print(c(data_m18$aic,data_m18$bic,data_m18$hq))
print(c(data_m19$aic,data_m19$bic,data_m19$hq))
print(c(data_m20$aic,data_m20$bic,data_m20$hq))
print(c(data_m21$aic,data_m21$bic,data_m21$hq))
print(c(data_m22$aic,data_m22$bic,data_m22$hq))
print(c(data_m23$aic,data_m23$bic,data_m23$hq))


data_m3_1 = refVAR(data_m3, thres=1.96)
data_m4_1 = refVAR(data_m4, thres=1.96)
data_m5_1 = refVAR(data_m5, thres=1.96)
data_m6_1 = refVAR(data_m6, thres=1.96)
data_m7_1 = refVAR(data_m7, thres=1.96)

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

print(c(data_m2_1$aic,data_m2_1$bic,data_m2_1$hq))
print(c(data_m3_1$aic,data_m3_1$bic,data_m3_1$hq))
print(c(data_m4_1$aic,data_m4_1$bic,data_m4_1$hq))
print(c(data_m5_1$aic,data_m5_1$bic,data_m5_1$hq))
print(c(data_m6_1$aic,data_m6_1$bic,data_m6_1$hq))
print(c(data_m7_1$aic,data_m7_1$bic,data_m7_1$hq))



summary(data_m1)
#Extracting residuals
resi = residuals(data_m3)
mq(resi, adj=4^2 *3)

obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = T))
x1 <- bekk_fit(obj_spec, resi, QML_t_ratios = T, crit = 1e-09,max_iter = 500)

x2 = virf(x1,index_series = 3,q = 0.99,n.ahead=100)

plot(x2)

summary(x1)

print(x1$A_t)

print(x1$B_t)

print(x1$G_t)


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


#(A)VIRF
#Configuring some parameters
#Horizon
symmetric = FALSE
N=dim(data)[2]
x=x1
n.ahead <- 100
#Trajectories
n.iterations <- 100000
#Percentile
q = 0.99
#1-4 denotes the national CETs, GECs, Hubei, and Guangdong CETs, respectively.
index_series = 1
#Starting time
time = 1
#Duration
span = dim(x1$H_t)[1]
span = 2
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
V <- avirf_bekk(x1,time,span-1,q,index_series,n.ahead,n.iterations,symmetric)
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

##co-volatility
eps_mean = apply(x1$e_t,MARGIN=2,mean)
x1$A[3,3] * x1$A[4,4] * eps_mean[3]
x1$A[3,3] * x1$A[4,4] * eps_mean[4]
