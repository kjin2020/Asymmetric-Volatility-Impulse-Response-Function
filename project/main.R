library(xts)
library(vars)
library(BEKKs)
library(ggplot2)
library(tidyr)
library(dplyr)
library(Rcpp)
library(aod)


your_wd = "/.../project"

Rcpp::sourceCpp(paste0(your_wd,"/src/avirf.cpp"))

data = read.csv(paste0(your_wd,"/data/data.csv"))

#VAR-BEKK-GARCH
#data_m1 = VAR(data,p=1)
data_m2 = VAR(data,p=2)
#data_m3 = VAR(data,p=3)
#data_m4 = VAR(data,p=4)

#data_m5 = VAR(data,p=5)

resi = residuals(data_m2)
#mq(resi, adj=4^2 * 2)

obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = T))

x1 <- bekk_fit(obj_spec, resi, QML_t_ratios = T, crit = 1e-09)

summary(x1)

#Wald test
coef_cov <- function(theta, r, signs) {
  l1 = BEKKs:::score_asymm_bekk(theta, r, signs)
  l1 = crossprod(l1)
  J1 = BEKKs:::hesse_asymm_bekk(theta, r, signs)
  cov_mat = solve(J1) %*% l1 %*% solve(J1)
  return(cov_mat)
}
cov_matrix = coef_cov(x1$theta,x1$data,matrix(rep(-1, N), ncol = 1))

A_index_1 = 10 + 12
A_index_2 = A_index_1 + 3

G_index_1 = 10 + 16 * 2 + 12
G_index_2 = G_index_1 + 3

B_index_1 = 10 + 16 * 1 + 12
B_index_2 = B_index_1 + 3

#H_0
wald.test(Sigma=cov_matrix, b=x1$theta, Terms = c(A_index_2,A_index_1,G_index_2,G_index_1))
#H_1
wald.test(Sigma=cov_matrix, b=x1$theta, Terms = c(B_index_1,B_index_2))


#(A)VIRF

n.ahead <- 100
n.iterations <- 2000
q = 0.99
index_series = 2
time = 1
span = dim(x1$H_t)[1]
#Generate a shocks matrix
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

df <- data.frame(
  periods   = 1:n.ahead, 
  CETs      = V[,1], 
  GECs      = V[,5],
  HBCM      = V[,8],
  GDCM      = V[,10]
)

#write.csv(df,file=paste0(your_wd,"/data/test.csv"),row.names = F)

df1 = read.csv(paste0(your_wd,"/data/99GECs.csv"))
df2 = read.csv(paste0(your_wd,"/data/01GECs.csv"))

df_long_1 <- tidyr::pivot_longer(df1, cols = c("HBCM","GDCM","CETs","GECs"), 
                                 names_to = "variable", values_to = "value")
df_long_2 <- tidyr::pivot_longer(df2, cols = c("HBCM","GDCM","CETs","GECs"), 
                                 names_to = "variable", values_to = "value")

my_labels_1 <- c(GECs = "Response in the variance \n of M2 to positive news in M1",CETs = "Response in the variance \n of M1 to positive news in M1",GDCM="Response in the variance \n of M4 to positive news in M1",HBCM="Response in the variance \n of M3 to positive news in M1")
my_labels_2 <- c(GECs = "Response in the variance \n of M2 to negative news in M1",CETs = "Response in the variance \n of M1 to negative news in M1",GDCM="Response in the variance \n of M4 to negative news in M1",HBCM="Response in the variance \n of M3 to negative news in M1")

p <- ggplot(df_long_1, aes(x = periods, y = value)) + 
  geom_line(color = "black") + 
  facet_wrap(~ variable, scales = "free_y", ncol = 1,labeller = labeller(variable = my_labels_1)) +
  scale_x_continuous(limits=c(0,100),expand=c(0,0))+
  scale_y_continuous(limits=c(min(df_long_1$value,df_long_2$value),max(df_long_1$value,df_long_2$value)),expand=c(0,0))+
  theme_minimal() + 
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
    aspect.ratio = 0.5,
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(face = "bold"), 
    plot.background = element_rect(fill = "white", color = "white"), 
    panel.border = element_rect(color = "black", fill = NA)  
  ) +
  labs(x = NULL, y = NULL)

p <- ggplot(df_long_2, aes(x = periods, y = value)) + 
  geom_line(color = "black") + 
  facet_wrap(~ variable, scales = "free_y", ncol = 1,labeller = labeller(variable = my_labels_2)) +
  scale_x_continuous(limits=c(0,100),expand=c(0,0))+
  scale_y_continuous(limits=c(min(df_long_1$value,df_long_2$value),max(df_long_1$value,df_long_2$value)),expand=c(0,0))+
  theme_minimal() + 
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
    aspect.ratio = 0.5,
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(face = "bold"), 
    plot.background = element_rect(fill = "white", color = "white"), 
    panel.border = element_rect(color = "black", fill = NA)  
  ) +
  labs(x = NULL, y = NULL)  
print(p)

