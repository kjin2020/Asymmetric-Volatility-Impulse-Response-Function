library(xts)
library(vars)
library(BEKKs)
library(ggplot2)
library(tidyr)
library(dplyr)
library(Rcpp)

your_wd = "/.../project"

Rcpp::sourceCpp(paste0(your_wd,"/src/avirf.cpp"))

data = read.csv(paste0(your_wd,"/data/data.csv"))

#VAR-BEKK-GARCH
#data_m1 = VAR(data,p=1)
data_m2 = VAR(data,p=2)
#data_m3 = VAR(data,p=3)
#data_m4 = VAR(data,p=4)

#data_m5 = VAR(data,p=5)

resi = data_m2$residuals
mq(resi, adj=4^2 * 2)

obj_spec <- bekk_spec(model = list(type = "bekk", asymmetric = T))

x1 <- bekk_fit(obj_spec, resi, QML_t_ratios = T, crit = 1e-09)

summary(x1)
#VIRF
N = 4
n = dim(x1$H_t)[1]
v <- 100
repl <- 20000
per = 0.01
index_series = 1
start = 1
end = n

sign_1 = matrix(rep(1, N), ncol = 1)
sign_2 = matrix(rep(0, N), ncol = 1)

C = x1$C0
A = x1$A
B = x1$B
G = x1$G
H_0 = x1$H_t

H_1 = matrix(NA,N,N*n)
for(i in 1:n){
  H_1[,(N*i-N+1):(N*i)] = matrix(H_0[i,],N,N)
}

CC = t(C)%*%C

shocks = matrix(c(rep(0,index_series-1),quantile(en[,index_series],per),rep(0,N-index_series)),nrow=1)

z1 = matrix(nrow=v,ncol=N*repl)
z2 = matrix(nrow=v,ncol=N*repl)

en = x1$e_t

z1[1,] = c(rep(0,(index_series-1)*repl),rep(quantile(en[,index_series],per),repl),rep(0,(N-index_series)*repl))
z1[2:v,1:repl] = matrix(sample(en[,1],size=(v-1)*repl,replace=T),nrow=v-1)
z1[2:v,(repl+1):(2*repl)]   = matrix(sample(en[,2],size=(v-1)*repl,replace=T),nrow=v-1)
z1[2:v,(2*repl+1):(3*repl)] = matrix(sample(en[,3],size=(v-1)*repl,replace=T),nrow=v-1)
z1[2:v,(3*repl+1):(4*repl)] = matrix(sample(en[,4],size=(v-1)*repl,replace=T),nrow=v-1)

z2[,1:repl] = matrix(sample(en[,1],size=v*repl,replace=T),nrow=v)
z2[,(repl+1):(2*repl)]   = matrix(sample(en[,2],size=v*repl,replace=T),nrow=v)
z2[,(2*repl+1):(3*repl)] = matrix(sample(en[,3],size=v*repl,replace=T),nrow=v)
z2[,(3*repl+1):(4*repl)] = matrix(sample(en[,4],size=v*repl,replace=T),nrow=v)
start_time <- Sys.time()
V = virf_bekka(start,end,H_1, A, B, G, CC, shocks, z1, z2, v , repl, N, sign_1,sign_2)
end_time <- Sys.time()
time_taken <- end_time - start_time
print(time_taken)

df <- data.frame(
  periods   = 1:v, 
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
