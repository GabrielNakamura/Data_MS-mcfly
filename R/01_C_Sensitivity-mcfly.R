
# Read packages and data for SEM ------------------------------------------

library(piecewiseSEM)

#Variable names
#Nspp_Metacommunity = Number of species in the metacommunity.
#W_Metacommunity = Real w slope of the metacommunity.
#Alpha_Metacommunity = Real OU's alpha of the metacommunity.
#Median_Posterior_Alpha = Median of the posterior distribution of OU's alpha.
#Prior_mode = Either 'uniform' or 'half-life'.
#log(Median_Posterior_W) - Log-transformed median of the posterior distribution of w slope.

res_full<-read.table(here::here("Data", "test_final.txt"),h=T)
unif<-as.data.frame(res_full[1001:2000,])
hl<-as.data.frame(res_full[1:1000,])

