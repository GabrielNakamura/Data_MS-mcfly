
# Analysis ----------------------------------------------------------------

#Structural equation models (SEM) used to evaluate the statistical performance of ABC analysis
#implemented in the package "mcfly"

#Variable names
#Nspp_Metacommunity = Number of species in the metacommunity.
#W_Metacommunity = Real w slope of the metacommunity.
#Alpha_Metacommunity = Real OU's alpha of the metacommunity.
#Median_Posterior_Alpha = Median of the posterior distribution of OU's alpha.
#Prior_mode = Either 'uniform' or 'half-life'.
#log(Median_Posterior_W) - Log-transformed median of the posterior distribution of w slope.

#SEM  - Fig. 3
#Fig. 3A
sensitivity.mod.alpha<-psem(lm(Nspp_Metacommunity~log(W_Metacommunity+1)+Alpha_Metacommunity,data=res_full),lm(Median_Posterior_Alpha~Nspp_Metacommunity+Prior_mode+log(W_Metacommunity+1)+Alpha_Metacommunity,data=res_full),
                            lm(log(Median_Posterior_W)~Nspp_Metacommunity+log(W_Metacommunity+1)+Alpha_Metacommunity,data=res_full),Median_Posterior_Alpha%~~%log(Median_Posterior_W))
res_mod.alpha<-summary(sensitivity.mod.alpha)

#Fig. 3B
sensitivity.mod.hl<-psem(lm(Nspp_Metacommunity~log(W_Metacommunity+1)+Log.HL_Metacommunity,data=res_full),
                         lm(Log.HL~Nspp_Metacommunity+log(W_Metacommunity+1)+Prior_mode+Log.HL_Metacommunity,data=res_full),
                         lm(log(Median_Posterior_W)~Nspp_Metacommunity+log(W_Metacommunity+1)+Log.HL_Metacommunity,data=res_full),Log.HL%~~%log(Median_Posterior_W))
res_mod_hl<-summary(sensitivity.mod.hl)

# Fig. 5A
sensitivity.mod.alpha.unif<-psem(lm(Nspp_Metacommunity~log(W_Metacommunity+1)+Alpha_Metacommunity,data=unif),lm(Median_Posterior_Alpha~log(W_Metacommunity+1)+Alpha_Metacommunity,data=unif),
                                 lm(log(Median_Posterior_W)~Nspp_Metacommunity+log(W_Metacommunity+1)+Alpha_Metacommunity,data=unif),Median_Posterior_Alpha%~~%log(Median_Posterior_W))
res_mod.unif<-summary(sensitivity.mod.alpha.unif)

# Fig. 5B
sensitivity.mod.hl.unif<-psem(lm(Nspp_Metacommunity~log(W_Metacommunity+1)+Log.HL_Metacommunity,data=unif),
                              lm(Log.HL~log(W_Metacommunity+1)+Log.HL_Metacommunity,data=unif),
                              lm(log(Median_Posterior_W)~Nspp_Metacommunity+log(W_Metacommunity+1)+Log.HL_Metacommunity,data=unif),Log.HL%~~%log(Median_Posterior_W))
res_mod.hl.unif<-summary(sensitivity.mod.hl.unif)

# Fig. 5C
sensitivity.mod.alpha.hl<-psem(lm(Nspp_Metacommunity~log(W_Metacommunity+1)+Alpha_Metacommunity,data=hl),lm(Median_Posterior_Alpha~log(W_Metacommunity+1)+Alpha_Metacommunity,data=hl),
                               lm(log(Median_Posterior_W)~Nspp_Metacommunity+log(W_Metacommunity+1)+Alpha_Metacommunity,data=hl),Median_Posterior_Alpha%~~%log(Median_Posterior_W))
res_mod.alpha.hl<-summary(sensitivity.mod.alpha.hl)

# Fig. 5D
sensitivity.mod.hl.hl<-psem(lm(Nspp_Metacommunity~log(W_Metacommunity+1)+Log.HL_Metacommunity,data=hl),
                            lm(Log.HL~log(W_Metacommunity+1)+Log.HL_Metacommunity,data=hl),
                            lm(log(Median_Posterior_W)~Nspp_Metacommunity+log(W_Metacommunity+1)+Log.HL_Metacommunity,data=hl),Log.HL%~~%log(Median_Posterior_W))
res_mod.hl.hl<-summary(sensitivity.mod.hl.hl)