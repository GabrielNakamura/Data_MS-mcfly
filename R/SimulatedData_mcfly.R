# Script for performing ABC analysis using simulated datasets. In these analyses,
# observed metacomunities are simulated using the package MCSim, and are hereafter called OSM.
# The function "RunSimTest_mcfly.R" simulated the OSM, and then perform ABC analysis using
# the package 'mcfly'.

source("RunSimTest_mcfly.R")
library(mcfly)

W0<-RunTest_mcfly(Nmeta=200,#Number of OSM.
                            Nspp.phy = 50,#Number of species in the phylogeny.
                            Nspp.comm = 50,#Number of species in each OSM.
                            Ncomm =50,#Number of assemblages in each OSM.
                            occurrence = TRUE,#TRUE if OSM is a binary matrix. FALSE if species abundances are available.
                            subset=FALSE,#TRUE if species in OSM are a subset of those in the phylogeny. FALSE otherwise.
                            W.r=0,#w slope parameter of the OSM.
                            W.r.prior = TRUE,#TRUE if w slope parameter is to be estimated.
                            alpha.mode="uniform",# OU's alpha parameter prior mode. Default option ("uniform") must be used.
                            alpha.interval=c(0.02310491,0.6931472),#OU's alpha intervals defined for the prior.
                            sample.size.posterior = 240,#Maximum sample size for the posterior distribution.
                            max.sample.size.prior = 14400,#Maximum sample size of the prior distribution.
                            parallel = 24)#Number of parallel processing cores.

W1<-RunTest_mcfly(Nmeta=200,
                            Nspp.phy = 50,
                            Nspp.comm = 50,
                            Ncomm =50,
                            occurrence = TRUE,
                            subset=FALSE,
                            W.r=1,
                            W.r.prior = TRUE,
                            alpha.mode="uniform",
                            alpha.interval=c(0.02310491,0.6931472),
                            sample.size.posterior = 240,
                            max.sample.size.prior = 14400,
                            parallel = 24)

W5<-RunTest_mcfly(Nmeta=200,
                                 Nspp.phy = 50,
                                 Nspp.comm = 50,
                                 Ncomm =50,
                                 occurrence = TRUE,
                                 subset=FALSE,
                                 W.r=5,
                                 W.r.prior = TRUE,
                                 alpha.mode="uniform",
                                 alpha.interval=c(0.02310491,0.6931472),
                                 sample.size.posterior = 240,
                                 max.sample.size.prior = 14400,
                                 parallel = 24)
W15<-RunTest_mcfly(Nmeta=200,
                                 Nspp.phy = 50,
                                 Nspp.comm = 50,
                                 Ncomm =50,
                                 occurrence = TRUE,
                                 subset=FALSE,
                                 W.r=15,
                                 W.r.prior = TRUE,
                                 alpha.mode="uniform",
                                 alpha.interval=c(0.02310491,0.6931472),
                                 sample.size.posterior = 240,
                                 max.sample.size.prior = 14400,
                                 parallel = 24)
W100<-RunTest_mcfly(Nmeta=200,
                                 Nspp.phy = 50,
                                 Nspp.comm = 50,
                                 Ncomm =50,
                                 occurrence = TRUE,
                                 subset=FALSE,
                                 W.r=100,
                                 W.r.prior = TRUE,
                                 alpha.mode="uniform",
                                 alpha.interval=c(0.02310491,0.6931472),
                                 sample.size.posterior = 240,
                                 max.sample.size.prior = 14400,
                                 parallel = 24)
