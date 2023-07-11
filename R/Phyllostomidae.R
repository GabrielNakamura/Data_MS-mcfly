devtools::install_github("sokole/MCSim")
devtools::install_github("GabrielNakamura/mcfly")
library(vegan)
library(ape)
library(MCSim)
library(picante)
library(scales)
library(parallel)
library(HDInterval)
library(mcfly)

test.phyllostomidae<-mcfly(comm=comm[1:50,1:50],occurrence=TRUE,entropy.order=1,
                     phylo=phy,envir=env[1:50,4],xy.coords=spa[1:50,],
                     niche.breadth=10,m=0.5,n.timestep=50,OU.alpha="uniform",W.r.prior=TRUE,
                     tol=1,sample.size.posterior=6,
                     max.sample.size.prior=36,HPD=0.9,
                     return.comm=TRUE,
                     parallel = 6,scenario.ID = "mcfly",output.dir.path = "delorean")


