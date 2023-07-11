library(ape)
library(here)
library(picante)
library(phangorn)

####Load general data and test with Phyllostomidae####
comm.full<- as.matrix(read.table(here::here("data", "comm_phyllostomidae_final.txt"),  h= T))
single.spp<-which(rowSums(comm.full)==0)
comm<-comm.full[-single.spp,]#community matrix.
grid.data<- read.table(here::here("data", "grid_data.txt"), h=T)
env<-grid.data[-single.spp,c(10:15)] #environmental variables.
spa<-as.matrix(grid.data[,3:4]) #spatial coordinates.
phy.raw<-ape::read.tree(here::here("data", "tree_phyllostomidae_final.tre")) #phylogenetic hypothesis for Phyllostomidae species.
phy<-phangorn::nnls.tree(cophenetic(phy.raw),phy.raw,rooted=TRUE)#make phylogeny ultrametric.
match<-picante::match.phylo.comm(phy,comm)
comm<-match$comm
dim(comm)
phy<-match$phy
rownames(comm)==rownames(env)

