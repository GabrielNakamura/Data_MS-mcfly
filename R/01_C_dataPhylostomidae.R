
# Read data and libraries -------------------------------------------------

library(mcfly)
library(robts)
library(here)

comm.full<-as.matrix(read.table(here::here("Data", "comm_phyllostomidae_final.txt"), h = T))
dim(comm.full)
single.spp<-which(rowSums(comm.full)<1)
comm<-comm.full[-single.spp,]
dim(comm)
data<-read.table(here::here("Data", "A_SAM_Grid_morrone.txt"), h=T)
env<-data[-single.spp,c(10:15)]
esp<-data[-single.spp,c(3:4)]
dim(esp)
plot(esp[,1],esp[,2])

jet.colors<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
colors2<-c("white",jet.colors(max(rowSums(comm))))
plot(esp, col=colors2[rowSums(comm)])

phy<-ape::read.tree(here::here("Data", "tree_phyllostomidae_final.tre"))
max(cophenetic(phy))/2
phy$tip.label
ape::is.ultrametric(phy)
phy.ultra<-phangorn::nnls.tree(cophenetic(phy),phy,rooted=TRUE)
ape::is.ultrametric(phy.ultra)
match<-picante::match.phylo.comm(phy.ultra,comm)
comm<-match$comm
dim(comm)
phy<-match$phy
rownames(comm)==rownames(env)

plot(phy.ultra)

div<-vegan::renyi(comm,scales=1)
hist(div)
envir<-scales::rescale(env$Temp_Seas,c(1,100))
plot(envir,div)
cor(envir,div)