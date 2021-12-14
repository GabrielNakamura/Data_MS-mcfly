devtools::install_github("GabrielNakamura/mcfly")

library(mcfly)
library(robts)
library(here)

####Load general data####
comm.full<-as.matrix(read.table("comm_phyllostomidae_final.txt",h=T))
dim(comm.full)
single.spp<-which(rowSums(comm.full)<1)
comm<-comm.full[-single.spp,]
dim(comm)
data<-read.table("A_SAM_Grid_morrone.txt",h=T)
env<-data[-single.spp,c(10:15)]
esp<-data[-single.spp,c(3:4)]
dim(esp)
plot(esp[,1],esp[,2])

jet.colors<-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
colors2<-c("white",jet.colors(max(rowSums(comm))))
plot(esp, col=colors2[rowSums(comm)])

phy<-ape::read.tree("tree_phyllostomidae_final.tre")
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

tol<-define_tolerance(comm=comm,
                      occurrence=TRUE, entropy.order = 1,
                      phylo=phy,
                      envir=envir,
                      xy.coords=esp,
                      parallel=24,
                      m = 0.5,
                      OU.alpha = "uniform",
                      n.timestep = 50,
                      max.sample.size.prior = 1200,
                      sample.size.posterior = 1200,
                      scenario.ID = "doc",
                      output.dir.path = "einstein")


HDInterval::hdi(1-tol$Posterior,0.9)
hist(1-tol$Posterior[which(1-tol$Posterior<=0.48)],breaks=15)
length(1-tol$Posterior[which(1-tol$Posterior<=0.48)])
1-quantile(tol$Posterior,c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99))
tol$Tolerance

length(which(1-tol$Posterior<=0.4))
hist(tol$Posterior,breaks=30)
tol$Tolerance
max(tol$Posterior)

test.phyllostomidae<-mcfly(comm=comm,occurrence=TRUE,entropy.order=1,
                     phylo=phy,envir=envir,xy.coords=esp,
                     niche.breadth=10,m=0.5,n.timestep=50,OU.alpha="uniform",W.r.prior=TRUE,
                     summary.stat="correlation",tol=0.1,sample.size.posterior=480,
                     max.sample.size.prior=36000,HPD=0.9,
                     return.comm=TRUE,return.w.priors=TRUE,return.alpha.priors = TRUE,
                     parallel = 24,scenario.ID = "mcfly",output.dir.path = "delorean")

test.phyllostomidae$Alpha_Limits
plot(density(log(log(2)/test.phyllostomidae$Alpha_Posterior_Distribution)))

#compute density function
dens.a.prior<-density(x=test.phyllostomidae$Alpha_Prior_Distribution,from=test.phyllostomidae$Alpha_Limits[2],to=test.phyllostomidae$Alpha_Limits[3])
alpha.prior<-dens.a.prior$x
dens.alpha.prior<-dens.a.prior$y
plot(alpha.prior,dens.alpha.prior)

dens.a.post<-density(x=test.phyllostomidae$Alpha_Posterior_Distribution,from=test.phyllostomidae$Alpha_Limits[2],to=test.phyllostomidae$Alpha_Limits[3])
alpha.post<-dens.a.post$x
dens.alpha.post<-dens.a.post$y
plot(alpha.post,dens.alpha.post)
hdi.09.alpha.post<-HDInterval::hdi(dens.a.post,allowSplit=T,0.9)#uncorrected
plot(dens.a.post)

est.alpha<-alpha.post[which(dens.alpha.post==max(dens.alpha.post))]


hist(log(2)/test.phyllostomidae$Alpha_Posterior_Distribution)
dens.hl<-density(x=log(log(2)/test.phyllostomidae$Alpha_Posterior_Distribution),
            to=log(log(2)/test.phyllostomidae$Alpha_Limits[2]),
            from=log(log(2)/test.phyllostomidae$Alpha_Limits[3]))
plot(dens.hl)
hdi.09.hl<-HDInterval::hdi(dens.hl,allowSplit=TRUE,0.9) 
b_b<-log(2)/exp(hdi.09.hl[1,1])
b_a<-log(2)/exp(hdi.09.hl[1,2])
a_b<-log(2)/exp(hdi.09.hl[2,1])
a_a<-log(2)/exp(hdi.09.hl[2,2])
corr_hdi.09.alpha<-matrix(c(a_a,b_a,a_b,b_b),2,2,
          dimnames=list(c(1,2),c("begin","end")))#corrected


dens.prior.w<-density(x=test.phyllostomidae$W_Prior_Distribution,
                      from=min(test.phyllostomidae$W_Prior_Distribution),
                      to=max(test.phyllostomidae$W_Prior_Distribution))
prior.w<-dens.prior.w$x
dens.w.prior<-dens.prior.w$y

dens.post.w<-density(x=test.phyllostomidae$W_Posterior_Distribution,from=min(test.phyllostomidae$W_Prior_Distribution),
                     to=max(test.phyllostomidae$W_Prior_Distribution))
hist(test.phyllostomidae$W_Posterior_Distribution)
post.w<-dens.post.w$x
dens.w.post<-dens.post.w$y
max(dens.post.w$y)
plot(dens.post.w)

hdi.09.w<-HDInterval::hdi(dens.post.w,allowSplit=T,0.9) # I used only the first interval. All others seem meaningless.
max(test.phyllostomidae$W_Posterior_Distribution) #value used to define maximum x-value in the plot

int.alpha.w<-which(test.phyllostomidae$Alpha_Posterior_Distribution>=hdi.09.alpha.post[1,1]&test.phyllostomidae$Alpha_Posterior_Distribution<=hdi.09.alpha.post[1,2]
                   &test.phyllostomidae$Alpha_Posterior_Distribution>=hdi.09.w[1,1]&test.phyllostomidae$Alpha_Posterior_Distribution<=hdi.09.w[1,2])
post.alpha.int<-test.phyllostomidae$Alpha_Posterior_Distribution[int.alpha.w]#mean.sim.entropy computed for hpd.
length(post.alpha.int)

sim.div<-matrix(NA,nrow=nrow(comm),ncol=length(test.phyllostomidae$Alpha_Posterior_Distribution[int.alpha.w]))
cor.sim.div<-numeric()
sim.comm<-test.phyllostomidae$COMM.sim[unlist(lapply(test.phyllostomidae$COMM.sim, function(x) !is.null(x)))]
length(sim.comm)
sim.comm<-sim.comm[int.alpha.w]

for(i in 1:length(test.phyllostomidae$Alpha_Posterior_Distribution[int.alpha.w])){
  sim.div[,i]<-vegan::renyi(sim.comm[[i]],scales=1)
  cor.sim.div[i]<-cor(div,sim.div[,i])
}

mean.sim.div<-rowMeans(sim.div)
length(mean.sim.div)
plot(envir,mean.sim.div)
cor(mean.sim.div,div)
plot(mean.sim.div,div)

dat<-cbind(div,mean.sim.div,esp,envir)
write.table(dat,"dados.txt",sep=" ")

test.phyllostomidae$Sample_Attributes

alpha.post<-test.phyllostomidae$Alpha_Posterior_Distribution
hl.post<-log(2)/test.phyllostomidae$Alpha_Posterior_Distribution
log.hl.post<-log(hl.post)
alpha.prior<-test.phyllostomidae$Alpha_Prior_Distribution
hl.prior<-log(2)/test.phyllostomidae$Alpha_Prior_Distribution
log.hl.prior<-log(hl.prior)

dat.alpha.post<-cbind(alpha.post,hl.post,log.hl.post)
dat.alpha.prior<-cbind(alpha.prior,hl.prior,log.hl.prior)
write.table(dat.alpha.post,"dados_alpha_post.txt",sep=" ")
write.table(dat.alpha.prior,"dados_alpha_prior.txt",sep=" ")

data.dens.prior<-cbind(alpha.prior,dens.alpha.prior,prior.w,dens.w.prior)
write.table(data.dens.prior,"dados_dens_alpha_prior.txt",sep=" ")
data.dens.post<-cbind(alpha.post,dens.alpha.post,post.w,dens.w.post)
write.table(data.dens.post,"dados_dens_alpha_post.txt",sep=" ")

w.prior.raw<-test.phyllostomidae$W_Prior_Distribution
write.table(w.prior.raw,"dados_w_prior.txt",sep=" ")
w.post.raw<-test.phyllostomidae$W_Posterior_Distribution
write.table(w.post.raw,"dados_w_post.txt",sep=" ")


# definir parâmetros de dispersão

est.w<-post.w[which(dens.w.post==max(dens.w.post))]

colnames(esp)<-c("Lon","Lat")
dist.km<-as.dist(geodist::geodist(x=esp,measure = "geodesic"),diag=T,upper=T)/1000
dist.km.vec<-as.vector(as.dist(dist.km,diag=F,upper=F))
r<-as.vector(as.dist(scales::rescale(dist.km,c(0,1)),diag=F,upper=F))

dlk<-matrix(NA,length(r),3,dimnames=list(1:length(r),c("r","dlk","km")))

  for (p in 1:length(r)){
    dlk[p,1]<-r[p]
    dlk[p,2]<-round(exp(1)^-(est.w*r[p]^2),3)
    dlk[p,3]<-dist.km.vec[p]
  }

dlk.sub<-dlk[sample(nrow(dlk),1000,replace=F),]


plot(dlk.sub[,3],dlk.sub[,2])
write.table(dlk.sub,"w_curve.txt",sep=" ")



