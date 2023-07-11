# Function to perform ABC analysis using OSM. Only objects not described in "SimulatedData_mcfly.R" are described here.

library(ape)
library(geiger)
library(adephylo)
library(scales)
library(MCSim)
library(HDInterval)
library(mcfly)

RunTest_mcfly<-function(Nmeta=2,
                        Nspp.phy = 50,
                        Nspp.comm = 50,
                        Ncomm =50,
                        occurrence = TRUE,
                        subset=FALSE,
                        n.timestep = 50, #Number of timesteps of species assembly in metacommunity simulation.
                        niche.breadth = 10, #Niche breadths of species.
                        W.r=0,
                        W.r.prior = FALSE,
                        alpha.mode=c("uniform","half-life"),
                        alpha.interval=c(0.01,0.05),
                        sample.size.posterior = 24,
                        max.sample.size.prior = 2400, 
                        HPD = 0.9,#Set the credibility level of HPD.
                        parallel = 24,
                        scenario.ID="simul_mcfly",
                        output.dir.path = "meta_output"){
  sum.res<-matrix(NA,Nmeta,19, dimnames=list(1:Nmeta,c("Total_Prior_Size",
                                                    "Final_Posterior_Size",
                                                    "Min_Alpha_Prior",
                                                    "Max_Alpha_Prior",
                                                    "Tol",
                                                    "Nspp_Metacommunity",
                                                    "Alpha_Metacommunity",
                                                    "K_L-Posterior",
                                                    "K_H-Posterior",
                                                    "Mean_Posterior_Alpha",
                                                    "SD_Posterior_Alpha",
                                                    "Median_Posterior_Alpha",
                                                    "Alpha_L-Posterior",
                                                    "Alpha_H-Posterior",
                                                    "Mean_Posterior_W",
                                                    "SD_Posterior_W",
                                                    "Median_Posterior_W",
                                                    "W_L-Posterior",
                                                    "W_H-Posterior")))
    for (i in 1:Nmeta){
      output.dir.path.1<-paste(output.dir.path,i,sep = ".") #Internal argument MCSim.
      sim.ID<-paste(scenario.ID,i,sep=".") #Internal argument MCSim.
      phylo.meta <- geiger::rescale(geiger::sim.bdtree(b= 0.1, d= 0, stop= "taxa",
                                n = Nspp.phy, extinct = FALSE),
                                model="depth",depth=30) # Simulate phylogeny.
      DRoot<-suppressWarnings(as.numeric(adephylo::distRoot(phylo.meta,
                                                      method = "patristic")[1])) #Define timespan of the tree.
      min.alpha<-log(2)/DRoot #Define minimum OU's alpha value for the niche position of species in the OSM.
      max.alpha<-log(2)/(0.03333333*DRoot)#Define maximum OU's alpha value for the niche position of species in the OSM.
      prior.alpha<-runif(10*max.sample.size.prior,
                           min=min(alpha.interval),max=max(alpha.interval)) #Simulate OU's alpha interval.
      env.meta <- runif(Ncomm, min=0, max=100) #Simulate the environmental gradient for the sites in the OSM.
      xy.coords.meta <- data.frame(x= runif(length(env.meta), 1, 100), 
                                   y= runif(length(env.meta), 1, 100))#Simulate the spatial coordinates for the sites in the OSM.
      comm.names.meta <- rownames(as.matrix(xy.coords.meta), FALSE, prefix= "comm") # Define labels for local assemblages in the OSM.
      names(env.meta)<-comm.names.meta # Define labels for local assemblages described by environment.
      alpha.meta<-sample(prior.alpha,1) # Sample a OU's alpha value for the OSM.
      sigma.meta<-sqrt(sd(env.meta)) # Define the 'sigma' parameter of niche simulation.
      theta.meta<-mean(env.meta) # Define the "theta' parameter of niche simulation.
      root.value.meta<-mean(env.meta) #Define the 'root.value' parameter of niche simulation.
      JL.meta <- rep(1000,Ncomm) #Define the number of individuals in each assemblage of the OSM.
      JM.meta <- sum(JL.meta) # Define the total number of individuals in the OSM.
      my.landscape.meta <- MCSim::make.landscape(site.coords= xy.coords.meta,
                                        Ef= env.meta,
                                        m= 0.5, JM= JM.meta, JL= JL.meta)# Define 'my.landscape' argument of function "metasim" of the package 'MCsim'.
      low.rich<-TRUE
      while(low.rich){#Set the minimum number of species in the OSM as 20 species.
      thresh <- TRUE
      while(thresh){
        niche.pos.meta <- ape::rTraitCont(phy= phylo.meta, model= "OU", 
                             alpha= alpha.meta,
                             sigma= sigma.meta, theta= theta.meta, 
                             root.value= root.value.meta) #Simulate niche positions of species in the OSM.
        fitOU.comm <- suppressWarnings(geiger::fitContinuous(phylo.meta,
                         niche.pos.meta,
                         model="OU",bounds=list(alpha=c(min.alpha,max.alpha)),
                         ncores=parallel)) # Estimate the OU's alpha parameter for simulated niche positions.
        thresh <- fitOU.comm$opt$alpha<(alpha.meta - alpha.meta * 0.2)|
                      fitOU.comm$opt$alpha>(alpha.meta + alpha.meta*0.2)#Bound the estimated OU's alpha parameter of niche positions between +- 20% of 
                                                                        #the nominal values (object 'alpha.meta').
      }
      
       if(subset == TRUE){ # Extract niche positions of species under subset == TRUE.
        niche.pos.meta <- sample(niche.pos.meta, size= Nspp.comm)
        niche.pos.meta[order(as.numeric(gsub("s","",names(niche.pos.meta))))]
      } else {niche.pos.meta<-niche.pos.meta}
      
      spp.names.meta <- names(niche.pos.meta) #Set the names of species.
      spp.freq.meta <- rmultinom(1, size= JM.meta, 
                          prob = scales::rescale(rnorm(length(niche.pos.meta),
                          0,10),to=c(0.4,1))) #Define the frequency of species in the OSM.
      spp.gamma<- spp.freq.meta/JM.meta # Define the 'spp.gamma' argument of function "metasim" of the package 'MCSim".
      sim.meta<- MCSim::metasim(landscape = my.landscape.meta, nu = 0,
                       speciation.limit = 0, JM.limit = JM.meta, 
                       n.timestep = n.timestep, 
                       W.r = W.r, save.sim = FALSE, trait.Ef = niche.pos.meta, 
                       trait.Ef.sd = niche.breadth, gamma.abund = spp.gamma, 
                       taxa.list = spp.names.meta,scenario.ID=scenario.ID, 
                       sim.ID=sim.ID,
                       output.dir.path = output.dir.path.1)# Simute a OSM.
      comm.out.meta <- sim.meta$J.long # Extract the complete OSM assembly steps from 'metasim' output.
      comm.frame.meta <-comm.out.meta[which(comm.out.meta[,1]==n.timestep),]# Keep only the last OSM assembly step as dataframe.
      comm.obs.meta <- tapply(comm.frame.meta$count,list(comm.frame.meta$site,
                                                         comm.frame.meta$spp),
                                                                          sum)# Reshape OSM as matrix.
      rownames(comm.obs.meta) <- comm.names.meta #Set labels for the assemblages in the OSM.
      if(occurrence){
        comm.obs.meta <- ifelse(comm.obs.meta > 0, yes = 1, no = 0)
      }# Transform species abundances to presences/absences.
      low.rich<-length(which(colSums(comm.obs.meta)>0))<20
      } #Evaluate whether the low.rich criteria was achieved (see L. 68)
      sum.res[i,6]<-length(which(colSums(comm.obs.meta)>0))
 ####initiate tolerance level definition####
  output.dir.path.2<-paste("einstein",output.dir.path.1,sep=".")
  training<-mcfly::mcfly(comm = comm.obs.meta,
              phylo = phylo.meta,
              envir = env.meta,
              xy.coords = xy.coords.meta,
              m=0.5,
              W.r.prior= W.r.prior,
              OU.alpha = alpha.mode,
              tol = 1,
              max.sample.size.prior = 100*parallel,
              sample.size.posterior = 20*parallel,
              parallel = parallel,
              output.dir.path=output.dir.path.2) # Run 'mcfly' with maximum tolerance value (tol = 1).
  tol.i<-1-quantile(training$distance.statistics, probs = 0.99) # Define the 99th percentile as argument 'tol' in ABC analysis.
  output.dir.path.3<-paste("delorean",output.dir.path.1,sep=".") #Internal argument of the package "MCSim".
  test<-mcfly::mcfly(comm = comm.obs.meta, # A 'i in 1:Nmeta' OMS.
              phylo = phylo.meta, # Phylogeny.
              envir = env.meta, # Environment.
              xy.coords = xy.coords.meta, # Spatial coordinates.
              m=0.5, # Probability of immigration.
              W.r.prior= W.r.prior,
              tol = tol.i, # Tolerance value estimated previously.
              OU.alpha = alpha.mode,
              max.sample.size.prior = max.sample.size.prior,
              sample.size.posterior = sample.size.posterior,
              parallel = parallel,
              output.dir.path=output.dir.path.3) #Run ABC analysis.
#Hereafter the results are extracted from 'test' object .
  sum.res[i,1]<-test$Sample_Attributes[2,]
  sum.res[i,2]<-test$Sample_Attributes[4,]
  sum.res[i,3]<-min.alpha
  sum.res[i,4]<-max.alpha
  sum.res[i,5]<-tol.i
  sum.res[i,7]<-alpha.meta
    if (is.numeric(test$K_niche)){
    sum.res[i,8:9]<-HDInterval::hdi(object = test$K_niche, credMass = HPD,
                                    allowSplit=TRUE)
    } else {sum.res[i,8:9]<-NA}
  sum.res[i,10]<-mean(test$Alpha_Posterior_Distribution)
  suppressWarnings(if(!is.numeric(test$Alpha_Posterior_Distribution)){
    sum.res[i,10] <- NA
  })
  sum.res[i,11]<-sd(test$Alpha_Posterior_Distribution)
  suppressWarnings(if(!is.numeric(test$Alpha_Posterior_Distribution)){
    sum.res[i,11] <- NA
  })
  sum.res[i,12]<-quantile(test$Alpha_Posterior_Distribution, probs = 0.5)
  suppressWarnings(if(!is.numeric(test$Alpha_Posterior_Distribution)){
    sum.res[i,12] <- NA
  })
  sum.res[i,13:14]<-test$HPD_Alpha
  suppressWarnings(if(!is.numeric(test$HPD_Alpha)){
    sum.res[i,13:14] <- NA
  })
  if(W.r.prior==TRUE){
  sum.res[i,15]<-mean(test$W_Posterior_Distribution)
  suppressWarnings(if(!is.numeric(test$W_Posterior_Distribution)){
    sum.res[i,15] <- NA
  })
  } else{sum.res[i,15]<-NA
    }
  if(W.r.prior==TRUE){
  sum.res[i,16]<-sd(test$W_Posterior_Distribution)
  suppressWarnings(if(!is.numeric(test$W_Posterior_Distribution)){
    sum.res[i,16] <- NA
  })
  } else {sum.res[i,16]<-NA
  }
  if(W.r.prior==TRUE){
    sum.res[i,17]<-quantile(test$W_Posterior_Distribution, probs = 0.5)
    suppressWarnings(if(!is.numeric(test$W_Posterior_Distribution)){
      sum.res[i,17] <- NA
    })
  } else {sum.res[i,17]<-NA
  }
   if(W.r.prior==TRUE){
    sum.res[i,18:19]<-test$HPD_w
      if(!is.numeric(test$HPD_w)){
        sum.res[i,18:19] <- NA
    }
  } else{sum.res[i,18:19]<-NA
  }
  print(i)
}
  RES<-list(Summary_Results=sum.res)
  return(RES)
}




