# Install the required R packages
install.packages("vegan")
install.packages("compiler")

library(vegan)    # Contains the function for biodiversity metrics 
library(compiler) # Contains the functions to compile the R-commands to byte code to speed up processing 


iterations <- 5  # Number of iterations (manuscript used 100 iterations)
ts.length <- 100 # Time series length (manuscript used 100,000 time-steps)
rep.freq <- 10   # Reporting frequency (manuscript reported biodiversity metrics every 1000 time-steps.)

Spec <- 0.0025        # Speciation rate, v (manuscript used these levels: 0.00075, 0.0025, 0.005, 0.01)
Dispersal <- 0.75     # Dispersal rate, m  (manuscript used these levels: 0.25, 0.5, 0.75)
Disp <- (Dispersal/(1-Spec))*1


# Reporters for null counter-factual (sr = species richness; ev = species evenness; beta = community dissimilarity; a = local diversity (alpha), g = regional diversity (gamma))
null.sr.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
null.sr.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
null.ev.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
null.ev.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
null.beta <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)

#Reporters for development (dev) scenario (sr = species richness; ev = species evenness; beta = community dissimilarity; a = local diversity (alpha), g = regional diversity (gamma))
dev.sr.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
dev.sr.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
dev.ev.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
dev.ev.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
dev.beta <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)

# Reporters for restoration (res) scenrio (sr = species richness; ev = species evenness; beta = community dissimilarity; a = local diversity (alpha), g = regional diversity (gamma))
res.sr.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
res.sr.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
res.ev.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
res.ev.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
res.beta <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)

# Reporters for translocation (tran) scenario (sr = species richness; ev = species evenness; beta = community dissimilarity; a = local diversity (alpha), g = regional diversity (gamma))
tran.sr.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
tran.sr.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
tran.ev.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
tran.ev.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
tran.beta <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)

################################################################################
################################################################################


# Function for the forward neutral simulation (death-immigration function)
# Based on: Buschke, F.T. et al. (2016) Ecology and Evolution, doi:10.1002/ece3.2379

death.imm.fun <- function (comm.t,new.K) {
species <- 1:dim(comm.t)[1]
if (runif(1) < Spec) {
  pres.ab <- ifelse(colSums(comm.t)== 0,0,1)
  N.id <- sample(samples,1,prob=(new.K*pres.ab))
  ind.k <- sample(species,1,prob=comm.t[,N.id])
  new.r <- rep(0,N)
  new.r[N.id] <- 1
  comm.t <- rbind(comm.t,new.r)
  if (sum(comm.t[,N.id]) == new.K[N.id]) {
    comm.t[ind.k,N.id] <- comm.t[ind.k,N.id]-1
  }
} else{
  if (runif(1)>Disp) {
      pres.ab <- ifelse(colSums(comm.t)== 0,0,1)
      N.id <- sample(samples,1,prob=(new.K*pres.ab))
      ind.k <- sample(species,1,prob=comm.t[,N.id])
      if (sum(comm.t[,N.id]) == new.K[N.id]) {
        comm.t[ind.k,N.id] <- comm.t[ind.k,N.id]-1
      }
      ind.rep <- sample(species,1,prob=comm.t[,N.id])
      comm.t[ind.rep,N.id] <- comm.t[ind.rep,N.id]+1
  } else {
      pres.ab <- ifelse(colSums(comm.t)== 0,0,1)
      N.id <- sample(samples,1,prob=(new.K*pres.ab))
      ind.k <- sample(species,1,prob=comm.t[,N.id])
      if (sum(comm.t[,N.id]) == new.K[N.id]) {
        comm.t[ind.k,N.id] <- comm.t[ind.k,N.id]-1
      }
          
      imm.cell <- sample(samples,1,replace=TRUE, prob = (as.vector(imm.prob[N.id,])*pres.ab))
      
      ind.imm <- sample(species,1,prob=comm.t[,imm.cell])
      comm.t[ind.imm,N.id] <- comm.t[ind.imm,N.id]+1
  }
}
comm.t
}

death.imm <- cmpfun(death.imm.fun)

####################################################################################################
# The follwing code sets up the landscape

for (k in 1:iterations) {
  M <- 179        # Number of patches
  Sites <- 1:M
  K <- rep(100,M) # Habitat capacity of each patch


# Setting up the dispersal kernel based on inverse-quared distance  
dist.mat <- matrix(NA,nrow=M,ncol=M)
  for (i in 1:M) {
    for (j in 1:M) {
      if (abs(Sites[i] - Sites[j])<=89) {
        dist.mat[i,j] <- abs(Sites[i] - Sites[j])
      } else {
        dist.mat[i,j] <- 179 - (abs(Sites[i] - Sites[j]))
       }
     }
  }

disp.mat <- 1/(dist.mat^2)
diag(disp.mat) <- 0


###############################################################################

# Create a vector of the site identification for each species unit.
sites <- rep(1:M,K)

J <- sum(K) #  total number of species units

sp.name <- 1:J  # Start by assigning each species unit a unique ID
ind <- 1:J

# This is a probability vector that determines whether a species unit has already coalesed.
ind.prob <- rep(1,J)

# This is simply to set the counter to display the number of time-steps (starts at 1)
x <- 1

  
################################################################################
# Algorithm for the coalescence process to set up the simulation
# See: Rosindell, J., et al. (2008) Ecological Informatics, 3, 259-271.  
################################################################################

# Start the loop and continue until the number of species in the simulation matches the observed species richness
while(sum(ind.prob)> 0) {

# Choose a species unit at random
ind.id <- sample(ind,1,prob=ind.prob)

# First, Speciation
if (runif(1)<Spec) {
  sp.name[which(sp.name==ind.id)] <- sp.name[ind.id]
  # Since the focal species unit already has an assigned species ID, it is excluded from the model from now on.
  ind.prob[ind.id] <- 0
    } else {
  # Second, dispersal
    if (runif(1)>Disp) {
      # Identify the site location of selected species unit
      site.id <- sites[ind.id]
      # Identify a subset of species units that occur in the source patch
      migr.id <- ind[which(sites==site.id)]
        if(length(table(sp.name[migr.id]))>1) {
          # Select a subset of species units that DO NOT already share the same species ID
          imm.id <- migr.id[which(sp.name[migr.id]!=sp.name[ind.id])]
          # Select a parent/immigrant from the source quadrat
          sist.sp <- sample(imm.id,1)
          # Assign the same species ID as the parent species unit to the randomly selected focal species unit
          # AND all other species units that already share that species ID
          sp.name[which(sp.name==ind.id)] <- sp.name[sist.sp]
        }
      # Since the focal bird unit already has an assigned species ID, it is excluded from the model from now on.
      ind.prob[ind.id] <- 0
    } else {
  # Identify the site location of selected species unit
  site.id <- sites[ind.id]
  # Identify the source of selected species unit based on the dispersal kernel
  source.id <- sample(1:M,1,prob=(disp.mat[site.id,]*(K/100)))
  # Identify a subset of species units that occur in the source quadrat
  migr.id <- ind[which(sites==source.id)]
    if(length(table(sp.name[migr.id]))>1) {
      # Select a subset of species units that DO NOT already share the same species ID
      imm.id <- migr.id[which(sp.name[migr.id]!=sp.name[ind.id])]
      # Select a parent/immigrant from the source quadrat
      sist.sp <- sample(imm.id,1)
      # Assign the same species ID as the parent bird unit to the randomly selected focal bird unit
      # AND all other bird units that already share that species ID
       sp.name[which(sp.name==ind.id)] <- sp.name[sist.sp]
     }
    # Since the focal bird unit already has an assigned species ID, it is excluded from the model from now on.
   ind.prob[ind.id] <- 0
   }
  }
}

# Rearrange the new assemblage into a species-by site matrix
Patches <- as.matrix(table(sp.name,sites))

  
# Remove the diversity at the offset-receiving site (i.e. recent, but unrelated, degraded habitat)
Patches[,81:100] <- 0
K[81:100] <- 0

  
#####################################################################################
#####################################################################################
# Run the forward siimulation for scenario 1 (counterfactual)
  
Patches.null <- Patches
imm.prob <- disp.mat

N <- dim(Patches.null)[2] 
S <- dim(Patches.null)[1]

species <- seq(1:S)    
samples <- seq(1:N) 

for (i in 1:ts.length){
  Patches.null <- death.imm(Patches.null,K)
    if (i/rep.freq == floor(i/rep.freq))
      {
      # Summarise all the biodiversity metrics
        null.sr.a[k,i/rep.freq] <- mean(colSums(decostand(Patches.null,"pa")))
        null.sr.g[k,i/rep.freq] <- length(which(rowSums(decostand(Patches.null,"pa"))>0))
        null.ev.a[k,i/rep.freq] <-  mean((diversity(t(Patches.null)))/(log(colSums(decostand(Patches.null,"pa")))),na.rm=TRUE)
        null.ev.g[k,i/rep.freq] <-  mean((diversity(colSums(Patches.null)))/(log(length(which(rowSums(decostand(Patches.null,"pa"))>0)))),na.rm=TRUE)
        null.beta[k,i/rep.freq] <- mean(betadisper(vegdist(t(Patches.null[,which(colSums(Patches.null)>0)])), 
        as.factor(rep(1,length(which(colSums(Patches.null)>0)))), type = "centroid")$distances)
      }
  }


#####################################################################################
#####################################################################################
# Run the forward siimulation for scenario 2 (No offsets)

Patches.dev <- Patches
Patches.dev[,c(1:10,170:179)] <- 0
K.dev <- K
K.dev[c(1:10,170:179)] <- 0 # Clear the development site

for (i in 1:ts.length){
  Patches.dev <- death.imm(Patches.dev,K.dev)
    if (i/rep.freq == floor(i/rep.freq))
      {
      # Summarise all the biodiversity metrics
        dev.sr.a[k,i/rep.freq] <- mean(colSums(decostand(Patches.dev,"pa")))
        dev.sr.g[k,i/rep.freq] <- length(which(rowSums(decostand(Patches.dev,"pa"))>0))
        dev.ev.a[k,i/rep.freq] <-  mean((diversity(t(Patches.dev)))/(log(colSums(decostand(Patches.dev,"pa")))),na.rm=TRUE)
        dev.ev.g[k,i/rep.freq] <-  mean((diversity(colSums(Patches.dev)))/(log(length(which(rowSums(decostand(Patches.dev,"pa"))>0)))),na.rm=TRUE)
        dev.beta[k,i/rep.freq] <- mean(betadisper(vegdist(t(Patches.dev[,which(colSums(Patches.dev)>0)])), 
        as.factor(rep(1,length(which(colSums(Patches.dev)>0)))), type = "centroid")$distances)
       }
  }


#####################################################################################
#####################################################################################
# Run the forward siimulation for scenario 3 (passive restoration)

Patches.res <- Patches
Patches.res[,c(1:10,170:179)] <- 0
K.res <- K
K.res[c(1:10,170:179)] <- 0   # Clear the development site
K.res[81:100] <- 100          # Increa the habitat capacity at the offset site

for (i in 1:ts.length){
  Patches.res <- death.imm(Patches.res,K.res)
    if (i/rep.freq == floor(i/rep.freq))
      {
      # Summarise all the biodiversity metrics
        res.sr.a[k,i/rep.freq] <- mean(colSums(decostand(Patches.res,"pa")))
        res.sr.g[k,i/rep.freq] <- length(which(rowSums(decostand(Patches.res,"pa"))>0))
        res.ev.a[k,i/rep.freq] <-  mean((diversity(t(Patches.res)))/(log(colSums(decostand(Patches.res,"pa")))),na.rm=TRUE)
        res.ev.g[k,i/rep.freq] <-  mean((diversity(colSums(Patches.res)))/(log(length(which(rowSums(decostand(Patches.res,"pa"))>0)))),na.rm=TRUE)
        res.beta[k,i/rep.freq] <- mean(betadisper(vegdist(t(Patches.res[,which(colSums(Patches.res)>0)])), 
        as.factor(rep(1,length(which(colSums(Patches.res)>0)))), type = "centroid")$distances)
      }
   }

#####################################################################################
#####################################################################################
# Run the forward siimulation for scenario 4 (Translocation)

Patches.tran <- Patches
Patches.tran[,81:100] <- Patches.tran[,c(1:10,170:179)] # Translocate the communities
Patches.tran[,c(1:10,170:179)] <- 0  # Clear the development site
K.tran <- K
K.tran[c(1:10,170:179)] <- 0
K.tran[81:100] <- 100

for (i in 1:ts.length){
  Patches.tran <- death.imm(Patches.tran,K.tran)
    if (i/rep.freq == floor(i/rep.freq))
      {
      # Summarise all the biodiversity metrics
        tran.sr.a[k,i/rep.freq] <- mean(colSums(decostand(Patches.tran,"pa")))
        tran.sr.g[k,i/rep.freq] <- length(which(rowSums(decostand(Patches.tran,"pa"))>0))
        tran.ev.a[k,i/rep.freq] <-  mean((diversity(t(Patches.tran)))/(log(colSums(decostand(Patches.tran,"pa")))),na.rm=TRUE)
        tran.ev.g[k,i/rep.freq] <-  mean((diversity(colSums(Patches.tran)))/(log(length(which(rowSums(decostand(Patches.tran,"pa"))>0)))),na.rm=TRUE)
        tran.beta[k,i/rep.freq] <- mean(betadisper(vegdist(t(Patches.tran[,which(colSums(Patches.tran)>0)])), 
        as.factor(rep(1,length(which(colSums(Patches.tran)>0)))), type = "centroid")$distances)
      }
  }

  print(k)
flush.console()
}



#################################################################################
# Average the biodiversity metric over the 100 iterations
#################################################################################

m.null.sr.a <- apply(null.sr.a,2,function(x){mean(x,na.rm=TRUE)})
m.null.sr.g <- apply(null.sr.g,2,function(x){mean(x,na.rm=TRUE)})
m.null.ev.a <- apply(null.ev.a,2,function(x){mean(x,na.rm=TRUE)})
m.null.ev.g <- apply(null.ev.g,2,function(x){mean(x,na.rm=TRUE)})
m.null.beta <- apply(null.beta,2,function(x){mean(x,na.rm=TRUE)})

m.dev.sr.a <- apply(dev.sr.a,2,function(x){mean(x,na.rm=TRUE)})
m.dev.sr.g <- apply(dev.sr.g,2,function(x){mean(x,na.rm=TRUE)})
m.dev.ev.a <- apply(dev.ev.a,2,function(x){mean(x,na.rm=TRUE)})
m.dev.ev.g <- apply(dev.ev.g,2,function(x){mean(x,na.rm=TRUE)})
m.dev.beta <- apply(dev.beta,2,function(x){mean(x,na.rm=TRUE)})

m.res.sr.a <- apply(res.sr.a,2,function(x){mean(x,na.rm=TRUE)})
m.res.sr.g <- apply(res.sr.g,2,function(x){mean(x,na.rm=TRUE)})
m.res.ev.a <- apply(res.ev.a,2,function(x){mean(x,na.rm=TRUE)})
m.res.ev.g <- apply(res.ev.g,2,function(x){mean(x,na.rm=TRUE)})
m.res.beta <- apply(res.beta,2,function(x){mean(x,na.rm=TRUE)})

m.tran.sr.a <- apply(tran.sr.a,2,function(x){mean(x,na.rm=TRUE)})
m.tran.sr.g <- apply(tran.sr.g,2,function(x){mean(x,na.rm=TRUE)})
m.tran.ev.a <- apply(tran.ev.a,2,function(x){mean(x,na.rm=TRUE)})
m.tran.ev.g <- apply(tran.ev.g,2,function(x){mean(x,na.rm=TRUE)})
m.tran.beta <- apply(tran.beta,2,function(x){mean(x,na.rm=TRUE)})


#################################################################################
# This plots the outputs from the simulation
# Each plot shows how the log response ratio between the offset scenario and the counterfactual varies through time.

par(mfrow=c(2,3))  # Set up a multi-panel plot with 2 rows and three columns.

# Top-left panel is for local species richness
yval <- seq(rep.freq,ts.length,by=rep.freq)
plot (log(m.dev.sr.a/m.null.sr.a)~yval,type="l",col="red",ylab="Alpha richness",ylim=c(-.2,.2), xlab="Time steps", main="Local species richness")
lines(log(m.res.sr.a/m.null.sr.a)~yval,col="blue")
lines(log(m.tran.sr.a/m.null.sr.a)~yval,col="green")
abline(h=0,lty=2)

# The top-middle panel is for regional species richness
plot (log(m.dev.sr.g/m.null.sr.g)~yval,type="l",col="red",ylab="Gamma richness",ylim=c(-.2,.2), xlab="Time steps", main="Regional species richness")
lines(log(m.res.sr.g/m.null.sr.g)~yval,col="blue")
lines(log(m.tran.sr.g/m.null.sr.g)~yval,col="green")
abline(h=0,lty=2)

# The top-right panel is for local species evenness
plot (log(m.dev.ev.a/m.null.ev.a)~yval,type="l",col="red",ylab="Alpha evenness",ylim=c(-.5,.5), xlab="Time steps", main="Local species evenness")
lines(log(m.res.ev.a/m.null.ev.a)~yval,col="blue")
lines(log(m.tran.ev.g/m.null.ev.g)~yval,col="green")
abline(h=0,lty=2)

# The bottom-left panel is for regional species evenness
plot (log(m.dev.ev.g/m.null.ev.g)~yval,type="l",col="red",ylab="Gamma evennes",ylim=c(-.2,.2), xlab="Time steps", main="Regional species evenness")
lines(log(m.res.ev.g/m.null.ev.g)~yval,col="blue")
lines(log(m.tran.ev.g/m.null.ev.g)~yval,col="green")
abline(h=0,lty=2)

# The bottom-middle panel is for community dissimilarity
plot (log(m.dev.beta/m.null.beta)~yval,type="l",col="red",ylab="Beta dispersion",ylim=c(-.05,.05), xlab="Time steps", main="Community dissimilarity")
lines(log(m.res.beta/m.null.beta)~yval,col="blue")
lines(log(m.tran.beta/m.null.beta)~yval,col="green")
abline(h=0,lty=2)

# This plots a makeshift legen in the bottom right panel
plot(0,0,type="n",xlim=c(0,10),ylim=c(0,10),axes=F, xlab="",ylab="")
lines (c(1,3),c(8,8),lty=2) 
text(3.5,8,"No Net loss",pos=4)
lines (c(1,3),c(6,6),col="red") 
text(3.5,6,"Scenario 2 (No offset)",pos=4)
lines (c(1,3),c(4,4),col="blue") 
text(3.5,4,"Scenario 3 (Passive restoration)",pos=4)
lines (c(1,3),c(2,2),col="green") 
text(3.5,2,"Scenario 4 (Translocation)",pos=4)

########################################################################################
########################################################################################

# Write all the outputs to file (will be saved in your default working directory)
write.table(null.sr.a,"Null_alpha.txt",quote=F,row.names=F,sep="\t")
write.table(null.sr.g,"Null_gamma.txt",quote=F,row.names=F,sep="\t")
write.table(null.ev.a,"Null_even_a.txt",quote=F,row.names=F,sep="\t")
write.table(null.ev.g,"Null_even_g.txt",quote=F,row.names=F,sep="\t")
write.table(null.beta,"Null_beta.txt",quote=F,row.names=F,sep="\t")

write.table(dev.sr.a,"Dev_alpha.txt",quote=F,row.names=F,sep="\t")
write.table(dev.sr.g,"Dev_gamma.txt",quote=F,row.names=F,sep="\t")
write.table(dev.ev.a,"Dev_even_a.txt",quote=F,row.names=F,sep="\t")
write.table(dev.ev.g,"Dev_even_g.txt",quote=F,row.names=F,sep="\t")
write.table(dev.beta,"Dev_beta.txt",quote=F,row.names=F,sep="\t")

write.table(res.sr.a,"Res_alpha.txt",quote=F,row.names=F,sep="\t")
write.table(res.sr.g,"Res_gamma.txt",quote=F,row.names=F,sep="\t")
write.table(res.ev.a,"Res_even_a.txt",quote=F,row.names=F,sep="\t")
write.table(res.ev.g,"Res_even_g.txt",quote=F,row.names=F,sep="\t")
write.table(res.beta,"Res_beta.txt",quote=F,row.names=F,sep="\t")

write.table(tran.sr.a,"Tran_alpha.txt",quote=F,row.names=F,sep="\t")
write.table(tran.sr.g,"Tran_gamma.txt",quote=F,row.names=F,sep="\t")
write.table(tran.ev.a,"Tran_even_a.txt",quote=F,row.names=F,sep="\t")
write.table(tran.ev.g,"Tran_even_g.txt",quote=F,row.names=F,sep="\t")
write.table(tran.beta,"Tran_beta.txt",quote=F,row.names=F,sep="\t")
