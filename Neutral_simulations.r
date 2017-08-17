install.packages("vegan")
install.packages("compiler")

library(vegan)
library(compiler)


iterations <- 100  # Number of iterations
ts.length <- 100000 # Time series length
rep.freq <- 1000   # Reporting frequency

Spec <- 0.0025 # Speciation rate     (0.00075, 0.0025, 0.005)
Dispersal <- 0.75 # Dispersal rate    (0.25, 0.5, 0.75)
Disp <- (Dispersal/(1-Spec))*1




# Reporters for null counter-factual
null.ts.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
null.ts.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
null.ev.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
null.ev.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
null.beta <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)

#Reporters for development scenario
dev.ts.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
dev.ts.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
dev.ev.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
dev.ev.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
dev.beta <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)

# Reporters for rehabilitation scenrio
res.ts.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
res.ts.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
res.ev.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
res.ev.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
res.beta <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)

# Reporters for translocation scenario
tran.ts.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
tran.ts.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
tran.ev.a <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
tran.ev.g <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)
tran.beta <- matrix(NA,ncol=(ts.length/rep.freq),nrow=iterations)

################################################################################

# Define the function for the forward neutral simulation


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

# Compile the function to byte code to make the algorithm run more efficiently

death.imm <- cmpfun(death.imm.fun)

####################################################################################################

# Set up the simulation landscape

for (k in 1:iterations) {

M <- 179 # number of sites


Sites <- 1:M

K <- rep(100,M)

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

# Initialising the model using coalescence


# Create a vector of the site identification for each organism unit.
sites <- rep(1:M,K)

J <- sum(K) #  total number of organism units

sp.name <- 1:J  # Start by assigning each organism unit a unique ID

ind <- 1:J

# This is a probability vector that determines whether a organism unit has already coalesed.
ind.prob <- rep(1,J)

# This is simply to set the counter to display the number of time-steps (starts at 1)
x <- 1


# Start the loop and continue until the number of species in the simulation matches the observed species richness
while(sum(ind.prob)> 0) {

# Choose a organism unit at random
ind.id <- sample(ind,1,prob=ind.prob)

# First, Speciation
if (runif(1)<Spec) {
sp.name[which(sp.name==ind.id)] <- sp.name[ind.id]
# Since the focal organism unit already has an assigned species ID, it is excluded from the model from now on.
ind.prob[ind.id] <- 0
} else {
# Second, dispersal
if (runif(1)>Disp) {
  # Identify the site location of selected bird unit
  site.id <- sites[ind.id]
  # Identify a subset of organism units that occur in the source quadrat
  migr.id <- ind[which(sites==site.id)]
  
  
if(length(table(sp.name[migr.id]))>1) {
  # Select a subset of organism units that DO NOT already share the same species ID
  imm.id <- migr.id[which(sp.name[migr.id]!=sp.name[ind.id])]

  # Select a parent/immigrant from the source quadrat
  sist.sp <- sample(imm.id,1)

  # Assign the same species ID as the parent organism unit to the randomly selected focal organism unit
  # AND all other organism units that already share that species ID
  sp.name[which(sp.name==ind.id)] <- sp.name[sist.sp]
}
  # Since the focal organism unit already has an assigned species ID, it is excluded from the model from now on.
  ind.prob[ind.id] <- 0
} else {
# Identify the site location of selected organism unit
site.id <- sites[ind.id]
# Identify the source of selected organism unit based on the dispersal kernel
source.id <- sample(1:M,1,prob=(disp.mat[site.id,]*(K/100)))

# Identify a subset of organism units that occur in the source quadrat
migr.id <- ind[which(sites==source.id)]

if(length(table(sp.name[migr.id]))>1) {

# Select a subset of organism units that DO NOT already share the same species ID
imm.id <- migr.id[which(sp.name[migr.id]!=sp.name[ind.id])]

# Select a parent/immigrant from the source quadrat
sist.sp <- sample(imm.id,1)

# Assign the same species ID as the parent organism unit to the randomly selected focal organism unit
# AND all other organism units that already share that species ID
sp.name[which(sp.name==ind.id)] <- sp.name[sist.sp]
}

# Since the focal organism unit already has an assigned species ID, it is excluded from the model from now on.
ind.prob[ind.id] <- 0
}
}
}


# This summarises the information as a site by species matrix and clear the offset receiving site of organisms

Patches <- as.matrix(table(sp.name,sites))
Patches[,81:100] <- 0
K[81:100] <- 0



#####################################################################################

# Set up the inputs for the simulation

Patches.null <- Patches
imm.prob <- disp.mat

N <- dim(Patches.null)[2] 
S <- dim(Patches.null)[1]

species <- seq(1:S)    # just a vector of species
samples <- seq(1:N) 



################################################################################
################################################################################

# Run the forward simulation for the counterfactual scenario (Scenario 1))

for (i in 1:ts.length){
Patches.null <- death.imm(Patches.null,K)
if (i/rep.freq == floor(i/rep.freq))
{
null.ts.a[k,i/rep.freq] <- mean(colSums(decostand(Patches.null,"pa")))
null.ts.g[k,i/rep.freq] <- length(which(rowSums(decostand(Patches.null,"pa"))>0))
null.ev.a[k,i/rep.freq] <-  mean((diversity(t(Patches.null)))/(log(colSums(decostand(Patches.null,"pa")))),na.rm=TRUE)
null.ev.g[k,i/rep.freq] <-  mean((diversity(colSums(Patches.null)))/(log(length(which(rowSums(decostand(Patches.null,"pa"))>0)))),na.rm=TRUE)
null.beta[k,i/rep.freq] <- mean(betadisper(vegdist(t(Patches.null[,which(colSums(Patches.null)>0)])), 
as.factor(rep(1,length(which(colSums(Patches.null)>0)))), type = "centroid")$distances)
}
}

################################################################################
################################################################################

# Run the forward simulation for the development without offsets scenario (Scenario 2))

Patches.dev <- Patches
Patches.dev[,c(1:10,170:179)] <- 0
K.dev <- K
K.dev[c(1:10,170:179)] <- 0

for (i in 1:ts.length){
Patches.dev <- death.imm(Patches.dev,K.dev)
if (i/rep.freq == floor(i/rep.freq))
{
dev.ts.a[k,i/rep.freq] <- mean(colSums(decostand(Patches.dev,"pa")))
dev.ts.g[k,i/rep.freq] <- length(which(rowSums(decostand(Patches.dev,"pa"))>0))
dev.ev.a[k,i/rep.freq] <-  mean((diversity(t(Patches.dev)))/(log(colSums(decostand(Patches.dev,"pa")))),na.rm=TRUE)
dev.ev.g[k,i/rep.freq] <-  mean((diversity(colSums(Patches.dev)))/(log(length(which(rowSums(decostand(Patches.dev,"pa"))>0)))),na.rm=TRUE)
dev.beta[k,i/rep.freq] <- mean(betadisper(vegdist(t(Patches.dev[,which(colSums(Patches.dev)>0)])), 
as.factor(rep(1,length(which(colSums(Patches.dev)>0)))), type = "centroid")$distances)
}
}


################################################################################
################################################################################

# Run the forward simulation for the passive restoration scenario (Scenario 3))

Patches.res <- Patches
Patches.res[,c(1:10,170:179)] <- 0
K.res <- K
K.res[c(1:10,170:179)] <- 0
K.res[81:100] <- 100

for (i in 1:ts.length){
Patches.res <- death.imm(Patches.res,K.res)
if (i/rep.freq == floor(i/rep.freq))
{
res.ts.a[k,i/rep.freq] <- mean(colSums(decostand(Patches.res,"pa")))
res.ts.g[k,i/rep.freq] <- length(which(rowSums(decostand(Patches.res,"pa"))>0))
res.ev.a[k,i/rep.freq] <-  mean((diversity(t(Patches.res)))/(log(colSums(decostand(Patches.res,"pa")))),na.rm=TRUE)
res.ev.g[k,i/rep.freq] <-  mean((diversity(colSums(Patches.res)))/(log(length(which(rowSums(decostand(Patches.res,"pa"))>0)))),na.rm=TRUE)
res.beta[k,i/rep.freq] <- mean(betadisper(vegdist(t(Patches.res[,which(colSums(Patches.res)>0)])), 
as.factor(rep(1,length(which(colSums(Patches.res)>0)))), type = "centroid")$distances)
}
}

################################################################################
################################################################################

# Run the forward simulation for the translocation scenario (Scenario 4))

Patches.tran <- Patches
Patches.tran[,81:100] <- Patches.tran[,c(1:10,170:179)]
Patches.tran[,c(1:10,170:179)] <- 0
K.tran <- K
K.tran[c(1:10,170:179)] <- 0
K.tran[81:100] <- 100

for (i in 1:ts.length){
Patches.tran <- death.imm(Patches.tran,K.tran)
if (i/rep.freq == floor(i/rep.freq))
{
tran.ts.a[k,i/rep.freq] <- mean(colSums(decostand(Patches.tran,"pa")))
tran.ts.g[k,i/rep.freq] <- length(which(rowSums(decostand(Patches.tran,"pa"))>0))
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

# Avergare simulation outputs across the 100 iterations

m.null.ts.a <- apply(null.ts.a,2,function(x){mean(x,na.rm=TRUE)})
m.null.ts.g <- apply(null.ts.g,2,function(x){mean(x,na.rm=TRUE)})
m.null.ev.a <- apply(null.ev.a,2,function(x){mean(x,na.rm=TRUE)})
m.null.ev.g <- apply(null.ev.g,2,function(x){mean(x,na.rm=TRUE)})
m.null.beta <- apply(null.beta,2,function(x){mean(x,na.rm=TRUE)})

m.dev.ts.a <- apply(dev.ts.a,2,function(x){mean(x,na.rm=TRUE)})
m.dev.ts.g <- apply(dev.ts.g,2,function(x){mean(x,na.rm=TRUE)})
m.dev.ev.a <- apply(dev.ev.a,2,function(x){mean(x,na.rm=TRUE)})
m.dev.ev.g <- apply(dev.ev.g,2,function(x){mean(x,na.rm=TRUE)})
m.dev.beta <- apply(dev.beta,2,function(x){mean(x,na.rm=TRUE)})

m.res.ts.a <- apply(res.ts.a,2,function(x){mean(x,na.rm=TRUE)})
m.res.ts.g <- apply(res.ts.g,2,function(x){mean(x,na.rm=TRUE)})
m.res.ev.a <- apply(res.ev.a,2,function(x){mean(x,na.rm=TRUE)})
m.res.ev.g <- apply(res.ev.g,2,function(x){mean(x,na.rm=TRUE)})
m.res.beta <- apply(res.beta,2,function(x){mean(x,na.rm=TRUE)})

m.tran.ts.a <- apply(tran.ts.a,2,function(x){mean(x,na.rm=TRUE)})
m.tran.ts.g <- apply(tran.ts.g,2,function(x){mean(x,na.rm=TRUE)})
m.tran.ev.a <- apply(tran.ev.a,2,function(x){mean(x,na.rm=TRUE)})
m.tran.ev.g <- apply(tran.ev.g,2,function(x){mean(x,na.rm=TRUE)})
m.tran.beta <- apply(tran.beta,2,function(x){mean(x,na.rm=TRUE)})


#################################################################################

# Make plots for all the simulation outputs

par(mfrow=c(2,3))
yval <- seq(rep.freq,ts.length,by=rep.freq)
plot (log(m.dev.ts.a/m.null.ts.a)~yval,type="l",col="red",ylab="Alpha richness",ylim=c(-.2,.2), xlab="Time steps")
lines(log(m.res.ts.a/m.null.ts.a)~yval,col="blue")
lines(log(m.tran.ts.a/m.null.ts.a)~yval,col="green")
abline(h=0,lty=2)

plot (log(m.dev.ts.g/m.null.ts.g)~yval,type="l",col="red",ylab="Gamma richness",ylim=c(-.2,.2), xlab="Time steps")
lines(log(m.res.ts.g/m.null.ts.g)~yval,col="blue")
lines(log(m.tran.ts.g/m.null.ts.g)~yval,col="green")
abline(h=0,lty=2)

plot (log(m.dev.ev.a/m.null.ev.a)~yval,type="l",col="red",ylab="Alpha evenness",ylim=c(-.5,.5), xlab="Time steps")
lines(log(m.res.ev.a/m.null.ev.a)~yval,col="blue")
lines(log(m.tran.ev.g/m.null.ev.g)~yval,col="green")
abline(h=0,lty=2)

plot (log(m.dev.ev.g/m.null.ev.g)~yval,type="l",col="red",ylab="Gamma evennes",ylim=c(-.2,.2), xlab="Time steps")
lines(log(m.res.ev.g/m.null.ev.g)~yval,col="blue")
lines(log(m.tran.ev.g/m.null.ev.g)~yval,col="green")
abline(h=0,lty=2)

plot (log(m.dev.beta/m.null.beta)~yval,type="l",col="red",ylab="Beta dispersion",ylim=c(-.05,.05), xlab="Time steps")
lines(log(m.res.beta/m.null.beta)~yval,col="blue")
lines(log(m.tran.beta/m.null.beta)~yval,col="green")
abline(h=0,lty=2)
 

# Write all the outputs to text files

write.table(null.ts.a,"Null_alpha.txt",quote=F,row.names=F,sep="\t")
write.table(null.ts.g,"Null_gamma.txt",quote=F,row.names=F,sep="\t")
write.table(null.ev.a,"Null_even_a.txt",quote=F,row.names=F,sep="\t")
write.table(null.ev.g,"Null_even_g.txt",quote=F,row.names=F,sep="\t")
write.table(null.beta,"Null_beta.txt",quote=F,row.names=F,sep="\t")

write.table(dev.ts.a,"Dev_alpha.txt",quote=F,row.names=F,sep="\t")
write.table(dev.ts.g,"Dev_gamma.txt",quote=F,row.names=F,sep="\t")
write.table(dev.ev.a,"Dev_even_a.txt",quote=F,row.names=F,sep="\t")
write.table(dev.ev.g,"Dev_even_g.txt",quote=F,row.names=F,sep="\t")
write.table(dev.beta,"Dev_beta.txt",quote=F,row.names=F,sep="\t")

write.table(res.ts.a,"Res_alpha.txt",quote=F,row.names=F,sep="\t")
write.table(res.ts.g,"Res_gamma.txt",quote=F,row.names=F,sep="\t")
write.table(res.ev.a,"Res_even_a.txt",quote=F,row.names=F,sep="\t")
write.table(res.ev.g,"Res_even_g.txt",quote=F,row.names=F,sep="\t")
write.table(res.beta,"Res_beta.txt",quote=F,row.names=F,sep="\t")

write.table(tran.ts.a,"Tran_alpha.txt",quote=F,row.names=F,sep="\t")
write.table(tran.ts.g,"Tran_gamma.txt",quote=F,row.names=F,sep="\t")
write.table(tran.ev.a,"Tran_even_a.txt",quote=F,row.names=F,sep="\t")
write.table(tran.ev.g,"Tran_even_g.txt",quote=F,row.names=F,sep="\t")
write.table(tran.beta,"Tran_beta.txt",quote=F,row.names=F,sep="\t")
