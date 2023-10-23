## L. Coggins, NOAA, SEFSC October 2023
## MetaPVA.R
### Builds on single stock Gulf Sturgeon PVA built by Brett Van Porten to 
### examine multi stock performance objectives 


rm(list = ls())
species <- "GulfSturgeon"


getLHpars <- function(LHdf=NULL){
  parNames <- LHdf$Parameter
  
  par <- list()
  for(pn in parNames){
    ind <- which(parNames == pn)
    val <- LHdf[ind,2]
    val <- as.numeric(unlist(strsplit(gsub(";", " ", val), " ")))
    par[[pn]] <- val
  }
  return(par)
}

"set.pars" <- function(Species=species){
  contrs <- read.csv(paste(Species,"controls_meta.csv",sep=" "),header=TRUE)
  controls <<- getLHpars(contrs)
  
  parm <- read.csv(paste(Species,"parameters_meta.csv",sep=" "),header=TRUE)
  parameters <<- getLHpars(parm)
  
}

controls <- list()
parameters <- list()
set.pars(species)

nT <- controls$nT
n.sim <- controls$n.sim
A <- controls$A
nS <- controls$nS
AR <- controls$AR
h.st <- controls$h.st
ep.fr <- controls$ep.fr
ep.J.M <- controls$ep.J.M
ep.A.M <- controls$ep.A.M
ep.spr <- controls$ep.spr
recfail <- controls$recfail
ep.fec <- controls$ep.fec
opt <- controls$opt
extir.threshold<-controls$extir.threshold

R0 <- parameters$R0                   # unfished equilibrium recruitment
reck <- parameters$reck               # recruitment compensation ratio
K <- parameters$K                     # von Bertalanffy metabolic parameter
Madult <- parameters$Madult           # vector of age-specific survival (after age at recruitment) [length=A-AR+1]
lh <- parameters$lh                   # length at 50% capture probability (logistic function)
afec <- parameters$afec               # slope of fecundity-weight relationship
Wmat <- parameters$Wmat               # weight at maturity
Ms <- parameters$Ms                   # maximum survival by stanza
Bs <- parameters$Bs                   # stanza-specific density effect (can be interpretted as amount of available habitat)
No <- parameters$No                   # initial vulnerable abundance (used to create initial population)
sd.S <- parameters$sd.S               # standard deviation of environmental effect on survival
n.pops<-parameters$n.pops             # number of populations 
moveProbs<-matrix(c(parameters$move.pear,
                    parameters$move.pasc,
                    parameters$move.esca,
                    parameters$move.yell,
                    parameters$move.choc,
                    parameters$move.apal,
                    parameters$move.suwa),n.pops,n.pops,byrow=T)

age <- AR:A                                     # ages
la <- (1-exp(-K*age))                           # unfished mean length-at-age
wa <- la^3                                      # unfished mean weight-at-age
fec <- (pmax(0,wa-Wmat)*afec)                   # unfished eggs at age
amat <- as.integer(-log(1-Wmat^(1/3))/K)        # age-at-maturity (used for differentiating subadults and adults)

Sa<-matrix(NA,n.pops,A)
for(j in 1:n.pops){Sa[j,] <- exp(-Madult[j]/la*0.66)}                           # length-based survival


lx<-matrix(NA,n.pops,A)
phie<-rep(NA,n.pops)
Rinit<-rep(NA,n.pops)              # initial recruitment given initial population size
PopN1<-rep(NA,n.pops)

for (j in 1:n.pops){
  lx[j,] <- c(1,Sa[j,1:(A-AR)])                         # incomplete survivorship to age
  lx[j,] <- cumprod(lx[j,])                               # survivorship
  phie[j] <- sum(lx[j,]*fec)                             # unfished eggs per recruit
  Rinit[j]<-No[j]/sum(lx[j,])
  PopN1[j]<-sum(Rinit[j]*lx[j,])             # initial abundance
}
  
R.A <- reck/phie                                # alpha of recruitment function 
R.B <- (reck-1)/(R0*phie)                       # beta of Beverton-Holt recruitment function for each population


moveProbs<-diag(n.pops)


Abun.Results=array(NA,c(nT,n.pops,n.sim))
Extir.Results=matrix(NA,n.sim,n.pops)

for(j in 1:n.sim){
#print(j)
start <- Sys.time()
rec.dev <- matrix(exp(rnorm((A-AR+1)*n.pops,0,sd.S)+0.5*sd.S^2),nrow=A-AR+1,ncol=n.pops)  # annual recruitment to initialize each population
anom <-    matrix(exp(rnorm((nT-1)*n.pops,0,sd.S)),nT-1,n.pops)


#Initialize Populations
Pops<-array(NA,c(nT,A,n.pops))
EtPops<-matrix(NA,nT,n.pops)
for (i in 1:n.pops){
  Pops[1,,i] <- round(Rinit[i]*lx[i,])  #;print(Pops[1,,i])
  Pops[1,AR:A,i] <- rmultinom(1,PopN1[i],lx[i,]*rec.dev[,i]/sum(lx[i,]*rec.dev[,i]))        # initial population 
  EtPops[,i] <-sum(Pops[1,AR:A,i]*fec)
  
}

#Episodic Definition
ep=0
if(ep.fr <= nT){
  n.samp <- as.integer(nT/ep.fr)
  ep <- replicate(1,sample(nT,n.samp,replace=FALSE))         # Years when episodes happen for each population
  if(n.samp>1){ep.p <- cbind(rep(NA,n.samp),ep[1])}else{ep.p <- cbind(rep(NA,n.samp),ep[1])}
  if(opt){
    ep <- rbind(ep,ep.p)
  }
  n.ep <- dim(ep.p)[1]
}

#Initialize Movement Array
moves<-array(NA,c(n.pops,A,n.pops))

#Compute the rest of the population dynamics
for(t in 2:nT){
    #episotic events
    ep.S <- rep(1,(A-AR+1))
    svec <- rep(1,(A-AR))
    if(t%in% ep){
      svec[(AR+1):amat-AR] <- 1-ep.J.M
      svec[amat:A-AR] <- 1-ep.A.M
      #print(svec)
    }
    
    for(i in 1:n.pops){
      #Pops[t,AR,i] <- rbinom(1,as.integer(pmax(0,EtPops[t-1,i])),pmin(R.A[i]*anom[t-1,i]/(1+R.B[i]*EtPops[t-1,i]),1))
      Pops[t,AR,i] <- as.integer((R.A[i]*EtPops[t-1,i]*anom[t-1,i])/(1+R.B[i]*EtPops[t-1,i]))
      Pops[t,AR,i] <- as.integer((R.A[i]*EtPops[t-1,i])/(1+R.B[i]*EtPops[t-1,i]))
      
      
      Pops[t,(AR+1):A,i] <- rbinom(A-AR,as.integer(pmax(0,Pops[t-1,AR:(A-1),i])),Sa[i,1:(A-AR)]*svec)  # survive fish to the next age
      Pops[t,(AR+1):A,i] <- (Pops[t-1,AR:(A-1),i]*Sa[i,1:(A-AR)]*svec)  # survive fish to the next age
      
      
      
      EtPops[t,i] <- as.integer(sum(Pops[t,AR:A,i]*fec));print(EtPops[t,i])
      moves[,,i]<-as.integer(sapply(Pops[t,,i],rmultinom,n=1,prob=moveProbs[i,]));sum(moves[,,i])
      #moves[,amat:A-AR,i]<-sapply(Pops[t,amat:A-AR,i],rmultinom,n=1,prob=moveProbs[i,]);sum(moves[,,i])
      
    }

#Redistribute the Fish after Movement    
    for(i in 1:n.pops){
      Pops[t,,i]<-rowSums(moves[i,,1:n.pops])
      #Pops[t,amat:A-AR,i]<-rowSums(moves[i,amat:A-AR,1:n.pops])
    }

    
  }



Abun.Results[,,j]<-apply(Pops,3,rowSums)
Extir.Results[j,]<-apply(Abun.Results[1:100,,j],2,min)<=extir.threshold
runtime <- Sys.time()-start
print(runtime)    


}

print(colSums(Extir.Results)/n.sim)

plot(1:200,rowSums(Pops[,,1]),type="l",ylim=c(0,8000),col=1)
lines(1:200,rowSums(Pops[,,2]),col=2)
lines(1:200,rowSums(Pops[,,3]),col=3)
lines(1:200,rowSums(Pops[,,4]),col=7)
lines(1:200,rowSums(Pops[,,5]),col=5)
lines(1:200,rowSums(Pops[,,6]),col=6)
lines(1:200,rowSums(Pops[,,7]),col=4)
legend("topleft",legend=c("Pearl","Pasc","Esca","Yell","Choc","Apal","Suwa"),lty=1,col=c(1:3,7,5,6,4),bty='n')

#matplot((Abun.Results[,1,]),type='l')








