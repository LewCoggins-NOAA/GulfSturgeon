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

R0 <- parameters$R0           # unfished equilibrium recruitment
reck <- parameters$reck       # recruitment compensation ratio
K <- parameters$K             # von Bertalanffy metabolic parameter
Madult <- parameters$Madult   # vector of age-specific survival (after age at recruitment) [length=A-AR+1]
lh <- parameters$lh           # length at 50% capture probability (logistic function)
afec <- parameters$afec       # slope of fecundity-weight relationship
Wmat <- parameters$Wmat       # weight at maturity
Ms <- parameters$Ms           # maximum survival by stanza
Bs <- parameters$Bs           # stanza-specific density effect (can be interpretted as amount of available habitat)
V1 <- parameters$V1           # initial vulnerable abundance (used to create initial population)
sd.S <- parameters$sd.S       # standard deviation of environmental effect on survival


age <- AR:A                                     # ages
la <- (1-exp(-K*age))                           # unfished mean length-at-age
wa <- la^3                                      # unfished mean weight-at-age
fec <- (pmax(0,wa-Wmat)*afec)                   # unfished eggs at age
amat <- as.integer(-log(1-Wmat^(1/3))/K)        # age-at-maturity (used for differentiating subadults and adults)
Sa <- exp(-Madult/la)                           # length-based survival
lx <- c(1,Sa[1:(A-AR)])                         # incomplete survivorship to age
lx <- cumprod(lx)                               # survivorship
phie <- sum(lx*fec)                             # unfished eggs per recruit
R.A <- reck/phie                                # alpha of recruitment function 
R.B <- (reck-1)/(R0*phie)                       # beta of Beverton-Holt recruitment function
n.pops<-5


Rinit=rep(NA,n.pops)              # initial recruitment given initial population
Rinit[1]=V1/sum(lx)
ScalePop2to5=c(1.5,2.2,2.4,2.6)
Rinit[2:5]=Rinit[1]*ScalePop2to5

PopN1=rep(NA,n.pops)
PopN1=colSums(sapply(Rinit,'*',lx))             # initial abundance

moveProbs=matrix(c(.3,.1,.1,.25,.25,
                   .0,.6,.1,.2,.1,
                   .6 ,.1 ,.1,.1,.1,
                   .1 ,.1 ,.1,.4,.3,
                   .1 ,.1 ,.1,.2,.5),5,5,byrow=T)

#moveProbs=diag(5)


Pops<-array(NA,c(nT,A,n.pops))

#Pop1<-matrix(NA,nT,A)
#Pop2<-matrix(NA,nT,A)
#Pop3<-matrix(NA,nT,A)
#Pop4<-matrix(NA,nT,A)
#Pop5<-matrix(NA,nT,A)


Pops[1,,1] <- round(Rinit[1]*lx)  ;Pops[1,,1]                         
Pops[1,,2] <- round(Rinit[2]*lx)  ;Pops[1,,2]
Pops[1,,3] <- round(Rinit[3]*lx)  ;Pops[1,,3]
Pops[1,,4] <- round(Rinit[4]*lx)  ;Pops[1,,4]
Pops[1,,5] <- round(Rinit[5]*lx)  ;Pops[1,,5]

rec.dev <- matrix(exp(rnorm((A-AR+1)*n.pops,0,sd.S)+0.5*sd.S^2),nrow=A-AR+1,ncol=n.pops)  # annual recruitment to initialize each population
anom <- matrix(exp(rnorm((nT-1)*n.pops,0,sd.S)),nT-1,n.pops)


Pops[1,AR:A,1] <- rmultinom(1,PopN1[1],lx*rec.dev[,1]/sum(lx*rec.dev[,1]))        # initial population 
Pops[1,AR:A,2] <- rmultinom(1,PopN1[2],lx*rec.dev[,2]/sum(lx*rec.dev[,2]))        # initial population 
Pops[1,AR:A,3] <- rmultinom(1,PopN1[3],lx*rec.dev[,3]/sum(lx*rec.dev[,3]))        # initial population 
Pops[1,AR:A,4] <- rmultinom(1,PopN1[4],lx*rec.dev[,4]/sum(lx*rec.dev[,4]))        # initial population 
Pops[1,AR:A,5] <- rmultinom(1,PopN1[5],lx*rec.dev[,5]/sum(lx*rec.dev[,5]))        # initial population 


EtPops<-matrix(NA,nT,n.pops)

#EtPop1 <- rep(NA,nT)                 # eggs for each year in each simulation
#EtPop2 <- rep(NA,nT)                 # eggs for each year in each simulation
#EtPop3 <- rep(NA,nT)                 # eggs for each year in each simulation
#EtPop4 <- rep(NA,nT)                 # eggs for each year in each simulation
#EtPop5 <- rep(NA,nT)                 # eggs for each year in each simulation

EtPops[,1] <-sum(Pops[1,AR:A,1]*fec)
EtPops[,2] <-sum(Pops[1,AR:A,2]*fec)
EtPops[,3] <-sum(Pops[1,AR:A,3]*fec)
EtPops[,4] <-sum(Pops[1,AR:A,4]*fec)
EtPops[,5] <-sum(Pops[1,AR:A,5]*fec)


moves=array(NA,c(n.pops,A,n.pops))

for(t in 2:nT){
    for(i in 1:n.pops){
      Pops[t,AR,i] <- rbinom(1,as.integer(pmax(0,EtPops[t-1,i])),pmin(R.A*anom[t-1,i]/(1+R.B*EtPops[t-1,i]),1))
      Pops[t,(AR+1):A,i] <- rbinom(A-AR,pmax(0,Pops[t-1,AR:(A-1),i]),Sa[1:(A-AR)])  # survive fish to the next age
      EtPops[t,i] <- as.integer(sum(Pops[t,AR:A,i]*fec))
      moves[,,i]=sapply(Pops[t,,i],rmultinom,n=1,prob=moveProbs[i,]);sum(moves[,,i])
    }
  
    for(i in 1:n.pops){
      Pops[t,,i]=rowSums(moves[i,,1:n.pops])
    }

  }

plot(1:200,rowSums(Pops[,,1]),type="l",ylim=c(0,50000),col=1)
lines(1:200,rowSums(Pops[,,2]),col=2)
lines(1:200,rowSums(Pops[,,3]),col=3)
lines(1:200,rowSums(Pops[,,4]),col=4)
lines(1:200,rowSums(Pops[,,5]),col=5)








