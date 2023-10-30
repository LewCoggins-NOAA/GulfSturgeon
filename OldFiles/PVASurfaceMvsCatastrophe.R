## Brett van Poorten, Modified by L. Coggins, NOAA, SEFSC September 2023
## University of British Columbia, Institute for the Oceans and Fisheries
## PVA.R
### R code recreates Pine et al. 2013 PVA to examine viability scenarios
### uses binomial probabilities instead of multiple Bernoulli trials to speed up computation

rm(list = ls())

library(MASS)
library(ggplot2)
library(grid)
library(plotly)
library(fields)


#------------------------------------------------------------------------------#
#    NOTES
#------------------------------------------------------------------------------#

# can put in parameters for multiple populations so you just choose the 
# population you will be evaluating

#------------------------------------------------------------------------------#
#    PROBLEMS
#------------------------------------------------------------------------------#

# ERROR!!: YOUR ADULT ABUNDANCE IS NOT INDIVIDUAL-BASED. NEED TO PUT RBINOM INTO
#          MORTALITY CALCULATION

#------------------------------------------------------------------------------#
#    PROCEDURES
#------------------------------------------------------------------------------#

"PVA" <- function(M,epi.fr=101){   
#  "PVA" <- function(controls,parameters,habitat,TP.st){   
    start <- Sys.time()
  cat("Calculating population projections ...\n")
  # controls is a list including number of years (nT), number of ages (A),
  #   number of juvenile stanzas (nS), age at recruitment (AR), 
  #   and the number of simulations (n.sim)
  # parameters is a list of all parameters

  nT <- controls$nT
  n.sim <- controls$n.sim
  A <- controls$A
  nS <- controls$nS
  AR <- controls$AR
  h.st <- controls$h.st
#  TP.st <- controls$TP.st
  ep.fr <-epi.fr
  #ep.fr <- controls$ep.fr
  ep.J.M <- controls$ep.J.M
  ep.A.M <- controls$ep.A.M
  ep.spr <- controls$ep.spr
  recfail <- controls$recfail
  ep.fec <- controls$ep.fec
  opt <- controls$opt
#  habitat <- controls$habitat

  R0 <- parameters$R0           # unfished equilibrium recruitment
  reck <- parameters$reck       # recruitment compensation ratio
  K <- parameters$K             # von Bertalanffy metabolic parameter
#  Madult <- parameters$Madult   # vector of age-specific survival (after age at recruitment) [length=A-AR+1]
  Madult <- M   # vector of age-specific survival (after age at recruitment) [length=A-AR+1]
  lh <- parameters$lh           # length at 50% capture probability (logistic function)
  lsd <- parameters$lsd         # slope of logistic function at 50% capture probability
  afec <- parameters$afec       # slope of fecundity-weight relationship
  Wmat <- parameters$Wmat       # weight at maturity
  Ms <- parameters$Ms           # maximum survival by stanza
  Bs <- parameters$Bs           # stanza-specific density effect (can be interpretted as amount of available habitat)
  V1 <- parameters$V1           # initial vulnerable abundance (used to create initial population)
  sd.S <- parameters$sd.S       # standard deviation of environmental effect on survival

  # define the relative change in spawning habitat
  lam <- mapply('/',habitat,habitat[1,1:nS])

  # define stanza-specific recruitment parameters
  age <- AR:A                                     # ages
  la <- (1-exp(-K*age))                           # unfished mean length-at-age
  wa <- la^3                                      # unfished mean weight-at-age
  fec <- (pmax(0,wa-Wmat)*afec)                   # unfished eggs at age
  amat <- as.integer(-log(1-Wmat^(1/3))/K)        # age-at-maturity (used for differentiating subadults and adults)
  Sa <- exp(-Madult/la)                           # length-based survival
  va <- 1/(1+exp(-(la-lh)/lsd))                   # vulnerability
  lx <- c(1,Sa[1:(A-AR)])                         # incomplete survivorship to age
  lx <- cumprod(lx)                               # survivorship
  phie <- sum(lx*fec)                             # unfished eggs per recruit
  R.A <- reck/phie                                # alpha of recruitment function 
  R.B <- (reck-1)/(R0*phie)                       # beta of Beverton-Holt recruitment function
  Rinit <- V1/sum(lx*va)                          # initial recruitment given initial vulnerable population
  N1 <- sum(Rinit*lx)                             # initial abundance
  
  # stanza-specific recruitment parameters
  A.s <- exp(log(R.A)*Ms/sum(Ms))                  # maximum survival for each stanza
  A.s.temp <- c(1,A.s)
  den <- vector()
  for(i in 1:nS)
    den[i] <- Bs[i]*prod(A.s.temp[1:i])
  B.s <- Bs*R.B/sum(den)                           # carrying capacity parameter for each stanza
 
  # initialize population
  Nt <- array(dim=c(nT,A,n.sim))                   # numbers at time and age for each simulation
  Et <- matrix(nrow=nT,ncol=n.sim)                 # eggs for each year in each simulation
  rec.dev <- matrix(exp(rnorm((A-AR+1)*n.sim,0,sd.S)+0.5*sd.S^2),nrow=A-AR+1,ncol=n.sim)  # annual recruitment deviate
  for(i in 1:n.sim)
    Nt[1,AR:A,i] <- rmultinom(1,N1,lx*rec.dev[,i])        # initial population 
  Et[1,] <- as.integer(colSums(sweep(Nt[1,AR:A,],MARGIN=1,fec,'*')))      # eggs produced in first year
  W.st <- rnorm(nS*nT*n.sim,0,sd.S)
  anom <- array(data=exp(W.st),dim=c(nT,nS,n.sim))  # survival anomolies for each stanza
  A.s <- A.s*(1-h.st)                       # remove fish for translocation
  B.st <- B.s/t(lam)

  if(ep.fr <= nT){
    n.samp <- as.integer(nT/ep.fr)
    ep <- replicate(n.sim,sample(nT,n.samp,replace=FALSE))
    if(n.samp>1){ep.p <- cbind(rep(NA,n.samp),ep[,1:n.sim-1])}else{ep.p <- cbind(rep(NA,n.samp),ep[1:n.sim-1])}
    if(opt){
      ep <- rbind(ep,ep.p)
    }
    n.ep <- dim(ep.p)[1]
  }
  
  # dynamics for subsequent years
  for(t in 2:nT){
    # episodic events
    ep.S <- matrix(rep(1,(A-AR+1)*n.sim),nrow=A-AR+1,ncol=n.sim)  # episodic survival
    svec <- rep(1,(A-AR))
    rec.mult <- rep(1,n.sim)  # recruitment multiplier during episodic spring flow reduction events
    fec.mult <- rep(1,n.sim)  # fecundity multiplier during episodic spring flow events
    if(ep.fr <= nT){
      if(recfail) ep.fec <- 0
      svec[(AR+1):amat-AR] <- 1-ep.J.M
      svec[amat:A-AR+1] <- 1-ep.A.M
      if(is.matrix(ep)){x <- which(ep==t,arr.ind=TRUE)[,2]} # which of the simulations has an episodic even this year
      if(is.vector(ep)){x <- which(ep==t,arr.ind=TRUE)} # which of the simulations has an episodic even this year
      if(length(x)!=0){for(i in 1:length(x)) ep.S[,x[i]] <- svec } # simulations with an episodic event have annual survival adjusted
      if(ep.spr) rec.mult[x] <- 0.2  # if an episodic reduction in coldspring flow, reduce recruitment to by 80%
      fec.mult[x] <- ep.fec  # simulations with an episodic event have fecundity adjusted (if it is different than normal)
    }

    N.st <- matrix(nrow=nS+1,ncol=n.sim)    # numbers surviving through each stanza in a year
    N.st[1,] <- Et[t-1,]
    for(st in 1:nS){
      N.st[st,] <- N.st[st,]*rec.mult + TP.st[t,st]  # add stocked fish 
      N.st[st+1,] <- rbinom(n.sim,as.integer(pmax(0,N.st[st,])),pmin(A.s[st]*anom[t,st,]/(1+B.st[st,t]*N.st[st,]),1))
    }
    Nt[t,AR,] <- as.integer(N.st[nS+1,]*ep.S[1,])       # numbers surviving through all stanzas become recruits to the population
    for(i in 1:n.sim)
      Nt[t,(AR+1):A,i] <- rbinom(A-AR,pmax(0,Nt[t-1,AR:(A-1),i]),Sa[1:(A-AR)]*ep.S[2:(A-AR+1),i])  # survive fish to the next age
    Et[t,] <- as.integer(colSums(sweep(Nt[t,AR:A,],MARGIN=1,fec,'*')))*fec.mult
  }
  p.extinct.50 <- NA
  p.extinct.100 <- NA
  p.extinct.200 <- NA
  if( nT >= 50 ) p.extinct.50 <- length(which(colSums(Nt[50,,],na.rm=TRUE)==0))/n.sim
  if( nT >= 100 ) p.extinct.100 <- length(which(colSums(Nt[100,,],na.rm=TRUE)==0))/n.sim
  if( nT >= 200 ) p.extinct.200 <- length(which(colSums(Nt[200,,],na.rm=TRUE)==0))/n.sim
  runtime <- Sys.time()-start
  print(runtime)

  #out <- list()
  #out$R.A <- R.A
  #out$R.B <- R.B
  #out$A.s <- A.s
  #out$B.s <- B.s
  #out$B.st <- B.st
  #out$Nt <- Nt
  #out$Et <- Et
  #out$p.extinct.50 <- p.extinct.50
  #out$p.extinct.100 <- p.extinct.100
  #out$p.extinct.200 <- p.extinct.200
  #out$runtime <- runtime
  #out<-c(p.extinct.50,p.extinct.100,p.extinct.200)
  out<-p.extinct.100
  
  return(out)
}
  

#------------------------------------------------------------------------------#
#  SET PARAMETERS AND CONTROLS
#------------------------------------------------------------------------------#

species <- readline(prompt="Please enter the species prefix on files: ")
GulfSturgeon
cat("To evaluate a scenario, type 'scenario.switch( scenario ), where scenario
    can include (so far) 'base', 'epis', 'recfail', 'springs', 'hab', 'dam', 'stock',
    'remove.brood' or 'supplementation'\n")
cat("Or, type 'heat.proj() to plot base results.\n")
cat("New scenarios can be run by altering available habitat or changing controls in appropriate *.csv files.\n")
cat("New species or populations can be evaluated by changing parameters in appropriate *.csv files.\n")
#
controls <- list()
parameters <- list()
habitat <- list()
TP.st <- list()

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
  contrs <- read.csv(paste(Species,"controls.csv",sep=" "),header=TRUE)
  controls <<- getLHpars(contrs)
  
  parm <- read.csv(paste(Species,"parameters.csv",sep=" "),header=TRUE)
  parameters <<- getLHpars(parm)
  
  hab <- read.csv(paste(Species, "habitat.csv",sep=" "),header=TRUE)
  habitat <<- hab[,1+1:controls$nS]
  
  TP <- read.csv(paste(Species, "stocking.csv",sep=" "),header=TRUE)
  TP.st <<- TP[,1+1:controls$nS]
}
set.pars(species)


set.pars(species)
numEpFreq=15
numM=12
M=rep(seq(.05,.16,by=.01),numEpFreq)
epi.fr=rep(seq(10,80,by=5),each=12)

# This is the call to compute the response surface
response=mapply(PVA,M,epi.fr)

results=matrix(response,numM,numEpFreq)
res.df=as.data.frame(cbind(M,epi.fr,response))
obj<- list(x=M[1:numM], y=seq(10,80,by=5),z=results)
plot.surface( obj, type="c", col="red",xlab="M",ylab="CatastropheFreq")

ProbabilityExtripation100Years=(results)
Mortality=M[1:numM]
CatastropheFreq=seq(10,80,by=5)
fig1=plot_ly(z=~ProbabilityExtripation100Years,y=~Mortality,x=~CatastropheFreq)
fig1 <- fig1 %>% add_surface()
fig1

fig2=plot_ly(z=~ProbabilityExtripation100Years,y=~Mortality,x=~CatastropheFreq,type="contour")
fig2

