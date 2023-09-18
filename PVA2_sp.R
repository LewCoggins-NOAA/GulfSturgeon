## Brett van Poorten
## University of British Columbia, Institute for the Oceans and Fisheries
## PVA.R
### R code recreates Pine et al. 2013 PVA to examine viability scenarios
### uses binomial probabilities instead of multiple Bernoulli trials to speed up computation

library(MASS)
library(ggplot2)
library(grid)

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

"PVA" <- function(controls,parameters,habitat,TP.st){   
  start <- Sys.time()
  cat("Calculating population projections ...\n")
  # controls is a list including number of years (nT), number of ages (A),
  #   number of juvenile stanzas (nS), age at recruitment (AR), 
  #   and the number of simulations (n.sim)
  # parameters is a list of all parameters

  nT <- controls$nT
  A <- controls$A
  nS <- controls$nS
  AR <- controls$AR
  n.sim <- controls$n.sim
  h.st <- controls$h.st
#  TP.st <- controls$TP.st
  ep.fr <- controls$ep.fr
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
  Madult <- parameters$Madult   # vector of age-specific survival (after age at recruitment) [length=A-AR+1]
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
    ep.p <- cbind(rep(NA,n.samp),ep[,1:n.sim-1])
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
      x <- which(ep==t,arr.ind=TRUE)[,2]  # which of the simulations has an episodic even this year
      for(i in 1:length(x))
        ep.S[,x[i]] <- svec  # simulations with an episodic event have annual survival adjusted
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

  out <- list()
  out$R.A <- R.A
  out$R.B <- R.B
  out$A.s <- A.s
  out$B.s <- B.s
  out$B.st <- B.st
  out$Nt <- Nt
  out$Et <- Et
  out$p.extinct.50 <- p.extinct.50
  out$p.extinct.100 <- p.extinct.100
  out$p.extinct.200 <- p.extinct.200
  out$runtime <- runtime
  
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

#------------------------------------------------------------------------------#
#  DEFINE SCENARIOS
#------------------------------------------------------------------------------#

scenario.switch <- function( 
  scenario = c("base", "epis", "recfail", "springs", "hab",
               "dam", "stock","remove.brood","supplementation"), 
  set.ymax = TRUE, #2000                                             ###### STEPHEN SWITCHED FOR GULFWIDE
  plot.switch=TRUE){  
  # First enter the scenario you want to evaluate; the function will ask for control parameters;
  # Set.ymax controls the maximum extent of the y-axis
  # plot.switch (logical) controls whether the time-series is plotted (TRUE) or if values are returned (FALSE)
  nT <- controls$nT
  set.pars()
  switch( scenario, 
          base = {
            parameters$Madult <- as.numeric(readline(prompt="Enter adult mortality: "))#rel.chg
            parameters$V1 <- as.numeric(readline(prompt="Enter initial vulnerable abundance: "))#V1
          },
          epis = {
            parameters$Madult <- as.numeric(readline(prompt="Enter adult mortality: "))#rel.chg
            parameters$V1 <- as.numeric(readline(prompt="Enter initial vulnerable abundance: "))#V1
            controls$ep.fr <- as.integer(readline(prompt="How many years between episodic mortality events: "))#freq
            controls$ep.J.M <- as.numeric(readline(prompt="Enter additional juvenile mortality during episodic events: "))#0.8
            controls$ep.A.M <- as.numeric(readline(prompt="Enter additional adult mortality during episodic events: "))#0.8
          },
          recfail = {
            parameters$Madult <- as.numeric(readline(prompt="Enter adult mortality: "))#rel.chg
            parameters$V1 <- as.numeric(readline(prompt="Enter initial vulnerable abundance: "))#V1
            controls$ep.fr <- as.integer(readline(prompt="How many years between episodic recruitment failure events: "))#freq
            controls$recfail <- TRUE
          },
          springs = {
            parameters$Madult <- as.numeric(readline(prompt="Enter adult mortality: "))#rel.chg
            parameters$V1 <- as.numeric(readline(prompt="Enter initial vulnerable abundance: "))#V1
            controls$ep.fr <- as.integer(readline(prompt="How many years between episodic coldspring flow reduction events: "))#freq
            controls$ep.spr <- TRUE
            controls$opt <- as.logical(readline(prompt="Do springs dry for two years (TRUE) or just one (FALSE):"))
            controls$ep.A.M <- as.numeric(readline(prompt="What is the additional adult mortality imposed during \n  coldspring flow reduction events: "))#0.001
            controls$ep.fec <- as.numeric(readline(prompt="What is the percent change in egg deposition \n  during coldspring flow reduction events (no change=100): "))/100#0.8
          },
          hab = {
            parameters$Madult <- as.numeric(readline(prompt="Enter adult mortality: "))#rel.chg
            parameters$V1 <- as.numeric(readline(prompt="Enter initial vulnerable abundance: "))#V1
            hab.chg <- as.numeric(readline(prompt="What is the percent change in habitat capacity (no change=100): "))/100
            st.imp <- as.integer(readline(prompt="Which stanza is impacted by this habitat change (between 1-nS): "))
            habitat[2:nT,st.imp] <- habitat[2:nT,st.imp]*hab.chg
          },
          dam = {
            parameters$Madult <- as.numeric(readline(prompt="Enter adult mortality: "))#rel.chg
            parameters$V1 <- as.numeric(readline(prompt="Enter initial vulnerable abundance: "))#V1
            dam.yr <- 1+as.integer(readline(prompt="What year is the dam installed/removed: "))
            dam.chg <- as.numeric(readline(prompt="What is the percent change in available habitat \n  with the addition/removal of a dam (no change=100): "))/100
            habitat[dam.yr:nT,] <- habitat[dam.yr:nT,]*dam.chg
          },
          stock = {
            parameters$Madult <- as.numeric(readline(prompt="Enter adult mortality: "))#rel.chg
            parameters$V1 <- as.numeric(readline(prompt="Enter initial vulnerable abundance: "))#V1
            stock.yr <- as.integer(readline(prompt="Which year does stocking start: "))
            stock.st <- as.integer(readline(prompt="Which stanza are fish stocked into: "))
            stock.num <- as.numeric(readline(prompt="How many fish are stocked annually: "))
            TP.st[stock.yr:nT,stock.st] <- stock.num
          },
          remove.brood = {
            parameters$Madult <- as.numeric(readline(prompt="Enter adult mortality: "))#rel.chg
            parameters$V1 <- as.numeric(readline(prompt="Enter initial vulnerable abundance: "))#V1
            ns <- controls$nS
            h.st <- vector()
            for(i in 1:ns)
              h.st[i] <- as.numeric(readline(prompt=paste("What proportion of eggs/larvae are removed from stanza ",i,": ",sep="")))
            controls$h.st <- h.st
          },
          supplementation = {
            parameters$Madult <- as.numeric(readline(prompt="Enter adult mortality: "))#rel.chg
            parameters$V1 <- as.numeric(readline(prompt="Enter initial vulnerable abundance: "))#V1
            ns <- controls$nS
            h.st <- vector()
            for(i in 1:ns)
              h.st[i] <- as.numeric(readline(prompt=paste("What proportion of eggs/larvae are removed from stanza ",i,": ",sep="")))
            controls$h.st <- h.st
            stock.yr <- as.integer(readline(prompt="Which year does stocking start: "))
            stock.st <- as.integer(readline(prompt="Which stanza are fish stocked into: "))
            stock.num <- as.numeric(readline(prompt="How many fish are stocked annually: "))
            TP.st[stock.yr:nT,stock.st] <- stock.num
          })
  heat.proj(controls,parameters,habitat,TP.st,set.ymax,plot.switch)
}


#  plots
#------------------------------------------------------------------------------#

"heat.proj" <- function( Controls=controls, Parameters=parameters, Habitat=habitat,
                         TP=TP.st, set.ymax = TRUE, plot.switch = TRUE ){            #### STEPHEN CHANGED YMAX FROM 2000
  ## set.ymax either set to TRUE (automatic calculation of plotted upper y value
  ##    or set as the desired upper limit (e.g. set.ymax=5000))
  n.sim <- Controls$n.sim
  nT <- Controls$nT
  pva <- PVA(Controls,Parameters,Habitat,TP)
  Na <- apply(pva$Nt,MARGIN=c(1,3),sum,na.rm=TRUE)
  Data <- data.frame(Na)
  names(Data) <- rep("Y",n.sim)
  Data$X <- 1:nT
  cat("Probability of extirpation after:\n")
  cat("   50 years - ",pva$p.extinct.50*100,"%\n")
  cat("   100 years - ",pva$p.extinct.100*100,"%\n")
  cat("   200 years - ",pva$p.extinct.200*100,"%\n")
  ymax <- set.ymax
  
  if(plot.switch)
    return(vwReg2(data=Data,controls=Controls,set.ymax=ymax))
  if( !plot.switch ) 
    return(pva)
}

vwReg2 <- function(data,controls,palette=colorRampPalette(c("purple4","blue","green","yellow","orange","red"), bias=2, space="rgb")(40),show.CI=FALSE,set.ymax=TRUE){
  nT <- controls$nT
  n.sim <- controls$n.sim
  
  # compute median and CI limits of sample
  library(plyr)
  library(reshape2)
  CI <- adply(data[,1:nT], 1, function(x) quantile(x, prob=c(.025, .5, .975, pnorm(c(-3, -2, -1, 0, 1, 2, 3))), na.rm=TRUE))[, -1]
  colnames(CI)[1:10] <- c("LL", "M", "UL", paste0("SD", 1:7))
  CI$X <- 1:nT
  CI$width <- CI$UL - CI$LL
  
  # scale the CI width to the range 0 to 1 and flip it (bigger numbers = narrower CI)
  CI$w2 <- (CI$width - min(CI$width))
  CI$w3 <- 1-(CI$w2/max(CI$w2))
  
  # convert spaghettis to long format
  b2 <- melt(as.matrix(data)[,1:n.sim])
  b2$X <- rep(1:nT,n.sim)
  colnames(b2) <- c("index", "rep", "value", "X")
  
  library(ggplot2)
  library(RColorBrewer)
  
  # Construct ggplot
  # All plot elements are constructed as a list, so they can be added to an existing ggplot
  
  p0 <- ggplot(data, aes(x=X, y=Y)) + theme_bw()
  
  # initialize elements with NULL (if they are defined, they are overwritten with something meaningful)
  gg.tiles <- gg.poly <- gg.spag <- gg.median <- gg.CI1 <- gg.CI2 <- gg.lm <- gg.points <- gg.title <- NULL
  
  cat("Computing density estimates for each vertical cut ...\n")
  flush.console()
  ymax <- as.integer(log10(max(data)))-1
  
  ifelse(set.ymax==TRUE,
         ylim <- c(0,as.integer(max(data)/ymax+1)*ymax),
         ylim <- c(0,set.ymax))
  d2 <- ddply(b2[, c("X", "value")], .(X), function(df) {
    res <- data.frame(density(df$value, na.rm=TRUE, n=n.sim, bw=ylim[2]/100,from=ylim[1], to=ylim[2])[c("x", "y")])
    colnames(res) <- c("Y", "dens")
    return(res)
  }, .progress="text")
  
  maxdens <- max(d2$dens,na.rm=TRUE)
  mindens <- min(d2$dens,na.rm=TRUE)
  d2$Density <- (d2$dens - mindens)/maxdens
  d2$X<-data$X[d2$X]
  
  ## Tile approach
  shade.alpha = 0.1
  d2$alpha.factor <- d2$Density^shade.alpha
  gg.tiles <-  list(geom_tile(data=d2, aes(x=X, y=Y, fill=Density, alpha=alpha.factor)), scale_fill_gradientn("Probability\ndensity\n", colours=palette), scale_alpha_continuous(range=c(0.001, 0.999),guide="none"))
  
  # Confidence limits
  if (show.CI == TRUE) {
    gg.CI1 <- geom_path(data=CI, aes(x=X, y=UL), size=1, color="red")
    gg.CI2 <- geom_path(data=CI, aes(x=X, y=LL), size=1, color="red")
  }
  
  cat("Build ggplot figure ...\n")
  flush.console()
  
  gg.elements <- list(gg.tiles, gg.poly, gg.spag, gg.median, gg.CI1, gg.CI2, gg.lm, gg.points, gg.title)#, gg.labs)#, theme(legend.position="none"))
  
  return(p0 + gg.elements + gg.CI1 + gg.CI2 + xlab("Year") + ylab("Recruited population numbers") + theme(text = element_text(size=18), legend.key.height= unit(2,"cm")))
}

# heat.proj()

# The scenario “base” allows evaluation of adult mortality and initial abundance 
#   to address uncertainty in these two critical parameters. The user is asked to 
#   supply these two values and then runs the model conditional on all other parameters 
#   and controls defined in set.pars(). 
scenario.switch("base")

# The scenario “epis” evaluates the impacts of episodic events on population persistence. 
# Depending on the nature of the event and the habitat in which it occurs, these events 
#   could include high and immediate additional mortality on either sub-adults or adults. 
# Sub-adults are recruited fish younger than the age at maturity.
# The user is asked to supply adult mortality and initial abundance, then the frequency of an episodic event occurring. 
# Finally, the user supplies the additional mortality on juvenile and adult fish. 
# Episodic mortality events occur in random years across simulations at the mean frequency determined by the user.
scenario.switch("epis")

# The scenario “recfail” evaluates the impact of occasional episodic recruitment failures on population persistence. 
# The user is asked to supply adult mortality and initial abundance, then the mean frequency of recruitment failure. 
# Episodic recruitment failure events occur in random years across simulations at the mean specified by the user.
scenario.switch("recfail")

# The scenario “springs” simulates an episodic event affecting consumption or metabolism of adults, 
#   which in turn may elevate adult mortality for that year as well as impact likelihood of spawning, 
#   and therefore aggregate fecundity for the population. 
# In Gulf sturgeon, this may occur as coldspring flows decline or cease in some years, 
#   which effects metabolism and energy reserves of adults resting during the summer. 
# Other species may have similar environmental bottlenecks. 
# The user is asked to supply adult mortality and initial abundance, then the frequency of an episodic event occurring.
# The user is then asked whether episodic events occur for a single year or two consecutive years, 
#   as well as the change to adult mortality and population fecundity. 
scenario.switch("springs")

# The scenario “hab” simulates the effects of changes to habitat quality or quantity. 
# Changes to habitat should not affect subadults or adults should not affect the population 
#   because there is no density dependent regulation at those stages: 
#   only pre-recruits experience density dependent survival. 
# Assuming pre-recruits experience some level of ontogenetic habitat shift 
#   (e.g. larvae settle in one habitat, juveniles move to another habitat to rear), 
#   habitat changes will affect one or more pre-recruit stanzas. 
# The user is asked to supply adult mortality and initial abundance, 
#   followed by an indication of the relative change (positive or negative) in habitat capacity as a percent. 
# Finally, the user is asked which pre-recruit stanza uses the changed habitat. 
scenario.switch("hab")






