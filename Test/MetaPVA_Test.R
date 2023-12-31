## L. Coggins, NOAA, SEFSC October 2023
## MetaPVA.R
### Builds on single stock Gulf Sturgeon PVA built by Brett Van Porten to 
### examine multi stock performance objectives
#This version allows specification of minimum age of movement movement

#explores No versus RecK

#this version includes parallel processing


rm(list = ls())

library(plotly)
library(snow)

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
  contrs <- read.csv(paste(Species,"controls_meta_test.csv",sep=" "),header=TRUE)
  controls <<- getLHpars(contrs)
  
  parm <- read.csv(paste(Species,"parameters_meta_test.csv",sep=" "),header=TRUE)
  parameters <<- getLHpars(parm)
  
}

#moveProbs<-diag(n.pops)


MetaPVA<-function(fact1,fact2){
  
  start <- Sys.time()
  cat("Calculating population projections ...\n")
 

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
  
  #R0 <- parameters$R0                  # NO LONGER NEEDED BECAUSE SPECIFY UNFISHED EQUILIBRIUM ADULT ABUNDANCE
  reck <- parameters$reck               # recruitment compensation ratio
  K <- parameters$K                     # von Bertalanffy metabolic parameter
  Madult <- parameters$Madult           # vector of age-specific survival (after age at recruitment) [length=A-AR+1]
  lh <- parameters$lh                   # length at 50% capture probability (logistic function)
  afec <- parameters$afec               # slope of fecundity-weight relationship
  Wmat <- parameters$Wmat               # weight at maturity
  Ms <- parameters$Ms                   # maximum survival by stanza
  Bs <- parameters$Bs                   # stanza-specific density effect (can be interpreted as amount of available habitat)
  n.pops<-parameters$n.pops             # number of populations 
  No <- parameters$No                   # initial vulnerable abundance (used to create initial population)
  sd.S <- parameters$sd.S               # standard deviation of environmental effect on survival
  A.move<-parameters$age.move           # minimum age of fish that move between rivers  pear<-parameters$move.pear
  pear<-parameters$move.pear            # movement probabilities from each river to the other rivers
  pasc<-parameters$move.pasc        
  esca<-parameters$move.esca
  yell<-parameters$move.yell
  choc<-parameters$move.choc
  apal<-parameters$move.apal
  suwa<-parameters$move.suwa
  moveProbs<-matrix(c(pear,
                      pasc,
                      esca,
                      yell,
                      choc,
                      apal,
                      suwa),n.pops,n.pops,byrow=T)  # Movement probability matrix
  
  
  Nadu.eq<-parameters$NaduEq
  
#HERE IS WHERE WE CHANGE EVERYTHING FROM THE INPUT FILE PARAMETERS  
  #No<-rep(2000,n.pops)
  #Nadu.eq<-rep(500,n.pops)
  #moveProbs=diag(n.pops)
  #Madult=rep(.25,n.pops)
  reck <- reck*fact2                   
  #Madult <-Madult*fact2              
  Nadu.eq<-Nadu.eq*fact1
  
  
   
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
    #Rinit[j]<-No[j]/sum(lx[j,])
    Rinit[j]<-No[j]/sum(lx[j,amat:A])           #this fixes the problem that No scales with Madult
    PopN1[j]<-sum(Rinit[j]*lx[j,])             # initial abundance
  }
  
  #This part computes R0 from Nadu.eq 
  R0<-Nadu.eq/rowSums(lx[,amat:A])
  
  R.A <- reck/phie                                # alpha of recruitment function 
  R.B <- (reck-1)/(R0*phie)                       # beta of Beverton-Holt recruitment function for each population
  
Abun.Results=array(NA,c(nT,n.pops,n.sim))
Adult.Results=array(NA,c(nT,n.pops,n.sim))
Extir.Results=matrix(NA,n.sim,n.pops)
numpops.ext=rep(NA,n.sim)

rec.dev<- array(exp(rnorm((A-AR+1)*n.pops*n.sim,0,sd.S)+0.5*sd.S^2),c((A-AR+1),n.pops,n.sim))  # annual recruitment to initialize each population
anom <-   array(exp(rnorm((nT-1)*n.pops*n.sim,0,sd.S)),c((nT-1),n.pops,n.sim))
Pops<-array(NA,c(nT,A,n.pops,n.sim))
EtPops<-array(NA,c(nT,n.pops,n.sim))
ep.S <- array(rep(1,nT*(A-AR)*n.sim),c(nT,(A-AR),n.sim)) #This contains age specific annual survival rates for each simulation



#Episodic Definition
ep=0
if(ep.fr <= nT){
  n.samp <- as.integer(nT/ep.fr)
  ep <- replicate(n.sim,sample(nT,n.samp,replace=FALSE))         # Years when episodes happen for each population
}

svec=rep(1,A-1)
svec[(AR+1):amat-AR] <- 1-ep.J.M
svec[amat:(A-AR)] <- 1-ep.A.M

#Initialize Populations and set up Episodic survival Rates 
for (j in 1:n.sim){
  for (i in 1:n.pops){
    Pops[1,,i,j] <- round(Rinit[i]*lx[i,])  
    Pops[1,AR:A,i,j] <- rmultinom(1,PopN1[i],lx[i,]*rec.dev[,i,j]/sum(lx[i,]*rec.dev[,i,j]))        # initial population 
    EtPops[,i,j] <-sum(Pops[1,AR:A,i,j]*fec)
  }
  if(ep.fr <= nT){ep.S[ep[,j],,j]=matrix(rep(svec,dim(ep)[1]),dim(ep)[1],(A-AR),byrow=T)}
}

#Initialize Movement Array
  moves<-array(NA,c(n.pops,A,n.pops,n.sim))

#Compute the rest of the population dynamics
for(t in 2:nT){

    for(i in 1:n.pops){
      #Pops[t,AR,i,] <- rbinom(n.sim,as.integer(pmax(0,EtPops[t-1,i,])),pmin(R.A[i]*anom[t-1,i,]/(1+R.B[i]*EtPops[t-1,i,]),1))
      Pops[t,AR,i,] <- as.integer((R.A[i]*EtPops[t-1,i,]*anom[t-1,i,])/(1+R.B[i]*EtPops[t-1,i,]))
      #Pops[t,AR,i] <- as.integer((R.A[i]*EtPops[t-1,i,j])/(1+R.B[i]*EtPops[t-1,i,j]))
      
      for(j in 1:n.sim){      
        Pops[t,(AR+1):A,i,j] <- rbinom(A-AR,as.integer(pmax(0,Pops[t-1,AR:(A-1),i,j])),ep.S[t-1,,j]*Sa[i,AR:(A-1)])  # survive fish to the next age
        EtPops[t,i,j] <- as.integer(sum(Pops[t,AR:A,i,j]*fec));#print(EtPops[t,i,j])
        #moves[,A.move:(A),i,j]<-sapply(Pops[t,A.move:(A),i,j],rmultinom,n=1,prob=moveProbs[i,])  #this is a redundant line that computes movement within simulation loop, I think doing it this way is slightly slower than computing movement outside the simulation loop as below
      }
      moves[,A.move:(A),i,]<-sapply(Pops[t,A.move:(A),i,],rmultinom,n=1,prob=moveProbs[i,])#;sum(moves[,,i,])
      
    }
# 
#Redistribute the Fish after Movement
  for(j in 1:n.sim){
    for(i in 1:n.pops){
    Pops[t,A.move:(A),i,j]<-rowSums(moves[i,A.move:(A),,j])
    }
  }

 }
  
    

for(j in 1:n.sim){      
  Abun.Results[,,j]<-apply(Pops[,,,j],3,rowSums)
  Adult.Results[,,j]<-apply(Pops[,amat:A,,j],3,rowSums)
  Extir.Results[j,]<-apply(Adult.Results[1:100,,j],2,min)<=extir.threshold
} 

#Compute probability of the number of populations becoming extirpated  
  numpops.ext=rowSums(Extir.Results)
  prop.extir=sapply(0:(n.pops),function(i)1-length(numpops.ext[numpops.ext<=i])/n.sim);prop.extir
  #prop.extir=sapply(0:(n.pops),function(i)length(numpops.ext[numpops.ext<=i]));prop.extir
  
  #proportion below adult starting abundance at 25 to get a probability of decline
  Nadult20to25=matrix(NA,n.pops,n.sim)
  Nadult20to25=apply(Adult.Results[20:25,,],3,colMeans)
  Declines=rowMeans(Nadult20to25<Adult.Results[1,,]*.75)

  #timeseries of mean age for each river

  meanAge=matrix(NA,nT,n.pops)
  temp=array(NA,c(nT,n.pops,n.sim))
  for (t in 1:nT){
    for(i in 1:n.pops){
      for (j in 1:n.sim){
        temp[t,i,j]=sum(Pops[t,,i,j]*1:A)/sum(Pops[t,,i,j])
      }
    }
    meanAge[t,]=rowMeans(temp[t,,])
  }
  
runtime <- Sys.time()-start
print(runtime)    

out<-list(Extir.Res=Extir.Results,Abun.Res=Abun.Results,Adult.Res=Adult.Results,Decline.Res=Declines,MeanAge.Res=meanAge,PropExt.Res=prop.extir)

return(out)

}


numfact1Vals=5
numfact2Vals=5
fact1=rep(seq(0.1,0.5,length=numfact1Vals),numfact2Vals)
fact2=rep(seq(0.5,1.5,length=numfact2Vals),each=numfact1Vals)

# set number of cores to use
library(parallel)
ncpu = detectCores()-2
# ncpu = 14

# makes each row a list element
# we parLapply over this
# will evaluate MetaPVA at each combination of epAM and epfreq
# but each call to MetaPVA happens on a different core
params_list = lapply(1:(numfact1Vals * numfact2Vals), function(i) c(fact1 = fact1[i], fact2 = fact2[i]))

# start a timer
starttime = Sys.time()

# initialize a cluster for parallel computing
my_cluster = snow::makeSOCKcluster(ncpu)

# use this if you ever need to send the packages to the cluster
# also and suppresses the output from being printed to console
# invisible({
#   snow::clusterEvalQ(my_cluster, {library("pkg1"); library("pkg2")})
# })

# send the necessary objects from the current session to the cluster
snow::clusterExport(my_cluster, c("set.pars", "getLHpars", "MetaPVA"))

# perform the apply in parallel rather than on one cpu
response = snow::parLapply(my_cluster, x = params_list, fun = function(params) {
  Species <- "GulfSturgeon"
  controls <- list()
  parameters <- list()
  set.pars(Species)
  MetaPVA(fact1 = params["fact1"], fact2 = params["fact2"])
})

# stop the cluster and timer
snow::stopCluster(my_cluster)
stoptime = Sys.time()
cat("\nElapsed:", format(stoptime - starttime, digits = 2))

# extract the extirpation outcomes at each combo of control parameters
Extir.Res = lapply(1:length(response), function(i) response[[i]]$Extir.Res)
Decline.Res = lapply(1:length(response), function(i) response[[i]]$Decline.Res)
MeanAge.Res = lapply(1:length(response), function(i) response[[i]]$MeanAge.Res)

#Get Control and Parameter values to main environment
Species <- "GulfSturgeon"
controls <- list()
parameters <- list()
set.pars(Species)
n.pops<-parameters$n.pops
Nadu.eq<-parameters$NaduEq

ProbExt1orMore=matrix(sapply(1:(numfact1Vals * numfact2Vals),function(i)response[[i]]$PropExt.Res[1]),numfact2Vals,numfact1Vals)
ProbExt2orMore=matrix(sapply(1:(numfact1Vals * numfact2Vals),function(i)response[[i]]$PropExt.Res[2]),numfact2Vals,numfact1Vals)
ProbExt3orMore=matrix(sapply(1:(numfact1Vals * numfact2Vals),function(i)response[[i]]$PropExt.Res[3]),numfact2Vals,numfact1Vals)
ProbExt4orMore=matrix(sapply(1:(numfact1Vals * numfact2Vals),function(i)response[[i]]$PropExt.Res[4]),numfact2Vals,numfact1Vals)
ProbExt5orMore=matrix(sapply(1:(numfact1Vals * numfact2Vals),function(i)response[[i]]$PropExt.Res[5]),numfact2Vals,numfact1Vals)
ProbExt6orMore=matrix(sapply(1:(numfact1Vals * numfact2Vals),function(i)response[[i]]$PropExt.Res[6]),numfact2Vals,numfact1Vals)
ProbExtEQ7=matrix(sapply(1:(numfact1Vals * numfact2Vals),function(i)response[[i]]$PropExt.Res[7]),numfact2Vals,numfact1Vals)
ProbExtEQ0=1-ProbExt1orMore

# calculate extirpation probabilities
# ready to be passed to the remainder of the summary prep/plotting code
# NOTE: colMeans is cleaner way to calculate proportion of TRUE elements than colSums/n.sim
extir = sapply(Extir.Res, colMeans)


fact1=matrix(fact1,numfact1Vals,numfact2Vals)[,1]
fact2=matrix(fact2,numfact1Vals,numfact2Vals)[1,]
#EqAduAbun=matrix(NA,n.pops,length(EquilibriumAdultAbundFactor))
#for(i in 1:7){EqAduAbun[i,]=round(matrix(EquilibriumAdultAbundFactor,n.pops,length(EquilibriumAdultAbundFactor),byrow=TRUE)[i,]*Nadu.eq[i])}

Pearl_ProbabilityExtripation100Years=matrix(extir[1,],numfact2Vals,numfact1Vals)
Pascagoula_ProbabilityExtripation100Years =matrix(extir[2,],numfact2Vals,numfact1Vals)
Escambia_ProbabilityExtripation100Years =matrix(extir[3,],numfact2Vals,numfact1Vals)
Yellow_ProbabilityExtripation100Years =matrix(extir[4,],numfact2Vals,numfact1Vals)
Choctawhatchee_ProbabilityExtripation100Years =matrix(extir[5,],numfact2Vals,numfact1Vals)
Apalachicola_ProbabilityExtripation100Years =matrix(extir[6,],numfact2Vals,numfact1Vals)
Suwannee_ProbabilityExtripation100Years =matrix(extir[7,],numfact2Vals,numfact1Vals)

contours = list(showlabels = TRUE,start = 0,end = 1,labelfont = list(size = 12, color = "lightgray"))
fig=plot_ly(z=~Pearl_ProbabilityExtripation100Years,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~Pascagoula_ProbabilityExtripation100Years,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~Escambia_ProbabilityExtripation100Years,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~Yellow_ProbabilityExtripation100Years,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~Choctawhatchee_ProbabilityExtripation100Years,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~Apalachicola_ProbabilityExtripation100Years,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~Suwannee_ProbabilityExtripation100Years,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig

fig=plot_ly(z=~ProbExtEQ0,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExt1orMore,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExt2orMore,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExt3orMore,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExt4orMore,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExt5orMore,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExt6orMore,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExtEQ7,y=~fact1,x=~fact2,type="contour",colors="Greys",contours=contours);fig




windows(record = T)
# matplot(response[[1]]$Adult.Res[,,1],type='l');abline(h=50,lwd=2,col='red')
# matplot(response[[12]]$Adult.Res[,,1],type='l');abline(h=50,lwd=2,col='red')
# matplot(response[[25]]$Adult.Res[,,1],type='l');abline(h=50,lwd=2,col='red')
# 
# 
 
# matplot(response[[1]]$Adult.Res[,1,],type='l');abline(h=50,lwd=2,col='red')
# matplot(response[[25]]$Adult.Res[,1,],type='l');abline(h=50,lwd=2,col='red')
# matplot(response[[50]]$Adult.Res[,1,],type='l');abline(h=50,lwd=2,col='red')

 for(i in 1:(numfact1Vals*numfact2Vals)){ matplot(response[[i]]$Adult.Res[,1,],type='l',ylim=c(0,1000));abline(h=50,lwd=2,col='red')}
 
 # 
# matplot(response[[1]]$Adult.Res[,7,],type='l');abline(h=50,lwd=2,col='red')
# matplot(response[[12]]$Adult.Res[,7,],type='l');abline(h=50,lwd=2,col='red')
# matplot(response[[25]]$Adult.Res[,7,],type='l');abline(h=50,lwd=2,col='red')
