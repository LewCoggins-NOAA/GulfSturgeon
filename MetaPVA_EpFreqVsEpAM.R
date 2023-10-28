## L. Coggins, NOAA, SEFSC October 2023
## MetaPVA.R
### Builds on single stock Gulf Sturgeon PVA built by Brett Van Porten to 
### examine multi stock performance objectives 

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
  contrs <- read.csv(paste(Species,"controls_meta.csv",sep=" "),header=TRUE)
  controls <<- getLHpars(contrs)
  
  parm <- read.csv(paste(Species,"parameters_meta.csv",sep=" "),header=TRUE)
  parameters <<- getLHpars(parm)
  
}

#moveProbs<-diag(n.pops)


MetaPVA<-function(epfreq,epAM){
  
  start <- Sys.time()
  cat("Calculating population projections ...\n")
 
  
  #ep.fr <- epfreq[4]
  #ep.A.M <- epAM[4]
   

  nT <- controls$nT
  n.sim <- controls$n.sim
  A <- controls$A
  nS <- controls$nS
  AR <- controls$AR
  h.st <- controls$h.st
  ep.fr <- epfreq
  #ep.fr <- controls$ep.fr
  ep.J.M <- controls$ep.J.M
  ep.A.M <- epAM
  #ep.A.M <- controls$ep.A.M
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
  Bs <- parameters$Bs                   # stanza-specific density effect (can be interpreted as amount of available habitat)
  No <- parameters$No                   # initial vulnerable abundance (used to create initial population)
  sd.S <- parameters$sd.S               # standard deviation of environmental effect on survival
  n.pops<-parameters$n.pops             # number of populations 
  pear<-parameters$move.pear
  pasc<-parameters$move.pasc
  esca<-parameters$move.esca
  yell<-parameters$move.yell
  choc<-parameters$move.choc
  apal<-parameters$move.apal
  suwa<- parameters$move.suwa
  moveProbs<-matrix(c(pear,
                     pasc,
                     esca,
                     yell,
                     choc,
                     apal,
                     suwa),n.pops,n.pops,byrow=T)
  
  
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
    #episotic events
#    ep.S <- rep(1,(A-AR+1))
    
    for(i in 1:n.pops){
      Pops[t,AR,i,] <- rbinom(n.sim,as.integer(pmax(0,EtPops[t-1,i,])),pmin(R.A[i]*anom[t-1,i,]/(1+R.B[i]*EtPops[t-1,i,]),1))
      #Pops[t,AR,i] <- as.integer((R.A[i]*EtPops[t-1,i,j]*anom[t-1,i,j])/(1+R.B[i]*EtPops[t-1,i,j]))
      #Pops[t,AR,i] <- as.integer((R.A[i]*EtPops[t-1,i,j])/(1+R.B[i]*EtPops[t-1,i,j]))
      
      for(j in 1:n.sim){      
        Pops[t,(AR+1):A,i,j] <- rbinom(A-AR,as.integer(pmax(0,Pops[t-1,AR:(A-1),i,j])),ep.S[t-1,,j]*Sa[i,AR:(A-1)])  # survive fish to the next age
        EtPops[t,i,j] <- as.integer(sum(Pops[t,AR:A,i,j]*fec));#print(EtPops[t,i,j])
        #moves[,,i,j]<-as.integer(sapply(Pops[t,,i,j],rmultinom,n=1,prob=moveProbs[i,]))
      }
      moves[,,i,]<-(sapply(Pops[t,,i,],rmultinom,n=1,prob=moveProbs[i,]))
      
    }
# 
#Redistribute the Fish after Movement
  for(j in 1:n.sim){
    for(i in 1:n.pops){
    Pops[t,,i,j]<-rowSums(moves[i,,1:n.pops,j])
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
        temp[t,i,j]=sum(Pops[t,,i,j]*1:60)/sum(Pops[t,,i,j])
      }
    }
    meanAge[t,]=rowMeans(temp[t,,])
  }
  
  
  
  
  
  
runtime <- Sys.time()-start
print(runtime)    

out<-list(Extir.Res=Extir.Results,Abun.Res=Abun.Results,Adult.Res=Adult.Results,Decline.Res=Declines,MeanAge.Res=meanAge,PropExt.Res=prop.extir)

return(out)

}


numepfreqVals=5
numepAM=11
epAM=rep(seq(.25,.75,length=numepAM),numepfreqVals)
epfreq=rep(seq(5,25,by=5),each=numepAM)


# set number of cores to use
ncpu = 14

# makes each row a list element
# we parLapply over this
# will evaluate MetaPVA at each combination of epAM and epfreq
# but each call to MetaPVA happens on a different core
ep_params_list = lapply(1:(numepfreqVals * numepAM), function(i) c(epAM = epAM[i], epfreq = epfreq[i]))

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
response = snow::parLapply(my_cluster, x = ep_params_list, fun = function(ep_params) {
  Species <- "GulfSturgeon"
  controls <- list()
  parameters <- list()
  set.pars(Species)
  MetaPVA(epfreq = ep_params["epfreq"], epAM = ep_params["epAM"])
})

# stop the cluster and timer
snow::stopCluster(my_cluster)
stoptime = Sys.time()
cat("\nElapsed:", format(stoptime - starttime, digits = 2))

# extract the extirpation outcomes at each combo of control parameters
Extir.Res = lapply(1:length(response), function(i) response[[i]]$Extir.Res)
Decline.Res = lapply(1:length(response), function(i) response[[i]]$Decline.Res)
MeanAge.Res = lapply(1:length(response), function(i) response[[i]]$MeanAge.Res)


ProbExt1orMore=matrix(sapply(1:(numepfreqVals * numepAM),function(i)response[[i]]$PropExt.Res[1]),numepAM,numepfreqVals)
ProbExt2orMore=matrix(sapply(1:(numepfreqVals * numepAM),function(i)response[[i]]$PropExt.Res[2]),numepAM,numepfreqVals)
ProbExt3orMore=matrix(sapply(1:(numepfreqVals * numepAM),function(i)response[[i]]$PropExt.Res[3]),numepAM,numepfreqVals)
ProbExt4orMore=matrix(sapply(1:(numepfreqVals * numepAM),function(i)response[[i]]$PropExt.Res[4]),numepAM,numepfreqVals)
ProbExt5orMore=matrix(sapply(1:(numepfreqVals * numepAM),function(i)response[[i]]$PropExt.Res[5]),numepAM,numepfreqVals)
ProbExt6orMore=matrix(sapply(1:(numepfreqVals * numepAM),function(i)response[[i]]$PropExt.Res[6]),numepAM,numepfreqVals)
ProbExtEQ7=matrix(sapply(1:(numepfreqVals * numepAM),function(i)response[[i]]$PropExt.Res[7]),numepAM,numepfreqVals)
ProbExtEQ0=1-ProbExt1orMore

# calculate extirpation probabilities
# ready to be passed to the remainder of the summary prep/plotting code
# NOTE: colMeans is cleaner way to calculate proportion of TRUE elements than colSums/n.sim
extir = sapply(Extir.Res, colMeans)


CatastropheFrequency=matrix(epfreq,numepAM,numepfreqVals)[1,]
CatastropheMortality=seq(.25,.75,length=numepAM)

Pearl_ProbabilityExtripation100Years=matrix(extir[1,],numepAM,numepfreqVals)
Pascagoula_ProbabilityExtripation100Years =matrix(extir[2,],numepAM,numepfreqVals)
Escambia_ProbabilityExtripation100Years =matrix(extir[3,],numepAM,numepfreqVals)
Yellow_ProbabilityExtripation100Years =matrix(extir[4,],numepAM,numepfreqVals)
Choctawhatchee_ProbabilityExtripation100Years =matrix(extir[5,],numepAM,numepfreqVals)
Apalachicola_ProbabilityExtripation100Years =matrix(extir[6,],numepAM,numepfreqVals)
Suwannee_ProbabilityExtripation100Years =matrix(extir[7,],numepAM,numepfreqVals)

contours = list(showlabels = TRUE,start = 0,end = 1,labelfont = list(size = 12, color = "lightgray"))
fig=plot_ly(z=~Pearl_ProbabilityExtripation100Years,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~Pascagoula_ProbabilityExtripation100Years,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~Escambia_ProbabilityExtripation100Years,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~Yellow_ProbabilityExtripation100Years,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~Choctawhatchee_ProbabilityExtripation100Years,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~Apalachicola_ProbabilityExtripation100Years,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~Suwannee_ProbabilityExtripation100Years,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig

fig=plot_ly(z=~ProbExtEQ0,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExt1orMore,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExt2orMore,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExt3orMore,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExt4orMore,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExt5orMore,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExt6orMore,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig
fig=plot_ly(z=~ProbExtEQ7,y=~CatastropheMortality,x=~CatastropheFrequency,type="contour",colors="Greys",contours=contours);fig



# for(i in 1:7){
#   matplot((results[[2]][,i,]),type='l',ylim=c(0,10000))
#   abline(h=50,lwd=2)
# }
# print(colSums(Extir.Results)/n.sim)
# 
# plot(1:200,rowSums(Pops[,,1]),type="l",ylim=c(0,16000),col=1)
# lines(1:200,rowSums(Pops[,,2]),col=2)
# lines(1:200,rowSums(Pops[,,3]),col=3)
# lines(1:200,rowSums(Pops[,,4]),col=7)
# lines(1:200,rowSums(Pops[,,5]),col=5)
# lines(1:200,rowSums(Pops[,,6]),col=6)
# lines(1:200,rowSums(Pops[,,7]),col=4)
# legend("topleft",legend=c("Pearl","Pasc","Esca","Yell","Choc","Apal","Suwa"),lty=1,col=c(1:3,7,5,6,4),bty='n')
# 
# matplot(results[[2]][,1,],type='l',ylim=c(0,10000))
# matplot(results[[2]][,2,],type='l',ylim=c(0,10000))
# matplot(results[[2]][,3,],type='l',ylim=c(0,4000))
# matplot(results[[2]][,4,],type='l',ylim=c(0,10000))
# matplot(results[[2]][,5,],type='l',ylim=c(0,4000))
# matplot(results[[2]][,6,],type='l',ylim=c(0,4000))
# matplot(results[[2]][,7,],type='l',ylim=c(0,10000))


# ProbabilityExtripation100Years=matrix(extir[1,],numepAM,numepfreqVals)
# CatastropheFrequency=matrix(epfreq,numepAM,numepfreqVals)[1,]
# CatastropheMortality=seq(.25,.75,length=numepAM)
# fig1=plot_ly(z=~ProbabilityExtripation100Years,y=~CatastropheMortality,x=~CatastropheFrequency,colors="Spectral")
# fig1 <- fig1 %>% add_surface()
# fig1





