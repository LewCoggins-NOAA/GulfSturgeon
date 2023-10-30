
# set number of cores to use
# I have 12, this allows doing *some* other stuff 
# more cores, more efficiency, so long as you have enough RAM too
ncpu = 10

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
  species <- "GulfSturgeon"
  controls <- list()
  parameters <- list()
  set.pars(species)
  MetaPVA(epfreq = ep_params["epfreq"], epAM = ep_params["epAM"])
})

# stop the cluster and timer
snow::stopCluster(my_cluster)
stoptime = Sys.time()
cat("\nElapsed:", format(stoptime - starttime, digits = 2))

# extract the extirpation outcomes at each combo of control parameters
Extir.Res = lapply(1:length(response), function(i) response[[i]]$Extir.Res)

# calculate extirpation probabilities
# ready to be passed to the remainder of the summary prep/plotting code
# NOTE: colMeans is cleaner way to calculate proportion of TRUE elements than colSums/n.sim
extir = sapply(Extir.Res, colMeans)
