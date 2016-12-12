library(tidyverse)
library(tidyverse)
library(jagsUI)

rm(list = ls()) 

# code to run the growth model function is at bottom

# truncated normal dist
rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

runGrowthModel <- function(iterToUse, firstNonBurnIter, chainToUse, simInfo, grDir, out){

  nInd <- simInfo$nInd
  nOcc <- simInfo$nOcc
  nRivers <- simInfo$nRivers
  nSeasons <- simInfo$nSeasons

# get river data
rivIn <- out[chainToUse]$sims.list$riverDATA[ firstNonBurnIter + iterToUse -1,]
# hold in matrix format for growth calculations
r <- matrix(rivIn, nrow=nInd, byrow=TRUE)

#proportion of NAs
pNA <- 0

# initial size data
initSize <- 65
initSizeSD <- 4

# residual growth error
gSD <- 1

# 1st column is mean growth, second col is size-dep coeff 
# rows are seasons
# third array index is river
growth <- 
  array(
    NA, c(nSeasons,2,nRivers)
  )

# setup for 3 rivers
growth[,,1] <- matrix( 
  c( 30, 5, 
     10, 0,
     8, -1,
     4, 0
  ), nrow=nSeasons, byrow=TRUE)
growth[,,2] <- matrix( 
  c( 30, 5, 
     10, 0,
     8, -1,
     4, 0
  ), nrow=nSeasons, byrow=TRUE)
growth[,,3] <- matrix( 
  c( 20, 5, 
     5, 0,
     4, -1,
     2, 0
  ), nrow=nSeasons, byrow=TRUE)


# dat holds body sizes  
dat <- matrix(NA,nInd,nOcc)#, dimnames=list(1:nOcc, 1:T)) 
# get means and sds across occasions for standardizing
meanLen <- matrix(NA,nOcc,nRivers)
sdLen <- matrix(NA,nOcc,nRivers)

dat[,1] <- rtnorm(nInd,initSize,initSizeSD,58,70) # preload occ 1 

# individual random slope multiplier on current body size
# it is proportional to initial size
ranSlope <- (dat[,1] - mean(dat[,1]))/sd(dat[,1]) * 0.75 

# fill matix with body sizes
for (t in 1:(nOcc-1)){
  seas <- (t - 2) %% nSeasons # (t-2) to offset season so occ1=seas3
  if(seas == 0) seas <- nSeasons
  for (i in 1:nInd){    
      dat[i,t+1] <- dat[i,t] + 
                   ( growth[seas,1,r[i,t]] +                                            # intercept   
                     growth[seas,2,r[i,t]] * (dat[i,t] - mean(dat[,t]))/sd(dat[,t]) +   #size-dep growth
                     ranSlope[i] *    (dat[i,t] - mean(dat[,t]))/sd(dat[,t]) ) + #individual random effect on growth
                    rnorm(1,0,gSD)
  }
  meanLen[t,] <- mean(dat[,t])
  sdLen[t,] <- sd(dat[,t]) 
}

# mean and sd for last occ
# placeholder for river-specific means if wanted. Using overall mean.
meanLen[nOcc,] <- mean(dat[,nOcc])
sdLen[nOcc,] <- sd(dat[,nOcc])

d <- gather(as.data.frame(t(dat))) # to long format
names(d) <- c("id","l")
d$occ <- rep(1:nOcc,nInd)
d$seas <- rep( rep((c(3,4,1,2)), nOcc/nSeasons, length.out=nOcc), nInd ) #start with season 3
d$riverDATA <- rivIn
d$ranSlope <- rep(ranSlope,each=nOcc)

d <- d %>% 
  group_by(id) %>%
  mutate( minOcc = min(occ), fOcc = (occ==minOcc)*1,
          maxOcc = max(occ), lOcc = (occ==maxOcc)*1,
          lagLen = lead(l),
          gr = lagLen - l
        )


evalRows <- which( d$lOcc == 0 )
firstObsRows <- which( d$fOcc == 1 )
lastObsRows <- which( d$lOcc == 1 )

nEvalRows <- length(evalRows)
nFirstObsRows <- length(firstObsRows)
nLastObsRows <- length(lastObsRows)

d$lNA <- ifelse(runif(nrow(d),0,1) < pNA & d$fOcc != 1, NA, d$l) #always have a l for fOcc

 ggplot(d, aes(occ,lNA, group=id)) + geom_line(aes(color=factor(riverDATA))) #+ scale_color_gradient2() #+ theme(legend.position="none")

 # mostly a placeholder for first version
 cutoffYOYDATA <- matrix(NA,nSeasons,nRivers) 
 cutoffYOYDATA[,] <- 110 

data <- list( lengthDATA=d$lNA,
              riverDATA=d$riverDATA, 
              nRivers=nRivers, nInd=nInd, nOcc=nOcc, 
              occ = d$occ, season = d$seas,
              nEvalRows=nEvalRows, evalRows = evalRows,
              nFirstObsRows = nFirstObsRows, firstObsRows = firstObsRows,
              nLastObsRows = nLastObsRows, lastObsRows = lastObsRows,
              nSeasons = nSeasons,
              lengthMean = meanLen, lengthSD = sdLen,
              cutoffYOYDATA = cutoffYOYDATA 
            )

inits <- function(){
  list(grBetaInt=array(rnorm(nSeasons*nRivers,0,2.25),c(nSeasons,nRivers)))
}

params <- c("grBetaInt","grBeta","grSigmaBeta")#, "length")

outGR <- jags(data = data,
            inits = inits,
            parameters.to.save = params,
            model.file = paste0(grDir,"grWMoveMod.jags"),
            n.chains = 3,
            n.adapt = 1000, #1000
            n.iter = 2000, # with pNA>0.25, need to run 50,000 iters, otherwise rhats are >1.1
            n.burnin = 1000,
            n.thin = 4,
            parallel = TRUE
           )

  outGR$movementModelIterUsed <- iter
  return(outGR)
}

## Run iterations from the movement model ##

moveDir <- "/home/ben/movementModel/"
grDir <- "/home/ben/growthModel/"
setwd(grDir)

# load output from movement model

load (paste0(moveDir,"simMoveOut.RData"))

chainToUse <- 1
firstNonBurnIter <- (out$mcmc.info$n.burn/out$mcmc.info$n.thin)+1
numItersToUse <- 3
itersToUse <- sort(sample(firstNonBurnIter:(out$mcmc.info$n.iter/out$mcmc.info$n.thin),numItersToUse))

runOverIters <- list()
ii <- 0
for (iter in itersToUse){
  ii <- ii+1
  print(c(ii,iter))
  # saving into a list for now, could also concat a dataframe with identifiers
  runOverIters[[ii]] <- runGrowthModel(iter, firstNonBurnIter, chainToUse, simInfo, grDir out)
}


