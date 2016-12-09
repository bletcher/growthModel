library(tidyverse)
library(jagsUI)

nInd <- 100
nOcc <- 10
nRivers <- 1
nSeasons <- 4

#proportion of NAs
pNA <- 0

# initial size data
initSize <- 65
initSizeSD <- 4

# residual growth error
gSD <- 2

growth <- array(NA,c(nSeasons,2))

# 1st column is mean growth, second col is size-dep coeff
growth <- matrix( 
  c( 25, 5, 
     10, 0,
     8, -2,
     4, -1
  ), nrow=nSeasons, byrow=TRUE)
  
# dat holds body sizes  
dat <- matrix(NA,nInd,nOcc)#, dimnames=list(1:nOcc, 1:T)) 
# get means and sds across occasions for standardizing
meanLen <- matrix(NA,nOcc,nRivers+1)
sdLen <- matrix(NA,nOcc,nRivers+1)

for (i in 1:nInd){ dat[i,1] <- rnorm(1,initSize,initSizeSD) } # preload occ 1 

for (t in 1:(nOcc-1)){
  seas <- (t - 2) %% nSeasons # (t-2) to offset season so occ1=seas3
  if(seas == 0) seas <- nSeasons
  for (i in 1:nInd){    
      dat[i,t+1] <- dat[i,t] + 
                   ( growth[seas,1] + growth[seas,2]*(dat[i,t] - mean(dat[,t]))/sd(dat[,t]) ) +
                    rnorm(1,0,gSD)
  }
  meanLen[t,1] <- mean(dat[,t])
  sdLen[t,1] <- sd(dat[,t]) 
}

# mean and sd for last occ
meanLen[nOcc,1] <- mean(dat[,nOcc])
sdLen[nOcc,1] <- sd(dat[,nOcc])

d <- gather(as.data.frame(t(dat))) # to long format
names(d) <- c("id","l")
d$occ <- rep(1:nOcc,nInd)
d$seas <- rep( rep((c(3,4,1,2)), nOcc/nSeasons, length.out=nOcc), nInd ) #start with season 3
d$riverDATA <- 1

d <- d %>% 
  group_by(id) %>%
  mutate( minOcc = min(occ), fOcc = (occ==minOcc)*1,
          maxOcc = max(occ), lOcc = (occ==maxOcc)*1
        )
evalRows <- which( d$lOcc == 0 )
firstObsRows <- which( d$fOcc == 1 )
lastObsRows <- which( d$lOcc == 1 )

nEvalRows <- length(evalRows)
nFirstObsRows <- length(firstObsRows)
nLastObsRows <- length(lastObsRows)

d$lNA <- ifelse(runif(nrow(d),0,1) < pNA & d$fOcc != 1, NA, d$l) #always have a r for fOcc

 ggplot(d, aes(occ,l)) + geom_line(aes(color=id))

sink("gr.jags")
cat("
    model {
    
    ############################
    # Variable standardization #
    ############################
  
    ############ standardized length for sizeBetas
    for( i in 1:(nEvalRows) ){ 
      stdLength[ evalRows[i] ] <-  ( length[ evalRows[i] ] - lengthMean[ season[ evalRows[i]],riverDATA[ evalRows[i] ] ] ) / 
        lengthSD[ season[ evalRows[i]],riverDATA[ evalRows[i] ] ]
  
    }

    for( i in 1:( nLastObsRows ) ){  
      stdLength[ lastObsRows[ i ] ] <-  ( length[ lastObsRows[i] ] - lengthMean[ season[ lastObsRows[i]],riverDATA[ lastObsRows[i] ] ] ) / 
        lengthSD[ season[ lastObsRows[i]],riverDATA[ lastObsRows[i] ] ]
    }



    ############################
    ###### growth model ########
    ############################

#    for( i in 1:nFirstObsRows ){
#      length[ firstObsRows[i] ] <- lengthDATA[ firstObsRows[i] ]
#    }
    for( i in 1:nFirstObsRows ){
      length[ firstObsRows[i] ] ~ dnorm( 65,0.001 )   
      lengthDATA[ firstObsRows[i] ] ~ dnorm( length[ firstObsRows[i] ],9 ) # length measurement error
    }

    for( i in 1:nEvalRows ){
      
      length[ evalRows[i]+1 ] <- length[ evalRows[i] ] + gr[ evalRows[i] ] 
      
      lengthDATA[ evalRows[i] + 1 ] ~ dnorm( length[ evalRows[i] + 1 ], 9/1 ) # length measurement error 
      
      gr[ evalRows[i] ] ~ dnorm( expectedGR[ evalRows[i] ], 
                                 1/expectedGRSigma[ evalRows[i] ]^2 
                               )
  
      expectedGR[ evalRows[i] ]  <- 
        grBetaInt[ season[ evalRows[i] ],riverDATA[ evalRows[i] ] ] 
      + grBeta[ 1, season[ evalRows[i] ],riverDATA[ evalRows[i] ] ] * stdLength[ evalRows[ i ] ]

      expectedGRSigma[ evalRows[i] ] <- grSigmaBeta[ season[ evalRows[i] ],riverDATA[ evalRows[i] ] ]   
    }

  ############## Growth Priors ###############

    for( r in 1:( nRivers ) ){

#      grBetaInt[ 1,r ] ~ dnorm( 8, 0.01 )
#      grBetaInt[ 2,r ] ~ dnorm( 4, 0.01 )
#      grBetaInt[ 3,r ] ~ dnorm( 25, 0.01 )
#      grBetaInt[ 4,r ] ~ dnorm( 10, 0.01 )

      for( s in 1:4 ){   

        grSigmaBeta[ s,r ] ~ dgamma(0.001,0.001)   
        grBetaInt[ s,r ] ~ dnorm( 0,0.0001 )
        grBeta[ 1,s,r ] ~ dnorm( 0,0.0001 )
      }
    }
  
  } # model
    ",fill = TRUE)
sink()

data <- list( lengthDATA=d$lNA,
              riverDATA=d$riverDATA, 
              nRivers=nRivers, nInd=nInd, nOcc=nOcc, 
              occ = d$occ, season = d$seas,
              nEvalRows=nEvalRows, evalRows = evalRows,
              nFirstObsRows = nFirstObsRows, firstObsRows = firstObsRows,
              nLastObsRows = nLastObsRows, lastObsRows = lastObsRows,
              nSeasons = nSeasons,
              lengthMean = meanLen, lengthSD = sdLen
            )

inits <- function(){
  list(grBetaInt=array(rnorm(nSeasons*nRivers,0,2.25),c(nSeasons,nRivers)))
}

params <- c("grBetaInt","grBeta","grSigmaBeta")#, "riverDATA")

out <- jags(data = data,
            inits = inits,
            parameters.to.save = params,
            model.file = "gr.jags",
            n.chains = 3,
            n.adapt = 100,
            n.iter = 1000,
            n.burnin = 500,
            n.thin = 2,
            parallel = TRUE
           )

