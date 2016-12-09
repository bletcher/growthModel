library(tidyverse)
library(jagsUI)

nInd <- 100
nOcc <- 10
nRivers <- 1
nSeasons <- 4

#proportion of NAs
pNA <- 0

# initial size SD
initSizeSD <- 0.025

# residual growth error
gSD <- 0.1

growth <- array(NA,c(nSeasons,2))

growth <- matrix( 
  c( 0.2, 0, 
     1.0, 0,
     0.4, 0,
     0.1, 0
  ), nrow=nSeasons, byrow=TRUE)
  
  
dat <- matrix(NA,nInd,nOcc)#, dimnames=list(1:nOcc, 1:T)) 

for (i in 1:nInd){
  dat[i,1] <- rnorm(1,-2,initSizeSD) 
  for (t in 1:(nOcc-1)){

      seas <- t %% nSeasons
      if(seas == 0) seas <- nSeasons

      dat[i,t+1] <- dat[i,t] + 
                   (growth[seas,1] + growth[seas,2]*dat[i,t]) +
                    rnorm(1,0,gSD)
  }
}

d <- gather(as.data.frame(t(dat))) # to long format
names(d) <- c("id","l")
d$occ <- rep(1:nOcc,nInd)
d$seas <- rep( rep((1:4), nOcc/nSeasons, length.out=nOcc), nInd )
d$riverDATA <- 1

d <- d %>% 
  group_by(id) %>%
  mutate( minOcc = min(occ), fOcc = (occ==minOcc)*1,
          maxOcc = max(occ), lOcc = (occ==maxOcc)*1
        )
evalRows <- which( d$lOcc == 0 )
firstObsRows <- which( d$fOcc == 1 )

nEvalRows <- length(evalRows)
nFirstObsRows <- length(firstObsRows)

d$lNA <- ifelse(runif(nrow(d),0,1) < pNA & d$fOcc != 1, NA, d$l) #always have a r for fOcc

#ggplot(d, aes(occ,l)) + geom_line(aes(color=id))

sink("gr.jags")
cat("
    model {
    
    ############################
    ###### growth model ########
    ############################

    for( i in 1:nFirstObsRows ){
      length[ firstObsRows[i] ] <- lengthDATA[ firstObsRows[i] ]
    }
    # for( i in 1:nFirstObsRows ){
    #   length[ firstObsRows[i] ] ~ dnorm( -2,0.001 )   
    #   lengthDATA[ firstObsRows[i] ] ~ dnorm( length[ firstObsRows[i] ],9 ) # length measurement error
    # }

    for( i in 1:nEvalRows ){
      
      length[ evalRows[i]+1 ] <- length[ evalRows[i] ] + gr[ evalRows[i] ] 
      
      lengthDATA[ evalRows[i] + 1 ] ~ dnorm( length[ evalRows[i] + 1 ], 9/1 ) # length measurement error 
      
      gr[ evalRows[i] ] ~ dnorm( expectedGR[ evalRows[i] ],#  0.001)   
                               #  1/expectedGRSigma[ evalRows[i] ]^2 )
                              #   expectedGRSigma[ evalRows[i] ] 
                               #  grSigmaBeta[ season[ evalRows[i] ],riverDATA[ evalRows[i] ] ]
                              #    tau
                               0.001
                               )
  
      expectedGR[ evalRows[i] ]  <- 
        grBetaInt[ season[ evalRows[i] ],riverDATA[ evalRows[i] ] ] 
  #    + grBeta[ 1, season[ evalRows[i] ],riverDATA[ evalRows[i] ] ] * length[ evalRows[i] ] 

   #   expectedGRSigma[ evalRows[i] ] ~ dnorm(0,0.001)T(0,)
   #   expectedGRSigma[ evalRows[i] ] <- tau #grSigmaBeta[ season[ evalRows[i] ],riverDATA[ evalRows[i] ] ]   
    }

  ############## Growth Priors ###############
  tau <- pow(sigma, -2)
  sigma ~ dunif(0,100) #dgamma(0.001,0.001)

    for( r in 1:( nRivers ) ){

#      grBetaInt[ 1,r ] ~ dnorm( 0.2, 0.01 )
#      grBetaInt[ 2,r ] ~ dnorm( 1, 0.01 )
#      grBetaInt[ 3,r ] ~ dnorm( 0.4, 0.01 )
#      grBetaInt[ 4,r ] ~ dnorm( 0.1, 0.01 )

      for( s in 1:4 ){   

  #      grSigmaBeta[ s,r ] ~ dgamma(0.001,0.001)#dnorm( 0,0.01 )T(0,)    
         grBetaInt[ s,r ] ~ dnorm( 0,0.0001 )
  #      grBeta[ 1,s,r ] ~ dnorm( 0,0.0001 )
      }
    }
 #   grSigmaBeta ~ dgamma(0.001,0.001) #dnorm( 0,0.001 )T(0,)    
  
  } # model
    ",fill = TRUE)
sink()

data <- list( lengthDATA=d$lNA,
              riverDATA=d$riverDATA, 
              nRivers=nRivers, nInd=nInd, nOcc=nOcc, 
              occ = d$occ, season = d$seas,
              nEvalRows=nEvalRows, evalRows = evalRows,
              nFirstObsRows = nFirstObsRows, firstObsRows = firstObsRows,
              nSeasons = nSeasons)

inits <- function(){
  list(grBetaInt=array(rnorm(nSeasons*nRivers,0,2.25),c(nSeasons,nRivers)))
}

params <- c("grBetaInt","grBeta","grSigmaBeta","tau","sigma","length")#, "riverDATA")

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


