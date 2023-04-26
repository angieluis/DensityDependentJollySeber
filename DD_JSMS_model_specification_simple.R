## Jolly Seber Multistate Model specification  
## Density Dependent K model. See SimulateDensityDependence_simple.R 
## Jolly-Seber models allow for individuals to enter the population
## from a superpopulation - can estimate the prob of entrance
## See DD_JSMS_model_run_simple.R for run code.


## To do: ----------------------------------------------------------##
## Trying Paul's suggestion of the 'ones trick'


## May run faster in NIMBLE instead of JAGS 
## https://r-nimble.org/quick-guide-for-converting-from-jags-or-bugs-to-nimble

## For N/K, do I need to add a 0.001 to the K so it can't be -Inf?

## -----------------------------------------------------------------##

sink("DD_JSMS_model_specification_simple.bug")
cat("
    model {

    # -------------------------------------------------
    # Parameters:
    # phi: survival probability 
    # p: recapture probability 
    # gamma: indiv prob of entry 
    # -------------------------------------------------
    # States (S):
    # 1 not yet entered
    # 2 alive
    # 3 dead
    
    # Observations (O):
    # 1 seen 
    # 2 not seen
    # -------------------------------------------------

          
      ##### PRIORS FOR PHI #####
      m0          ~ dnorm(0, 1)T(-10, 10)    # prior for mortality intercept
      me          ~ dnorm(0, 0.4)T(-10, 10)    # prior for (logit) mortality at equil 

      #### PRIORS for birth ####
      b0          ~ dnorm(0, 0.4)T(-10, 10)   # prior for birth intercept (log at N=2)
      
      #### PRIORS for K  ######
      k.0         ~ dnorm(0, 1)T(-10, 10)    # prior for k coef on intercept for k
      k.ndvi      ~ dnorm(0, 0.4)T(-10, 10)    # prior for k coef on ndvi12
      


	  
      ##### PRIORS FOR RECAPTURE #####
      p     ~ dunif(0, 1)  # prior for p
 

     ##### PRIOR FOR initial population size #####
     N1.est ~ dnorm(20, 5)T(0, 100)   


       
    ##### MODEL FOR PHI #####
        for(m in 1:(n.months - 1)) {              
          
          # calculate K for that month based on ndvi and season
          K[m] <- exp(k.0  + k.ndvi*ndvi[m])

          for(i in 1:n.inds) {
 
            ### Phi 
            logit(mort[i, m]) <-            
              m0 + (me-m0) * N[m]/(K[m]+0.001) 
              
            phi[i, m] <- 1 - mort[i,m]

            } # i for individual 
          } # m for months
      
        
        



 
  ######### Define state-transition and observation matrices --------------#
  
    for (i in 1:n.inds){
      # Define probabilities of State(t+1) given State(t)
      for (m in 1:(n.months-1)){
        ps[1,i,m,1] <- 1-gamma[m]                 # still not yet entered
        ps[1,i,m,2] <- gamma[m]                   # just entered 
        ps[1,i,m,3] <- 0                          # not yet entered to dead
        ps[2,i,m,1] <- 0                          # alive to not yet entered
        ps[2,i,m,2] <- phi[i,m]                   # alive to alive (survival)
        ps[2,i,m,3] <- 1-phi[i,m]                 # alive to dead
        ps[3,i,m,1] <- 0                          # dead to not yet entered
        ps[3,i,m,2] <- 0                          # dead to alive
        ps[3,i,m,3] <- 1                          # dead stay dead
      } #m
     } # i

    
    
    
    # Define probabilities of Observation(m) given State(m)
    # first index is state, last index is observation
    # could potentially include observation as wrong state (false neg or pos)
        po[1,1] <- 0              # not yet entered and observed 
        po[1,2] <- 1              # not yet entered and not observed
        po[2,1] <- p              # alive and observed 
        po[2,2] <- 1-p            # alive and not observed
        po[3,1] <- 0              # dead and observed 
        po[3,2] <- 1              # dead and not observed



    ############### Likelihood ---------------------------------------#
    
      for (i in 1:n.inds){
        
        # STATE PROCESS
        # Define latent state at first dummy occasion
         z[i,1] <- 1   # Make sure that all individuals are in state 1 at t=1 (dummy occasion)
         # No one has entered yet (state 1) at t=1, because of dummy occasion

        for (m in 2:n.months){  #
          # State process: draw S(t) given S(t-1)
          z[i, m] ~ dcat(ps[z[i, m-1], i, m-1, ])


        # OBSERVATION PROCESS
        # draw O(t) given S(t)

          for(d in 1:n.sec.occ){
            y[i, m, d] ~ dcat(po[z[i, m], ])   
          } #d
          
        } #m
      } #i
  

 




    ########################## Calculate derived population parameters
    
    be <- log((exp(me)/(1+exp(me))))   #log(rev.logit(me))
  
            
      for (i in 1:n.inds){
        for (m in 1:n.months){
          alive[i, m] <- equals(z[i, m], 2)
          not.yet.entered[i, m] <- equals(z[i, m], 1)
        } #m
        for (m in 2:n.months){
          just.entered[i, m] <- equals(z[i, m-1]-alive[i, m],0)
        } #m
      } #i

      N[1] <- N1.est                        # number alive now (estimate what it needs to have been for first time step - not actual, but what needs to be to get entries)
      B[1] <- 0                             # Number entries in z from CMR mod
      B.predicted[1] <- 0                   # predicted new this months (based on DD mod)
      Bpossible[1] <- sum(not.yet.entered[ , 1]) # total number that could possibly enter next time (haven't yet entered this time step)
    
      for (m in 2:n.months){
        N[m] <- sum(alive[ , m]) # number alive now
        B[m] <- sum(just.entered[ , m]) # number new entries in z from CMR mod
        f[m-1] <- exp(b0 + (be-b0) * N[m-1]/(K[m-1]+0.001))  # total per capita recruitment rate predicted from DD model
        B.predicted[m]   <-  f[m-1] * N[m-1]         # predicted new recruits = total per capita recruitment rate  *N
        Bpossible[m] <- sum(not.yet.entered[ , m])
        gamma[m-1] <- B.predicted[m]/Bpossible[m-1]
        
        # Want to minimize difference between B.predicted (from N/K model) and B (predicted based on JS, which is new entries to z)
        # define likelihood - here must lie between 0 and 1 
        pr.B[m] <- exp( -0.5*(B.predicted[m] - B[m])^2 )  # this assumes difference fits normal distr with variance = 1, I don't think it matters? because it's hump shaped and max at 0
        ones[m] ~ dbern( pr.B[m] )   #------------------------------------ the 'ones trick' - want to maximize pr.B by fitting to 1
       


      } #m
      
  

}
    ",fill = TRUE)
sink()














