## Jolly Seber Multistate Model run code, using 
## Density Dependent K model - see "SimulateDensityDependence_simple.R"
## Jolly-Seber models allow for individuals to enter the population
## from an augmented population - can estimate the prob of entrance & pop size


## To Do-------------------------------------------------------------------

## error with the ones trick. indexing?

## ------------------------------------------------------------------------


load("SimDensDepData_simple.RData")
source("DD_JSMS_model_specification_simple.R")
#source("~/Documents/JAGS/BayesianMarkRecapSNV/RealData/Code/01_sw_data_functions_more.R")

library(R2jags)
library(abind)
## for this version, the robust design capture histories are in an  array with dimensions
# y[i, m, d] #individual, month (primary occasion), day (secondary occasion)


## need to augment the capture histories to allow for more individuals than 
##   those detected
n.aug <- 200 # number of individuals added to capture histories (not observed)
y.aug <- abind(y, array(0, dim = c(n.aug, dim(y)[2], dim(y)[3])), along = 1)

## Add a dummy occasion to beginning when no one has entered yet
y.aug.du <- abind(array(0, dim = c(dim(y.aug)[1], 1, dim(y.aug)[3])), y.aug, along = 2)

## collapse to a monthly (primary) capture history for known state and init functions
## this one needs dummy occasion but not augmented
y.primary <- apply(y, 1:2, sum)
y.primary <- replace(y.primary, y.primary>1, 1)
y.primary.du <- cbind(0, y.primary) # add dummy occasion
y.primary.aug.du <- rbind(y.primary.du, matrix(0, nrow = n.aug, ncol=dim(y.primary.du)[2]))# augemented with dummy occasion

# function to create matrix with info about known latent state z for Jolly-Seber
# models 
known.state.js <- function(ch) {
  state <- ch
  for (i in 1:dim(ch)[1]) {
    n1 <- min(which(ch[i, ] == 1))
    n2 <- max(which(ch[i, ] == 1))
    state[i, n1:n2] <- 1
    # only filling in those that were 0s but we know were alive
    # because caught before and after
   }
  state[state == 0] <- NA
  return(state)
}
knownz <- known.state.js(y.primary.du) 
## this has 1=alive. 
## But in the model, state 1=not yet entered, 2=alive, 3=dead
## so add 1
knownz <- knownz + 1
#add the augmented rows 
knownz <- rbind(knownz, matrix(NA, nrow=n.aug, ncol=dim(knownz)[2]))



##### Function to give initial values for z ----------------------------#

################################################################## 
## function to specify initial values for Jolly-Seber multistate
# 1 not yet entered
# 2 alive 
# 3 dead
# this just puts in an initial value of alive (2) if ever seen
# or 1 if not yet seen & has NA for dummy occasion & known state
##################################################################
MSJS.init.z <- function(ch=y.primary.aug.du, knownz){ # assumes ch is a array [i,m,d] of primary occasions already augmented with dummy occasion if appropriate; 
  # knownz is a matrix [i,m] - assumes dummy occasion for first time point, and augmented, if appropriate
  initz <- matrix(NA, ncol = dim(knownz)[2], nrow=dim(knownz)[1]) 
  suppressWarnings(f <- apply(ch,1,function(x){min(which(x > 0))}))
  ever.f <- which(is.finite(f))
  state <- matrix(1, nrow = dim(knownz)[1], ncol = dim(knownz)[2]) # default value of not yet entered
  for(i in ever.f){
    state[i,f[ever.f[i]]:dim(state)[2]] <- 2 # fill in 2 for alive after enter  
  }
  # remove those that are in the known state
  state <- replace(state,!is.na(knownz),NA)
  initz <- state
  initz[, 1] <- NA # dummy occasion is in likelihood
  return(initz)
}


#initz <- MSJS.init.z(ch=y.primary.aug.du, knownz)


jags.data <- list(
  n.inds =  dim(y.aug.du)[1],
  n.months = 41, # total months in z : length w +1 for dummy occasion
  n.sec.occ = 3, # here simplified but could be dimensions:[w, Prim[i,m]]
  ndvi = temp.data$ndvi, 
  #sex = indiv.data$sex, # mat with dimensions: [i,w]
  y = replace(y.aug.du,y.aug.du==0,2), # replace 0's with 2's 
  primary.ch = y.primary.aug.du, # still with 0's for init function
  ones = rep(1,  dim(y.aug.du)[2]), #[w,m] for the 'ones' trick
  z = knownz, # for known states
  MSJS.init.z = MSJS.init.z
  
)

inits <- function(){list(m0 = runif(1, 0, 1), 
                         me = runif(1, 0, 1), 
                         b0 = runif(1, 0, 1), 
                         k.0 = runif(1, 0, 1), 
                         k.ndvi = runif(1, 0, 1), 
                         p = runif(1, 0, 1),
                         N1.est = rnorm(1, 20, 2), 
                         z = MSJS.init.z(primary.ch, z)
                           )}


parameters <- c("m0",  #<- -2.944439 #logit(0.05)  #survival is 0.95 at N=0
              "me", # <- -0.8472979 #logit(0.3) # survival at equilibrium is 0.7
              "b0", # <- 0.6931472 #log(2) # birth rate at N=0 is 2
              "be", # <- log((exp(me)/(1+exp(me))))   #log(rev.logit(me))
              "k.0", # <- 3.8        # K is exp(k.0) at mean ndvi in winter 
              "k.ndvi", # <- 0.35
              "N1.est", # 
              "p", #            <-  0.58 
              "N",
              "K")



date() #date and time before run
DD.JSMS.simulation.output <- jags.parallel(data     = jags.data,
                                            inits,
                                            parameters,
                                            "DD_JSMS_model_specification_simple.bug",
                                            n.chains = 3,
                                            n.thin   = 3,
                                            n.iter   = 20000,
                                            n.burnin = 5000)

date() #date and time after run
save(DD.JSMS.simulation.output,file="DDJSMSsimoutput.RData")

# took 41 hours to say:
# Error in checkForRemoteErrors(val) : 
#   2 nodes produced errors; first error: Error in node z[1,18]
# Node inconsistent with parents

# had previously run 1 chain 1 iteration with no errors