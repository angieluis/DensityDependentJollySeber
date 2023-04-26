#############################################################################
## Simulate Data and explore the model for density-dependent IPM
##  Robust design model, where environmental drivers (NDVI)
## affect K and that drives survival and birth rates. 
#############################################################################


######################################################
## Set parameter Values for simulation
######################################################

## Mortality parameters
m0 <- -2.944439 #logit(0.05)  #survival is 0.95 at N=0
me <- -0.8472979 #logit(0.3) # survival at equilibrium is 0.7


# Recruitment parameters
b0 <- 0.6931472 #log(2) # birth rate at N=0 is 2
# be <- log((exp(me)/(1+exp(me))))   #log(rev.logit(me))
# birth rate at equilibrium is 0.3 = death rate

# K parameters
k.0 <- 3.5        # K is exp(k.0) at mean ndvi  
k.ndvi <- 0.35

# detection parameters
p <- 0.58 



#############################################################################
############## Explore relationships
#############################################################################

NK<-seq(0,2,length=50) # N/K, so equilibrium is at 1.

rev.logit <- function(x) {
  exp(x) / (1 + exp(x))
}

ddbirthfun <- function(NK, b0, be){
  exp(b0 + (be-b0)*NK)
}

ddmortfun <- function(NK,m0, me){
  rev.logit(m0 + (me-m0)*NK)
}

birth.rate <- exp(b0 + (be-b0)*NK)
mortality.rate  <- rev.logit(m0 + (me-m0)*NK)

plot(NK,birth.rate,type="l",col="slateblue",ylim=c(0,1.5),
     ylab="Rates",xlab="N/K",main="Density-Dependent Rates")
lines(NK,mortality.rate,col="tomato")
abline(v=1,lty=3)
text(1,1,"equilibrium")
legend("topright",c("recruitment rate","mortality rate"),lty=1,
       col=c("slateblue","tomato"),bty="n",cex=0.8)


### Make K time varying as a function of ndvi --------------------#

K <- exp(k.0  + k.ndvi*seq(-2,2,by=0.1))
plot(seq(-2,2,by=0.1),K,type="l",col="blue",main="carrying capacity by ndvi",xlab="ndvi") 


##########################################################################
## Simulate data 
# Use site grandcanyon.e as a template
# so use covariate data for that sites for the first n.months
##########################################################################


source("RealData/Code/01_sw_data_functions_more.R")
temp.data <- read.csv("RealData/Data/southwest_covariates_norm.csv")
temp.data$date2 <- lubridate::dmy(paste("1", temp.data$date))
temp.data$site.web <- paste(temp.data$site,temp.data$web,sep=".")
web <- "grandcanyon.e"
n.start <- 30 ## numebr of individuals to start simulation with at first time step
n.months <- 40   # months to simulate
start.date <- temp.data$date2[20] # start covariate data at Aug 1994
end.date <- start.date + months(n.months-1)
n.sec.occ <- 3 # secondary occasions
prim.dates <- seq(start.date,end.date,by="month")
sessions <- character()
for(i in 1:length(prim.dates)){
  sessions[i] <- ifelse(month(prim.dates[i])<10 ,paste(year(prim.dates[i]),"0",month(prim.dates[i]),sep=""),paste(year(prim.dates[i]),month(prim.dates[i]),sep=""))
}
## Cut down temporal covariate data to just this web and dates
temp.data <- temp.data[which(temp.data$site.web==web),]
s.ind <- which(temp.data$date2==start.date) 
temp.data <- temp.data[s.ind:(s.ind+(n.months-1)),] 

plot.ts(exp(k.0  + k.ndvi*temp.data$ndvi), ylab="K",main="K(t) driven by NDVI(t)",
        col="blue")



## Simulate data ------------------------------------------------------------#


## Simulate Z -----------------------------------------------------#
z <- matrix(NA,n.start,n.months)
zsex <- rep(0:1,n.start/2) # start each web with 1/2 males and 1/2 females
z[,1] <- 1 # all n.start are alive at time 1
N <- rep(NA,n.months)
K <- rep(NA,n.months)
Phi.f <-rep(NA,n.months) 
  
for(m in 2: n.months){
    
    # calculate K for that month based on ndvi and season
    K[m-1] <- exp(k.0  + k.ndvi*temp.data$ndvi[m-1])
    
    # calculate population size 
    N[m-1] <- sum(length(which(z[,m-1]>0)) )
    
    
    for(i in 1:dim(z)[1]) {
      
      # calculate survival based on N/K (last month)
      mort <- rev.logit(m0 + (me-m0)* N[m-1]/ K[m-1])
      phi <- 1 - mort
      
      last.state <- z[i, m-1]
      
      #if last.state dead
      if(last.state==0){
        z[i, m] <- 0  #stay dead
      }
      
      # if last.state alive
      if(last.state==1){
        # stay alive (1) with prob phi
        # die (0) with prob 1-phi
        
        z[i,m] <- rbinom(1, 1, phi)
      }
      
     } #i
    
    
    # keep track of Phis
    Phi.f[m-1] <- phi
 
    
    ########## Add individuals (rows) to z based on births
  
    be <- log((exp(me)/(1+exp(me))))   #log(rev.logit(me))
    birth.rate <- exp(b0 + (be-b0) * N[m-1]/ K[m-1])
    new.indivs <- rpois(1, birth.rate*sum(length(which(z[,m]>0)))) 
    z.new <- matrix(rep(c(rep(0, m-1), # before now, the animal was not here
                          1, # now at month m, animal is here
                          rep(NA, n.months-m)), # rest of row is NA
                        new.indivs), byrow=TRUE, nrow=new.indivs, ncol=n.months)
    z <- rbind(z, z.new)
    # assume equal chance of males and females
    zsex <- c(zsex, sample(c(0,1),size=new.indivs,replace=TRUE)) 
    
    
    ############## Save individual data
      indiv.data <- data.frame(ID = 1:dim(z)[1], site.web = web, 
                                    sex = zsex)
      
} #m
  


# take a look at simulated numbers
colSums(z)





# Simulate observations in a robust design array format ------------------#

y <- array(NA, dim = c(dim(z),n.sec.occ))

for(i in 1:dim(y)[1]){
  for(m in 1:dim(y)[2]){
    y[i, m, ] <- rbinom(n.sec.occ, 1, p) * z[i, m]
  }
}

# Remove rows from observation data for which animals were never observed ---#
rm.ind <- which(apply(y, 1, sum)==0)
y <- y[-rm.ind, , ]
indiv.data <- indiv.data[-rm.ind, ]


save(y, temp.data, indiv.data, file="SimDensDepData_simple.RData")

