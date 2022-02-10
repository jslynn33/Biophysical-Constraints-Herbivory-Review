#####################################################
### Code for analyses of herbivory across studies ###
### Lynn et al. 2022 Ecography                    ###
### DOI: 10.1111/ecog.06114                       ###
#####################################################

# load required packages
require(rjags);require(R2jags);require(coda)
require(car);library(tidyverse);require(grid)


# load data 
herbdat <- read.csv("LynnETAL_2022_Ecography_sodCaseStudy.csv")

# load in covarience matrix for incorporating phylogeny
source("SpeciesCovMatirx.R") # covariance is "vcvdat"

# logit transform herbivroy data
herbdat$l_herbsite <- logit(herbdat$mean_herb/100, adjust=0.01)
hist(herbdat$l_herbsite)

# functions for standardization 
standard <- function(x) (x - mean(x, na.rm=T)) / (2* (sd(x, na.rm=T)))

# begin data preparation for model

# random effects of 
study <- factor(herbdat$study)

# fixed effects 
sod <- standard(log(herbdat$na_wdcpc)) # natural log transform
lat <- standard(abs(herbdat$latitude)) # take absolute value
mat <- standard(herbdat$mat)
aet <- standard(herbdat$aet)

# phylogentic variance covariance matrix 
spmat <- as.matrix(solve(vcvdat))
# match species names to matrix
spp <- herbdat %>% dplyr::select(species_name) %>% 
  mutate(species_name= str_replace(species_name, " ", "_")) 
spp <- factor(spp$species_name)
Nspp <- as.numeric(length(levels(spp)))

# dependent variable
herb <- herbdat$l_herbsite
N <- as.numeric(length(herb))

# data and parameters to save for the models model
jags.data <- list("herb", "N", "spp", "study", "spmat", "Nspp", "sod")
jags.param <- c("a", "b", "sig1", "sig2","tau","lambda", "sigma", "zlog.lik",
                "rss", "rss.new") 

herbmod <- function() {
  # make a zero for the dmnorm
  for(i in 1:Nspp){zero[i] <- 0} 
  # specification of the random effect for spp with phylo signal
  arter[1:Nspp] ~ dmnorm(zero[],TAU[,])
  
  # study random effect 
  for(j in 1:14){stud[j]~dnorm(0, prec2)}
  for(i in 1:N){
    #likelihood
    herb[i] ~ dnorm(mu[i], prec1)
    mu[i] <- a + b*sod[i] + arter[spp[i]]+ stud[study[i]] # arter is spp random effect
    # 1) form of linear models - just swich sod for mat, aet, and lat
    # 2) drop "b*sod[i]" for null model
    # 3) add "bmat*mat[i] + bna*sod[i]+bnm*sod[i]*mat[i]" for the interaction model
    
    # for model comparison
    zlog.lik[i] <- logdensity.norm(herb[i],mu[i],prec1)
    
    #for model checks
    ssq[i] <- pow((herb[i]-mu[i]),2)
    herb.new[i]~dnorm(mu[i], prec1)
    ssq.new[i] <- pow((herb.new[i]-mu[i]),2)
  }
  
  #priors
  a ~dnorm(0,1.0E-6) # intercept
  b ~dnorm(0,1.0E-6) # slope 
  prec1~dgamma(0.001,0.001)
  sig1 <- 1/sqrt(prec1)
  prec2~dgamma(0.001,0.001)
  sig2 <- 1/sqrt(prec2)
  #derived parameters
  rss <- sum(ssq[])
  rss.new <- sum(ssq.new[])
  
  #Prior for phylo correction
  lambda ~ dunif(0,1) # phylo signal
  tau ~ dgamma(1,1) # prior for precision of phylo autocorrelation
  sigma <- 1/sqrt(tau) # overall SD
  
  for(i in 1:Nspp){
    for(j in 1:Nspp){ # number of speceis
      #LAMBDA=matrix with off-diagonal lambda value and 1 in diagonal
      LAMBDA[i,j] <- 1+(lambda-1)*(1-equals(i,j))
      #sppmat is the inverse variance covariance matrix
      TAU[i,j] <- tau*LAMBDA[i,j]*spmat[i,j]
    }
  }
}

herbres <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                          n.iter=100000,model.file=herbmod, n.thin=5, n.chains=3)

herbres
