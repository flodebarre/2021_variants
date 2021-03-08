#rm(list = ls())
library(MASS)
library(plyr)
library(VGAM)
library(plyr)

#load("inference.RData")
source("functions_v2.R")

# parameters and priors common to all:
gtpar <- c(7.8, 0.833); meangt <- 6.5 # generation time parameter mean 6.5, sd 2.3
gtpar <- c(6, 0.833); meangt <- 5 # mean 5, sd 2

solve_rtoR(param = gtpar, r = 0.0775)
upper_limit_rho <- 0.001 #1/(getbetaparams(0.03, 0.02^2)$alpha + getbetaparams(0.03, 0.02^2)$beta + 1) # motivated by the spatial heterogeneity at first point
sd_normal <- 0.3
niter <- 100000

# CREATE A FAST FUNCTION TO CONVERT FROM R to r
Rseq <- seq(0.5, 6, 0.1)
plot(rseq <- sapply(Rseq, function(myR) solve_Rtor(param = gtpar, myR)), Rseq, type = "l")
ss_rtoR <- smooth.spline(rseq, Rseq)
ss_Rtor <- smooth.spline(Rseq, rseq)
lines(ss_rtoR, col = "red")
fast_solve_rtoR <- function(r) if(all(r > min(rseq) & r < max(rseq))) predict(ss_rtoR, r)$y else stop("r value not in range")
fast_solve_Rtor <- function(R) if(all(R >= 0.5 & R <= 6)) predict(ss_Rtor, R)$y else {print(R); stop("R value not in range")}

# FD Compare formulas
solve_Rtor(param = gtpar, R = 2)
new_Rtor(param = gtpar, R = 2)
solve_rtoR(param = gtpar, r = 0.09)
new_rtoR(param = gtpar, r = 0.09)

new_fast_solve_rtoR <- function(r){new_rtoR(param = gtpar, r = r)}
new_fast_solve_Rtor <- function(R){new_Rtor(param = gtpar, R = R)}
# endFD

regional_lik <- function(par, region, epitable, freqdatatable, return_full = FALSE){

  stopifnot(length(par)==8)
  # region specific parameters:
  p0 <- par[1]
  Rinit <- par[2]
  Rfinal <- par[3]
  k <- par[4] # > 0, not too large
  inflection_date <- par[5] # must be in minday:maxday
  
  # generic parameters
  transmission_advantage <- par[6]
  size_binomial <- par[7] # dispersion of negative binomial
  rho <- par[8] # dispersion of beta binomial
  
  RWT <- Rinit + (Rfinal - Rinit) / (1 + exp(- k * (minday:(maxday-1) - inflection_date)))
  RVOC <-  RWT * (1 + transmission_advantage)
  
  rWT <- new_fast_solve_Rtor(RWT) # FD
  rVOC <- new_fast_solve_Rtor(RVOC) # FD
  
  stopifnot(length(rWT) == n-1)
  stopifnot(length(RVOC) == n-1)
  
  epitable$PVOCsim <- NA
  epitable$PWTsim <- NA

  # initialise prevalence of VOC and WT using frequency parameter
  epitable[1, "PVOCsim"] <- p0 * epitable[1, "Psmoothed"] 
  epitable[1, "PWTsim"] <- (1-p0) * epitable[1, "Psmoothed"]
  
  # extrapolate with exp growth from minday+1 to maxday
  if(!return_full){
    epitable[2:n, "PWTsim"]  <- epitable[1 , "PWTsim"] * exp(cumsum(rWT))
    epitable[2:n, "PVOCsim"] <- epitable[1 , "PVOCsim"] * exp(cumsum(rVOC))
  } else {
    # add rs values equal to the last one:
    finalrWT <- rWT[n-1]
    finalrVOC <- rVOC[n-1]
    rWT <- c(rWT, rep(finalrWT, 90))
    rVOC <- c(rVOC, rep(finalrVOC, 90))
    epitable <- rbind.fill(epitable, expand.grid(reg2 = region, dateday = (maxday+1):(maxday+90)))
    epitable[2:nrow(epitable), "PWTsim"]  <- epitable[1 , "PWTsim"] * exp(cumsum(rWT))
    epitable[2:nrow(epitable), "PVOCsim"] <- epitable[1 , "PVOCsim"] * exp(cumsum(rVOC))
  }
  
  epitable$PTOTsim <- epitable$PVOCsim+epitable$PWTsim
  
  # component of the likelihood based on number of Psmoothed daily compared to inital prediction
  ll0 <- - sum(dnbinom(x = epitable$Psmoothed, size = size_binomial, mu = epitable$PTOTsim, log = T), na.rm = T)
  
  # LIKELIHOOD COMPONENT CORRESPONDING TO FREQUENCIES
  
  if(nrow(freqdatatable) > 0){ # if frequency data is available
    lik_freqs <- sum(
      apply(freqdatatable, 1, function(vec) getlik_frequency(NVOCdata = as.numeric(vec["NVOC"]), Ndata = as.numeric(vec["N"]), daymin = as.numeric(vec["daymin"]), daymax = as.numeric(vec["daymax"]), simtable = epitable, myrho = rho))
    )
  }

  if(return_full){
    return(list(ll = ll0+sum(lik_freqs), epitable = epitable))
  } else {
    return(ll0+sum(lik_freqs))
  }
}

########################################### TABLE OF FREQUENCY DATA ###########################################

# Flash data
freqs <- freqs[-match(c("N_SNEG","N_SEQ","NSEQ_OK","NVOC", "f_neg", "fVOC_in_neg"), names(freqs))]
freqs$daymin <- minday
freqs$daymax <- minday + 1

# Flash data2 from 27/01
freqs2$daymin <- 393
freqs2$daymax <- 393

# now add some data
# PITIE IDF 11-21 JAN
# https://www.nouvelobs.com/coronavirus-de-wuhan/20210126.OBS39397/le-variant-anglais-represente-desormais-environ-10-des-cas-en-ile-de-france-s-inquiete-l-ap-hp.html
freqs3 <- data.frame(region = "ILE DE FRANCE", N = 1080, fVOC = 0.094, daymin = 377, daymax = 387)

freqs4 <- rbind(freqs, freqs2, freqs3)

for(mycol in c("N", "fVOC", "daymin", "daymax")) freqs4[, mycol] <- as.numeric(freqs4[, mycol])
freqs4$NVOC <- round(freqs4$N * freqs4$fVOC)
freqs4$lower <- freqs4$fVOC - 1.96 * sqrt(freqs4$fVOC * (1-freqs4$fVOC) / freqs4$N)
freqs4$upper <- freqs4$fVOC + 1.96 * sqrt(freqs4$fVOC * (1-freqs4$fVOC) / freqs4$N)
freqs4$daymid <- 0.5 * (freqs4$daymin + freqs4$daymax)

###########################################  OPTIMISE ONE REGION  ###########################################

## OPTIMISE ONE REGION ONLY
myregion <- "ILE DE FRANCE"
regional_lik(
  par = c(
    0.069, 1.1, 0.9, 0.1, 390, # region-specific pars
    0.5, 10, 0.1 # generic pars
  ),
  region = myregion,
  epitable = a2[which(a2$reg2 == myregion), ],
  freqdatatable = freqs4[which(freqs4$region == myregion), ]
)
# init frequency; Rinit; Rfin; inflection speed; inflection mid-point; selection advantage; dispersion negbin
lower0 <- c(0, 0.5, 0.5, 0.001, minday, 0, 0.001, 0)
upper0 <- c(0.1, 2, 2, 10     , maxday, 1, 1000, upper_limit_rho)
npar0  <- 8
init_fun0 <- function(x) sapply(1:npar0, function(i) runif(1, lower0[i], upper0[i]))
opt0 <- optim.fun.repeated(n.repeats = 10, lik.fun = regional_lik, init.fun = init_fun0, lower = lower0, upper = upper0,
                           region = myregion,
                           epitable = a2[which(a2$reg2 == myregion), ],
                           freqdatatable = freqs4[which(freqs4$region == myregion), ])

betaprior <- getbetaparams(mu = 0.033, var = 0.02^2) # prior beta with mean 0.033 and sd 0.02
proposal.sd0 <- c(0.1, 0.01, 0.01, 1, 0.01, 0.05, 1, 0.1)
prior0 <- function(par){
  sd_normal <- 0.3
  return(
    (dbeta(x = par[1], shape1 = betaprior$alpha, shape2 = betaprior$beta, log = T)) + # beta prior for VOC frequency
      (dnorm(x = par[2], mean = 1, sd = sd_normal, log = T)) + # normal prior for R wild type init
      (dnorm(x = par[3], mean = 1, sd = sd_normal, log = T))   # normal prior for R wild type final
    # uniform prior fork, inflection dates, transmission advantage, size of the binomial, dispersion beta biniomial (do not code them)
  )
}


posterior_idf <- function(par) prior0(par) - regional_lik(par, region = myregion, epitable = a2[which(a2$reg2 == myregion), ], freqdatatable = freqs4[which(freqs4$region == myregion), ])
mcmc_idf <- run_MCMC_Gibbs(start.value = opt0$pars, N.iter = niter, proposal.sd = proposal.sd0, posterior.fun = posterior_idf, npar = npar0, lower = lower0, upper = upper0)

########################################################################### OPTIMISE ALL REGIONS TOGETHER  ###########################################################################

epitable0 <- ddply(a2, .(dateday), summarise, Psmoothed = sum(Psmoothed), reg2="ALL")
freqdatatable0 <- ddply(freqs4, .(daymin, daymax), summarise, N = sum(N), NVOC = sum(NVOC), region ="ALL") 
freqdatatable0 <- freqdatatable0[-2, ]# remove row 2 which corresponds to IDF only
freqdatatable0 <- rbind(freqdatatable0, c(408, 408, 70498, 26063, "ALL")) # ADD 37% MENTIONED IN SPF REPORT 18/02/2021, ANNEX 19/02/2020 (405-411 or 408)
freqdatatable0$N <- as.numeric(freqdatatable0$N)
freqdatatable0$NVOC <- as.numeric(freqdatatable0$NVOC)
epitable0$Psmoothed <- as.numeric(epitable0$Psmoothed)

posterior0 <- function(par) prior0(par) - regional_lik(par, region = "all", epitable = epitable0, freqdatatable = freqdatatable0)
mcmc0 <- run_MCMC_Gibbs(start.value = opt0$pars, N.iter = niter, proposal.sd = proposal.sd0, posterior.fun = posterior0, npar = npar0, lower = lower0, upper = upper0)

write.csv(x = cbind(mcmc0$chain, mcmc0$all.lik), file ="results/mcmc_sample_VOC_ALL_v1_binomial.csv", row.names = F)

##### sensitivity analysis with only 90% of variants at points 2 and 3 #####
if(do_sensitivity_analysis <- T){
  
  freqdatatable0_SA <- freqdatatable0
  freqdatatable0_SA$NVOC[2:3] <- round(0.9 * freqdatatable0_SA$NVOC[2:3]) # 90% of NVOC only
  
  posterior0_SA <- function(par) prior0(par) - regional_lik(par, region = "all", epitable = epitable0, freqdatatable = freqdatatable0_SA)
  mcmc0_SA <- run_MCMC_Gibbs(start.value = opt0$pars, N.iter = niter, proposal.sd = proposal.sd0, posterior.fun = posterior0_SA, npar = npar0, lower = lower0, upper = upper0)
  write.csv(x = cbind(mcmc0_SA$chain, mcmc0_SA$all.lik), file ="results/mcmc_sample_VOC_ALL_SA_v1.csv", row.names = F)
  
}

########################################### OPTIMISE ALL REGIONS SIMULTANEOUSLY ###########################################

national_lik <- function(par, regionset = 1:nregions, return_all_regions = F){
  stopifnot(length(par) == 5*nregions+3)
  lls <- c()
  for(i in regionset){
    par_i <- par[(5*(i-1)+1):(5*i)] # parameters of region i
    par_general <- par[(5*nregions+1):(5*nregions+3)]
    myregion <- allregions[i]
    lls <- c(lls,
             regional_lik(
               par = c(
                 par_i, # region-specific pars
                 par_general # generic pars
               ),
               region = myregion,
               epitable = a2[which(a2$reg2 == myregion), ],
               freqdatatable = freqs4[which(freqs4$region == myregion), ]
             ))
  }
  if(return_all_regions){
    return(lls)
  } else {
    return(sum(lls))
  }
}
par <- c(rep(c(0.069, 1.1, 0.9, 0.1, 390), nregions), 0.5, 10, 0.1)
national_lik(par)

lower1 <- c(rep(c(0, 0.5, 0.5, 0.001, minday), nregions), 0, 0.001, 0.)
upper1 <- c(rep(c(0.1, 2, 2, 10     , maxday), nregions), 1, 1000, upper_limit_rho)
npar1  <- 5*nregions+3
init_fun1 <- function(x) sapply(1:npar1, function(i) runif(1, lower1[i], upper1[i]))
opt1 <- optim.fun.repeated(
  n.repeats = 5, lik.fun = national_lik, init.fun = init_fun1, verbose = T, lower = lower1, upper = upper1
)

# Bayesian approach:
prior <- function(par, regionset = 1:nregions, return_all_regions = F){
  
  if(return_all_regions){
    return(
      (dbeta(x = par[5*(regionset-1) + 1], shape1 = betaprior$alpha, shape2 = betaprior$beta, log = T)) + # beta prior for VOC frequency
      (dnorm(x = par[5*(regionset-1) + 2], mean = 1, sd = sd_normal, log = T)) + # normal prior for R wild type init
      (dnorm(x = par[5*(regionset-1) + 3], mean = 1, sd = sd_normal, log = T))   # normal prior for R wild type final
      # uniform prior for transmission advantage, k, inflection dates, size of the binomial (do not code them)
    )
  } else {
    return(
        sum(dbeta(x = par[5*(regionset-1) + 1], shape1 = betaprior$alpha, shape2 = betaprior$beta, log = T)) + # beta prior for VOC frequency
        sum(dnorm(x = par[5*(regionset-1) + 2], mean = 1, sd = sd_normal, log = T)) + # normal prior for R wild type init
        sum(dnorm(x = par[5*(regionset-1) + 3], mean = 1, sd = sd_normal, log = T))   # normal prior for R wild type final
      # uniform prior for transmission advantage, k, inflection dates, size of the binomial (do not code them)
    )
  }
}

posterior <- function(par, regionset = 1:nregions, return_all_regions = F) -national_lik(par = par, regionset = regionset, return_all_regions = return_all_regions) + prior(par = par, regionset = regionset, return_all_regions = return_all_regions)

# run the MCMC chain but do not update all parameters, only region by region
start.value <- opt1$pars
N.iter <- niter
proposal.sd1 <- c(rep(c(0.1, 0.01, 0.01, 1, 0.01), nregions), 0.05, 0.1, 0.1)
posterior.fun <- posterior
npar <- npar1
lower <- lower1
upper <- upper1

  # Gibbs sampling
  # log-normal proposal
  # code adapted from Henrik Salje, Cécile Tran Kiem  on https://zenodo.org/record/3813815#.YBx5MC1Q0Ut
  # explanations here for the log-normal proposal  https://umbertopicchini.wordpress.com/2017/12/18/tips-for-coding-a-metropolis-hastings-sampler/
  
  chain = array(dim = c(N.iter, npar))
  all.lik <- rep(NA, N.iter)
  all.lik.regions <- matrix(NA, nrow = N.iter, ncol = nregions)
  
  acceptance <- matrix(0, nrow = N.iter, ncol = npar)
  chain[1,] = start.value
  all.lik.regions[1, ] <- posterior.fun(par = start.value, return_all_regions = T)
  all.lik[1] <- posterior.fun(start.value, return_all_regions = F)
  stopifnot(abs(sum(all.lik.regions[1, ]) - all.lik[1]) < 1e-9)
  
  for (i in 2:N.iter){
    # print(i)
    if(i / 1000 == floor(i / 1000)){
      print(i)
      acceptance_byparam <- colMeans(acceptance[1:i, ]) # assess acceptance rate
      cat("acceptance rate in last 1000 generations: ", acceptance_byparam, "\n")
      cat("likelihood: ", all.lik[i-1], "\n")
      for(iparam in 1:npar) plot(chain[, iparam], type = "l", main = iparam)
    }
    
    # Getting parameters and posterior obtained at the previous iteration
    chain[i,]  <-  chain[i-1,]
    all.lik[i] <- all.lik[i-1]
    all.lik.regions[i, ] <- all.lik.regions[i-1, ]
    
    # Parameter updates
    for(iparam in 1:npar){
      
      old_param <- chain[i - 1, iparam]
      
      # sampling a new candidate in log-normal distribution with sd proposal.sd
      new_param <- old_param * exp(proposal.sd1[iparam] * rnorm(1))
      
      if(new_param > upper[iparam] || new_param < lower[iparam]){ # if not in the limits, just keep old parameter
        chain[i, iparam] <- old_param
      } else { # if in the limits
        
        chain[i, iparam] <- as.numeric(new_param)
        
        # posterior function for this 
        if(iparam <= 5*nregions){ # if parameter is region-specific, just update this region
          region_index <- ceiling(iparam/5)
          new_lik_singleregion <- posterior.fun(chain[i,], regionset = region_index,  return_all_regions = F)
          newlik <- all.lik[i] + new_lik_singleregion - all.lik.regions[i, region_index]
          #newlik_slow <- posterior.fun(chain[i,], return_all_regions = F)
          #stopifnot(abs(newlik-newlik_slow) < 1e-9)
        } else { # else update all regions (for the two last parameters)
          newlik_allregions <-  posterior.fun(chain[i,], return_all_regions = T)
          newlik <- sum(newlik_allregions)
        }
      
        # acceptance of rejection
        log_acceptance_ratio <- newlik - all.lik[i] + log(new_param) - log(old_param) # the latter factor is a correction for the log-normal proposal
        
        if(log(runif(1)) < log_acceptance_ratio){
          all.lik[i] <- newlik
          if(iparam <= 5*nregions) all.lik.regions[i, region_index] <- new_lik_singleregion else all.lik.regions[i,] <- newlik_allregions # update regional likelihoods
          stopifnot(abs(all.lik[i]-sum(all.lik.regions[i,])) < 1e-9)
          acceptance[i, iparam] <- 1
        } else{
          chain[i, iparam] <- old_param
        }
      }
    }
  }
mcmc <- list(chain = chain, all.lik = all.lik)
write.csv(x = cbind(mcmc$chain, mcmc$all.lik), file ="results/mcmc_sample_VOC_France_v4_binomial.csv", row.names = F)


########################################################################### ENVELOPPE FOR TRAJECTORIES ###########################################################################

regionset <- 1:nregions
burnin_fraction <- 0.2
samp <- sample((burnin_fraction*niter):niter, 1000)

# create enveloppe for trajectories
env_TOT0 <- c()
env_VOC0 <- c()
date_50_0 <- c()
env_VOC_f_0 <- c()
env_RWT0 <- c()
env_RVOC0 <- c()
for(i in samp){
  mypar <- mcmc0$chain[i, ]
  truc <- regional_lik(par = unlist(mypar), region = "ALL",
                       epitable = epitable0,
                       freqdatatable = freqdatatable0, return_full = T)$epitable
  Rinit <- mypar[2]
  Rfinal <- mypar[3]
  k <- mypar[4] # > 0, not too large
  inflection_date <- mypar[5] # must be in minday:maxday
  transmission_advantage <- mypar[6]
  
  RWT <- Rinit + (Rfinal - Rinit) / (1 + exp(- k * (minday:(maxday+90-1) - inflection_date)))
  RVOC <-  RWT * (1 + transmission_advantage)
  
  truc$freq <- truc$PVOCsim/(truc$PVOCsim + truc$PWTsim); date_50_0 <- c(date_50_0, truc$dateday[which(truc$freq>0.5)[1]])
  truc <- truc[!is.na(truc$PTOTsim), ]
  env_TOT0 <- rbind(env_TOT0, truc$PTOTsim)
  env_VOC0 <- rbind(env_VOC0, truc$PVOCsim)
  env_VOC_f_0 <- rbind(env_VOC_f_0, truc$PVOCsim/truc$PTOTsim)
  env_RWT0 <- rbind(env_RWT0, RWT)
  env_RVOC0 <- rbind(env_RVOC0, RVOC)
  
}
env_TOT0 <- apply(env_TOT0, 2, quantile, c(0.025, 0.5, 0.975))
env_VOC0 <- apply(env_VOC0, 2, quantile, c(0.025, 0.5, 0.975))
date_50_0 <- quantile(date_50_0, c(0.025, 0.5, 0.975), na.rm = T)
env_VOC_f_0 <- apply(env_VOC_f_0, 2, quantile, c(0.025, 0.5, 0.975))
env_RWT0 <- apply(env_RWT0, 2, quantile, c(0.025, 0.5, 0.975))
env_RVOC0 <- apply(env_RVOC0, 2, quantile, c(0.025, 0.5, 0.975))


# function to create enveloppe for trajectories
create_enveloppe <- function(rr){
  print("Warning: need to make absolutely sure that region_index is well defined")
  region_index <- which(rr==allregions)
    
  env_TOT <- c()
  env_VOC <- c()
  env_WT <- c()
  env_VOC_f <- c()
  date_50 <- c()
  env_RWT <- c()
  env_RVOC <- c()
  
  for(i in samp){
    
    truc <- regional_lik(par = unlist(mcmc$chain[i, c(5*(region_index-1) + (1:5), 5*nregions+1:3)]), region = rr,
                         epitable = a2[which(a2$reg2 == rr), ],
                         freqdatatable = freqs4[which(freqs4$region == rr), ], return_full = T)$epitable
    
    truc$freq <- truc$PVOCsim/(truc$PVOCsim + truc$PWTsim); date_50 <- c(date_50, truc$dateday[which(truc$freq>0.5)[1]])
    truc <- truc[!is.na(truc$PTOTsim), ]
    env_TOT <- rbind(env_TOT, truc$PTOTsim)
    env_VOC <- rbind(env_VOC, truc$PVOCsim)
    env_WT <- rbind(env_WT, truc$PTOTsim - truc$PVOCsim)
    env_VOC_f <- rbind(env_VOC_f, truc$PVOCsim/truc$PTOTsim)
    
    # sigmoids:
    mypar <- mcmc$chain[i, c(5*(region_index-1) + (1:5), 5*nregions+1:2)]
    Rinit <- mypar[2]
    Rfinal <- mypar[3]
    k <- mypar[4] # > 0, not too large
    inflection_date <- mypar[5] # must be in minday:maxday
    transmission_advantage <- mypar[6]
    
    RWT <- Rinit + (Rfinal - Rinit) / (1 + exp(- k * (minday:(maxday-1) - inflection_date)))
    RVOC <-  RWT * (1 + transmission_advantage)
    
    env_RWT <- rbind(env_RWT, RWT)
    env_RVOC <- rbind(env_RVOC, RVOC)
  }
  return(list(env_TOT = env_TOT,
              env_VOC = env_VOC,
              env_WT  = env_WT,
              env_VOC_f = env_VOC_f,
              date_50  = date_50,
              env_RWT = env_RWT,
              env_RVOC = env_RVOC)
         )
}

rm(region_index)
save.image(paste0("results/inference_", Sys.Date(), "_", upper_limit_rho, "_", meangt, ".RData"))

# PLOT FIGURES:
source("plot_all_figures.R")







