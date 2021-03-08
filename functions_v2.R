myintegrate <- function(f, tmin, tmax, deltat, param, ...){ # my function to simply integrate a function with parameters "param"
  if(tmin == tmax) return(0)
  tmp <- sapply(seq(tmin, tmax, deltat), f, param = param, ...)
  n <- length(tmp)
  weights <- rep(2, n); weights[1] <- 1; weights[n] <- 1
  return(sum(0.5 * (tmax-tmin)/(n-1) * weights * tmp))
}

solve_rtoR <- function(param, r){
  f <- function(t, param, r) dgamma(t, shape = param[1], scale = param[2]) * exp(-r * t)
  return(
    1 /  myintegrate(f, tmin = 0, tmax = 40, delta = 0.1, param = param, r = r)
  )
}

solve_Rtor <- function(param, R) uniroot(f = function(r) {solve_rtoR(param, r) - R}, interval = c(-1, 1))$root

# FD New version: with Gamma function, explicit formula
new_rtoR <- function(param, r){tmp <- (r + 1/param[2])^(-param[1]) * param[2]^(-param[1]); 1/tmp}
new_Rtor <- function(param, R){-1/param[2] + (param[2]^param[1]/R)^(-1/param[1])}
# endFD
  
# weekly_to_r <- function(mf){
#   log(mf)/7
# }

getbetaparams <- function(mu, var) {
  # gets parameter for a beta distrib with mean mu and variance var
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

optim.fun.repeated <- function(n.repeats, lik.fun, init.fun, verbose = F, lower = NULL, upper = NULL, ...){
  # a function to MINIMISE properly the negative log-likelihood function
  # using repeated use of function optim
  generate.init.fun <- match.fun(init.fun)
  lik.fun <- match.fun(lik.fun)
  all.opt <- list()
  for(ii in 1:n.repeats){
    init <- generate.init.fun()
    
    if(length(init)==1){ # one parameter function, use optimise
      if(is.null(upper) | is.null(lower)){ # if unspecified bounds for parameters, use (0,1)
        interval <- c(0, 1)
      } else {
        interval <- c(lower[1], upper[1])
      }
      all.opt[[ii]] <- optimise(lik.fun, interval = interval, ...)
      names(all.opt[[ii]]) <- c("par", "value")
      all.opt[[ii]]$message <- "Successful convergence"
    } else { # otherwise use nmk
      
      if(is.null(upper) | is.null(lower)){ # if unspecified bounds for parameters, use nmk
        all.opt[[ii]] <- nmk(init, lik.fun, ...)
        flag <- F
      } else {
        all.opt[[ii]] <- NA
        while(is.na(all.opt[[ii]])){
          init <- generate.init.fun()
          tryCatch(
            all.opt[[ii]] <- nmkb(init, lik.fun, lower = lower, upper = upper, ...),
            error = function(e){
              print(e)
            }
          )
        }
      }
    }
    
    if(verbose)print(ii)
  }
  if(verbose) cat(sum(unlist(lapply(all.opt, function(x)x$message == "Successful convergence"))), "optimisations converge\n")
  all.ll <- unlist(lapply(all.opt, function(x)x$value))
  min(all.ll) -> minll.pow
  if(verbose) cat(length(which(all.ll < minll.pow + 1 & all.ll > minll.pow)), " optim within 1 log lik of minimum likelihood\n")
  output <- list(minll.pow, all.opt[[which.min(all.ll)]]$par)
  names(output) <- c("lik", "pars")
  return(output)
}
run_MCMC_Gibbs <- function(start.value, N.iter, proposal.sd, posterior.fun, npar, lower, upper, verbose = T, ...){
  # Gibbs sampling
  # log-normal proposal
  # code adapted from Henrik Salje, Cécile Tran Kiem  on https://zenodo.org/record/3813815#.YBx5MC1Q0Ut
  # explanations here for the log-normal proposal  https://umbertopicchini.wordpress.com/2017/12/18/tips-for-coding-a-metropolis-hastings-sampler/

  chain = array(dim = c(N.iter, npar))
  all.lik <- rep(NA, N.iter)
  acceptance <- matrix(0, nrow = N.iter, ncol = npar)
  chain[1,] = start.value
  all.lik[1] <- posterior.fun(start.value, ...)

  for (i in 2:N.iter){
    
    if(i / 1000 == floor(i / 1000) & verbose){
      print(i)
      acceptance_byparam <- colMeans(acceptance[1:i, ]) # assess acceptance rate
      cat("acceptance rate in last 1000 generations: ", acceptance_byparam, "\n")
      cat("likelihood: ", all.lik[i-1], "\n")
      for(iparam in 1:npar) plot(chain[, iparam], type = "l", main = iparam)
    }
    
    # Getting parameters and posterior obtained at the previous iteration
    chain[i,]  <-  chain[i-1,]
    all.lik[i] <- all.lik[i-1]
    
    # Parameter updates
    for(iparam in 1:npar){
      
      old_param <- chain[i - 1, iparam]
      
      # sampling a new candidate in log-normal distribution with sd proposal.sd
      new_param <- old_param * exp(proposal.sd[iparam] * rnorm(1))
      
      if(new_param > upper[iparam] || new_param < lower[iparam]){ # if not in the limits, just keep old parameter
        chain[i, iparam] <- old_param
      } else { # if in the limits
        
        chain[i, iparam] <- as.numeric(new_param)
        
        # posterior function for this 
        newlik <- posterior.fun(chain[i,], ...)
        
        # acceptance of rejection
        log_acceptance_ratio <- newlik - all.lik[i] + log(new_param) - log(old_param) # the latter factor is a correction for the log-normal proposal
        
        if(log(runif(1)) < log_acceptance_ratio){
          all.lik[i] <- newlik
          acceptance[i, iparam] <- 1
        } else{
          chain[i, iparam] <- old_param
        }
      }
    }
  }
  return(list(chain = chain, all.lik = all.lik))
}
# testing the new MCMC function:
# fakedata <- rnorm(1000, mean = 1, sd = 3)
# mypost <- function(par) sum(dnorm(x = fakedata, mean = par[1], sd = par[2], log = T))
# mcmc <- run_MCMC_Gibbs(start.value = c(2, 1), N.iter = 10000, proposal.sd = c(0.1, 0.1), posterior.fun = mypost, npar = 2, lower=  c(0,0), upper = c(10, 10), verbose = T)
# plot(mcmc$all.lik)
# plot(mcmc$chain[,1])
# plot(mcmc$chain[,2])
run_MCMC_MH <- function(start.value, N.iter, proposal.cov, posterior.fun, npar, lower, upper, verbose = T, ...){ # same as run_MCMC but tuning the proposal distribution
  # THIS IS THE LATEST AND MOST CORRECT VERSION OF THE MCMC ALGO
  # Metropolis-Hastings algorithm (requires symmetric proposal distribution -> multivariate Gaussian OK)
  init.proposal.cov <- proposal.cov
  factor.cov <- 10
  
  posterior.fun <- match.fun(posterior.fun)
  stopifnot(dim(proposal.cov) == c(npar, npar))
  
  chain = array(dim = c(N.iter + 1, npar))
  all.lik <- rep(NA, N.iter + 1)
  chain[1,] = start.value
  all.lik[1] <- posterior.fun(start.value, ...)
  
  for (i in 1:N.iter){
    if(i / 1000 == floor(i / 1000)){
      if(verbose) print(i)
      acceptance <- 1 - mean(duplicated(chain[(i - 1000 + 1):i, ])) # assess acceptance rate
      if(verbose) cat("acceptance rate in last 1000 generations: ", acceptance, "\n")
      if(verbose) cat("likelihood: ", all.lik[i], "\n")
      df.mcmc <- data.frame(chain[1:i, ])
      df.mcmc <- df.mcmc[!duplicated(df.mcmc), ]    # unique values
      if(nrow(df.mcmc) == 1) {
        print("no new parameter accepted; reducing the proposal distribution")
        proposal.cov <- init.proposal.cov / factor.cov
      } else {
        if(npar > 1){
          proposal.cov <- cov(df.mcmc)                  # update proposal cov
        } else {
          proposal.cov <- var(df.mcmc)
        }
      }
      if(verbose) print("updated proposal covariance:")
      if(verbose) print(proposal.cov / factor.cov)
    }
    # draw new parameters
    if(npar > 1){
      proposal <- chain[i,] + mvrnorm(1, mu = rep(0, npar), Sigma = proposal.cov / factor.cov) # multivariate normal distribution
    } else {
      proposal <- chain[i,] + rnorm(1, mean = 0, sd = sqrt(proposal.cov)) # univariate normal
    }
    # NOTE HERE WE CANNOT SIMPLY DRAW NEW PROPOSAL UNTIL WE FIND SUITABLE ONE - THIS WOULD BIAS THE POSTERIOR SAMPLE https://darrenjw.wordpress.com/2012/06/04/metropolis-hastings-mcmc-when-the-proposal-and-target-have-differing-support/
    # see also https://stats.stackexchange.com/questions/73885/mcmc-on-a-bounded-parameter-space
    
    if(all(proposal > lower & proposal < upper)) {
      all.lik[i + 1] <- posterior.fun(proposal, ...)
    } else {
      all.lik[i + 1] <- NA # set to NA if unfeasible value
    }
    
    if(!is.na(all.lik[i+1])) probab <- exp(all.lik[i + 1] - all.lik[i]) else probab <- 0 # and set acceptance proba to 0 if NA...
    
    if(runif(1) < probab){
      chain[i+1,] <- proposal
    }else{
      chain[i+1,] <- chain[i,]
      all.lik[i + 1] <- all.lik[i]
    }
  }
  return(list(chain = chain, all.lik = all.lik))
}
nmkb <- function (par, fn, lower = -Inf, upper = Inf, control = list(), mytol = 1e-2, ...)  # from dfoptim package
{
  ctrl <- list(tol = mytol, maxfeval = min(5000, max(1500, 
                                                     20 * length(par)^2)), regsimp = TRUE, maximize = FALSE, 
               restarts.max = 3, trace = FALSE)
  namc <- match.arg(names(control), choices = names(ctrl), 
                    several.ok = TRUE)
  if (!all(namc %in% names(ctrl))) 
    stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
  if (!is.null(names(control))) 
    ctrl[namc] <- control
  ftol <- ctrl$tol
  maxfeval <- ctrl$maxfeval
  regsimp <- ctrl$regsimp
  restarts.max <- ctrl$restarts.max
  maximize <- ctrl$maximize
  trace <- ctrl$trace
  n <- length(par)
  g <- function(x) {
    gx <- x
    gx[c1] <- atanh(2 * (x[c1] - lower[c1])/(upper[c1] - 
                                               lower[c1]) - 1)
    gx[c3] <- log(x[c3] - lower[c3])
    gx[c4] <- log(upper[c4] - x[c4])
    gx
  }
  ginv <- function(x) {
    gix <- x
    gix[c1] <- lower[c1] + (upper[c1] - lower[c1])/2 * (1 + 
                                                          tanh(x[c1]))
    gix[c3] <- lower[c3] + exp(x[c3])
    gix[c4] <- upper[c4] - exp(x[c4])
    gix
  }
  if (length(lower) == 1) 
    lower <- rep(lower, n)
  if (length(upper) == 1) 
    upper <- rep(upper, n)
  if (any(c(par <= lower, upper <= par))) 
    stop("Infeasible starting values!", call. = FALSE)
  low.finite <- is.finite(lower)
  upp.finite <- is.finite(upper)
  c1 <- low.finite & upp.finite
  c2 <- !(low.finite | upp.finite)
  c3 <- !(c1 | c2) & low.finite
  c4 <- !(c1 | c2) & upp.finite
  if (all(c2)) 
    stop("Use `nmk()' for unconstrained optimization!", call. = FALSE)
  if (maximize) 
    fnmb <- function(par) -fn(ginv(par), ...)
  else fnmb <- function(par) fn(ginv(par), ...)
  x0 <- g(par)
  if (n == 1) 
    stop(call. = FALSE, "Use `optimize' for univariate optimization")
  if (n > 30) 
    warning("Nelder-Mead should not be used for high-dimensional optimization")
  V <- cbind(rep(0, n), diag(n))
  f <- rep(0, n + 1)
  f[1] <- fnmb(x0)
  V[, 1] <- x0
  scale <- max(1, sqrt(sum(x0^2)))
  if (regsimp) {
    alpha <- scale/(n * sqrt(2)) * c(sqrt(n + 1) + n - 1, 
                                     sqrt(n + 1) - 1)
    V[, -1] <- (x0 + alpha[2])
    diag(V[, -1]) <- x0[1:n] + alpha[1]
    for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
  }
  else {
    V[, -1] <- x0 + scale * V[, -1]
    for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
  }
  f[is.nan(f)] <- Inf
  nf <- n + 1
  ord <- order(f)
  f <- f[ord]
  V <- V[, ord]
  rho <- 1
  gamma <- 0.5
  chi <- 2
  sigma <- 0.5
  conv <- 1
  oshrink <- 1
  restarts <- 0
  orth <- 0
  dist <- f[n + 1] - f[1]
  v <- V[, -1] - V[, 1]
  delf <- f[-1] - f[1]
  diam <- sqrt(colSums(v^2))
  sgrad <- c(solve(t(v), delf))
  alpha <- 1e-04 * max(diam)/sqrt(sum(sgrad^2))
  simplex.size <- sum(abs(V[, -1] - V[, 1]))/max(1, sum(abs(V[, 
                                                              1])))
  itc <- 0
  conv <- 0
  message <- "Succesful convergence"
  while (nf < maxfeval & restarts < restarts.max & dist > ftol & 
         simplex.size > 1e-06) {
    fbc <- mean(f)
    happy <- 0
    itc <- itc + 1
    xbar <- rowMeans(V[, 1:n])
    xr <- (1 + rho) * xbar - rho * V[, n + 1]
    fr <- fnmb(xr)
    nf <- nf + 1
    if (is.nan(fr)) 
      fr <- Inf
    if (fr >= f[1] & fr < f[n]) {
      happy <- 1
      xnew <- xr
      fnew <- fr
    }
    else if (fr < f[1]) {
      xe <- (1 + rho * chi) * xbar - rho * chi * V[, n + 
                                                     1]
      fe <- fnmb(xe)
      if (is.nan(fe)) 
        fe <- Inf
      nf <- nf + 1
      if (fe < fr) {
        xnew <- xe
        fnew <- fe
        happy <- 1
      }
      else {
        xnew <- xr
        fnew <- fr
        happy <- 1
      }
    }
    else if (fr >= f[n] & fr < f[n + 1]) {
      xc <- (1 + rho * gamma) * xbar - rho * gamma * V[, 
                                                       n + 1]
      fc <- fnmb(xc)
      if (is.nan(fc)) 
        fc <- Inf
      nf <- nf + 1
      if (fc <= fr) {
        xnew <- xc
        fnew <- fc
        happy <- 1
      }
    }
    else if (fr >= f[n + 1]) {
      xc <- (1 - gamma) * xbar + gamma * V[, n + 1]
      fc <- fnmb(xc)
      if (is.nan(fc)) 
        fc <- Inf
      nf <- nf + 1
      if (fc < f[n + 1]) {
        xnew <- xc
        fnew <- fc
        happy <- 1
      }
    }
    if (happy == 1 & oshrink == 1) {
      fbt <- mean(c(f[1:n], fnew))
      delfb <- fbt - fbc
      armtst <- alpha * sum(sgrad^2)
      if (delfb > -armtst/n) {
        if (trace) 
          cat("Trouble - restarting: \n")
        restarts <- restarts + 1
        orth <- 1
        diams <- min(diam)
        sx <- sign(0.5 * sign(sgrad))
        happy <- 0
        V[, -1] <- V[, 1]
        diag(V[, -1]) <- diag(V[, -1]) - diams * sx[1:n]
      }
    }
    if (happy == 1) {
      V[, n + 1] <- xnew
      f[n + 1] <- fnew
      ord <- order(f)
      V <- V[, ord]
      f <- f[ord]
    }
    else if (happy == 0 & restarts < restarts.max) {
      if (orth == 0) 
        orth <- 1
      V[, -1] <- V[, 1] - sigma * (V[, -1] - V[, 1])
      for (j in 2:ncol(V)) f[j] <- fnmb(V[, j])
      nf <- nf + n
      ord <- order(f)
      V <- V[, ord]
      f <- f[ord]
    }
    v <- V[, -1] - V[, 1]
    delf <- f[-1] - f[1]
    diam <- sqrt(colSums(v^2))
    simplex.size <- sum(abs(v))/max(1, sum(abs(V[, 1])))
    f[is.nan(f)] <- Inf
    dist <- f[n + 1] - f[1]
    sgrad <- c(solve(t(v), delf))
    if (trace & !(itc%%2)) 
      cat("iter: ", itc, "\n", "value: ", f[1], "\n")
  }
  if (dist <= ftol | simplex.size <= 1e-06) {
    conv <- 0
    message <- "Successful convergence"
  }
  else if (nf >= maxfeval) {
    conv <- 1
    message <- "Maximum number of fevals exceeded"
  }
  else if (restarts >= restarts.max) {
    conv <- 2
    message <- "Stagnation in Nelder-Mead"
  }
  return(list(par = ginv(V[, 1]), value = f[1] * (-1)^maximize, 
              feval = nf, restarts = restarts, convergence = conv, 
              message = message))
}
myattach <- function(mylist){
  stopifnot(is.list(mylist))
  mynames <- names(mylist)
  n <- length(mylist)
  stopifnot(length(mynames) == n)
  cat("attaching: ", mynames, "\n")
  for(i in 1:n) assign(x = mynames[i], value = mylist[[i]], envir = .GlobalEnv)
}

###############                                   LOAD AND CLEAN DATA                                     ############### 

# frequency of the variant in cases 7-8 Jan according to SPF report 28/01/2021 UPDATED
# total frequency 3.3%
freqs <- read.csv("data/freqs_data_flash.csv", sep = ";")
freqs$f_neg <- freqs$N_SNEG / freqs$N
freqs$fVOC_in_neg <- freqs$NVOC / freqs$NSEQ_OK
freqs$fVOC <- freqs$fVOC_in_neg * freqs$f_neg
plot((1/freqs$fVOC), (1/freqs$fVOC_in_neg), pch = 20)
# if frequency of VOC among SGTF y = x / (a + x) where a is background frequency and x is frequency of VOC,
# then (1/y) = (a/x) + 1
x <- (1/freqs$fVOC)
y <- (1/freqs$fVOC_in_neg)
lm0 <- lm(y ~ x) # not working much

freqs2 <- read.csv("data/freqs_data_flash2_v4.csv", sep = ";")
freqs2$fVOC <- freqs2$NVOC/freqs2$N
freqs2 <- freqs2[-3]
freqs2 <- freqs2[c("region", "N", "fVOC")]

# https://www.data.gouv.fr/fr/datasets/donnees-relatives-aux-resultats-des-tests-virologiques-covid-19/
# https://www.data.gouv.fr/fr/datasets/taux-dincidence-de-lepidemie-de-covid-19/#_

a <- read.csv("data/sp-pe-tb-quot-reg-2021-02-23-19h20.csv", sep = ";") # daily cases per region
regioncode <- read.csv("data/region2020.csv")
a$reg2 <- regioncode$ncc[match(a$reg, regioncode$reg)] # add regions
table(a$reg2, useNA = "ifany") # 3 regions are NA -> 5, 7, 8
unique(a$reg[is.na(a$reg2)]) # the three regions 5, 7, 8 are not present in regioncode
# according to https://www.data.gouv.fr/en/datasets/donnees-relatives-aux-resultats-des-tests-virologiques-covid-19/ this is: 975 977 978
a <- a[a$reg2 %in% freqs$region, ]

# convert dates in day of year:
a$year <- as.numeric(strftime(a$jour, format = "%y"))
a$dateday <- as.numeric(strftime(a$jour, format = "%j"))
a$dateday[a$year == 21] <- a$dateday[a$year == 21] + 366

a2 <- a[a$dateday >366  & a$cl_age90 == 0,]
maxday <- max(a2$dateday)
# 7-days smoothed cases?
a2$logP <- log(a2$P)

allregions <- unique(a2$reg2); allregions <- as.character(allregions[!is.na(allregions)])
nregions <- length(allregions)
ndates <- nrow(a2) / nregions
stopifnot(all(a2$reg2 == rep(allregions, ndates)))
a2$cases <- NA
for(rr in allregions){
  suba2 <- a2[which(a2$reg2==rr),]
  for(dd in (370):(maxday-3)){ # 370 is min day where smoothed cases can be computed

    a2[which(a2$reg2 == rr & a2$dateday == dd), "Psmoothed"] <- round(exp(sum(suba2$logP[suba2$dateday %in% (dd-3):(dd+3)])/7)) # number of cases smootheda2_withpred[a2_withpred$reg2==rr, ]
  }
}

freqs <- freqs[match(allregions, freqs$region), ]; stopifnot(all(freqs$region==allregions))

initdate <- "2021-01-07" # date of initial "flash" investigation
minday <- unique(a$dateday[a$jour==initdate])
a2 <- a2[a2$dateday>=minday,]
n <- maxday-minday+1
stopifnot(nrow(a2) == n * nregions) # check we do have n * nregions in epi datatable

#View(a2)

getlik_frequency <- function(NVOCdata, Ndata, daymin, daymax, simtable, myrho){ # get likelihood component
  
  # get simulated frequency of VOC in sim:
  idx <- which(simtable$dateday >= daymin  & simtable$dateday <= daymax)
  if(length(idx)<1) stop("data not present")
  
  fsim <- sum(simtable$PVOCsim[idx]) / sum(simtable$PTOTsim[idx]) 
  
  # probability to get NVOCdata VOC if true frequency is fsim
  return(
    - sum(dbetabinom(x = NVOCdata,  prob = fsim, rho = myrho, size = Ndata, log = T))
    # note to myself: the maximum likelihood estimator of beta binomial will not be fsim = NVOCdata/Ndata but higher
    #- sum(dbinom(x = NVOCdata,  prob = fsim, size = Ndata, log = T))
  )
}


############################################################ OLD STUFF ############################################################

if(F){ # DO NOT EVALUATE
  
process_params <- function(par, nr){
  if(length(par)!=nregions+nr+2) stop()
  f_init <- par[1:nregions] # initial frequency in each region
  RWT <- par[(nregions+1):(nregions+nr)] # region-specific R
  transmission_advantage <- par[nregions+nr+1]
  size <- par[nregions+nr+2]
  return(list(f_init = f_init, RWT = RWT, transmission_advantage = transmission_advantage, size = size))
}


mylik <- function(par, return_full = F, nr, regions = allregions, freqdatatable = freqs4){ # original function with the option to have 1 parameter for all regions or many
  
  stopifnot(nr==1 | nr == nregions)
  myattach(process_params(par, nr)) # process parameters
  
  RVOC <- RWT * transmission_advantage
  rVOC <- new_fast_solve_Rtor(RVOC) # FD
  rWT <- new_fast_solve_Rtor(RWT) # FD
  
  if(nr==nregions){
    names(rVOC) <- regions
    names(rWT) <- regions
  }
  names(f_init) <- regions
  
  a2$PVOCsim <- NA
  a2$PWTsim <- NA
  
  if(return_full){
    a2 <- rbind.fill(a2 , expand.grid(reg2 = allregions, dateday = (maxday+1):(maxday+90)))
  }
  
  # NOW EXTRAPOLATE Psmoothed VARIANT / WT AT EACH DAY IN EACH REGION GIVEN RWT AND TRANSMISSION ADVANTAGE
  for(rr in regions){
    
    if(nr==nregions){
      myrVOC <- rVOC[rr]
      myrWT <- rWT[rr]
    } else {
      myrVOC <- rVOC
      myrWT <- rWT
    }
    
    # initialise prevalence of VOC and WT using frequency parameters
    initPVOC <- f_init[rr] * a2[which(a2$jour==initdate & a2$reg2==rr), "Psmoothed"] 
    initPWT <- a2[which(a2$jour==initdate & a2$reg2==rr), "Psmoothed"] - initPVOC
    
    a2[which(a2$jour==initdate & a2$reg2==rr), "PVOCsim"] <- initPVOC
    a2[which(a2$jour==initdate & a2$reg2==rr), "PWTsim"] <- initPWT
    
    # extrapolate with exp growth from minday+1 to maxday
    if(!return_full){
      a2[which(a2$reg2==rr)[match((minday+1):maxday, a2[which(a2$reg2==rr), "dateday"])], "PVOCsim"] <- initPVOC * exp(myrVOC * (1:(maxday-minday)))
      a2[which(a2$reg2==rr)[match((minday+1):maxday, a2[which(a2$reg2==rr), "dateday"])], "PWTsim"] <- initPWT * exp(myrWT * (1:(maxday-minday)))
    } else {
      a2[which(a2$reg2==rr)[match((minday+1):(maxday+90), a2[which(a2$reg2==rr), "dateday"])], "PVOCsim"] <- initPVOC * exp(myrVOC * (1:(maxday+90-minday)))
      a2[which(a2$reg2==rr)[match((minday+1):(maxday+90), a2[which(a2$reg2==rr), "dateday"])], "PWTsim"] <- initPWT * exp(myrWT * (1:(maxday+90-minday)))
    }
  }
  a2$PTOTsim <- a2$PVOCsim+a2$PWTsim
  #a2$jour_region <- paste(a2$jour, a2$reg2, sep = "_")
  #a2$fsim <- a2$PVOCsim / a2$PTOTsim
  
  # component of the likelihood based on number of Psmoothed daily compared to inital prediction
  ll0 <- - sum(dnbinom(x = a2$Psmoothed, size = size, mu = a2$PTOTsim, log = T), na.rm = T)
  
  # LIKELIHOOD COMPONENT CORRESPONDING TO FREQUENCIES
  
  # FLASH 7-8 JAN BINOMIAL LIKELIHOOD FREQUENCY
  lik_freqs <- sum(
    apply(freqdatatable, 1, function(vec) getlik_frequency(NVOCdata = as.numeric(vec["NVOC"]), Ndata = as.numeric(vec["N"]), region = vec["region"], daymin = as.numeric(vec["daymin"]), daymax = as.numeric(vec["daymax"]), simtable = a2))
  )
  
  if(return_full){
    return(list(ll = ll0+sum(lik_freqs), a2 = a2))
  } else {
    return(ll0+sum(lik_freqs))
  }
}

# OBSOLETE!!!
getlik_frequency <- function(NVOCdata, Ndata, region, daymin, daymax, simtable, myrho){ # get likelihood component
  # get simulated frequency of VOC in sim:
  print("this version is obsolete!!!")
  idx <- which(simtable$dateday >= daymin  & simtable$dateday <= daymax & simtable$reg2==region)
  if(length(idx)<1) stop("data not present")
  
  fsim <- sum(simtable$PVOCsim[idx]) / sum(simtable$PTOTsim[idx]) 
  
  # probability to get NVOCdata VOC if true frequency is fsim
  return(
    - sum(dbetabinom(x = NVOCdata,  prob = fsim, rho = myrho, size = Ndata, log = T))
  )
}

# playing around with betabinomial
library(VGAM)
NVOC <- 148
N <- 2149
plot(probs <- seq(0, 0.2, 0.001),
     y3 <- dbetabinom(x = NVOC,  prob = probs, rho = 0.1, size = N, log = T), type = "l", ylim = c(-10,0)
)
points(probs,
       y2 <- dbetabinom(x = NVOC,  prob = probs, rho = 0.05, size = N, log = T), type = "l", col = "blue"
)
points(probs,
       y1 <- dbetabinom(x = NVOC,  prob = probs, rho = 0.01, size = N, log = T), type = "l", col = "blue"
)
points(probs,
       y0 <- dbetabinom(x = NVOC,  prob = probs, rho = 0.0, size = N, log = T), type = "l", col = "red"
)

probs[which.max(y3)]
probs[which.max(y2)]
probs[which.max(y1)]
probs[which.max(y0)]


}


