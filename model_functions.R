###########################################################################
###########################################################################
###########################################################################
#### All functions required to run sociohydrologic model  ####
###########################################################################
###########################################################################
###########################################################################

## TO DO:

# search: ## !! ## to look for to-do's within script, such as edits to be made.

## NOTES:

#  "Qsw.use.S" in out.S is surface water used (reported by the S = social model) in a given time step (t). "Qsw.store.H" in out.H is surface water that is available (reported by the H = hydrology model) from the previous time step (t-1); the water remaining at the end of a given time step (t-1) from the hydrology model is the water available at the start of the subsequent time step (t).
# l indexes current time step; in a given iteration, l-1 is (t-1); l = 1 is the initialization step.

## load packages for analysis

require(reshape2)
require(ggplot2)
require(hydromad) # loads zoo

###########################################################################
###########################################################################
###########################################################################
#### Initialization Functions  ####
###########################################################################
###########################################################################
###########################################################################

f.init <- function( # all defaults set
  
  #### Specify number of iterations
  iters = 10,
  
  #### Specify groups (names)
  J = c("a", "u"), # j indexes groups (ag = a, urban = u) in social model (only)
  
  # Note: In the methodology description, i indexes water type (surface water = sw, groundwater = gw) or product (crop = c, hours = h), but this index is only used for methods explanation, not in the functions to this model.
  # all water quantities are assumed basin- or multi-basin averaged depths [mm]
  
  #### Specify state variables
  vars.out.H = c("P", "E", "Quse.H", "Qsw.store.H", "Qgw.store.H", "Qstream.H", "Qwaste.H"), # rain, ET, sw store, gw store, streamflow (environmental flows), wastewater flows
  vars.out.S = c("X", "N", "Pi", "Qsw.access", "Qgw.access", "Qsw.use.S", "Qgw.use.S", "Q.i", "p.i"), # along J; note: Q.i and p.i are total quantity of products (h,c) and per-unit price of products (h,c) -> h (hours) corresponds to "u" group, c (crops) corresponds to "a" group.
  
  #### Specify starting values
  
  ## Hydrology: singe value applies to all groups
  start.Qsw.store.H = 10^2, # surface water available
  start.Qgw.store.H = 10^3, # groundwater available
  
  ## Social: in J (a,u) order
  start.Qh = 50, # total hours worked by "u" group
  start.Qc = 100, # total crops produced by "a" group
  
  #### Specify static parameters
  
  ## Social:
  # for population and group switching functions
  params.pop.S = list(N = c(100,100), YN.delta.X = "N" , delta.X.limit = NA),
  # for water use avail. and access functions
  params.water.S = list(Qsw.limit.a = 2, Qgw.limit.a = 5, Qsw.limit.u = 1, Qgw.limit.u = 1, Qsw.fraction.u = 0.8),
  # for price (crops, hours) function and set water prices
  params.price.S = list(k.a = 1, l.a = 1, m.a = 1, k.u = 1, l.u = 1, m.u = 1),
  # for production (crops, hours) function
  params.prod.S = list(alpha.a=1, beta.sw.a = 0.3, beta.gw.a = 0.3, hm.u = 5, alpha.u = 1, beta.w.u = 0.1),
  params.profit.S = list(psw.a = 0.1, pgw.a = 0.2, pw.u = 0.5),
  
  ## Hydrology:
  
  params.H = list(Qd.frac.E.a = 0.85,  Qd.frac.E.u = 0.85,  Qd.frac.W.u = 0.5, Q.A.frac = 0.05, Q.S.frac = 0.7, s1.use.frac.a = 0.2, s1.use.frac.u = 0.4),
  
  ## Modifications:
  
  params.M = list(Qsw.slope.trigger_psw.a_percent = NULL, Qgw.slope.trigger_pgw.a_percent = NULL),
  
  ... 
  
){
  
  #### Iterations/time steps
  counter <- 1:iters
  
  #### Output Arrays
  out <- list(
    H = array(NA, dim=c(length(counter),length(vars.out.H)), dimnames = list(counter, vars.out.H)), # 2-dimensiona array; indexed along iteration (rows), variables (cols)
    S = array(NA, dim=c(length(counter),length(vars.out.S),length(J)), dimnames = list(counter, vars.out.S, J)), # 3-dimensional array; indexed along iteration (rows), variables (cols), and groups (layers)
    log = array(NA, dim=c(length(counter),1), dimnames = list(counter, "Notes"))
  )
  
  out$log[1,"Notes"] <- "Initial values."
  
  #### Enter starting values into output arrays
  
  ## Hydrology
  
  out$H[1,"Qsw.store.H"] <- start.Qsw.store.H;
  out$H[1,"Qgw.store.H"] <- start.Qgw.store.H;
  
  ## Social
  
  #out$S[1,"Pi",] <- Pi;
  out$S[1,"N",] <- params.pop.S$N;
  out$S[1,"Q.i","u"] <- start.Qh;
  out$S[1,"Q.i","a"] <- start.Qc;
  
  # proportion of total population in group
  out$S[1,"X",] <- c(a = out$S[1,"N","a"]/sum(out$S[1,"N",]), u = out$S[1,"N","u"]/sum(out$S[1,"N",])); 
  if(sum(out$S[1,"X",]) != 1) stop("Sum of proportions of total population in each group does not sum to 1.") 
  
  ## Modifications
  ## !! ##
  
  ## Elements to Return
  
  returns <- list(
    out = out, 
    J = J,
    vars.out.H = vars.out.H, 
    vars.out.S = vars.out.S, 
    params.H = params.H, 
    params.pop.S = params.pop.S,
    params.water.S = params.water.S,
    params.price.S = params.price.S,
    params.prod.S = params.prod.S,
    params.profit.S = params.profit.S,
    params.M = params.M
  )
  
  return(returns)
  
}


###########################################################################
###########################################################################
###########################################################################
#### Functions Internal to Submodels ####
###########################################################################
###########################################################################
###########################################################################

#############################################################################
########### Social ###########
#############################################################################

########### Urban ###########

#### Optimization ####

## objective/profit function, with embedded production function  
# X == quantities of water to be used by individuals (here, just one X here: qw.used.u)
f.pi.u <-function(X,x,l){(x$out$S[l,"p.i","u"] * (x$params.prod.S$hm.u + (x$params.prod.S$alpha.u * X^x$params.prod.S$beta.w.u)) ) - (x$params.profit.S$pw.u * X)}

## gradient of objective function (profit): derivative wrt X; returns a value
f.pi.u.grad <-function(X,x,l){( x$out$S[l,"p.i","u"] * x$params.prod.S$alpha.u * x$params.prod.S$beta.w.u * (X^(x$params.prod.S$beta.w.u - 1)) ) - x$params.profit.S$pw.u}

########### Agricultural ###########

#### Optimization ####

## objective/profit function, with embedded production function 

# when sw and gw are available. X == quantities of water to be used by individuals (here, two X: qsw.used.a [1], qgw.used.a [2])
f.pi.a <-function(X,x,l){ (x$out$S[l,"p.i","a"] * ( x$params.prod.S$alpha.a * (X[1]^x$params.prod.S$beta.sw.a) * (X[2]^x$params.prod.S$beta.gw.a) )) - (x$params.profit.S$psw.a * X[1]) - (x$params.profit.S$pgw.a * X[2]) }

# when only sw is available. X = sw
f.pi.a_sw <-function(X,x,l){ (x$out$S[l,"p.i","a"] * ( x$params.prod.S$alpha.a * (X^x$params.prod.S$beta.sw.a) )) - (x$params.profit.S$psw.a * X) }

# when only gw is available. X = gw
f.pi.a_gw <-function(X,x,l){ (x$out$S[l,"p.i","a"] * ( x$params.prod.S$alpha.a * (X^x$params.prod.S$beta.gw.a) )) - (x$params.profit.S$pgw.a * X) }

## sw and gw: gradient of objective function (profit): derivative wrt X[1] and X[2]; returns a vector
f.pi.a.grad <-function(X,x,l){
  # dpi/dsw
  dsw <- (x$out$S[l,"p.i","a"] * x$params.prod.S$alpha.a * x$params.prod.S$beta.sw.a * (X[1]^(x$params.prod.S$beta.sw.a - 1)) * (X[2]^x$params.prod.S$beta.gw.a)) - x$params.profit.S$psw.a
  # dpi/dgw
  dgw <- (x$out$S[l,"p.i","a"] * x$params.prod.S$alpha.a * x$params.prod.S$beta.gw.a * (X[1]^x$params.prod.S$beta.sw.a) * (X[2]^(x$params.prod.S$beta.gw.a-1))) - x$params.profit.S$pgw.a
  # return both
  return(c(dsw,dgw))
}

## sw only: derivative wrt X[1]
f.pi.a.grad_sw <-function(X,x,l){
  (x$out$S[l,"p.i","a"] * x$params.prod.S$alpha.a * x$params.prod.S$beta.sw.a * (X^(x$params.prod.S$beta.sw.a - 1))) - x$params.profit.S$psw.a
}

## gw only: derivative wrt X[2]
f.pi.a.grad_gw <-function(X,x,l){
  (x$out$S[l,"p.i","a"] * x$params.prod.S$alpha.a * x$params.prod.S$beta.gw.a * (X^(x$params.prod.S$beta.gw.a-1))) - x$params.profit.S$pgw.a
}

#############################################################################
########### Hydrology ###########
#############################################################################

##### Data Simulation

## simulates P,ET, returns zoo ts object with index (not date format) times

f.PE <- function(cor = "neg", pmean.s1 = 8, pmean.s2 = 3, yrs = 10, s1.days = 200, s2.days = 165){
  
  # cor is correlation between rainfall and ET
  # pmean.s1 and pmean.s2 are daily average rainfall intensities by season
  # yrs is total number of yrs to simulate, should be same as iterations (and first year will be used to initialize the hydrologic model)
  # s1.days is total number of days in first (e.g. wet) season, which is split (beginning and end of year starting in January); s2.days is total number in 2nd (dry) season (middle of year)
  # assuming everything in mm
  
  if (!(s1.days + s2.days == 365)) stop ("Number of days in seasons must sum to 365.")
  
  cor.set <- cor
  
  # sample from rain mean distribution
  pmean1 <- pmean.s1
  pmean2 <- pmean.s2
  
  #yrs <- (yrs + 1) #add extra year for trimming
  
  # time series of rainfall; starts in January.
  for (i in 1:yrs){
    if (i == 1){
      rain <- c(rexp( s1.days/2, rate = 1/pmean1), rexp(s2.days, rate = 1/pmean2), rexp( s1.days/2,rate = 1/pmean1))
    } else {
      rain <- c(rain, rexp( s1.days/2, rate = 1/pmean1), rexp(s2.days, rate = 1/pmean2), rexp( s1.days/2,rate = 1/pmean1))
    }
  }
  #plot(rain, type="l")
  
  ## et
  x<-1:(365*yrs)
  A <- (pmean1-3)/pmean2 # amplitude
  B <- (2*pi)/365 # wavelength, for one year = (2*pi)/length(x) 
  D <- 2 + A # vertical offset
  # horizontal offset
  if (cor.set == "neg"){
    C <- 181 # neg cor, was 181
  } else if (cor.set == "pos"){
    C <- 0 # pos cor
  } else if (cor.set == "zero"){
    C <- 91 # zero cor, was 91 = 365/4
  }
  
  et = A * cos(B*(x-C)) + D + rnorm(length(x), mean=0, sd = 0.5)
  et[et < 0] <- 0
  
  #plot(rain, type="l", col=4)
  #lines(et, col=2)
  
  #   ## (fit model for simulated flow)
  #   
  #   flow = filter(0.2*((rain*0.3) - (et*0.1)), filter=c(0.7,0.2), method="recursive")
  #   flow[flow < 0] <- 0
  #   #plot(flow, type="l")
  #   dates <- seq(as.Date("2001/01/01"), as.Date("2003/12/31"), by="days")
  #   dat <- data.frame(date = dates, P = rain, E = et, Q = flow)
  #   dat <- zoo(dat[,c("P", "E", "Q")], order.by = dat$date, frequency = 1)
  #   
  #   mod0 <- hydromad(dat, sma = "bucket", routing = "expuh", delay = 1, tau_s = 0.02, tau_q = 8, v_s = 0.8)
  #   
  #   fit0 <- fitByOptim(mod0, samples = 100, method="PORT", sampletype = c("latin.hypercube"))
  #   #xyplot(fit0)
  #   
  #   #fit0pos <- unlist(fit0$parlist)
  #   #fit0neg <- unlist(fit0$parlist)
  #   #fit0zero <- unlist(fit0$parlist)
  #   # colMeans(rbind(fit0pos, fit0neg, fit0zero)) # used for params below
  #   
  #   # sims <- predict(mod0, newdata = dat[,-3])
  #   #plot(as.numeric(sims[-c(1:50, (length(sims)-49):length(sims))]), type="l")
  #   #lines(flow[-c(1:50, (length(sims)-49):length(sims))], col=2)
  #   
  #   ##
  
  dat <- data.frame( date = 1:length(rain), P = rain, E = et, year = rep(1:yrs,each=365), season = rep(c(rep(1,s1.days/2), rep(2,s2.days), rep(1,s1.days/2)), yrs) )
  dat <- zoo(dat[,c("P", "E", "year", "season")], order.by = dat$date, frequency = 1)
#   
#   # pre-fitted model (based on fit to AR filtered data above)
#   mod0 <- hydromad(dat, sma = "bucket", routing = "expuh",
#                    # routing
#                    delay = 1, tau_s = 0.02, tau_q = 8, v_s = 0.8,
#                    # sma
#                    Sb = 600, fc = 0.04, a.ei = 0.3, M = 0.2, a.ss = 0.09
#   )
#
#   # gw.frac (would need to enter into model inputs) is percentage of daily rainfall that infiltrates to groundwater aquifer (gw.frac = 0.05)
#
#   dat$Q <- predict(mod0) - (gw.frac * predict(mod0)) # streamflow
#   dat$A <- gw.frac * dat$Q # aquifer
#   # xyplot(dat)
  
#   dat <- as.data.frame(dat[-c(1:365), ], row.names = 1:(365*(yrs-1))) # trim edges
  
  return (dat)
}

## simulates Q
f.Q <- function(dat){
  
#   dat is PE input as a zoo object
#   Q.S.frac is fraction of simulated Q partitioned to stream: environmental (min) flow
#   Q.A.frac is raction of simulated Q that infiltrates to aquifer
#   assuming everything in mm
  
  # arbitrary model
    mod0 <- hydromad(dat[,c("P","E")], sma = "bucket", routing = "expuh",
                     # routing
                     delay = 1, tau_s = 0.02, tau_q = 8, v_s = 0.8,
                     # sma
                     Sb = 600, fc = 0.04, a.ei = 0.3, M = 0.2, a.ss = 0.09
    )
  
    dat$Q <- predict(mod0)
#     Qtot <- predict(mod0)
#     dat$Q.A <- Q.A.frac*Qtot
#     dat$Q.S <- Q.S.frac*Qtot
#     dat$Q.R <- (1-Q.A.frac-Q.S.frac)*Qtot
#     
#     if (!all.equal( rowSums(cbind(dat$Q.A,dat$Q.S,dat$Q.R)), as.numeric(Qtot))) stop(paste0("Sum of flow partitions is not equal to total flow in iteration: ",l,"."))
#     
    # re-order
    # dat <- dat[,c("P", "E", "Q.A", "Q.S", "Q.R", "year", "season")]
    dat <- dat[,c("P", "E", "Q", "year", "season")]
    # xyplot(dat)
  
    #dat <- as.data.frame(dat[-c(1:365), ], row.names = 1:(365*(yrs-1))) # trim edges
  
  return (dat)
}

## Creates a daily vector of seasonally-summed variables (according to season indicators in supplied data)
f.season.vec <- function(dat, Qd.val){
  
  # Qd.val is a vector of length two w/ names == "s1" and "s2" of daily values
  # dat is data frame with column labeled "season"; uses only values from year = l
  
  if (any(names(Qd.val) != c("s1", "s2"))) stop ("Names of input Qd.val to f.season.vec() function must be c('s1', 's2').")
  
  df.season <- data.frame(season = as.numeric(dat[,"season"])) # season indicator vector
  s.val <- rep(NA, nrow(df.season)) # empty seasonal value vector
  s.val[which(df.season == 1)] <- Qd.val["s1"]
  s.val[which(df.season == 2)] <- Qd.val["s2"]
  
  return(s.val)
}

###########################################################################
###########################################################################
###########################################################################
#### Submodel Functions ####
###########################################################################
###########################################################################
###########################################################################

#############################################################################
########### Social ###########
#############################################################################

submodel.S <- function(x, d, l){
  
  ## !! ## d not currently used in this function, but is included for future modifications - e.g. rules or functions based on climate (although this would most likely end up in the "modifications" submodel, not here.)
  
  #############################################################################
  ########### Group populations and group switching ###########
  #############################################################################
  
  if (l == 2){ ## if first iteration, use initialization values as is.
    
    x$out$S[(l),"X",] <- x$out$S[(l-1),"X",]
    x$out$S[(l),"N",] <- x$out$S[(l-1),"N",]
    
  } else { ## calculate popoulation changes
    
    # average population profit (not recorded)
    Pi.avg <- (x$out$S[(l-1),"Pi",] %*% x$out$S[(l-1),"N",])/sum(x$out$S[(l-1),"N",])
    
    # changes in group proportions
    
    if (x$params.pop.S$YN.delta.X == "Y"){
      # if group switching is desired
      delta.X <- x$out$S[(l-1),"X",] * ( ( x$out$S[(l-1),"Pi", ] / Pi.avg) - 1 )
      
      if (!all.equal( as.numeric(delta.X["a"]) , as.numeric(-delta.X["u"]) ) ) stop("Proportional changes in groups are not matched between groups.")
      
      # population change (proportion) limit
      delta.X <- sign(delta.X) * apply( rbind(abs(delta.X), rep(x$params.pop.S$delta.X.limit,2)), 2, min, na.rm=T) # if limit is NA, it is ignored and returns the same delta.X
      
    } else if (x$params.pop.S$YN.delta.X == "N"){
      # if group switching is not desired (static group populations)
      delta.X <- 0
    }
    
    # enter group proportions in current time step
    x$out$S[(l),"X",] <- x$out$S[(l-1),"X",] + delta.X
    
    if (!all.equal(sum(x$out$S[(l),"X",]), 1)) stop("Proportions of population in each group do not sum to 1")
    
    # enter group populations in current time step
    x$out$S[(l),"N",] <- x$out$S[(l),"X",] * sum(x$out$S[(l-1),"N",])
    if (!all.equal( sum(x$out$S[(l),"N",]), sum(x$out$S[(l-1),"N",]) )) stop("Total population count has changed.")
  }
  
  #############################################################################
  ########### Urban ###########
  #############################################################################
  
  ########### Water availability subject to physical constraints ###########
  
  # available water (from previous time step)
  Qsw.avail.u <- min(x$out$H[(l-1),"Qsw.store.H"], x$params.water.S$Qsw.limit.u)
  Qgw.avail.u <- min(x$out$H[(l-1),"Qgw.store.H"], x$params.water.S$Qgw.limit.u)
  Qw.avail.u <- Qsw.avail.u + Qgw.avail.u
  
  # accessible surface water
  if ( Qsw.avail.u >= (x$params.water.S$Qsw.fraction.u * Qw.avail.u) ) {
    x$out$S[l,"Qsw.access","u"] <- x$params.water.S$Qsw.fraction.u * Qw.avail.u
  } else if ( Qsw.avail.u >= 0 & Qsw.avail.u < (x$params.water.S$Qsw.fraction.u * Qw.avail.u) ){
    x$out$S[l,"Qsw.access","u"] <- Qsw.avail.u
  } else {
    x$out$S[l,"Qsw.access","u"] <- NA
  }
  
  if (is.na(x$out$S[l,"Qsw.access","u"])) stop("Calculation of accessible surface water failed.")
  
  # accessible total water
  Qw.access.u <- x$out$S[l,"Qsw.access","u"]/x$params.water.S$Qsw.fraction.u
  
  # accessible groundwater
  if ( Qgw.avail.u >= ((1-x$params.water.S$Qsw.fraction.u)*Qw.access.u) ) {
    x$out$S[l,"Qgw.access","u"] <- (1-x$params.water.S$Qsw.fraction.u)*Qw.access.u
  } else if ( Qgw.avail.u >= 0 & Qgw.avail.u < ((1-x$params.water.S$Qsw.fraction.u)*Qw.access.u) ) {
    x$out$S[l,"Qgw.access","u"] <- Qgw.avail.u
  } else {
    x$out$S[l,"Qgw.access","u"] <- NA
  }
  
  if (is.na(x$out$S[l,"Qgw.access","u"])) stop("Calculation of accessible groundwater failed.")
  
  # check
  if (Qw.access.u != x$out$S[l,"Qsw.access","u"] + x$out$S[l,"Qgw.access","u"] ) stop("Sum of accessible surface water and groundwater does not equal total accessible water (for the urban group).")
  
  # Accessible water (sw + gw) for individuals = (individual level) constraint
  qw.access.u <- Qw.access.u / x$out$S[l,"N","u"]
  
  ########### Price ###########
  
  # price for hours in current time step is a function of total quantity of hours worked in previous time step. If quantity produed (Q.i) in last time step is less than 1, the price is set as the maximum price.
  x$out$S[l,"p.i","u"] <- ifelse( (x$out$S[(l-1),"Q.i","u"] < 1), ((x$params.price.S$l.u * x$params.price.S$m.u) + x$params.price.S$k.u), (x$params.price.S$k.u + (x$params.price.S$l.u * (x$params.price.S$m.u / x$out$S[(l-1),"Q.i","u"])) ) )
  
  ## Visulize price function over different quantities (of hours)
  # k.u_test <- x$params.price.S$k.u
  # l.u_test <- x$params.price.S$l.u
  # m.u_test <- x$params.price.S$m.u
  # Qh.u_test <- 1:100
  # ph_test <- k.u_test + (l.u_test * (m.u_test/Qh.u_test))
  # plot(Qh.u_test, ph_test, type="l")
  
  ########### Optimization ###########
  
  # objective/profit function is f.pi.u() and gradient function is f.pi.u.grad() (see `model_functions.R')
  
  # if population is non-zero, and water access > 0, proceed with optimization; else, set water use to zero
  
  if (x$out$S[l,"N","u"] > 1e-05 & qw.access.u > 1e-05){ # if non-zero population in urban and water available > 0
    
    ## constraints: 
    # note: qw.used.u <= qw.access.u and qw.used.u >= 0 (the second is duplicative since qw.access.u >= 0 always); form is (- qw.used.u >= - qw.access.u) -> (-x1 >= -c1).
    ui_mat <-  rbind(c(-1)) # constraint matrix (k x p); dim=1x1
    ci_vec <- c(-qw.access.u) # constraint vector (k); dim=1
    theta_vec <- (solve(ui_mat) %*% ci_vec ) - 0.00001 # starting value; ui %*% theta - ci >= 0 defines the feasible region for starting values (theta).
    if (!(ui_mat %*% theta_vec - ci_vec >= 0)) stop("Optimization starting values (`theta') are not in the feasible region for group 'u'. See ?constrOptim.")
    
    ## profit maximization
    sols <- constrOptim(theta = theta_vec, f = f.pi.u, ui = ui_mat, ci = ci_vec, control = list(fnscale = -1), grad=f.pi.u.grad, x=x, l=l)
    # control = list(fnscale = -1) is maximization
    # method = "BFGS" default; one-dimensional fails with "Nelder-Mead"
    # if you do not want to supply gradient, use "Nelder-Mead" method and for one-dimensional, just use optimize(), see commented out section below - "Nelder-Mead" will fail on one-dimensional.
    
    # check convergence
    if (!(sols$convergence == 0)) paste0("Optimization did not converge on iteration: ",l,". See error code: ",sols$convergence, " in ?optim.")
    
    #   ## single variable optim - same answer as above constrOptim()
    #   sols <- optimize(f = f.pi.u, lower = 0, upper = qw.access.u, maximum = TRUE, x=x)
    #   # sols$maximum # optimal water use
    #   # sols$objective # profit, given optimal use
    #   names(sols) <- c("par", "value") # rename to match output to constrOptim()
    
    ## returns (used in commands below)
    qw.used.u <- as.numeric(sols$par) # optimal water use by inidividual
    pi.u <- as.numeric(sols$value) # individual profit, given optimal use
    
  } else { # zero population in urban or zero water available
    
    ## returns (used in commands below)
    qw.used.u <- 0 # optimal water use by inidividual
    pi.u <- 0 # individual profit, given optimal use
    
  }
  
  ########### Water Use ###########
  
  # water use by group
  Qw.used.u <- qw.used.u * x$out$S[l,"N","u"]
  x$out$S[l,"Qsw.use.S","u"] <- x$params.water.S$Qsw.fraction.u * Qw.used.u
  x$out$S[l,"Qgw.use.S","u"] <- (1-x$params.water.S$Qsw.fraction.u) * Qw.used.u
  ## !! ## see footnote about there potentially not beign enough groundwater to fulfil this - and make if/else statement to cover this.
  
  ########### Production ###########
  
  # quantity of product (hours) produced (production function)
  qh <- (ifelse(x$out$S[l,"N","u"] > 1e-05,1,0) * x$params.prod.S$hm.u) + (x$params.prod.S$alpha.u * (qw.used.u^x$params.prod.S$beta.w.u)) # individual production
  x$out$S[l,"Q.i","u"] <- qh * x$out$S[l,"N","u"] # total production
  
  ########### Profit ###########
  
  # total group profit (individual profit * number of individuals in group)
  x$out$S[l,"Pi","u"] <- pi.u * x$out$S[l,"N","u"]
  
  ########### Checks ###########
  
  # head(x$out$S[,,"u"])
  
  # check current time step values are completed
  if (any(is.na(x$out$S[l,,"u"]))) stop(paste0("There are NA values in social model output remaining at the end of iteration: ", l, " for group 'u'."))
  
  # optimal price is same as recorded price
  # ph_optimal = (pi.u + (x$params.profit.S$pw.u * qw.used.u))/qh
  
  # profit calculation (individual profit * number of individuals) is same as when calculated using straightforward formula for profit function with optimal values.
  # x$out$S[l,"Pi","u"] == (x$out$S[l,"p.i","u"] * qh - (x$params.profit.S$pw.u * qw.used.u))* x$out$S[l,"N","u"]
  
  ## Revenue: profit + costs
  # Note: in out: Q.i * p.i is total revenue, and will not equal Pi.
  # x$out$S[l,"p.i","u"]*x$out$S[l,"Q.i","u"] # == x$out$S[l,"p.i","u"] * qh * x$out$S[l,"N","u"]
  
  #############################################################################
  ########### Agricultural ###########
  #############################################################################
  
  ########### Water availability subject to physical constraints ###########
  
  # available water (from previous time step and after urban group use)
  Qsw.avail.a <- min(x$out$H[(l-1),"Qsw.store.H"] - x$out$S[l,"Qsw.use.S","u"], x$params.water.S$Qsw.limit.a)
  Qgw.avail.a <- min(x$out$H[(l-1),"Qgw.store.H"] - x$out$S[l,"Qgw.use.S","u"], x$params.water.S$Qgw.limit.a)
  Qw.avail.a <- Qsw.avail.a + Qgw.avail.a
  
  # accessible water (equal to available water in agriculture)
  x$out$S[l,"Qsw.access","a"] <- Qsw.avail.a
  x$out$S[l,"Qgw.access","a"] <- Qgw.avail.a
  Qw.access.a <- Qw.avail.a
  
  # check
  if (!all.equal(Qw.access.a, (x$out$S[l,"Qsw.access","a"] + x$out$S[l,"Qgw.access","a"])) ) stop("Sum of accessible surface water and groundwater does not equal total accessible water (for the agricultural group).")
  
  # Accessible water (sw + gw) for individuals = (individual level) constraint
  qw.access.a <- Qw.access.a / x$out$S[l,"N","a"]
  qsw.access.a <- x$out$S[l,"Qsw.access","a"] / x$out$S[l,"N","a"]
  qgw.access.a <- x$out$S[l,"Qgw.access","a"] / x$out$S[l,"N","a"]
  
  ########### Price ###########
  
  # price for crops in current time step is a function of total quantity of crops produced in previous time step. If the quantity produce (Q.i) in previous time step is less than 1, then the price is set to the maximum price: (l*m + k)
  
  x$out$S[l,"p.i","a"] <- ifelse( (x$out$S[(l-1),"Q.i","a"] < 1), ((x$params.price.S$l.a * x$params.price.S$m.a) + x$params.price.S$k.a), (x$params.price.S$k.a + (x$params.price.S$l.a * (x$params.price.S$m.a / x$out$S[(l-1),"Q.i","a"]))) )
  
  ## Visulize price function over different quantities (of hours)
  # k.a_test <- x$params.price.S$k.a
  # l.a_test <- x$params.price.S$l.a
  # m.a_test <- x$params.price.S$m.a
  # Qc.a_test <- 1:100
  # pc_test <- k.a_test + (l.a_test * (m.a_test/Qc.a_test))
  # plot(Qc.a_test, pc_test, type="l")
  
  ########### Optimization ###########
  
  ## objective/profit function is f.pi.a() and gradient function is f.pi.a.grad() (see `model_functions.R')
  
  # if population is non-zero, proceed with optimization; else, set water use to zero
  if (x$out$S[l,"N","a"] > 1e-05){ # if non-zero population in agriculture
    
    # if sw and gw are > 0, proceed with two-dimensional optimization; else, proceed with one-dimensional optimization of the source that has water
    
    ## !! ## it seems like the two-dimensional optimization should work, even if one of the variables is set to zero - potentially I have to set up the constraints explicity for this, which I have not done, but should consider doing later; right now, I just have individual functions written out for each case: two-dimensional for gw + sw > 0, and two one-dimensional profit and gradient functions for sw and gw - which are used according the if statements below.
    
    if (qsw.access.a > 1e-05 & qgw.access.a > 1e-05){ # if sw and gw are both > 0
      
      ## constraints: 
      # note: qsw.used.a <= qsw.access.a and qsw.used.u >= 0 (the second is duplicative since qsw.access.a >= 0 always; same for gw); form is (- qsw.used.a >= - qsw.access.a) -> (-x1 >= -c1 and -x2 >= -c2).
      ui_mat <-  rbind(c(-1,0), c(0,-1)) # constraint matrix (k x p); dim=2x2
      ci_vec <- c(-qsw.access.a, -qgw.access.a) # constraint vector (k); dim=1
      theta_vec <- (solve(ui_mat) %*% ci_vec ) - 0.00001 # starting value; ui %*% theta - ci >= 0 defines the feasible region for starting values (theta).
      if (!all(ui_mat %*% theta_vec - ci_vec >= 0)) stop("Optimization starting values (`theta') are not in the feasible region for group 'a' (optimization with both SW and GW). See ?constrOptim.")
      
      ## profit maximization
      sols <- constrOptim(theta = theta_vec, f = f.pi.a, ui = ui_mat, ci = ci_vec, control = list(fnscale = -1), grad=f.pi.a.grad, x=x, l=l)
      # control = list(fnscale = -1) is maximization
      # method = "BFGS" default
      # if you do not want to supply gradient, use "Nelder-Mead" method
      
      # check convergence
      if (!(sols$convergence == 0)) paste0("Optimization did not converge on iteration: ",l,". See error code: ",sols$convergence, " in ?optim.")
      
      ## returns (used in commands below)
      qsw.used.a <- as.numeric(sols$par)[1] # optimal surface water use by inidividual
      qgw.used.a <- as.numeric(sols$par)[2] # optimal groundwater use by inidividual
      pi.a <- as.numeric(sols$value) # individual profit, given optimal use
    
    } else if(qsw.access.a < 1e-05 & qgw.access.a > 1e-05){ # sw < 0; only gw available
      
      ## constraints: 
      ui_mat <-  rbind(c(-1)) # constraint matrix (k x p); dim=1x1
      ci_vec <- c(-qgw.access.a) # constraint vector (k); dim=1
      theta_vec <- (solve(ui_mat) %*% ci_vec ) - 0.00001 # starting value; ui %*% theta - ci >= 0 defines the feasible region for starting values (theta).
      if (!(ui_mat %*% theta_vec - ci_vec >= 0)) stop("Optimization starting values (`theta') are not in the feasible region for group 'a' (optimization with only GW). See ?constrOptim.")
      
      ## profit maximization
      sols <- constrOptim(theta = theta_vec, f = f.pi.a_gw, ui = ui_mat, ci = ci_vec, control = list(fnscale = -1), grad=f.pi.a.grad_gw, x=x, l=l)
      # control = list(fnscale = -1) is maximization
      # method = "BFGS" default; one-dimensional fails with "Nelder-Mead"
      # if you do not want to supply gradient, use "Nelder-Mead" method and for one-dimensional, just use optimize(), see commented out section below - "Nelder-Mead" will fail on one-dimensional.
      
      # check convergence
      if (!(sols$convergence == 0)) paste0("Optimization did not converge on iteration: ",l,". See error code: ",sols$convergence, " in ?optim.")
      
      ## returns (used in commands below)
      qsw.used.a <- 0 # optimal surface water use by inidividual
      qgw.used.a <- as.numeric(sols$par) # optimal groundwater use by inidividual
      pi.a <- as.numeric(sols$value) # individual profit, given optimal use
      
    } else if (qsw.access.a > 1e-05 & qgw.access.a < 1e-05){ # gw < 0; only sw available
      
      ## constraints: 
      ui_mat <-  rbind(c(-1)) # constraint matrix (k x p); dim=1x1
      ci_vec <- c(-qsw.access.a) # constraint vector (k); dim=1
      theta_vec <- (solve(ui_mat) %*% ci_vec ) - 0.00001 # starting value; ui %*% theta - ci >= 0 defines the feasible region for starting values (theta).
      if (!(ui_mat %*% theta_vec - ci_vec >= 0)) stop("Optimization starting values (`theta') are not in the feasible region for group 'a' (optimization with only SW). See ?constrOptim.")
      
      ## profit maximization
      sols <- constrOptim(theta = theta_vec, f = f.pi.a_sw, ui = ui_mat, ci = ci_vec, control = list(fnscale = -1), grad=f.pi.a.grad_sw, x=x, l=l)
      # control = list(fnscale = -1) is maximization
      # method = "BFGS" default; one-dimensional fails with "Nelder-Mead"
      # if you do not want to supply gradient, use "Nelder-Mead" method and for one-dimensional, just use optimize(), see commented out section below - "Nelder-Mead" will fail on one-dimensional.
      
      # check convergence
      if (!(sols$convergence == 0)) paste0("Optimization did not converge on iteration: ",l,". See error code: ",sols$convergence, " in ?optim.")
      
      ## returns (used in commands below)
      qsw.used.a <- as.numeric(sols$par) # optimal surface water use by inidividual
      qgw.used.a <- 0 # optimal groundwater use by inidividual
      pi.a <- as.numeric(sols$value) # individual profit, given optimal use
      
    }
    
  } else { # zero population in agriculture
    
    ## returns (used in commands below)
    qsw.used.a <- 0 # optimal surface water use by inidividual
    qgw.used.a <- 0 # optimal groundwater use by inidividual
    pi.a <- 0 # individual profit, given optimal use
    
  }
    
  ########### Water Use ###########
  
  # water use by group
  x$out$S[l,"Qsw.use.S","a"] <- qsw.used.a * x$out$S[l,"N","a"]
  x$out$S[l,"Qgw.use.S","a"] <- qgw.used.a * x$out$S[l,"N","a"]
  
  ########### Production ###########
  
  # quantity of product (hours) produced (production function)
  qc <- x$params.prod.S$alpha.a * (qsw.used.a^x$params.prod.S$beta.sw.a) *(qgw.used.a^x$params.prod.S$beta.gw.a) # individual production
  x$out$S[l,"Q.i","a"] <- qc * x$out$S[l,"N","a"] # total production
  
  ########### Profit ###########
  
  # total group profit (individual profit * number of individuals in group)
  x$out$S[l,"Pi","a"] <- pi.a * x$out$S[l,"N","a"]
  
  ########### Checks ###########
  
  # head(x$out$S[,,"a"])
  
  # check current time step values are completed
  if (any(is.na(x$out$S[l,,"a"]))) stop(paste0("There are NA values in social model output remaining at the end of iteration: ", l, " for group 'a'."))
  
  # optimal price is same as recorded price
  # pc_optimal = (pi.a + (x$params.profit.S$psw.a * qsw.used.a) + (x$params.profit.S$pgw.a * qgw.used.a))/qc
  
  # profit calculation (individual profit * number of individuals) is same as when calculated using straightforward formula for profit function with optimal values.
  # x$out$S[l,"Pi","a"] == (x$out$S[l,"p.i","a"] * qc - (x$params.profit.S$psw.a * qsw.used.a) - (x$params.profit.S$pgw.a * qgw.used.a) ) * x$out$S[l,"N","a"]
  
  ## Revenue: profit + costs
  # Note: in out: Q.i * p.i is total revenue, and will not equal Pi.
  # x$out$S[l,"p.i","a"]*x$out$S[l,"Q.i","a"] # == x$out$S[l,"p.i","a"] * qc * x$out$S[l,"N","a"]
  
  ########### Model Ouput ###########
  
  return(list(x = x, d = d))
  
}

#############################################################################
########### Hydrology ###########
#############################################################################

#### NOTES ####

# this section also indexed by l

# demand in this step = use. Demand was constratined with the previous iteration's (year's) water stores (Qsw.store.H and Qgw.store.H).

# note: demand/use (Qd) fractions are the same across sources - sources don't affect how they filter thorugh the soil (in a) or through infrastructure (in u) and return to environment in form of E,A,Q

## params: see x$params.H

## agriculture
# Qd.frac.E.a == fraction of water used (applied to crops) in a that goes to E
# the remaining amount goes to P

## urban
# Qd.frac.E.u == fraction of APPLIED water (total water minus wastewater) that goes to ET

# Qd.frac.W.u == fraction of used water that goes to wastewater; exits the system (could add in a fraction of this that is recycled ## !! ##)

## seasonal use fractions
# s1.use.frac.a == 1st season (wet) fraction of water use (sw + gw) - agriculture
# s1.use.frac.u == # 1st season (wet) fraction of water use (sw + gw) - urban

## flow partitioning
# Q.A.frac == fraction of flows that infiltrate to groundwater aquifer
# Q.S.frac == fraction of flows for environment (minimum flows)
# remaining amount goes to sw reservoir

submodel.H <- function(x, d, l){
  
  #### Adjustments ####
  
  # agriculture
  Qd.a <- x$out$S[l,"Qsw.use.S","a"] + x$out$S[l,"Qgw.use.S","a"] # total water used
  Qd.2E.a <- Qd.a * x$params.H$Qd.frac.E.a # water used to E
  Qd.2P.a <- Qd.a - Qd.2E.a # water remaining ~ applied water; same as rainfall (a water input)
  
  if (sum(Qd.2E.a,Qd.2P.a) != Qd.a) stop(paste0("Sum of water demand partitions do not equal total water demand in group 'a', iteration:", l,"."))
  
  # urban
  Qd.u <- x$out$S[l,"Qsw.use.S","u"] + x$out$S[l,"Qgw.use.S","u"]
  Qd.2W.u <- Qd.u * x$params.H$Qd.frac.W.u # water used to wastewater (exit system)
  Qd.u_applied <- Qd.u - Qd.2W.u # applied water; remaining for i.e. landscaping
  Qd.2E.u <- Qd.u_applied * x$params.H$Qd.frac.E.u # applied water used to E
  Qd.2P.u <- Qd.u - Qd.2E.u - Qd.2W.u # water remaining ~ same as rainfall (a water input)
  
  if (!all.equal(sum(Qd.2E.u,Qd.2W.u,Qd.2P.u), Qd.u)) stop(paste0("Sum of water demand partitions do not equal total water demand in group 'u', iteration:", l,"."))
  
  #### Daily Values ####
  ## calculate daily values, varied according to season (fractions of use by season)
  
  Qd.2P.a_daily <- c(s1 = (Qd.2P.a*x$params.H$s1.use.frac.a), s2 = (Qd.2P.a*(1-x$params.H$s1.use.frac.a)))/365
  Qd.2P.u_daily <- c(s1 = (Qd.2P.u*x$params.H$s1.use.frac.u), s2 = (Qd.2P.u*(1-x$params.H$s1.use.frac.u)))/365
  
  Qd.2E.a_daily <- c(s1 = (Qd.2E.a*x$params.H$s1.use.frac.a), s2 = (Qd.2E.a*(1-x$params.H$s1.use.frac.a)))/365
  Qd.2E.u_daily <- c(s1 = (Qd.2E.u*x$params.H$s1.use.frac.u), s2 = (Qd.2E.u*(1-x$params.H$s1.use.frac.u)))/365
  
  Qd.2W.u_daily <- c(s1 = (Qd.2W.u*x$params.H$s1.use.frac.a), s2 = (Qd.2W.u*(1-x$params.H$s1.use.frac.a)))/365
  # no wastewater for agriculture
  
  #### Modify PE data ####
  ## add demand/use daily value returns to PE data
  
  dfi <- which(d$year == l) # index of days in year l
  
  # P
  d[dfi,"P"] <- d[dfi,"P"] + f.season.vec(d[dfi,], Qd.2P.a_daily) + f.season.vec(d[dfi,], Qd.2P.u_daily)
  
  # E
  d[dfi,"E"] <- d[dfi,"E"] + f.season.vec(d[dfi,], Qd.2E.a_daily) + f.season.vec(d[dfi,], Qd.2E.u_daily)
  
  #### Simulate and Partition Flow ####
  # note: re-simulates all years flows from start to finish, with adjusted PE data that is updated at each year based on water use.
  
  df.PEQ <- f.Q(d)
  # xyplot(df.PEQ)
  
  # subtract wastewater from total flows
  df.PEQ[dfi,"Q"] <- df.PEQ[dfi,"Q"] - f.season.vec(df.PEQ[dfi,],Qd.2W.u_daily) # just replace current year's Q values; all else remain the same
  
  df.PEQ$Q.A <- x$params.H$Q.A.frac * df.PEQ$Q
  df.PEQ$Q.S <- x$params.H$Q.S.frac * df.PEQ$Q
  df.PEQ$Q.R <- (1- x$params.H$Q.A.frac - x$params.H$Q.S.frac) * df.PEQ$Q
  
  if (!all.equal( rowSums(cbind(df.PEQ$Q.A,df.PEQ$Q.S,df.PEQ$Q.R)), as.numeric(df.PEQ$Q))) stop(paste0("Sum of flow partitions is not equal to total flow in iteration: ",l,"."))
  
  # re order
  df.PEQ <- df.PEQ[,c("P", "E", "Q", "Q.R", "Q.A", "Q.S", "year", "season")]
  # Q.R is flow (depth/day) to reservoir (to Qsw.store)
  # Q.A is flow (depth/day) to aquifer (to Qgw.store)
  # Q.S is flow (depth/day) left to stream = environmental flows
  
  #### Format H output ####
  ## Sum variables over year, add return flows to A and R stores, and subtract demand/use from stores
  
  # sw
  x$out$H[l, "Qsw.store.H"] <- x$out$H[(l-1), "Qsw.store.H"] + sum(df.PEQ[dfi,"Q.R"]) - x$out$S[l,"Qsw.use.S","a"] - x$out$S[l,"Qsw.use.S","u"]
  
  # gw
  x$out$H[l, "Qgw.store.H"] <- x$out$H[(l-1), "Qgw.store.H"] + sum(df.PEQ[dfi,"Q.A"]) - x$out$S[l,"Qgw.use.S","a"] - x$out$S[l,"Qgw.use.S","u"]
  
  # stream (environment)
  x$out$H[l, "Qstream.H"] <- sum(df.PEQ[dfi,"Q.S"])
  
  # wastewater
  x$out$H[l, "Qwaste.H"] <- Qd.2W.u
  
  # use (applied sw and gw in both groups)
  x$out$H[l, "Quse.H"] <- (x$out$S[l,"Qsw.use.S","a"] + x$out$S[l,"Qgw.use.S","a"] + x$out$S[l,"Qsw.use.S","u"] + x$out$S[l,"Qgw.use.S","u"])
  
  # P (minus applied water = actual rainfall)
  x$out$H[l, "P"] <- sum(df.PEQ[dfi,"P"]) - x$out$H[l, "Quse.H"]
  
  # E
  x$out$H[l, "E"] <- sum(df.PEQ[dfi,"E"])
  
  # returns full P,E,Q data set on last iteration, otherwise passes just P,E data so that it can be used in the next iteration.
  # note: data values for P are not split into P and use (applied) becuase it might potentially affect flow calculations at each step (via soil moisture).
  
  if (l == n){
    datareturn <- df.PEQ
  } else {
    # P is NOT adjusted
    datareturn <- d
  }
  
  #### Model Outputs ####
  return(list(x = x, d = datareturn))
}

#############################################################################
########### Modifications ###########
#############################################################################

## !! ## e.g... list of if/else statements that change parameters - e.g. water use limits, subject to some decision; this is where `management' would come in. A management decision would be represented by there existing a trigger that changes parameters.

submodel.M <- function(x, d, l){
  
  ## !! ## To be completed.
  
  # if Qsw.slope.trigger_psw.a_percent is entered into input (is not NULL)
  if (!is.null(x$params.M$Qsw.slope.trigger_psw.a_percent)){
    # if slope of SW store is negative
    if ( sign(lm(x$out$H[1:l,"Qsw.store.H"] ~ c(1:l))$coefficients[2]) < 0 ){
      # then, increase price of sw by the designated percent
      x$params.profit.S$psw.a <- x$params.profit.S$psw.a + (x$params.profit.S$psw.a * x$params.M$Qsw.slope.trigger_psw.a_percent)
      x$out$log[l,"Notes"] <- ifelse( is.na(x$out$log[l,"Notes"]),  paste0("Qsw.slope.trigger_psw.a_percent: Slope of the Qsw.store was negative over all previous time steps, psw.a increased by: ",x$params.M$Qsw.slope.trigger_psw.a_percent, " to: ", round(x$params.profit.S$psw.a,3)), paste0(x$out$log[l,"Notes"], " // " , paste0("Qsw.slope.trigger_psw.a_percent: Slope of the Qsw.store was negative over all previous time steps, psw.a increased by: ",x$params.M$Qsw.slope.trigger_psw.a_percent, " to: ", round(x$params.profit.S$psw.a,3)) )) 
    }
  }
  
  # if Qgw.slope.trigger_pgw.a_percent is entered into input (is not NULL)
  if (!is.null(x$params.M$Qgw.slope.trigger_pgw.a_percent)){
    # if slope of gw store is negative
    if ( sign(lm(x$out$H[1:l,"Qgw.store.H"] ~ c(1:l))$coefficients[2]) < 0 ){
      # then, increase price of gw by the designated percent
      x$params.profit.S$pgw.a <- x$params.profit.S$pgw.a + (x$params.profit.S$pgw.a * x$params.M$Qgw.slope.trigger_pgw.a_percent)
      x$out$log[l,"Notes"] <- ifelse( is.na(x$out$log[l,"Notes"]),  paste0("Qgw.slope.trigger_pgw.a_percent: Slope of the Qgw.store was negative over all previous time steps, pgw.a increased by: ",x$params.M$Qgw.slope.trigger_pgw.a_percent, " to: ", round(x$params.profit.S$pgw.a,3)), paste0(x$out$log[l,"Notes"], " // " , paste0("Qgw.slope.trigger_pgw.a_percent: Slope of the Qgw.store was negative over all previous time steps, pgw.a increased by: ",x$params.M$Qgw.slope.trigger_pgw.a_percent, " to: ", round(x$params.profit.S$pgw.a,3)) )) 
      
    }
  }
  
  #### Model Outputs ####
  return(list(x = x, d = d))
}

###########################################################################
###########################################################################
###########################################################################
#### Master Model Functions ####
###########################################################################
###########################################################################
###########################################################################

# x is initialization object above
# n is number of iterations (years)
# d is data, which is PE data only (currently ## !! ##)

shmod <- function(x,d,n){
  
  ## checks
  if (!((length(unique(d$year)) == n) & (nrow(d)/365 == n))) stop("Length of data time series in years must match the number of iterations 'n'.")
  if (!(nrow(d) %% 365 == 0)) stop("Data time series must be full years with 365 days per year.")
  if (!(nrow(x$out$S) == n)) stop("In f.init(), 'iters' must match the length of data time series.")
  
  ## iterations loop
  
  for (l in 2:n){
    
    if (l == 2){
      
      cat("Iteration: ", l)
      
      SS <- submodel.S(x, d, l)
      SH <- submodel.H(SS$x, SS$d, l)
      SM <- submodel.M(SH$x, SH$d, l) ## !! ## right now: all.equal(SM,SH)
      
      ## check
      if (!any( (as.numeric(SM$d$P) != as.numeric(d$P)) )) stop(paste0("There are no differences the between initially supplied data (P,E) and the data updated by the hydrologic model, indicating an issue in either the social or hydrologic submodel at iteration: ",l,"."))
      
    } else {
      
      cat(",",l)
      
      SS <- submodel.S(SM$x, SM$d, l)
      SH <- submodel.H(SS$x, SS$d, l)
      SM <- submodel.M(SH$x, SH$d, l) ## !! ## right now: all.equal(SM,SH)
      
      ## check
      if (!any( (as.numeric(SM$d$P) != as.numeric(d$P)) )) stop(paste0("There are no differences the between initially supplied data (P,E) and the data updated by the hydrologic model, indicating an issue in either the social or hydrologic submodel at iteration: ",l,"."))
      
    }
  }
  
  ## Output plots
  
  # melt
  dfm.S <- melt(SM$x$out$S[-1,,])
  dfm.H <- melt(SM$x$out$H[-1,])
  
  # social
  SM$plot.S <- ggplot(dfm.S, aes(x=X1,y=value, colour = X3)) + geom_line() + facet_wrap(~X2, scales = "free") + theme_bw() + xlab("Years")
  # hydrology
  SM$plot.H <-ggplot(dfm.H, aes(x=X1,y=value)) + geom_line() + facet_wrap(~X2, scales = "free") + theme_bw() + xlab("Years")
  
  # data:
  # note: SM$d is modified daily climate and hydrologic data at end of run; d is original daily climate data
  
  # Split rainfall (P) and used/applied (U) water:
  SM$d$U <- SM$d$P - d$P
  SM$d$P <- SM$d$P - SM$d$U
  if (!all.equal(SM$d$P,d$P)) stop("Adjusted P data (with applied water (U) separated) does not match original P data.")
  SM$plot.d <- xyplot(SM$d)
  
  ## Returns
  
  return(SM)
  
}

