###########################################################################
###########################################################################
#### This script loads model functions, and gives examples of use.
###########################################################################
###########################################################################

###########################################################################
#### Set Up ####
###########################################################################

rm(list= ls()); # clear workspace

#### install required R packages ####

## packages on CRAN
# install.packages(c("reshape2", "ggplot2", "zoo", "latticeExtra", "polynom", "car", "Hmisc"))
## packages not on CRAN
# install.packages("hydromad", repos="http://hydromad.catchment.org", type="source")

#### load all model functions ####

## NOTE: Navigate to local workspace (uncomment below)
# master <- setwd('/Users/.../[file where model_functions.R is saved]') # mac
# master <- setwd('C:\\Users\\...\\[file where model_functions.R is saved]') # PC
# setwd(master)
# source('model_functions.R')

source('/Users/mcl/Dropbox/SESYNC/model/model_functions.R') # Morgan's computer

###########################################################################
#### Input Parameter Initialization Example ####
###########################################################################
 
## NOTE: 
# See 'model_paramaters_key.csv' for full list of parameters, their descriptions, value types and ranges, and default values.
# Parameter naming conventions: .S refers to the social submodel, .H to the hydrologic submodel, and .M to the modifications submodel; .a refers to the agricultural group and .u refers to the urban group

## default initialization of starting values and parameters;
x <- f.init()

## look at some of the initialization objects

names(x) # all init objects that can be called by appending to x$
x$out$S[,,"a"] # output table for agriculture group
names(x$params.water.S) # water parameters, i.e. use limits
x$params.water.S$Qsw.limit.a # a specific water parameter
names(x$params.price.S) # parameters to the price functions
names(x$params.H) # parameters to the hydrologic model

###########################################################################
#### Input Data Example: rainfall (P) and ET (E) ####
###########################################################################

# Simulate data; this example is the same as df.PE <- f.PE(), all values are the default values.

df.PE <- f.PE(cor = "neg", # relationship between P and E
              pmean.s1 = 8, # mean daily rainfall in season 1 (wet) 
              pmean.s2 = 3, # mean daily rainfall in season 2 (dry) 
              yrs = 10, # number of years to simulate
              s1.days = 200, # number of days in season 1 (wet)
              s2.days = 165 # number of days in season 2 (dry)
              )

xyplot(df.PE) # plot

## NOTE:
# The total number of days must be 365; season 1 is split into beginning and end of the year (years start in January).
# 'cor' is the relationship (correlation) desired between P and E; if P and E are negatively correlated ("neg"), this is similar to a climate with warmer dry seasons, such as Mediterranean and/or tropical savanna climates; positively correlated ("pos") would be for warmer wet seasons, such monsoon climates; "zero" correlation is in between.
# daily data is currently P and E only (E is assumed ET, but could be modified to be temperature, with a coefficient that converts to ET.)

###########################################################################
#### Testing ####
###########################################################################

## Example of how to set up inputs, run model, and look at outputs.

#### Inputs ####

# years
n <- 10

# simulated data
# df.PE <- f.PE(cor="neg", yrs=n, pmean.s1 = 4, pmean.s2 = 1)

# testing data (static)
df.PE <- f.PE_static(cor="neg", yrs=n, pmean = 5, etmean = 3)

# parameters
x <- f.init(iters=n,
            start.Qsw.store.H = 2,
            # start.Qgw.store.H = 10,
            # start.Qh = 50,
            # start.Qc = 1000,
            params.pop.S = list(
              N = c(100,100), 
              YN.delta.X = "N", # "Y"
              delta.X.limit = 0.25),
            params.water.S = list(
              Qsw.limit.a = 50, 
              Qgw.limit.a = 100, 
              Qsw.limit.u = 25, 
              Qgw.limit.u = 25, 
              Qsw.fraction.u = 0.8),
            params.profit.S = list(
              psw.a = 0.1, 
              pgw.a = 0.2, 
              pw.u = 0.5),
              params.M = list(
                Qsw.slope.trigger_psw.a_percent = NULL, # 0.05
                Qgw.slope.trigger_pgw.a_percent = NULL) # 0.1
            )

## NOTE:
## the parameters included above are not necessarily defaults, although some are; these are a selection of parameters and groups of parameters that are probably the most interesting to start experimenting with... and some examples of things to start fiddling with (see '#' next to entered parameters...)
## If you edit grouped parameters that are entered in list() format, include the full list, even if you're not changing anything (e.g. copy/paste from the defaults in 'model_paramaters_key.csv').

#### Run Model ####

test <- shmod(
  x = x, # parameters
  d = df.PE, # data
  n=n # years
  )

#### Look at Outputs ####

## Output objects

# output object names
names(test)
names(test$x)

# plots
test$plot.S # social output plot (annual)
test$plot.H # hydrology output plot (annual)
test$plot.d # hydrology data (daily)

# data
test$x$out # output matrices: H, S(a), S(u), log
test$x$out$H[-1, ] # hydrologic output
test$x$out$S[-1,,"a"] # social output for group "a"
test$x$out$S[-1,,"u"] # social output for group "u"
## NOTE: the first row in all outputs is the initialization data, and should be excluded from analysis.
head(test$d) ; tail(test$d) # daily hydrologic data

# record of modifications
test$x$out$log

###########################################################################
#### Extra ####
###########################################################################

#### For debugging:

# debug(shmod)
# debug(submodel.S)
# debug(submodel.H)
# 
# undebug(shmod)
# undebug(submodel.S)
# undebug(submodel.H)

# # debug (manual); go into functions manually in 'model_functions.R'
# n = 10
# df.PE <- f.PE()
# x <- f.init()
# d = df.PE
# l = 2
#
#
# ## all defaults run
# 
# rm(list= ls()); # clear workspace
# source('/Users/mcl/Dropbox/SESYNC/model_functions.R')
# 
# n = 10
# set.seed(1208)
# df.PE <- f.PE()
# x <- f.init()
# test <- shmod(x = x, d = df.PE, n=n)
