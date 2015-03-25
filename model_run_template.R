###########################################################################
#### Set Up ####
###########################################################################

rm(list= ls()); # clear workspace

#### load all model functions ####

## NOTE: Navigate to local workspace (uncomment below)
master <- setwd('/Users/.../[file where model_functions.R is saved]') # mac
#master <- setwd('C:\\Users\\...\\[file where model_functions.R is saved]') # PC
setwd(master)
source('model_functions.R')

###########################################################################
#### Testing ####
###########################################################################

#### Inputs ####

# years
n <- 10

# simulation data
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

# record of modifications
test$x$out$log