## extras
### Survival

# Survival of post-larval urchins $\textbf{p}$ is an age-specific function of fishing pressure, selectivity-at-age $q$ (assumed to be constant over time), and natural mortality, assumed to be constant across all adult age classes. Survival of age-0 urchins is a separate parameter. If fishing pressure is zero, survival will be determined solely by natural mortality.
# 
# ```{r surv}
# surv.fxn <- function(fish, q.vec, mort,surv0) {
#   p.vec <- exp(-fish*q.vec-mort)
#   p.vec[1] <- surv0
#   p.vec[length(p.vec)] <- 1/(1-surv) # survival of plus group
#   # where fish is the fishing intensity, q.vec is the selectivity vector, and mort is natural mortality
#   return(p.vec)
# }
# ## with our parameters:
# p.vec <- surv.fxn (fish=fish,q.vec=q.vec,mort=mort,surv0=surv0)
# plot(0:max.age,p.vec,main="Age vs. Survival",type='l',xlab='Age, years',ylab='Proportional Survival',
#      xaxs='i')
# ```

## Fecundity

# This function relates length to fecundity, and is a power function that can be parameterized on a species-specific basis:
#   
#   $$f_i = f_1D_i^{f_2}$$
#     
#     where $f_i$ is fecundity at age $i$, and $f_1$ and $f_2$ are parameters. For now, these parameters were chosen to given a slightly positive population growth rate.
#   
#   ```{r fertility function}
#   fert <- function(age.vec, g.params, fert1,fert2) {
#     diams <- Dt(Dinf=g.params[1],K=g.params[2],n=g.params[3],age.vec=age.vec)
#     f <- sapply(diams, function(x) fert1*x^fert2)
#     return(f)
#   }
#   
#   # with our parameters:
#   fert.vec <- fert(age.vec=ages,g.params=growth.params,fert1=fert1,fert2=fert2)
#   plot(0:max.age,fert.vec,main="Age vs. Fecundity",type='l',xlab='Age, years',ylab='Fertility (# eggs)',
#        xaxs='i')
#   ```

## Population Projection Matrix

# The density-independent model for a local population takes the form
# 
# $$\textbf{N}(\textit{t} + 1) = \textbf{MN}(\textit{t})$$
#   
#   where **M** is a Leslie matrix, and **N** is the vector of numbers-at-age.
# 
# We can now construct the Leslie matrix **M** using the age vector, fertility vector, and proportional survival vector:
#   
#   ```{r Leslie matrix density independent}
# calcM <- function (max.age,fert.vec, p.vec) {
#   M=matrix(nrow=max.age+1,ncol=max.age+1)
#   M[1,]  <- fert.vec #first row
#   for (i in 2:(max.age+1)) { #rows 2 to max age
#     row <- numeric(max.age+1)
#     row[i-1] <- p.vec[i-1]
#     M[i,] <- row
#   }
#   return(M)
# }
# 
# # test that it's working (different values for the fertility params)
# f1 <- seq(3,8,by=0.5)
# f2 <- seq(2,8,by=0.5)
# fs<-expand.grid(f1,f2)
# names(fs)=c("f1","f2")
# fs <- as.data.frame(fs)
# fs$eigen <- NA
# for(i in 1:nrow(fs)) {
#   f <- fert(age.vec=ages,g.params=growth.params,fert1=fs[i,1],fert2=fs[i,2])
#   out <- eigen(calcM(max.age=max.age,fert.vec=f,p.vec=p.vec))$value[1]
#   fs$eigen[i] <- out
# }
# 
# fert(age.vec=ages,g.params=growth.params,fert1=fert1,fert2=fert2)
# test<-calcM(max.age=max.age,fert.vec=fert.vec,p.vec=p.vec)
# eigen(test)$value
# 
# ## a vector of stable age distribution (at equilibrium) can be calculated from the right eigenvector of the M matrix
# ## we will use this later to 'seed' simulations with appropriate abundances
# stableage <- Re(eigen(calcM(max.age=max.age,fert.vec=fert.vec,p.vec=p.vec))$vector)[,1]
# stableage <- stableage/sum(stableage)
# ```
# 
# ## Fishing Mortality Vs. Yield
# 
# To show a relationship between fishing mortality and yield, we simulate the population to "equilibrium" at each level of fishing mortality.
# 
# ```{r fishing mortality vs yield}
# sim.pop <- function(nyears=200,max.age=30,surv0=5e-8, surv=0.9, mort=1-0.9,fish=0, 
#                     Dinf=9.31,growthrate=0.0926,shape=-0.3208,
#                     fert1=4.5,fert2=7, 
#                     dw.a=0.0012659, dw.b=2.7068, size.limit=8.3) {
#   ## age vector
#   ages <- 0:max.age
#   
#   ## growth parameters as a vector
#   g.params <- c(Dinf,growthrate,shape)
#   
#   ## diameters vector (size)
#   diams <- Dt(Dinf = g.params[1],K=g.params[2],n=g.params[3],age.vec=ages)
#   
#   ## fertility vector
#   fert.vec <- fert(age.vec=ages,g.params=g.params,fert1=fert1,fert2=fert2)
#   
#   ## selectivity-at-age
#   q.vec <- sel(g.params=g.params,age.vec=ages,size.limit = size.limit)
#   
#   ## survival vector
#   p.vec <- surv.fxn(fish=fish, q.vec=q.vec,mort=mort,surv0=surv0)
#   
#   ## weight-at-age
#   w.vec <- weight.at.age(age.vec=ages,g.params=g.params,dw.a=dw.a,dw.b=dw.b)
#   
#   ## stable age distribution (for seeding simulation)
#   stableage <- Re(eigen(calcM(max.age=max.age,fert.vec=fert.vec,p.vec=p.vec))$vector)[,1]
#   stableage <- stableage/sum(stableage)
#   
#   ## the simulation output matrix
#   sim <- matrix(nrow=(max.age+1), ncol=nyears)
#   simtot <- numeric() # total population
#   Ninit <- rep(250,max.age+1)*stableage # starting abundance
#   sim[,1] <- Ninit
#   simtot <- sum(Ninit)
#   eigens <- numeric() # for holding eigenvalues
#   SSB <- numeric() # for holding value of SSB in each year
#   Yield <- numeric() #for holding value of yield in each year
#   
#   ## run the simulation
#   for (i in 1:(nyears-1)) {
#     M <- calcM(max.age=max.age,fert.vec=fert.vec,p.vec=p.vec)
#     eigens[i] <- eigen(M)[[1]][1] #store the dominant eigen value
#     N.vec <- M%*%sim[,i]
#     sim[,i+1] <- N.vec
#     simtot[i+1] <- sum(N.vec)
#     SSB[i] <- ssb(w.vec=w.vec,n.vec=sim[,i])
#     Yield[i] <- yield(fish=fish,q.vec=q.vec,mort=mort,w.vec=w.vec,n.vec=sim[,i])
#   }
#   ## for the last year
#   M <- calcM(max.age=max.age,fert.vec=fert.vec,p.vec=p.vec)
#   eigens[nyears] <- eigen(M)[[1]][1] #store the dominant eigen value
#   SSB[nyears] <- ssb(w.vec=w.vec,n.vec=sim[,nyears])
#   Yield[nyears] <- yield(fish=fish,q.vec=q.vec,mort=mort,w.vec=w.vec,n.vec=sim[,nyears])
#   
#   ## Returns a list with the last 100 years of the simulation, including population (every age class), total population summed over age classes, proportional population in each age class, eigenvalues of each year's matrix, and Reproductive Value for the average individual each year, adjusted for settler survival (gamma) and self-recruitment rate a.
#   return(list(pop=sim[,(nyears-99):nyears],tot=simtot[(nyears-99):nyears],props=apply(sim,
#                                                                                       MARGIN=2,FUN=function(x) x/sum(x))[,(nyears-99):nyears], eigens=Re(eigens[(nyears-99):nyears]),
#               Yield=Yield[(nyears-99):nyears]))
# }
# simtest<- sim.pop()
# ```
