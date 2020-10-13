
# start by finding probability of not being 0 power!


alpha=1.5
beta=0.1
nu=seq(0,50,by=0.1)
gammadist<-dgamma(nu,shape=alpha,rate=beta)
plot(nu,gammadist)

# now how Do I get the probablity of getting nonzero? find valeus of area between nu=13 and 30

pgamma(13,shape=alpha,rate=beta) # comparing the cumulative between these 2
pgamma(30,shape=alpha,rate=beta)

# probability of some power is 
alpha=1.5
beta=0.1
Probability_of_any_power<-pgamma(30,shape=alpha,rate=beta)-pgamma(13,shape=alpha,rate=beta) # 0.34%

# Here we ave to integrate over possible values of nu! Begin b defing the function

powerfunction<-function(nu){
  if (nu < 13 | nu >=30){
    value=0
  }
  else if (13<=nu & nu<23){
  value=cos(0.1*pi*(nu-3))+1
  }
  else 
    value=(1031+46*nu-nu^2)/780
  return(value)
}

N=10000 # defines how many simulated to use
sims_of_nu<-rgamma(N,shape=alpha,rate=beta) # the simulated values to evaluate function at
#sum(sapply(sims_of_nu,powerfunction))/N
# combine with sapply to make function take vectors and nor just single values
plot(nu,sapply(nu,powerfunction)) # seems to work!
plot(nu,sapply(nu,powerfunction)*dgamma(nu,shape=alpha,rate=beta)) # interesting!
powerfunction<-Vectorize(powerfunction) # TO avoid using sapply

# simulate a ton of values, then sum these together
hbar<-mean(powerfunction(sims_of_nu)) # the simulated expectedvalue!
#sum((powerfunction(sims_of_nu)-hbar)^2)/N^2
var_hbar<-var(powerfunction(sims_of_nu))/(N^2)*(N-1) # The variance,
# I use built in function but then have to divide by N and multiply by N-1 to get the right formula
confint=c(hbar-1.96*sqrt(var_hbar),hbar+1.96*sqrt(var_hbar))
hbar
var_hbar
confint

h_times_f<-function(x){dgamma(x,shape=alpha,rate=beta)*Vectorize(powerfunction)(x)} # the relevant function!
points_to_evaluate_function_at<-seq(from=13,to=30,length.out = N)
sum(h_times_f(points_to_evaluate_function_at[-N])*
      diff(points_to_evaluate_function_at))



integrate(h_times_f,lower=13,upper=30) # the discretization of the integral value!

curve(h_times_f,
      from=0, to=40, main = "Product of gamma distribution and power function 
      in black, rescaled uniform distribution Unif(13,30) in red")
lines(seq(0,40,by=1),0.5*dunif(seq(0,40,by=1),min=13,max=30),col="red")
### F, using importance sampling..
# In this case we might want to sample from a distribution which has large values about where our function has large values.


importance_function<-function(x){
  1/(30-13)
}
h_times_f_over_g<-function(x){h_times_f(x)/Vectorize(importance_function)(x)}

sims_of_important_nu<-runif(N,min=13,max=30)

h_bar_important<-mean(h_times_f_over_g(sims_of_important_nu))
var_importance<-var(h_times_f_over_g(sims_of_important_nu))/(N^2)*(N-1) 
confint_importance<-c(h_bar_important-1.96*sqrt(var_importance),
                      h_bar_important+1.96*sqrt(var_importance))
h_bar_important
var_importance
confint_importance


##########################################################################
# Part 2 -----------------------------------------------------------------
##########################################################################


days<-c(162,267,271,185,111,61,27,8,3,1)
deaths<-seq(0,9)
data<-rbind(deaths,days)
num_days<-sum(days)
mean_death<-sum(days*deaths)/num_days
Likelihoodfun<-function(days,deaths,lambda){
  num_days<-sum(days)
  mean_death<-sum(days*deaths)/num_days
  
  denominator<-prod(factorial(deaths)^days) # really cumbersome computation, generally does not work...
  
  Lh=lambda^(num_days*mean_death)*exp(-num_days*lambda)/denominator
  return(Lh)
}

MLE_lambda<- sum(days*deaths)/num_days # the mean

plot(deaths,days,type="b",cex=2, main = " Observed count in black and \n predicted count using poisson model in red",
     ylab="Days/Count",
     xlab="Deaths on day")
points(num_days*dpois(0:9,lambda=MLE_lambda),col=2,cex=2,pch=2)


#Part 2g

# To implement a gibbs sampler, we need to define a few things first
# we need initial values of the distribution of all lambda parameters, p, Z0-Z9 and 
# We also have to store the sequence of new values somehow, preferably in some relevant vector.
Number_of_samples<- 100000
Z_matrix<-matrix(0,nrow=length(deaths),ncol=Number_of_samples)
lambda_matrix<-matrix(0,nrow=2,ncol=Number_of_samples) # includes both lambdas
p_sequence<-vector(mode="double",length=Number_of_samples)

# initialiation step, somewhat arbitrary
Z_matrix[,1]=floor(days*0.7) 
lambda_matrix[,1]=4
p_sequence[1]=0.7

for (j in 2:Number_of_samples){
  
  # simulate ate from the posterior for lambda1 first, then lambda2, then p, then z
  
  lambda_matrix[1,j]<-rgamma(1,
                             shape=(1+sum(Z_matrix[,j-1]*(0:9))),
                             rate=1+sum(Z_matrix[,j-1]))
  lambda_matrix[2,j]<-rgamma(1,
                             shape=1+sum((days-Z_matrix[,j-1])*(0:9)),
                             rate=(1+sum(days-Z_matrix[,j-1]))
  )
  p_sequence[j]=rbeta(1,shape1=1+sum(Z_matrix[,j-1]),shape2=1+sum(days-Z_matrix[,j-1]))
  
  # For the 10 Z_i variables
  
  for (i in 1:10){
    a= p_sequence[j]*(lambda_matrix[1,j])^(i-1)*exp(-lambda_matrix[1,j])
    b=  (1-p_sequence[j])*(lambda_matrix[2,j])^(i-1)*exp(-lambda_matrix[2,j])
    Z_matrix[i,j]<-rbinom(1,days[i],a/(a+b))  
  }
  # That concludes all simulations for gibbs sampling!
}


mean(p_sequence)
mean(lambda_matrix[1,])
mean(lambda_matrix[2,])
apply(Z_matrix,1,mean)
hist(Z_matrix[1,])
hist(lambda_matrix[2,])


# The proposal function

par(mfrow=c(2,2))
hist(p_sequence, main= "Histogram of simulated p:s")
hist(lambda_matrix[1,],main="Histogram of simulated lambda1:s")
hist(lambda_matrix[2,],main="Histogram of simulated lambda2:s")
hist(Z_matrix[2,],main= "Histogram of simulated Z2:s")
# A Burn in has not been used, however It might be useful!

# To compute expected value, # note E(h(x| the rest))=integral(h(x|the rest)f(the rest)dx) for

# new distribution for deaths is sum of 2 poisson distributions
#dev.off()
dist_improved<-num_days*((mean(p_sequence)*dpois(0:9,lambda=mean(lambda_matrix[1,]))+
                            (1-mean(p_sequence))*dpois(0:9,lambda=mean(lambda_matrix[2,]))))
plot(deaths,days,type="b",cex=2, main = " Observed count in black, predicted \n  count using poisson model in red and 
     predicted count using extended model in blue",
     ylab="Days/Count",
     xlab="Deaths on day")
points(num_days*dpois(0:9,lambda=MLE_lambda),col=2,cex=2,pch=2)
points(deaths,dist_improved, col= 4, cex=3,pch=18)


#relevant distributions are the gamma, poisson, 



