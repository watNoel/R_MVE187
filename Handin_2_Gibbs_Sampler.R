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



