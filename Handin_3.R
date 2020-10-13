# Handin 3, based on handin 2

# Old model!
days<-c(162,267,271,185,111,61,27,8,3,1)
deaths<-seq(0,9)
data<-rbind(deaths,days)
num_days=sum(days)
MLE_lambda<- sum(days*deaths)/num_days # the mean

plot(deaths,days,type="b",cex=2, main = " Observed count in black and \n predicted count using poisson model in red",
     ylab="Days/Count",
     xlab="Deaths on day")
points(num_days*dpois(0:9,lambda=MLE_lambda),col=2,cex=2,pch=2)

# Defining initial guesses for MLE:s, based on means from previous gibbs sampling


# Function to calculate value of expected value E_z[loglik]

#Expectation step
Expectation_over_Zis<- function(days,p_old,lambda_old){
  i=seq(1:10)
  a= p_old*(lambda_old[[1]]^(i-1)*exp(-lambda_old[[1]]))
  b=  (1-p_old)*(lambda_old[[2]])^(i-1)*exp(-lambda_old[[2]])
  EZi=days*a/(a+b)
  
  return(EZi)
}

p_old=0.1
lambda_old=c(l1=6,l2=1)
lambda_new=lambda_old
counter=1
while(counter < 200){
  EZis<-Expectation_over_Zis(days,p_old,lambda_old)
  #Maximization step
  p_new<-(1+(sum(EZis)/sum(days-EZis)))^-1
  lambda_new[[1]]=sum(0:9*EZis)/(1+sum(EZis))
  lambda_new[[2]]=sum(0:9*(days-EZis))/(1+sum(days-EZis))
  # Stopping criteria
  if (abs(p_new-p_old)<0.01 & 
      abs(lambda_new[[1]]-lambda_old[[1]]) <0.01 &
      abs(lambda_new[[2]]-lambda_old[[2]]) <0.01){
      break
      } 

  p_old<-p_new
  lambda_old<-lambda_new
  counter=counter+1 # To check steps until convergence
}

# Works GREAT!

 
# Computing probabilities from improved model with MAP parameters.
dist_improved<-num_days*(p_new*dpois(0:9,lambda=lambda_new[[1]])+
                            (1-p_new)*dpois(0:9,lambda=lambda_new[[2]]))
plot(deaths,days,type="b",cex=2, main = " Observed count in black, predicted \n  count using poisson model in red and 
     predicted count using extended model in blue",
     ylab="Days/Count",
     xlab="Deaths on day")
points(num_days*dpois(0:9,lambda=MLE_lambda),col=2,cex=2,pch=2)
points(deaths,dist_improved, col= 4, cex=3,pch=18)

sum()