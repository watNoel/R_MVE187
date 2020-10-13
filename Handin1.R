
theta<-seq(1.5,1.55,length.out = 1000)
fractions<-c(0.33,0.57,0.10)
expectations<-c(1.5163,1.5197,1.5203)
standardErrors<-c(0.001,0.001,0.005)
standardErrorOfX<-0.001
X<-seq(1.498,1.54,length.out = 1000) # sequence for X in plots except last one with lilelihood ratio

jointDistributionGenerator<-
  function(theta,fractions,
           expectations,standardErrors) {
    temp<-rep(0,length(theta))
    for (i in 1:length(fractions)){
      temp= temp+dnorm(
        theta,expectations[i],standardErrors[i])*fractions[i]
    }
  
    return(Distribution=temp)
  }
   
plot(theta,jointDistributionGenerator(theta,fractions,
                                      expectations,
                                      standardErrors),type="l",main="Prior distribution for theta", ylab= "Probability density")                                

# 2b ----------------------------------------------------------------------
# prior predictive density.
# use result given in lecture 2.3 slide 5 for prior density given a normal density for theta
# sum the 3 priors to a total prior using formula on 3.1 slide 4, that being the law of total variance and expectations




priorX=rep(0,length(X))

for (index in 1:length(fractions))
  {
    priorX<-priorX+fractions[index]*dnorm(X,
                            expectations[index],
                            sqrt(standardErrorOfX^2+standardErrors[index]^2)
                            )
} # the desired cuantity



plot(X,priorX,type="l")


# Posterior predicitve ----------------------------------------------------

# for this we need the normalized weights, and the posteror for each individual disribution


priorXpdf<-function(valuesOfX,fractions,expectations,standardErrors,standardErrorOfX){
  priorX=rep(0,length(valuesOfX))
  for (index in 1:length(fractions)){
  priorX<-priorX+fractions[index]*dnorm(valuesOfX,
                                        expectations[index],
                                        sqrt(standardErrorOfX^2+standardErrors[index]^2)
  )
  }
  return(priorX)
}


plot(X,priorXpdf(X,fractions,expectations,standardErrors,standardErrorOfX),type="l", main = "Prior predictive distribution",ylab="Probability density")

rm(standardErrorOfX)

newWeightsGenerator<-function(valuesOfX,fractions,expectations,standardErrors,standardErrorOfX){
  # only compatible with single value input I think
  newWeights<-rep(0,length(fractions))
  denominator<-priorXpdf(valuesOfX,fractions,expectations,standardErrors,standardErrorOfX)
  for (i in 1:length(fractions)){
    newWeights[i]<-priorXpdf(valuesOfX,fractions[i],
                             expectations[i],
                             standardErrors[i],
                             standardErrorOfX)/denominator
  }
  return(newWeights)
}

posteriorPdfForThetaGenerator<-function(theta,
                                        singleXValue,
                                        fractions,
                                        expectations,
                                        standardErrors,
                                        standardErrorOfX){
  
  newWeights<-newWeightsGenerator(singleXValue,fractions,expectations,standardErrors,standardErrorOfX)
  
  tau<-1/standardErrors^2 
  tauX<-1/standardErrorOfX^2 # useful for updated parameters for dnorm
  
    
  
  posteriorExpectations<-(tauX*singleXValue+tau*expectations)/(tauX+tau)
  posteriorStandarderrors<-sqrt(1/(tau+tauX))
  
  posteriorPdfForTheta<-jointDistributionGenerator(theta,newWeights,
                                                   posteriorExpectations,
                                                   posteriorStandarderrors) # same as the old one
  return(posteriorPdfForTheta)

  }
  

### Next step

standardErrorOfX=0.001  
singleXValue= 1.52083

plot(theta,
  posteriorPdfForThetaGenerator(theta,
                                   singleXValue,
                                   fractions,
                                   expectations,
                                   standardErrors,
                                   standardErrorOfX),
     
     type="l", main="Posterior distribution for theta", ylab= "Probability density") # IT works!!!

##### NOw how do I get a credible interval? Simulate?


thetaRelativeProbs<-posteriorPdfForThetaGenerator(theta,
                              singleXValue,
                              fractions,
                              expectations,
                              standardErrors,
                              standardErrorOfX)
test_simulation<-sample(theta,100000,replace=TRUE,prob=thetaRelativeProbs)


lowerb<-test_simulation[order(test_simulation)[2500]]
upperb<-test_simulation[order(test_simulation)[97500]]

print(paste("Lower bound =",lowerb,"Upper bound =",upperb))
# these give the approximate lower and upper bounds on theta for a 95% CI
integrate(function(theta){posteriorPdfForThetaGenerator(theta,
                                                        singleXValue,
                                                        fractions,
                                                        expectations,
                                                        standardErrors,
                                                        standardErrorOfX)},lowerb,upperb) #  a sanity check, it seems to work out fine
# approximate bounds 

#### The posterior predicitve:
# we make a linear combination of three posterior predictives, each corresponding to using the three separate components of the prior for 
## ditribution for theta. The We may write the integral as pi(xnew/x)=integral(pi(xnew/theta)*pi(theta/xold))dtheta
## and note that the latter of the two factors may be split into a sum of the three components as shown in lecture 3.1
# each of the three contributions to the posteriorpredicitve may the be calculated as follows:

# not for a normal normal pair that the posterior predictive takes the form normal(theta,sigmax^2+sigma_i^2)
# and the arguments theta and sigmai^2 are the parameters of the posterior distribution of theta!


posteriorPredictiveGenerator<-function(valuesOfX,
                                                singleXValue,
                                                fractions,
                                                expectations,
                                                standardErrors,
                                                standardErrorOfX){
  
  newWeights<-newWeightsGenerator(singleXValue,fractions,expectations,standardErrors,standardErrorOfX)
  
  tau<-1/standardErrors^2 
  tauX<-1/standardErrorOfX^2 # useful for updated parameters for dnorm
  
  posteriorExpectations<-(tauX*singleXValue+tau*expectations)/(tauX+tau)
  posteriorStandarderrors<-sqrt(1/(tau+tauX))
  
  
  posteriorPredictivePdf<-priorXpdf(valuesOfX,fractions=newWeights,expectations=posteriorExpectations,standardErrors=posteriorStandarderrors,standardErrorOfX)
  
  return(posteriorPredictivePdf)

}

plot(X,posteriorPredictiveGenerator(X, singleXValue, fractions, expectations,standardErrors,standardErrorOfX),
     type="l", main="Posterior predictive for X", ylab="Probability density") 
## Seems legit 




# Part e ------------------------------------------------------------------

# Likelihood ratios

XValuesForLikelihoodratio<-seq(1.51,1.53,length.out=1000)

numerator<-posteriorPredictiveGenerator(valuesOfX=XValuesForLikelihoodratio,
                                        singleXValue,
                                        fractions,
                                        expectations,
                                        standardErrors,
                                        standardErrorOfX)
denominator<-priorXpdf(valuesOfX=XValuesForLikelihoodratio,
                       fractions,expectations,
                       standardErrors,standardErrorOfX)
plot(XValuesForLikelihoodratio,numerator/denominator,ylab="Likelihood Ratio",xlab=bquote(x[c]), cex.lab=1.5)
abline(v=singleXValue, col="red",lty=2)
legend("topright",legend=bquote(x[s]),col="red",lty=2,cex=1.5)


#plot(X,posteriorPredictiveGenerator(X, singleXValue, fractions, expectations,standardErrors,standardErrorOfX),
#type="l", main="Posterior predictive for X", ylab="Probability density")
#lines(X,priorXpdf(X,fractions,expectations,standardErrors,standardErrorOfX))
#points(XValuesForLikelihoodratio,300/2.5*numerator/denominator,ylab="Something")
