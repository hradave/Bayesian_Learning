### 2a How can you also base your conclusion about NSmallChild on
### the fact that the interval does not cover zero?
# Because the interval does not cover zero, it means that NSmallChild is relevant to determine wether the woman works.
# If it covered zero it'd mean that has neither potivie nor negative impact.
  



### 2b This is not really the predictive distribution but very much related. It is the predictive
### probability of y=1 for different posterior draws of beta. The predictive distribution is a distribution
### over only two values, 0 and 1, so for each probability you should also simulate a binary response variable.
  
# in the log_reg, generate binary variable and plot over that! 
  
  
### 2c This is not correct. You should sample from the binomial distribution once for every draw of your
### "pred" variable and not based on the mean of pred as you do. You can plot a discrete distribution in a 
### nicer way using barplot(table(…))
