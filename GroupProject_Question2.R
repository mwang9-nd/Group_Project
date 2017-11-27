library(deSolve)
library(ggplot2)

#Function for Rosenzweig-MacArthur Predator/Prey Equations
ddRM <- function (t,y,p){
  
  #Unpack State Variables
  H = y[1] #Prey
  P = y[2] #Predator
  
  #Unpack Model Parameters
  b = p[1] #Prey birth rate
  e = p[2] #Conversion efficency of prey to predator
  s = p[3] #Predator death rate
  w = p[4] #Predator attack rate 
  d = p[5] #Density of prey when predator's kill rate is 1/2*max
  alpha = p[6] #Prey self-limiting term
  
  #Defining Rates of Change
  dHdt = b*H*(1 - alpha*H) - w*H/(d + H)*P
  dPdt = e*w*H/(d + H)*P- s*P
  
  return(list(c(dHdt, dPdt)))
}

#Defining Initial Parameters
b = 0.8
e = 0.07
s = 0.2
w = 5
d = 400
alpha = 0.001
params = c(b, e, s, w, d, alpha)

H0 = 500
P0 = 120
Y0 = c(H0, P0)

times = 1:100

modelSim = ode(y=Y0, times = times, func = ddRM, parms = params)

#Subsetting the Results
modelOutput = data.frame(time=modelSim[,1],H=modelSim[,2],P = modelSim[,3])
attach(modelOutput)

#Graphing the Results
ggplot(modelOutput) + 
  geom_line(aes(x = time, y = H), colour = "blue") + 
  geom_line(aes(x = time, y = P), colour = "red")
