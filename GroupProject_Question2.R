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

times = 1:400

modelSim = ode(y=Y0, times = times, func = ddRM, parms = params)

#Subsetting the Results
modelOutput = data.frame(time=modelSim[,1],H=modelSim[,2],P = modelSim[,3])
attach(modelOutput)

#Graphing the Results
ggplot(modelOutput) + 
  geom_line(aes(x = time, y = H), colour = "blue") + 
  geom_line(aes(x = time, y = P), colour = "red") + 
  ggtitle("Initial Parameters") +
  labs(y="Population Size", x = "Time")
ggplot(modelOutput, aes(H,P)) + geom_point() + 
  ggtitle("Predator vs. Prey Population") +
  labs(y="Prey Population Size", x = "Predator Population Size")

###Changing Parameters 

##W - the Predator Attack Rate
w1 = 2*w
w2 = 1/2*w

#w1
params = c(b, e, s, w1, d, alpha)
modelSim = ode(y=Y0, times = times, func = ddRM, parms = params)

#Subsetting the Results
modelOutput = data.frame(time=modelSim[,1],H=modelSim[,2],P = modelSim[,3])
attach(modelOutput)

#Graphing the Results
ggplot(modelOutput) + 
  geom_line(aes(x = time, y = H), colour = "blue") + 
  geom_line(aes(x = time, y = P), colour = "red") + 
  ggtitle("Initial Parameters \nwith w = 2*w") +
  labs(y="Population Size", x = "Time")
ggplot(modelOutput, aes(H,P)) + geom_point() + 
  ggtitle("Predator vs. Prey Population") +
  labs(y="Prey Population Size", x = "Predator Population Size")

#w2
params = c(b, e, s, w2, d, alpha)
modelSim = ode(y=Y0, times = times, func = ddRM, parms = params)

#Subsetting the Results
modelOutput = data.frame(time=modelSim[,1],H=modelSim[,2],P = modelSim[,3])
attach(modelOutput)

#Graphing the Results
ggplot(modelOutput) + 
  geom_line(aes(x = time, y = H), colour = "blue") + 
  geom_line(aes(x = time, y = P), colour = "red") +
  ggtitle("Initial Parameters \nwith w = 1/2*w") +
  labs(y="Population Size", x = "Time")
ggplot(modelOutput, aes(H,P)) + geom_point() + 
  ggtitle("Predator vs. Prey Population") +
  labs(y="Prey Population Size", x = "Predator Population Size")

##d - Prey population when predator kill rate is 1/2*max
d1 = 2*d
d2 = 1/2*d


#d1
params = c(b, e, s, w, d1, alpha)
modelSim = ode(y=Y0, times = times, func = ddRM, parms = params)

#Subsetting the Results
modelOutput = data.frame(time=modelSim[,1],H=modelSim[,2],P = modelSim[,3])
attach(modelOutput)

#Graphing the Results
ggplot(modelOutput) + 
  geom_line(aes(x = time, y = H), colour = "blue") + 
  geom_line(aes(x = time, y = P), colour = "red") + 
  ggtitle("Initial Parameters \nwith d = 2*d") +
  labs(y="Population Size", x = "Time")
ggplot(modelOutput, aes(H,P)) + geom_point() + 
  ggtitle("Predator vs. Prey Population") +
  labs(y="Prey Population Size", x = "Predator Population Size")

#d2
params = c(b, e, s, w, d2, alpha)
modelSim = ode(y=Y0, times = times, func = ddRM, parms = params)

#Subsetting the Results
modelOutput = data.frame(time=modelSim[,1],H=modelSim[,2],P = modelSim[,3])
attach(modelOutput)

#Graphing the Results
ggplot(modelOutput) + 
  geom_line(aes(x = time, y = H), colour = "blue") + 
  geom_line(aes(x = time, y = P), colour = "red") +
  ggtitle("Initial Parameters \nwith d = 1/2*d") +
  labs(y="Population Size", x = "Time")
ggplot(modelOutput, aes(H,P)) + geom_point() + 
  ggtitle("Predator vs. Prey Population") +
  labs(y="Prey Population Size", x = "Predator Population Size")



##alpha - Prey self-limiting factor 
alpha1 = 2*alpha
alpha2 = 1/2*alpha


#alpha1
params = c(b, e, s, w, d, alpha1)
modelSim = ode(y=Y0, times = times, func = ddRM, parms = params)

#Subsetting the Results
modelOutput = data.frame(time=modelSim[,1],H=modelSim[,2],P = modelSim[,3])
attach(modelOutput)

#Graphing the Results
ggplot(modelOutput) + 
  geom_line(aes(x = time, y = H), colour = "blue") + 
  geom_line(aes(x = time, y = P), colour = "red") + 
  ggtitle("Initial Parameters \nwith alpha = 2*alpha") +
  labs(y="Population Size", x = "Time")
ggplot(modelOutput, aes(H,P)) + geom_point() + 
  ggtitle("Predator vs. Prey Population") +
  labs(y="Prey Population Size", x = "Predator Population Size")

#alpha2
params = c(b, e, s, w, d, alpha2)
modelSim = ode(y=Y0, times = times, func = ddRM, parms = params)

#Subsetting the Results
modelOutput = data.frame(time=modelSim[,1],H=modelSim[,2],P = modelSim[,3])
attach(modelOutput)

#Graphing the Results
ggplot(modelOutput) + 
  geom_line(aes(x = time, y = H), colour = "blue") + 
  geom_line(aes(x = time, y = P), colour = "red") +
  ggtitle("Initial Parameters \nwith alpha = 1/2*alpha") +
  labs(y="Population Size", x = "Time") +
ggplot(modelOutput, aes(H,P)) + geom_point() + 
  ggtitle("Predator vs. Prey Population") +
  labs(y="Prey Population Size", x = "Predator Population Size")