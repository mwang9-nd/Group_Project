library(deSolve)
library(ggplot2)

### Function for Rosenzweig-MacArthur Predator/Prey Equations ###
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
  dHdt = b*H*(1 - alpha*H) - w*H/(d + H)*P #Rate of Prey population over time
  dPdt = e*w*H/(d + H)*P- s*P #Rate of Predator population over time
  
  return(list(c(dHdt, dPdt)))
}

### Defining Initial Parameters ###
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

### Simulating model with Inital conditions
modelSim = ode(y=Y0, times = times, func = ddRM, parms = params)

#Subsetting the Results
modelOutput = data.frame(time=modelSim[,1],H=modelSim[,2],P = modelSim[,3])

#Graphing the Results
ggplot(modelOutput) + 
  geom_line(aes(x = time, y = H), colour = "blue") + 
  geom_line(aes(x = time, y = P), colour = "red") + 
  ggtitle("Initial Parameters") +
  labs(y="Population Size", x = "Time")
ggplot(modelOutput, aes(H,P)) + geom_point() + 
  ggtitle("Predator vs. Prey Population") +
  labs(y="Prey Population Size", x = "Predator Population Size")

### Changing Parameters ### 

##W - the Predator Attack Rate
w1 = 2*w
w2 = 1/2*w

#w1
params = c(b, e, s, w1, d, alpha)
modelSim = ode(y=Y0, times = times, func = ddRM, parms = params)

#Subsetting the Results
modelOutput = data.frame(time=modelSim[,1],H=modelSim[,2],P = modelSim[,3])

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

#Graphing the Results
ggplot(modelOutput) + 
  geom_line(aes(x = time, y = H), colour = "blue") + 
  geom_line(aes(x = time, y = P), colour = "red") +
  ggtitle("Initial Parameters \nwith alpha = 1/2*alpha") +
  labs(y="Population Size", x = "Time")
ggplot(modelOutput, aes(H,P)) + geom_point() + 
  ggtitle("Predator vs. Prey Population") +
  labs(y="Prey Population Size", x = "Predator Population Size")


### Explaining Parameters ###
  
# ---------------
# Prey Birth Rate
# ---------------

# At low birth rate (0.2-1.6)
# Prey and Predator population will decrease initially before increase in prey then increase in predator to equilibrium
  # Predator death due to starvation
  # Starvation period allows Prey population to increase to self-limiting density
  # Increase in Prey population allows Predator population to increase to capacity
  # Prey population decreases to equilibrium capacity due to Predator Population

# At high birth rate (2.0-3.2)
# Prey birthrate increases to self-limiting density before Predator decreases Prey rate to equilibrium
  # Initial Predator population does not require as much Prey to sustain itself, therefore
  # Prey population reaches self-limiting density
  # Predator population will continue to increase until capacity is reached
  # Capacity for Prey and Predator will increase as birth rate of prey increases

  
# -----------------------------------------
# Conversion efficiency of Prey to Predator
# -----------------------------------------

# Low efficiency conversion (0.0175-0.0525)
  # Rate of Prey will initially decrease because Predator is present in system
  # Rate of Predator decreases because of efficiency to convert energy to more Predator
  # Population of Predator will decrease to 0 as Prey population increases to self-limiting density 
      # Rate to capacity determined by other parameters in the function

# Medium efficiency conversion and mutual equilibrium (0.07-0.0875)
  # Rate of Prey and Predator will oscillate
      # Prey population drops as Predator starts consuming
      # Predator reaches a population size that is not sustained by the amount of Prey
          # Predator population will decrease
      # Prey will bounce back and cycle will begin again
  # Prey and Predator will eventually reach a stable equilibrium

# High efficiency conversion (0.1050-0.28)
  # Rate of exchange from Prey to Predator occurs too quickly
  # System will change from decreasing Prey and increasing Predator until Predator density is reached
  # Predators will decrease from starvation or lack of Prey to support system until Prey numbers bounch back
  # Increase in efficiency will increase maximum capacity of predator
 
  
# -------------------
# Predator Death Rate
# -------------------
  
# Low Death Rate (0.05-0.2)
  # Low death rate leads to depletion of Prey population
  # Decrease in Prey population decreases Predator density limit
  # Oscillation of both population until equilibrium is reached, or will continue to oscillate between capacity and 0
  
# High Death Rate (0.25-0.8)
  # Initial dip in Prey population results from initial Predator population
  # High death rate results in eradication of Predator and increase of Prey Population to self-limiting density
  

# -------------------
# Predator Attack Rate
# -------------------

# Low Attack Rate 
  # Low attack rate means that predators cannot sustain themselves, and their population goes to zero
  # Prey population reach self-limiting equilibrium

# High Attack Rate
  # Predators are able to survive
  # High attack rate implies that prey will be driven to a small population. Predators can't survive if prey population is small, so they decreases
  # As predator population decreases, prey population is able to increase. The cyle repeats
  # Long-term behaivor is stable equilibrium

# ---------------
# Prey Denisty to Maximize Predator Attack
# ---------------
  
# Low d values
  # Predator population can maximize attack rate more easily
  # Maximizing attack rate drives prey population down
  # Small prey population hurts the predators, so predators decrease
  # The system oscillates 
  
# High d values 
  # Predators cannot reach an attack rate that sustains their population
  # Predator population goes to zero; prey go to carrying capacity
  
# -----------------------------------------
# Carrying Capacity/Alpha (alpha = 1/k)
# -----------------------------------------
  
# Low alpha
  # Low alpha implies higher carrying capacity of prey
  # Higher carrying capacity implies higher population; this sustains predators
  # System oscillates
  
# High alpha
  # Carrying Capacity is lower, and predators cannot sustain themseleves when there are less prey
  # Prey go to their carrying capacity 
