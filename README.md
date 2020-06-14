# Simulation-Tool-TFE-OMAS
Simulation tool used during my master thesis at UCLouvain: "Data Fusion in Open Multi-Agent Systems for Decentralized Estimation"
(Advisor: Julien Hendrickx)

 ## Simulation
This repository contains a matlab simulation function of the system, the problem and different algorithms to solve it (see file `OMAS_simu.m`. The algorithm has been validated using theoreitcal results about the gossip averaging algorithm (see file `validation_simu.m`). There also is a framework for plotting the resulting performance of different simulated algorithms (see file `run_simu.m`)

See the content of my master thesis in order to have details about the problem, its performance limitations, the algorithms and their in-depth analysis. A synthetic statement of the problem considered is exposed below.

## Synthetic Problem Statement
**The system** is an open multi-agent system, subject to replacements of agents during its operation. The agents present in the system at time t are described by the set of indices N(t). The size of the system, denoted n = |N(t)|, is constant thanks to the replacement assumption: the departure of an agent from the system is immediately followed by the arrival of a new agent.

**The agents** are supposed to be identical, to have a bounded memory, to be capable of local computations and capable of communicating with each other, without having access to any kind of identifier. Each agent ![i \in N(t)](https://render.githubusercontent.com/render/math?math=i%20%5Cin%20N(t))  holds an initial value ![z_i](https://render.githubusercontent.com/render/math?math=z_i) (also called initial measurement) independently drawn from an identical distribution of constant mean ![\mu](https://render.githubusercontent.com/render/math?math=%5Cmu) and variance ![\sigma^2](https://render.githubusercontent.com/render/math?math=%5Csigma%5E2). The mean ![\mu](https://render.githubusercontent.com/render/math?math=%5Cmu) is also called the value of interest. Each agent ![i \in N(t)](https://render.githubusercontent.com/render/math?math=i%20%5Cin%20N(t)) is supposed to execute the same algorithm in order to achieve a common desired goal for which it has its own solution estimate ![y_i(t)](https://render.githubusercontent.com/render/math?math=y_i(t)).

**The goal** of each agent is to estimate \emph{as accurately as possible} the mean distribution ![\mu](https://render.githubusercontent.com/render/math?math=%5Cmu):

![y_i(t) \to \mu, \qquad \text{for all i} ](https://render.githubusercontent.com/render/math?math=y_i(t)%20%5Cto%20%5Cmu%2C%20%5Cqquad%20%5Ctext%7Bfor%20all%20i%7D%20)

Since the value of interest ![\mu](https://render.githubusercontent.com/render/math?math=%5Cmu): is constant, the best way of estimating it is to compute the external average of the system. The external average of the system is the average of the initial values ![z_k](https://render.githubusercontent.com/render/math?math=z_k) of all agents k that have ever been in the system until time t. The goal of each agent ![i \in N(t)](https://render.githubusercontent.com/render/math?math=i%20%5Cin%20N(t)) can then be summarized as follows:

![y_i(t) \to \mu \approx \frac{1}{|\cup_{s\le t}N(s)|}\sum_{k\in\cup_{s\le t}N(s)}z_k \quad \text{for all i} \qquad \text{(External averaging problem)}](https://render.githubusercontent.com/render/math?math=y_i(t)%20%5Cto%20%5Cmu%20%5Capprox%20%5Cfrac%7B1%7D%7B%7C%5Ccup_%7Bs%5Cle%20t%7DN(s)%7C%7D%5Csum_%7Bk%5Cin%5Ccup_%7Bs%5Cle%20t%7DN(s)%7Dz_k%20%5Cquad%20%5Ctext%7Bfor%20all%20i%7D%20%5Cqquad%20%5Ctext%7B(External%20averaging%20problem)%7D)
   
**The metric** used to quantify the accuracy and the success of an algorithm to achieve this goal is the mean squared error criterion (MSE) computed between the agents estimates ![y_i(t)](https://render.githubusercontent.com/render/math?math=y_i(t)) and the value of interest ![\mu](https://render.githubusercontent.com/render/math?math=%5Cmu). It is defined as:

![\text{MSE}(t, \omega) = \frac{1}{n} \sum_{i\in N(t, \omega)} (y_i(t, \omega)-\mu)^2](https://render.githubusercontent.com/render/math?math=%5Ctext%7BMSE%7D(t%2C%20%5Comega)%20%3D%20%5Cfrac%7B1%7D%7Bn%7D%20%5Csum_%7Bi%5Cin%20N(t%2C%20%5Comega)%7D%20(y_i(t%2C%20%5Comega)-%5Cmu)%5E2)

It depends on the realization ![\omega \in \Omega](https://render.githubusercontent.com/render/math?math=%5Comega%20%5Cin%20%5COmega) of the randomness into the system (the realization of the initial values of the agents and of the sequence of communication and replacement events that occur in the system). The expected mean square error ![E\[\text{MSE}(t, \omega)\]](https://render.githubusercontent.com/render/math?math=E%5B%5Ctext%7BMSE%7D(t%2C%20%5Comega)%5D) is used to evaluate the expected performance of an algorithm.
    
**The communications** occur in the system according to a Poisson clock of global rate ![n\lambda_c](https://render.githubusercontent.com/render/math?math=n%5Clambda_c). Whenever a communication occurs, two randomly uniformly and independently selected agents ![i,j \in N(t)](https://render.githubusercontent.com/render/math?math=i%2Cj%20%5Cin%20N(t)) (possibly twice the same agent) exchange information with each other. In a nutshell, communications are said to be asynchronous, pairwise, symmetric and random. 

 **The replacements** occur in the system according to a Poisson clock of global rate ![n\lambda_r](https://render.githubusercontent.com/render/math?math=n%5Clambda_r). The value of ![\lambda_r](https://render.githubusercontent.com/render/math?math=%5Clambda_r) therefore the individual replacement rate. A replacement corresponds to the departure of one randomly and uniformly chosen agent from the system, instantaneously followed by the arrival of a new agent into the system.
 

