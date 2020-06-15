% Simulate and show an animation of a specific system for a specific realization of the randomness.
% The animation show the behavior of the agents estimates with the given algorithm.

% System configuration
% can be MODIFIED
n = 10;     % number of agents
rho = 10;   % rate ratio
ass_arrcom = false; % assumption of arrival with communication

% supposed always identical
mu = 0;     % distribution mean
sigma = 1;  % standart deviation

%fixed global com rate (n*lc = 10)
lci = (10/n); lri = lci/rho; % individual rates

% Simulation parameters
type = 'event'; % type of simulation -- can be MODIFIED
nevents = 500;  % number of events (only for event-based simulation) -- can be MODIFIED

tf =100;        % final time (only for time-based simulation) --  can be MODIFIED
ts = 0.1;       % time step (only for time-based simulation) --  can be MODIFIED
t = 0:ts:tf;    % time vector
[~,nt] = size(t);

anim = 2;       % 2 for dynamic and 1 for static animation --  can be MODIFIED

% Algorithm(s)
%  can be MODIFIED
algo = 'sympushsum';
param = [false,false]; 
seed = 'shuffle'; % random seed

[mu_est,MSE,Bias,Var,MSE_ext] = OMAS_simu(n,mu,sigma,lri,lci,nevents,tf,ts,anim,type,algo,param,ass_arrcom,seed);

