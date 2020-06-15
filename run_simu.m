% Script to simulate and plot the performance of different algorithms designed in my Master Thesis.
% These algorithms operate in open multi-agent systems subject to replacement of agents. 
% The purpose of each agent is to estimate the value of interest mu from which each agent holds a measurement z.
% To achieve this goal (estimating mu), the agents aims to compute the 'external' average of the system
%(that average the measurements of all the agents that have ever been in the system)
% The metric used to evaluate the performance is the mean squared error (MSE). 
% Its values are averaged over multiple realizations of the same simulation in order to get an approximation of its expectation.

date = '15_06_2020';

%%%%%%%% INITIALIZATION %%%%%%%%%%
% System configuration
% can be MODIFIED
n = 50; % number of agents
rho = 10; % rate ratio
ass_arrcom = false; % assumption of arrival with communication

% supposed always identical
mu = 0;     % distribution mean
sigma = 1;  % standart deviation

%fixed global com rate (n*lc = 10)
lci = (10/n); lri = lci/rho; % individual rates

% Simulation parameters 
type = 'time'; % event-based simulation is not suited for these plots
tf = 300;   % final time  -  can be MODIFIED
ts = 1;     % time step
t = 0:ts:tf; % time vector
[~,nt] = size(t);
nevents = 500;  % number of events (not used)
anim = 0;       % no animation
nreal = 1000;     % number of realizations -  can be MODIFIED

% Algorithm(s) to use : 'gossip', 'gossip-dt', 'sympushsum', 'sps-dt', 'kw-gossip'
% can be MODIFIED
algo = 'sympushsum';
algo2 = 'sps-dt'; 

gossip = true; % upper comparison
best = true;    % lower comparison

% param for the algo (sps-based):
param = [false,false]; 
param2 = [false,false]; % param(1): estimate Ntot, param(2): bias correction

%%%%%%%%%% LOAD results %%%%%%%%%%
%           [algo, algo2, gossip, best]
%load_res = [false, false, false, false]; % simulate all the results
load_res = [true, true, true, true];      % load all the results


file1 = 'data/simulation-sympushsum-15_06_2020-0';
file2 = 'data/simulation-sps-dt-15_06_2020-0';
file_g = 'data/simulation-gossip-15_06_2020-0'; %10, 9, %6
file_b = 'data/simulation-best-15_06_2020-0'; %27_04_2020-10 10, 9, %6 // 06_05_2020-19

if load_res(1)
    res1 = load(file1);
    MSE = res1.MSE(1:nt); Bias=res1.Bias(1:nt); Var = res1.Var(1:nt);
else
    MSE = zeros(nreal,tf/ts+1);  Bias = zeros(nreal,tf/ts+1); Var = zeros(nreal,tf/ts+1);
    mu_est = zeros(nreal,tf/ts+1,n);
end

if load_res(2)
    res2 = load(file2);
    MSE2 = res2.MSE(1:nt); Bias2=res2.Bias(1:nt); Var2 = res2.Var(1:nt);
else
    MSE2 = zeros(nreal,tf/ts+1);  Bias2 = zeros(nreal,tf/ts+1); Var2 = zeros(nreal,tf/ts+1);
    mu_est2 = zeros(nreal,tf/ts+1,n);
end

% Gossip
if load_res(3)
    res = load(file_g);
    MSEg = res.MSE(1:nt); Biasg=res.Bias(1:nt); Varg = res.Var(1:nt);
else
    MSEg = zeros(nreal,tf/ts+1);  Biasg = zeros(nreal,tf/ts+1); Varg = zeros(nreal,tf/ts+1);
    mu_estg = zeros(nreal,tf/ts+1,n);
end

% Best
if load_res(4)
    res = load(file_b);
    MSEbest = res.MSE(1:nt); Biasbest=res.Bias(1:nt); Varbest = res.Var(1:nt);
else
    MSEbest = zeros(nreal,tf/ts+1); Biasbest = zeros(nreal,tf/ts+1); Varbest = zeros(nreal,tf/ts+1);
    mu_estbest = zeros(nreal,tf/ts+1,n);
end
%%
%%%%%%%%% RUN simulations %%%%%%%%%%%
fprintf("n = %d, rho = %1.2f \n",n,rho);
fprintf("running...\t");
MSE_ext = zeros(nreal,tf/ts+1); % MSE external average
for i=1:nreal
   if mod(i,50) == 0, fprintf("%d\t",i); end 
   seed = i;
   if ~strcmp(algo,'false') && ~load_res(1)
       [mu_est(i,:,:),MSE(i,:),Bias(i,:),Var(i,:),MSE_ext] = OMAS_simu(n,mu,sigma,lri,lci,nevents,tf,ts,anim,type,algo,param,ass_arrcom,seed);
   end
   if ~strcmp(algo2,'false') && ~load_res(2)
       [mu_est2(i,:,:),MSE2(i,:),Bias2(i,:),Var2(i,:)] = OMAS_simu(n,mu,sigma,lri,lci,nevents,tf,ts,anim,type,algo2,param2,ass_arrcom,seed);
   end
   if gossip && ~load_res(3)
       [mu_estg(i,:,:),MSEg(i,:),Biasg(i,:),Varg(i,:)] = OMAS_simu(n,mu,sigma,lri,lci,nevents,tf,ts,anim,type,'gossip',param,ass_arrcom,seed);
   end
   if best && ~load_res(4)
       [mu_estbest(i,:,:),MSEbest(i,:),Biasbest(i,:),Varbest(i,:),] = OMAS_simu(n,mu,sigma,lri,lci,nevents,tf,ts,anim,type,'best',param,ass_arrcom,seed);
   end
end
fprintf('\n'); 

%%%%% SAVE results %%%%%%%%
if ~strcmp(algo,'false') && ~load_res(1)
    file = save_in_file(date,algo,ass_arrcom,param,n,rho,lri,lci,mu,sigma,type,tf,ts,nevents,nreal,mu_est,MSE,Bias,Var,MSE_ext);
end
if ~strcmp(algo2,'false') && ~load_res(2)
    file = save_in_file(date,algo2,ass_arrcom,param2,n,rho,lri,lci,mu,sigma,type,tf,ts,nevents,nreal,mu_est2,MSE2,Bias2,Var2,MSE_ext);
end
if ~strcmp(gossip,'false') && ~load_res(3)
    file = save_in_file(date,'gossip',ass_arrcom,param,n,rho,lri,lci,mu,sigma,type,tf,ts,nevents,nreal,mu_estg,MSEg,Biasg,Varg,MSE_ext);
end
if ~strcmp(best,'false') && ~load_res(4)
    file = save_in_file(date,'best',ass_arrcom,param,n,rho,lri,lci,mu,sigma,type,tf,ts,nevents,nreal,mu_estbest,MSEbest,Biasbest,Varbest,MSE_ext);
end
    
    
%%%%%%%%% FIGURES - plot results %%%%%%%%
f1=figure('Position', [50 200 600 400]);    
s1 = subplot(1,1,1);
if ~strcmp(algo,'false')
    semilogy(s1,t, mean(MSE,1),'r','LineWidth', 3,'DisplayName',sprintf('%s',algo)); hold on
end
if ~strcmp(algo2,'false')
    semilogy(s1,t, mean(MSE2,1),'b','LineWidth', 3,'DisplayName',sprintf('%s ',algo2)); hold on
end
if gossip
    semilogy(s1,t, mean(MSEg,1),'k','LineWidth', 3,'DisplayName','Classical averaging gossip'); hold on
end
if best
    semilogy(s1,t, mean(MSEbest,1),'Color',[0 0.7 0],'LineWidth', 3,'DisplayName','Empirical lower bound');
end
grid on;

title(s1,sprintf("Expectation (%d realizations) of the mean squared error (MSE)\n for estimation of $\\mu$ by an open multi-agent system.\n $n=%d$, $\\lambda_c =$ %1.2f and $\\lambda_r =$%1.2f",nreal,n,lci,lri),'Interpreter','Latex');
leg1 = legend(s1);
set(leg1,'Location','northeast','Interpreter','latex','FontSize',13);
xlabel(s1,"Time (s.)",'Interpreter','latex','FontSize',14);
ylabel(s1,"E[MSE]",'Interpreter','latex','FontSize',14);
%ylim(s1,[0.01,1]);
    
f2=figure('Position', [650 200 600 400]);
s2 = subplot(1,1,1); 
if gossip
    semilogy(s2,t, mean(Biasg,1),'k','LineWidth', 3,'DisplayName',sprintf('Gossip : Bias$^2$')); hold on;
    semilogy(s2,t, mean(Varg,1),':k','LineWidth', 2.5,'DisplayName',sprintf('Gossip : Variance'));
end
semilogy(s2,t, mean(Bias,1),'r','LineWidth', 3,'DisplayName',sprintf('theorical %s : Bias$^2$',algo)); hold on;
semilogy(s2,t, mean(Var,1),':r','LineWidth', 2.5,'DisplayName',sprintf('theorical %s : Variance',algo));
if ~strcmp(algo2,'false')
    semilogy(s2,t, mean(Bias2,1),'b','LineWidth', 3,'DisplayName',sprintf('implementable %s : Bias$^2$ ',algo2));
    semilogy(s2,t, mean(Var2,1),':b','LineWidth', 2.5,'DisplayName',sprintf('implementable %s : Variance',algo2));
end
if best
    semilogy(s2,t, mean(Biasbest,1),'Color',[0 0.7 0],'Linewidth', 3,'DisplayName','Empirical lower bound : Bias$^2$');
    semilogy(s2,t, mean(Varbest,1),':','Color',[0 0.7 0],'Linewidth',2.5,'DisplayName','Empirical lower bound : Variance');
end
title(s2,sprintf("Expectation (%d iterations) of the variance and the squared bias\n for estimation of $\\mu$  by an open multi-agent system.\n $n = %d$, $\\lambda_c =$ %1.2f and $\\lambda_r =$%1.2f",nreal,n,lci,lri),'Interpreter','Latex');
leg2 = legend(s2,'E[Bias$^2$]','E[Variance]');
set(leg2,'Location','northeast','Interpreter','latex','FontSize',13);
xlabel(s2,"Time (s.)",'Interpreter','latex','FontSize',14);
grid on;
%%
% SAVE PDF
% set(f1,'PaperSize',[15 11]); %set the paper size to what you want  
% print(f1,'plots/perf_mse_evol_example','-dpdf'); % then print it
% set(f2,'PaperSize',[15 11]); %set the paper size to what you want  
% print(f2,'plots/perf_bv_evol_example','-dpdf'); % then print it

