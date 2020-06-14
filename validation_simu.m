% Validation of the simulation with the theoretical evolution of the Gossip algorithm in open system

% System configuration
sigma = 1;  % standart deviation
mu = 0;     % distribution mean
rho = 19;   % rate ratio (lc/lr)
n = 50;     % number of agents
nevents = 3000; % number of events

% Load Simulation results for the given system
res = load('data/simulation-gossip-02_06_2020-1.mat');
MSE = res.MSE(1:nevents+1); Bias=res.Bias(1:nevents+1); Var=res.Var(1:nevents+1);

% Compute Theoretical results for the same system
p = 1/(rho+1);
A = [1-2*p/n, p/n^2; (1-p)/n, 1-1/n]; b = sigma^2*[p/n^2;p/n];

% descriptor system evolution
EX = zeros(2,nevents+1);
EX(:,1) = [sigma^2/n;sigma^2];
for e=1:nevents
    EX(:,e+1) = A*EX(:,e) + b;
end

% Resulting theoretical expected performance
th_Bias = EX(1,:);
th_Var = EX(2,:) - EX(1,:);
th_MSE = th_Bias + th_Var;

% PLOT simulation vs theory
f1=figure('units','normalized','outerposition',[0.05 0.2 0.9 0.5]);
s1=subplot(1,3,1);
semilogy(0:nevents, MSE,'b','LineWidth',2.5); hold on
semilogy(0:nevents, th_MSE,'--r','LineWidth',3); 
title("Expected MSE")
l=legend('Average of realizations', 'Theoretical expectation');
set(l,'FontSize',11);
xlabel("Events",'FontSize',11);
%ylim([0,1]);
xlim([0,nevents]);

s2=subplot(1,3,3);
semilogy(0:nevents, Bias,'b','LineWidth',2.5); hold on
semilogy(0:nevents, th_Bias,'--r','LineWidth',3); 
title("Expected Squared Bias");
l=legend('Average of realizations', 'Theoretical expectation','Location','northeast');
set(l,'FontSize',11);
xlabel("Events",'FontSize',11);
ylim([0.01,1]);
xlim([0,nevents]);

s3=subplot(1,3,2);
semilogy(0:nevents, Var,'b','LineWidth',2.5); hold on
semilogy(0:nevents, th_Var,'--r','LineWidth',3); 
title("Expected Variance");
l=legend('Average of realizations', 'Theoretical expectation','Location','northeast');
set(l,'FontSize',11);
xlabel("Events",'FontSize',11);
%ylim([0,1]);
xlim([0,nevents]);

% SAVE PDF fig
%set(f1,'PaperSize',[26 8]); %set the paper size to what you want  
%print(f1,'plots/Gossip/th_vs_simu','-dpdf'); % then print it
