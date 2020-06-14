function file_name = save_in_file(date,algo,ass_arrcom,param,n,rho,lri,lci,mu,sigma,type,tf,ts,nevents,nreal,mu_est,MSE,Bias,Var,MSE_ext)
% Save results of a set of simulations from 'OMAS_simu' function in a .mat document
% Error information (MSE, Bias, Var and MSE_ext) are averaged over all the realizations.
%
% INPUT arguments 
%   - date      :   date for the file_name
%   - algo      :   specify the algorithm used in order to solve the problem. 
%                   It can be 'gossip', 'sympushsum','gossip-dt','sps-dt',
%                   'kw-gossip' (for knowledge weighted gossip) or 'best' (for the best ideal algorithm)
%   - ass_arrcom:   boolean to indicate if the assumption of arrival with direct communication holds or not. (Assumption 5.1 in the document)  
%   - param     :   array of parameters to use for the specified algorithm.
%   - n         :   number of agents in the system
%   - rho       :   rate ratio in the system (lci/lri)
%   - lri       :   individual replacement rate (n*lr corresponds to the global replacement rate)
%   - lci       :   individual communication rate (n*lr corresponds to the global communication rate)
%   - mu        :   true distribution mean of the measurements of the agents
%   - sigma     :   standard deviation of the distribution of the measurements of the agents
%   - type      :   'time' for a time step simulation and 'event' for an event step simulation
%   - tf        :   final time of the simulation (only for time-base simulation)
%   - ts        :   time step (only for time-base simulation)
%   - nevents   :   number of events to simulate(only for event-base simulation)
%   - nreal     :   number of realization csimulated
%   - mu_est    :   estimate evolution of mu for each agent
%   - MSE       :   evolution of the MSE of the system for each run
%   - Bias      :   evolution of the Bias of the system for each run
%   - Var       :   evolution of the Var of the system for each run
%   - MSE_ext   :   evolution of the error of the MSE of the external average for each run
%
%   OUTPUT:
%       - file_name:  file_name where data have been saved

    % File name
    version = 0;
    file_name = sprintf('data/simulation-%s-%s-%d.mat',algo,date,version);
    while isfile(file_name) % File already exists.
        version = version + 1;
        file_name = sprintf('data/simulation-%s-%s-%d.mat',algo,date,version);
    end 
    
    % Mean of the collected data
    MSE = mean(MSE,1); Bias = mean(Bias,1);
    Var = mean(Var,1); MSE_ext = mean(MSE_ext,1);
    mu_est = mean(mu_est,1);
    
    % save parameters and data in a .mat file
    save(file_name,'algo','ass_arrcom','param','n','lri','lci','rho','mu','sigma','type','tf','ts','nevents','nreal','mu_est','MSE','Bias','Var','MSE_ext');
end