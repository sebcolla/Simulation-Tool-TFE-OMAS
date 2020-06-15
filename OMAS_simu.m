function [yevol,MSE,Bias2,Var,MSE_ext,N,t] = OMAS_simu(n,mu,sigma,lr,lc,nevents,tf,ts,anim,type,algo,param,ass_arrcom,seed)
% Simulation an Open Multi-Agent System (OMAS) for solving the external averaging estimation problem.
% System of n agents, subject to replacements (occuring according to a Poisson clock of global rate n*lr)
% Each agent holds a measurement indentically and independently dranw from a normal distribution of mean mu and variance sigma squared
% Each agents aims to estimate the value of mu.
% Communication: Two randomly, uniformly and independently selected agents communicates in order to update their estimate of mu. 
%                Communications occurs according to a Poisson clock of global rate n*lc
% The simulation is either event-base and run until 'nevents' have occured or time-based and run until time tf is reach (using a time-step ts).
% The simulation return mainly the estimates of the agents and their performances.

% INPUT arguments 
%   - n         :   number of agents in the system
%   - mu        :   true distribution mean of the measurements of the agents
%   - sigma     :   standard deviation of the distribution of the measurements of the agents
%   - lr        :   individual replacement rate (n*lr corresponds to the global replacement rate)
%   - lc        :   individual communication rate (n*lr corresponds to the global communication rate)
%   - nevents   :   number of events to simulate(only for event-base simulation)
%   - tf        :   final time of the simulation (only for time-base simulation)
%   - ts        :   time step (only for time-base simulation)
%   - anim      :   0 no figures produced, 1 bar plot before/after, 2 evolutive barplot
%   - type      :   'time' for a time step simulation and 'event' for an event step simulation
%   - algo      :   specify the algorithm to use in order to solve the problem. 
%                   It can be 'gossip', 'sympushsum','gossip-dt','sps-dt',
%                   'kw-gossip' (for knowledge weighted gossip) or 'best' (for the best ideal algorithm)
%   - ass_arrcom:   boolean to indicate if the assumption of arrival with direct communication holds or not. (Assumption 5.1 in the document)  
%   - param     :   array of parameters to use for the specified algorithm.
%                   Only required for sympushsum or sps-dt. 
%                       param(1): boolean to indicate if the agents should estimate Ntot in order to initialize their weight. 
%                       (if false, the true value of Ntot is used)
%                       param(2): boolean to indicate if a bias correction is performed for the estimate of Ntot 
%                       (should not be true if param(1) is false)
%                   Can be used for algorithmic design purposes (e.g. changing methods parameters easily).
%   - seed      :   random seed. 'shuffle' to chose it randomly.
%
% OUTPUT arguments - TO DO
%   - yevol         :   Matrix with the estimate of each agent after each
%                       event or time step (in function of the parameter
%                       'type'). Each column correspond to one agent and
%                       each line to one event or time step. The first line
%                       correspond to the initial estimate of the agents.
%   - MSE           :   Vector of mean square error for the estimates after each
%                       event or time step (in function of the parameter 'type'). 
%                       First element is the MSE before all the events. 
%                       Vector length is nevents+1 for 'event' based simulation and tf/ts+2 for 'time' based simulation.
%   - Bias2 and Var :   Vectors similar to MSE. We have that MSE = Bias + Var.
%   - MSE_ext       :   Vector of error between the true external mean and the best mean agents can reach, i.e the mean
%                       computed with all the data that have been in the system
%   - N             :   Vector with the value of Ntot for each time step (only for 'time' based simulation)
%   - t             :   Vector recording time corresponding to each value of MSE. Identical dimensions than MSE.

    % Seed of the simulation
    rng(seed);
    
    % INITIALIZATION
    %=================
    % System data
    p = lr/(lr+lc);     % probability of replacement
    %ag_sys = 1:n;       % index of agents in the system

    % Agents data
    z = sigma*randn(n,1) + mu;   % initial measurement of each agent in the system
    ninter = ones(n,1);          % number of interactions of each agents

    % general parameters
    unsucc_com = 0;
    q = 1;
    
    % Algo-dependent initialization
    switch algo
        case 'sympushsum'
            w = ones(n,1);          % weight value for each agent (used in pushsum)
            x = z.*w;               % scale initial intrinsic value
            Nest = n*ones(n,1);
            estim_param = param(1); % agents should estimate the system param (Ntot, SumW,...) or not
            if param(2)
                q = 0.8*n*p+1;
            end
        case 'sps-dt'
            w = ones(n,1);          % weight value for each agent (used in pushsum)
            x = z.*w;               % scale initial intrinsic value
            delta = zeros(n,1);  % delta value of each agent
            deltaw = zeros(n,1); % delta value of each agent on their weight
            beta = 2*p;
            Nest = n*ones(n,1);
            estim_param = param(1); % agents should estimate the system param (Ntot, SumW,...) or not
            if param(2)
                q = 0.8*n*p+1;
            end
        case 'gossip-dt'
            delta = zeros(n,1);  % delta value of each agent
            w = ones(n,1);       % weight value for each agent (used in pushsum)
            x = z.*w;            % scale initial intrinsic value
            beta = p;
        case 'kw-gossip'
            k_set = num2cell(1:n);% knowledge sets of each agent. 'j' is in knowledge set of 'i' if i 'knows' j (maybe not directly) 
            x = z;
            w_set = num2cell(ones(n));% knowledge weights of each agent
            ninter = ones(n,1);
        case 'best'
            k_set = num2cell(1:n);% knowledge sets of each agent. 'j' is in knowledge set of 'i' if i 'knows' j (maybe not directly) 
        otherwise
            x = z;
    end
    
    % Simulation data
    nr = 0;                 % total number of replacements
    
    % Global information
    sumZ = sum(z);          % sum of all initial measurements
    Ntot = n;               % total number of agents that have ever been present in the system
    extav = sumZ/Ntot;    % external average of these agents.
    
    % Initial estimates
    y = evalMean(); 
    
    %figures for ANIMATIONS
    if anim > 0
        %close all;
        figAnim=figure();
        b = bar(1:n,y,'FaceColor','flat'); hold on
        title(sprintf("Open Multi-Agent system - objective: $\\mu = %1.2f $",mu),'Interpreter','latex');
        xlabel("Agents",'Interpreter','latex');
        ylabel("Estimate of $\mu$",'Interpreter','latex');
        ylim([mu-sigma mu+sigma]);
        if anim == 2 % dynamic animation
            l = plot(xlim,[extav extav],'g','LineWidth',1.5);
            an1 = annotation('textbox',[0.2 0.8 0.15 0.08],'String','t = 0 s','FitBoxToText','on','BackgroundColor','white');
            an2 = annotation('textbox',[0.4 0.8 0.15 0.08],'String',sprintf('%d events including %d replacements',0, nr),'FitBoxToText','on','BackgroundColor','white');
            lgd = legend(l,{sprintf('External Average = %.2f',extav)},'Location','southeast');
            set(lgd,'Interpreter','latex');
        end
    end

    if strcmp(type,'time')
    %% ========== TIME STEP SIMULATION ============
        % Time INIT
        t = 0:ts:tf;        % time vector
        [~,nt] = size(t);   % nt stand for the number of time steps
        N = n*ones(nt,1);   % number of agents that have been in the system at each time.
        yevol = zeros(nt,n);% time evolution of the estimates of each agents
        yevol(1,:) = y;

        nev = 0;            % total number of events

        % Error analysis for each time step. 
        MSE = zeros(nt,1);
        Bias2 = zeros(nt,1);
        Var = zeros(nt,1);
        [MSE(1), Bias2(1), Var(1)] = mse(y,mu);
        
        % MSE of the external average
        MSE_ext = zeros(nt,1);
        MSE_ext(1) = (extav-mu)^2;
        
        % Loop on each time-step
        for i = 1:nt-1 
            nevents_t = poissrnd(n*(lr + lc)*ts);   % number of events that occurs at time t (during time step ts)
            for e = 1:nevents_t  % Resolution of each events
                rnb = rand;
                if rnb < p  % event is a REPLACEMENT
                    ag = randi(n);  % random agent to replace
                    repl(ag);       % replacement protocol function
                    N(i+1:end) = Ntot;
                    if(anim == 2)   % dynamic animation
                        b.CData(ag,:) = [1,0,0];
                    end
                    %fprintf("Event %d - Replacement of agent %d : time = %.1fs - Ntot = %d - new x = %.2f - total mean = %.4f\n",nev+e,ag,t(i), Ntot,x(ag),meanAll);
                else        % Event is a COMMUNICATION
                    commu();      % communication protocol function
                end
                if anim == 2 % dynamic animation
                    if ~ishghandle(figAnim)
                        return
                    end
                    y = evalMean(); % estimates updates
                    if mod(i,ceil(nt/(n))) == 0 % update color
                        b.CData(:,:) = b.CData(:,:) + 0.6*(repmat([0,0.447,0.741],n,1) - b.CData(:,:));
                    end
                    set(b,'YDATA',y);
                    set(an1,'String',sprintf("t = %.1f s",t(i)));
                    set(an2,'String',sprintf('%d events including %d replacements',nev, nr));
                    %set(l,'YDATA',[mu(t(i+1)),mu(t(i+1))]);
                    set(l,'YDATA',[extav,extav]);
                    set(lgd,'String',sprintf('External Average = %.2f',extav))%mu(t(i+1))));
                    drawnow;
                end
            end
            nev = nev + nevents_t; %update total nber of events
            y = evalMean(); % estimates updates
            yevol(i+1,:) = y;

            % Error analysis
            [MSE(i+1), Bias2(i+1), Var(i+1)] = mse(y,mu);
            MSE_ext(i+1) = (extav-mu)^2;
        end
        if anim == 1 % static ANIMATION
            bar(1:n,y);
            plot(xlim,[extav, extav],'-g','LineWidth',2);
            %plot(xlim,[mu(tf) mu(tf)],'-g','LineWidth',2);
            annotation('textbox',[0.2 0.8 0.15 0.08],'String',sprintf('%d replacements',nr),'FitBoxToText','on','BackgroundColor','white');
            lgd = legend('0 events - 0.0 s',sprintf('%d events - %.1f s.',nev,t(i+1)),sprintf('External Average = %.2f',extav));
            set(lgd,'Interpreter','latex');
        end
        %fprintf('%d events inwhich %d repl\n',nev,Ntot-n);
    elseif strcmp(type,'event')
    %% ========== EVENT STEP SIMULATION ===========

        % Error analysis for each events. 
        MSE = zeros(nevents+1,1);
        Bias2 = zeros(nevents+1,1);
        Var = zeros(nevents+1,1);
        [MSE(1),Bias2(1),Var(1)] = mse(y,mu);
        
        % MSE of the external average
        MSE_ext = zeros(nevents+1,1);
        MSE_ext(1) = (extav-mu)^2;
        
        % Evolution of the estimates of the agents
        yevol = zeros(nevents+1,n);
        yevol(1,:) = y;
        
        % time between each events
        events = exprnd(1/(n*(lr+lc)),nevents,1); 

        t=zeros(nevents+1,1);   % time vector for the time at which each events occurs
        N = n*ones(nevents+1,1); % agents that have ever been in the system at each event.

        % Loop on all the events
        for e = 1:nevents
            rnb = rand;
            replacement = rnb < p;
            if replacement % event if a replacement 
                    ag = randi(n);  % random agent to replace
                    repl(ag);       % replacement protocol function
                    N(e+1:end) = Ntot;
                if(anim == 2)
                    b.CData(ag,:) = [1,0,0];
                end
                %fprintf("Event %d - Replacement of agent %d : time = %.1fs - Ntot = %d - new x = %.2f - total mean = %.4f - error = %.4f\n",e,ag,t, Ntot,x(ag),meanTot,abs(mu-meanTot));
            else %communication
                commu();          % communication protocol function
            end
            y = evalMean();         % update the estimates
            yevol(e+1,:) = y;
            t(e+1)=t(e)+ events(e); % time of the next event
            
            % Error analysis after this event
            [MSE(e+1),Bias2(e+1),Var(e+1)] = mse(y,mu);
            MSE_ext(e+1) = (extav-mu)^2;
            
            if anim == 2 % dynamic animation
                if ~ishghandle(figAnim)
                    return
                end
                if mod(e,ceil(nevents/(n))) == 0
                    b.CData(:,:) = b.CData(:,:) + 0.6*(repmat([0,0.447,0.741],n,1) - b.CData(:,:));
                end
                set(b,'YDATA',y);
                set(an1,'String',sprintf("t = %.1f s",t(e+1)));
                set(an2,'String',sprintf('%d events including %d replacements',e, nr));
                set(l,'YDATA',[extav,extav]);
                set(lgd,'String',sprintf('External Average = %.2f',extav));
                drawnow;
            end
        end
        if anim == 1 % static animation
            bar(1:n,y);
            plot(xlim,[extav extav],'-g','LineWidth',1.5);
            annotation('textbox',[0.2 0.8 0.15 0.08],'String',sprintf('%d replacements',nr),'FitBoxToText','on','BackgroundColor','white');
            legend('0 events - 0s',sprintf('%d events - %.2f s.',nevents,t(e+1)),sprintf('External Average = %.2f',extav));
        end
    end


 %% =========== AUXIALIARY FUNCTIONS ====================
    function commu(ag1,ag2)
    % Communication of two agents according to the algorithm of the simulation.
    % If no agents indices are provided, it selects independently and uniformly two agents (including twice the same agent)
    % If only one agent is provided, it selects a different second agent.
        
    % Selection of communicating agents
        if nargin == 0
            ags = randi(n,2,1);
        elseif nargin==1
            ags = [ag1,randi(n,1,1)];
            while ags(1) == ags(2)
                ags = [ag1,randi(n,1,1)];
            end
        else
            ags=[ag1,ag2];
        end
        
    % Algo dependent update
        switch algo
            case 'gossip'
                gossipMean = (x(ags(1))+x(ags(2)))/2;
                x(ags(1)) = gossipMean;
                x(ags(2)) = gossipMean;
                %fprintf("%d and %d have %1.3f and %1.3f \n",ags(1),ags(2),x(ags(1)),x(ags(2)));
            case 'sympushsum'
                if estim_param % Estimation of Ntot via max-consensus
                    if all(Nest(ags)==1)
                        unsucc_com = unsucc_com + 1;
                    else
                        if any(Nest(ags)==1) 
                            ag = ags(Nest(ags)==1);
                            Nest(ag) = Nest(ags(Nest(ags)~=1));
                            ntoti = q*(Nest(ag)-n)+n;
                            nri = ntoti - n + 1;
                            w(ag) = compute_weights(n,nri); % new weight
                            x(ag) = x(ag).*w(ag); % scale the new value to get a relevant initial guess
                            Nest(ags) = (Nest(ag)+1)*ones(2,1);
                        else
                           % maximum consensus about Nest
                           Nest(ags) = (max(Nest(ags)))*ones(2,1);                           
                        end
                    end
                end
                % update the x
                inter_x = (x(ags(1)) + x(ags(2)))/2;
                x(ags(1)) = inter_x; x(ags(2)) = inter_x;
                
                % update the w
                inter_w = (w(ags(1)) + w(ags(2)))/2;
                w(ags(1)) = inter_w; w(ags(2)) = inter_w;
            case 'gossip-dt'
                if all(ninter(ags)==0) % Communication between two 'new' agents
                    unsucc_com = unsucc_com+1;    
                elseif any(ninter(ags)==0) % FIRST communication of a 'new' agent with an 'old' agent
                    a1 = ags(ninter(ags)==0); % new agent
                    a2 = ags(ninter(ags)>0);  % old agent
                    
                    % initialization of the pretening value of a1
                    delta(a1) = x(a1)-x(a2);
                    x(a1)=x(a2);
                end
                % update the x
                x(ags) = x(ags) + beta.*delta(ags); % incorporation of the differences
                inter_x = (x(ags(1)) + x(ags(2)))/2;% average of their estimates
                x(ags(1)) = inter_x; x(ags(2)) = inter_x;
                
                % update the delta
                inter_d = 1/2*sum(delta(ags).*(1-beta)); % average of their difference
                delta(ags(1)) = inter_d; delta(ags(2)) = inter_d;
                
            case 'sps-dt'
                if all(ninter(ags)==0) % Communication between two 'new' agents
                    unsucc_com = unsucc_com+1;    
                elseif any(ninter(ags)==0) % FIRST communication of a 'new' agent with an 'old' agent
                    a1 = ags(ninter(ags)==0); % new agent
                    a2 = ags(ninter(ags)>0);  % old agent
                    
                    % initialization of the pretening value of a1 (for x)
                    delta(a1) = x(a1)-x(a2);
                    x(a1)=x(a2);
                    
                    % Initialization of the weight if estim_param is true
                    if estim_param && any(Nest(ags)==1) 
                        ag = ags(Nest(ags)==1);     % new agent
                        Nest(ag) = Nest(ags(Nest(ags)~=1));
                        ntoti = q*(Nest(ag)-n)+n;   % corrected estimation of Ntot
                        nri = ntoti - n + 1;        % corrected estimatof of nbr replacements
                        w(ag) = compute_weights(n,nri); % initial weight
                        x(ag) = x(ag).*w(ag);       % scale the new value to get a relevant initial guess
                        Nest(ags) = (Nest(ag)+1)*ones(2,1); % update their Ntot estimate
                    end
                    
                    % initialization of the pretening value of a1 (for w)
                    deltaw(a1) = w(a1)-w(a2);
                    w(a1)=w(a2);
                end
                % max consensus for Ntot estimation
                if estim_param
                    Nest(ags) = (max(Nest(ags)))*ones(2,1);                           
                end
                
                % update the x
                x(ags) = x(ags) + beta.*delta(ags); % incorporation of the differences
                inter_x = (x(ags(1)) + x(ags(2)))/2;% average of their estimates
                x(ags(1)) = inter_x; x(ags(2)) = inter_x;
                
                % update the w
                w(ags) = w(ags) + beta.*deltaw(ags);% incorporation of the differences
                inter_w = (w(ags(1)) + w(ags(2)))/2;% average of their weights
                w(ags(1)) = inter_w; w(ags(2)) = inter_w;
                
                % update the delta_x
                inter_d = 1/2*sum(delta(ags).*(1-beta)); % average of their difference
                delta(ags(1)) = inter_d; delta(ags(2)) = inter_d;
                
                % update the delta_w
                inter_dw = 1/2*sum(deltaw(ags))*(1-beta);
                deltaw(ags(1)) = inter_dw; deltaw(ags(2)) = inter_dw;
            case 'kw-gossip'
                % knowledge set and knowledge weight of both communicating agents
                K1 = k_set{ags(1)};  W1 = w_set{ags(1)};
                K2 = k_set{ags(2)};  W2 = w_set{ags(2)};
                L1 = length(K1);     L2 = length(K2);
                u_set = union(K1,K2,'sorted'); % union of knowledge sets
                [~,i1,i2] = intersect(K1,K2);  % intersection of the knowledge sets
                Lu = length(u_set);
                
                % setting the weighting between both agents
                if (sum(W1.^2)+sum(W2.^2) - 2*sum(W1(i1).*W2(i2)))==0 % if denominator is 0
                    alpha = 1/2;
                else
                    % weighting that minimize the variance
                    alpha = min(1,max(0,    (sum(W2.^2) - sum(W1(i1).*W2(i2)))/(sum(W1.^2)+sum(W2.^2) - 2*sum(W1(i1).*W2(i2)))  ));
                end

                % updating the knowledge weights of both agents
                uw_set = zeros(Lu,1);
                k1 = 1; k2 = 1; ku = 1;
                for u = u_set
                    if k1 <= L1 && u==K1(k1) && k2 <= L2 && u==K2(k2)
                        uw_set(ku) = alpha*W1(k1) + (1-alpha)*W2(k2);
                        ku = ku + 1; k1 = k1 + 1; k2 = k2 + 1;
                    elseif k1 <= L1 && u==K1(k1)
                        uw_set(ku) = alpha*W1(k1);
                        ku = ku + 1; k1 = k1 + 1;
                    elseif k2 <= L2 && u==K2(k2)
                        uw_set(ku) = (1-alpha)*W2(k2);
                        ku = ku + 1; k2 = k2 + 1;
                    end
                end
                
                % update the x
                wy = [alpha, 1-alpha];
                gossipMean = wy*x(ags);
                x(ags(1)) = gossipMean;
                x(ags(2)) = gossipMean;
                
                % update the knowledge sets and weights
                k_set(ags(1)) = {u_set};
                k_set(ags(2)) = {u_set};
                w_set(ags(1)) = {uw_set};
                w_set(ags(2)) = {uw_set};
            case 'best'
                % union of the knowledge sets
                u_set = union(k_set{ags(1)},k_set{ags(2)});
                k_set(ags(1)) = {u_set};
                k_set(ags(2)) = {u_set};
        end
        % Update global informations
        if ~all(ninter(ags)==0) || (~strcmp(algo,'gossip-dt') && ~strcmp(algo,'sps-dt'))
            ninter(ags) = ninter(ags)+1;
        end
    end

    function repl(ag)
    % replacement of agent 'ag' in the system
        
        % new agent data
        ag_id = Ntot+1;
        z(ag_id) = sigma*randn(1,1) + mu; % initial measurement of the new agent
        ninter(ag) = 0; 
        
        % update global information and simulation data
        nr = nr + 1; % index of this replacement
        sumZ = sumZ+z(ag_id);
        
        % Algo dependent initialization
       switch algo
            case 'sympushsum'
                if estim_param
                    Nest(ag) = 1;
                    x(ag) = z(ag_id);
                    w(ag) = 1;
                else
                    w(ag) = compute_weights(n,nr);  % new weight
                    x(ag) = z(ag_id)*w(ag);         % scale the new value to get a relevant initial guess
                end
            case 'gossip-dt'
                delta(ag) = 0;
                x(ag) = z(ag_id);
            case 'sps-dt'
                % new weight
                if estim_param
                    Nest(ag) = 1;
                    w(ag) = 1;
                else
                    w(ag) = compute_weights(n,nr);
                end
                delta(ag) = 0;
                x(ag) = z(ag_id);
                deltaw(ag) = 0;
                x(ag) = z(ag_id)*w(ag);
            case 'best'
                k_set(ag)={ag_id};
            case 'kw-gossip'
                k_set(ag)={ag_id};
                w_set(ag)={1};
                x(ag) = z(ag_id);
            otherwise
                x(ag) = z(ag_id);
       end
        % Communication of the arrival agent
        if ass_arrcom
            commu(ag);
        end
        
        % Update of global information
        Ntot = Ntot+1;      % update the total number of agent that have ever been in the system
        extav = sumZ/Ntot;  % update the sum of the external average
    end

    function y = evalMean()
    % Evaluation of the estimate of each agent (depending of the algorithm)
        switch algo
            case {'sympushsum','sps-dt'}
                    y = x./w;
            case 'best'
                y = zeros(n,1);
                for id=1:n 
                    % average of the measurements of all the agents in the knowledge set of id
                    y(id) = sum(z(k_set{id}))/length(k_set{id});
                end
            otherwise % for gossip, gossip-dt and kw-gossip
                y = x;
        end
    end

    function [mse,bias2,var] = mse(y,mu)
    % compute the MSE, squared bias and variance of a given vector of
    % estimate, for a given value of mu.
        yb = 1/length(y)*sum(y);
        y2b = 1/length(y)*sum(y.^2);
        yb2 = yb^2;
        bias2 = (yb-mu)^2;
        var = y2b - yb2;
        if var < 1e-16
            var = 0;
        end
        mse = bias2 + var;
    end

    function new_w = compute_weights(n,nr)
    % initial weight of the incoming agents (for sympushsum or sps-dt)
        new_w = ((n-1)/n)^(nr);
    end
end