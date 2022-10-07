% parameter recovery: determine trial number

% This script generates several fake datasets of different number of
% trials of a ternary TOJ task for pre and post test, bootstrap for 1000
% trials to obtain the confidence interval of parameters, then repeat the
% whole process for 20 times (20 fake datasets)

clear all; close all; clc; rng(1); parpool(2);

%% Experimental parameters
%the levels of SOA
ExpInfo.ifi        = 1000/60;
ExpInfo.s_unique   = [-500, -300, -200, -100, -3*ExpInfo.ifi, -2*ExpInfo.ifi,...
    -ExpInfo.ifi, 0, ExpInfo.ifi, 2*ExpInfo.ifi, 3*ExpInfo.ifi, 100, 200, 300, 500];
% s_unique = [-500, -300, -200, -100, -1000/60*4, -1000/60*2, 0, 1000/60*2, 1000/60*4, 100, 200, 300, 500];
%the total number of levels
ExpInfo.len_deltaT = length(ExpInfo.s_unique);

%% ground truth
%The measurement distribution has a width of sigma_deltaT
Param.sigma_deltaT = 63;
% lapse rate
Param.lapse        = 0.06;
% The absolute criterion for switching decision from ‘A-coincides-V’ to ‘A-precedes-V’ or ‘A-follows-V’
Param.c            = 147.76;
% set true parameters
Param.muPrePost    = [-30, 0];

%% functions 
%the psychometric function after taking lapses into account
Ptilde             = @(P, lambda) lambda/3 + (1-lambda).*P;
P_Afirst           = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*...
                        normcdf(-c, SOA - mu, sig);
P_Vfirst           = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*...
                        (1 - normcdf(c, SOA-mu, sig));
P_simultaneous     = @(SOA, mu, sig, lambda, c) 1 - (lambda/3 + (1-lambda).*...
                        normcdf(-c, SOA - mu, sig)) - (lambda/3 + (1-lambda).*...
                        (1 - normcdf(c, SOA-mu, sig)));
lb                 = [ -150, -150,  10, 1e-2, 10];
ub                 = [  150,  150, 150, 0.06, 300];
options            = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');

%% variables we can manipulate
SimInfo.nDatasets  = 1e2;
SimInfo.numBtst    = 1e3;
SimInfo.nTrials    = 10:5:30;%the number of trial for each audiovisual pair
lenTestTrials      = length(SimInfo.nTrials);
%3 mins for 1 datasets and 30 trials per level
%3 mins x 100 x 5 = 1500 mins = 25 hours (if we pick 5 cores, it can be
%done in less than 6 hours)

%init
[estP_btst, CI_btst_lb, CI_btst_ub] = deal(cell(1,length(SimInfo.nTrials)));

parfor t = 1:length(SimInfo.nTrials)
    disp(['Tested nTrials: ', num2str(t)]);
    simTrial = SimInfo.nTrials(t);
    %[estP_btst_temp, CI_btst_lb_temp, CI_btst_ub_temp] = deal(cell(1,SimInfo.nDatasets));
    %to use parfor, we need to intialize cells here
    for n = 1:SimInfo.nDatasets
        disp(['Tested #dataset: ', num2str(n)]);
        
        r_org_temp = cell(1,2);
        for s  = 1:2 % for pre and post test
            mu  = Param.muPrePost(s);

            % function for creating a cumulative Gausian, the probability of reporting
            % 'A-precedes-V','A-follows-V','A-coincides-V' is:
            P_A_follows_V = 1 - normcdf(Param.c, ExpInfo.s_unique - mu, Param.sigma_deltaT);
            P_A_precedes_V = normcdf(-Param.c, ExpInfo.s_unique - mu, Param.sigma_deltaT);
            P_A_coincides_V = 1 - P_A_precedes_V - P_A_follows_V;

            % on \lambda of the trials they ignore the stimulus and gues, with equal
            % numbers of gueses for each response
            Ptilde_A_follows_V = Ptilde(P_A_follows_V, Param.lapse);
            Ptilde_A_precedes_V = Ptilde(P_A_precedes_V, Param.lapse);
            Ptilde_A_coincides_V = Ptilde(P_A_coincides_V, Param.lapse);
            % simulate fake data
            [bool_Vfirst, ~, bool_Afirst, ~, bool_simul, sim_prob_simul] ...
                = simTernaryTOJ(ExpInfo.len_deltaT, simTrial, Ptilde_A_follows_V,...
                Ptilde_A_precedes_V, Ptilde_A_coincides_V);
            % nT_Vfirst + nT_Afirst + nT_simul should add up to trial number
            % organize raw simulated response into a single matrix (Vfirst = 1,
            % simultaneous = 2, Afirst = 3)
            r_org_temp{1,s} = bool_Vfirst + bool_simul* 2 + bool_Afirst*3;
        end

%         [estP_btst{1,n}, CI_btst_lb{1,n}, CI_btst_ub{1,n}]  = ...
%             BootstrapPrePosttestTOJ(ExpInfo.s_unique, ...
%             r_org_temp{1,1}, r_org_temp{1,2}, simTrial, SimInfo.numBtst, P_Afirst,...
%             P_Vfirst, P_simultaneous, lb, ub, options);  
        
        % bootstrap fake data of pre and post tests
        [estP_btst_temp{1,n}, ~, CI_btst_lb_temp{1,n}, CI_btst_ub_temp{1,n}] = ...
            BootstrapPrePosttestTOJ(ExpInfo.s_unique, ...
            r_org_temp{1,1}, r_org_temp{1,2}, simTrial, SimInfo.numBtst, P_Afirst,...
            P_Vfirst, P_simultaneous, lb, ub, options);       
    end
    estP_btst{t}  = estP_btst_temp;
    CI_btst_lb{t} = CI_btst_lb_temp;
    CI_btst_ub{t} = CI_btst_ub_temp;
end

SimD.estP_btst  = estP_btst;  %
SimD.CI_btst_lb = CI_btst_lb; %
SimD.CI_btst_ub = CI_btst_ub; %

%% save all the data to one file
Data = {ExpInfo, Param, SimInfo, SimD};
save(['PR_determineNumTrials','_',datestr(datetime('now')), '.mat'],'Data');


