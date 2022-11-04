% This script fits Bayesian causal inference model to the pre- and post-
% TOJ data, for each sesion, each participant

% corresponding pathway on hpc: #cd $/scratch/ll3981/project1/CI_model_fit

clear all; close all; clc; rng('Shuffle'); parpool(8);

%% manage paths
currentDir = pwd;
exptDir = currentDir(1:regexp(pwd,'01cau_inf')-1); % project folder
addpath(genpath([exptDir 'data'])); % data folder, necessary to run organize_data
addpath(genpath([exptDir 'bads']));% add bads to path, to be changed later

all_sub = 1:5;
all_ses = 1:9;

for sub = all_sub
    for  ses = all_ses

        %% organize data
        data(sub, ses)= organize_data(sub, ses);

        %% define model parameters

        % set fixed & set-up parameters
        model(sub, ses).expo_num_sim = 1e3; % number of simulation for exposure phase
        model(sub, ses).expo_num_trial = 250; % number of *real* trials in exposure phase
        model(sub, ses).num_runs = 100; %fit the model 100 times, each with a different initialization
        model(sub, ses).num_bin = 100; % num bin to approximate mu_shift
        model(sub, ses).thre_r2 = 0.95; %if R2<0.95, we use KDE

        % define the grid for free parameters
        model(sub, ses).paraID = {'\mu_{pre}','\sigma_{pre}','c_{pre}','\lambda','\sigma_{post}',...
            'c_{post}','p_{common}','\sigma_{soa}','\sigma_{c1}','\sigma_{c2}','alpha'};
        model(sub, ses).numPara = length(model(sub, ses).paraID);

        % hard bounds, the range for LB, UB, larger than soft bounds
        paraH.mu1 = [-0.3, 0.3]; % s
        paraH.sigma1 = [0.01, 0.35]; %s
        paraH.c1 = [0.01, 0.35]; % s
        paraH.lambda = [1e-4, 0.06]; % percentage
        paraH.sigma2 = [0.01, 0.35]; %s
        paraH.c2 = [0.01, 0.35]; % s
        paraH.p_common = [1e-4, 1-1e-4]; % weight
        paraH.sigma_soa  = [0.01, 0.2]; % s
        paraH.sigma_c1 = [1e-4, 0.05]; % s
%         paraH.sigma_c2  = [0.2, 1]; % s
        paraH.sigma_c2  = [0.8, 0.8]; % s
%         paraH.alpha = [1e-3, 0.1]; % percentage
        paraH.alpha = [0.01, 0.01];

        % soft bounds, the range for PLB, PUB
        paraS.mu1 = [-0.25, 0.25]; % s
        paraS.sigma1 = [0.01, 0.3]; %s
        paraS.c1 = [0.01, 0.3]; % s
        paraS.lambda = [1e-4, 0.06]; % percentage
        paraS.sigma2 = [0.01, 0.3]; %s
        paraS.c2 = [0.01, 0.3]; % s
        paraS.p_common = [0.2, 0.5]; % weight
        paraS.sigma_soa  = [0.01, 0.18]; % s
        paraS.sigma_c1 = [1e-4, 0.03]; % s
%         paraS.sigma_c2  = [0.5, 0.8]; % s
        paraS.sigma_c2  = [0.8, 0.8]; % s
%         paraS.alpha = [1e-3, 0.01]; % percentage
        paraS.alpha = [0.01, 0.01];

        % reorganize parameter bounds to feed to bads
        fn = fieldnames(paraH);
        for k = 1:numel(fn)
            model(sub, ses).lb(:,k) = paraH.(fn{k})(1);
            model(sub, ses).ub(:,k) = paraH.(fn{k})(2);
            model(sub, ses).plb(:,k) = paraS.(fn{k})(1);
            model(sub, ses).pub(:,k) = paraS.(fn{k})(2);
        end
        model(sub, ses).paraS = paraS; model(sub, ses).paraH = paraH;

        % set OPTIONS to tell bads that my objective function is noisy
        OPTIONS.UncertaintyHandling = 1;

        %% run bads for optimization

        % define function calculating nll
        funcNLL = @(p) cal_nLL_CI(p(1), p(2), p(3), p(4), p(5), p(6), p(7), ...
            p(8), p(9), p(10), p(11), model(sub, ses), data(sub, ses));

        % get grid initializations
        numREachInit    = 1;
        numSections     = 3;
        model(sub, ses).init      = gridInitializations_11d(model(sub, ses).plb, model(sub, ses).pub, numSections, ...
            model(sub, ses).num_runs, numREachInit);

        %initialize matrices that store negative log likelihood and best-fit paramters
        minNLL          = NaN(1, model(sub, ses).num_runs);
        estimatedP      = NaN(model(sub, ses).num_runs, length(model(sub, ses).lb));

        parfor i = 1:model(sub, ses).num_runs
            disp(i);
            try
                [estimatedP(i,:),minNLL(i)] = bads(funcNLL, model(sub, ses).init(i,:), model(sub, ses).lb,...
                    model(sub, ses).ub, model(sub, ses).plb, model(sub, ses).pub, [], OPTIONS);
                disp(estimatedP(i,:));
                disp(round(minNLL(i),4));
            catch e
                disp('Error!')
            end
        end

        model(sub, ses).estimatedP = estimatedP;
        model(sub, ses).minNLL     = minNLL;

    end

    %% save the data for each participant
    save(['a1ModelFitResults_',datestr(datetime('now'))],'data','model')
    fprintf('sub %.0f sesion %.0f saved', sub, ses)

end