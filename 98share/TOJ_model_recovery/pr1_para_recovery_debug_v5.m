% This scripts is based on 6 models specified in v5

% The purpose is to simulate data from 6 models with a large number
% of trials (ExpInfo.nTrials = 1e4, numTrials per SOA) and counts = 1 to
% debug

% The simulated parameters are drawn from distributions from the fitted
% parameters from 'mc_results_v5_n9.mat'

clear all; close all; clc; rng('shuffle');

%% set fixed parameters

% set path
currentDir = pwd;
exptDir = currentDir(1:regexp(pwd,'TOJ_model_recovery')-1);
outDir = [pwd '/pr1_figures'];
addpath(genpath([exptDir 'TOJ_model_recovery']));

if ~exist(outDir, 'dir')
    mkdir(outDir)
end

% simulation parameters
nModel = 6;
CM = zeros(nModel);
counts = 1; % num of simulation

% experimental parameters
ExpInfo.ifi              = 1000/60;
ExpInfo.SOA              = [-0.5000, -0.3000, -0.2500, -0.2000, -0.1500, -0.1000,...
    -0.0500, 0, 0.0500, 0.1000, 0.1500, 0.2000, 0.2500, 0.3000, 0.5000]*1000; % in ms
ExpInfo.lenS             = length(ExpInfo.SOA);
ExpInfo.nTrials          = 1e4; % num of trials per SOA
paraLabel = {{'\mu_{pre}','\mu_{post}','\sigma_{pre}','\sigma_{post}','criterion_{pre}','criterion_{post}','\lambda_{pre}','\lambda_{post}'},...
    {'\mu_{pre}','\mu_{post}','\sigma','criterion_{pre}','criterion_{post}','\lambda_{pre}','\lambda_{post}'},...
    {'\mu_{pre}','\mu_{post}','\sigma_{pre}','\sigma_{post}','criterion','\lambda_{pre}','\lambda_{post}'},...
    {'\mu_{pre}','\mu_{post}','\sigma','criterion','\lambda_{pre}','\lambda_{post}'},...
    {'\mu','\sigma','criterion_{pre}','criterion_{post}','\lambda_{pre}','\lambda_{post}'},...
    {'\mu','\sigma','criterion','\lambda_{pre}','\lambda_{post}'}};

%% summarize fitted parameters
% load model fitted parameters
load('mc_results_v5_n9.mat')

% for each model, extract each parameter from ava_sub, 9 sessions
nPara = [8,7,7,6,6,5];
sub2use = [3,4,6,7,8,9,10];
nSub = numel(sub2use);
nSes = 9;
realP = cell(1, nModel);

% only mu can be negative. Create a matrix for the
% index of mu in each model
mu_idx = [1,1;1,2;2,1;2,2;3,1;3,2;4,1;4,2;5,1;6,1];
sig_idx = [1,3;1,4;2,3;3,3;3,4;4,3;5,2;6,2];
c_idx = [1,5;1,6;2,4;2,5;3,5;4,4;5,3;5,4;6,3];
lambda_idx = [1,7;1,8;2,6;2,7;3,6;3,7;4,5;4,6;5,5;5,6;6,4;6,5];

for m = 1:nModel

    % number of parameters for this model
    m_nPara = nPara(m);

    % initiate parameter matrix (i_para x session x ava_sub)
    m_p = NaN(m_nPara, nSes, numel(sub2use));

    for p = 1:m_nPara
        for s = 1:nSub
            % for each subject, store a matrix (i_para x session)
            m_p(:,:,s) = reshape([estP{sub2use(s)}{:, m}], [], 9);
        end

        % reorganize matrix
        % for each parameter, get a matrix (session x sub)
        realP{m}{p} = squeeze(m_p(p, :, :));

    end
end

% calculate mean and sd for each model, each parameter
[realP_mean, realP_SD] = deal(cell(1, nModel));

for m = 1:nModel

    % number of parameters for this model
    m_nPara = nPara(m);

    for p = 1:m_nPara
        if ismember([m,p], mu_idx,'rows') % if current parameter is mu, which can be negative, summarize them assuming they are from a normal dist
            realP_mean{m}{p} = mean(realP{m}{p},'all');
            realP_SD{m}{p} = std(realP{m}{p},[],'all');

        else % if current parameter is sigma/criterion/lambda, firstly take its log, assume the distribution is normal, take the mean and sd
            realP_mean{m}{p} = mean(log(realP{m}{p}),'all');
            realP_SD{m}{p} = std(log(realP{m}{p}),[],'all');

        end
    end

end


%% simulate and fit for each model

% preallocate
[sampledP, estP, AIC, NLL] = deal(cell(counts, nModel));

parfor count = 1:counts
    disp(count)

    for m = 1:nModel

        %% sample simulated parameters for each model from fitted parameters

        % number of parameters for this model
        m_nPara = nPara(m);

        for p = 1:m_nPara
            % extract mean and sd from fitted parameters
            temp_mean = realP_mean{m}{p};
            temp_sd = realP_SD{m}{p};

            % if current parameter is mu, which can be negative, sample from a normal distribution
            if ismember([m,p], mu_idx,'rows')
                temp_sampledP = randn * temp_sd + temp_mean;

                % if current parameter is sigma/criterion/lambda, sample log(p)
                % from a normal dist, and transfer back p = exp(log(p)), so
                % that p is always positive
            else
                temp_sampledP = exp(randn * temp_sd + temp_mean);

                % if  current parameter is sigma/criterion
                if ismember([m,p], sig_idx,'rows')||ismember([m,p], c_idx,'rows')
                    % if it's larger than 350, resample it
                    while temp_sampledP > 350
                        temp_sampledP = exp(randn * temp_sd + temp_mean);
                    end
                end

                % if  current parameter is sigma/criterion
                if ismember([m,p], lambda_idx,'rows')
                    % if it's larger than 0.06, resample it
                    while temp_sampledP > 0.06
                        temp_sampledP = exp(randn * temp_sd + temp_mean);
                    end
                end

            end
            % store sampled sampled parameters
            sampledP{count, m}(1,p) = temp_sampledP;
        end

        %% simulate fake data

        % input: i_model, model-specific parameter, sim_trial
        sim_r_org = sim_6M_v5(m, sampledP{count, m}, ExpInfo);

        %% fit fake data to all the models

        [AIC{count, m}, NLL{count, m}, estP{count, m}] = fit_6M_v5(sim_r_org, ExpInfo);
    end

end


%% plotting parameter recovery results: sim vs fit

for m = 1:6

    % number of parameters for this model
    m_nPara = nPara(m);

    figure; hold on
    sgtitle(['M' num2str(m)])
    set(gcf, 'Position',[10 10 1500 1200])

    % initiate idx for 'bad' simulation trials
    bad_count = [];
    [sim_p, fit_p] = deal(NaN(m_nPara, counts));

    for p = 1:m_nPara

        for count = 1:counts
            % plot fit parameter against simulated parameter
            sim_p(p, count) = sampledP{count,m}(p);
            fit_p(p, count)= estP{count, m}{m}(p);

        end

        %% plot
        subplot(3,3,p); hold on
        set(gca, 'LineWidth', 1, 'FontSize', 10)
        axis equal
        scatter(sim_p(p,:), fit_p(p,:),40)
        xl = get(gca, 'xlim');
        plot(xl, xl, 'k--')
        title(paraLabel{m}{p})

    end

%     %% plot 'bad' counts in red for all parameters
%     % 'bad' as in simulated values exceeds reasonable range
% 
%     % find 'bad' criterion values > 350
%     [~,idx] = find(sim_p>350);
%     bad_count = [bad_count; idx];
% 
%     % find 'bad' lambda values > 0.06
%     [~,idx] = find(sim_p((m_nPara-1):m_nPara,:)>0.06);
%     bad_count = [bad_count; idx];
%     bad_count = unique(bad_count);
% 
%     for p = 1:m_nPara
%         subplot(3,3,p); hold on
%         scatter(sim_p(p, bad_count), fit_p(p, bad_count),40,'red')
%     end

% save figure
flnm = ['pr_sim_vs_fit_para_M' num2str(m) '_v5'];
saveas(gca, fullfile(outDir, flnm),'epsc')

end
