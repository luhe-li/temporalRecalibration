% This scripts is based on 6 models specified in v5

% The purpose is to obtain a distribution of fitted parameters.

% simulate data with real number of trials (ExpInfo.nTrials = 20, numTrials
% per SOA) and counts = 1000

% for sim_m = 1:n_models
% 	sample the true parameters from fitted parameters
% 	for c = 1:counts
% 		simulate fake dataset (where noise is added because of few trials)
% 		fit with the same model
% 	end
% end

clear all; close all; clc; %rng(1);

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
counts = 1000; % num of simulation

% experimental parameters
ExpInfo.ifi              = 1000/60;
ExpInfo.SOA              = [-0.5000, -0.3000, -0.2500, -0.2000, -0.1500, -0.1000,...
    -0.0500, 0, 0.0500, 0.1000, 0.1500, 0.2000, 0.2500, 0.3000, 0.5000]*1000; % in ms
ExpInfo.lenS             = length(ExpInfo.SOA);
ExpInfo.nTrials          = 20; % num of trials per SOA
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
[realP_mean, realP_SD, sampledP] = deal(cell(1, nModel));

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


%% start parameter recovery each model

% preallocate
[estP, AIC, NLL, deltaAIC] = deal(cell(counts, nModel));

parfor m = 1:nModel

    %% sample ONCE true parameters for each model from fitted parameters

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
        sampledP{m}(1,p) = temp_sampledP;
    end

    %% simulate and fit for each model

    for count = 1:counts

        disp(count)

        %% simulate fake data

        % input: i_model, model-specific parameter, sim_trial
        sim_r_org = sim_6M_v5(m, sampledP{m}, ExpInfo);

        %% fit fake data to all the models

        [AIC{count, m}, NLL{count, m}, estP{count, m}] = fit_6M_v5(sim_r_org, ExpInfo);

        %% update confusion matrix
        % iBest is the index of model that best fits this specific fake data set
        deltaAIC{count, m} = AIC{count, m} - min(AIC{count, m});
        [M iBEST] = min(deltaAIC{count, m});
        BEST = deltaAIC{count, m} == M;
        BEST = BEST / sum(BEST);
        CM(m,:) = CM(m,:) + BEST;

    end

end

save('pr2_results_2')

%% plot results

clear all; close all; clc; 

load('pr2_results_2.mat')

%% plot the estimated parameter against the true parameters

for m = 1:6 % for each model
    true_para = sampledP{m};
    n_para = numel(true_para);
    fit_para = NaN(size(true_para));
    for c = 1:counts
        fit_para(c,:) = estP{c,m}{m}; %c: trial of fits; first m: generating model; second m: fitted model
    end
    figure; clf; hold on % create a figure for each model
    set(gcf, 'Position', [0,0,200*n_para,200])
    sgtitle(['M' num2str(m)])
    for p = 1:n_para % create subplot for each parameter
        subplot(1,n_para,p)
        histogram(fit_para(:,p),'EdgeAlpha',0,'FaceColor','k','FaceAlpha',0.3);
        xline(true_para(p),'r','LineWidth',2)
        title(paraLabel{m}{p})
    end

    flnm = ['pr2_M' num2str(m)];
    saveas(gca, fullfile(outDir, flnm),'epsc')
end

%% plot parameters against each other to check trade-off

% calculated subplot dimension for each model

for m = 1:6
    true_para = sampledP{m};
    n_para = numel(true_para);
    fit_para = NaN(size(true_para));
    for c = 1:counts
        fit_para(c,:) = estP{c,m}{m}; %c: trial of fits; first m: generating model; second m: fitted model
    end

    figure; clf; hold on % create a figure for each model
    set(gcf, 'Position', get(0, 'Screensize'));
    sgtitle(['M' num2str(m)])

    combi = nchoosek([1:(n_para-2)], 2);
    n_combi = size(combi,1);
    n_col = max(factor(n_combi));
    n_row = n_combi/n_col;
    for cc = 1:size(combi,1)
        subplot(n_row, n_col, cc)
        axis square
        x_para = combi(cc,1);
        y_para = combi(cc,2);

        plot(fit_para(:, x_para), fit_para(:, y_para),'o','Color',[repmat(0.3, 1, 3)])

        xlabel(paraLabel{m}{x_para})
        ylabel(paraLabel{m}{y_para})

        % calculate correlation
        [R P_value] = corrcoef(fit_para);
        i_r = R(x_para, y_para);
        i_p = P_value(x_para, y_para);
        
        ttl = sprintf('corr = %.2f, p = %.2f', i_r, i_p);
        title(ttl)
    end
    
    flnm = ['pr2_corr_M' num2str(m)];
    saveas(gca, fullfile(outDir, flnm),'epsc')
end

%% plot confusion matrix because we actually fit with all models

figure; clf;
FM = round(100*CM/sum(CM(1,:)))/100;
t = imageTextMatrix(FM);
colormap("bone")
set(t(FM'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(CM);
set(t, 'fontsize', 22)
title(['count = ' num2str(counts)]);
set(gca, 'xtick', [1:6], 'ytick', [1:6], 'fontsize', 20, ...
    'xaxislocation', 'top', 'tickdir', 'out')
xlabel('fit model')
ylabel('sim model')

flnm = 'pr2_CM';
saveas(gca, fullfile(outDir, flnm),'epsc')
