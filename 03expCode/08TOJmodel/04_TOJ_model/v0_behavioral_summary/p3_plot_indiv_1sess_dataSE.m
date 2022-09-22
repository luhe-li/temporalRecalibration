%% a05_p2_plot_behavioral_results

% part1: plot psychometric curve fitting for one subject, one session, with
% bootstrapping

% part2: plot delta_mu, errorbar = 1SE (68%CI) across subjects for each
% session with subj5 excluded

% part3: plot individual PSS/slope, errorbar = 1SE (68%CI) across sessions
% for each subj with subj5 excluded


clear all; close all; clc; rng(1);

% set data path
currentDir            = pwd;
exptDir               = currentDir(1:regexp(pwd,'03expCode')-1);
outDir                = [exptDir '02figures/fig2_behavioral'];
addpath(genpath([exptDir '03ExpCode/05functions']));
addpath(genpath([exptDir '03ExpCode/06helperFunctions']));
addpath(genpath([exptDir '03ExpCode/01pretest/data']));
addpath(genpath([exptDir '03ExpCode/04posttest/data']));

% best fitting parameters are orgnized by subject x session
load('best_para_2.mat')

% define subject and session to use
isub                  = 8;
isess                 = 9;
bestP                 = best_para{isub, isess};

%% obtain raw data

%%%%% pre-test
load(['pretest_sub' num2str(isub) '_session' num2str(isess) '.mat'])
pre_s_unique          = ExpInfo.SOA; % unique SOA levels, in s
pre_ms_unique         = pre_s_unique * 1e3; % unique SOA levels, in ms
pre_numTrials         = ExpInfo.nTrials; % num of trials per SOA
% inititate
pre_r_org             = NaN(length(pre_s_unique), pre_numTrials);
pre_respCount         = NaN(3, length(pre_s_unique));
for i                 = 1:length(pre_s_unique)
    iSOA                  = pre_s_unique(i);
    iResp                 = Response.order(ExpInfo.trialSOA == iSOA);
    pre_r_org(i,:)        = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j                 = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        pre_respCount(j,i)    = sum(iResp == j);
    end
end
pre_pResp             = pre_respCount/pre_numTrials;

%%%%% post-test
% load data and define key parameters
load(['posttest_sub' num2str(isub) '_session' num2str(isess) '.mat'])
post_s_unique         = ExpInfo.SOA; % unique SOA levels, in ms
post_ms_unique        = post_s_unique * 1e3; % unique SOA levels, in s
post_numTrials        = ExpInfo.nTrials; % num of trials per SOA
% inititate
post_r_org            = NaN(length(post_s_unique), post_numTrials);
post_respCount        = NaN(3, length(post_s_unique));
for i                 = 1:length(post_s_unique)
    iSOA                  = post_s_unique(i);
    iResp                 = Response.order(ExpInfo.trialSOA == iSOA);
    post_r_org(i,:)       = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j                 = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        post_respCount(j,i)   = sum(iResp == j);
    end
end
post_pResp            = post_respCount/post_numTrials;

%% fit psychometric curves

%make a finer grid for the timing difference between the auditory and the
%visual stimulus
SOA_finer             = linspace(pre_ms_unique(1), pre_ms_unique(end), 1000);

P_Afirst              = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst              = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous        = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

pre_fit(1,:) = P_Vfirst(SOA_finer, bestP(1), bestP(3), bestP(5), bestP(7));
pre_fit(2,:) = P_simultaneous(SOA_finer, bestP(1), bestP(3), bestP(5), bestP(7));
pre_fit(3,:) = P_Afirst(SOA_finer, bestP(1), bestP(3), bestP(5), bestP(7));
post_fit(1,:) = P_Vfirst(SOA_finer, bestP(2), bestP(4), bestP(6), bestP(8));
post_fit(2,:) = P_simultaneous(SOA_finer, bestP(2), bestP(4), bestP(6), bestP(8));
post_fit(3,:) = P_Afirst(SOA_finer, bestP(2), bestP(4), bestP(6), bestP(8));

%% plotting

% set plotting parameters
cmp  = [240,99,146; 175, 213, 128; 88,193,238; ... % pre-test
    250,207,226; 226,236, 175; 171,223,235]./255; % post-test
jitter = 10;

% calculate se of data
num_trials = pre_numTrials; % trial number per soa is the same for pre and post;
calSE = @(pv) sqrt(pv.*(1-pv)./num_trials);
pre_SE = arrayfun(@(i) calSE(pre_pResp(i,:)), 1:3, 'UniformOutput',false);
post_SE = arrayfun(@(i) calSE(post_pResp(i,:)), 1:3, 'UniformOutput',false);


% draw error bars on the data points
figure; hold on
for i = 1:3
    e = errorbar(pre_ms_unique + jitter + i*10, pre_pResp(i,:), pre_SE{i}, ...
        'o','MarkerSize',8,'MarkerFaceColor', cmp(i,:), 'MarkerEdgeColor','none');
    e.Color = cmp(i,:); e.LineWidth = 1; e.CapSize = 0;
    e = errorbar(post_ms_unique - jitter - i*10, post_pResp(i,:), post_SE{i}, ...
        'o','MarkerSize',8,'MarkerFaceColor', cmp(i+3,:), 'MarkerEdgeColor','none');
    e.Color = cmp(i+3,:); e.LineWidth = 1; e.CapSize = 0;
end


% initiate figure
figure; hold on
set(gca, 'LineWidth', 2, 'FontSize', 15)
set(gcf, 'Position',[10 10 500 400])

for i = 1:3

    %% pre
    % plot shaded regions (1SE) around data points 
    [l p] = boundedline(pre_ms_unique + jitter, pre_pResp(i,:), pre_SE{i}, 'o','cmap', cmp(i+3,:),'lineWidth', 0.1,'alpha');
    l.MarkerSize = 6; l.MarkerFaceColor = cmp(i+3,:); l.MarkerEdgeColor = 'none';
    % plot best-fitting PMF
    plot(SOA_finer, pre_fit(i,:),'Color', cmp(i+3,:),'LineWidth',2);

    %% post
    % plot shaded regions (1SE) around data points 
    [l p] = boundedline(post_ms_unique - jitter, post_pResp(i,:), post_SE{i}, 'o','cmap', cmp(i,:),'lineWidth', 0.1,'alpha');
    l.MarkerSize = 6; l.MarkerFaceColor = cmp(i,:); l.MarkerEdgeColor = 'none';
    % plot best-fitting PMF
    plot(SOA_finer, post_fit(i,:),'Color', cmp(i,:),'LineWidth',2);
end


