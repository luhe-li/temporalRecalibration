% this script plots the exemplar data, for one subject, one session

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
load('best_para_1.mat')

% define subject and session to use
isub                  = 7;
isess                 = 7;
bestP                 = best_para{isub, isess};
i_estP_btst           = estP_btst{isub, isess};
i_lb_68CI = lb_68CI{isub, isess};
i_ub_68CI = ub_68CI{isub, isess};
numBtst = 1000;


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

%% fit PMF

%make a finer grid for the timing difference between the auditory and the
%visual stimulus
SOA_finer             = linspace(pre_ms_unique(1), pre_ms_unique(end), 1000);

% define functions
P_Afirst              = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst              = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous        = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));
psychFunctions = {P_Vfirst; P_simultaneous; P_Afirst;};

%plug in t_diff_finer and the best-fitting mu, sigma and lapse rate to
%function P

for i = 1:3
   pmf(i,:) = psychFunctions{i} (SOA_finer, bestP(1), bestP(3), bestP(5), bestP(7)); % pre, best fitting
   pmf(i+3,:) = psychFunctions{i} (SOA_finer, bestP(2), bestP(4), bestP(6), bestP(8)); % post, best fitting

%    % a) calculate errorbar by using 68CI of parameters 
%    % a1) find the lb and ub of parameters;
%    % a2) plug in psychFunction'
%    errorbar_lb(i,:) = psychFunctions{i} (SOA_finer,i_lb_68CI(1), i_lb_68CI(3), i_lb_68CI(5), i_lb_68CI(7)); % pre, lower
%    errorbar_lb(i+3,:) = psychFunctions{i} (SOA_finer,i_lb_68CI(2), i_lb_68CI(4), i_lb_68CI(6), i_lb_68CI(8)); % post, lower
%    errorbar_ub(i,:) = psychFunctions{i} (SOA_finer,i_ub_68CI(1), i_ub_68CI(3), i_ub_68CI(5), i_ub_68CI(7)); % pre, higher
%    errorbar_ub(i+3,:) = psychFunctions{i} (SOA_finer,i_ub_68CI(2), i_ub_68CI(4), i_ub_68CI(6), i_ub_68CI(8)); % post, higher
end


% b) calculate errorbar by using 68CI of the curves at each time point
% b1) plug in all btst params in psychFunction;
% b2) find out 68/95 CI of the boundaries at each time point;
% pmf_btst size: numBtst x numTimepoint
for j = 1: numBtst
    for i = 1:3
     pmf_btst{i}(j,:) = psychFunctions{i} (SOA_finer, i_estP_btst(j,1), i_estP_btst(j,3), i_estP_btst(j,5), i_estP_btst(j,7));  
     pmf_btst{i+3}(j,:) = psychFunctions{i} (SOA_finer, i_estP_btst(j,2), i_estP_btst(j,4), i_estP_btst(j,6), i_estP_btst(j,8));
    end
end

for i = 1:6
    for p = 1:1000 % for each time point
        [CI_lb{i}(1,p), CI_ub{i}(1,p)] = get68CI(pmf_btst{i}(:,p));
    end
end

%% plotting

% set plotting parameters
% cmp = [250,207,226; 226,236, 175; 171,223,235;...
%     240,99,146; 175, 213, 128; 88,193,238]./255; % post-test
% cmp = [242, 185, 197; 137, 208, 236; 207, 234, 142;...% lighter
%     221, 60, 96; 52, 176, 229;  166, 211, 30;]./255; % darker
% cmp =  [229, 158, 168; 226,236, 175; 171,223,235;...
cmp =  [229, 158, 168; 203, 227, 172; 171,223,235;...
    216, 49, 91; 175, 213, 128; 88,193,238]./255;

% decide jitter direction by adaptor soa
if ExpInfo.adaptor <= 0
    jitter = 10;
else
    jitter = -10;
end
alpha_value = [repmat(0.6, 1, 3), repmat(0.6, 1, 3)];
line_pattern = {'--','--','--','-','-','-'};

figure; hold on
set(gca, 'LineWidth', 2, 'FontSize', 20)
set(gcf, 'Position',[10 10 700 400])

% plot pmf
for i = [4:6, 1:3] % so that lighter lines are at the front
%     % a) parameter CI
%     patch([SOA_finer, fliplr(SOA_finer)], [errorbar_lb(i,:), fliplr(errorbar_ub(i,:))],cmp(i,:), 'EdgeColor','none','FaceAlpha', alpha_value(i));
    % b) curve CI
    patch([SOA_finer, fliplr(SOA_finer)], [CI_lb{i}, fliplr(CI_ub{i})],cmp(i,:), 'EdgeColor','none','FaceAlpha', alpha_value(i));
    lines(i) = plot(SOA_finer, pmf(i,:), line_pattern{i},'Color', cmp(i,:), 'LineWidth', 1.5);
end

% plot raw data
for i = 1:3
    dots(i) = plot(pre_ms_unique + jitter, pre_pResp(i,:), 'o','MarkerSize',7,'MarkerFaceColor', 'none','MarkerEdgeColor',cmp(i+3,:),'LineWidth',1.2);
    dots(i+3) = plot(post_ms_unique - jitter, post_pResp(i,:), 'o','MarkerSize',7,'MarkerFaceColor', cmp(i+3,:),'MarkerEdgeColor',cmp(i+3,:),'LineWidth',1.2);
end

% plot mu and criterion
h1    = xline(bestP(1),':','LineWidth',2,'Color', repmat(0.3, 1, 3));
% h2                    = xline(bestP(1) - bestP(5),':','LineWidth',1,'Color', repmat(0.5, 1, 3));
% h3                    = xline(bestP(1) + bestP(5),':','LineWidth',1,'Color', repmat(0.5, 1, 3));

h4   = xline(bestP(2),'LineWidth',2,'Color', repmat(0.3, 1, 3));
% h5                    = xline(bestP(2) - bestP(6), 'LineWidth',1,'Color', repmat(0.5, 1, 3));
% h6                    = xline(bestP(2) + bestP(6),'LineWidth',1,'Color', repmat(0.5, 1, 3));

% look better
tick_soa = [-500, -300:100:300,500];
tick_p = [0:0.25:1];
xticks(tick_soa)
xticklabels(strsplit(num2str(tick_soa)))
yticks(tick_p)
yticklabels(strsplit(num2str(tick_p)))
xlabel('SOA (ms)')
ylabel('probability')

% add legend
legendlabels = {'p_{pre}(vision lead)', 'p_{pre}(simultaneous)', 'p_{pre}(audition lead)',...
    'p_{post}(vision lead)', 'p_{post}(simultaneous)','p_{post}(audition lead)',...
    'PSS_{pre}','PSS_{post}'};
legend([dots, h1, h4],legendlabels,'Location','bestoutside')

% save figure
flnm = 'fig2_indiv_pmf';
saveas(gca, fullfile(outDir, flnm), 'epsc')
