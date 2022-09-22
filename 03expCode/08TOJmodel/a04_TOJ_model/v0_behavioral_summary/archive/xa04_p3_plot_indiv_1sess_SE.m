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
isub                  = 8;
isess                 = 9;
bestP                 = best_para{isub, isess};
i_estP_btst           = estP_btst{isub, isess};


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

% pre-test

%plug in t_diff_finer and the best-fitting mu, sigma and lapse rate to
%function P
pmf(1,:)       = P_Vfirst(SOA_finer, bestP(1), bestP(3), bestP(5), bestP(7));
pmf(2,:)       = P_simultaneous(SOA_finer, bestP(1), bestP(3), bestP(5), bestP(7));
pmf(3,:)       = P_Afirst(SOA_finer, bestP(1), bestP(3), bestP(5), bestP(7));

% errorbar of PMF
numBtst               = 100;

for i                 = 1:numBtst
    Pre_Vfirst_btst(i,:)  = P_Vfirst(SOA_finer, i_estP_btst(i,1), i_estP_btst(i,3), i_estP_btst(i,5), i_estP_btst(i,7));
    Pre_Afirst_btst(i,:)  = P_Afirst(SOA_finer, i_estP_btst(i,1), i_estP_btst(i,3), i_estP_btst(i,5), i_estP_btst(i,7));
    Pre_simul_btst(i,:)   = P_simultaneous(SOA_finer, i_estP_btst(i,1), i_estP_btst(i,3), i_estP_btst(i,5), i_estP_btst(i,7));
end

bounds{1}          = [max(Pre_Vfirst_btst, [], 1) - pmf(1,:); pmf(1,:) - min(Pre_Vfirst_btst, [], 1)];
bounds{2}      = [max(Pre_simul_btst, [], 1) - pmf(2,:); pmf(2,:) - min(Pre_simul_btst, [], 1)];
bounds{3}         = [max(Pre_Afirst_btst, [], 1) - pmf(3,:); pmf(3,:) - min(Pre_Afirst_btst, [], 1)];

% post-test

%plug in t_diff_finer and the best-fitting mu, sigma and lapse rate to
%function P
pmf(4,:)      = P_Vfirst(SOA_finer, bestP(2), bestP(4), bestP(6), bestP(8));
pmf(5,:)       = P_simultaneous(SOA_finer, bestP(2), bestP(4), bestP(6), bestP(8));
pmf(6,:)       = P_Afirst(SOA_finer, bestP(2), bestP(4), bestP(6), bestP(8));

% errorbar of PMF
for i                 = 1:numBtst
    Post_Vfirst_btst(i,:) = P_Vfirst(SOA_finer, i_estP_btst(i,2), i_estP_btst(i,4), i_estP_btst(i,6), i_estP_btst(i,8));
    Post_Afirst_btst(i,:) = P_Afirst(SOA_finer, i_estP_btst(i,2), i_estP_btst(i,4), i_estP_btst(i,6), i_estP_btst(i,8));
    Post_simul_btst(i,:)  = P_simultaneous(SOA_finer, i_estP_btst(i,2), i_estP_btst(i,4), i_estP_btst(i,6), i_estP_btst(i,8));
end

bounds{4}         = [max(Post_Vfirst_btst, [], 1) - pmf(4,:); pmf(4,:) - min(Post_Vfirst_btst, [], 1)];
bounds{5}   = [max(Post_simul_btst, [], 1) - pmf(5,:); pmf(5,:) - min(Post_simul_btst, [], 1)];
bounds{6}         = [max(Post_Afirst_btst, [], 1) - pmf(6,:); pmf(6,:) - min(Post_Afirst_btst, [], 1)];

%% plotting

% set plotting parameters
% cmp = [250,207,226; 226,236, 175; 171,223,235;...
%     240,99,146; 175, 213, 128; 88,193,238]./255; % post-test
% cmp = [242, 185, 197; 137, 208, 236; 207, 234, 142;...% lighter
%     221, 60, 96; 52, 176, 229;  166, 211, 30;]./255; % darker
cmp =  [229, 158, 168; 226,236, 175; 171,223,235;...
    216, 49, 91; 175, 213, 128; 88,193,238]./255;

jitter                = 10;
alpha_value = [repmat(0.6, 1, 3), repmat(0.4, 1, 3)];
line_pattern = {'--','--','--','-','-','-'};


figure; hold on
set(gca, 'LineWidth', 2, 'FontSize', 15)
set(gcf, 'Position',[10 10 500 400])

% plot pmf
for i = [4:6,1:3]
ls(i) = shadedErrorBar(SOA_finer, pmf(i,:), bounds{i}, ...
    'lineProps',{line_pattern{i}, 'Color', cmp(i,:),'lineWidth', 2},'patchSaturation', alpha_value(i));
% ls(i).mainLine.Color = 
% ls(i).patch.FaceColor = patch_cmp(i,:); ls(i).patch.FaceAlpha = 0.6;
end

% plot raw data
for i = 1:3
    plot(pre_ms_unique + jitter, pre_pResp(i,:), 'o','MarkerSize',7,'MarkerFaceColor', 'none','MarkerEdgeColor',cmp(i+3,:),'LineWidth',1.2);
    plot(post_ms_unique - jitter, post_pResp(i,:), 'o','MarkerSize',7,'MarkerFaceColor', cmp(i+3,:),'MarkerEdgeColor',cmp(i+3,:),'LineWidth',1.2);
end

% plot mu and criterion
h1                    = xline(bestP(1),':','LineWidth',1,'Color', repmat(0.3, 1, 3));
% h2                    = xline(bestP(1) - bestP(5),':','LineWidth',1,'Color', repmat(0.5, 1, 3));
% h3                    = xline(bestP(1) + bestP(5),':','LineWidth',1,'Color', repmat(0.5, 1, 3));

h4                    = xline(bestP(2),'LineWidth',1,'Color', repmat(0.3, 1, 3));
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
set(gca, 'FontSize', 20)

% save figure
flnm = 'fig2_indiv_pmf';
saveas(gca, fullfile(outDir, flnm), 'svg')
