% This function plots the raw data against fitting with specified model and
% best-fitting parameters

function model_comp_plotting_v2(subjID, sess, idxBestModel, bestP)

%% organize data
%%%%% pre-test
load(['pretest_sub' num2str(subjID) '_session' num2str(sess) '.mat'])
pre_s_unique = ExpInfo.SOA; % unique SOA levels, in s
pre_ms_unique = pre_s_unique * 1e3; % unique SOA levels, in ms
pre_numTrials = ExpInfo.nTrials; % num of trials per SOA
% inititate
pre_r_org = NaN(length(pre_s_unique), pre_numTrials);
pre_respCount = NaN(3, length(pre_s_unique));
for i = 1:length(pre_s_unique)
    iSOA = pre_s_unique(i);
    iResp = Response.order(ExpInfo.trialSOA == iSOA);
    pre_r_org(i,:) = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        pre_respCount(j,i) = sum(iResp == j);
    end
end
pre_pResp = pre_respCount/pre_numTrials;
pre_nT_A1st = pre_respCount(3,:);
pre_nT_V1st = pre_respCount(1,:);
pre_nT_simul = pre_respCount(2,:);

%%%%% post-test
% load data and define key parameters
load(['posttest_sub' num2str(subjID) '_session' num2str(sess) '.mat'])
post_s_unique = ExpInfo.SOA; % unique SOA levels, in ms
post_ms_unique = post_s_unique * 1e3; % unique SOA levels, in s
post_numTrials = ExpInfo.nTrials; % num of trials per SOA
% inititate
post_r_org = NaN(length(post_s_unique), post_numTrials);
post_respCount = NaN(3, length(post_s_unique));
for i = 1:length(post_s_unique)
    iSOA = post_s_unique(i);
    iResp = Response.order(ExpInfo.trialSOA == iSOA);
    post_r_org(i,:) = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        post_respCount(j,i) = sum(iResp == j);
    end
end
post_pResp = post_respCount/post_numTrials;
post_nT_A1st = post_respCount(3,:);
post_nT_V1st = post_respCount(1,:);
post_nT_simul = post_respCount(2,:);

%% define and select best model

%make a finer grid for the timing difference between the auditory and the
%visual stimulus
SOA_finer = linspace(pre_ms_unique(1), pre_ms_unique(end), 1000);

% define PMF
P_Afirst = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

M1 = {@(p) P_Afirst(SOA_finer, p(1), p(3), p(5), p(7));...
    @(p) P_Vfirst(SOA_finer, p(1), p(3), p(5), p(7));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(5), p(7));...
    @(p) P_Afirst(SOA_finer, p(2), p(4), p(6), p(8));...
    @(p) P_Vfirst(SOA_finer, p(2), p(4), p(6), p(8));...
    @(p) P_simultaneous(SOA_finer, p(2), p(4), p(6), p(8))};

M2 = {@(p) P_Afirst(SOA_finer, p(1), p(3), p(4), p(6));...
    @(p) P_Vfirst(SOA_finer, p(1), p(3), p(4), p(6));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(4), p(6));...
    @(p) P_Afirst(SOA_finer, p(2), p(3), p(5), p(7));...
    @(p) P_Vfirst(SOA_finer, p(2), p(3), p(5), p(7));...
    @(p) P_simultaneous(SOA_finer, p(2), p(3), p(5), p(7))};

M3 = {@(p) P_Afirst(SOA_finer, p(1), p(3), p(5), p(6));...
    @(p) P_Vfirst(SOA_finer, p(1), p(3), p(5), p(6));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(5), p(6));...
    @(p) P_Afirst(SOA_finer, p(2), p(4), p(5), p(7));...
    @(p) P_Vfirst(SOA_finer, p(2), p(4), p(5), p(7));...
    @(p) P_simultaneous(SOA_finer, p(2), p(4), p(5), p(7))};

M4 = {@(p) P_Afirst(SOA_finer, p(1), p(3), p(4), p(5));...
    @(p) P_Vfirst(SOA_finer, p(1), p(3), p(4), p(5));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(4), p(5));...
    @(p) P_Afirst(SOA_finer, p(2), p(3), p(4), p(6));...
    @(p) P_Vfirst(SOA_finer, p(2), p(3), p(4), p(6));...
    @(p) P_simultaneous(SOA_finer, p(2), p(3), p(4), p(6))};

M5 = {@(p) P_Afirst(SOA_finer, p(1), p(2), p(3), p(4));...
    @(p) P_Vfirst(SOA_finer, p(1), p(2), p(3), p(4));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(3), p(4));...
    @(p) P_Afirst(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_Vfirst(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(3), p(5))};

models = {M1; M2; M3; M4; M5};

%% get djusted R^2
[~, aR2] = model_comp_gof_v2(subjID, sess, idxBestModel, bestP);

%% plug in best-fitting parameters and plot

% plot raw data
cMAP1 = [158,202,225; 161,217,155; 252,174,145]./255;
cMAP2 = [33,113,181; 35,139,69; 203,24,29]./255;
jitter = 5;

% plot raw data
figure; hold on
set(gca, 'LineWidth', 2, 'FontSize', 15)
set(gcf, 'Position',[10 10 900 600])
for i = 1:3
    scatter(pre_ms_unique + jitter, pre_pResp(i,:),'MarkerFaceColor', cMAP1(i,:),...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.7)
    scatter(post_ms_unique - jitter, post_pResp(i,:),'MarkerFaceColor', cMAP2(i,:),...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.7)
end

% plot fitting curves
Pre_Afirst_fit = models{idxBestModel}{1}(bestP);
Pre_Vfirst_fit = models{idxBestModel}{2}(bestP);
Pre_simul_fit = models{idxBestModel}{3}(bestP);
Post_Afirst_fit = models{idxBestModel}{4}(bestP);
Post_Vfirst_fit = models{idxBestModel}{5}(bestP);
Post_simul_fit = models{idxBestModel}{6}(bestP);

plot(SOA_finer, Pre_Afirst_fit, '--','Color', cMAP1(3,:), 'lineWidth', 3)
plot(SOA_finer, Pre_Vfirst_fit, '--','Color', cMAP1(1,:), 'lineWidth', 3)
plot(SOA_finer, Pre_simul_fit, '--','Color', cMAP1(2,:), 'lineWidth', 3)

plot(SOA_finer, Post_Afirst_fit, '-','Color', cMAP2(3,:), 'lineWidth', 3)
plot(SOA_finer, Post_Vfirst_fit, '-','Color', cMAP2(1,:), 'lineWidth', 3)
plot(SOA_finer, Post_simul_fit, '-','Color', cMAP2(2,:), 'lineWidth', 3)

% display best fitting parameters
annotation( 'textbox', 'String', round(bestP,3), ...
        'FontSize', 14, 'Units', 'normalized', 'EdgeColor', 'none', ...
        'Position', [0.9,0.9,0.03,0.03], 'Interpreter','Latex')

% look better
xlabel('SOA (ms)')
ylabel('probability')
title(['sub' num2str(subjID) ' adaptor SOA = ' num2str(ExpInfo.adaptor*1000) ' ms'],['adjusted R^2 = ' num2str(aR2)])
% title(['sub' num2str(subjID) ' adaptor SOA = ' num2str(ExpInfo.adaptor*1000) ' ms'])
end