%% TOJ fitting by cumulative gaussian model v4

% 2020/04

% This script fits **cumulative Gaussian PMF** for the ternary TOJ data.

% In data-cleaning, it extracts two key outputs:
% r_org :  this matrix has a size of length(s_unique) x numTrials
% respCount: this matrix has the size of 3 (type of response: 1 = V first,
% 2 = simultaneous, 3 = A first) x length(s_unique)

% In fitting, we fitted the PMF use fmincon for 1e3 times and chose the
% best-fitting parameters based on the smallest nLL.

% Data were then bootsrapped for 1e3 times to calculate error bars for
% each data point. No fitting in bootstrapping.

% In plotting, we plotted the raw data, the fitted PMF with the
% best-fitting parameters, and errorbars of data, and saved the figure.

clear all; close all; clc;

subjID                  = 3;
sess                    = 2;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% data-cleaning %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean raw data
addpath('/Volumes/GoogleDrive-107780329657346345479/My Drive/2021A-22/01project/04ExperimentalCode/03ExpCpde/01pretest/data')
load(['pretest_sub' num2str(subjID) '_session' num2str(sess) '.mat'])

% define variables used in MLE
s_unique                = ExpInfo.SOA; % unique SOA levels, in ms
numTrials               = ExpInfo.nTrials; % number of trials per SOA

% organize response
for i                   = 1:length(s_unique)
    iSOA                    = s_unique(i);
    iResp                   = Response.order(ExpInfo.trialSOA == iSOA);
    r_org(i,:)              = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j                   = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        respCount(j,i)          = sum(iResp == j);
    end
end

% convert to probability
pResp                   = respCount/numTrials;

% convert s_unique from s to ms for the following fitting and plotting
s_unique                = s_unique*1000;

% plot raw data
cMAP1                   = [189,215,231; 186,228,179; 252,174,145]./255;
cMAP2                   = [33,113,181; 35,139,69; 203,24,29]./255;
f1                      = figure; hold on
for i                   = 1:3
    scatter(s_unique, pResp(i,:),'MarkerFaceColor', cMAP1(i,:),...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.5)
end
xlabel('SOA (ms)')
ylabel('probability')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define the scaled psychometric function

% define PMF
P_Afirst                = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst                = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous          = @(SOA, mu, sig, lambda, c) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

% % check PMF with arbiturary parameters
% figure
% t_diff_finer = linspace(s_unique(1), s_unique(end), 1000);
% P = [-50, 50, 0.06, 100];
% P_Afirst_finer = P_Afirst(t_diff_finer, P(1), P(2), P(3), P(4));
% P_Vfirst_finer = P_Vfirst(t_diff_finer, P(1), P(2), P(3), P(4));
% P_simultaneous_finer = P_simultaneous(t_diff_finer, P(1), P(2), P(3), P(4));
% 
% figure
% plot(t_diff_finer, P_Afirst_finer, 'Color', cMAP(3,:), ...
%     'lineWidth', 3); hold on;
% plot(t_diff_finer, P_Vfirst_finer, 'Color', cMAP(1,:), ...
%     'lineWidth', 3); hold on;
% plot(t_diff_finer, P_simultaneous_finer, 'Color', cMAP(2,:), ...
%     'lineWidth', 3); hold on;
% legend({'$P(A-precedes-V)$','$P(A-follows-V)$',...
%     '$P(A-coincides-V)$'},'Location','bestoutside', 'Interpreter','Latex');

%% define cost function
% set parameters
nT_A1st                 = respCount(3,:);% the number of A-first responses for each SOA
nT_V1st                 = respCount(1,:); % the number of V-first responses for each SOA

%nLL cost function
nLogL                   = @(p) -nT_A1st*log(P_Afirst(s_unique, p(1), p(2), p(3), p(4)))' ...
    -nT_V1st*log(P_Vfirst(s_unique, p(1), p(2), p(3), p(4)))'...
    -(repmat(numTrials,size(nT_A1st)) - nT_A1st - nT_V1st)...
    * log(P_simultaneous(s_unique, p(1), p(2), p(3), p(4)))';

% set lower and upper bounds
lb                      = [-150, 10, 1e-2, 50];
ub                      = [150, 200, 0.06, 250];

% choose random initial values for 1e3 times
for i                   = 1:1e3
    init                    = rand(1,length(lb)).*(ub-lb) + lb;
    %You can also define how many times you want MATLAB to search
    options                 = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');
    
    %fmincon returns best-fitting parameters that minimize the cost function as
    %well as the corresponding value for the cost function (in this case, the
    %negative log likelihood)
    [estP(i,:), min_NLL(i)] = fmincon(nLogL, init,[],[],[],[],lb,ub,[],options);
end

% use the best-fitting parameters with the smallest NLL among 1e3 fittings
[value idx]             = min(min_NLL);
bestP                   = estP(idx,:);
fprintf('mu = %4.2f sig = %4.2f lambda = %4.2f criterion = %4.2f\n', bestP)

%% bootstrap to add error bars
numBtst                 = 1e3;
for i                   = 1:numBtst
    %initialize resampled responses
    r_slc                   = NaN(length(s_unique), numTrials);
    %for each stimulus location, we resample the responses
    for j                   = 1:length(s_unique)
        %randomly select trial indices (indices are allowed to occur more
        %than once, since we resample with replacement).
        idx                     = randi([1 numTrials],[1 numTrials]);
        %store the resampled responses
        r_slc(j,:)              = r_org(j,idx);
    end
    %compute the total number of comparison-more responses given each
    %intensity level
    nT_Vfirst_slc(i,:)      = sum(r_slc == 1,2)';
    nT_simul_slc(i,:)       = sum(r_slc == 2,2)';
    nT_Afirst_slc(i,:)      = sum(r_slc == 3,2)';
end

% convert into bootstrapped p at each SOA
btst_p                  = {nT_Vfirst_slc./numTrials;  nT_simul_slc./numTrials; nT_Afirst_slc./numTrials;};
for c                   = 1:3
    btst_resp_p             = btst_p{c};
    for x                   = 1:length(s_unique)
        [lb(c,x), ub(c,x)]      = get68CI(btst_resp_p(:,x));
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot with the best-fitting parametrs
%make a finer grid for the timing difference between the auditory and the
%visual stimulus
SOA_finer               = linspace(s_unique(1), s_unique(end), 1000);

%plug in t_diff_finer and the best-fitting mu, sigma and lapse rate to
%function P
P_Vfirst_fit            = P_Vfirst(SOA_finer, bestP(1), bestP(2), bestP(3), bestP(4));
P_Afirst_fit            = P_Afirst(SOA_finer, bestP(1), bestP(2), bestP(3), bestP(4));
P_simultaneous_fit      = P_simultaneous(SOA_finer, bestP(1), bestP(2), bestP(3), bestP(4));

f2                      = figure; hold on
plot(SOA_finer, P_Afirst_fit, 'Color', cMAP1(3,:), ...
    'lineWidth', 3); hold on;
plot(SOA_finer, P_Vfirst_fit, 'Color', cMAP1(1,:), ...
    'lineWidth', 3); hold on;
plot(SOA_finer, P_simultaneous_fit, 'Color', cMAP1(2,:), ...
    'lineWidth', 3); hold on;
h1                      = xline(bestP(1),'LineWidth',2);
h2                      = xline(bestP(1) - bestP(4),'--','LineWidth',2);
h3                      = xline(bestP(1) + bestP(4),'--','LineWidth',2);

% errorbar
for c                   = 1:3
    [hl, hp]                = boundedline(s_unique, pResp(c,:), ...
        [pResp(c,:) - lb(c,:); ub(c,:) - pResp(c,:)]',...
        '.','Color', cMAP1(c,:), 'alpha');
    hl.MarkerSize           = 20;
    hl.HandleVisibility     = 'off';
end

% display best-fitting parameters
bestP                   = round(bestP,2); 
text                    = {['$\mu = $' num2str(bestP(1))], ...
    ['$\sigma               = $' num2str(bestP(2))],...
    ['$\lambda              = $' num2str(bestP(3))],...
    ['$criterion            = $' num2str(bestP(4))]};
annotation( 'textbox', 'String', text, ...
            'FontSize', 14, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.72,0.5,0.4,0.3], 'Interpreter','Latex')

% look better
legend( {'$P_{pre}(A-precedes-V)$','$P_{pre}(A-follows-V)$',...
    '$P_{pre}(A-coincides-V)$'},'Location','bestoutside', 'Interpreter','Latex');
xlabel('SOA (ms)')
ylabel('probability')
set(f2,'Position',[10 10 900 600])
set(gca, 'FontSize', 15)
title(['sub' num2str(subjID)])

% save figure
fignm                   = ['v4_pretest_sub' num2str(subjID) '_session' num2str(sess)];
saveas(gca,fignm,'epsc')
