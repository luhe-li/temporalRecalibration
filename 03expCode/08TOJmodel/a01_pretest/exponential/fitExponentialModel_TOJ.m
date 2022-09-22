%% TOJ fitting by exponential model v4

% 2020/04

% This script fits ternary TOJ data by the exponential model in García-Pérez &
% Alcalá-Quintana (2012) 

clear all; close all; clc;
subjID                  = 3;
sess                    = 1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% data-cleaning %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean raw data
addpath('/Volumes/GoogleDrive/My Drive/2021A-22/01project/04ExperimentalCode/03ExpCpde/01pretest/data')
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
% cMAP = [200, 40, 40; 255, 128, 0; 13, 183, 200]./255;
% set plotting parameters
cMAP1                   = [189,215,231; 186,228,179; 252,174,145]./255;
cMAP2                   = [33,113,181; 35,139,69; 203,24,29]./255;

f1                      = figure;
hold on
for i                   = 1:3
    scatter(s_unique, pResp(i,:),'MarkerFaceColor', cMAP1(i,:),...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.5)
end
xlabel('SOA (ms)')
ylabel('probability')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define cost function
% set parameters
nT_A1st                 = respCount(3,:);% the number of A-first responses for each SOA
nT_V1st                 = respCount(1,:); % the number of V-first responses for each SOA

% define PMF and cost function
% nLL = exponentialModelNLL(p, s_unique, nT_A1st, nT_V1st, numTrials)

% set lower and upper bounds for
% tau, criterion, lambda_a, lambda_v, epislon, kappa, where
lb                      = [-100, 0, 0, 0, 1e-2, 0.5];
ub                      = [100, 200, 100, 100, 0.06, 0.5];

% choose random initial values for 1e3 times
for i                   = 1:1e3
    init                    = rand(1,length(lb)).*(ub-lb) + lb;
    %You can also define how many times you want MATLAB to search
    options                 = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');
    
    %fmincon returns best-fitting parameters that minimize the cost function as
    %well as the corresponding value for the cost function (in this case, the
    %negative log likelihood)
    [estP(i,:), min_NLL(i)] = fmincon(@(p) exponentialModelNLL(p, s_unique, nT_A1st, nT_V1st, numTrials),init,[],[],[],[],lb,ub,[],options);
end

% use the best-fitting parameters with the smallest NLL among 1e3 fittings
[value idx]             = min(min_NLL);
bestP                   = estP(idx,:);
disp(bestP) % check if bestP lies in the bounds reasonably
fprintf('tau            = %4.2f criterion = %4.2f mu_a = %4.2f mu_v = %4.2f epsilon = %4.2f\n kappa = %4.2f\n', bestP)

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
%% plotting with the best-fitting variable

% set plotting parameters
cMAP1                   = [189,215,231; 186,228,179; 252,174,145]./255;
cMAP2                   = [33,113,181; 35,139,69; 203,24,29]./255;

%make a finer grid for the timing difference between the auditory and the
%visual stimulus
SOA_finer               = linspace(s_unique(1), s_unique(end), 1000);

% best-fitting parameters
bestP                   = round(bestP,2); 
tau                     = bestP(1);
criterion               = bestP(2);
mu_a                = bestP(3);
mu_v                = bestP(4);
epsilon                 = bestP(5);
kappa                   = bestP(6);

% probability of each response
p_afirst                = cumulativeD(SOA_finer, tau, -criterion, 1/mu_a, 1/mu_v);
p_vfirst                = 1 - cumulativeD(SOA_finer, tau, criterion, 1/mu_a, 1/mu_v);
p_simul                 = cumulativeD(SOA_finer, tau, criterion, 1/mu_a, 1/mu_v)...
    - cumulativeD(SOA_finer, tau, -criterion, 1/mu_a, 1/mu_v);

% add lapse to each response probability
p_afirst_fit            = add_lapse(p_afirst, p_vfirst, p_simul, epsilon, kappa);
p_vfirst_fit            = add_lapse(p_vfirst, p_afirst, p_simul, epsilon, kappa);
p_simul_fit             = add_lapse(p_simul, p_afirst, p_vfirst, epsilon, kappa);

f2                      = figure; hold on
l1                      = plot(SOA_finer, p_afirst_fit, 'Color', cMAP1(3,:), ...
    'lineWidth', 3); hold on;
l2                      = plot(SOA_finer, p_vfirst_fit, 'Color', cMAP1(1,:), ...
    'lineWidth', 3); hold on;
l3                      = plot(SOA_finer, p_simul_fit, 'Color', cMAP1(2,:), ...
    'lineWidth', 3); hold on;

% add vertical lines
mu = mu_a - mu_v + tau;
h1                      = xline(mu,'LineWidth',2); % mean
h2                      = xline(mu - criterion,'--','LineWidth',2); % lower criterion
h3                      = xline(mu + criterion,'--','LineWidth',2); % higher criterion


% errorbar
for c                   = 1:3
    [hl, hp]                = boundedline(s_unique, pResp(c,:), ...
        [pResp(c,:) - lb(c,:); ub(c,:) - pResp(c,:)]',...
        '.','Color', cMAP1(c,:), 'alpha');
    hl.MarkerSize           = 20;
    hl.HandleVisibility     = 'off';
end

% display best-fitting parameters
text                    = {['$\tau = $' num2str(tau)], ...
    ['$criterion            = $' num2str(criterion)],...
    ['$\mu_{a}          = $' num2str(mu_a)],...
    ['$\mu_{v}          = $' num2str(mu_v)],...
    ['$\epsilon             = $' num2str(epsilon)],...
    ['$\kappa               = $' num2str(kappa)]};
annotation( 'textbox', 'String', text, ...
            'FontSize', 16, 'Units', 'normalized', 'EdgeColor', 'none', ...
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
fignm                   = ['pretest_sub' num2str(subjID) '_session' num2str(sess) '_model2'];
saveas(gca,fignm,'epsc')