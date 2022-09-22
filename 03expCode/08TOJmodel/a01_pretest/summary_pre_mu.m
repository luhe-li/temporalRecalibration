%% TOJ fitting by cumulative gaussian model v5 summary

% 2020/05

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
% best-fitting parameters, and errorbars of PMF, and saved the figure.

clear all; close all; clc;


sess   = 1;
load('pre_data_summary.mat')

for subjID = 7:7
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% data-cleaning %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean raw data
% addpath(genpath(/Volumes/GoogleDrive-107780329657346345479/My Drive/2021A-22/01project/04ExperimentalCode/03ExpCpde/01pretest/data'))
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
% f1                      = figure; hold on
% for i                   = 1:3
%     scatter(s_unique, pResp(i,:),'MarkerFaceColor', cMAP1(i,:),...
%         'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.5)
% end
% xlabel('SOA (ms)')
% ylabel('probability')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define the scaled psychometric function

% define PMF
P_Afirst                = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst                = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous          = @(SOA, mu, sig, lambda, c) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

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
all_bestP(subjID,:) = bestP;
fprintf('mu = %4.2f sig = %4.2f lambda = %4.2f criterion = %4.2f\n', bestP)

%% bootstrap to add error bars

% set parameter for btst
 numBtst                 = 1e3;
% lb and ub is the same as previous fmincon

[estP_btst, minNLL, lb_68CI(subjID,:), ub_68CI(subjID,:)] = bootstrapPre(s_unique, r_org, numTrials,...
    numBtst, P_Afirst, P_Vfirst, P_simultaneous, lb, ub, options);

clearvars -except all_bestP lb_68CI ub_68CI subjID sess

end

% save data
save('pre_data_summary.mat', 'all_bestP', 'lb_68CI','ub_68CI')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure;hold on
set(f1,'Position',[0 0 300 300])
set(gca, 'FontSize', 20, 'LineWidth', 1.5)
markers = ['o','s','d','^','v','>','<'];
yneg = all_bestP(:,2) - lb_68CI(:,2);
ypos = ub_68CI(:,2) - all_bestP(:,2);
xneg = all_bestP(:,1) - lb_68CI(:,1);
xpos = ub_68CI(:,1) - all_bestP(:,1);
for sub = 1:size(all_bestP, 1)
    errorbar(all_bestP(sub,1), all_bestP(sub,2), yneg(sub), ypos(sub), xneg(sub), xpos(sub),...
        markers(sub), 'MarkerSize',15, 'Color','k')
end
% look better
ylim([0 200])
xlim([-50 100])
xlabel('PSS')
ylabel('slope')

% save data
fignm = 'summary_pre';
saveas(gca,fignm,'epsc')
