% This script fits ternary TOJ data by the exponential model in García-Pérez &
% Alcalá-Quintana (2012) 

clear all; close all; clc;

%% declare global variable
global s_unique nT_A1st nT_V1st numTrials

%% clean raw data
subjID = 1;
sess = 1;
addpath('/Volumes/GoogleDrive/My Drive/2021A-22/01project/04ExperimentalCode/03ExpCpde/01pretest/data')
load(['pretest_sub' num2str(subjID) '_session' num2str(sess) '.mat'])

s_unique = ExpInfo.SOA; % unique SOA levels, in ms
for i = 1:length(s_unique)
    iSOA = s_unique(i);
    iResp = Response.order(ExpInfo.trialSOA == iSOA);
    r_org(i,:) = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        respCount(j,i) = sum(iResp == j);
    end
end
pResp = respCount/ExpInfo.numTrials;

% convert s_unique from s to ms for the following fitting and plotting
s_unique = s_unique*1000;

% plot raw data
% cMAP = [200, 40, 40; 255, 128, 0; 13, 183, 200]./255;
% set plotting parameters
cMAP1 = [189,215,231; 186,228,179; 252,174,145]./255;
cMAP2 = [33,113,181; 35,139,69; 203,24,29]./255;

f1 = figure;
hold on
for i = 1:3
    scatter(s_unique, pResp(i,:),'MarkerFaceColor', cMAP1(i,:),...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.5)
end
xlabel('SOA (ms)')
ylabel('probability')

%% define cost function
% set parameters
nT_A1st = respCount(3,:);% the number of A-first responses for each SOA
nT_V1st = respCount(1,:); % the number of V-first responses for each SOA
numTrials = ExpInfo.numTrials;

% define PMF and cost function
% nLL = exponentialModelNLL(p, s_unique, nT_A1st, nT_V1st, numTrials)
fun = @exponentialModelNLL;

% set lower and upper bounds for
% tau, criterion, lambda_a, lambda_v, epislon, kappa
lb      = [-100, 0, 0, 0, 1e-2, 0.5];
ub      = [100, 200, 1, 1, 0.06, 0.5];

% choose random initial values for 1e3 times
for i = 1:1e3
    init    = rand(1,length(lb)).*(ub-lb) + lb;
    %You can also define how many times you want MATLAB to search
    options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');
    
    %fmincon returns best-fitting parameters that minimize the cost function as
    %well as the corresponding value for the cost function (in this case, the
    %negative log likelihood)
    [estP(i,:), min_NLL(i)] = fmincon(fun,init,[],[],[],[],lb,ub,[],options);
end

% use the best-fitting parameters with the smallest NLL among 1e3 fittings
[value idx] = min(min_NLL);
bestP = estP(idx,:);
disp(bestP) % check if bestP lies in the bounds reasonably

%% plotting with the best-fitting variable

% set plotting parameters
cMAP1 = [189,215,231; 186,228,179; 252,174,145]./255;
cMAP2 = [33,113,181; 35,139,69; 203,24,29]./255;

%make a finer grid for the timing difference between the auditory and the
%visual stimulus
SOA_finer = linspace(s_unique(1), s_unique(end), 1000);

% best-fitting parameters
tau = bestP(1);
criterion = bestP(2);
lambda_a = bestP(3);
lambda_v = bestP(4);
epsilon = bestP(5);
kappa = bestP(6);

% probability of each response
p_afirst = cumulativeD(SOA_finer, tau, -criterion, lambda_a, lambda_v);
p_vfirst = 1 - cumulativeD(SOA_finer, tau, criterion, lambda_a, lambda_v);
p_simul = cumulativeD(SOA_finer, tau, criterion, lambda_a, lambda_v)...
    - cumulativeD(SOA_finer, tau, -criterion, lambda_a, lambda_v);

% add lapse to each response probability
p_afirst_fit = add_lapse(p_afirst, p_vfirst, p_simul, epsilon, kappa);
p_vfirst_fit = add_lapse(p_vfirst, p_afirst, p_simul, epsilon, kappa);
p_simul_fit = add_lapse(p_simul, p_afirst, p_vfirst, epsilon, kappa);

figure(f1); hold on
l1 = plot(SOA_finer, p_afirst_fit, 'Color', cMAP1(3,:), ...
    'lineWidth', 3); hold on;
l2 = plot(SOA_finer, p_vfirst_fit, 'Color', cMAP1(1,:), ...
    'lineWidth', 3); hold on;
l3 = plot(SOA_finer, p_simul_fit, 'Color', cMAP1(2,:), ...
    'lineWidth', 3); hold on;

