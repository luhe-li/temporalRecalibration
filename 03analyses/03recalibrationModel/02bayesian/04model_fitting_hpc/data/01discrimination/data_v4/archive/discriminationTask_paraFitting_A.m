%% discrimination task: parameter fitting
% This script obtains auditory JND at 75% accuracy.

% In data-cleaning, it extracts two key outputs:
% r_org :  this matrix has a size of length(s_unique) x numTrials
% respCount: this matrix has the size of 2(type of response, 1 =
% standard-more, 2 = comparison more) x length(s_unique)

% In fitting, we used log_weibull to fit volume in log scale, fit the PMF use
% fmincon for 1e3 times and chose the best-fitting parameters based on the
% smallest nLL.

% Data were then bootsrapped for 1e3 times to calculate the CI of
% parameters.

% In plotting, we plotted the raw data, the fitted PMF with the
% best-fitting parameters, and errorbars(?)

clear all; close all; clc;

%% define psychometric function

% define PMF
% PF_cumGaussian = @(x, gamma, lambda, mu, sig) gamma + (1 - lambda - gamma).* normcdf(x, mu, sig);
PF_log_weibull = @(x, gamma, lambda, alpha, beta) gamma + (1 - lambda - gamma).* (1 - exp(-10.^(beta.*(x-alpha))));
PF_weibull = @(x, gamma, lambda, alpha, beta) gamma + (1 - lambda - gamma).* (1 - exp( -(x./alpha).^beta));

% check PMF with arbiturary parameters
% int_finer = linspace(s_unique(1), s_unique(end), 1000);
% P =[0.5, 0.05, 0.7, 4];
% P_comp_more = PF_log_weibull(int_finer, P(1), P(2), P(3), P(4));
% figure; hold on
% plot(int_finer, P_comp_more)
% ylim([0.5 1.0])
% ylabel('proportion of perceiving comparison stimulus longer')
% xlabel('comparison volumn')
% xticks(round(levels,2))

% set parameter for both modalities to loop
M = {'discrimination_sub99_A.mat', 'discrimination_sub99_V.mat';
    PF_log_weibull, PF_weibull;
    'auditory JND', 'visual JND'};

%% clean raw data

load('discrimination_sub99_A.mat')

% % fake response
% Response.more = Shuffle([ones(1,70), 2*ones(1,70)]);
% ExpInfo.VolumeRange = logspace(log10(0.5), log10(1), 7);

levels = ExpInfo.VolumeRange;
s_unique = levels; % unique intensity levels

% convert response from "which order is more" to "standard(1) or comparison(2) is
% more"
trialType = {[1,2];[2,1]}; % 1 = standard first; 2 = comparison first
cbTrials = reshape([trialType{ExpInfo.counterbalance}],2, [])';
for t = 1:ExpInfo.numTotalTrials
    RespComp(t) = cbTrials(t, Response.more(t)); % 1 = chose standard more; 2 = chose comparison more
end

% organize response based on intensity level
for i = 1:length(levels)
    iLevel = levels(i);
    iResp = RespComp(AudInfo.shuffledIntensity == iLevel);
    r_org(i,:) = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j = unique(RespComp) % 1 = chose standard more; 2 = chose comparison more
        respCount(j,i) = sum(iResp == j);
    end
end
pResp = respCount/ExpInfo.numTrials;

% plot raw data
cMAP = [215,48,39; 252,141,89; 254,224,144; 24,243,248; 145,191,219; 69,117,180]./255;
f(m) = figure; hold on
l1 = scatter(s_unique, pResp(2,:),40,'MarkerFaceColor', cMAP(6,:),...
    'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.5);

%% define psychometric function

% define PMF
% PF_cumGaussian = @(x, gamma, lambda, mu, sig) gamma + (1 - lambda - gamma).* normcdf(x, mu, sig);

% check PMF with arbiturary parameters
% int_finer = linspace(s_unique(1), s_unique(end), 1000);
% P =[0.5, 0.05, 0.7, 4];
% P_comp_more = PF_log_weibull(int_finer, P(1), P(2), P(3), P(4));
% figure; hold on
% plot(int_finer, P_comp_more)
% ylim([0.5 1.0])
% ylabel('proportion of perceiving comparison stimulus longer')
% xlabel('comparison volumn')
% xticks(round(levels,2))

%% define cost function

% set parameters
nT_compMore = respCount(2,:); % number of comparison stimulus more responses for each level
nT_standardMore = respCount(1,:); % number of standard stimulus more responses for each level
numTrials = ExpInfo.numTrials;

% nLL cost function
nLL = @(p) -nT_compMore * log(M{2,m}(s_unique, p(1), p(2), p(3), p(4)))'...
    -nT_standardMore * log(1 - M{2,m}(s_unique, p(1), p(2), p(3), p(4)))';

% set lower and upper bounds
lb      = [0.5, 0, 0.5, 0];
ub      = [0.5, 0.06, 1, 4];

% choose random initial values for 1e3 times
for i = 1:1e3
    init(i,:)    = rand(1,length(lb)).*(ub-lb) + lb;
    %You can also define how many times you want MATLAB to search
    options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');

    %fmincon returns best-fitting parameters that minimize the cost function as
    %well as the corresponding value for the cost function (in this case, the
    %negative log likelihood)
    [estP(i,:), min_NLL(i)] = fmincon(nLL, init(i,:),[],[],[],[],lb,ub,[],options);
end

% use the best-fitting parameters with the smallest NLL among 1e3 fittings
[value idx] = min(min_NLL);
bestP = estP(idx,:);
disp(bestP) % check if bestP lies in the bounds reasonably

%% bootstrap to obtain error bars

% definee bounds for each parameter: gamma, lambda, alpha, beta
b_lb = [ 0, 0, 0, 0];
b_ub = [1, 0.06, 1, 10];
b_init = rand(1,length(b_lb)).*(b_ub-b_lb) + b_lb;
b_options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');

numBtst = 1e3;
[estP_btst, minNLL, lb_95CI, ub_95CI, nT_compMore_slc] = BootstrapDiscrimination...
    (s_unique, r_org, numTrials, numBtst, M{2,m}, b_lb, b_ub, b_options);

% % compute 95% CI for bootstrapped comparison-more response at each intensity
% % level
% for j = 1:length(s_unique)
%     [CI_lb(j), CI_ub(j), eb(j)] = get95CI(nT_compMore_slc(:,j));
% end

%% plot with best-fitting parameter
figure(f(m)); hold on
set(gca,'FontSize',15,'linewidth',2)
int_finer = linspace(s_unique(1), s_unique(end), 1000);
P = bestP;
P_comp_more = PF_log_weibull(int_finer, P(1), P(2), P(3), P(4));
l2 = plot(int_finer, P_comp_more,  'Color', cMAP(6,:), 'LineWidth',3);
% errorbar(s_unique, pResp(2,:), eb, 'Color', cMAP(6,:),'MarkerSize',4, 'LineWidth',1);
%  look better

ylim([0.4 1.0])
ylabel('p(comparison longer)')
xlabel('comparison volume')
title(M{3,m})
xticks(round(levels,2))
legend([l1, l2],{'data','fit'},'Location','northwest')
%     saveas(f(m), M{3,m}, 'eps')

%% obtain JND at x.75
[value idx]  = min(abs(P_comp_more - 0.75));

x75(sub, m) = int_finer(idx);
save('JND', x75)
save(M{3,m})