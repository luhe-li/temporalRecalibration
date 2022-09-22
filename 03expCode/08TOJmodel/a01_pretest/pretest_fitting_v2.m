%% TOJ data analysis2
% fit multiple times to find the best-fitting parameters with the min NLL

clear all; close all; clc;

%% clean raw data
load('pretest_sub99.mat')
SOA = ExpInfo.SOA;
s_unique = SOA; % unique SOA levels
for i = 1:length(SOA)
    iSOA = SOA(i);
    iResp = Response.order(ExpInfo.trialSOA == iSOA);
    for j = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        respCount(j,i) = sum(iResp == j);
    end
end
pResp = respCount/ExpInfo.numTrials;

cMAP = [200, 40, 40; 255, 128, 0; 13, 183, 200]./255;
f1 = figure;
hold on
for i = [3,1,2]
    %     plot(SOA,pResp(i,:),'-o','Color',cMAP(i,:),...
    %         'MarkerEdgeColor',cMAP(i,:),'MarkerFaceColor', cMAP(i,:));
    scatter(SOA, pResp(i,:),'MarkerFaceColor', cMAP(i,:),...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.5)
    
end
xlabel('SOA (ms)')
ylabel('probability')

%% define the (unscaled/scaled) psychometric function and the cost function
%P is a cumulative gaussian distribution. It takes three inputs
%1. x: the timing difference between the auditory and the visual stimulus
%2. mu: the center of the psychometric function (a.k.a., PSS,
%       perceived subjective simultaneity). Past studies have shown that
%       PSS is not equal to 0, as sounds are normally perceived faster than
%       visual stimuli by ~60ms. In other words, participants perceive
%       an auditory and a visual stimulus as simultaneous when the auditory
%       stimulus is delayed by 60ms.
%3. sig: the measurement noise regarding the temporal difference

%the scaled cumulative gaussian distribution by a lapse rate
P_Afirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous = @(SOA, mu, sig, lambda, c) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

% % plot models with random set of parameters to check
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

%% fitting the scaled psychometric function by minimizing the cost function
% set parameters
nT_A1st = respCount(3,:);% the number of A-first responses for each SOA
nT_V1st = respCount(1,:); % the number of V-first responses for each SOA
numTrials = ExpInfo.numTrials;

%nLL cost function
nLogL = @(p) -nT_A1st*log(P_Afirst(s_unique, p(1), p(2), p(3), p(4)))' ...
    -nT_V1st*log(P_Vfirst(s_unique, p(1), p(2), p(3), p(4)))'...
    -(repmat(numTrials,size(nT_A1st)) - nT_A1st - nT_V1st)...
    * log(P_simultaneous(s_unique, p(1), p(2), p(3), p(4)))';

%To find the combination of mu, sigma and lapse rate that minimizes the
%negative log likelihood, we will use fmincon.m in MATLAB. To use this
%function, we need to define lower and upper bounds for each parameter
%(i.e., search space) as well as an initial point for MATLAB to start
%searching.
lb      = [-150, 10, 1e-2, 50];
ub      = [150, 200, 0.06, 250];

for i = 1:1e3
    init    = rand(1,length(lb)).*(ub-lb) + lb;
    %You can also define how many times you want MATLAB to search
    options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');
    
    %fmincon returns best-fitting parameters that minimize the cost function as
    %well as the corresponding value for the cost function (in this case, the
    %negative log likelihood)
    [estP(i,:), min_NLL(i)] = fmincon(nLogL, init,[],[],[],[],lb,ub,[],options);
end

[value idx] = min(min_NLL);
bestP = estP(idx,:);
disp(bestP)

%% Plot with the best-fitting parametrs
%make a finer grid for the timing difference between the auditory and the
%visual stimulus
t_diff_finer = linspace(s_unique(1), s_unique(end), 1000);
%plug in t_diff_finer and the underlying true mu, sigma and lapse rate to
%function Ptilde_finer
P_Vfirst_fit = P_Vfirst(t_diff_finer, bestP(1), bestP(2), bestP(3), bestP(4));
P_Afirst_fit = P_Afirst(t_diff_finer, bestP(1), bestP(2), bestP(3), bestP(4));
P_simultaneous_fit = P_simultaneous(t_diff_finer, bestP(1), bestP(2), bestP(3), bestP(4));

figure(f1)
plot(t_diff_finer, P_Afirst_fit, 'Color', cMAP(3,:), ...
    'lineWidth', 3); hold on;
plot(t_diff_finer, P_Vfirst_fit, 'Color', cMAP(1,:), ...
    'lineWidth', 3); hold on;
plot(t_diff_finer, P_simultaneous_fit, 'Color', cMAP(2,:), ...
    'lineWidth', 3); hold on;
h1= xline(bestP(1),'LineWidth',2);
h2 = xline(bestP(1) - bestP(4),'--','LineWidth',2);
h3 = xline(bestP(1) + bestP(4),'--','LineWidth',2);
legend({'$P(A-precedes-V)$','$P(A-follows-V)$',...
    '$P(A-coincides-V)$'},'Location','bestoutside', 'Interpreter','Latex');

fignm = ['pretest_sub', num2str(ExpInfo.subjID)];
saveas(gca,fignm,'png')
