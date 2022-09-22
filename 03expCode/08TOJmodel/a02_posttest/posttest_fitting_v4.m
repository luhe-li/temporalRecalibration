%% TOJ fitting v4

% 2020/04

% This script fits **cumulative Gaussian PMF** for the ternary TOJ data.

clear all; close all; clc;

%% input

subjID = input('Please enter participant ID#: ') ;
sess = input('Please enter participant session#: ') ;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% data-cleaning %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pretest: raw data
% clean raw data
addpath('/Volumes/GoogleDrive-107780329657346345479/My Drive/2021A-22/01project/04ExperimentalCode/03ExpCpde/01pretest/data')
load(['pretest_sub' num2str(subjID) '_session' num2str(sess) '.mat'])

pre_s_unique = ExpInfo.SOA; % unique SOA levels, in ms
pre_numTrials = ExpInfo.nTrials;
for i = 1:length(pre_s_unique)
    iSOA = pre_s_unique(i);
    iResp = Response.order(ExpInfo.trialSOA == iSOA);
    pre_r_org(i,:) = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        pre_respCount(j,i) = sum(iResp == j);
    end
end
pre_pResp = pre_respCount/pre_numTrials;

% convert s_unique from s to ms for the following fitting and plotting
pre_s_unique = pre_s_unique*1000;

% plot raw data
cMAP1 = [189,215,231; 186,228,179; 252,174,145]./255;
cMAP2 = [33,113,181; 35,139,69; 203,24,29]./255;
alpha = 0.5; % transparency
f1 = figure;
hold on
for i = 1:3
    scatter(pre_s_unique+5, pre_pResp(i,:),60,'MarkerFaceColor', cMAP1(i,:),...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',alpha)
end
xlabel('SOA (ms)')
ylabel('probability')

%% posttest: raw data
addpath('/Volumes/GoogleDrive-107780329657346345479/My Drive/2021A-22/01project/04ExperimentalCode/03ExpCpde/04posttest/data')
load(['posttest_sub' num2str(subjID) '_session' num2str(sess) '.mat'])

post_s_unique = ExpInfo.SOA; % unique SOA levels, in ms
post_numTrials = ExpInfo.nTrials;
for i = 1:length(post_s_unique)
    iSOA = post_s_unique(i);
    iResp = Response.order(ExpInfo.trialSOA == iSOA);
    post_r_org(i,:) = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        post_respCount(j,i) = sum(iResp == j);
    end
end
post_pResp = post_respCount/post_numTrials;

% convert s_unique from s to ms for the following fitting and plotting
post_s_unique = post_s_unique*1000;

% plot raw data
figure(f1)
hold on
for i = 1:3
    scatter(post_s_unique-5, post_pResp(i,:),60,'MarkerFaceColor', cMAP2(i,:),...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',alpha)
end
xlabel('SOA (ms)')
ylabel('probability')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% fittings %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define the scaled psychometric function

% define PMF
P_Afirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous = @(SOA, mu, sig, lambda, c) ...
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
pre_nT_A1st = pre_respCount(3,:);% the number of A-first responses for each SOA
pre_nT_V1st = pre_respCount(1,:); % the number of V-first responses for each SOA

post_nT_A1st = post_respCount(3,:);% the number of A-first responses for each SOA
post_nT_V1st = post_respCount(1,:); % the number of V-first responses for each SOA

%nLL cost function
nLogL = @(p) -pre_nT_A1st*log(P_Afirst(pre_s_unique, p(1), p(3), p(4), p(5)))' ...
    -pre_nT_V1st*log(P_Vfirst(pre_s_unique, p(1), p(3), p(4), p(5)))'...
    -(repmat(pre_numTrials,size(pre_nT_A1st)) - pre_nT_A1st - pre_nT_V1st)...
    * log(P_simultaneous(pre_s_unique, p(1), p(3), p(4), p(5)))'...
    -post_nT_A1st*log(P_Afirst(post_s_unique, p(2), p(3), p(4), p(5)))' ...
    -post_nT_V1st*log(P_Vfirst(post_s_unique, p(2), p(3), p(4), p(5)))'...
    -(repmat(post_numTrials,size(post_nT_A1st)) - post_nT_A1st - post_nT_V1st)...
    * log(P_simultaneous(post_s_unique, p(2), p(3), p(4), p(5)))';

% set lower and upper bounds
lb      = [-150, -150, 10, 1e-2, 50];
ub      = [150, 150, 200, 0.06, 250];

% choose random initial values for 1e3 times
for i = 1:1e3
    init    = rand(1,length(lb)).*(ub-lb) + lb;
    %You can also define how many times you want MATLAB to search
    options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');

    %fmincon returns best-fitting parameters that minimize the cost function as
    %well as the corresponding value for the cost function (in this case, the
    %negative log likelihood)
    [estP(i,:), min_NLL(i)] = fmincon(nLogL, init,[],[],[],[],lb,ub,[],options);
end

% use the best-fitting parameters with the smallest NLL among 1e3 fittings
[value idx] = min(min_NLL);
bestP = estP(idx,:);
fprintf('pre mu = %4.2f post mu = %4.2f sig = %4.2f lambda = %4.2f criterion = %4.2f\n', bestP)

%% pretest: bootstrap to add error bars
numBtst = 1e3;
for i = 1:numBtst
    %initialize resampled responses
    pre_r_slc = NaN(length(pre_s_unique), pre_numTrials);
    %for each stimulus location, we resample the responses
    for j = 1:length(pre_s_unique)
        %randomly select trial indices (indices are allowed to occur more
        %than once, since we resample with replacement).
        idx        = randi([1 pre_numTrials],[1 pre_numTrials]);
        %store the resampled responses
        pre_r_slc(j,:) = pre_r_org(j,idx);
    end
    %compute the total number of comparison-more responses given each
    %intensity level
    pre_nT_Vfirst_slc(i,:) = sum(pre_r_slc == 1,2)';
    pre_nT_simultaneous_slc(i,:) = sum(pre_r_slc == 2,2)';
    pre_nT_Afirst_slc(i,:) = sum(pre_r_slc == 3,2)';
end

% convert into p
pre_btst_p = {pre_nT_Vfirst_slc./pre_numTrials;  pre_nT_simultaneous_slc./pre_numTrials; pre_nT_Afirst_slc./pre_numTrials};
for c = 1:3
    post_p = pre_btst_p{c};
    for x = 1:length(pre_s_unique)
        [pre_lb(c,x), pre_ub(c,x)] = get68CI(post_p(:,x));
    end
end

%% posttest: bootstrap to add error bars
numBtst = 1e3;
for i = 1:numBtst
    %initialize resampled responses
    post_r_slc = NaN(length(post_s_unique), post_numTrials);
    %for each stimulus location, we resample the responses
    for j = 1:length(post_s_unique)
        %randomly select trial indices (indices are allowed to occur more
        %than once, since we resample with replacement).
        idx        = randi([1 post_numTrials],[1 post_numTrials]);
        %store the resampled responses
        post_r_slc(j,:) = post_r_org(j,idx);
    end
    %compute the total number of comparison-more responses given each
    %intensity level
    post_nT_Vfirst_slc(i,:) = sum(post_r_slc == 1,2)';
    post_nT_simultaneous_slc(i,:) = sum(post_r_slc == 2,2)';
    post_nT_Afirst_slc(i,:) = sum(post_r_slc == 3,2)';
end

% convert to p
post_btst_p = {post_nT_Vfirst_slc./post_numTrials;  post_nT_simultaneous_slc./post_numTrials; post_nT_Afirst_slc./post_numTrials};
for c = 1:3
    post_p = post_btst_p{c};
    for x = 1:length(post_s_unique)
        [post_lb(c,x), post_ub(c,x)]  = get68CI(post_p(:,x));
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pre-test: plot with the best-fitting parametrs

% set plotting parameters
cMAP1 = [189,215,231; 186,228,179; 252,174,145]./255;
cMAP2 = [33,113,181; 35,139,69; 203,24,29]./255;
if ExpInfo.adaptor <= 0
    jitter = 5;
else
    jitter = -5;
end

%make a finer grid for the timing difference between the auditory and the
%visual stimulus
t_diff_finer = linspace(pre_s_unique(1), pre_s_unique(end), 1000);

%plug in t_diff_finer and the best-fitting mu, sigma and lapse rate to
%function P
P_Vfirst_fit = P_Vfirst(t_diff_finer, bestP(1), bestP(3), bestP(4), bestP(5));
P_Afirst_fit = P_Afirst(t_diff_finer, bestP(1), bestP(3), bestP(4), bestP(5));
P_simultaneous_fit = P_simultaneous(t_diff_finer, bestP(1), bestP(3), bestP(4), bestP(5));

f2 = figure; hold on;
l1 = plot(t_diff_finer, P_Afirst_fit, 'Color', [cMAP1(3,:)], ...
    'lineWidth', 3); hold on;
l2 = plot(t_diff_finer, P_Vfirst_fit, 'Color', [cMAP1(1,:)], ...
    'lineWidth', 3); hold on;
l3 = plot(t_diff_finer, P_simultaneous_fit, 'Color', [cMAP1(2,:)], ...
    'lineWidth', 3); hold on;
h1 = xline(bestP(1),'LineWidth',2);
h2 = xline(bestP(1) - bestP(5),'-','LineWidth',2);
h3 = xline(bestP(1) + bestP(5),'-','LineWidth',2);
% errorbar
for c = 1:3
    [hl hp]= boundedline(pre_s_unique + jitter, pre_pResp(c,:), ...
        [pre_pResp(c,:) - pre_lb(c,:); pre_ub(c,:) - pre_pResp(c,:)]',...
        '.','Color', cMAP1(c,:), 'alpha');
    hl.MarkerSize = 20;
    hl.HandleVisibility = 'off';
end

%% post-test: plot with the best-fitting parametrs
%make a finer grid for the timing difference between the auditory and the
%visual stimulus
t_diff_finer = linspace(post_s_unique(1), post_s_unique(end), 1000);

%plug in t_diff_finer and the best-fitting mu, sigma and lapse rate to
%function P
P_Vfirst_fit = P_Vfirst(t_diff_finer, bestP(2), bestP(3), bestP(4), bestP(5));
P_Afirst_fit = P_Afirst(t_diff_finer, bestP(2), bestP(3), bestP(4), bestP(5));
P_simultaneous_fit = P_simultaneous(t_diff_finer, bestP(2), bestP(3), bestP(4), bestP(5));

figure(f2)
l4 = plot(t_diff_finer, P_Afirst_fit, '-', 'Color', cMAP2(3,:), ...
    'lineWidth', 3); hold on;
l5 = plot(t_diff_finer, P_Vfirst_fit, '-','Color', cMAP2(1,:), ...
    'lineWidth', 3); hold on;
l6 = plot(t_diff_finer, P_simultaneous_fit, '-','Color', cMAP2(2,:), ...
    'lineWidth', 3); hold on;
h4 = xline(bestP(2),'--','LineWidth',2);
h5 = xline(bestP(2) - bestP(5),'--','LineWidth',2);
h6 = xline(bestP(2) + bestP(5),'--','LineWidth',2);
% errorbar
for c = 1:3
    [hl, hp] = boundedline(post_s_unique-jitter, post_pResp(c,:), ...
        [post_pResp(c,:) - post_lb(c,:); post_ub(c,:) - post_pResp(c,:)]',...
        '.','Color', cMAP2(c,:), 'alpha');
    hl.MarkerSize = 20;
    hl.HandleVisibility = 'off';
end

% display best-fitting parameters
bestP = round(bestP,2); 
text = {['$\mu_{pre} = $' num2str(bestP(1))], ...
    ['$\mu_{post} = $' num2str(bestP(2))],...
    ['$\sigma = $' num2str(bestP(3))],...
    ['$\lambda = $' num2str(bestP(4))],...
    ['$criterion = $' num2str(bestP(5))]};
annotation( 'textbox', 'String', text, ...
            'FontSize', 14, 'Units', 'normalized', 'EdgeColor', 'none', ...
            'Position', [0.72,0.5,0.4,0.2], 'Interpreter','Latex')

% look better
legend([l1, l2, l3, l4, l5, l6], {'$P_{pre}(A-precedes-V)$','$P_{pre}(A-follows-V)$',...
    '$P_{pre}(A-coincides-V)$', '$P_{post}(A-precedes-V)$','$P_{post}(A-follows-V)$',...
    '$P_{post}(A-coincides-V)$'},'Location','bestoutside', 'Interpreter','Latex');

xlabel('SOA (ms)')
ylabel('probability')
set(f2,'Position',[10 10 900 600])
set(gca, 'FontSize', 15)
title(['sub' num2str(subjID) ' adaptor SOA = ' num2str(ExpInfo.adaptor*1000) ' ms'])
fignm = ['posttest_sub', num2str(subjID) '_session' num2str(sess)];

saveas(gca,fignm,'epsc')
