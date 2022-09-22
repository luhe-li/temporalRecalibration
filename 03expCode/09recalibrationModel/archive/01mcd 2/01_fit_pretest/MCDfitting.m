%% TOJ fitting by MCD model

clear all; clc; close all;

%% organize data

pre_filenm = 'pretest_sub1_session1.mat';
load(pre_filenm);

% extract key parameters
s_unique_sec = ExpInfo.SOA; % unique SOA levels, in sec
s_unique = s_unique_sec*1000; % in ms
num_soa = length(s_unique); % number of SOA
numTrials = ExpInfo.nTrials; % number of trials per SOA level

% initiate response vectors
r_org = NaN(num_soa , numTrials); % organized responses, col = SOA, row = response at each trial
respCount = NaN(3, num_soa); % number of first-responses for V, simul, A at each SOA level

for i = 1:length(s_unique_sec)
    iSOA = s_unique_sec(i);
    iResp = Response.order(ExpInfo.trialSOA == iSOA);
    r_org(i,:) = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        respCount(j,i) = sum(iResp == j);
    end
end

% convert number of first-responses to proportion
pResp = respCount/numTrials;

% extract count and proportion of first-responses by conditions
nT_afirst = respCount(3,:); nT_vfirst = respCount(1,:); nT_simul = respCount(2,:); 
p_afirst = pResp(3,:); p_vfirst = pResp(1,:); p_simul = pResp(2,:);

%% play around with random parameters
% 
% % define fixed parameters, see parameter details in MCD_prob.m
% duration = 16; % in s
% fs = 1e3; % hz
% onset = 9; % in sec
% stim_dura = 0.033; % in sec
% 
% % randonly select free parameters
% bestP = [0.0873, 0.0684, 0.7859, 0.1, 0.1, 0.1];
% [p_afirst, p_vfirst, p_simul] = MCD_prob (s_unique_sec, duration, fs, onset, stim_dura, ...
%     bestP(1), bestP(2), bestP(3), bestP(4), bestP(5), bestP(6));
% 
% figure; hold on
% plot(s_unique, [p_afirst; p_vfirst; p_simul] , '-ok')

%% fit model

% define fixed parameters, see parameter details in MCD_prob.m
duration = 16; % in s
fs = 1e3; % hz
onset = 9; % in sec
stim_dura = 0.033; % in sec

% define nLL cost function
nll = @(p) MCD_nll(s_unique_sec, duration, fs, onset, stim_dura, ... % fixed parameters
    nT_vfirst, nT_simul, nT_afirst,... % behavioral responses to be fitted
    p); % free parameters

% set lower and upper bounds for
% tv, ta, tav, mu, sigma, criterion
lb      = [0.001, 0.001, 0.1,  -1, 0.01, 0.01];
ub      = [  0.1,   0.1,   1,   1,    1,    1];

% choose random initial values for 1e3 times
for i = 1:1 % this should be 1e3 once we figured out the best bounds
    init    = rand(1,length(lb)).*(ub-lb) + lb;
    %You can also define how many times you want MATLAB to search
    options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');

    %fmincon returns best-fitting parameters that minimize the cost function as
    %well as the corresponding value for the cost function (in this case, the
    %negative log likelihood)
    [estP(i,:), min_NLL(i)] = fmincon(nll, init,[],[],[],[],lb,ub,[],options);
end

% use the best-fitting parameters with the smallest NLL among 1e3 fittings
[value idx] = min(min_NLL);
bestP = estP(idx,:);
fprintf('tau_v = %4.2f tau_a = %4.2f tau_av = %4.2f mu = %4.2f sigma = %4.2f\n criterion = %4.2f\n', bestP)

%% plotting

% fit with best-fitting parameters and finer grids

soa_finer = round(linspace(s_unique_sec(1), s_unique_sec(end), 1e3),3);
[p_afirst_fit, p_vfirst_fit, p_simul_fit] = MCD_prob(-soa_finer, duration, fs, onset, stim_dura, ...
    bestP(1), bestP(2), bestP(3), bestP(4), bestP(5), bestP(6));

% color parameters
cMAP1 = [158,202,225; 161,217,155; 252,174,145]./255;
alpha = 0.5;

% plot fitting curve
figure; hold on
set(gca,'FontSize',15,'linewidth',1.5)
box off
plot(-soa_finer, p_afirst_fit, 'Color',cMAP1(3,:),'LineWidth',2);
plot(-soa_finer, p_vfirst_fit,'Color',cMAP1(1,:),'LineWidth',2);
plot(-soa_finer, p_simul_fit,'Color',cMAP1(2,:),'LineWidth',2);

% plot raw data points
for i = 1:3
scatter(-s_unique_sec, pResp(i,:),60,'MarkerFaceColor', cMAP1(i,:),'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',alpha)
end
xlabel('SOA (s)'); xticks([-s_unique_sec(1), 0, -s_unique_sec(end)]);
ylabel('probability')