% This function returns the raw percentage of each response and fitted percentage of each response
% with specified model and corresponding best-fitting parameters

function  [pre_pResp, post_pResp, pre_fit, post_fit] = model_comp_plot_each_session_v4(subjID, sess, idxModel, bestP)

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

% %make a finer grid for the timing difference between the auditory and the
% %visual stimulus
% SOA_finer = linspace(pre_ms_unique(1), pre_ms_unique(end), 1000);

% define PMF
P_Afirst = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

M1 = {@(p) P_Afirst(pre_ms_unique, p(1), p(3), p(5), p(7));...
    @(p) P_Vfirst(pre_ms_unique, p(1), p(3), p(5), p(7));...
    @(p) P_simultaneous(pre_ms_unique, p(1), p(3), p(5), p(7));...
    @(p) P_Afirst(post_ms_unique, p(2), p(4), p(6), p(8));...
    @(p) P_Vfirst(post_ms_unique, p(2), p(4), p(6), p(8));...
    @(p) P_simultaneous(post_ms_unique, p(2), p(4), p(6), p(8))};

M2 = {@(p) P_Afirst(pre_ms_unique, p(1), p(3), p(4), p(6));...
    @(p) P_Vfirst(pre_ms_unique, p(1), p(3), p(4), p(6));...
    @(p) P_simultaneous(pre_ms_unique, p(1), p(3), p(4), p(6));...
    @(p) P_Afirst(post_ms_unique, p(2), p(3), p(5), p(7));...
    @(p) P_Vfirst(post_ms_unique, p(2), p(3), p(5), p(7));...
    @(p) P_simultaneous(post_ms_unique, p(2), p(3), p(5), p(7))};

M3 = {@(p) P_Afirst(pre_ms_unique, p(1), p(3), p(5), p(6));...
    @(p) P_Vfirst(pre_ms_unique, p(1), p(3), p(5), p(6));...
    @(p) P_simultaneous(pre_ms_unique, p(1), p(3), p(5), p(6));...
    @(p) P_Afirst(post_ms_unique, p(2), p(4), p(5), p(7));...
    @(p) P_Vfirst(post_ms_unique, p(2), p(4), p(5), p(7));...
    @(p) P_simultaneous(post_ms_unique, p(2), p(4), p(5), p(7))};

M4 = {@(p) P_Afirst(pre_ms_unique, p(1), p(3), p(4), p(5));...
    @(p) P_Vfirst(pre_ms_unique, p(1), p(3), p(4), p(5));...
    @(p) P_simultaneous(pre_ms_unique, p(1), p(3), p(4), p(5));...
    @(p) P_Afirst(post_ms_unique, p(2), p(3), p(4), p(6));...
    @(p) P_Vfirst(post_ms_unique, p(2), p(3), p(4), p(6));...
    @(p) P_simultaneous(post_ms_unique, p(2), p(3), p(4), p(6))};

M5 = {@(p) P_Afirst(pre_ms_unique, p(1), p(2), p(3), p(5));...
    @(p) P_Vfirst(pre_ms_unique, p(1), p(2), p(3), p(5));...
    @(p) P_simultaneous(pre_ms_unique, p(1), p(2), p(3), p(5));...
    @(p) P_Afirst(post_ms_unique, p(1), p(2), p(4), p(6));...
    @(p) P_Vfirst(post_ms_unique, p(1), p(2), p(4), p(6));...
    @(p) P_simultaneous(post_ms_unique, p(1), p(2), p(4), p(6))};

M6 = {@(p) P_Afirst(pre_ms_unique, p(1), p(2), p(3), p(4));...
    @(p) P_Vfirst(pre_ms_unique, p(1), p(2), p(3), p(4));...
    @(p) P_simultaneous(pre_ms_unique, p(1), p(2), p(3), p(4));...
    @(p) P_Afirst(post_ms_unique, p(1), p(2), p(3), p(5));...
    @(p) P_Vfirst(post_ms_unique, p(1), p(2), p(3), p(5));...
    @(p) P_simultaneous(post_ms_unique, p(1), p(2), p(3), p(5))};

models = {M1; M2; M3; M4; M5; M6};

%% plug in best-fitting parameters and plot

% get fitting by the model
for c = 1:3
    pre_fit(i,:) = models{idxModel}{i}(bestP); % 1: A; 2: V; 3: simul
    post_fit(i+3, :) = models{idxModel}{i+3}(bestP); % 1: A; 2: V; 3: simul
end

end