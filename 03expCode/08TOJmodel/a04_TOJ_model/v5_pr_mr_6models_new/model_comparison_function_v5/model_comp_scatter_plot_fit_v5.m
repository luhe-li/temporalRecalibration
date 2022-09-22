% This function plots the raw data with model fitting with a specified
% model its corresponding best-fitting parameters

% inputs:
% sub       : a scalar, subject id (e.g., 8)
% sess      : a scalar, session number (e.g., 9)
% imodel    : a scalar, index of one specific model (1-6)
% bestP     : a vector (size: 1 x num_parameters), must correspond with the
%           specified model (M1-M6 have different number of parameters)

function [data, pmf_fit] = model_comp_scatter_plot_fit_v5(sub, ses, iModel, bestP)

%% organize data

%%%%% pre-test
% load data and define key parameters
load(['pretest_sub' num2str(sub) '_session' num2str(ses) '.mat'])
pre_s_unique             = ExpInfo.SOA; % unique SOA levels, in s
pre_ms_unique            = pre_s_unique * 1e3; % unique SOA levels, in ms
pre_numTrials            = ExpInfo.nTrials; % num of trials per SOA
% inititate
pre_r_org                = NaN(length(pre_s_unique), pre_numTrials);
pre_respCount            = NaN(3, length(pre_s_unique));
for i                    = 1:length(pre_s_unique)
    iSOA                     = pre_s_unique(i);
    iResp                    = Response.order(ExpInfo.trialSOA == iSOA);
    pre_r_org(i,:)           = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j                    = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        pre_respCount(j,i)       = sum(iResp == j);
    end
end
pre_pResp                = pre_respCount/pre_numTrials;

%%%%% post-test
% load data and define key parameters
load(['posttest_sub' num2str(sub) '_session' num2str(ses) '.mat'])
post_s_unique            = ExpInfo.SOA; % unique SOA levels, in ms
post_ms_unique           = post_s_unique * 1e3; % unique SOA levels, in s
post_numTrials           = ExpInfo.nTrials; % num of trials per SOA
% inititate
post_r_org               = NaN(length(post_s_unique), post_numTrials);
post_respCount           = NaN(3, length(post_s_unique));
for i                    = 1:length(post_s_unique)
    iSOA                     = post_s_unique(i);
    iResp                    = Response.order(ExpInfo.trialSOA == iSOA);
    post_r_org(i,:)          = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j                    = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        post_respCount(j,i)      = sum(iResp == j);
    end
end
post_pResp               = post_respCount/post_numTrials;

data                     = [pre_pResp; post_pResp];

%% define and select best model

% grid is evaluated at 15 test_SOAs
SOA_finer                = pre_ms_unique;

% define PMF
P_Afirst                 = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst                 = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous           = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

M1                       = {@(p) P_Vfirst(SOA_finer, p(1), p(3), p(5), p(7));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(5), p(7));...
    @(p) P_Afirst(SOA_finer, p(1), p(3), p(5), p(7));...
    @(p) P_Vfirst(SOA_finer, p(2), p(4), p(6), p(8));...
    @(p) P_simultaneous(SOA_finer, p(2), p(4), p(6), p(8));...
    @(p) P_Afirst(SOA_finer, p(2), p(4), p(6), p(8))};

M2                       = {@(p) P_Vfirst(SOA_finer, p(1), p(3), p(4), p(6));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(4), p(6));...
    @(p) P_Afirst(SOA_finer, p(1), p(3), p(4), p(6));...
    @(p) P_Vfirst(SOA_finer, p(2), p(3), p(5), p(7));...
    @(p) P_simultaneous(SOA_finer, p(2), p(3), p(5), p(7));...
    @(p) P_Afirst(SOA_finer, p(2), p(3), p(5), p(7))};

M3                       = {@(p) P_Vfirst(SOA_finer, p(1), p(3), p(5), p(6));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(5), p(6));...
    @(p) P_Afirst(SOA_finer, p(1), p(3), p(5), p(6));...
    @(p) P_Vfirst(SOA_finer, p(2), p(4), p(5), p(7));...
    @(p) P_simultaneous(SOA_finer, p(2), p(4), p(5), p(7));...
    @(p) P_Afirst(SOA_finer, p(2), p(4), p(5), p(7))};

M4                       = {@(p) P_Vfirst(SOA_finer, p(1), p(3), p(4), p(5));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(4), p(5));...
    @(p) P_Afirst(SOA_finer, p(1), p(3), p(4), p(5));...
    @(p) P_Vfirst(SOA_finer, p(2), p(3), p(4), p(6));...
    @(p) P_simultaneous(SOA_finer, p(2), p(3), p(4), p(6));...
    @(p) P_Afirst(SOA_finer, p(2), p(3), p(4), p(6))};

M5                       = {@(p) P_Vfirst(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_Afirst(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_Vfirst(SOA_finer, p(1), p(2), p(4), p(6));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(4), p(6));...
    @(p) P_Afirst(SOA_finer, p(1), p(2), p(4), p(6))};

M6                       = {@(p) P_Vfirst(SOA_finer, p(1), p(2), p(3), p(4));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(3), p(4));...
    @(p) P_Afirst(SOA_finer, p(1), p(2), p(3), p(4));...
    @(p) P_Vfirst(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_Afirst(SOA_finer, p(1), p(2), p(3), p(5))};

models                   = {M1; M2; M3; M4; M5; M6};

%% fit PMF

pmf_fit                  = arrayfun(@(i) models{iModel}{i}(bestP), 1:6,'UniformOutput',false);

end