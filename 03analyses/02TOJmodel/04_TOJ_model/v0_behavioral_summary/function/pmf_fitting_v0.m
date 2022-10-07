% This function fits each session with the full model, assuming all 8
% parameters.

% b) to include bootstrapping, use:
function [bestP, idx_adaptor_soa, fit_nll, estP_btst, btst_nll, lb_68CI, ub_68CI] = ...
    pmf_fitting_v0(pre_filenm, post_filenm)

% %a) to exclude bootstrapping, use:
% function [bestP, fit_nll, idx_adaptor_soa] = ...
%     pmf_fitting_v0(pre_filenm, post_filenm)

%% clean data

%%%%% pre-test
load(pre_filenm);
pre_s_unique                  = ExpInfo.SOA; % unique SOA levels, in s
pre_ms_unique                 = pre_s_unique * 1e3; % unique SOA levels, in ms
pre_numTrials                 = ExpInfo.nTrials; % num of trials per SOA
% inititate
pre_r_org                     = NaN(length(pre_s_unique), pre_numTrials);
pre_respCount                 = NaN(3, length(pre_s_unique));
for i                         = 1:length(pre_s_unique)
    iSOA                          = pre_s_unique(i);
    iResp                         = Response.order(ExpInfo.trialSOA == iSOA);
    pre_r_org(i,:)                = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j                         = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        pre_respCount(j,i)            = sum(iResp == j);
    end
end
% pre_pResp                             = pre_respCount/pre_numTrials;
pre_nT_V1st                   = pre_respCount(1,:);
pre_nT_simul                  = pre_respCount(2,:);
pre_nT_A1st                   = pre_respCount(3,:);


%%%%% post-test
% load data and define key parameters
load(post_filenm);
post_s_unique                 = ExpInfo.SOA; % unique SOA levels, in ms
post_ms_unique                = post_s_unique * 1e3; % unique SOA levels, in s
post_numTrials                = ExpInfo.nTrials; % num of trials per SOA
adaptor_soa                   = ExpInfo.adaptor; % in s
adaptors = [-700, -300:100:300, 700]/1000;
idx_adaptor_soa = find(adaptors == adaptor_soa);

% inititate
post_r_org                    = NaN(length(post_s_unique), post_numTrials);
post_respCount                = NaN(3, length(post_s_unique));
for i                         = 1:length(post_s_unique)
    iSOA                          = post_s_unique(i);
    iResp                         = Response.order(ExpInfo.trialSOA == iSOA);
    post_r_org(i,:)               = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j                         = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        post_respCount(j,i)           = sum(iResp == j);
    end
end
% post_pResp                            = post_respCount/post_numTrials;
post_nT_V1st                  = post_respCount(1,:);
post_nT_simul                 = post_respCount(2,:);
post_nT_A1st                  = post_respCount(3,:);

%% PMF fitting

% define the scaled psychometric function
P_Afirst                      = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst                      = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous                = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

% define cost function using the full model, assuming 8 parameters
nLogL                         = @(p) -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(3), p(5), p(7)))'...
    -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(3), p(5), p(7)))'...
    -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(3), p(5), p(7)))'...
    -post_nT_A1st*log(P_Afirst(post_ms_unique, p(2), p(4), p(6), p(8)))' ...
    -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(2), p(4), p(6), p(8)))'...
    -post_nT_simul*log(P_simultaneous(post_ms_unique, p(2), p(4), p(6), p(8)))';

% set lower and upper bounds
% mu_pre, mu_post, sigma_pre, sigma_post, criterion_pre, criterion_post,
% lambda_pre, lambda_post
lb                            = [-250, -250, 10, 10, 10, 10, 1e-4, 1e-4];
ub                            = [250, 250, 350, 350, 350, 350, 0.06, 0.06];

% choose random initial values for 100 times
for i                         = 1:100
    init                          = rand(1,length(lb)).*(ub-lb) + lb;
    %You can also define how many times you want MATLAB to search
    options                       = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');

    %fmincon returns best-fitting parameters that minimize the cost function as
    %well as the corresponding value for the cost function (in this case, the
    %negative log likelihood)
    [estP(i,:), min_NLL(i)]       = fmincon(nLogL, init,[],[],[],[],lb,ub,[],options);
    %     disp(i)
end

% use the best-fitting parameters with the smallest NLL among 1e3 fittings
[fit_nll idx]                   = min(min_NLL);
bestP                         = estP(idx, :);

% to include bootstrapping, use:
%% bootstrap

numBtst = 1e3;
[estP_btst, btst_nll, lb_68CI, ub_68CI] = bootstrap_pmf_fitting(pre_ms_unique, post_ms_unique, ...
    pre_r_org, post_r_org, pre_numTrials, post_numTrials,...
    numBtst, P_Afirst, P_Vfirst, P_simultaneous, lb, ub, options) ;

end