function [AIC, min_NLL, estP] = model_comp_fit_all_sessions_v4(subjID, fix_lambda)

% v4

% Fit M7: all sessions are fitted together, delta_AIC is summed across
% sessions for each model, each subject

% output
% AIC       : a scalar
% min_NLL   : a scalar
% estP      : a cell (size = num_session x num_model), each element is fitted
%             parameter (a vector)

%% model parameters outside the loop
% define PMF
P_Afirst                                   = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst                                   = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous                             = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

numM                                       = 1; %number of models
n_init                                     = 10;
numP                                       = 12; %number of free parameters
num_sess                                   = 9; %number of sessions

%define upper and lower bounds
lb                                         = {[repmat(-200,[1,10]), 10, 10]}; %M7
ub                                         = {[repmat(200,[1,10]), 350, 400]}; %M7
init_fun                                   = @(a,b) rand(1,length(a)).*(b-a) + a;
options                                    = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');
estP                                       = cell(1,numM); % fit all sessions together
[min_NLL, AIC]                             = deal(NaN(1, numM)); % fit all sessions together

%% organize data inside the loop

% create an empty cell (size: 1x9) for the number of response, pre, post,
% in each session
[pre_nT_A1st, pre_nT_V1st, pre_nT_simul, ...
    post_nT_A1st, post_nT_V1st, post_nT_simul] = deal(cell(1, num_sess)); 

for s                                      = 1:num_sess
    %%%%% pre-test
    load(['pretest_sub' num2str(subjID) '_session' num2str(s) '.mat'])
    pre_s_unique                               = ExpInfo.SOA; % unique SOA levels, in s
    pre_ms_unique                              = pre_s_unique * 1e3; % unique SOA levels, in ms
    pre_numTrials                              = ExpInfo.nTrials; % num of trials per SOA
    % inititate
    pre_r_org                                  = NaN(length(pre_s_unique), pre_numTrials);
    pre_respCount                              = NaN(3, length(pre_s_unique));
    for i                                      = 1:length(pre_s_unique)
        iSOA                                       = pre_s_unique(i);
        iResp                                      = Response.order(ExpInfo.trialSOA == iSOA);
        pre_r_org(i,:)                             = iResp; % this matrix has a size of length(s_unique) x numTrials
        for j                                      = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
            pre_respCount(j,i)                         = sum(iResp == j);
        end
    end
    %     pre_pResp                     = pre_respCount/pre_numTrials;
    pre_nT_A1st{s}                             = pre_respCount(3,:);
    pre_nT_V1st{s}                             = pre_respCount(1,:);
    pre_nT_simul{s}                            = pre_respCount(2,:);

    %%%%% post-test
    % load data and define key parameters
    load(['posttest_sub' num2str(subjID) '_session' num2str(s) '.mat'])
    post_s_unique                              = ExpInfo.SOA; % unique SOA levels, in ms
    post_ms_unique                             = post_s_unique * 1e3; % unique SOA levels, in s
    post_numTrials                             = ExpInfo.nTrials; % num of trials per SOA
    % inititate
    post_r_org                                 = NaN(length(post_s_unique), post_numTrials);
    post_respCount                             = NaN(3, length(post_s_unique));
    for i                                      = 1:length(post_s_unique)
        iSOA                                       = post_s_unique(i);
        iResp                                      = Response.order(ExpInfo.trialSOA == iSOA);
        post_r_org(i,:)                            = iResp; % this matrix has a size of length(s_unique) x numTrials
        for j                                      = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
            post_respCount(j,i)                        = sum(iResp == j);
        end
    end
    %     post_pResp                    = post_respCount/post_numTrials;
    post_nT_A1st{s}                            = post_respCount(3,:);
    post_nT_V1st{s}                            = post_respCount(1,:);
    post_nT_simul{s}                           = post_respCount(2,:);

end

%% specify models outside the loop

NLL                                        = @(r1, r2, r3, r4, r5, r6, p) ...
    -r1*log(P_Vfirst(pre_ms_unique, p(1), p(3), p(4), fix_lambda))'...
    -r2*log(P_simultaneous(pre_ms_unique, p(1), p(3), p(4), fix_lambda))'...
    -r3*log(P_Afirst(pre_ms_unique, p(1), p(3), p(4), fix_lambda))'...
    -r4*log(P_Vfirst(post_ms_unique, p(2), p(3), p(4), fix_lambda))'...
    -r5*log(P_simultaneous(post_ms_unique, p(2), p(3), p(4), fix_lambda))'...
    -r6*log(P_Afirst(post_ms_unique, p(2), p(3), p(4), fix_lambda))';


M7_NLL = @(p) sum(arrayfun(@(idx) NLL(pre_nT_V1st{idx}, pre_nT_simul{idx}, pre_nT_A1st{idx}, ...
    post_nT_V1st{idx}, post_nT_simul{idx}, post_nT_A1st{idx}, ...
    [p(1), p(idx+1), p(11), p(12)]), 1:num_sess));

%loop through models
for m                                      = 1:numM
    for ii                                     = 1:n_init
        %the initial point for matlab to start searching
        init                                       = init_fun(lb{m}, ub{m});
        %use fmincon.m to fit
        [ii_estP{ii}, ii_min_NLL(ii)]              = fmincon(M7_NLL, init,[],[],[],[],...
            lb{m}, ub{m},[],options);
    end
    % % use the best-fitting parameters with the smallest NLL among ninit fittings
    [min_NLL idx]                              = min(ii_min_NLL);
    estP                                       = ii_estP{idx};
    %compute the AIC/BIC
    AIC                                        = 2*min_NLL + 2*numP;
    %     BIC(m) = 2*min_NLL(m) + numP(m)*log(pre_numTrials);
end

end