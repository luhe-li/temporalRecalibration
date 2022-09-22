function [sum_AIC, sum_NLL, estP]     = model_comp_fit_each_session_v4(subjID, fix_lambda)

% v4

% This function fits 6 models.
% Each session is fitted independently, delta_AIC is summed across sessions for each model, each subject

% output
% sum_AIC       : a vector
% sum_NLL       : a vector
% estP          : a cell (size = num_session x num_model), each element is fitted
% parameter (a vector)

%% model parameters outside the loop
% define PMF
P_Afirst                      = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst                      = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous                = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

numM                          = 6; %number of models
n_init                        = 10;
numP                          = [6,5,5,4,4,3]; %number of free parameters for M1, M2, M3, M4, M5 respectively
num_sess                      = 9; %number of sessions

%define upper and lower bounds
lb                            = {[-200, -200, 10, 10, 10, 10], ...%M1
    [-200, -200, 10, 10, 10] ...%M2
    [-200, -200, 10, 10, 10], ...%M3
    [-200, -200, 10, 10], ...%M4
    [-200, 10, 10, 10], ...%M5
    [-200, 10, 10]};%M6
ub                            = {[200, 200, 350, 350, 450, 450], ...%M1
    [200, 200, 350, 450, 450], ...%M2
    [200, 200, 350, 350, 450], ...%M3
    [200, 200, 350, 450], ...%M4
    [200, 350, 450, 450], ...%M5
    [200, 350, 450]};%M6
init_fun                      = @(a,b) rand(1,length(a)).*(b-a) + a;
options                       = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');
estP                          = cell(num_sess,numM); % 9 sessions x 6 models
[min_NLL, AIC]                = deal(NaN(num_sess, numM)); % 9 sessions x 6 models

%% organize data inside the loop
for sess                      = 1:9
    %%%%% pre-test
    load(['pretest_sub' num2str(subjID) '_session' num2str(sess) '.mat'])
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
    pre_pResp                     = pre_respCount/pre_numTrials;
    pre_nT_V1st                   = pre_respCount(1,:);
    pre_nT_simul                  = pre_respCount(2,:);
    pre_nT_A1st                   = pre_respCount(3,:);

    %%%%% post-test
    % load data and define key parameters
    load(['posttest_sub' num2str(subjID) '_session' num2str(sess) '.mat'])
    post_s_unique                 = ExpInfo.SOA; % unique SOA levels, in ms
    post_ms_unique                = post_s_unique * 1e3; % unique SOA levels, in s
    post_numTrials                = ExpInfo.nTrials; % num of trials per SOA
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
    post_pResp                    = post_respCount/post_numTrials;
    post_nT_V1st                  = post_respCount(1,:);
    post_nT_simul                 = post_respCount(2,:);
    post_nT_A1st                  = post_respCount(3,:);

    %% specify models inside the loop
    M1                            = @(p)  -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(3), p(5), fix_lambda))'...
        -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(3), p(5), fix_lambda))'...
        -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(3), p(5), fix_lambda))'...
        -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(2), p(4), p(6), fix_lambda))'...
        -post_nT_simul*log(P_simultaneous(post_ms_unique, p(2), p(4), p(6), fix_lambda))'...
        -post_nT_A1st*log(P_Afirst(post_ms_unique, p(2), p(4), p(6), fix_lambda))' ;

    M2                            = @(p) -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(3), p(4), fix_lambda))'...
        -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(3), p(4), fix_lambda))'...
        -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(3), p(4), fix_lambda))' ...
        -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(2), p(3), p(5), fix_lambda))'...
        -post_nT_simul*log(P_simultaneous(post_ms_unique, p(2), p(3), p(5), fix_lambda))'...
        -post_nT_A1st*log(P_Afirst(post_ms_unique, p(2), p(3), p(5), fix_lambda))' ;

    M3                            = @(p) -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(3), p(5), fix_lambda))'...
        -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(3), p(5), fix_lambda))'...
        -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(3), p(5), fix_lambda))' ...
        -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(2), p(4), p(5), fix_lambda))'...
        -post_nT_simul*log(P_simultaneous(post_ms_unique, p(2), p(4), p(5), fix_lambda))'...
        -post_nT_A1st*log(P_Afirst(post_ms_unique, p(2), p(4), p(5), fix_lambda))';

    M4                            = @(p) -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(3), p(4), fix_lambda))'...
        -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(3), p(4), fix_lambda))'...
        -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(3), p(4), fix_lambda))' ...
        -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(2), p(3), p(4), fix_lambda))'...
        -post_nT_simul*log(P_simultaneous(post_ms_unique, p(2), p(3), p(4), fix_lambda))'...
        -post_nT_A1st*log(P_Afirst(post_ms_unique, p(2), p(3), p(4), fix_lambda))' ;

    M5                            = @(p) -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(2), p(3), fix_lambda))'...
        -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(2), p(3), fix_lambda))'...
        -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(2), p(3), fix_lambda))' ...
        -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(1), p(2), p(4), fix_lambda))'...
        -post_nT_simul*log(P_simultaneous(post_ms_unique, p(1), p(2), p(4), fix_lambda))'...
        -post_nT_A1st*log(P_Afirst(post_ms_unique, p(1), p(2), p(4), fix_lambda))';

    M6                            = @(p) -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(2), p(3), fix_lambda))'...
        -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(2), p(3), fix_lambda))'...
        -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(2), p(3), fix_lambda))' ...
        -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(1), p(2), p(3), fix_lambda))'...
        -post_nT_simul*log(P_simultaneous(post_ms_unique, p(1), p(2), p(3), fix_lambda))'...
        -post_nT_A1st*log(P_Afirst(post_ms_unique, p(1), p(2), p(3), fix_lambda))' ;

    NLLs                          = {M1; M2; M3; M4; M5; M6};

    %loop through models
    for m                         = 1:numM
        for ii                        = 1:n_init
            %the initial point for matlab to start searching
            init                          = init_fun(lb{m}, ub{m});
            %use fmincon.m to fit
            [ii_estP{ii}, ii_min_NLL(ii)] = fmincon(NLLs{m}, init,[],[],[],[],...
                lb{m}, ub{m},[],options);
        end
        % % use the best-fitting parameters with the smallest NLL among ninit fittings
        [min_NLL(sess, m) idx]        = min(ii_min_NLL);
        estP{sess, m}                 = ii_estP{idx};
        %compute the AIC/BIC
        AIC(sess, m)                  = 2*min_NLL(sess, m) + 2*numP(m);
        %     BIC(m) = 2*min_NLL(m) + numP(m)*log(pre_numTrials);
    end
end

%% add AIC across sessions and plot
sum_AIC                         = sum(AIC, 1);
sum_NLL = sum(min_NLL, 1);
% deltaAIC                      = s_AIC - min(s_AIC);

end