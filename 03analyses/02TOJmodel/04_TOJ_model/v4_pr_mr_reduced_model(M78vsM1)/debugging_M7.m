%% debug M7_NLL codes. Parameters are set the same so that the NLL shouldn't differ much between M1 and M7

subjID = 3; fix_lambda = 0.3; num_sess = 9;

% define PMF
P_Afirst                                   = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst                                   = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous                             = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

p1 = [-10, 40, 70, 70, 200, 200];
M1_nll=0;

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
    pre_nT_A1st                   = pre_respCount(3,:);
    pre_nT_V1st                   = pre_respCount(1,:);
    pre_nT_simul                  = pre_respCount(2,:);

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
    post_nT_A1st                  = post_respCount(3,:);
    post_nT_V1st                  = post_respCount(1,:);
    post_nT_simul                 = post_respCount(2,:);

    %% specify models inside the loop
    M1                            = @(p) -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(3), p(5), fix_lambda))'...
        -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(3), p(5), fix_lambda))'...
        -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(3), p(5), fix_lambda))'...
        -post_nT_A1st*log(P_Afirst(post_ms_unique, p(2), p(4), p(6), fix_lambda))' ...
        -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(2), p(4), p(6), fix_lambda))'...
        -post_nT_simul*log(P_simultaneous(post_ms_unique, p(2), p(4), p(6), fix_lambda))';
    
    M1_nll = M1_nll+M1(p1);


end

display(M1_nll)

clearvars post_nT_A1st post_nT_V1st post_nT_simul pre_nT_A1st pre_nT_V1st pre_nT_simul

subjID = 3; fix_lambda = 0.3; num_sess = 9;
p2 = [-10, 40, 40, 40, 40, 40, 40, 40, 40, 40, 70, 200];

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

post_mu_idx                                = [2:10];

NLL                                        = @(r1, r2, r3, r4, r5, r6, p) ...
    -r1*log(P_Afirst(pre_ms_unique, p(1), p(3), p(4), fix_lambda))'...
    -r2*log(P_Vfirst(pre_ms_unique, p(1), p(3), p(4), fix_lambda))'...
    -r3*log(P_simultaneous(pre_ms_unique, p(1), p(3), p(4), fix_lambda))'...
    -r4*log(P_Afirst(post_ms_unique, p(2), p(3), p(4), fix_lambda))'...
    -r5*log(P_Vfirst(post_ms_unique, p(2), p(3), p(4), fix_lambda))'...
    -r6*log(P_simultaneous(post_ms_unique, p(2), p(3), p(4), fix_lambda))';

M7_NLL = @(p) sum(arrayfun(@(idx) NLL(pre_nT_A1st{idx}, pre_nT_V1st{idx}, pre_nT_simul{idx}, ...
    post_nT_A1st{idx}, post_nT_V1st{idx}, post_nT_simul{idx}, ...
    [p(1), p(idx+1), p(11), p(12)]), 1:num_sess));

disp(M7_NLL(p2))