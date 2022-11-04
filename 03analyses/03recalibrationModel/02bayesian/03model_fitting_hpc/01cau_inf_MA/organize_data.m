% This function loads the pretest and posttest TOJ data for one session,
% one subject.

function data = organize_data(subjID, sess)

%--------------------------------------------------------------------------
% Outputs:
% pResp     :proportion of v_first, simultaneous, a_first responses for
%           15 SOA levels (3x15)
% nT_V1st  :number/counts of v_first responses for 15 SOA levels
%           (1x15)
% nT_simul  :number/counts of simultaneous responses for 15 SOA levels
%           (1x15)
% nT_A1st  :number/counts of a_first responses for 15 SOA levels
%           (1x15)
%--------------------------------------------------------------------------
    
    %% %%% pre-test
    load(['pretest_sub' num2str(subjID) '_session' num2str(sess) '.mat'])
    data.pre_s_unique                  = ExpInfo.SOA; % unique SOA levels, in s
    data.pre_ms_unique                 = data.pre_s_unique * 1e3; % unique SOA levels, in ms
    data.pre_numTrials                 = ExpInfo.nTrials; % num of trials per SOA
    % inititate
    pre_r_org                     = NaN(length(data.pre_s_unique), data.pre_numTrials);
    pre_respCount                 = NaN(3, length(data.pre_s_unique));
    for i                         = 1:length(data.pre_s_unique)
        iSOA                          = data.pre_s_unique(i);
        iResp                         = Response.order(ExpInfo.trialSOA == iSOA);
        pre_r_org(i,:)                = iResp; % this matrix has a size of length(s_unique) x numTrials
        for j                         = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
            pre_respCount(j,i)            = sum(iResp == j);
        end
    end
    data.pre_pResp                     = pre_respCount/data.pre_numTrials;
    data.pre_nT_V1st                   = pre_respCount(1,:);
    data.pre_nT_simul                  = pre_respCount(2,:);
    data.pre_nT_A1st                   = pre_respCount(3,:);

    %% %%% post-test
    % load data and define key parameters
    load(['posttest_sub' num2str(subjID) '_session' num2str(sess) '.mat'])
    data.post_s_unique                 = ExpInfo.SOA; % unique SOA levels, in ms
    data.post_ms_unique                = data.post_s_unique * 1e3; % unique SOA levels, in s
    data.post_numTrials                = ExpInfo.nTrials; % num of trials per SOA
    % inititate
    post_r_org                    = NaN(length(data.post_s_unique), data.post_numTrials);
    post_respCount                = NaN(3, length(data.post_s_unique));
    for i                         = 1:length(data.post_s_unique)
        iSOA                          = data.post_s_unique(i);
        iResp                         = Response.order(ExpInfo.trialSOA == iSOA);
        post_r_org(i,:)               = iResp; % this matrix has a size of length(s_unique) x numTrials
        for j                         = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
            post_respCount(j,i)           = sum(iResp == j);
        end
    end
    data.post_pResp                    = post_respCount/data.post_numTrials;
    data.post_nT_V1st                  = post_respCount(1,:);
    data.post_nT_simul                 = post_respCount(2,:);
    data.post_nT_A1st                  = post_respCount(3,:);

    %% exposure phase
    data.adaptor_soa = ExpInfo.adaptor; % in s

end