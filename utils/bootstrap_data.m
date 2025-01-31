function btstData = bootstrap_data(sub)

for ses = 1:9

    data = organize_data(sub, ses);

    %% pre-test

    %initialize resampled responses
    pre_r_btst = NaN(length(data.pre_ms_unique), data.pre_numTrials);

    %for each stimulus location, we resample the responses
    for j = 1:length(data.pre_ms_unique)
        %randomly select trial indices (indices are allowed to occur more
        %than once, since we resample with replacement).
        idx        = randi([1 data.pre_numTrials],[1 data.pre_numTrials]);
        %store the resampled responses
        pre_r_btst(j,:) = data.pre_r_org(j,idx);
    end
    %compute the total number of comparison-more responses given each
    %intensity level
    b_Data.pre_nT_V1st=  sum(pre_r_btst == 1,2)';
    b_Data.pre_nT_simul = sum(pre_r_btst == 2,2)';
    b_Data.pre_nT_A1st = sum(pre_r_btst == 3,2)';
    b_Data.pre_r_org = pre_r_btst;
    b_Data.pre_s_unique = data.pre_s_unique;
    b_Data.pre_ms_unique = data.pre_ms_unique;
    b_Data.pre_numTrials = data.pre_numTrials;
    b_Data.pre_pResp = [b_Data.pre_nT_V1st; b_Data.pre_nT_simul; b_Data.pre_nT_A1st]./b_Data.pre_numTrials;

    %% post-test

    %initialize resampled responses
    post_r_btst = NaN(length(data.post_ms_unique), data.post_numTrials);
    %for each stimulus location, we resample the responses
    for j = 1:length(data.post_ms_unique)
        %randomly select trial indices (indices are allowed to occur more
        %than once, since we resample with replacement).
        idx        = randi([1 data.post_numTrials],[1 data.post_numTrials]);
        %store the resampled responses
        post_r_btst(j,:) = data.post_r_org(j,idx);
    end
    %compute the total number of comparison-more responses given each
    %intensity level
    b_Data.post_nT_V1st=  sum(post_r_btst == 1,2)';
    b_Data.post_nT_simul = sum(post_r_btst == 2,2)';
    b_Data.post_nT_A1st = sum(post_r_btst == 3,2)';
    b_Data.post_r_org = post_r_btst;
    b_Data.post_s_unique = data.post_s_unique;
    b_Data.post_ms_unique = data.post_ms_unique;
    b_Data.post_numTrials = data.post_numTrials;
    b_Data.post_pResp = [b_Data.post_nT_V1st; b_Data.post_nT_simul; b_Data.post_nT_A1st]./b_Data.post_numTrials;
    b_Data.adaptor_soa = data.adaptor_soa;

    btstData(ses) = b_Data;

end