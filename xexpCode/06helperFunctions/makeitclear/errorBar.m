% calculate errorbar for p in psychometric function

%% 1 SD
% closed form: sd = sqrt(p*(1-p)/N)
% p: 1-d array, p at different levels of x
% N: number of trials at x
sd = arrayfun(@(p) sqrt(p*(1-p)/post_numTrials), post_pResp);

%% bootstrap 1 SD
% bootstrap with replacement, obtain 68% of confidence interval
% r_org: a matrix of N x X, raw response of N trials at X levels
numBtst = 1e3;
for i = 1:numBtst
    %initialize resampled responses
    r_slc = NaN(length(X), N);
    %for each stimulus location, we resample the responses
    for j = 1:length(X)
        %randomly select trial indices (indices are allowed to occur more
        %than once, since we resample with replacement).
        idx        = randi([1 N],[1 N]);
        %store the resampled responses
        r_slc(j,:) = r_org(j,idx);
    end
    %compute the total number of comparison-more responses given each
    %intensity level
    nT_Vfirst_slc(i,:) = sum(r_slc == 1,2)';
    nT_simultaneous_slc(i,:) = sum(r_slc == 2,2)';
    nT_Afirst_slc(i,:) = sum(r_slc == 3,2)';
end

% convert into p
btst_p = {nT_Vfirst_slc./N;  nT_simultaneous_slc./N; nT_Afirst_slc./N};
for c = 1:3 % loop through 3 responses
    each_p = btst_p{c};
    for x = 1:length(X)
        [pre_lb(c,x), pre_ub(c,x)] = get68CI(each_p(:,x));
    end
end
