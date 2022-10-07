% boostrap p at each SOA to plot data error bar with 68% confidence
% interval

%% pre-test

numBtst                 = 1e3;
for i                   = 1:numBtst
    %initialize resampled responses
    r_slc                   = NaN(length(s_unique), numTrials);
    %for each stimulus location, we resample the responses
    for j                   = 1:length(s_unique)
        %randomly select trial indices (indices are allowed to occur more
        %than once, since we resample with replacement).
        idx                     = randi([1 numTrials],[1 numTrials]);
        %store the resampled responses
        r_slc(j,:)              = r_org(j,idx);
    end
    %compute the total number of comparison-more responses given each
    %intensity level
    nT_Vfirst_slc(i,:)      = sum(r_slc == 1,2)';
    nT_simul_slc(i,:)       = sum(r_slc == 2,2)';
    nT_Afirst_slc(i,:)      = sum(r_slc == 3,2)';
end

% convert into bootstrapped p at each SOA
btst_p                  = {nT_Vfirst_slc./numTrials;  nT_simul_slc./numTrials; nT_Afirst_slc./numTrials;};
for c                   = 1:3
    btst_resp_p             = btst_p{c};
    for x                   = 1:length(s_unique)
        [lb(c,x), ub(c,x)]      = get68CI(btst_resp_p(:,x));
    end
end

% plotting
%make a finer grid for the timing difference between the auditory and the
%visual stimulus
SOA_finer               = linspace(s_unique(1), s_unique(end), 1000);

%plug in t_diff_finer and the best-fitting mu, sigma and lapse rate to
%function P
P_Vfirst_fit            = P_Vfirst(SOA_finer, bestP(1), bestP(2), bestP(3), bestP(4));
P_Afirst_fit            = P_Afirst(SOA_finer, bestP(1), bestP(2), bestP(3), bestP(4));
P_simultaneous_fit      = P_simultaneous(SOA_finer, bestP(1), bestP(2), bestP(3), bestP(4));

f2                      = figure; hold on
plot(SOA_finer, P_Afirst_fit, 'Color', cMAP1(3,:), ...
    'lineWidth', 3); hold on;
plot(SOA_finer, P_Vfirst_fit, 'Color', cMAP1(1,:), ...
    'lineWidth', 3); hold on;
plot(SOA_finer, P_simultaneous_fit, 'Color', cMAP1(2,:), ...
    'lineWidth', 3); hold on;
h1                      = xline(bestP(1),'LineWidth',2);
h2                      = xline(bestP(1) - bestP(4),'--','LineWidth',2);
h3                      = xline(bestP(1) + bestP(4),'--','LineWidth',2);

for c                   = 1:3
    [hl, hp]                = boundedline(s_unique, pResp(c,:), ...
        [pResp(c,:) - lb(c,:); ub(c,:) - pResp(c,:)]',...
        'Color', cMAP1(c,:), 'alpha');
    hl.MarkerSize           = 20;
    hl.HandleVisibility     = 'off';
end
