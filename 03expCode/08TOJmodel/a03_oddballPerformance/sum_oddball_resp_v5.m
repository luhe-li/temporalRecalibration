%% sum oddball response v5

% This script summarize dprime averaged between exposure phase and
% post-test phase, for each session, for each individual in expo_d and
% post_d.

% It then plots group averaged auditory and visual dprime by session, for
% both the exposure and post-phase.

% This version uses function **summarizeBlockPerformance2.m**

clear all; clc; close all;

%% extract individual dprime for each session for a/v
subjIDs = 3:10;
sessions = 1:9;

[expo_d, post_d] = deal(NaN(length(subjIDs), length(sessions), 2));

for i = 1:length(subjIDs)
    subjID = subjIDs(i);

    for j = 1:length(sessions)
        sess = sessions(j);

        load(['exposure_sub' num2str(subjID) '_session' num2str(sess) '.mat']);
        [expo_hitRate, expo_FA, expo_btst_d, expo_lb_68CI, expo_ub_68CI] = ...
            summarizeBlockPerformance2(1, ExpInfo.nTTrials, ExpInfo.idxOddballA, ExpInfo.idxOddballV, Response.oddball);

        load(['posttest_sub' num2str(subjID) '_session' num2str(sess) '.mat']);
        [post_hitRate, post_FA, post_btst_d, post_lb_68CI, post_ub_68CI] = ...
            summarizeBlockPerformance2(1, ExpInfo.expoNTTrials, ExpInfo.idxOddballA, ExpInfo.idxOddballV, ExpoResponse.oddball);

        for k = 1:2 % 1 = A; 2 = V
            expo_d(i, j, k) = calculate_dprime(expo_hitRate(k), expo_FA(k));
            post_d(i, j, k) = calculate_dprime(post_hitRate(k), post_FA(k));
        end

    end
end


%% summarize group data

% replace infinite data point by the max value for later averaging/std
max_d1 = max(expo_d(isfinite(expo_d)),[],'all');
expo_d(expo_d == inf) = max_d1;
max_d2 = max(post_d(isfinite(post_d)),[],'all');
post_d(post_d == inf) = max_d2;

% mean and sd
g_expo_d = squeeze(mean(expo_d, 1));
g_se_expo_d = squeeze(std(expo_d, 1))./sqrt(length(subjIDs));
g_post_d = squeeze(mean(post_d, 1));
g_se_post_d = squeeze(std(post_d, 1))./sqrt(length(subjIDs));

%% plotting

label = {'auditory'; 'visual'};

% exposure phase
figure; hold on
set(gca, 'FontSize', 20, 'LineWidth', 1.5)
box off
for k = 1:2
    errorbar([sessions + 0.2 *(k - 1.5)], g_expo_d(:,k), g_se_expo_d(:,k),'-o','LineWidth', 1.5)
    xticks(sessions)
    xticklabels(sessions)
    ylabel('d''')
    xlabel('sessions')
    ylim([0 5])
    title('exposure phase')
    legend({'auditory oddball','visual oddball'})
    expo_outlier{k} = isoutlier(expo_d(:, :, k), 'mean');
end
fignm = ['expo_nsub' num2str(length(subjIDs))];
saveas(gca,fignm,'epsc')

% post-test phase
figure; hold on
set(gca, 'FontSize', 20, 'LineWidth', 1.5)
box off
for k = 1:2
    errorbar([sessions + 0.2 *(k - 1.5)], g_post_d(:,k), g_se_expo_d(:,k),'-o','LineWidth', 1.5)
    xticks(sessions)
    xticklabels(sessions)
    ylabel('d''')
    xlabel('sessions')
    ylim([0 5])
    title('post-test phase')
    legend({'auditory oddball','visual oddball'})
    post_outlier{k} = isoutlier(expo_d(:, :, k), 'mean');
end
fignm = ['post_nsub' num2str(length(subjIDs))];
saveas(gca,fignm,'epsc')


%% function part
function d = calculate_dprime(iHR, iFA)
if iHR == 1
    iHR = 1 - (1/2*(iHR+iFA));
elseif iHR == 0
    iHR = 1/2*(iHR+iFA);
elseif iFA == 1
    iFA = 1 - (1/2*(iHR+iFA));
elseif iFA == 0
    iFA = 1/2*(iHR+iFA);
end
d = norminv(iHR,0,1) - norminv(iFA,0,1);
end
