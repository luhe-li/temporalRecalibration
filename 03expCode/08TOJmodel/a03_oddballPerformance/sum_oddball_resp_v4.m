%% sum oddball response v3

% This script summarize dprime averaged between exposure phase and
% post-test phase, for each session, for each individual

% It plots auditory and visual dprime by session

% This version uses function **summarizeBlockPerformance2.m**

clear all; clc; close all;

subjIDs = 3:8;
sessions = 1:9;

[expo_d, post_d] = deal(NaN(length(subjIDs), length(sessions), 2));

for i = 1:length(subjIDs)
    subjID = subjIDs(i);
    for j = 1:length(sessions)
        sess = sessions(j);

        load(['exposure_sub' num2str(subjID) '_session' num2str(sess) '.mat']);
        [expo_hitRate, expo_FA, expo_btst_d, expo_lb_68CI, expo_ub_68CI] = summarizeBlockPerformance2(1, ExpInfo.nTTrials, ExpInfo.idxOddballA, ExpInfo.idxOddballV, Response.oddball);

        load(['posttest_sub' num2str(subjID) '_session' num2str(sess) '.mat']);
        [post_hitRate, post_FA, post_btst_d, post_lb_68CI, post_ub_68CI] = summarizeBlockPerformance2(1, ExpInfo.expoNTTrials, ExpInfo.idxOddballA, ExpInfo.idxOddballV, ExpoResponse.oddball);

        for k = 1:2 % 1 = A; 2 = V
            expo_d(i, j, k) = calculate_dprime(expo_hitRate(k), expo_FA(k));
            post_d(i, j, k) = calculate_dprime(post_hitRate(k), post_FA(k));
        end

    end
end

%% plotting

label = {'auditory'; 'visual'};

% exposure phase
for k = 1:2
    figure; hold on
    set(gca, 'FontSize', 20, 'LineWidth', 1.5)
    box off
    for i = 1:length(subjIDs)
    plot(sessions, expo_d(i, :, k),'-o','LineWidth', 1.5,'DisplayName',['sub' num2str(subjIDs(i))])
    end
    xticks(sessions)
    xticklabels(sessions)
    title([label{k} ' oddball'])
    ylabel('d''')
    xlabel('exposure sessions')
    ylim([0 5])
    expo_outlier{k} = isoutlier(expo_d(:, :, k), 'mean');
    legend('Location','bestoutside')
    fignm = ['expo_sub' num2str(subjID) label{k}];
    saveas(gca,fignm,'epsc')
end

% post-test phase
for k = 1:2
    figure; hold on
    set(gca, 'FontSize', 20, 'LineWidth', 1.5)
    box off
    for i = 1:length(subjIDs)
        plot(sessions, post_d(i, :, k),'-o','LineWidth', 1.5,'DisplayName',['sub' num2str(subjIDs(i))])
    end
    xticks(sessions)
    xticklabels(sessions)
    title([label{k} ' oddball'])
    ylabel('d''')
    xlabel('post-test sessions')
    ylim([0 5])
    post_outlier{k} = isoutlier(expo_d(:, :, k), 'mean');
    legend('Location','bestoutside')
    fignm = ['post_sub' num2str(subjID) label{k}];
    saveas(gca,fignm,'epsc')
end


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
