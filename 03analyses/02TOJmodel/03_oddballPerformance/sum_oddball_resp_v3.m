%% sum oddball response v3

% This script summarize hit rate and false alarm in exposure phase and
% post-test phase, also calculates d' across individual blocks

% This version uses function **summarizeBlockPerformance.m**

clear all; close all;

%% input

subjID = input('Please enter participant ID#: ') ;
sess = input('Please enter participant session#: ') ;

%% Exposure phase
load(['exposure_sub' num2str(subjID) '_session' num2str(sess) '.mat']);
[expo_hitRate, expo_FA, expo_btst_d, expo_lb_68CI, expo_ub_68CI] = summarizeBlockPerformance2(1, ExpInfo.nTTrials, ExpInfo.idxOddballA, ExpInfo.idxOddballV, Response.oddball);

% print
fprintf('Exposure: auditory oddball hit rate  = %4.2f, visual oddball hit rate =  %4.2f\n', expo_hitRate)
fprintf('Exposure: auditory oddball FA rate   = %4.2f, visual oddball FA rate =  %4.2f\n', expo_FA)
fprintf('Visual standard = %4.2f, comparison = %4.2f\n', unique(ExpInfo.VIntensity))
fprintf('Auditory standard = %4.2f, comparison = %4.2f\n', unique(ExpInfo.AIntensity))

%% Post-test phase
load(['posttest_sub' num2str(subjID) '_session' num2str(sess) '.mat']);
[post_hitRate, post_FA, post_btst_d, post_lb_68CI, post_ub_68CI] = summarizeBlockPerformance2(1, ExpInfo.expoNTTrials, ExpInfo.idxOddballA, ExpInfo.idxOddballV, ExpoResponse.oddball);

% print
fprintf('Post-test: auditory oddball hit rate = %4.2f, visual oddball hit rate =  %4.2f\n', post_hitRate)
fprintf('Post-test: auditory oddball FA rate  = %4.2f, visual oddball FA rate =  %4.2f\n', post_FA)
fprintf('Visual standard = %4.2f, comparison = %4.2f\n', unique(ExpInfo.VIntensity))
fprintf('Auditory standard = %4.2f, comparison = %4.2f\n', unique(ExpInfo.AIntensity))

%% individual block post-test phase
n = ExpInfo.nTrialsinBlock * ExpInfo.nRep; % number of exposure trials in each block
endtrial = 0:n:ExpInfo.expoNTTrials;
if rem(ExpInfo.nTTrials, ExpInfo.nTrialsinBlock) ~= 0
    endtrial = [endtrial(2:end), ExpInfo.expoNTTrials];
else
    endtrial = endtrial(2:end);
end
starttrial = endtrial - n + 1;

for i = 1:length(endtrial)
    [hitRate(i,:), FA(i,:), btst_d{i}, lb_68CI(i,:), ub_68CI(i,:)] = ...
        summarizeBlockPerformance2(starttrial(i), endtrial(i), ...
        ExpInfo.idxOddballA, ExpInfo.idxOddballV, ExpoResponse.oddball);
end

% calculate dprime for each session
all_HR = [expo_hitRate; hitRate];
all_FA = [expo_FA; FA];
all_lb = [expo_lb_68CI; lb_68CI]';
all_ub = [expo_ub_68CI; ub_68CI]';

for s = 1:size(all_HR, 1)
    for m = 1:2
        iHR = all_HR(s,m);
        iFA = all_FA(s,m);
        if iHR == 1
            iHR = 1 - (1/2*(iHR+iFA));
        elseif iHR == 0
            iHR = 1/2*(iHR+iFA);
        elseif iFA == 1
            iFA = 1 - (1/2*(iHR+iFA));
        elseif iFA == 0
            iFA = 1/2*(iHR+iFA);
        end
        d(s,m) = norminv(iHR,0,1) - norminv(iFA,0,1);
    end
end

clt = [189,0,38;65,182,196]./255;
figure; hold on
set(gca,'FontSize',15,'linewidth',2)
set(gcf,'Position', [0, 0, 600, 400]);
nblock = 1:size(all_HR, 1); 
% plot(nblock, d(:,1)','-o','Color',clt(1,:),'linewidth',2); 
% plot(nblock, d(:,2)','-o','Color',clt(2,:),'linewidth',2); 
errorbar(nblock, d(:,1)',d(:,1)' - all_lb(1,:), all_ub(1,:) - d(:,1)',...
    '-o','Color',clt(1,:),'linewidth',2)
errorbar(nblock-0.1, d(:,2)',d(:,2)' - all_lb(2,:), all_ub(2,:) - d(:,2)',...
    '-o','Color',clt(2,:),'linewidth',2)
xlim([0, length(nblock)+1])
ylim([0 5])
xlabel('exposure + posttest blocks')
ylabel('d''')
xticks(nblock)
title(['sub' num2str(subjID)])
% xticklabels({'exposure','post1','post2','post3','post4','post5'})
legend({'A','V'})

fignm = ['dprime_sub', num2str(subjID) '_session' num2str(sess)];
saveas(gca,fignm,'epsc')

% save dprime averaged across blocks
ad = mean(d(:,1));
vd = mean(d(:,2));

load("all_dprime.mat")
adprime(subjID, sess) = ad;
vdprime(subjID, sess) = vd;
save('all_dprime','adprime','vdprime')