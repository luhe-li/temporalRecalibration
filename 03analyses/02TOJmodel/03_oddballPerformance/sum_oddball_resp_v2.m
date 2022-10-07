%% sum oddball response v2

% This script summarize hit rate and false alarm in exposure phase and
% post-test phase, also provides HR and FAR across individual blocks

% This version uses function **summarizeBlockPerformance.m**

clear all; close all;

%% input

subjID = input('Please enter participant ID#: ') ;
sess = input('Please enter participant session#: ') ;

%% Exposure phase
load(['exposure_sub' num2str(subjID) '_session' num2str(sess) '.mat']);
[hitRate1, FA1]                               = summarizeBlockPerformance(1, ExpInfo.nTTrials, ExpInfo.idxOddballA, ExpInfo.idxOddballV, Response.oddball);

% print
fprintf('Exposure: auditory oddball hit rate  = %4.2f, visual oddball hit rate =  %4.2f\n', hitRate1)
fprintf('Exposure: auditory oddball FA rate   = %4.2f, visual oddball FA rate =  %4.2f\n', FA1)
fprintf('Visual standard = %4.2f, comparison = %4.2f\n', unique(ExpInfo.VIntensity))
fprintf('Auditory standard = %4.2f, comparison = %4.2f\n', unique(ExpInfo.AIntensity))

%% Post-test phase
load(['posttest_sub' num2str(subjID) '_session' num2str(sess) '.mat']);
[hitRate2, FA2]                               = summarizeBlockPerformance(1, ExpInfo.expoNTTrials, ExpInfo.idxOddballA, ExpInfo.idxOddballV, ExpoResponse.oddball);

% print
fprintf('Post-test: auditory oddball hit rate = %4.2f, visual oddball hit rate =  %4.2f\n', hitRate2)
fprintf('Post-test: auditory oddball FA rate  = %4.2f, visual oddball FA rate =  %4.2f\n', FA2)
fprintf('Visual standard = %4.2f, comparison = %4.2f\n', unique(ExpInfo.VIntensity))
fprintf('Auditory standard = %4.2f, comparison = %4.2f\n', unique(ExpInfo.AIntensity))

%% individual block post-test phase
n = ExpInfo.nTrialsinBlock * 4; % number of exposure trials in each block
endtrial = 0:n:ExpInfo.expoNTTrials;
endtrial = [endtrial(2:end), ExpInfo.expoNTTrials];
starttrial = endtrial - n + 1;

for i = 1:length(endtrial)
    [hitRate(i,:), FA(i,:)] = summarizeBlockPerformance(starttrial(i), endtrial(i), ExpInfo.idxOddballA, ExpInfo.idxOddballV, ExpoResponse.oddball);
end

figure; 
set(gca,'FontSize',15,'linewidth',2)
set(gcf,'Position', [10, 10, 800, 300]);

f1 = subplot(1,2,1);hold on
set(f1, 'FontSize',15,'linewidth',2)
nblock = 0:length(endtrial); 
plot(nblock, [hitRate1; hitRate]','-o','linewidth',2); 
xlim([0, 7])
ylim([0 1])
xlabel('block')
ylabel('hit rate')
xticks(nblock)
xticklabels({'exposure','post1','post2','post3','post4','post5','post6','post7'})


f2 = subplot(1,2,2);
set(f2, 'FontSize',15,'linewidth',2); hold on
plot(nblock, [FA1; FA]','-o','linewidth',2)
xlim([0, 7])
ylim([0 1])
xlabel('block')
ylabel('FA rate')
legend({'auditory','visual'})
xticks(nblock)
xticklabels({'exposure','post1','post2','post3','post4','post5','post6','post7'})

sgtitle('Oddball task performance')

fignm = ['oddball_sub', num2str(subjID) '_session' num2str(sess)];
saveas(gca,fignm,'epsc')