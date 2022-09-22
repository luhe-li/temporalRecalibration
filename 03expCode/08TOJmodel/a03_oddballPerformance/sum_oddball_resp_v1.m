% This script summarize hit rate and false alarm in exposure phase and
% post-test phase

clear all; close all;

%% input

subjID = input('Please enter participant ID#: ') ;
session = input('Please enter participant session#: ') ;

%% Exposure phase
addpath('/Volumes/GoogleDrive/My Drive/2021A-22/01project/04ExperimentalCode/03ExpCpde/03exposure/data')
load(['exposure_sub' num2str(subjID) '_session' num2str(session) '.mat']);
[hitRate1, FA1]                               = summarizePerformance(ExpInfo.nTTrials, ExpInfo.nTTrials, ExpInfo.idxOddballA, ExpInfo.idxOddballV, Response.oddball);

% print
fprintf('Exposure: auditory oddball hit rate  = %4.2f, visual oddball hit rate =  %4.2f\n', hitRate1)
fprintf('Exposure: auditory oddball FA rate   = %4.2f, visual oddball FA rate =  %4.2f\n', FA1)

%% Post-test phase
load(['posttest_sub' num2str(subjID) '_session' num2str(session) '.mat']);
[hitRate2, FA2]                               = summarizePerformance(ExpInfo.expoNTTrials, ExpInfo.expoNTTrials, ExpInfo.idxOddballA, ExpInfo.idxOddballV, ExpoResponse.oddball);

% print
fprintf('Post-test: auditory oddball hit rate = %4.2f, visual oddball hit rate =  %4.2f\n', hitRate2)
fprintf('Post-test: auditory oddball FA rate  = %4.2f, visual oddball FA rate =  %4.2f\n', FA2)

%% cumulative block-wise post-test phase
% i = ExpInfo.nTrialsinBlock * 4; % number of exposure trials in each block
% endtrial = 0:i:ExpInfo.expoNTTrials;
% endtrial = [endtrial(2:end), ExpInfo.expoNTTrials];
% for i = 1:length(endtrial)
%     [hitRate(i,:), FA(i,:)] = summarizePerformance(endtrial(i), ExpInfo.expoNTTrials, ExpInfo.idxOddballA, ExpInfo.idxOddballV, ExpoResponse.oddball);
% end
% 
% figure; hold on
% nblock =1:length(endtrial); 
% plot(nblock, hitRate','-o')
% ylim([0 1])
% xlabel('1:n block')
% ylabel('cumulative hit rate')
% legend({'auditory','visual'})

