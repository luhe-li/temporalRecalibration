% a05_pmf_behavioral_results
% This script fits the psychometric function (full model with 8 parameters) to the pre and post data.
% It also returns bootsrapped parameters for each session

%% fitting
% Fit raw data with psychometric function
clear all; close all; clc; rng(1);

% add data folder to the path
currentDir = pwd;
exptDir = currentDir(1:regexp(pwd,'03analyses')-1); % project folder
addpath(genpath([exptDir '02data'])); % data folder, necessary to run organize_data
addpath(genpath([exptDir '/03analyses/01analysesFunctions'])); % function folder
addpath(genpath('function')) % function specific to this analysis
folder1 = [exptDir '02data/02pretest'];
folder2 = [exptDir '02data/04posttest'];

filePattern1 = fullfile(folder1, 'pretest_sub*.mat');
files1 = dir(filePattern1);
filePattern2 = fullfile(folder2, 'posttest_sub*.mat');
files2 = dir(filePattern2);

% load('best_para_2.mat')

% define subjects and sessions to use
sub2add = 1;
sess2add = 1;

% fit each session and each subject individually
for i = sub2add

    for j = sess2add

    pre_filenm = ['pretest_sub' num2str(i) '_session' num2str(j) '.mat'];
    post_filenm = ['posttest_sub' num2str(i) '_session' num2str(j) '.mat'];

    %     % a) fit and return best-fitting parameters
    %     [best_para{i,j}, fit_nll(i,j), idx_adaptor_soa(i,j)] = pmf_fit_sub_sess(pre_filenm, post_filenm);

    % b) to use bootstrapping, use:
    [best_para{i,j}, idx_adaptor_soa(i,j), fit_nll(i,j),...
        estP_btst{i,j}, btst_nll{i,j}, lb_68CI{i,j}, ub_68CI{i,j}] = pmf_fitting_v0(pre_filenm, post_filenm);

    end

end

save('best_para_2.mat')
