% a05_pmf_behavioral_results
% This script fits the psychometric function (full model with 8 parameters) to the pre and post data.
% It also returns bootsrapped parameters for each session

%% fitting
% Fit raw data with psychometric function
clear all; close all; clc; rng(1);

% add data folder to the path
currentDir = pwd;
exptDir = currentDir(1:regexp(pwd,'03expCode')-1);
folder1 = [exptDir '03ExpCode/01pretest/data'];
folder2 = [exptDir '03ExpCode/04posttest/data'];

addpath(genpath([exptDir '03ExpCode/05functions']));
addpath(genpath([exptDir '03ExpCode/06helperFunctions']));
addpath(genpath('function'))
addpath(genpath(folder1));
addpath(genpath(folder2));

filePattern1 = fullfile(folder1, 'pretest_sub*.mat');
files1 = dir(filePattern1);
filePattern2 = fullfile(folder2, 'posttest_sub*.mat');
files2 = dir(filePattern2);

% define subjects and sessions to use
all_sub = 3:10;
all_sess = 1:9;

% % a) initiate
% [best_para, estP_btst]= deal(cell(all_sub(end), all_sess(end)));

% b) initiate: to use bootstrapping, use
[best_para, estP_btst, btst_nll, lb_68CI, ub_68CI]= deal(cell(all_sub(end), all_sess(end)));

[fit_nll, idx_adaptor_soa] = deal(NaN(all_sub(end), all_sess(end)));

% fit each session and each subject individually
for dd = 1:size(files1,1)

    pre_filenm = files1(dd).name;
    post_filenm = files2(dd).name;

    % extract subj and sess

    strSplit    = split(pre_filenm, ["_", "."]);
    % subject_id: i
    [st,ed] = regexp(strSplit{2}, 'sub');
    i = str2num(strSplit{2}(ed+1:end));
    % session: j
    [st,ed] = regexp(strSplit{3}, 'session');
    j = str2num(strSplit{3}(ed+1:end));

%     % a) fit and return best-fitting parameters
%     [best_para{i,j}, fit_nll(i,j), idx_adaptor_soa(i,j)] = pmf_fit_sub_sess(pre_filenm, post_filenm);

    % b) to use bootstrapping, use:
    [best_para{i,j}, idx_adaptor_soa(i,j), fit_nll(i,j),...
        estP_btst{i,j}, btst_nll{i,j}, lb_68CI{i,j}, ub_68CI{i,j}] = pmf_fitting_v0(pre_filenm, post_filenm);

end

save('best_para_2')
