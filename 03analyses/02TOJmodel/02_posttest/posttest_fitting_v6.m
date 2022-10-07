%% TOJ fitting v6

% 2020/05

% This script fits **cumulative Gaussian PMF** for the ternary TOJ data.

clear all; close all; clc;

%% loop through all the files

folder1 = '/Volumes/GoogleDrive-107780329657346345479/My Drive/2021A-22/01project/04experimentalCode/03ExpCpde/01pretest/data';
filePattern1 = fullfile(folder1, 'pretest_sub*.mat');
files1 = dir(filePattern1);

folder2 = '/Volumes/GoogleDrive-107780329657346345479/My Drive/2021A-22/01project/04experimentalCode/03ExpCpde/04posttest/data';
filePattern2 = fullfile(folder2, 'posttest_sub*.mat');
files2 = dir(filePattern2);

for dd = 1:size(files1,1)
    disp(dd)

    pre_filenm = files1(dd).name;
    post_filenm = files2(dd).name;

    posttest_fitting_ind(pre_filenm, post_filenm, true, true)

end
