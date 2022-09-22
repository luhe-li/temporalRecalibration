%% TOJ fitting v5

% 2020/04

% This script fits **cumulative Gaussian PMF** for the ternary TOJ data.

clear all; close all; clc;

%% individual input

subjID = input('Please enter participant ID#: ') ;
sess = input('Please enter participant session#: ') ;
    

pre_filenm = ['pretest_sub' num2str(subjID) '_session' num2str(sess) '.mat'];
post_filenm = ['posttest_sub' num2str(subjID) '_session' num2str(sess) '.mat'];

posttest_fitting_ind(pre_filenm, post_filenm, true, true)
