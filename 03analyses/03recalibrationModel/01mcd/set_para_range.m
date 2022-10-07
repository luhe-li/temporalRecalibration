% This script takes the fitted parameters from Parise et al., 2021 to find
% out the best range for tau parameters. The lower and upper bound was set
% by 95% CI given parameter samples

clear all; close all; clc; rng('shuffle');

%% extract fitted parameters

% keep the order: a, v, av
fit_ta = [84, 80, 58, 110, 92, 73, 50, 53, 52, 52, 54, 64, 74];
fit_tv = [130, 137, 120, 152, 142, 141, 122, 123, 114, 109, 103, 119, 122];
fit_tav = [519, 490, 562, 1147, 882, 1087, 748, 372, 764, 779, 799, 779, 770];

fit_para = [fit_ta; fit_tv; fit_tav];

%% calculate bounds for 95%CI
m = mean(fit_para, 2);
se = std(fit_para, [], 2)./sqrt(size(fit_para, 2));
ub = m + 1.96.*se;
lb = m - 1.96.*se;

% check where the bounds are
figure; hold on; histogram(fit_ta,10);xline(lb(1),'r'); xline(ub(1),'r')

% doesn't look right... maybe we should use min and max +- se

%% use min and max +- se

min_p = min(fit_para, [], 2);
max_p = max(fit_para, [], 2);
ub2 = (max_p + se)./1e3; % s
lb2 = (min_p - se)./1e3; % s
