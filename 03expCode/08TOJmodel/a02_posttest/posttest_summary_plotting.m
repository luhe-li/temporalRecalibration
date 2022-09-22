% This script plots mu shift and gain across adaptor SOA for each subject
% based on data_summary.mat

clear all; clc; close all;

subjID = input('Please enter participant ID#: ') ;
load('data_summary.mat')

%% plot
adaptor = [adaptor{subjID,:}];
delta_mu = [delta_mu{subjID,:}];
lb_68 = [delta_mu_lb68{subjID,:}];
ub_68 = [delta_mu_ub68{subjID,:}];

% plotting parameters
cMAP1 = [158,202,225; 161,217,155; 252,174,145]./255;
cMAP2 = [33,113,181; 35,139,69; 203,24,29]./255;
alpha = 0.6; % transparency
cMAPgray = [0, 0, 0; alpha.*ones(1,3)];

% open figure
f1 = figure; hold on;
set(gca, 'FontSize', 20, 'LineWidth', 1.5)
box off
colororder(cMAPgray);

% sort sessions
[B I] = sort(adaptor);
sort_adaptor = adaptor(I);
sort_delta_mu = delta_mu(I);
sort_lb68 = lb_68(I);
sort_ub68 = ub_68(I);

% plot delta_mu in ms
yyaxis left
errorbar(sort_adaptor, sort_delta_mu, ...
    sort_delta_mu - sort_lb68, sort_ub68 - sort_delta_mu, ...
    'o','LineWidth', 1.5)
ylim([-200 200])
yline(0)
ylabel('\mu_{post} - \mu_{pre}')

% plot sorted gain excluding SOA=0
non_zero_idx =  sort_adaptor ~= 0;
non_zero_adaptor = sort_adaptor(non_zero_idx);
gain = sort_delta_mu(non_zero_idx)./non_zero_adaptor;
gain_lb68 = sort_lb68(non_zero_idx)./non_zero_adaptor;
gain_ub68 = sort_ub68(non_zero_idx)./non_zero_adaptor;

yyaxis right
errorbar(non_zero_adaptor - 5, gain, ...
    gain - gain_lb68, gain_ub68 - gain, ...
    'o','LineWidth', 1.5)
ylim([-1.5 1.5])
ylabel('gain')

% look better
adaptors = [-700, -300, -200, -100,  0,  100, 200, 300, 700];
xticks(adaptors)
xlim([min(adaptors)-50, max(adaptors)+50])
xlabel('adaptor SOA')
title(['sub' num2str(subjID)])
fignm = ['summary_sub', num2str(subjID) ];
saveas(gca,fignm,'epsc')
