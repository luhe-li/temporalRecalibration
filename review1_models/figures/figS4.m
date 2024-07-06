% fig S4: outlier behavior

clear; clc; close all;

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% set up

sub_slc = 5;
adp_slc = 1;
save_fig = 1;

athe_path = fullfile(projectDir, 'atheoretical_models_VBMC','exp_shiftMu');
files = dir(fullfile(athe_path, 'sub-*'));

ss = sub_slc;
D = load(fullfile(athe_path, files(ss).name));

%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure set up

cmp1 = [229, 158, 168; 203, 227, 172; 171,223,235;]./255;
cmp2 = [216, 49, 91; 175, 213, 128; 88,193,238]./255;

lw = 0.5;
fontSZ = 7;
titleSZ = 9;
dotSZ = 15;
adaptor_soa = D.pred.adaptor_soa;
test_soa = D.data.pre_ms_unique;

figure;
set(gcf, 'Position',[0,0,420,130]);

subplot(1,2,1); 
set(gca, 'Position', [0.1, 0.2, 0.45, 0.8]);
set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
set(gca, 'FontName', 'Helvetica Neue');
set(gca, 'FontWeight', 'Light');
hold on

% pre
set(gca, 'ColorOrder', cmp1)
scatter(test_soa, D.data.pre_pResp, dotSZ,'LineWidth',lw);

% post
set(gca, 'ColorOrder', cmp2)
scatter(test_soa, D.data.post_pResp, dotSZ, 'filled');

% adaptor
xline(adaptor_soa(adp_slc))

% look better
xlim([-550, 550])
x_ticks = [-500, -300:100:300, 500];
y_ticks = [0:0.25:1];
xticks(x_ticks)
xticklabels(strsplit(num2str(x_ticks./1000)))
yticks(y_ticks)
xtickangle(45)
yticklabels(strsplit(num2str(y_ticks)))
xlabel('Test test_soa (s)', 'FontName', 'Helvetica Neue', 'FontWeight', 'Light');
ylabel('Probability', 'FontName', 'Helvetica Neue', 'FontWeight', 'Light');
title(sprintf('S%i, Adaptor test_soa = -0.7 s', sub_slc))
