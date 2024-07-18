% A cute exponential TOJ demo with gif output

clear all; clc; close all;


%% path
% 
% cur_dir                            = pwd;
% [project_dir, ~]                   = fileparts(cur_dir);
% out_dir = fullfile(cur_dir, 'fig6_exponential_demo');
% addpath(genpath(fullfile(project_dir, 'data'))); 
% addpath(genpath(fullfile(project_dir, 'recalibrationModel'))); 
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% set free param

% set free parameters for measurement distribution
% tau = 20;
sigma_a = 60;
sigma_v = 90;

tau_a         = 80;
tau_v         = 20;
tau = tau_a - tau_v;
soa = 0;

%% set up

lw = 3;
fontSz = 7;
titleFontSz = 10;
mksz = 3;
cmp =  [216, 49, 91; 175, 213, 128; 88,193,238]./255;

%% subplot1: measurement pdf

figure; hold on
set(gcf, 'Position',[0 0 120 200])
set(gcf, 'color', 'white')

subplot 211;
left = 0.2;      % Left position of the subplot
bottom = 0.65;    % Bottom position of the subplot
width = 0.7;     % Width of the subplot
height = 0.25;    % Height of the subplot
subplot('Position', [left, bottom, width, height]);
set(gca, 'LineWidth', lw, 'FontSize', fontSz); hold on
set(gca,'TickDir','out');

x_axis        = -0:200; % ms

g_v           = 1/sigma_v .* exp(-(1/sigma_v).*(x_axis - (soa + tau_v)));
g_a           = 1/sigma_a .* exp(-(1/sigma_a).*(x_axis - (soa + tau_a)));

% g_v = g_v./sum(g_v);
% g_a = g_a./sum(g_a);

idx_tau_v     = find(x_axis == soa + tau_v);
idx_tau_a     = find(x_axis == soa + tau_a);

v = plot(x_axis(idx_tau_v:end), g_v(idx_tau_v:end),'Color',cmp(1,:),'LineWidth',lw);
plot([idx_tau_v,idx_tau_v], [0, g_v(idx_tau_v)],'Color',cmp(1,:),'LineWidth',lw)

a = plot(x_axis(idx_tau_a:end), g_a(idx_tau_a:end),'Color',cmp(3,:),'LineWidth',lw);
plot([idx_tau_a,idx_tau_a], [0, g_a(idx_tau_a)],'Color',cmp(3,:),'LineWidth',lw);

lgd = legend([v,a],{'v', 'a'},'Location','best');
legend boxoff 
lgd.ItemTokenSize = [10,10];

xlim([x_axis(1), x_axis(end)])
xlabel('Arrival latency (s)')
xticks([0:100:200])
xticklabels({'0','0.1','0.2'})

% Make the Y-axis and Y-axis label invisible
ax = gca;
ax.YAxis.Visible = 'off';
ylabel(ax, 'Probability');

%% subplot2: likelihood 

subplot 212;
left = 0.2;      % Left position of the subplot
bottom = 0.2;    % Bottom position of the subplot
width = 0.7;     % Width of the subplot
height = 0.25;    % Height of the subplot
subplot('Position', [left, bottom, width, height]);
set(gca, 'LineWidth', lw, 'FontSize', fontSz); hold on
set(gca,'TickDir','out');

% pre-recalibration

soa = -400;
min_x = -600;
max_x = 200;
x_axis = min_x:max_x;

idx_peak = find(x_axis == (soa + tau));
lx_axis = x_axis(1:idx_peak);
rx_axis = x_axis(idx_peak+1:end);

lf = (1/(sigma_a + sigma_v)).* exp(1/sigma_v * (lx_axis - ((soa + tau))));
rf = (1/(sigma_a + sigma_v)).* exp(-1/sigma_a * (rx_axis - ((soa + tau))));
likelihood = [lf, rf];

plot(lx_axis, lf,'--','Color',cmp(1,:),'LineWidth',lw);
plot(rx_axis, rf,'--','Color',cmp(3,:),'LineWidth',lw);

% post-recalibration

soa = -200;
min_x = -600;
max_x = 200;
x_axis = min_x:max_x;

idx_peak = find(x_axis == (soa + tau));
lx_axis = x_axis(1:idx_peak);
rx_axis = x_axis(idx_peak+1:end);

lf = (1/(sigma_a + sigma_v)).* exp(1/sigma_v * (lx_axis - ((soa + tau))));
rf = (1/(sigma_a + sigma_v)).* exp(-1/sigma_a * (rx_axis - ((soa + tau))));
likelihood = [lf, rf];

plot(lx_axis, lf,'Color',cmp(1,:),'LineWidth',lw);
plot(rx_axis, rf,'Color',cmp(3,:),'LineWidth',lw);

% look better
xlabel('SOA (s)')
xticks([-500, 0, 500])
xticklabels({'-0.5','0','0.5'})

ax = gca;
ax.YAxis.Visible = 'off';
ylabel(ax, 'Likelihood');



flnm = 'exponential_demo';
saveas(gca,fullfile(out_dir,flnm),'pdf')