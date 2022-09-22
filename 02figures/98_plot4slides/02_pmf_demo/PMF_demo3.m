
clear all; clc; close all;

% plot parameters
cmp =  [216, 49, 91; 175, 213, 128; 88,193,238]./255;
mksz = 10;

% set parameters
min_x = -1000;
max_x = 1000;
x_axis = [min_x:max_x];
soas = [-500:100:500];
sig = 200;
mu = 0;
c = 150;
lc = mu-c;
rc = mu+c;

%% physical space

figure;
set(gca, 'LineWidth', 2, 'FontSize', 30); hold on

plot(x_axis, zeros(1,numel(x_axis)),'Color',[1, 1, 1])
% xline(-700, 'LineWidth',2)
% xline(-200,':','LineWidth',2)

xticks([-700 -200])
xticklabels({'lag','0'})
set(gca,'TickDir','out');

ax = gca;
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.Visible = 'off';

saveas(gca, 'pmf1', 'epsc')

%% internal space

figure;
set(gca, 'LineWidth', 2, 'FontSize', 30); hold on

plot(x_axis, zeros(1,numel(x_axis)),'Color',[1, 1, 1])
% xline(-500, 'LineWidth',2)
% xline(0,':','LineWidth',2)

xticks([-500 0])
xticklabels({'lag','PSS'})
set(gca,'TickDir','out');

ax = gca;
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.Visible = 'off';

saveas(gca, 'pmf2', 'epsc')

%% add noise, plot the distribution of SOA

figure;
set(gca, 'LineWidth', 2, 'FontSize', 30); hold on

soa = -500;
gauss = normpdf(x_axis, soa, sig);

plot(x_axis, gauss, 'k','LineWidth',3)

% look better
ylim([0 max(gauss)*1.5])
% tick options
xticks([soa 0])
xticklabels({'lag','PSS'})
set(gca,'TickDir','out');
ax = gca;
ax.XAxis.LineWidth = 2;
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.Visible = 'off';

saveas(gca, 'pmf3', 'epsc')

%% criterion: no color

figure;
set(gca, 'LineWidth', 2, 'FontSize', 30); hold on

% dist
soa = -500;
gauss = normpdf(x_axis, soa, sig);

% lines
xline(0,':','LineWidth',2)
xline(lc,'LineWidth',2)
xline(rc,'LineWidth',2)

plot(x_axis, gauss, 'k','LineWidth',3)
% look better
ylim([0 max(gauss)*1.5])
% tick options
xticks([soa 0])
xticklabels({'lag','PSS'})
set(gca,'TickDir','out');
ax = gca;
ax.XAxis.LineWidth = 2;
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.Visible = 'off';

saveas(gca, 'pmf4', 'epsc')


%% criterion: a-lead

figure;
set(gca, 'LineWidth', 2, 'FontSize', 30); hold on

% dist
soa = -500;
gauss = normpdf(x_axis, soa, sig);

% lines
xline(0,':','LineWidth',2)
xline(lc,'LineWidth',2)
xline(rc,'LineWidth',2)

% add colors
yl = ylim;
max_y = yl(2);

% % areas left to -c
area1 = patch([min_x:lc, fliplr(min_x:lc)], ...
    [zeros(1, numel(min_x:lc)), repmat(max_y, 1, numel(min_x:lc))],cmp(3,:), ...
    'FaceAlpha',.5,'EdgeColor','none');

plot(x_axis, gauss, 'k','LineWidth',3)
% look better
ylim([0 max(gauss)*1.5])
% tick options
xticks([soa 0])
xticklabels({'lag','PSS'})
set(gca,'TickDir','out');
ax = gca;
ax.XAxis.LineWidth = 2;
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.Visible = 'off';

saveas(gca, 'pmf5', 'epsc')

%% criterion: simul

figure;
set(gca, 'LineWidth', 2, 'FontSize', 30); hold on

% dist
soa = -500;
gauss = normpdf(x_axis, soa, sig);

% lines
xline(0,':','LineWidth',2)
xline(lc,'LineWidth',2)
xline(rc,'LineWidth',2)

% add colors
yl = ylim;
max_y = yl(2);

% plot area between -c and c
area2 = patch([lc:rc, fliplr(lc:rc)], ...
    [zeros(1, numel(lc:rc)), repmat(max_y, 1, numel(lc:rc))],cmp(2,:), ...
    'FaceAlpha',.5,'EdgeColor','none');

plot(x_axis, gauss, 'k','LineWidth',3)
% look better
ylim([0 max(gauss)*1.5])
% tick options
xticks([soa 0])
xticklabels({'lag','PSS'})
set(gca,'TickDir','out');
ax = gca;
ax.XAxis.LineWidth = 2;
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.Visible = 'off';

saveas(gca, 'pmf6', 'epsc')

%% criterion: v-lead

figure;
set(gca, 'LineWidth', 2, 'FontSize', 30); hold on

% dist
soa = -500;
gauss = normpdf(x_axis, soa, sig);

% lines
xline(0,':','LineWidth',2)
xline(lc,'LineWidth',2)
xline(rc,'LineWidth',2)

% add colors
yl = ylim;
max_y = yl(2);

% plot area right to c
area2 = patch([rc:max_x, fliplr(rc:max_x)], ...
    [zeros(1, numel(rc:max_x)), repmat(max_y, 1, numel(rc:max_x))],cmp(1,:), ...
    'FaceAlpha',.5,'EdgeColor','none');

plot(x_axis, gauss, 'k','LineWidth',3)
% look better
ylim([0 max(gauss)*1.5])
% tick options
xticks([soa 0])
xticklabels({'lag','PSS'})
set(gca,'TickDir','out');
ax = gca;
ax.XAxis.LineWidth = 2;
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.Visible = 'off';

saveas(gca, 'pmf7', 'epsc')