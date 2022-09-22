% a simple gaussian for slides
clear all; clc; close all;

% plot parameters
cmp =  [216, 49, 91; 175, 213, 128; 88,193,238]./255;
mksz = 10;

% set parameters
min_x = -600;
max_x = 600;
x_axis = [min_x:max_x];
sig = 200;
mu = 0;
soa = -500;
gauss = normpdf(x_axis, soa, sig);

%% PSS shift
figure;
set(gca, 'LineWidth', 2, 'FontSize', 30); hold on

soa = 0;
gauss = normpdf(x_axis, soa, sig);

plot(x_axis, gauss, 'Color',cmp(2,:),'LineWidth',12)

% look better
ylim([0 max(gauss)*1.5])

% tick options
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

saveas(gca,'wide_gaussian','png')

%% sigma shift: large sigma
figure;
set(gca, 'LineWidth', 2, 'FontSize', 30); hold on

gauss = normcdf(x_axis, soa, sig);

plot(x_axis, gauss, 'Color',cmp(1,:),'LineWidth',12)

% look better
ylim([0 max(gauss)*1.5])

% tick options
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
saveas(gca,'slope_demo1','png')

%% sigma shift: small sigma
figure;
set(gca, 'LineWidth', 2, 'FontSize', 30); hold on

sig = 50;
gauss = normcdf(x_axis, soa, sig);

plot(x_axis, gauss, 'Color',cmp(1,:),'LineWidth',12)

% look better
ylim([0 max(gauss)*1.5])

% tick options
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
saveas(gca,'slope_demo2','png')



%% criterion shift: narrow gaussian
figure;
set(gca, 'LineWidth', 2, 'FontSize', 30); hold on

soa = 0;
sig = 60;
gauss = normpdf(x_axis, soa, sig);

plot(x_axis, gauss, 'Color',cmp(2,:),'LineWidth',12)

% look better
ylim([0 max(gauss)*1.5])

% tick options
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

saveas(gca,'narrow','png')