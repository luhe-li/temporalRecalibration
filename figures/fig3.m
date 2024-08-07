% fig 3: illustration of recalibration models
% A. distribution of auditory and visual arrival latencty
% B. Recalibration of measurement distribution
% C. recalibration models

% out_dir = fullfile(cur_dir, 'fig6_exponential_demo');

clear; close all;

out_dir = fullfile(pwd, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% set free param

sigma_a = 60;
sigma_v = 90;

tau_a = 80;
tau_v = 20;
tau = tau_a - tau_v;
soa = 0;

lw = 1;
fontSZ = 7;
titleSZ = 9;
dotSZ = 10;
cmp =  [216, 49, 91; 175, 213, 128; 88,193,238]./255;

%% A. distribution of auditory and visual arrival latencty

figure('Position', [0, 0, 420, 60])

subplot 121;
set(gca, 'LineWidth', 0.5, 'FontSize', fontSZ, 'TickDir', 'out', 'FontName', 'Helvetica');
hold on

x_axis        = 0:200; % ms

g_v           = 1/sigma_v .* exp(-(1/sigma_v).*(x_axis - (soa + tau_v)));
g_a           = 1/sigma_a .* exp(-(1/sigma_a).*(x_axis - (soa + tau_a)));

idx_tau_v     = find(x_axis == soa + tau_v);
idx_tau_a     = find(x_axis == soa + tau_a);

v = plot(x_axis(idx_tau_v:end), g_v(idx_tau_v:end),'Color',cmp(1,:),'LineWidth',lw);
plot([idx_tau_v,idx_tau_v], [0, g_v(idx_tau_v)],'Color',cmp(1,:),'LineWidth',lw)

a = plot(x_axis(idx_tau_a:end), g_a(idx_tau_a:end),'Color',cmp(3,:),'LineWidth',lw);
plot([idx_tau_a,idx_tau_a], [0, g_a(idx_tau_a)],'Color',cmp(3,:),'LineWidth',lw);

lgd = legend([v,a],{'Visual signal', 'Auditory Signal'},'Location','best');
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

%% B. recalibration of likelihood

subplot 122;
set(gca, 'LineWidth', 0.5, 'FontSize', fontSZ, 'TickDir', 'out', 'FontName', 'Helvetica');
hold on

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
xlabel('Measurement of SOA (s)')
xticks([-500, 0, 500])
xticklabels({'-0.5','0','0.5'})

ax = gca;
ax.YAxis.Visible = 'off';
ylabel(ax, 'Probability');

flnm = 'AB_exponential';
saveas(gca,fullfile(out_dir,flnm),'pdf')

%% C. Illustration of models

%% simulation

tau = 0;
sigma_a = 60;
sigma_v = 80;
criterion = 77.23;
lambda = 0.018;
p_common = 0.5;
alpha = 0.0052;
sigma_C1 = 51.9;
sigma_C2 = 261.39;

lw = 0.5;
adapter = [-0.7:0.01:0.7];
xtks = [-0.7, 0, 0.7];

% heuristic model
p = measurementGiven0(adapter*1e3, tau, sigma_a, sigma_v);

% causal inference model
model.num_ses = 9;
model.thres_R2 = 0.95;
model.expo_num_sim = 1e2; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_bin  = 100; % numer of bin to approximate tau_shift distribution
model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.4*1e3; % in second, where measurements are likely to reside
model.num_sample = 1e3; % number of samples for simulating psychometric function with causal inference, only used in pmf_exp_CI
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = adapter*1e3;
model.toj_axis_finer = 0; % simulate pmf with finer axis
model.adaptor_axis_finer = 0; % simulate with more adpators
model_str = 'cauInf_asym';
[projectDir, ~]= fileparts(pwd);
addpath(genpath(fullfile(projectDir, 'recalibration_models_VBMC',model_str)))
currModel = str2func(['nll_' model_str]);
model.mode       = 'predict';
pred =  currModel([tau, sigma_a, sigma_v, criterion, lambda, p_common, alpha, sigma_C1, sigma_C2], model, []);
recal = mean(pred.pss_shift, 2);

%% plot
figure('Position', [0, 0, 420, 90]);

subplot 131
set(gca, 'LineWidth', 0.5, 'FontSize', fontSZ, 'TickDir', 'out', 'FontName', 'Helvetica'); hold on
plot(adapter, 0.3*adapter,'k','LineWidth',lw)
xticks(xtks)
xlim([-0.7, 0.7])
ylim([-0.7, 0.7])
yline(0,'--')
yticks(xtks)
xlabel('Adapter SOA (s)')

subplot 132
set(gca, 'LineWidth', 0.5, 'FontSize', fontSZ, 'TickDir', 'out', 'FontName', 'Helvetica'); hold on
plot(adapter, adapter.* p * 3.4,'k','LineWidth',lw);
yline(0,'--')
ylim([-0.001, 0.001])
xlim([-0.7, 0.7])
yticks([])
xticks(xtks)
xlabel('Adapter SOA (s)')

subplot 133
set(gca, 'LineWidth', 0.5, 'FontSize', fontSZ, 'TickDir', 'out', 'FontName', 'Helvetica'); hold on
plot(adapter, recal', 'k','LineWidth',lw);
yline(0,'--')
ylim([-100, 100])
xlim([-0.7, 0.7])
yticks([])
xticks(xtks)
xlabel('Adapter SOA (s)')

% save
flnm = 'C_model';
saveas(gca,fullfile(out_dir,flnm),'pdf')

%% utility function

function p = measurementGiven0(soa_m, tau, sigma_a, sigma_v)

p = zeros(size(soa_m));
% right side of measurement distribution
bool_m = soa_m > tau;
p(bool_m) = (1/(sigma_a + sigma_v)).* exp(-1/sigma_v .* (soa_m(bool_m) - tau));
% left side of measurement distribution
p(~bool_m) = (1/(sigma_a + sigma_v)).* exp(1/sigma_a .* (soa_m(~bool_m) - tau));

end