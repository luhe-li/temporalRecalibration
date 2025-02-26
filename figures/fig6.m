% fig 6: demonstration of the normative model of psychometric functions of
% ternary TOJ task

clear; clc; close all;

%% manage paths

restoredefaultpath;
[projectDir, ~]= fileparts(pwd);
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(pwd, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% set free param

% set free parameters for TOJ
tau = 60;
sigma_a = 60;
sigma_v = 90;
criterion = 70;
lambda = 0.02;
soa = 0;

% set free parameters for causal inference
p_common = 0.2;
sigma_C1 = 10;
sigma_C2 = 500;

% calculate PSS
if sigma_a <= sigma_v
mu  = -tau - sigma_v * log((sigma_a + sigma_v)/(2*sigma_v));
else
mu  = -tau + sigma_a * log((sigma_a + sigma_v)/(2*sigma_a));
end

%% set fixed param

% set fixed parameters for TOJ
min_x = -500;
max_x = 500;
x_axis = min_x:max_x;
soas = min_x:50:max_x;

% set fixed param for causal inference
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.bound_full = 10*1e3; % in ms, the bound for prior axis
model.bound_int = 2*1e3; % in ms, where measurements are likely to reside
model.sim_adaptor_soa = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;

% prior
fixP.x_axis     = -model.bound_full:1:model.bound_full;
fixP.x_axis_int = -model.bound_int:1:model.bound_int;
fixP.l_window   = find(fixP.x_axis == fixP.x_axis_int(1));
fixP.r_window   = find(fixP.x_axis == fixP.x_axis_int(end));
fixP.prior_C1   = normpdf(fixP.x_axis_int, 0, sigma_C1);
fixP.prior_C2   = normpdf(fixP.x_axis_int, 0, sigma_C2);

% likelihood that centers around 0
v_peak                       = 0;
idx_peak  = ceil(length(fixP.x_axis)/2);

lf = (1/(sigma_a + sigma_v)).* exp(1/sigma_a .* (fixP.x_axis(1:idx_peak)));
rf = (1/(sigma_a + sigma_v)).* exp(-1/sigma_v .* (fixP.x_axis(idx_peak+1:end)));
fixP.df_likelihood           = [lf, rf] + 1e-40;
fixP.df_likelihood_int       = fixP.df_likelihood(fixP.l_window:fixP.r_window);

fixP.y_criterion             = sigma_v/(sigma_a + sigma_v);

% set up for PMF with causal inference by simulation (pmf_exp_CI)
fixP.num_sample             = 1e4; % number of simulation samples of measurements
model.test_soa              = 0;

% pre-calculate possible likelihoods by shifting default likelihood
fixP.shift_bound = model.bound_int;
shiftRange = -fixP.shift_bound:1:fixP.shift_bound;
numShifts = numel(shiftRange);
indices = (fixP.l_window - shiftRange(:)) + (0:(fixP.r_window - fixP.l_window));
indices = max(min(indices, length(fixP.df_likelihood)), 1);
fixP.likelihoods = fixP.df_likelihood(indices);

% pre-calculate two posteriors
fixP.protopost_C1s = fixP.likelihoods .* repmat(fixP.prior_C1, numShifts, 1);
fixP.protopost_C2s = fixP.likelihoods .* repmat(fixP.prior_C2, numShifts, 1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 0.5;
fontSz = 7;
titleFontSz = 10;
mksz = 3;
cmp = [216, 49, 91; 175, 213, 128; 88,193,238]./255;

%% subplot1: measurement pdf

figure; hold on
set(gcf, 'Position',[10 10 200 150])
set(gcf, 'color', 'white')

subplot 211;
left = 0.2;      % Left position of the subplot
bottom = 0.65;    % Bottom position of the subplot
width = 0.7;     % Width of the subplot
height = 0.25;    % Height of the subplot
subplot('Position', [left, bottom, width, height]);
set(gca, 'LineWidth', lw, 'FontSize', fontSz); hold on
set(gca,'TickDir','out');

idx_peak = find(x_axis == (soa + tau));
lx_axis = x_axis(1:idx_peak);
rx_axis = x_axis(idx_peak+1:end);

lf = (1/(sigma_a + sigma_v)).* exp(1/sigma_v * (lx_axis - ((soa + tau))));
rf = (1/(sigma_a + sigma_v)).* exp(-1/sigma_a * (rx_axis - ((soa + tau))));
measurement_dist = [lf, rf];

l3 = plot(x_axis, measurement_dist, 'k','LineWidth',lw); 
l1 = xline(tau,'LineWidth',lw);
l2 = xline(-mu,':','LineWidth',lw);

% lgd = legend([l1, l2, l3],{'s + \tau','\mu','\Psi_A(s = 0)'},'Location','bestoutside');
% legend boxoff 
% lgd.ItemTokenSize = [10,10];

% look better
ylabel('Probability density')
xlabel('Measurement  of a zero SOA (s)')
% yticks([yl(1), yl(2)])
xticks([-500, 0, 500])
xticklabels({'-0.5','0','0.5'})

%% subplot2: estimate pdf

subplot 212;
left = 0.2;      % Left position of the subplot
bottom = 0.2;    % Bottom position of the subplot
width = 0.7;     % Width of the subplot
height = 0.25;    % Height of the subplot
subplot('Position', [left, bottom, width, height]);
set(gca, 'LineWidth', lw, 'FontSize', fontSz); hold on
set(gca,'TickDir','out');

[p_afirst, p_simul, p_vfirst, shat_sample] = pmf_exp_CI(model.test_soa, fixP,...
    tau, sigma_a, sigma_v, criterion, lambda, p_common);

% Perform kernel smoothing density estimation
[f, xi] = ksdensity(shat_sample,'NumPoints',1000);

indices = xi<=-criterion;
fill([xi(indices), fliplr(xi(indices))], [f(indices), zeros(size(f(indices)))], cmp(3,:), 'EdgeColor', 'none');

indices = xi >= -criterion & xi <= criterion;
fill([xi(indices), fliplr(xi(indices))], [f(indices), zeros(size(f(indices)))], cmp(2,:), 'EdgeColor', 'none');

indices = xi >= criterion;
fill([xi(indices), fliplr(xi(indices))], [f(indices), zeros(size(f(indices)))], cmp(1,:), 'EdgeColor', 'none');

ylabel('Probability density')
xlabel('Estimate of a zero SOA (s)')

xlim([min_x, max_x])
xticks([min_x, 0, max_x])
xticklabels({'-0.5','0','0.5'})

flnm = 'TOJdemo1';
saveas(gca,fullfile(out_dir,flnm),'pdf')

%% subplot3: pmf with causal inference, derived from the estimate distribution

figure; hold on
set(gcf, 'Position',[10 10 230 150])
set(gcf, 'color', 'white')

left = 0.2;      % Left position of the subplot
bottom = 0.2;    % Bottom position of the subplot
width = 0.7;     % Width of the subplot
height = 0.7;    % Height of the subplot
subplot('Position', [left, bottom, width, height]);
set(gca, 'LineWidth', lw, 'FontSize', fontSz); hold on
model.test_soa              = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;

[pre_afirst, pre_simul, pre_vfirst] = pmf_exp_CI(model.test_soa, fixP,...
        tau, sigma_a, sigma_v, criterion, lambda, p_common);
p_resp(1,:) = pre_vfirst; p_resp(3,:) = pre_afirst; p_resp(2,:) = pre_simul;

for ii = 1:3
    plot(model.test_soa, p_resp(ii,:),'-o','LineWidth',lw,'Color',cmp(ii,:),...
        'MarkerSize',mksz,'MarkerFaceColor',cmp(ii,:),'MarkerEdgeColor','none')
end

lgd = legend({'\Psi_A(s)','\Psi_S(s)','\Psi_V(s)'},'Location','bestoutside');
legend boxoff
lgd.ItemTokenSize = [10,10];

% look better
yticks([0 1])
ylabel('Probability')
xlim([min(model.test_soa),max(model.test_soa)])
xlabel('Test SOA (sec)')
xticks([min(model.test_soa),0,max(model.test_soa)])
xticklabels({'-0.5','0','0.5'})
set(gca,'TickDir','out');

flnm = 'TOJdemo2';
saveas(gca,fullfile(out_dir,flnm),'pdf')