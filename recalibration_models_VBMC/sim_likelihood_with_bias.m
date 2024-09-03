% simulate the update of bias when the likelihood does not include the
% bias, the update, or neither.

clear; close all
%% manage paths
[projectDir, ~]= fileparts(pwd);
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(pwd, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% set up

exp_trial = 250;
adapter_soa = [-700, -300:100:300, 700];
idx_zero = find(adapter_soa == 0);
tau = -30;
sigma_a = 68;
sigma_v = 60;
p_common = 0.864;
alpha = 0.002;
sigma_C1 = 20;
sigma_C2 = 250;
num_adapter_soa = numel(adapter_soa);

model_names = {'Correct for bias and update','Correct for bias only','No correction'};
%% fixed parameters

model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.4*1e3; % in second, where measurements are likely to reside

% prior
fixP.x_axis     = -model.bound_full:1:model.bound_full;
fixP.x_axis_int = -model.bound_int:1:model.bound_int;
fixP.l_window   = find(fixP.x_axis == fixP.x_axis_int(1));
fixP.r_window   = find(fixP.x_axis == fixP.x_axis_int(end));
fixP.prior_C1   = normpdf(fixP.x_axis_int, 0, sigma_C1);
fixP.prior_C2   = normpdf(fixP.x_axis_int, 0, sigma_C2);
fixP.bound_int  = model.bound_int;

% likelihood that centers around 0
idx_peak   = ceil(length(fixP.x_axis)/2);
lf = (1/(sigma_a + sigma_v)).* exp(1/sigma_a .* (fixP.x_axis(1:idx_peak)));
rf = (1/(sigma_a + sigma_v)).* exp(-1/sigma_v .* (fixP.x_axis(idx_peak+1:end)));
fixP.df_likelihood     = [lf, rf] + realmin;

% iCDF
fixP.y_criterion = sigma_v/(sigma_a + sigma_v);

%% M1: correct for bias and update

delta_tau                    = zeros(num_adapter_soa, exp_trial+1);
for tt                       = 1:exp_trial

    %% sample a noisy measurement

    % update tau
    i_tau                        = tau + delta_tau(:,tt);

    % draw y from U[0, 1] for each adaptor soa
    y                            = rand(num_adapter_soa, 1);
    bool_y                       = y < fixP.y_criterion;
    soa_m                        = NaN(num_adapter_soa,1);

    % sample from iCDF of exponential measurement distribution
    % for y smaller than threshold:
    soa_m(bool_y)                = adapter_soa(bool_y)' + i_tau(bool_y) + sigma_v .* log((sigma_a + sigma_v)/ sigma_v.* y(bool_y));
    % for y larger than threshod:
    soa_m(~bool_y)               = adapter_soa(~bool_y)' + i_tau(~bool_y) - sigma_a .* log((sigma_a + sigma_v)/sigma_a .* (1 - y(~bool_y)));
    
    %% now we turn to the observer's perspective to make inference

    % perform causal inference to map measurements to estimates
    [shat, post_C1] = cau_inf_map(soa_m - i_tau, p_common, fixP);

    % fail-safe: update the mean of the measurement if shat is not nan. nan can happen if L_C1/L_C2 = 0/0
    idx_update                   = ~isnan(shat);
    delta_tau(idx_update, tt+1)  = delta_tau(idx_update, tt) + alpha .* (shat(idx_update) - soa_m(idx_update));
    delta_tau(~idx_update, tt+1) = delta_tau(~idx_update, tt); % carry over mu from the last trial

end

% model x trial
all_recal{1} = -delta_tau;

%% M2: correct for bias, but not aware of the update

delta_tau                    = zeros(num_adapter_soa, exp_trial+1);
for tt                       = 1:exp_trial

    %% sample a noisy measurement

    % update tau
    i_tau                        = tau + delta_tau(:,tt);

    % draw y from U[0, 1] for each adaptor soa
    y                            = rand(num_adapter_soa, 1);
    bool_y                       = y < fixP.y_criterion;
    soa_m                        = NaN(num_adapter_soa,1);

    % sample from iCDF of exponential measurement distribution
    % for y smaller than threshold:
    soa_m(bool_y)                = adapter_soa(bool_y)' + i_tau(bool_y) + sigma_v .* log((sigma_a + sigma_v)/ sigma_v.* y(bool_y));
    % for y larger than threshod:
    soa_m(~bool_y)               = adapter_soa(~bool_y)' + i_tau(~bool_y) - sigma_a .* log((sigma_a + sigma_v)/sigma_a .* (1 - y(~bool_y)));
    
    %% now we turn to the observer's perspective to make inference

    % perform causal inference to map measurements to estimates
    [shat, post_C1] = cau_inf_map(soa_m - tau, p_common, fixP);

    % fail-safe: update the mean of the measurement if shat is not nan. nan can happen if L_C1/L_C2 = 0/0
    idx_update                   = ~isnan(shat);
    delta_tau(idx_update, tt+1)  = delta_tau(idx_update, tt) + alpha .* (shat(idx_update) - soa_m(idx_update));
    delta_tau(~idx_update, tt+1) = delta_tau(~idx_update, tt); % carry over mu from the last trial

end

% model x trial
all_recal{2} = -delta_tau;

%% M3: does not correct for bias nor update

delta_tau                    = zeros(num_adapter_soa, exp_trial+1);
for tt                       = 1:exp_trial

    %% sample a noisy measurement

    % update tau
    i_tau                        = tau + delta_tau(:,tt);

    % draw y from U[0, 1] for each adaptor soa
    y                            = rand(num_adapter_soa, 1);
    bool_y                       = y < fixP.y_criterion;
    soa_m                        = NaN(num_adapter_soa,1);

    % sample from iCDF of exponential measurement distribution
    % for y smaller than threshold:
    soa_m(bool_y)                = adapter_soa(bool_y)' + i_tau(bool_y) + sigma_v .* log((sigma_a + sigma_v)/ sigma_v.* y(bool_y));
    % for y larger than threshod:
    soa_m(~bool_y)               = adapter_soa(~bool_y)' + i_tau(~bool_y) - sigma_a .* log((sigma_a + sigma_v)/sigma_a .* (1 - y(~bool_y)));
    
    %% now we turn to the observer's perspective to make inference

    % perform causal inference to map measurements to estimates
    [shat, post_C1] = cau_inf_map(soa_m, p_common, fixP);

    % fail-safe: update the mean of the measurement if shat is not nan. nan can happen if L_C1/L_C2 = 0/0
    idx_update                   = ~isnan(shat);
    delta_tau(idx_update, tt+1)  = delta_tau(idx_update, tt) + alpha .* (shat(idx_update) - soa_m(idx_update));
    delta_tau(~idx_update, tt+1) = delta_tau(~idx_update, tt); % carry over mu from the last trial

end

% model x trial
all_recal{3} = -delta_tau;

%% plot zero adapter

figure;
set(gcf, 'Position',[0 0 800 300])
sgtitle(sprintf('Recalibration at zero adapter, %i trials', exp_trial))
ymax = max(cellfun(@(x) max(x(idx_zero,:)), all_recal));
ymin = min(cellfun(@(x) min(x(idx_zero,:)), all_recal));
for m = 1:3
    subplot(1,3,m)
    plot(0:exp_trial, all_recal{m}(idx_zero,:))
    xlabel('Trial')
    ylabel('Shift of PSS')
    ylim([ymin, ymax])
    xlim([0, exp_trial])
    title(model_names{m})
end

saveas(gca,fullfile(out_dir,sprintf('%i trials', exp_trial)),'png')

%% plot all adapters

figure;
set(gcf, 'Position',[0 0 800 300])
sgtitle(sprintf('Recalibration across adapters, %i trials', exp_trial))
ymax = max(cellfun(@(x) max(x(:)), all_recal));
ymin = min(cellfun(@(x) min(x(:)), all_recal));
for m = 1:3
    subplot(1,3,m)
    final_pss = all_recal{m}(:,end);
    plot(adapter_soa, final_pss)
    xlabel('Adapter SOA')
    ylabel('Shift of PSS')
    ylim([ymin, ymax])
    title(model_names{m})
    yline(0,'--')
end

saveas(gca,fullfile(out_dir,sprintf('%i trials', exp_trial)),'png')