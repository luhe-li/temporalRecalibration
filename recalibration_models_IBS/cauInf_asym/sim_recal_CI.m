function tau_shift = sim_recal_CI(tau, sigma_a, sigma_v, p_common, alpha, model)

n_expo_trial = model.n_expo_trial;
adaptor_soa = model.adaptor_soa;

% We model the recalibration process as the shift of measurement after
% encountering the adaptor SOA in each exposure trial. 
checkPlot                    = 0;

% Initiate delta_mu as zero before the exposure phase
n_adaptor_soa              = numel(adaptor_soa);
delta_tau                    = zeros(n_adaptor_soa, n_expo_trial+1);

for tt                       = 1:n_expo_trial

    %% sample a noisy measurement

    % update tau
    i_tau                        = tau + delta_tau(:,tt);

    % draw y from U[0, 1] for each adaptor soa
    y                            = rand(n_adaptor_soa, 1);
    bool_y                       = y < model.y_criterion;
    soa_m                        = NaN(n_adaptor_soa,1);

    % sample from iCDF of exponential measurement distribution
    % for y smaller than threshold:
    soa_m(bool_y)                = adaptor_soa(bool_y)' + i_tau(bool_y) + sigma_v .* log((sigma_a + sigma_v)/ sigma_v.* y(bool_y));
    % for y larger than threshod:
    soa_m(~bool_y)               = adaptor_soa(~bool_y)' + i_tau(~bool_y) - sigma_a .* log((sigma_a + sigma_v)/sigma_a .* (1 - y(~bool_y)));
    
    %% now we turn to the observer's perspective to make inference

    % perform causal inference to map measurements to estimates
    [shat, ~] = cau_inf_map(soa_m, p_common, model);

    % fail-safe: update the mean of the measurement if shat is not nan. nan can happen if L_C1/L_C2 = 0/0
    idx_update                   = ~isnan(shat);
    delta_tau(idx_update, tt+1)  = delta_tau(idx_update, tt) + alpha .* (shat(idx_update) - soa_m(idx_update));
    delta_tau(~idx_update, tt+1) = delta_tau(~idx_update, tt); % carry over mu from the last trial

    if checkPlot
        figure;hold on
        subplot(1,2,1)
        plot(model.x_axis_int, protopost_C1);
        title('protopost_c1')
        subplot(1,2,2)
        plot(model.x_axis_int, protopost_C1./sum(protopost_C1, 2));
        title('normalized posterior_c1')
        legend

        figure;hold on
        subplot(1,2,1)
        plot(model.x_axis_int, protopost_C2);
        title('protopost_c2')
        subplot(1,2,2)
        plot(model.x_axis_int, protopost_C2./sum(protopost_C2, 2));
        title('normalized posterior_c2')
        legend
    end

end

% last tau of recalibration phase is the recalibrated tau in the post-test
% phase
tau_shift                    = delta_tau(:, end);

end