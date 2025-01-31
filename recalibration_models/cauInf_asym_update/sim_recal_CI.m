function [tau_shift, shat, post_C1] = sim_recal_CI(exp_trial, adaptor_soa, fixP, ...
    tau, sigma_a, sigma_v, p_common, alpha)

% We model the recalibration process as the shift of measurement after
% encountering the adaptor SOA in each exposure trial. 

checkPlot                    = 0;

% Initiate delta_mu as zero before the exposure phase
num_adaptor_soa              = numel(adaptor_soa);
delta_tau                    = zeros(num_adaptor_soa, exp_trial+1);

for tt                       = 1:exp_trial

    %% sample a noisy measurement

    % update tau
    i_tau                        = tau + delta_tau(:,tt);

    % draw y from U[0, 1] for each adaptor soa
    y                            = rand(num_adaptor_soa, 1);
    bool_y                       = y < fixP.y_criterion;
    soa_m                        = NaN(num_adaptor_soa,1);

    % sample from iCDF of exponential measurement distribution
    % for y smaller than threshold:
    soa_m(bool_y)                = adaptor_soa(bool_y)' + i_tau(bool_y) + sigma_v .* log((sigma_a + sigma_v)/ sigma_v.* y(bool_y));
    % for y larger than threshod:
    soa_m(~bool_y)               = adaptor_soa(~bool_y)' + i_tau(~bool_y) - sigma_a .* log((sigma_a + sigma_v)/sigma_a .* (1 - y(~bool_y)));
    
    %% now we turn to the observer's perspective to make inference

    % perform causal inference to map measurements to estimates
    [shat, post_C1] = cau_inf_map(soa_m - i_tau, p_common, fixP);

    % fail-safe: update the mean of the measurement if shat is not nan. nan can happen if L_C1/L_C2 = 0/0
    % update to measurement but by the rate of post_C1
    idx_update                   = ~isnan(soa_m);
    delta_tau(idx_update, tt+1)  = delta_tau(idx_update, tt) - alpha .* soa_m(idx_update) .* post_C1;
    delta_tau(~idx_update, tt+1) = delta_tau(~idx_update, tt); % carry over mu from the last trial

    if checkPlot
        figure;hold on
        subplot(1,2,1)
        plot(fixP.x_axis_int, protopost_C1);
        title('protopost_c1')
        subplot(1,2,2)
        plot(fixP.x_axis_int, protopost_C1./sum(protopost_C1, 2));
        title('normalized posterior_c1')
        legend

        figure;hold on
        subplot(1,2,1)
        plot(fixP.x_axis_int, protopost_C2);
        title('protopost_c2')
        subplot(1,2,2)
        plot(fixP.x_axis_int, protopost_C2./sum(protopost_C2, 2));
        title('normalized posterior_c2')
        legend
    end

end

% last tau of recalibration phase is the recalibrated tau in the post-test
% phase
tau_shift                    = delta_tau(:, end);

end