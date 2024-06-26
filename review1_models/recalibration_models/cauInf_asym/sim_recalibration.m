function [tau_shift, post_C1, shat] = sim_recalibration(exp_trial, adaptor_soa, fixP, ...
    tau, sigma_a, sigma_v, p_common, alpha)

checkPlot                    = 0;

% We model the recalibration process as the shift of measurement after
% encountering the adaptor SOA in each exposure trial. Initiate delta_mu as
% zero before the exposure phase.
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
    % right side of the equation differs by y. For y smaller than threshold:
    soa_m(bool_y)                = adaptor_soa(bool_y)' + i_tau(bool_y) + sigma_v .* log((sigma_a + sigma_v)/ sigma_v.* y(bool_y));
    % for y larger than threshod:
    soa_m(~bool_y)               = adaptor_soa(~bool_y)' + i_tau(~bool_y) - sigma_a .* log((sigma_a + sigma_v)/sigma_a .* (1 - y(~bool_y)));
    
    %% now we turn to the observer's perspective to make inference

    % find closest index of shift to use the closest propoposterior
    clamped_soa_m = max(min(round(soa_m), fixP.shift_bound), -fixP.shift_bound);
    idx_soa_m = clamped_soa_m  + fixP.shift_bound + 1;

    % shift default likelihood and calculates posterior for each adaptor
    % soa
    protopost_C1 = fixP.protopost_C1s(idx_soa_m, :);
    protopost_C2 = fixP.protopost_C2s(idx_soa_m, :);

    [~,idx_C1] = max(protopost_C1,[],2);
    [~,idx_C2] = max(protopost_C2,[],2);

    % likelihood of common cause/separate causes by numerically integrating
    % the unnormalized posterior
    L_C1 = sum(protopost_C1,2);
    L_C2 = sum(protopost_C2,2);

    % the posterior probability of a common source for the SOA measurement
    post_C1                      = L_C1.* p_common ./ (L_C1 .* p_common +  L_C2 .* (1 - p_common));
    post_C2                      = 1 - post_C1;

    % mode as the intermediate estimates for common cause and separate causes
    shat_C1                      = fixP.x_axis_int(idx_C1)';
    shat_C2                      = fixP.x_axis_int(idx_C2)';

    % compute the final estimates: assume model averaging
    % MA: the final estimate is the sum of two intermediate estimates weighted
    % by the posterior of the correponding causal structure
    % Eq. 4 in Wozny et al., 2010
    shat                         = shat_C1 .* post_C1 + shat_C2 .* post_C2;

    % fail-safe: update the mean of the measurement if shat is not nan. nan can happen if L_C1/L_C2 = 0/0
    idx_update                   = ~isnan(shat);
    delta_tau(idx_update, tt+1)  = delta_tau(idx_update, tt) + alpha .* (shat(idx_update) - soa_m(idx_update));
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

tau_shift                    = delta_tau(:, end);

end