function tau_shift = sim_recal_heu(exp_trial, adaptor_soa, ...
    tau, sigma_a, sigma_v, criterion, alpha)

% initiate delta_mu as zero before the exposure phase
num_adaptor_soa              = numel(adaptor_soa);
delta_tau                    = zeros(num_adaptor_soa, exp_trial+1);
y_criterion = sigma_v/(sigma_a + sigma_v);

for tt = 1:exp_trial

    % update tau
    i_tau                        = tau + delta_tau(:,tt);

    % draw y from U[0, 1] for each adaptor soa
    y                            = rand(num_adaptor_soa, 1);
    bool_y                       = y < y_criterion;
    soa_m                        = NaN(num_adaptor_soa,1);

    % sample measurements from iCDF of exponential distribution
    % for y smaller than threshold:
    soa_m(bool_y)                = adaptor_soa(bool_y)' + i_tau(bool_y) + sigma_v .* log((sigma_a + sigma_v)/ sigma_v.* y(bool_y));
    % for y larger than threshod:
    soa_m(~bool_y)               = adaptor_soa(~bool_y)' + i_tau(~bool_y) - sigma_a .* log((sigma_a + sigma_v)/sigma_a .* (1 - y(~bool_y)));

    % only update measurement if it's within a window
    idx_update = (i_tau - criterion < soa_m) & (soa_m < i_tau + criterion);
    delta_tau(idx_update,tt+1) = delta_tau(idx_update,tt) - alpha .* soa_m(idx_update);

end

% last tau of recalibration phase is the recalibrated tau in the post-test
% phase
tau_shift = delta_tau(:, end);

end

function p = measurementGiven0(soa_m, tau, sigma_a, sigma_v)

p = zeros(size(soa_m));
% right side of measurement distribution
bool_m = soa_m > tau;
p(bool_m) = (1/(sigma_a + sigma_v)).* exp(-1/sigma_v .* (soa_m(bool_m) - tau(bool_m)));
% left side of measurement distribution
p(~bool_m) = (1/(sigma_a + sigma_v)).* exp(1/sigma_a .* (soa_m(~bool_m) - tau(~bool_m)));

end