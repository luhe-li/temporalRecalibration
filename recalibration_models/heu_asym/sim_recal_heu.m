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

    % probability of perceiving simultaneity given measurement
    p_simul = - expCDF(soa_m, i_tau, -criterion, sigma_a, sigma_v) ...
        + expCDF(soa_m, i_tau, criterion, sigma_a, sigma_v);

    % alpha is modulated by probability of perceiving simultaneity. If it
    % is high, the system should recalibrate more since the
    % measurement is likely due to the perepcetual misalignment
    delta_tau(:,tt+1) = delta_tau(:,tt) - alpha .* p_simul .* soa_m;

end

% last tau of recalibration phase is the recalibrated tau in the post-test
% phase
tau_shift                    = delta_tau(:, end);

end


% CDF of double exponential distribution
% Eq.3 in García-Pérez & Alcalá-Quintana (2012)
% note that SOAs and tau are both vectors
function p_resp = expCDF(SOAs, tau, m, sigma_a, sigma_v)
p_resp = NaN(size(SOAs));
for i = 1:length(SOAs)
    SOA = SOAs(i); % evaluate function at each level of SOA
    current_tau = tau(i); % corresponding tau for each SOA
    if m <= SOA + current_tau
        p_resp(i) = sigma_v / (sigma_a + sigma_v) * exp(1/sigma_v * (m - SOA - current_tau));
    else
        p_resp(i) = 1 - sigma_a / (sigma_a + sigma_v) * exp(-1/sigma_a * (m - SOA - current_tau));
    end
end

end