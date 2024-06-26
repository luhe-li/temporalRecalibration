function tau_shift = sim_recal_heu(exp_trial, adaptor_soa, ...
    tau, sigma_a, sigma_v, criterion, alpha)

checkPlot                    = 0;

% Initiate delta_mu as zero before the exposure phase
num_adaptor_soa              = numel(adaptor_soa);
delta_tau                    = zeros(num_adaptor_soa, exp_trial+1);

for tt = 1:exp_trial

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
    
    % probability of perceiving simultaneity given measurement
    p_simul = 1 - expCDF(soa_m, i_tau, -criterion, sigma_a, sigma_v) ...
        + expCDF(soa_m, i_tau, criterion, sigma_a, sigma_v);

    delta_tau(:,tt+1) = delta_tau(:,tt) - alpha .* p_simul .* soa_m;

end







% last tau of recalibration phase is the recalibrated tau in the post-test
% phase
tau_shift                    = delta_tau(:, end);

end

% CDF of doueble exponential distribution
% Eq.3 in García-Pérez & Alcalá-Quintana (2012)
function p_resp = expCDF(SOAs, tau, m, sigma_a, sigma_v)
    p_resp = NaN(size(SOAs));
    for i = 1:length(SOAs)
        SOA = SOAs(i); % evaluate function at each level of SOA
        if m <= SOA + tau
            p_resp(i) =   sigma_v/(sigma_a + sigma_v) * exp(1/sigma_v * (m - SOA - tau));
        else
            p_resp(i) = 1 -  sigma_a/(sigma_a + sigma_v) * exp(- 1/sigma_a * (m - SOA - tau));
        end
    end
end