function [p_afirst_lapse, p_simul_lapse, p_vfirst_lapse] = pmf_exp_CI_noSim(test_soa, fixP,...
    tau, sigma_a, sigma_v, criterion, lambda, p_common)

% Assuming that the mapping function of causal-inference from measurement
% to estiamte is monotonic, i.e., larger the measurement, larger the
% estimate. This allows us to find the criteria in the measurement space
% (m_lc, m_uc), which can be used to determine response probabilities using
% closed-form CDF without simulation.

% cau_inf_map works best on a limit of shift range. We chose a reasonable
% range where measurements could happen, which is twice the stimulus SOA
% range, fixP.bound_int. Check if the bounds are on the one side of either criterion.
real_m_min = -fixP.bound_int;
real_m_max = fixP.bound_int;
est_min         = cau_inf_map(real_m_min, p_common, fixP);
est_max         = cau_inf_map(real_m_max, p_common, fixP);

if est_min < -criterion && est_max < -criterion
    % estimates lie completely below -criterion, report "A" only
    p_vfirst_lapse  = zeros(size(test_soa))+realmin;
    p_simul_lapse   = zeros(size(test_soa))+realmin;
    p_afirst_lapse  = 1 - p_vfirst_lapse - p_simul_lapse;

elseif est_min > criterion && est_max > criterion
    % estimates lie completely above criterion, report "V" only
    p_afirst_lapse  = zeros(size(test_soa))+realmin;
    p_simul_lapse   = zeros(size(test_soa))+realmin;
    p_vfirst_lapse  = 1 - p_afirst_lapse - p_simul_lapse;

elseif est_min < -criterion && est_max > criterion
    % the function of estimates intersects both criteria
    % find the measurement that corresponds to upper and lower criterion
    m_uc = fzero(@(soa_m) cau_inf_map(soa_m, p_common, fixP) - criterion, [real_m_min, real_m_max]);
    m_lc = fzero(@(soa_m) cau_inf_map(soa_m, p_common, fixP) + criterion, [real_m_min, real_m_max]);

    % the probability of each response is cdf of measurement distribution
    % evaluated at criterion
    p_afirst_lapse  = lambda/3 + (1-lambda).* expCDF(test_soa, tau, m_lc, sigma_a, sigma_v);
    p_vfirst_lapse  = lambda/3 + (1-lambda).* (1 - expCDF(test_soa, tau, m_uc, sigma_a, sigma_v));
    p_simul_lapse   = 1 - p_afirst_lapse - p_vfirst_lapse;

elseif est_min < -criterion && est_max < criterion
    % the function of estimates intersects the lower criterion -c only,
    % only report "A" or "simul"
    m_lc            = fzero(@(soa_m) cau_inf_map(soa_m, p_common, fixP) + criterion, [real_m_min, real_m_max]);

    p_afirst_lapse  = lambda/3 + (1-lambda).* expCDF(test_soa, tau, m_lc, sigma_a, sigma_v);
    p_vfirst_lapse  = zeros(size(test_soa))+realmin;
    p_simul_lapse   = 1 - p_afirst_lapse - p_vfirst_lapse;

elseif est_min > -criterion && est_max > criterion
    % the function of estimates intersects the upper criterion c only, only
    % report "V" or "simul"
    m_uc            = fzero(@(soa_m) cau_inf_map(soa_m, p_common, fixP) - criterion, [real_m_min, real_m_max]);

    p_afirst_lapse  = zeros(size(test_soa))+realmin;
    p_vfirst_lapse  = lambda/3 + (1-lambda).* (1 - expCDF(test_soa, tau, m_uc, sigma_a, sigma_v));
    p_simul_lapse   = 1 - p_afirst_lapse - p_vfirst_lapse;

elseif est_min > -criterion && est_max < criterion
    % the function of estimates lie within two criteria, only report "simul"
    p_afirst_lapse  = zeros(size(test_soa))+realmin;
    p_vfirst_lapse  = zeros(size(test_soa))+realmin;
    p_simul_lapse   = 1 - p_afirst_lapse - p_vfirst_lapse;

end

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

