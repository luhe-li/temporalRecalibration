
function [p_afirst_lapse, p_simul_lapse, p_vfirst_lapse] = pmf_exp(test_soa,...
    tau, sigma_a, sigma_v, lc, uc, lambda)

% the probability of each response is cdf of measurement distribution
% evaluated at criterion
p_afirst_lapse        = lambda/3 + (1-lambda).* expCDF(test_soa, tau, lc, sigma_a, sigma_v) + realmin;
p_vfirst_lapse        = lambda/3 + (1-lambda).* (1 - expCDF(test_soa, tau, uc, sigma_a, sigma_v))+ realmin; 
p_simul_lapse         = 1 - p_afirst_lapse - p_vfirst_lapse;

% check PMF if needed
% figure; hold on
% plot(test_soa, [p_afirst_lapse; p_simul_lapse; p_vfirst_lapse])

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