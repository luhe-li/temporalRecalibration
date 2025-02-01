function [p_afirst_lapse, p_simul_lapse, p_vfirst_lapse] = pmf_exp(test_soa,...
    tau, sigma_a, sigma_v, lc, uc, lambda)

% the probability of each response is cdf of measurement distribution
% evaluated at criterion
p_afirst_lapse        = lambda/3 + (1-lambda).* expCDF(test_soa, tau, lc, sigma_a, sigma_v) + realmin;
p_vfirst_lapse        = lambda/3 + (1-lambda).* (1 - expCDF(test_soa, tau, uc, sigma_a, sigma_v)) + realmin; 
p_simul_lapse         = 1 - p_afirst_lapse - p_vfirst_lapse;

end

% CDF of doueble exponential distribution
% Eq.3 in García-Pérez & Alcalá-Quintana (2012)
function p_resp = expCDF(SOAs, tau, m, sigma_a, sigma_v)
    p_resp = zeros(size(SOAs));
    for i = 1:length(SOAs)
        SOA = SOAs(i);
        delta = m - SOA - tau;
        if delta <= 0
            p_resp(i) = sigma_v / (sigma_a + sigma_v) * exp(delta / sigma_v);
        else
            p_resp(i) = 1 - sigma_a / (sigma_a + sigma_v) * exp(-delta / sigma_a);
        end
    end
    % ensure p_resp is within [0, 1]
    p_resp = max(min(p_resp, 1), 0);
end