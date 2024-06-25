
function [p_afirst_lapse, p_simul_lapse, p_vfirst_lapse] = pmf_gauss(test_soa,...
    mu, sigma, lc, uc, lambda)

% the probability of each response is cdf of measurement distribution
% evaluated at criterion
p_afirst_lapse  = lambda/3 + (1-lambda).*normcdf(lc, test_soa - mu, sigma) + realmin;
p_vfirst_lapse  = lambda/3 + (1-lambda).*(1 - normcdf(uc, test_soa - mu, sigma)) + realmin;
p_simul_lapse   = ones(size(p_afirst_lapse)) - p_afirst_lapse - p_vfirst_lapse;

% check PMF if needed
% figure; hold on
% plot(test_soa, [p_afirst_lapse; p_simul_lapse; p_vfirst_lapse])

end

