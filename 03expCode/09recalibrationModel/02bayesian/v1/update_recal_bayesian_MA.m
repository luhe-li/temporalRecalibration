function mu = update_recal_bayesian_MA(expTrial, adaptor_soa, mu_pre, ...
    p_common, sigma_soa, sigma_c1, sigma_c2, alpha)

% compute constants (these will come in handy when writing the equations
% for the likelihood/posterior of a common cause and separate causes.
JS                    = sigma_soa^-2;
JP1                   = sigma_c1^-2;
JP2                   = sigma_c2^-2;
const1                = sigma_soa^2 + sigma_c1^2;
const2                = sigma_soa^2 + sigma_c2^2;

% We assume that ∆ s are updated at the end of each
% exposure trial i.
mu = zeros(1, expTrial + 1);

% mu_1 is mu_pre
mu(1) = mu_pre;

for tt = 1:expTrial

    % In each trial (1 ≤ i ≤ 250), observer makes a noisy sensory measurement
    % of SOA, soa_m, with a standard deviation sigma_soa centered on adaptor_soa + delta_s
    soa_m = randn * sigma_soa + (adaptor_soa + mu(tt));

    % The likelihood of a common source of SOA measurement in a trial i is:
    L_C1 = 1 / (2*pi*sqrt(const1)) * exp( -0.5 * soa_m^2 / (const1));
    L_C2 = 1 / (2*pi*sqrt(const2)) * exp( -0.5 * soa_m^2 / (const2));

    % The posterior probability of a common source for the SOA measurement in trial i is:
    post_C1 = L_C1 * p_common / (L_C1 .* p_common +  L_C2 * (1 - p_common));
    post_C2 = 1 - post_C1;

    % compute the intermediate SOA estimates
    % In the case of a common source, the MAP estimates of SOA is:
    shat_C1 = soa_m * JS / (JS + JP1);
    % In the case of two sources, the MAP estimates of SOA is:
    shat_C2 = soa_m * JS / (JS + JP2);

    %compute the final location estimates if we assume model averaging.
    %Based on this strategy, the final location estimate is the sum of the
    %two intermediate location estimates, weighted by the corresponding
    %causal structure.
    %Eq. 4 in Wozny et al., 2010
    shat = shat_C1 * post_C1 + shat_C2 * post_C2;

    % update the mean of the measurement using three methods

    mu(tt+1) = mu(tt) - alpha * (shat - soa_m);
end
end