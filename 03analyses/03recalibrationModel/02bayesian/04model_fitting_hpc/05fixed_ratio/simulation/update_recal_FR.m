function mu_shift = update_recal_FR(expTrial, adaptor_soa, mu_pre, ...
    sigma_soa,  alpha)

% We assume that delta_s are the cumulative shift of s' and is updated at
% the end of each exposure trial i. delta_s(1) is delta_s before adaptation
% (i.e., trial 0).
delta_s = zeros(1, expTrial + 1);

for tt = 1:expTrial

    % In each trial (1 ≤ i ≤ 250), the observer makes a noisy sensory
    % measurement of SOA (soa_m) with a standard deviation (sigma_soa),
    % centered on s', where s' = s + delta_s - mu (shifted and remapped)
    soa_m = randn * sigma_soa + adaptor_soa + delta_s(tt) - mu_pre;

    % update the mean of the measurement towards 0
    delta_s(tt+1) = delta_s(tt) + alpha * (0 - soa_m);

end

% after the exposure phase, the shift of s' is delta_s(end). delta_s(250) =
% s'_250 - s'_0 = (s - mu_250) - (s - mu_0) = mu_0 - mu_250. Since we are
% interested in the shift of mu:
mu_shift = - delta_s(end);

end