% outputs

% delta_s   : a matrix (3 x exposure_trial + 1), row represents different strate
function delta_s = update_recal_bayesian(exposure_trial, soa, ...
    p_common, sig_soa, sig_C1, sig_C2, alpha)

% compute constants (these will come in handy when writing the equations
% for the likelihood/posterior of a common cause and separate causes.
JS                    = sig_soa^-2;
JP1                   = sig_C1^-2;
JP2                   = sig_C2^-2;
const1                = sig_soa^2 + sig_C1^2;
const2                = sig_soa^2 + sig_C2^2;

% For each physical stimulus onset asynchrony (SOA; s) between an
% audiovisual stimuli pair, the observer remaps it into the internal space
% with a bias ∆ s . We assume that ∆ s are updated at the end of each
% exposure trial i. Before the exposure phase (i.e., i = 0), the value of ∆
% s,0 reﬂects the default bias of SOA mapping. Thus, the internal value
% corresponding to the stimulus SOA is
delta_s = zeros(3, exposure_trial + 1); % initiate for3 models

for tt = 1:exposure_trial

    % In each trial (1 ≤ i ≤ 250), observer makes a noisy sensory measurement
    % of SOA, soa_m, which is corrupted by Gaussian-distributed sensory noise
    % with a standard deviation sigma_soa centered on soa_.
    soa_m = randn * sig_soa + soa + delta_s(:,tt);

    % The likelihood of a common source of SOA measurement in a trial i is:
    L_C1 = 1 / (2*pi*sqrt(const1)) .* exp( -0.5 * soa_m.^2 ./ (const1));
    L_C2 = 1 / (2*pi*sqrt(const2)) .* exp( -0.5 * soa_m.^2 ./ (const2));

    % The posterior probability of a common source for the SOA measurement in trial i is:
    post_C1 = L_C1 .* p_common ./ (L_C1 .* p_common +  L_C2 .* (1 - p_common));
    post_C2 = 1 - post_C1;

    % compute the intermediate SOA estimates
    % In the case of a common source, the MAP estimates of SOA is:
    shat_C1 = soa_m * JS ./ (JS + JP1);
    % In the case of two sources, the MAP estimates of SOA is:
    shat_C2 = soa_m * JS ./ (JS + JP2);

    %compute the final location estimates if we assume model averaging.
    %Based on this strategy, the final location estimate is the sum of the
    %two intermediate location estimates, weighted by the corresponding
    %causal structure.
    %Eq. 4 in Wozny et al., 2010
    shat(1) = shat_C1(1) * post_C1(1) + shat_C2(1) * post_C2(1);

    %compute the final location estimates if we assume model selection.
    %Based on this strategy, the final location estimate depends purely on
    %the causal structure that is more probable.
    %Eq. 5 in Wozny et al., 2010
    if post_C1 <=0.5
        shat(2) = shat_C1(2);
    else
        shat(2) = shat_C2(2);
    end

    %compute the final location estimates if we assume probability matching.
    %Based on this strategy, the final location estimate is the integrated
    %one with a probability of post_C1, and is the segregated one with a
    %probability of post_C2.
    %Eq. 6 in Wozny et al., 2010
    if post_C1 > rand
        shat(3) = shat_C1(3);
    else
        shat(3) = shat_C2(3);
    end

    % update the mean of the measurement using three methods
    for m = 1:3
        delta_s(m, tt+1) = delta_s(m, tt) + alpha * (shat(m) - soa_m(m));
    end

end

end