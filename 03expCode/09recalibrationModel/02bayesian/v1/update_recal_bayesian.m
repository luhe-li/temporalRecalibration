% outputs

% mu   : 1 x exposure_trial+1, cumulative shift of mu from pretest (mu_1)
%       to posttest (mu_{exposure_trial+1})
function mu = update_recal_bayesian(exposure_trial, adptor_soa, mu_pre,...
    p_common, sig_soa, sig_C1, sig_C2, alpha)

% compute constants (these will come in handy when writing the equations
% for the likelihood/posterior of a common cause and separate causes.
JS                    = sig_soa^-2;
JP1                   = sig_C1^-2;
JP2                   = sig_C2^-2;
const1                = sig_soa^2 + sig_C1^2;
const2                = sig_soa^2 + sig_C2^2;

% For each adaptor SOA (s) that is fixed within each exposure phase, the
% observer remaps it into the internal space relative to the PSS in each
% trial. We use μi to denote the cumulative shift of PSS by the end of each
% trial, thus μ1 = μpre, and μ251 = μpost.
mu = zeros(3, exposure_trial + 1); % initiate for3 models

% mu_1 is mu_pre
mu(:,1) = mu_pre;

for tt = 1:exposure_trial

    %    Observer then makes a noisy sensory measurement of SOA, m′, which is
    %    corrupted by Gaussian-distributed sensory noise with a standard
    %    deviation σs centered on shifted adaptor SOA, s + μ
    soa_m = randn * sig_soa + adptor_soa + mu(:,tt);

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
    if post_C1 > 0.5
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
        mu(m, tt+1) = mu(m, tt) - alpha * (shat(m) - soa_m(m));
    end

end

end