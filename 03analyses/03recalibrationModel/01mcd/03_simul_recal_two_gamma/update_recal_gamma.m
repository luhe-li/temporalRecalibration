function recal_effect = update_recal_gamma(exposure_trial, soa, fs, stim_dura, s_pad,...
    ka, theta_a, kv, theta_v, kav, theta_av, bias, learning_rate)

% initiate recal_effect assuming the participants are calibrated in the
% first trial, thus the recal_effect(1) = 0;
recal_effect = zeros(1, exposure_trial + 1);

for t = 1:exposure_trial

    % remap visual and auditory stimulus onsets
    t_v_ = 0;
    t_a_ = 0 + soa + bias;

    % make a noisy measurement of visual and auditory stimulus onsets. We
    % assume that measurements of auditory and visual onsets are
    % distributed as shifted Gamma distributions centered around the
    % remapped onsets, t_v' and t_a'. Moreover, 
    m_v = gamrnd(kv, theta_v) + t_v_ - kv*theta_v;
    m_a = gamrnd(ka, theta_a) + t_a_ - ka*theta_a;

    % measurement of soa = m_a - m_v
    soa_m = m_a - m_v;

    % obtain MCD_corr and MCD_lag based on soa_m
    [MCD_corr, MCD_lag]   = MCD_corr_lag_gamma (soa_m, fs, stim_dura, s_pad, ...
    ka, theta_a, kv, theta_v, kav, theta_av);

    % update recal_effect based on MCD_corr and soa_m
    recal_effect(t+1) = recal_effect(t) + learning_rate * MCD_corr * (m_a - m_v);

    % update soa so that next soa_m is sampled from a shifted mu
    soa = soa + recal_effect(t+1);
end


end