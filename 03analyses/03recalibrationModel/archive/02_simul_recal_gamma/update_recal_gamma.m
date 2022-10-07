function recal_effect = update_recal_gamma(exposure_trial, soa, bias,...
    duration, fs, onset, stim_dura,...
    ka, kv, kav, theta_av,...
    learning_rate)

% initiate recal_effect assuming the participants are calibrated in the
% first trial, thus the recal_effect(1) = 0;
recal_effect = zeros(1, exposure_trial + 1);

for t = 1:exposure_trial

    % remap visual and auditory stimulus onsets
    t_v_ = onset;
    t_a_ = onset + soa + bias;

    % calculate other parameters
    theta_v = t_v_/kv;
    theta_a = t_a_/ka;

    % make a noisy measureoment of visual and auditory stimulus onsets
    m_v = gamrnd(kv, t_v_/kv);
    m_a = gamrnd(ka, t_a_/ka);

    % obtain MCD_corr and MCD_lag based on soa_m
    [MCD_corr, MCD_lag] = MCD_corr_lag_gamma (duration, fs, stim_dura, m_v, m_a,...
        ka, theta_a, kv, theta_v, kav, theta_av);

    % update recal_effect based on MCD_corr and soa_m
    recal_effect(t+1) = recal_effect(t) + learning_rate * MCD_corr * MCD_lag;

    % update soa so that next soa_m is sampled from a shifted mu
    soa = soa + recal_effect(t+1);
end


end