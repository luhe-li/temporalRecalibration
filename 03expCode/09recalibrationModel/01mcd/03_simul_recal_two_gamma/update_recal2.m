function recal_effect = update_recal_gamma(exposure_trial, soa, bias, sigma_soa,...
            duration, fs, onset, stim_dura,...
            ka, kv, kav, theta_av,...
            learning_rate)

% remap visual and auditory stimulus onsets
t_v_ = onset;
t_a_ = onset + soa + bias;

for i = 1:exposure_trial

    % make a noisy measureoment of visual and auditory stimulus onsets
    m_v = gamrnd(kv, t_v_/kv);
    m_a = gamrnd(ka, t_a_/ka);

    
end



end