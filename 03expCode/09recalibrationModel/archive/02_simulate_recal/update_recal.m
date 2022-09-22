function recal_effect = update_recal(exposure_trial, soa, bias, sigma_soa,...
    duration, fs, onset, stim_dura,...
    ta, tv, tav,...
    learning_rate)

% input arguments:
% exposure_trial: trial number in the exposure phase
% soa           : adaptor soa, fixed in one session of the exposure phase
% bias          : in sec, the initial bias that shifts soa to internally remapped soa', could be the same value as pre_mu?
% sigma_soa     : the variance of soa measurement
% learning rate : a scalar
% duration, fs, onset, stim_dura, ta, tv, tav: see MCD_corr_lag.m

% output:
% recal_effect  : a vector of delta_soa shift for each exposure trial

% -------
% initiate recal_effect assuming the participants are calibrated in the
% first trial, thus the recal_effect(1) = 0; 
recal_effect = zeros(1, exposure_trial + 1);

% participants remap physical soa in the internal space, soa' = soa + bias
soa_ = soa + bias;

for t = 1:exposure_trial

    % Participants formed a noisy sensory measurement soa_m. We assume that soa_m is sampled from N(soa', sigma_soa)
    soa_m = randn * sigma_soa + soa_;
    soa_m = round(soa_m, 3); % round to interger ms

    % obtain MCD_corr and MCD_lag based on soa_m
    [MCD_corr, MCD_lag]   = MCD_corr_lag (soa_m, duration, fs, onset, stim_dura, ...
        ta, tv, tav);

%     if t == 250
%     disp(soa_m)
%     disp([MCD_corr, MCD_lag])
%     end
    
    % update recal_effect based on MCD_corr and soa_m
    recal_effect(t+1)       = recal_effect(t) + learning_rate * MCD_corr * soa_m;

    % update soa' so that next soa_m is sampled from a shifted mu
    soa_ = soa_ + recal_effect(t+1);
end

end