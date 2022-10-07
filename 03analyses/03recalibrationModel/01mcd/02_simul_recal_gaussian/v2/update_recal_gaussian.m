function recal_effect = update_recal_gaussian(exposure_trial, soa, ...
    fs, stim_dura, s_pad, ...
    bias, sigma_soa, learning_rate, ...
    ta, tv, tav)

% input arguments:
%___fixed paramaters___%
% exposure_trial: trial number in the exposure phase
% soa           : adaptor soa, fixed in one session of the exposure phase
% fs, stim_dura, s_pad: see  MCD_corr_lag1.m

%___free paramaters___%
% bias          : in sec, the initial bias that shifts soa to internally remapped soa', could be the same value as pre_mu?
% sigma_soa     : the variance of soa measurement
% learning rate : a scalar
% duration, fs, onset, stim_dura, ta, tv, tav: see MCD_corr_lag1.m

% output:
% recal_effect  : a vector of delta_soa shift for each exposure trial

% -------
% initiate recal_effect assuming the participants are calibrated in the
% first trial, thus the recal_effect(1) = 0;
recal_effect = zeros(1, exposure_trial + 1);

% participants remap physical soa in the internal space, soa' = soa + bias
soa_ = soa + bias;

for tt = 1:exposure_trial

    % Participants formed a noisy sensory measurement soa_m. We assume that soa_m is sampled from N(soa', sigma_soa)
    soa_m = randn * sigma_soa + soa_;
    soa_m = round(soa_m, 3); % round to interger ms

    % flip the sign because MCD uses SOA = t_v - t_a, the contrary of our experiment
    soa_m = -soa_m;

    %% reconstruct signals
    duration = stim_dura + abs(soa);
    nsample = round(duration * fs + 1); % round to integer as index !!
    stim_pair_signals = zeros([2, nsample]); % reset signals

    if soa > 0 % auditory first (positive soa)
        t_a_onset = 1;
        t_a_offset = 1 + stim_dura * fs;

        t_v_onset = t_a_onset + soa * fs;
        t_v_offset = t_v_onset + stim_dura * fs;

    elseif soa < 0 % visual first (negative soa)
        t_v_onset = 1;
        t_v_offset = 1 + stim_dura * fs;

        t_a_onset = t_v_onset - soa * fs;
        t_a_offset = t_a_onset + stim_dura * fs;

    elseif soa == 0 % simultaneous
        t_v_onset = 1;
        t_v_offset = 1 + stim_dura * fs;

        t_a_onset = 1;
        t_a_offset = 1 + stim_dura * fs;
    end

    % round to integer to be used as index
    t_v_onset = round(t_v_onset); t_v_offset = round(t_v_offset);
    t_a_onset = round(t_a_onset); t_a_offset = round(t_a_offset);

    % assign impulse
    stim_pair_signals(1, t_v_onset:t_v_offset) = 1;
    stim_pair_signals(2, t_a_onset:t_a_offset) = 1;

    % add padding
    padding = zeros(2, s_pad * fs); %!!
    signals = [padding, stim_pair_signals, padding]; %!!

    % pass signals to MCD model
    nsamp=length(signals); %!!
    duration=(nsamp-1)/fs;%!!
    t=linspace(0,duration,nsamp);%!!

    stimV=signals(1,:);
    stimA=signals(2,:);

    % temporal constants
    param=[ta, tv, tav];

    % low-pass filters (Equation 1)
    fv=fft(t.*exp(-t./param(1)));%!!
    fa=fft(t.*exp(-t./param(2)));%!!
    fav=fft(t.*exp(-t./param(3)));%!!

    % early filtering
    st_v=fft(stimV).*fv;
    st_a=fft(stimA).*fa;

    % late filtering
    st_v_av=ifft(st_v.*fav);
    st_a_av=ifft(st_a.*fav);

    % xcorrelate
    u1=st_a_av.*ifft(st_v);         % Equation 2
    u2=st_v_av.*ifft(st_a);         % Equation 3

    % MCD correlation detector output
    MCD_corr_signal=u2.*u1;         % Equation 4
    MCD_corr=mean(MCD_corr_signal); % Equation 6

    % MCD lag detector output
    MCD_lag_signal=u2-u1;           % Equation 5
    MCD_lag=mean(MCD_lag_signal);   % Equation 7

    % update recal_effect based on MCD_corr and soa_m
    recal_effect(tt+1)       = recal_effect(tt) + learning_rate * MCD_corr * soa_m;

    % update soa' so that next soa_m is sampled from a shifted mu
    soa_ = soa_ + recal_effect(tt+1);
end

end