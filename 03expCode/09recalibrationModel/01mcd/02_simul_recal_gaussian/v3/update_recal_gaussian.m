function recal_effect = update_recal_gaussian(exposure_trial, soa, ...
    ts_duration, fs, stim_dura, ...
    fa, fv, fav,...
    bias, sigma_soa, learning_rate)

%% input arguments:
%---fixed paramaters---%
% exposure_trial: trial number in the exposure phase
% soa           : adaptor soa, fixed in one session of the exposure phase
% ts_duration   : in s, duration of the reconstructed time series
% fs            : sampling rate
% stim_dura     : in s, stimulus duration
% fa            : unimodal filter after fft
% fv            : unimodal filter after fft
% fav           : supramodal filter after fft

%---free paramaters---%
% bias          : in sec, the initial bias that shifts soa to internally remapped soa', could be the same value as pre_mu?
% sigma_soa     : the variance of soa measurement
% learning rate : a scalar

% output:
% recal_effect  : a vector of delta_soa shift for each exposure trial

%% ---function starts---

% initiate recal_effect assuming the participants are calibrated in the
% first trial, thus the recal_effect(1) = 0;
recal_effect = zeros(1, exposure_trial + 1);

% participants remap physical soa in the internal space, soa' = soa + bias
soa_ = soa + bias;

for tt = 1:exposure_trial

    % Participants formed a noisy sensory measurement soa_m. We assume that soa_m is sampled from N(soa', sigma_soa)
    soa_m = randn * sigma_soa + soa_;

    % flip the sign because MCD uses SOA = t_v - t_a, the contrary of our experiment
    MCD_soa_m = -soa_m;

    %% reconstruct signals

    % initiate empty time seris
    nsample = ts_duration * fs + 1;
    stimV = zeros([1, nsample]); % reset v signals
    stimA = zeros([1, nsample]); % reset v signals
    midpoint = ts_duration/2 * fs + 1; % soa are centered around this midpoint

    % note that soa_m are rounded to be split into half and be centered around the midpoint
    half_soa = round(MCD_soa_m/2 * fs);

    % find out index of a/v stimulus onset and offset
    t_v_onset = midpoint + half_soa;
    t_v_offset = t_v_onset + stim_dura * fs;

    t_a_onset = midpoint - half_soa;
    t_a_offset = t_a_onset + stim_dura * fs;

    % assign 1 to stimulus inpulse
    stimV(t_v_onset:t_v_offset) = 1;
    stimA(t_a_onset:t_a_offset) = 1;

    %% pass to MCD model

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

    %% update the mean of soa measurement

    % update recal_effect based on MCD_corr and soa_m
    recal_effect(tt+1)       = recal_effect(tt) + learning_rate * MCD_corr * soa_m;

    % update soa' so that next soa_m is sampled from a shifted mu
    soa_ = soa_ + recal_effect(tt+1);
    disp(soa_)
end

end