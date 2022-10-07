function [MCD_corr, MCD_lag] = MCD_corr_lag1 (soas, fs, stim_dura, s_pad,...
    ta, tv, tav)
% ------updated 2022/06-----

% This function (1) takes in physical SOAs and creates time series signals,
% (2) pass it to MCD1.m to obtain MCD_corr and MCD_lag

%input arguments:
%___fixed paramaters___%
%soas     : can be a scalar or a vector, unique levels of SOA, in sec
%fs       : Hz
%stim_dura: stimulus duration in sec
%s_pad    : the time points of zero padding before and after the
%           audiovisual stimulus pair. Suggesed value is 7.5 s.

%___free paramaters___%
%ta       : a scalar, tau_a, the parameter for auditory filter in MCD1.m
%tv       : a scalar, tau_v, the parameter for visual filter in MCD1.m
%tav      : a scalar, tau_av, the parameter for modality-independent filter in MCD1.m

%outputs:
% MCD_corr: a vector, MCD_corr for each
% MCD_lag : a vector, probability of reporting vision-first for each level of SOA

%% (1) takes in physical SOAs and creates time series signals
% flip the sign because MCD uses SOA = t_v - t_a, the contrary of our experiment
soas = -soas;

% initialize MCD_corr and MCD_lag
[MCD_corr, MCD_lag] = deal(NaN(1, length(soas)));

%% (2) pass it to MCD1.m to compute MCD_corr and MCD_lag
for i = 1:length(soas)
    soa = soas(i); % for each individual SOA

    % stimulus duration
    duration = stim_dura + abs(soa);
    nsample = round(duration * fs + 1); % round to integer as index
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
    padding = zeros(2, s_pad * fs);
    signals = [padding, stim_pair_signals, padding];

    % pass signals to MCD model
    [MCD_corr(i), MCD_lag(i), ~] = MCD1(signals, fs, ta, tv, tav);
end

end