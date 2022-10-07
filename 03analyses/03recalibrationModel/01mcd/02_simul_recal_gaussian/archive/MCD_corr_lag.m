function [MCD_corr, MCD_lag] = MCD_corr_lag (soas, duration, fs, onset, stim_dura, ...
    ta, tv, tav)
% ------updated 2022/06-----

% This function (1) takes in physical SOAs and creates time series signals,
% (2) pass it to MCD1.m to obtain MCD_corr and MCD_lag

%input arguments:
%___fixed paramaters___%
%soas     : can be a scalar or a vector, unique levels of SOA, in sec
%duration : in sec
%fs       : Hz
%onset    : the onset of the visual stimulus
%stim_dura: stimulus duration in sec

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

% create empty signals based on duration and fs
nsample = duration * fs + 1;

% initialize MCD_corr and MCD_lag
[MCD_corr, MCD_lag] = deal(NaN(1, length(soas))); 

%% (2) pass it to MCD1.m to compute MCD_corr and MCD_lag
for i = 1:length(soas)
    soa = soas(i); % for each individual SOA
    signals = zeros(2, nsample); % reset signals

    t_v_onset = onset * fs + 1; % timepoint/sample of visual impulse onset
    t_v_offset = t_v_onset + stim_dura * fs;  % timepoint/sample of visual impulse offset
    
    t_a_onset = t_v_onset - soa * fs; % timepoint/sample of auditory impulse onset
    t_a_offset = t_a_onset + stim_dura * fs; % timepoint/sample of auditory impulse offset

    signals(1, t_v_onset:t_v_offset) = 1;
    signals(2, t_a_onset:t_a_offset) = 1;
    [MCD_corr(i), MCD_lag(i), ~] = MCD1(signals, fs, ta, tv, tav);
end

end