function [p_afirst, p_vfirst, p_simul] = MCD_prob (soas, duration, fs, onset, stim_dura, ...
    ta, tv, tav, mu, sigma, criterion)
% ------updated 2022/06-----

% This function (1) takes in physical SOAs and creates time series signals,
% (2) pass it to MCD1.m to obtain MCD_lag, and (3) convert MCD_lag to the
% perceived probability of a_first, v_first, and simultaneous responses in
% a ternary TOJ task.

%input arguments:
%___fixed paramaters___%
%soas     : a vector, unique levels of SOA, in sec
%duration : in sec
%fs       : Hz
%onset    : the onset of the visual stimulus
%stim_dura: stimulus duration in sec

%___free paramaters___%
%ta       : a scalar, tau_a, the parameter for auditory filter in MCD1.m
%tv       : a scalar, tau_v, the parameter for visual filter in MCD1.m
%tav      : a scalar, tau_av, the parameter for modality-independent filter in MCD1.m
%mu       : mean of cumulative Gaussian 
%sigma    : sigma of cumulative Gaussian 
%criterion: a scalar, the criterion that's compared with mu to determine the perceived
%probabilities for 3 responses.

%outputs:
% p_afirst : a vector, probability of reporting audition-first for each level of SOA
% p_vfirst : a vector, probability of reporting vision-first for each level of SOA
% p_simul : a vector, probability of reporting simultaneous for each level of SOA


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
    [MCD_corr(i), MCD_lag(i), ~] = MCD1(signals,fs, ta, tv, tav);
end

%% (3) convert MCD_lag to the perceived probability 
% convert MCD_lag into the proportion of correct responses via a general
% linear model with a probit link function (assuming additive Gaussian
% noise)

p_vfirst = normcdf(mu - criterion, MCD_lag, sigma);
p_afirst = 1 - normcdf(mu + criterion, MCD_lag, sigma);
p_simul = normcdf(mu + criterion, MCD_lag, sigma) - normcdf(mu - ...
    criterion, MCD_lag, sigma);

end