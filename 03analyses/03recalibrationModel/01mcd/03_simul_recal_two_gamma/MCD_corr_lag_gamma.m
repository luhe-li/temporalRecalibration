function [MCD_corr, MCD_lag]   = MCD_corr_lag_gamma (soa, fs, stim_dura, s_pad,...
    ka, theta_a, kv, theta_v, kav, theta_av)

%% (1) takes in physical SOAs and creates time series signals
% flip the sign because MCD uses SOA = t_v - t_a, the contrary of our experiment
soa = -soa;

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

elseif soa == 0
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

%% (2) pass signals to adapted GAMMA MCD model

nsamp = length(signals);
duration = (nsamp-1)/fs;
t=linspace(0,duration,nsamp);

stimV=signals(1,:);
stimA=signals(2,:);

% low-pass filters (Equation 1)
fv=fft(gampdf(t, kv, theta_v));
fa=fft(gampdf(t, ka, theta_a));
fav=fft(gampdf(t, kav, theta_av));

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

end