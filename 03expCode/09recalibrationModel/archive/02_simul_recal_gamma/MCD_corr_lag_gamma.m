function [MCD_corr, MCD_lag]   = MCD_corr_lag_gamma (duration, fs, stim_dura, m_v, m_a,...
    ka, theta_a, kv, theta_v, kav, theta_av)

%% reconstruct signals

% create empty signals based on duration and fs
nsample = duration * fs + 1;
signals = zeros(2, nsample); % initiate signals

% fill in 1 when the stimuli are on
t_v_onset = m_v * fs + 1; % timepoint/sample of visual impulse onset
t_v_offset = t_v_onset + stim_dura * fs;  % timepoint/sample of visual impulse offset

t_a_onset = m_a * fs + 1; % timepoint/sample of auditory impulse onset
t_a_offset = t_a_onset + stim_dura * fs; % timepoint/sample of auditory impulse offset

signals(1, t_v_onset:t_v_offset) = 1;
signals(2, t_a_onset:t_a_offset) = 1;

%% pass signals to MCD_gamma model

t=linspace(0,duration,nsample);

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