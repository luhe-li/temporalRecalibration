% test_MCD_simulate_recalibration

clear all; close all; clc; rng(1);

%% set parameters
fs          = 1e3; % hz
stim_dura   = 0.033; % in sec
ts_duration = 16; % in sec, 7 x 2 padding + 2 stimulus duration

%% for a single soa -> can test with a bunch of soas later
soa = 0.2; % in sec

% temporal constants
param=[0.0873, 0.0684, 0.7859, 0.7859]; % tv, ta, tav
mcd_c = 0.1;

% initiate
nsample      = ts_duration * fs + 1;
t            = linspace(0, ts_duration, nsample);

% low-pass filters (Equation 1)
fv                          = fft(t.*exp(-t./param(1)));
fa                          = fft(t.*exp(-t./param(2)));
fav1                         = fft(t.*exp(-t./param(3)));
fav2                         = fft(t.*exp(-t./param(4)));

%filpped soa cuz definition of soa is reversed in MCD
mcd_soa = -soa; 

%% reconstruct signals

% initiate empty time seris
stimV = zeros([1, nsample]); % reset v signals
stimA = zeros([1, nsample]); % reset v signals
midpoint = ts_duration/2 * fs + 1; % soa are centered around this midpoint

% note that soa_m are rounded to be split into half and be centered around the midpoint
half_soa = round(mcd_soa/2 * fs);

% find out index of a/v stimulus onset and offset
t_v_onset = midpoint + half_soa;
t_v_offset = t_v_onset + stim_dura * fs;

t_a_onset = midpoint - half_soa;
t_a_offset = t_a_onset + stim_dura * fs;

% assign 1 to stimulus inpulse
stimV(t_v_onset:t_v_offset) = 1;
stimA(t_a_onset:t_a_offset) = 1;

%     % plot to check signals
%     figure; hold on
%     plot(stimA); plot(stimV);

% early filtering
st_v=fft(stimV).*fv;
st_a=fft(stimA).*fa;

% late filtering
st_v_av=ifft(st_v.*fav1);
st_a_av=ifft(st_a.*fav2);

% xcorrelate
u1=st_a_av.*ifft(st_v);         % Equation 2
u2=st_v_av.*ifft(st_a);         % Equation 3

% MCD correlation detector output
MCD_corr_signal=u2.*u1;         % Equation 4
MCD_corr=mean(MCD_corr_signal); % Equation 6

% MCD lag detector output
MCD_lag_signal=u2-u1;           % Equation 5
MCD_lag=mean(MCD_lag_signal);   % Equation 7

%% recalibration stage

if MCD_lag > mcd_c
   
    
end
