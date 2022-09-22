function [MCD_corr, MCD_lag, MCD_output_signals] = MCD(signals,fs,offset)
%MCD computes the outputs of the MCD model
%
%   [MCD_corr, MCD_lag] = MCD(signals,fs) calculates the mean of the
%   output of an MCD unit (Equations 6 and 7, respectively) for audiovisual
%   time-varying signals. Signals is a 2-by-m matrix, whose first row
%   contains the visual signal, and the second row contains the auditory
%   signal. Fs is the sampling frequency of the signals (Hz).
%
%   [MCD_corr, MCD_lag] = MCD(signals,fs,offset) adds an offset to the
%   signals. Default is 0.
%
%   [MCD_corr, MCD_lag, MCD_output_signals] = MCD(signals,fs) and
%   [MCD_corr, MCD_lag, MCD_output_signals] = MCD(signals,fs,offset)
%   additionally returns a 2-by-m matrix with the time-varying outputs of
%   the MCD unit. The first row returns the output of the correlation detector
%   (Equation 4), the second row returns the output of the lag detector
%   (Equation 5).
%
%   To minimize artifacts due the temporal filtering, it is recommended to
%   embed the signals in a matrix so that there will be at least 7 seconds
%   of void (i.e. zeroes) at the beginning and the end.
%
%   Written by Cesare V. Parise (2016). A full description of the model and
%   its equations can be found in the paper:
%
%   Parise, CV and Ernst, MO (2016) Correlation detection as a general
%   mechanism for multisensory integration. Nature Communications.


if nargin==2
    offset=0;
end

nsamp=length(signals);
duration=nsamp/fs;
t=linspace(0,duration,nsamp);

stimV=signals(1,:)+offset;
stimA=signals(2,:)+offset;

% temporal constants
param=[0.0873, 0.0684, 0.7859];

% low-pass filters (Equation 1)
fv=fft(t.*exp(-t./param(1)));
fa=fft(t.*exp(-t./param(2)));
fav=fft(t.*exp(-t./param(3)));

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

if nargout==3
    MCD_output_signals = [MCD_corr_signal;MCD_lag_signal];
end