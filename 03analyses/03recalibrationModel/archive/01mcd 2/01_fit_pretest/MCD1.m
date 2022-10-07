function [MCD_corr, MCD_lag, MCD_output_signals] = MCD1(signals, fs, ta, tv, tav)
% ------updated 2022/06-----

% Major changes to the original MCD function: MCD1 sets the offset to the
% signal to 0. It also takes ta, tv, tav as free parameters for later data
% fitting.

% MCD1 also corrects the fencepost problem on the third line: 
% duration = (nsample - 1)/fs.

% MCD parameters are in seconds.

% There's a potential issue with the model. There impulse response
% functions are not normalized. They are all t*exp(-t/tau) with no constant
% out front. As such, they all have different integrals and different
% maximum values. That probably means that the computed lag variable will
% be biased. Maybe that's a feature, but it sounds like a bug in that the
% bias is a consequence of the different tau values, but the bias likely
% will not equal the difference in the taus. Maybe that's fine.

% ------original comments-----
% MCD computes the outputs of the MCD model
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


offset=0;
nsamp=length(signals);
duration=(nsamp-1)/fs;
t=linspace(0,duration,nsamp);

stimV=signals(1,:)+offset;
stimA=signals(2,:)+offset;

% temporal constants
param=[ta, tv, tav];

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