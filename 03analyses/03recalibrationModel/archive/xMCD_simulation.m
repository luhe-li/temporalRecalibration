
% obtain MCD_lag with different SOAs
duration = 17; %s
fs = 1e3; % hz
origin = 8000; % ms, offset of the visual stimulus
stim_dura = 33; % ms, stimulus duration
soas = s_unique; % ms, SOA = t_a - t_v, negative SOA means A-lead signals

for i = 1:length(soas)
    soa = soas(i);
    signals = zeros(2, duration*1000);
    signals(1, origin:(origin + stim_dura)) = 1;
    signals(2, (origin + soa):(origin + soa + stim_dura)) = 1;
    [MCD_corr(i), MCD_lag(i), ~] = MCD(signals,fs, 0.0873, 0.0684, 0.7859);
end

figure
plot(soas, MCD_lag, 'o')
xlabel('SOA')
ylabel('MCD_{lag}')

% convert MCD_lag into the proportion of correct responses via a general
% linear model with a probit link function (assuming additive Gaussian
% noise)

p = normcdf(MCD_lag, 0, 1);
figure
plot(soas, p,'o')