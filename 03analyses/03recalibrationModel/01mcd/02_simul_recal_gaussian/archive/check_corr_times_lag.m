a = [-700, -300:100:300, 700]./1000;
for i = 1:length(a)

    duration = 16; % in s
    fs = 1e3; % hz
    onset = 9; % in sec
    stim_dura = 0.033; % in sec

    ta = 0.078;
    tv = 0.068;
    tav = 0.786;
%     learning_rate = 0.005;

    % For each stimulus with an fixed soa (soa = t_a - t_v), participants formed a
    % noisy sensory measurement m_soa with a bias in the head. We assume that
    % m_soa is sampled from N(soa + bias, sigma).
    bias = 0.06; % in sec
    sigma_soa = 0.05; % in sec

    soa_m = a(i);
    [MCD_corr(i), MCD_lag(i)]  = MCD_corr_lag (soa_m, duration, fs, onset, stim_dura, ...
        ta, tv, tav);

end


figure; hold on
plot(MCD_lag)
plot(MCD_corr)
figure
plot(MCD_corr.*a*1000)