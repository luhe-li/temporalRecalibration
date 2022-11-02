function [conca_samples, conca_amp, ExpInfo] = generateSound(t, ExpInfo)

% extract durations
stanDura = ExpInfo.duraStandard; % in sec
compDura = ExpInfo.trialCompDura(t); % in sec
sigma = ExpInfo.trialSigma(t); 
blankDura = (ExpInfo.trialDura - ExpInfo.jFixation(t) - ExpInfo.jISI(t)...
    - stanDura - compDura)./2; 
blankJitter = rand/10; % jitter the pre/post-stimulus blank by at most 0.1 sec

% counterbalance stimulus
switch ExpInfo.counterbalance(t)
    case 1 % standard first
        stim1 = stanDura;
        stim2 = compDura;
    case 2 % comparison first
        stim1 = compDura;
        stim2 = stanDura;
end

% generate each interval and record the actual duration for each interval
% (duration might be jittered and floored)
[noiseD1, amp1, real_dura(1)] = generateOneSound(blankDura - blankJitter, ExpInfo.bin, ExpInfo.sampRate, ExpInfo.min_noise, sigma);
[signal1, amp2, real_dura(2)] = generateOneSound(stim1, ExpInfo.bin, ExpInfo.sampRate, ExpInfo.min_signal, sigma);
[noiseISI, amp3, real_dura(3)] = generateOneSound(ExpInfo.jISI(t), ExpInfo.bin, ExpInfo.sampRate, ExpInfo.min_noise, sigma);
[signal2, amp4, real_dura(4)] = generateOneSound(stim2, ExpInfo.bin, ExpInfo.sampRate, ExpInfo.min_signal, sigma);
[noiseD2, amp5, real_dura(5)] = generateOneSound(blankDura + blankJitter, ExpInfo.bin, ExpInfo.sampRate, ExpInfo.min_noise, sigma);

% concatenate samples
conca_samples = [noiseD1, signal1, noiseISI, signal2, noiseD2];
conca_amp = [amp1, amp2, amp3, amp4, amp5];
ExpInfo.realDura(t,:) = real_dura; % record the real duration of each trial component for each trial
% figure; plot(conca_samples) 

% % normalize amplitude
% range = max(conca_samples) - min(conca_samples);
% conca_samples = conca_samples./range;

%% function section

    function [samples, amp_noise, real_duration] = generateOneSound(dur, bin_dura, sampRate, min_amp, sigma)

        % calculate time stamp in each bin
        % when the whole samples are not divisible by bin number, floor it
        % to take the integer part

        n_bin = fix(dur / bin_dura);
        bin_time_stamp = bin_dura * sampRate;
        real_duration = bin_dura * n_bin;

        % generate orginal signal (v1) for each bin
        % which is sampled from a normal dist
        bin_v1 = randn(n_bin, bin_time_stamp);

        % sample amplitude from a uniform distribution
        amp_noise = rand(1, n_bin) * sigma + min_amp;

        % multiple amplitude to each bin
        for i = 1:n_bin
            bin_v2(i,:) = bin_v1(i,:) * amp_noise(i);
        end

        % reshape to get a 1D signal sample
        samples = reshape(bin_v2', [], numel(bin_v2));

    end
end