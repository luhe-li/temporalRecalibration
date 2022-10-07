function sim_sample = sampleVector(p, nT)

% inputs
% p: the real distribution where we sample from
% nT: the number of trials to simulate

% outputs
% sim_sample: a vector (1 x nT) of simulated raw responses in each trial

s_rand = rand(1,nT);
cum_p = [0, cumsum(p)]; % added zero as the lowest boundary to create the range
sim_sample = NaN(size(s_rand));
for i = 2:length(cum_p)
    % when sample value falls into the range, sample is assigned as i-1 (because we added zero as the first index)
    sim_sample(cum_p(i - 1)<= s_rand & s_rand < cum_p(i)) = i - 1;
end