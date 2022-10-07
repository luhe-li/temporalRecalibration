% make sense of exponential function in the MCD model

clear all; close all; clc; rng(1);

x = 0:0.1:5;
taus = 0.01:0.01:0.5;
for i = 1:numel(taus)
    tau = taus(i);
    y1(i,:) = x.*exp(-x./tau); % compared with the exponential distribution, there is no scaling facture (which should be 1/tau before exp)
    y2(i,:) = exppdf(x, tau); % a standard exponential pdf for comparison
end

figure
subplot 121
plot(x, y1)
title('MCD exponential function')
subplot 122
plot(x, y2)
title('standard exponential distribution')

% MCD exponential function's peak vary with tau, which shouldn't because we
% would expect that changing tau only affect its shift and delay. 

%% Lets' try to normalize it in the next step.

y1_norm = y1./sum(y1,2);

figure
subplot 121
plot(x, y1)
title('MCD exponential function')
subplot 122
plot(x, y1_norm)
title('normalized MCD exponential function')

% looks better, seems that exponential function only looks okay when tau is
% within a reasonable range.

%% Let's try to figure out the range for tau
clm_i = floor(linspace(1,200,numel(taus))); 
clm = gray;

% when tau gets larger, color gets lighter
figure; 
for i = 1:numel(taus)

    subplot 121; hold on
    plot(x, y1(i,:), 'Color', clm(clm_i(i),:));
    title('MCD exponential function')

    subplot 122; hold on
    plot(x, y1_norm(i,:), 'Color', clm(clm_i(i),:))
    title('normalized MCD exponential function')

end

% still looks not good enough even when tau reached its exemplar value.

%% Maybe we should scale the exponential function by its peak instead of normalization
y1_same_peak = y1./max(y1, [], 2);

figure; 
for i = 1:numel(taus)

    subplot 121; hold on
    plot(x, y1(i,:), 'Color', clm(clm_i(i),:));
    title('MCD exponential function')

    subplot 122; hold on
    plot(x, y1_same_peak(i,:), 'Color', clm(clm_i(i),:))
    title('scaled MCD exponential function')

end

% Yay!! In the right figure, when tau gets larger, the scaled exponential
% distribution shifts to the right and has a longer decay, which is
% consistent with supplementary figure 8(a) in Parise2016-wg. In other
% words, to realize recalibration we should have a knob factor k that
% scales up the visual/auditory filter and scales down the other.
