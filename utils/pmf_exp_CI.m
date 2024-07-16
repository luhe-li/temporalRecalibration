function [p_afirst_lapse, p_simul_lapse, p_vfirst_lapse, shat] = pmf_exp_CI(test_soa, fixP,...
    tau, sigma_a, sigma_v, criterion, lambda, p_common)

checkPlot = 0;

n_test_soa = numel(test_soa);
n_sample = fixP.num_sample;

%% sample noisy measurements

% draw y from U[0, 1] for each test soa
y = rand(n_test_soa, n_sample);
bool_y = y < fixP.y_criterion;
soa_m = NaN(size(y));
m_test_soa = repmat(test_soa', 1, n_sample);
m_tau = repmat(tau, n_test_soa, n_sample);

% sample from iCDF of exponential measurement distribution
% right side of the equation differs by y. For y smaller than threshold:
soa_m(bool_y)  = m_test_soa(bool_y) + m_tau(bool_y) + sigma_v .* log((sigma_a + sigma_v)./ sigma_v.* y(bool_y));
% for y larger than threshod:
soa_m(~bool_y) = m_test_soa(~bool_y) + m_tau(~bool_y) - sigma_a .* log((sigma_a + sigma_v)./ sigma_a .* (1 - y(~bool_y)));

if checkPlot
    figure; hold on
    for i = 1:15
        histogram(soa_m(i,:),50)
    end
    legend
end

%% now we turn to the observer's perspective to make inference

% find closest index of shift to use the closest propoposterior
clamped_soa_m = max(min(round(soa_m), fixP.shift_bound), -fixP.shift_bound);
idx_soa_m = clamped_soa_m + fixP.shift_bound + 1;

[protopost_C1, protopost_C2] = deal(NaN(n_test_soa, n_sample, length(fixP.x_axis_int)));
for mm = 1:n_test_soa
    
    i_idx = idx_soa_m(mm, :);
    protopost_C1(mm,:,:) = fixP.protopost_C1s(i_idx, :);
    protopost_C2(mm,:,:) = fixP.protopost_C2s(i_idx, :);
    
end

[~,idx_C1] = max(protopost_C1,[],3);
[~,idx_C2] = max(protopost_C2,[],3);

% likelihood of common cause/separate causes by numerically integrating
% the unnormalized posterior
L_C1 = sum(protopost_C1,3);
L_C2 = sum(protopost_C2,3);

% the posterior probability of a common source for the SOA measurement
post_C1    = L_C1.* p_common ./ (L_C1 .* p_common +  L_C2 .* (1 - p_common));
post_C2    = 1 - post_C1;

% mode as the intermediate estimates for common cause and separate causes
shat_C1    = fixP.x_axis_int(idx_C1);
shat_C2    = fixP.x_axis_int(idx_C2);

% compute the final estimates: assume model averaging
% MA: the final estimate is the sum of two intermediate estimates weighted
% by the posterior of the correponding causal structure
% Eq. 4 in Wozny et al., 2010
shat  = shat_C1 .* post_C1 + shat_C2 .* post_C2;

%% compare the shat distribution to the criterion

p_vfirst_lapse = lambda/3 + (1-lambda).* sum(shat >= repmat(criterion, size(shat)),2)'./n_sample + realmin;
p_afirst_lapse = lambda/3 + (1-lambda).* sum(shat < -repmat(criterion, size(shat)),2)'./n_sample + realmin;
p_simul_lapse = 1 - p_afirst_lapse - p_vfirst_lapse;

end