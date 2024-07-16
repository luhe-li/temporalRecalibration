function [shat, post_C1] = cau_inf_map(soa_m, p_common, fixP)

% This function takes in a measurement, and a set of free parameters,
% carries out the inference process as the recalibration phase and outpurs
% a final estimate. Thus it is a function of estimate given measurement.
% Works for soa_m as scaler or vector.

% shift default likelihood and calculates posterior
shift = round(soa_m);
numShifts = numel(shift);
clamp_shift = max(min(shift, fixP.l_window - 1),fixP.r_window - length(fixP.df_likelihood));

likelihood = zeros(numShifts, length(fixP.x_axis_int));
for i = 1:numShifts
    indices = (fixP.l_window - clamp_shift(i)):(fixP.r_window - clamp_shift(i));
    likelihood(i,:) = fixP.df_likelihood(indices);
end

protopost_C1 = bsxfun(@times, fixP.prior_C1, likelihood);
protopost_C2 = bsxfun(@times, fixP.prior_C2, likelihood);

% likelihood of common cause/separate causes by numerically integrating
% the unnormalized posterior
L_C1 = sum(protopost_C1,2);
L_C2 = sum(protopost_C2,2);

% the posterior probability of a common source for the SOA measurement
post_C1   = L_C1.* p_common ./ (L_C1 .* p_common +  L_C2 .* (1 - p_common));
post_C2   = 1 - post_C1;

% mode as the intermediate estimates for common cause and separate causes
[~,idx_C1]= max(protopost_C1,[],2);
[~,idx_C2]= max(protopost_C2,[],2);
shat_C1   = fixP.x_axis_int(idx_C1)';
shat_C2   = fixP.x_axis_int(idx_C2)';

% compute the final estimates: assume model averaging
% MA: the final estimate is the sum of two intermediate estimates weighted
% by the posterior of the correponding causal structure
% Eq. 4 in Wozny et al., 2010
shat      = shat_C1 .* post_C1 + shat_C2 .* post_C2;

end