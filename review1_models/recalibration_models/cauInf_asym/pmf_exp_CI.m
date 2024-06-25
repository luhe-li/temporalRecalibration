% This function predicts the probability of three responses by performing
% causal inference on the measurement distribution. For each possible test
% SOA, it produces a measurement. The criteria were then compared with the
% measurement distribution to obtain the probability of three responses.

% Functionally it's the same as predPMFv1, but uses simulation and
% thus is more computationally consuming.

function [p_afirst_lapse, p_simul_lapse, p_vfirst_lapse, shat_sample] = pmf_exp_CI(test_soa, fixP,...
    tau, sigma_a, sigma_v, criterion, lambda, p_common)

% --------------------------------------------------------------------------
% Outputs:
% p_afirst      :probability of reporting auditory first, a vector (1 x
%                num_test_soa); same for other two responses;
% shat_sample          : a matrix, test_soa(stimulus) x sample_count
% --------------------------------------------------------------------------

checkPlot       = 0;

%% for each possible SOA, define a measurement distribution

x_axis = fixP.x_axis_int;
% [p_afirst_lapse, p_simul_lapse, p_vfirst_lapse] = deal(NaN());

for i = 1:numel(test_soa)

    %% define measurement distribution for this soa
    soa = test_soa(i);
    [~, idx_peak] = min(abs(x_axis - (soa + tau)));
    
    lx_axis = x_axis(1:idx_peak);
    rx_axis = x_axis(idx_peak+1:end);
    lf = (1/(sigma_a + sigma_v)).* exp(1/sigma_v * (lx_axis - (soa + tau)));
    rf = (1/(sigma_a + sigma_v)).* exp(-1/sigma_a * (rx_axis - (soa + tau)));
    measurement_pdf = [lf, rf];

    if  checkPlot
        figure; plot(x_axis, measurement_pdf);
    end

    %% draw samples from the measurement distirbution
    cdf = cumsum(measurement_pdf) / sum(measurement_pdf);

    % generate random samples using inverse transform sampling
    m_samples = zeros(1, fixP.sample_count);
    for j = 1:fixP.sample_count
        u = rand(); % Uniform random number between 0 and 1
        idx = find(cdf >= u, 1); % Find the index where CDF >= u
        m_samples(j) = x_axis(idx); % Map the index back to the x_axis values
    end

    if  checkPlot
        figure;
        histogram(m_samples);
    end

    %% for each measurement sample, carry out causal inference to get a final estimate
    for j = 1:fixP.sample_count
        soa_m = m_samples(j);

        % likelihood
        [~, idx_peak] = min(abs(x_axis - soa_m));
        lx_axis = x_axis(1:idx_peak);
        rx_axis = x_axis(idx_peak+1:end);

        lf = (1/(sigma_a + sigma_v)).* exp(1/sigma_a .* (lx_axis - soa_m));
        rf = (1/(sigma_a + sigma_v)).* exp(-1/sigma_v .* (rx_axis - soa_m));
        likelihood = [lf, rf];

        protopost_C1  = fixP.prior_C1.*likelihood;
        protopost_C2  = fixP.prior_C2.*likelihood;

        % likelihood of common cause/separate causes by numerically integrating
        % the unnormalized posterior
        L_C1 = sum(protopost_C1,2);
        L_C2 = sum(protopost_C2,2);

        % the posterior probability of a common source for the SOA measurement
        post_C1                      = L_C1.* p_common ./ (L_C1 .* p_common +  L_C2 .* (1 - p_common));
        post_C2                      = 1 - post_C1;

        % mode as the intermediate estimates for common cause and separate causes
        [~,idx_C1]                   = max(protopost_C1);
        [~,idx_C2]                   = max(protopost_C2);
        shat_C1                      = fixP.x_axis_int(idx_C1)';
        shat_C2                      = fixP.x_axis_int(idx_C2)';

        % compute the final estimates: assume model averaging
        % MA: the final estimate is the sum of two intermediate estimates weighted
        % by the posterior of the correponding causal structure
        % Eq. 4 in Wozny et al., 2010
        shat_sample(i,j)                         = shat_C1 .* post_C1 + shat_C2 .* post_C2;

        if  checkPlot
            figure
            subplot(1,3,1)
            plot(likelihood)
            subplot(1,3,2); hold on
            plot(fixP.prior_C1); plot(fixP.prior_C2)
            subplot(1,3,3); hold on
            plot(protopost_C1./sum(protopost_C1)); plot(protopost_C2./sum(protopost_C2))
        end

    end

    %% compare the shat distribution to the criterion
    p_vfirst_lapse(i,:) = lambda/3 + (1-lambda).* sum(shat_sample(i,:) >= criterion)./fixP.sample_count;
    p_afirst_lapse(i,:) = lambda/3 + (1-lambda).* sum(shat_sample(i,:) < -criterion)./fixP.sample_count;
    p_simul_lapse = 1 - p_afirst_lapse - p_vfirst_lapse;


end
