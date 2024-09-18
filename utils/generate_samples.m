function samples = generate_samples(Val, mean_values, sd_values, num_sample)
% Generate samples of parameter sets given the mean and sd of group
% parameter estimates from model fitting for all models.

% Extract the parameter IDs
paraIds = Val.paraID;

% Initialize the samples matrix
num_parameters = length(paraIds);
samples = zeros(num_parameters, num_sample);

% Loop through each parameter
for i = 1:num_parameters
    paraId = paraIds{i};

    switch paraId
        case {'\tau','p_{common}', '\lambda', '\alpha'}
            % Case 1: For 'tau'
            param_samples = normrnd(mean_values(i), sd_values(i), [1, num_sample]);

        case {'\sigma_{A}', '\sigma_{V}', '\sigma', 'c', '\sigma_{C=1}', '\sigma_{C=2}'}
            % Case 2: For 'sigma_a', 'sigma_v', 'sigma', 'criterion', 'sigma_C1', 'sigma_C2'
            v = sd_values(i).^2;
            log_mu = log((mean_values(i).^2) ./ sqrt(v + mean_values(i).^2));
            log_sigma = sqrt(log(v ./ (mean_values(i).^2) + 1));
            param_samples = lognrnd(log_mu, log_sigma, [1, num_sample]);

        otherwise
            error('Unknown parameter ID: %s', paraId);
    end

    % Ensure all samples are within the bounds
    for j = 1:num_sample
        while param_samples(j) < Val.lb(i) || param_samples(j) > Val.ub(i)
            switch paraId
                case {'\tau','p_{common}', '\lambda', '\alpha'}
                    param_samples(j) = normrnd(mean_values(i), sd_values(i));

                case {'\sigma_{A}', '\sigma_{V}', '\sigma', 'c', '\sigma_{C=1}', '\sigma_{C=2}'}
                    param_samples(j) = lognrnd(log_mu, log_sigma);

                otherwise
                    error('Unknown parameter ID: %s', paraId);
            end
        end
    end

    samples(i, :) = param_samples;
end
end
