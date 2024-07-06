% Parameters
N = 100; % Number of neurons
sigma = 10; % Width of tuning function
sigma_a = 5; % Width of gain field
G0 = 1; % Unadapted response gain
alpha = 0.5; % Maximal proportional gain reduction
SOA_a = 0; % Adapted SOA
physical_SOA = -50:10:50; % Physical SOAs

% Preferred SOAs of neurons
SOAi = linspace(-50, 50, N);

% Tuning function
fi = @(SOA, SOAi) G0 * exp(-(SOA - SOAi).^2 / (2 * sigma^2));

% Adaptation model
Gi = @(SOAi) G0 * (1 - alpha * exp(-(SOAi - SOA_a).^2 / (2 * sigma_a^2)));

% Poisson noise
poisson_noise = @(lambda) poissrnd(lambda);

% Response of each neuron
R = zeros(N, length(physical_SOA));
for i = 1:N
    for j = 1:length(physical_SOA)
        SOA = physical_SOA(j);
        R(i, j) = poisson_noise(fi(SOA, SOAi(i)) * Gi(SOAi(i)));
    end
end

% Log likelihood calculation
logL = @(SOA, R) sum(R .* log(fi(SOA, SOAi))) - sum(fi(SOA, SOAi)) - sum(log(factorial(R)));

% Estimate SOA
t = 50;% soa = 50
SOA_range = -100:1:100;
log_likelihoods = arrayfun(@(SOA) logL(SOA, R(:, 11)), SOA_range, 'UniformOutput', false);


estimated_SOA = zeros(1, length(physical_SOA));
for j = 1:length(physical_SOA)
    SOA_range = -100:1:100;
    log_likelihoods = arrayfun(@(SOA) logL(SOA, R(:, j)), SOA_range, 'UniformOutput', false);
    [~, max_idx] = max(log_likelihoods);
    estimated_SOA(j) = SOA_range(max_idx);
end

% Plot
figure;
scatter(physical_SOA, estimated_SOA);
xlabel('Physical SOA (ms)');
ylabel('Estimated SOA (ms)');
title('Distribution of Estimated SOA given Physical SOA');
grid on;

% Display the plot
disp('Plot of distribution of estimated SOA given the physical SOA');