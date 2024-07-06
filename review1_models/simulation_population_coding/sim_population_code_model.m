% Parameters
N = 29; % number of neurons
x_max = 500;
SOAs = linspace(-x_max, x_max, N); % preferred SOAs of neurons (in ms)
sigma = 220.6; % width of tuning curves
alpha = 0.41; % maximal proportional gain reduction
sigma_a = 122.6; % breadth of the gain field for adaptation
G0 = 1; % unadapted response gain
adapted_SOA = 50; % SOA to which adaptation occurs

% Define the Gaussian tuning functions
f = @(SOA, SOAi) G0 * exp(-(SOA - SOAi).^2 / (2 * sigma^2));

% Define the adaptation function
G = @(SOAi, SOAa) G0 * (1 - alpha * exp(-(SOAi - SOAa).^2 / (2 * sigma_a^2)));

% SOA range for the stimulus
stimulus_SOAs = linspace(-x_max, x_max, 100);

% Initialize response matrices
responses_unadapted = zeros(N, length(stimulus_SOAs));
responses_adapted = zeros(N, length(stimulus_SOAs));

% Compute the responses for each neuron to each stimulus SOA
for i = 1:N
    for j = 1:length(stimulus_SOAs)
        responses_unadapted(i, j) = f(stimulus_SOAs(j), SOAs(i));
        responses_adapted(i, j) = G(SOAs(i), adapted_SOA) * f(stimulus_SOAs(j), SOAs(i));
    end
end

% Plot the results
figure;

% (a) Unadapted responses
subplot(1, 4, 1);
plot(stimulus_SOAs, responses_unadapted);
title('(a) Unadapted Responses');
xlabel('physical SOA (ms)');
ylabel('response');
xlim([-x_max, x_max]);

% (b) Adapted responses
subplot(1, 4, 2);
plot(stimulus_SOAs, responses_adapted);
title('(b) Adapted Responses');
xlabel('physical SOA (ms)');
ylabel('response');
xlim([-x_max, x_max]);
hold on;
plot([adapted_SOA adapted_SOA], [0 G0], 'r--'); % Mark the adapted SOA

% (c) Difference between adapted and unadapted responses
subplot(1, 4, 3);
responses_difference = responses_adapted - responses_unadapted;
plot(stimulus_SOAs, sum(responses_unadapted, 1), 'k', 'LineWidth', 2);
hold on;
plot(stimulus_SOAs, sum(responses_adapted, 1), 'r', 'LineWidth', 2);
title('(c) Population Response');
xlabel('physical SOA (ms)');
ylabel('response');
xlim([-x_max, x_max]);
legend('Unadapted', 'Adapted');
hold on;
plot([adapted_SOA adapted_SOA], [0 max(sum(responses_adapted, 1))], 'r--'); % Mark the adapted SOA

%%
% Parameters for adapted conditions
adapted_SOAs = [-100, 0, 100]; % a-lead, simultaneuous, v-lead
colors = ['b', 'k', 'r']; % Colors for the different adapted conditions

% Calculate probabilities for each adapted condition
for k = 1:length(adapted_SOAs)
    adapted_SOA = adapted_SOAs(k);
    
    % Compute adapted responses
    responses_adapted = zeros(N, length(stimulus_SOAs));
    for i = 1:N
        for j = 1:length(stimulus_SOAs)
            responses_adapted(i, j) = G(SOAs(i), adapted_SOA) * f(stimulus_SOAs(j), SOAs(i));
        end
    end
    
    % Sum the population response for adapted condition
    population_response_adapted = sum(responses_adapted, 1);
    
    % Calculate the probability of "sound first" response
    % Assuming that the probability is proportional to the population response
    p_sound_first(k, :) = population_response_adapted / max(population_response_adapted);
end

% (d) Probability of "sound first" response
subplot(1, 4, 4);
hold on;
for k = 1:length(adapted_SOAs)
    plot(stimulus_SOAs, p_sound_first(k, :), 'Color', colors(k), 'LineWidth', 2);
    scatter(stimulus_SOAs, p_sound_first(k, :), 'filled', 'MarkerFaceColor', colors(k));
end
title('(d) p(sound first)');
xlabel('physical SOA (ms)');
ylabel('p(sound first)');
xlim([-x_max, x_max]);
legend('Adapt to visual lead', 'No adaptation', 'Adapt to auditory lead', 'Location', 'Northwest');
plot([-100 -100], [0 1], 'b--');
plot([100 100], [0 1], 'r--');
hold off;