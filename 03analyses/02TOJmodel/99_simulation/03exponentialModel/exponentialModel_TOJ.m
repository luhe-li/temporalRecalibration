%% exponential model fitting v1

% This script follows the exponential model in García-Pérez &
% Alcalá-Quintana (2012) to simulate a ternary TOJ task

clear all; clc; close all;

%% Equation 1, fig. 1a

% Let the arrival latencies T_v and T_a of visual and auditory
% signals be random variables with respective densities g_v and g_a given
% by the shifted exponential distributions

% x is the axis of arrival latency, where exppdf is be evaluated at each
% point.
x = [-100:400];

% SOA is delta_t, the actual/physical onset difference. SOA = s_a - s_v. We
% use vision as the reference and audition as indepdendent avriable, thus a
% minus SOA means audition-leading and a positive SOA means vision-leading.
% note that soa is a single level of soa whereas SOA is an array of soa-s.
soa = 50;

% tau is the systematic processing delay.
tau_a = 20;
tau_v = 40;

% mu is the rate parameter that defines exppdf. mu = 1/lambda.
mu_a = 30;
mu_v = 60;

% define arrival latency density function by Equation 1
g_a = exppdf(x - (soa + tau_a), mu_a);
g_v = exppdf(x - (0 + tau_v), mu_v);

% plot
figure; hold on;
plot(x, [g_a; g_v], 'LineWidth', 2);
legend({'Auditory','Visual'})
xlabel('arrival latency')
ylabel('probability')
set(gca, 'FontSize', 15)

%% Equation 2, fig. 1b

% leave only one tau, where τ = τ_a − τ_v is a processing advantage such
% that τ < 0 indicates faster auditory processing and τ > 0 indicates
% faster visual processing.
tau = tau_a - tau_v;

% define lambda
lambda_a = 1/mu_a; % e.g., 1/30
lambda_v = 1/mu_v;

% let D be the difference between these two random variables, D = T_a -
% T_v, the density distribution of D is defined as:
% if d <= SOA + tau
x = [-200:200];
x1 = x(x <= soa + tau);
d1 = lambda_a * lambda_v / (lambda_a + lambda_v) * exp(lambda_v * (x1 - soa - tau));
% if d > SOA + tau
x2 = x(x > soa + tau);
d2 = lambda_a * lambda_v / (lambda_a + lambda_v) * exp(-lambda_a * (x2 - soa - tau));

figure; hold on;
plot(x, [d1, d2], 'k', 'LineWidth', 2);
xlabel('arrival latency difference, D = t_a - t_v')
ylabel('probability')
set(gca, 'FontSize', 15)

%% Equation 3, cumulative distribution for D, F(x;soa)

% This section fixes soa and vary x. 
% define x-axis
SIGMA = [-200:200];

% the cumulative distribution for D at a specific SOA is:
for i = 1:length(SIGMA)
    sigma = SIGMA(i); % evaluate function at each level of SOA
    if sigma <= soa + tau
        p(i) =   lambda_a / (lambda_a + lambda_v) * exp(lambda_v * (sigma - soa - tau));
    else
        p(i) = 1 - lambda_v / (lambda_a + lambda_v) * exp( - lambda_a * (sigma - soa - tau));
    end
end

% plot
figure; hold on;
plot(x, p, 'k', 'LineWidth', 2);
xlabel('arrival latency difference, D = t_a - t_v')
ylabel('cumulative probability')
set(gca, 'FontSize', 15)

%% Fig 2a. Define PMF, the probability of each judgment as a function of SOA 

% In the following sections, the funcion cumulativeD.m fixes sigma and varies soa.
SOA_finer = [-300:300];

% sigma/criterion is always positive
criterion = 40; 

% probability of a-first responses are the cumulative d evaluated at a
% specific negative criterion, evaluated across each SOA; simialr for other
% two p
p_afirst = cumulativeD(SOA_finer, tau, -criterion, lambda_a, lambda_v);
p_vfirst = 1 - cumulativeD(SOA_finer, tau, criterion, lambda_a, lambda_v);
p_simul = cumulativeD(SOA_finer, tau, criterion, lambda_a, lambda_v)...
    - cumulativeD(SOA_finer, tau, -criterion, lambda_a, lambda_v);

figure; hold on;
plot(SOA_finer, [p_afirst; p_vfirst; p_simul], 'LineWidth', 2);
xlabel('arrival latency difference, D = t_a - t_v')
ylabel('cumulative probability')
legend({'A first','simultaneous', 'V first'})
set(gca, 'FontSize', 15)

%% add lapse
% if we use 1 epislon for all conditions, and define kappa as 0.5
kappa = 0.5;
epislon = 0.06; 

add_lapse = @(p1, p2, p3) (1-epislon) * p1 + epislon * kappa * p2 + epislon * kappa * p3;
p_afirst_lampse = add_lapse(p_afirst, p_vfirst, p_simul);
p_vfirst_lampse = add_lapse(p_vfirst, p_afirst, p_simul);
p_simul_lampse = add_lapse(p_simul, p_afirst, p_vfirst);

figure; hold on;
plot(SOA_finer, [p_afirst_lampse; p_vfirst_lampse; p_simul_lampse], 'LineWidth', 2);
xlabel('arrival latency difference, D = t_a - t_v')
ylabel('cumulative probability')
legend({'A first','simultaneous', 'V first'})
set(gca, 'FontSize', 15)