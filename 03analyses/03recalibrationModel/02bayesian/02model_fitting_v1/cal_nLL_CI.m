% This function calculates and adds up the negative log likelihood (nLL)
% of the pretest and posttest

function nLL = cal_nLL_CI(mu1, sigma1, c1, lambda,... % pretest free para
    sigma2, c2, ... % posttest free para
    p_common, sigma_soa, sigma_c1, sigma_c2, alpha, ... % recal simul para
    model, data) % fixed para, data

%--------------------------------------------------------------------------
% Inputs:

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Outputs:

%--------------------------------------------------------------------------

%% 
% define psychometric function (PMF)
P_Afirst = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

%% calculate pretest nLL

pre_nLL = -data.pre_nT_A1st*log(P_Afirst(data.pre_ms_unique, mu1, sigma1, c1, lambda))'...
        -data.pre_nT_V1st*log(P_Vfirst(data.pre_ms_unique, mu1, sigma1, c1, lambda))'...
        -data.pre_nT_simul*log(P_simultaneous(data.pre_ms_unique, mu1, sigma1, c1, lambda))';

%% calculate posttest nLL
clear all; clc;

model.num_bin=100;
model.expo_num_sim = 1000;
model.expo_num_trial = 250;
data.adaptor_soa = 0.3;
mu1 = 0.06;
p_common            = 0.3;
sigma_c1            = 0.01;
sigma_c2            = 1;
sigma_soa           = 0.16;
alpha               = 0.01;

% simulate shift_mu through the full exposure phase and return the final
% shift of mu
for t   = 1:model.expo_num_sim 
    mu_shift(t)   = simulate_mu_shift_MA(model.expo_num_trial, data.adaptor_soa, mu1,...
        p_common, sigma_soa, sigma_c1, sigma_c2, alpha);
end

% approximate the probability of shift_mu by either Gaussian or KDE,
% computes the nLL given responses in the post TOJ task

% find the max and min of mu_shift
shift_min = min(mu_shift);
shift_max = max(mu_shift);
shift_range = shift_max - shift_min;
% define the lower and upper boundaries (i.e., min-(max-min), max+(max-min))
shift_lb = shift_min - shift_range;
shift_ub = shift_max + shift_range;
delta_shift_mu = linspace(shift_lb, shift_ub, model.num_bin);
binsize = diff(delta_shift_mu(1:2));
% fit a Gaussian and calculate R squared
gauss_mu = mean(mu_shift);
gauss_sigma = sqrt(sum((mu_shift - gauss_mu).^2)./numel(mu_shift));
gauss_pdf = normpdf(delta_shift_mu, gauss_mu, gauss_sigma);
predicted_y = gauss_pdf./sum(gauss_pdf);
delta_shift_mu_edges = [delta_shift_mu, delta_shift_mu(end) + binsize] - binsize/2;
observed_y = histcounts(mu_shift, delta_shift_mu_edges)./numel(mu_shift);
R = corr(predicted_y(:), observed_y(:));
rSquare = R^2;

figure;hold on;
bar(delta_shift_mu, observed_y)
plot(delta_shift_mu, predicted_y)


end