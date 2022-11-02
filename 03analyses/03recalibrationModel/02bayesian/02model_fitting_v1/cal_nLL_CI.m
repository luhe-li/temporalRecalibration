% This function calculates and adds up the negative log likelihood (nLL)
% of the pretest and posttest

function nLL = cal_nLL_CI(mu1, sigma1, c1, lambda,... % pretest free para
    sigma2, c2, ... % posttest free para
    p_common, sigma_soa, sigma_c1, sigma_c2, alpha, ... % recal simul para
    model, data) % fixed para, data

%--------------------------------------------------------------------------
% Inputs:
% mu1       : mu in pretest % in s
% sigma1    : sigma in pretest % in s
% c1        : criterion in pretest % in s
% lambda    : lapse rate, shared in pretest & posttest % in s
% sigma2    : sigma in posttest % in s
% c2        : criterion in posttest % in s
% p_common  : common-cause prior % between [0 1]
% sigma_soa : sigma of perceiving SOA % in s
% sigma_c1  : width of the prior_C1 % in s
% sigma_c2  : width of the prior_C2s % in s
% alpha     : learning rate
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Outputs:
% nLL       : total nLL of pre and posttest
%--------------------------------------------------------------------------

%% define psychometric function (PMF)
P_Afirst = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

%% calculate pretest nLL

pre_LL = data.pre_nT_A1st*log(P_Afirst(data.pre_s_unique, mu1, sigma1, c1, lambda))'...
    + data.pre_nT_V1st*log(P_Vfirst(data.pre_s_unique, mu1, sigma1, c1, lambda))'...
    + data.pre_nT_simul*log(P_simultaneous(data.pre_s_unique, mu1, sigma1, c1, lambda))';

%% calculate posttest nLL

%% simulate shift_mu through the full exposure phase and return the final
% shift of mu
mu_shift = NaN(1, model.expo_num_sim);
for t   = 1:model.expo_num_sim
    mu_shift(t) = simulate_mu_shift_MA(model.expo_num_trial, data.adaptor_soa, mu1,...
        p_common, sigma_soa, sigma_c1, sigma_c2, alpha);
end

%% approximate the probability of shift_mu by either Gaussian or KDE,
% computes the nLL given responses in the post TOJ task

% find the max and min of mu_shift
shift_min = min(mu_shift);
shift_max = max(mu_shift);
shift_range = shift_max - shift_min;

% define the lower and upper boundaries (i.e., min-(max-min), max+(max-min))
shift_lb = shift_min - shift_range;
shift_ub = shift_max + shift_range;
delta_mu_shift = linspace(shift_lb, shift_ub, model.num_bin);
binsize = diff(delta_mu_shift(1:2));

% fit a Gaussian and calculate R squared
gauss_mu = mean(mu_shift);
gauss_sigma = sqrt(sum((mu_shift - gauss_mu).^2)./numel(mu_shift)); % denominator is N instead of N-1
gauss_pdf = normpdf(delta_mu_shift, gauss_mu, gauss_sigma); % approximated gaussian pdf
predicted_y = gauss_pdf./sum(gauss_pdf); % normalize
delta_shift_mu_edges = [delta_mu_shift, delta_mu_shift(end) + binsize] - binsize/2; % create edges around delta_mu_shift
observed_y = histcounts(mu_shift, delta_shift_mu_edges)./numel(mu_shift); % manually normalize counts to probability
R = corr(predicted_y(:), observed_y(:));
r2 = R^2;

% if r2 is larger than criterion, use this Gaussian to fit
if r2 > model.thre_r2
    pdf_delta = gauss_pdf;
else % if not, use ksd
    [pdf_delta, ~] = ksdensity(mu_shift, delta_mu_shift); % ksdensity evaluated at selected delta_mu_shift values
end
pdf_delta = pdf_delta./sum(pdf_delta);

% % plot to check Gaussian fitting
% figure;hold on;
% bar(delta_mu_shift, observed_y)
% plot(delta_mu_shift, predicted_y)
% title('approximated by Gaussian')

% % plot to check ksd fitting
% figure;hold on;
% histogram(mu_shift,'BinEdges', delta_shift_mu_edges, Normalization='probability');
% % bar(delta_mu_shift, observed_y) % bar plot with manual normaliztion
% % should be the same as the histogram plot
% plot(delta_mu_shift, pdf_delta)
% title('approximated by ksdensity')

% compute the likelihood of approxiamated delta: P(resp|delta_mu_shift, M,
% \theta)
LL_delta = NaN(1, length(delta_mu_shift));
for i = 1:numel(delta_mu_shift)
    mu2 = mu1 + delta_mu_shift(i);
    LL_delta(i) = data.post_nT_A1st*log(P_Afirst(data.post_s_unique, mu2, sigma2, c2, lambda))'...
        + data.post_nT_V1st*log(P_Vfirst(data.post_s_unique, mu2, sigma2, c2, lambda))'...
        + data.post_nT_simul*log(P_simultaneous(data.post_s_unique, mu2, sigma2, c2, lambda))';
end

% post LL is the log sum of (likelihood of approximated delta x
% probability of approximated delta)
% post_LL = log(sum(exp(LL_delta) .* pdf_delta)); 

% re-written to avoid underflow of likelihood. Note that const can be
% subtracted to the exponent and added later because it is NOT summed, and
% log(exp(const)) = const
const = max(LL_delta + log(pdf_delta)); 
post_LL = log(sum(exp(LL_delta + log(pdf_delta) - const))) + const;

%% sum the negative likelihood of pre and post test 
nLL = - pre_LL - post_LL;

end