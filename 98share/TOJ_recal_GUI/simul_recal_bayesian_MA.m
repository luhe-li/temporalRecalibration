clear all; close all; clc

%% set up parameters
simTrial            = 1000;
expTrial            = 250;
adaptor_soas        = [-700,-300:100:300, 700]./1000;
mu_pre              = 0.06; % s

%% get free parameters
p_common            = 0.3;
sigma_c1            = 0.01;
sigma_c2            = 0.8;
sigma_soa           = 0.16;
alpha               = 0.01;

%% initiation
mu                  = NaN(numel(adaptor_soas), simTrial, expTrial+1);
mu_shift            = NaN(numel(adaptor_soas), simTrial);

%% running simulation
for i               = 1:numel(adaptor_soas)

    adaptor_soa         = adaptor_soas(i);% in s, adaptor_soa, fixed in session

    for t               = 1:simTrial

        mu_shift(i,t)           = update_recal_bayesian_MA(expTrial, adaptor_soa, mu_pre,...
            p_common, sigma_soa, sigma_c1, sigma_c2, alpha);

    end
end

%% summarize recalibration effetct
mu_shift_iSOA_error = std(mu_shift, [], 2);
mu_shift_iSOA       = mean(mu_shift, 2);

%% plotting
figure; hold on; box off
set(gca,'FontSize',15,'linewidth',2)
set(gcf,'Position', [0, 0, 500, 400]);
errorbar(adaptor_soas, mu_shift_iSOA, mu_shift_iSOA_error,'.','LineWidth',2)
yline(0)
xticks(adaptor_soas)
xticklabels(adaptor_soas)
xlabel('adaptor SOA (s)')
ylabel('recalibration (s)')

%% calculate r^2 for each soa
model.num_bin = 100;

for i = 1:numel(adaptor_soas)
    adaptor_soa         = adaptor_soas(i);% in s, adaptor_soa, fixed in session
    i_mu_shift = mu_shift(i,:);

    % find the max and min of i_mu_shift
    shift_min = min(i_mu_shift);
    shift_max = max(i_mu_shift);
    shift_range = shift_max - shift_min;

    % define the lower and upper boundaries (i.e., min-(max-min), max+(max-min))
    shift_lb = shift_min - shift_range;
    shift_ub = shift_max + shift_range;
    delta_i_mu_shift = linspace(shift_lb, shift_ub, model.num_bin);
    binsize = diff(delta_i_mu_shift(1:2));

    % fit a Gaussian and calculate R squared
    gauss_mu = mean(i_mu_shift);
    gauss_sigma = sqrt(sum((i_mu_shift - gauss_mu).^2)./numel(i_mu_shift)); % denominator is N
    gauss_pdf = normpdf(delta_i_mu_shift, gauss_mu, gauss_sigma); % approximated gaussian pdf
    predicted_y = gauss_pdf./sum(gauss_pdf); % normalize
    delta_shift_mu_edges = [delta_i_mu_shift, delta_i_mu_shift(end) + binsize] - binsize/2; % create edges around delta_i_mu_shift
    observed_y = histcounts(i_mu_shift, delta_shift_mu_edges)./numel(i_mu_shift); % manually normalize counts to probability
    R = corr(predicted_y(:), observed_y(:));
    r2(i) = R^2;

end

figure; hold on
set(gca,'FontSize',15,'linewidth',2)
plot(adaptor_soas, r2,'ko','MarkerFaceColor','k','MarkerSize',7)
title(['r^2, sim trial = ' num2str(simTrial)])
ylim([0.9, 1])