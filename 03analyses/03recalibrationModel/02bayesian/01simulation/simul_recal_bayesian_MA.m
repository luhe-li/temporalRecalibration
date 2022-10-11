clear all; close all; clc

%% set up parameters
simTrial            = 1000;
expTrial            = 250;
adaptor_soas        = [-700,-300:100:300, 700]./1000;
mu_pre              = 0.06; % s

%% get free parameters
p_common            = 0.3;
sigma_c1            = 0.01;
sigma_c2            = 1;
sigma_soa           = 0.16;
alpha               = 0.01;

%% initiation
mu                  = NaN(numel(adaptor_soas), simTrial, expTrial+1);
mu_shift            = NaN(numel(adaptor_soas), simTrial);

%% running simulation
for i               = 1:numel(adaptor_soas)

    adaptor_soa         = adaptor_soas(i);% in s, adaptor_soa, fixed in session

    for t               = 1:simTrial

        mu(i,t,:)           = update_recal_bayesian_MA(expTrial, adaptor_soa, mu_pre,...
            p_common, sigma_soa, sigma_c1, sigma_c2, alpha);

        mu_shift(i,t)       = mu(i, t, end) - mu(i, t, 1);

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
