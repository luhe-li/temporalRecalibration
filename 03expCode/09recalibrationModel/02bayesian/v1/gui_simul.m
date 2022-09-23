%% set up parameters
simTrial = 10;
expTrial = 250;
adaptor_soas = [-700,-300:100:300, 700]./1000;
last_recal  = NaN(numel(adaptor_soas), simTrial);

%% get free parameters
p_common = 0.3;
sigma_c1 = 0.01;
sigma_c2 = 1;
sigma_soa = 0.16;
alpha = 0.01;

%% running simulation
for i                 = 1:numel(adaptor_soas)

    adaptor_soa                   = adaptor_soas(i);% in s, adaptor_soa, fixed in session

    for t                 = 1:simTrial

        delta_s = update_recal_bayesian_MA(expTrial, adaptor_soa, ...
            p_common, sigma_soa, sigma_c1, sigma_c2, alpha);

        last_recal     = - delta_s(:, end);

    end
end

%% summarize recalibration effetct
last_recal_iSOA_error = std(last_recal, [], 2);
last_recal_iSOA       = mean(last_recal, 2);

%% plotting
figure; hold on; box off
set(gca,'FontSize',15,'linewidth',2)
set(gcf,'Position', [0, 0, 500, 400]);
errorbar(adaptor_soas, last_recal_iSOA, last_recal_iSOA_error,'.','LineWidth',2)
yline(0)
xticks(adaptor_soas)
xticklabels(adaptor_soas)
xlabel('adaptor SOA (s)')
ylabel('recalibration (s)')
title(titles{m})