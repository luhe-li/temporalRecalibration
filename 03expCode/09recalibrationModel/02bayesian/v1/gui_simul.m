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