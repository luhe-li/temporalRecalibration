% This script simulates the recalibration effect (i.e., the shift of mu in
% pre/post-tests) in the exposure phase by the Bayesian inference model.

clear all; clc; close all;

%% fixed parameters

% set parameters
sim_trial             = 100;
adaptor_soas          = [-700,-300:100:300, 700]./1000; % s
n_soas                = length(adaptor_soas);
n_strategy             = 3;

% set session parameters
exposure_trial        = 250;

%% free parameters

%In the Baysian causal inference model, participants have a prior belief
%about how often an auditory and a visual stimulus comes from a common
%cause. This common-cause prior is learned from past experience. A small
%common-cause prior means participants believe the two stimuli rarely come
%from the same source; a large common-cause prior means participants
%believe the two stimuli always come together from a common source.
p_common              = 0.3;

% In the case of a common source $(C = 1)$, the prior distiribution is a
% narrow Gaussian distribution centered on zero.
sig_C1                = 0.001; % in sec

% When there are two separate sources (C = 2), the prior distribution is a
% broad Gaussian distribution.
sig_C2                = 1; % in sec

% define measurement noise
sig_soa               = 0.13;

% The updates of the measurement shifts are scaled by a learning rate α
alpha                 = 0.009;

%% initiate recalibration effect for all soas

delta_s               = NaN(n_soas, sim_trial, n_strategy, exposure_trial + 1);
last_recal            = NaN(n_soas, sim_trial, n_strategy); % summarize the last recalibration effect

%% run simulation

for i                 = 1:n_soas

    soa                   = adaptor_soas(i);% in s, adaptor_soa, fixed in session

    for t                 = 1:sim_trial

        delta_s(i, t, :, :)       = update_recal_bayesian(exposure_trial, soa, ...
            p_common, sig_soa, sig_C1, sig_C2, alpha);

%         figure; hold on
%         subplot 131
%         plot(squeeze(delta_s(i, t, 1, :)))
%         subplot 132
%         plot(squeeze(delta_s(i, t, 2, :)))
%         subplot 133
%         plot(squeeze(delta_s(i, t, 3, :)))

        % delta_s on the last trial corresponds to the simulated delta_PSS,
        % therefore we only look at delta_s(end) for each simulation trial (t),
        % each adaptor SOA (i)
        last_recal(i, t, :)      = - delta_s(i, t, :, end);
        
    end
end

%% summarize and plot recalibration effect by adaptor SOA
titles={'model averaging','model selection','probability matching'};

for m = 1:3

    last_recal_iSOA_error = std(last_recal(:,:,m), [], 2);
    last_recal_iSOA       = mean(last_recal(:,:,m), 2);

    figure; hold on; box off
    set(gca,'FontSize',15,'linewidth',2)
    set(gcf,'Position', [0, 0, 500, 400]);
    errorbar(adaptor_soas, last_recal_iSOA, last_recal_iSOA_error,'.','LineWidth',2)
%     errorbar(adaptor_soas, last_recal_iSOA'./adaptor_soas, last_recal_iSOA_error,'.','LineWidth',2)
    yline(0)
    xticks(adaptor_soas)
    xticklabels(adaptor_soas)
    xlabel('adaptor SOA (s)')
    ylabel('recalibration (s)')
    title(titles{m})

%     %% plot the histogram of recalibration effect for each adaptor SOA

%     figure; hold on; box off
%     set(gca,'FontSize',15,'linewidth',2)
%     set(gcf,'Position', [0, 0, 1000, 800]);
%     for i                 = 1:n_soas
%         subplot(3,3,i)
%         histogram(last_recal(i,:,m),20,'FaceColor','k','FaceAlpha',0.3,'EdgeColor','w');
%         title(['' num2str(adaptor_soas(i)*1e3) ' ms'])
%     end
%     sgtitle('recalibration effect distribution by adaptor SOA')

end