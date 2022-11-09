% This script plots the recalibration effect from the data and the model
% fitting.

clear all; close all; clc;

%% manage path
currentDir = pwd;
exptDir = currentDir(1:regexp(pwd,'03analyses')-1); % project folder
addpath(genpath([exptDir '/03analyses/02TOJmodel/04_TOJ_model/v0_behavioral_summary'])); % add behavioral data
outDir = [currentDir '/a6_figures'];
addpath(genpath([exptDir '/02data']));  % to load data
addpath(genpath([exptDir '/03analyses/01analysesFunctions']))
addpath(genpath([exptDir '/01expCodes/05expFunctions'])); % to load randomize file

load('best_para_2.mat')

%% set sub/ses to run

all_sub = 1:10;
all_ses = 1:9;
sub_slc = [1:4, 6:10];
adaptor_soa = [-0.7, -0.3: 0.1: 0.3, 0.7];
adaptor_soa_finer = linspace(min(adaptor_soa), max(adaptor_soa), 1000);

%% extract behavioral results

% initiate parameters sorted by soa
best_para_            = cell(all_sub(end), all_ses(end));

% sort by adaptor soa order
for fit_mu_shift = all_sub
    [B, I]                = sort(idx_adaptor_soa(fit_mu_shift,:));
    [best_para_{fit_mu_shift,:}]     = deal(best_para{fit_mu_shift,I});
end

% preallocate
data_para = cell(1,3);

% extract parameters
for i = sub_slc
    for j = all_ses
        % delta_mu = post_mu - pre_mu
        data_para{1}(i,j)  = (best_para_{i,j}(2) - best_para_{i,j}(1))./1000;
        % delta_sigma = post_sigma - pre_sigma
        data_para{2}(i,j)  = (best_para_{i,j}(4) - best_para_{i,j}(3))./1000;
        % delta_c = post_c - pre_c
        data_para{3}(i,j)  = (best_para_{i,j}(6) - best_para_{i,j}(5))./1000;
    end
end

%% extract model fitting results

load("a1_modelFittingResults.mat")
load('RandomAdaptorOrder.mat' );

for sub = sub_slc

    % sort by adaptor soa
    i_adaptor_soa = randomAdaptorOrder(sub,:);
    [sorted_adaptor_soa, order] = sort(i_adaptor_soa);

    for  ses = order

        %% extract the best-fitting parameters based on minNLL

        %find the index that corresponds to the minimum nLL
        [min_val, min_idx] = min(model(sub, ses).minNLL);
        %find the corresponding estimated parameters
        p = model(sub, ses).estimatedP(min_idx,:);
        % assign it to the struct to store
        fits(sub, adaptor).p = p;
        fits(sub, adaptor).min_nll = min_val;

        % extract temporary variables for this script
        p_common(sub,adaptor) = p(7);
        mu1(sub,adaptor) = p(1);

        % simulate mu_shift to obtain mu_post
        mu_shift = NaN(1, model(sub, ses).expo_num_sim);
        for t   = 1: model(sub, ses).expo_num_sim
            mu_shift(t)   = simulate_mu_shift_MA(model(sub, ses).expo_num_trial, ...
                data(sub, ses).adaptor_soa, p(1),...
                p(7), p(8), p(9), p(10), p(11));
        end
        fits(sub, adaptor).mu_shift = mu_shift;

        % obtain 95CI for mu_shift
        [fits(sub, adaptor).fit_lb, fits(sub, adaptor).fit_ub] = get95CI(mu_shift);

        % average simulated shift_mu
        delta_mu(sub, adaptor) = mean(mu_shift);

        % mu_post = mu_pre + simulated mean(mu_shift)
        mu2(sub, adaptor) = mu1(sub,adaptor) + delta_mu(sub, adaptor);

        % manually switch to the next adaptor soa
        adaptor = adaptor + 1;
    end
end


%% plot individual behavioral results

% color parameter for model fitting
clt = [0.6, 0.6, 0.6];

f1 = figure; hold on
set(gcf, 'Position', get(0, 'Screensize'));

idx_p = 1; % only look at delta_mu, the first parameter
% plot delta_mu from data
for i = 1:numel(sub_slc)
    sub = sub_slc(i);
    plots{i} = subplot(3,3,i); hold on
    set(gca, 'LineWidth', 2, 'FontSize', 15);

    % plot
    plot(adaptor_soa, data_para{idx_p}(sub,:),'ko','LineWidth',1, 'MarkerFaceColor','k');
    yline(0,'-');

end

%% plot individual model fitting results

figure(f1)

for i = 1:numel(sub_slc)
    sub = sub_slc(i);
    subplot(plots{i});

    % plot model prediction with error bars
%     patch([adaptor_soa, fliplr(adaptor_soa)], ...
%         [[fits(sub, :).fit_lb], fliplr([fits(sub, :).fit_ub])], ...
%         clt, 'FaceAlpha', 1, 'EdgeColor','none')
    plot(adaptor_soa, delta_mu(sub, :), '-','Color',clt,'LineWidth',1)

    % look better
    title(['S' num2str(sub)])
    xticks(adaptor_soa)
    xlim([min(adaptor_soa)-0.050, max(adaptor_soa)+0.050])
%     xlabel('adaptor SOA')

end
linkaxes([plots{:}],'xy')

% save figure
flnm = 'indiv_recal_fits';
saveas(gca, fullfile(outDir, flnm),'epsc')

%% plot group behavioral results

% calculate group mean
mean_delta_mu = mean(data_para{idx_p}(sub_slc, :), 1, 'omitnan');
se_delta_mu = std(data_para{idx_p}(sub_slc, :), [], 1, 'omitnan')./sqrt(numel(sub_slc));

f2 = figure; hold on; box off;
set(gca, 'FontSize', 20, 'LineWidth', 1.5)
set(gcf, 'Position',[10 10 500 400])
e = errorbar(adaptor_soa, mean_delta_mu, se_delta_mu, ...
    'o','LineWidth', 1.5);
e.CapSize = 0; e.Color = 'k'; e.MarkerFaceColor = 'k';

%% plot group model fitting results

% calculate group mean
mean_delta_mu = mean(delta_mu(sub_slc, :), 1, 'omitnan');
se_delta_mu = std(delta_mu(sub_slc, :), [], 1, 'omitnan')./sqrt(numel(sub_slc));

% calculate lower and upper bound using se
lb_delta_mu = mean_delta_mu - se_delta_mu;
ub_delta_mu = mean_delta_mu + se_delta_mu;

figure(f2)
p = patch([adaptor_soa, fliplr(adaptor_soa)], [lb_delta_mu, fliplr(ub_delta_mu)], ...
    clt, 'FaceAlpha', 0.2, 'EdgeColor','none');

yticks(linspace(-150, 150, 5))
yline(0)
ylabel('recalibration effect (ms)')
xticks(adaptor_soa)
xlim([min(adaptor_soa)-0.050, max(adaptor_soa)+0.050])
xlabel('adaptor SOA (ms)')
set(gca,'TickDir','out');

% save figure
flnm = 'group_recal_fits';
saveas(gca, fullfile(outDir, flnm),'epsc')

