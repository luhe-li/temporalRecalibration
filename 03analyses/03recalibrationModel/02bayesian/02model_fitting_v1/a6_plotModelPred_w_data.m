% This script plots the recalibration effect from the data and the model
% fitting.

clear all; close all; clc;

%% manage path
currentDir = pwd;
exptDir = currentDir(1:regexp(pwd,'03analyses')-1); % project folder
addpath(genpath([exptDir '/03analyses/02TOJmodel/04_TOJ_model/v0_behavioral_summary'])); % add behavioral data
outDir = [currentDir '/a6_figures'];
addpath(genpath([exptDir '/02data']));  % to load data
addpath(genpath([exptDir '/01expCodes/05expFunctions'])); % to load randomize file

load('best_para_2.mat')

%% set sub/ses to run

all_sub = 1:10;
all_ses = 1:9;
sub_slc = [1:4, 6:10];
adaptor_soa = [-0.7, -0.3: 0.1: 0.3, 0.7];

%% behavioral results

%% extract parameters

% initiate parameters sorted by soa
best_para_            = cell(all_sub(end), all_ses(end));

% sort by adaptor soa order
for s = all_sub
    [B, I]                = sort(idx_adaptor_soa(s,:));
    [best_para_{s,:}]     = deal(best_para{s,I});
end

% initiate delta parameters
delta_para = cell(1,3);

% extract parameters
for i = sub_slc
    for j = all_ses
        % delta_mu = post_mu - pre_mu
        delta_para{1}(i,j)  = (best_para_{i,j}(2) - best_para_{i,j}(1))./1000;
        % delta_sigma = post_sigma - pre_sigma
        delta_para{2}(i,j)  = (best_para_{i,j}(4) - best_para_{i,j}(3))./1000;
        % delta_c = post_c - pre_c
        delta_para{3}(i,j)  = (best_para_{i,j}(6) - best_para_{i,j}(5))./1000;
    end
end

p = 1; % only look at delta_mu, the first parameter

%% plot individual behavioral results
f1 = figure; hold on
set(gcf, 'Position', get(0, 'Screensize'));

% plot delta_mu from data
for i = 1:numel(sub_slc)
    sub = sub_slc(i);
    plots{i} = subplot(3,3,i); hold on

    % plot
    plot(all_ses, delta_para{p}(sub,:),'-ko','LineWidth',1.5);
    yline(0,'--');
%     yline(mean(delta_para{p}(sub,:)),'-','LineWidth',1.5);% mean across sessions

end
linkaxes([plots{:}],'xy')

%% plot group behavioral results

% calculate group mean
mean_delta_criterion = mean(delta_para{p}(sub_slc, :), 1, 'omitnan');
se_delta_criterion = std(delta_para{p}(sub_slc, :), [], 1, 'omitnan')./sqrt(numel(sub_slc));

f2 = figure; hold on; box off;
set(gca, 'FontSize', 20, 'LineWidth', 1.5)
set(gcf, 'Position',[10 10 500 400])
alpha = 0.6; cMAPgray = [0, 0, 0; alpha.*ones(1,3)]; colororder(cMAPgray);
e = errorbar(adaptor_soa, mean_delta_criterion, se_delta_criterion, ...
    'o','LineWidth', 1.5);
e.CapSize = 0; e.MarkerFaceColor = 'k';
yticks(linspace(-150, 150, 5))
yline(0)
ylabel('recalibration effect (ms)')

%% model fitting results

%% extract parameters

load("a1_modelFittingResults.mat")
load('RandomAdaptorOrder.mat' );

for sub = sub_slc

    % sort by adaptor soa
    i_adaptor_soa = randomAdaptorOrder(sub,:);
    [sorted_adaptor_soa, order] = sort(i_adaptor_soa);

    % use counter in ALL OUTPUT VARIABLE to order them by adaptor SOA
    adaptor = 1; 

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

        % average simulated shift_mu
        delta_mu(sub, adaptor) = mean(mu_shift);

        % mu_post = mu_pre + simulated mean(mu_shift)
        mu2(sub, adaptor) = mu1(sub,adaptor) + delta_mu(sub, adaptor);

        % manually switch to the next adaptor soa
        adaptor = adaptor + 1; 
    end
end

%% plot individual model fitting results

figure(f1)

for i = 1:numel(sub_slc)
    sub = sub_slc(i);
    axes(plots{i}); hold on
    set(gca, 'LineWidth', 2, 'FontSize', 10)

    % plot
    plot(all_ses, delta_mu(sub,:),'--o','LineWidth',1.5,'Color',[0.5, 0.5, 0.5]);

    % look better
    title(['sub' num2str(sub)])
    tick_adaptor_soa = [-700, -300:100:300, 700]./1000;
    xticks(all_sess)
    xticklabels(strsplit(num2str(tick_adaptor_soa)))
    xlabel('adaptor soa')

end
linkaxes([plots{:}],'xy')