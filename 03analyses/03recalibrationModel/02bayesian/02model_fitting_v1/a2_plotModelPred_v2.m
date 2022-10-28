% This script takes the fitting results from a1_modelFitting.m and plots
% the model prediction to check whether the best-fitting parameters are
% reasonable

clear all; clc; close all;

%% manage path
currentDir = pwd;
exptDir = currentDir(1:regexp(pwd,'03analyses')-1); % project folder
outDir = [currentDir '/a2_figures'];
addpath(genpath([exptDir '/02data']));  % to load data
addpath(genpath([exptDir '/01expCodes/05expFunctions'])); % to load randomize file

load("a1_modelFittingResults.mat")
load('RandomAdaptorOrder.mat' );

%% set sub/ses to run

all_sub = 1:10;
sub_slc = [1:4, 6:10];

[p_common, delta_mu, mu1, mu2] = deal(NaN(length(all_sub), 9)); % all sub, 9 ses

% set a single figure for delta_mu
f2 = figure; 
set(gcf, 'Position', get(0, 'Screensize'));

% set variables
adaptor_soa = [-0.7, -0.3: 0.1: 0.3, 0.7];

for sub = sub_slc

    f1 = figure; % for PMF
    set(gcf, 'Position', get(0, 'Screensize'));

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

        % average simulated shift_mu
        delta_mu(sub, adaptor) = mean(mu_shift);

        % mu_post = mu_pre + simulated mean(mu_shift)
        mu2(sub, adaptor) = mu1(sub,adaptor) + delta_mu(sub, adaptor);

        %% plot PMF predictions for each session, each subject

        %make a finer grid for the timing difference between the auditory and the
        %visual stimulus
        SOA_finer  = linspace(data(sub, ses).pre_s_unique(1), data(sub, ses).pre_s_unique(end), 1000);

        % define psychometric function (PMF)
        P_Afirst = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
        P_Vfirst = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
        P_simultaneous = @(SOA, mu, sig, c, lambda) ...
            1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
            - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

        % fit PMF with best-fitting parameters
        pmf = {P_Vfirst(SOA_finer, p(1), p(2), p(3), p(4)), ...
            P_simultaneous(SOA_finer, p(1), p(2), p(3), p(4)), ...
            P_Afirst(SOA_finer, p(1), p(2), p(3), p(4)), ...
            P_Vfirst(SOA_finer, mu2(sub, adaptor), p(5), p(6), p(4)), ...
            P_simultaneous(SOA_finer, mu2(sub, adaptor), p(5), p(6), p(4)), ...
            P_Afirst(SOA_finer, mu2(sub, adaptor), p(5), p(6), p(4))};

        % set plotting parameters
        cmp  = [229, 158, 168; 203, 227, 172; 171,223,235;...
            216, 49, 91; 175, 213, 128; 88,193,238]./255;
        alpha_value  = [repmat(0.6, 1, 3), repmat(0.6, 1, 3)];
        line_pattern   = {'--','--','--','-','-','-'};

        % decide jitter direction by adaptor soa
        if data(sub, ses).adaptor_soa <= 0
            jitter  = 0.01;
        else
            jitter = -0.01;
        end

        % initiate figure
        figure(f1)
        subplot(3,3,ses); hold on
        set(gca, 'LineWidth', 2, 'FontSize', 20)

        % plot raw data
        for i = 1:3
            dots(i) = plot(data(sub, ses).pre_s_unique + jitter, data(sub, ses).pre_pResp(i,:), 'o',...
                'MarkerSize',7,'MarkerFaceColor', 'none','MarkerEdgeColor',cmp(i+3,:),'LineWidth',1.2);
            dots(i+3) = plot(data(sub, ses).pre_s_unique - jitter, data(sub, ses).post_pResp(i,:), 'o',...
                'MarkerSize',7,'MarkerFaceColor', cmp(i+3,:),'MarkerEdgeColor',cmp(i+3,:),'LineWidth',1.2);
        end

        % plot pmf (without error bars)
        for i  = [4:6, 1:3] % lighter lines (pretest) are at the front
            lines(i)  = plot(SOA_finer, pmf{i}, line_pattern{i},'Color', cmp(i,:), 'LineWidth', 2);
        end

        % plot mu
        h1 = xline(mu1(sub,adaptor),':','LineWidth',2,'Color', repmat(0.3, 1, 3));
        h2 = xline(mu2(sub, adaptor),'LineWidth',2,'Color', repmat(0.3, 1, 3));

        % look better
        tick_soa = [-500, -300:100:300,500]./1000;
        tick_p = [0:0.25:1];
        xticks(tick_soa)
        xticklabels(strsplit(num2str(tick_soa)))
        yticks(tick_p)
        yticklabels(strsplit(num2str(tick_p)))
        xlabel('SOA (s)')
        ylabel('proportion of responses')
        title(['adaptor SOA = ' num2str(data(sub, ses).adaptor_soa) ' s'])

        % manually switch to the next adaptor soa
        adaptor = adaptor + 1;

        % save figure
        flnm = ['sub' num2str(sub)];
        saveas(gca, fullfile(outDir, flnm),'epsc')

    end

    %% plot individual recalibration effect
    figure(f2); hold on

    subplot(3, 4, sub)
    plot(adaptor_soa, delta_mu(sub,:),'-ko')
    ylim([-0.3 0.3])
    xlim([-1 1])
    yline(0)
    title(['sub' num2str(sub)])

    fignm = 'indiv_modelpred_delta_mu';
    saveas(gca, fullfile(outDir, fignm), 'epsc')

end

%% plot group recalibration effect

adaptor_soa = [-0.7, -0.3: 0.1: 0.3, 0.7];

% calculate group mean
mean_delta_criterion = mean(delta_mu(sub_slc, :), 1, 'omitnan');
se_delta_criterion = std(delta_mu(sub_slc, :), [], 1, 'omitnan')./sqrt(numel(sub_slc));

% plot
figure; hold on; box off;
set(gca, 'FontSize', 20, 'LineWidth', 1.5)
set(gcf, 'Position',[10 10 500 400])
alpha = 0.6; cMAPgray = [0, 0, 0; alpha.*ones(1,3)]; colororder(cMAPgray);

yyaxis left
e = errorbar(adaptor_soa, mean_delta_criterion, se_delta_criterion, ...
    'o','LineWidth', 1.5);
e.CapSize = 0; e.MarkerFaceColor = 'k';
% ylim([-150 150 ])
yticks(linspace(-150, 150, 5))
yline(0)
ylabel('recalibration effect (ms)')

% plot sorted gain excluding SOA=0
nnz_idx =  adaptor_soa ~= 0;
nnz_adaptor_soa = adaptor_soa(nnz_idx);
gain = delta_mu(:,nnz_idx)./nnz_adaptor_soa;
g_gain = mean(gain, 1, 'omitnan');
sem_gain = std(gain,[],1,'omitnan')./sqrt(numel(sub_slc));

% yyaxis right
% e = errorbar(nnz_adaptor_soa + 0.020, g_gain, sem_gain, ...
%     'o','LineWidth', 1.5);
% e.CapSize = 0; e.MarkerFaceColor = 'auto';
% ylim([-1.5 1.5])
% yticks(linspace(-1.5, 1.5, 5))
% ylabel('gain')

% look better
xticks(adaptor_soa)
xlim([min(adaptor_soa)-0.050, max(adaptor_soa)+0.050])
xlabel('adaptor SOA (ms)')
set(gca,'TickDir','out');

% save figure
fignm = 'group_modelpred_delta_mu';
saveas(gca, fullfile(outDir, fignm), 'epsc')
