%% plotting for slides

% exclude sub5
% extract key parameters from the best fitting model

clear all; close all; clc; rng(1);

% set data path
currentDir = pwd;
exptDir = currentDir(1:regexp(pwd,'03expCode')-1);
outDir = [exptDir '/02figures/02_minicon_plot4slides/07_parameter_estimates_results'];
addpath(genpath([exptDir '03ExpCode/05functions']));
addpath(genpath([exptDir '03ExpCode/01pretest/data']));
addpath(genpath([exptDir '03ExpCode/04posttest/data']));
addpath('./model_comparison_function_v5')

load('mc_results_v5_n10.mat')
load('RandomAdaptorOrder.mat' );

% extract the best model across 9 sessions based on delta_AIC
[row col] = find(deltaAIC == 0);
perBestModel = NaN(1,n_all_sub);
perBestModel(row) = col;

% define subjects and sessions to use
sub_slc = [1:4,6:10]; % exclude sub5
all_adaptor = 1:9;
adaptor_soa = [-700, -300, -200, -100,  0,  100, 200, 300, 700];

% pre-allocate
[delta_mu,mu_pre,sigma_pre,delta_sigma, delta_c] = deal(NaN(sub_slc(end), all_adaptor(end)));

for sub = sub_slc
    bestModel = perBestModel(sub);

    % sort by adaptor soa
    i_adaptor_soa = randomAdaptorOrder(sub,:);
    [sorted_adaptor_soa, order] = sort(i_adaptor_soa);

    counter = 1;
    for sess = order
        switch bestModel
            case 1 % M1
                % delta_mu = post_mu - pre_mu
                delta_mu(sub,counter) = estP{sub}{sess,bestModel}(2) - estP{sub}{sess,bestModel}(1);
                % pre_mu, pre_sigma
                mu_pre(sub,counter) = estP{sub}{sess,bestModel}(1);
                sigma_pre(sub,counter) = estP{sub}{sess,bestModel}(3);
                sigma_post(sub,counter) = estP{sub}{sess,bestModel}(4);
                sigma_ave(sub,counter) = (sigma_pre(sub,counter) + sigma_post(sub,counter))/2;

                % delta_sigma
                delta_sigma(sub,counter) = estP{sub}{sess,bestModel}(4) - estP{sub}{sess,bestModel}(3);
                % delta_criterion
                delta_c(sub,counter) = estP{sub}{sess,bestModel}(6) - estP{sub}{sess,bestModel}(5);
            case 2 % M2
                % delta_mu = post_mu - pre_mu
                delta_mu(sub,counter) = estP{sub}{sess,bestModel}(2) - estP{sub}{sess,bestModel}(1);
                % pre_mu, pre_sigma
                mu_pre(sub,counter) = estP{sub}{sess,bestModel}(1);
                sigma_pre(sub,counter) = estP{sub}{sess,bestModel}(3); % only 1 sigma
                % take only one sigma as the average sigma
                sigma_ave(sub,counter) = estP{sub}{sess,bestModel}(3); % only 1 sigma

                % delta_sigma
                delta_sigma(sub,counter) = estP{sub}{sess,bestModel}(3) - estP{sub}{sess,bestModel}(3); % only 1 sigma so should be 0
                % delta_criterion
                delta_c(sub,counter) = estP{sub}{sess,bestModel}(5) - estP{sub}{sess,bestModel}(4);
        end

        counter = counter + 1;
    end

end

%% plot causal inference 
% calculate group mean
mean_delta_criterion = mean(delta_mu, 1, 'omitnan');
se_delta_criterion = std(delta_mu, [], 1, 'omitnan')./sqrt(numel(sub_slc));

% plot
figure; hold on; box off;
set(gca, 'FontSize', 20, 'LineWidth', 2)
set(gcf, 'Position',[10 10 500 400])
alpha = 0.6; cMAPgray = [0, 0, 0; alpha.*ones(1,3)]; colororder(cMAPgray);

e = errorbar(adaptor_soa, mean_delta_criterion, se_delta_criterion, ...
    'o','MarkerSize',7, 'LineWidth', 1.5);
e.CapSize = 0;
e.MarkerFaceColor = 'k';
ylim([-150 150])
yticks(linspace(-150, 150, 5))
yline(0)
ylabel('recalibration effect (ms)')

% look better
xticks(adaptor_soa)
xlim([min(adaptor_soa)-50, max(adaptor_soa)+50])
xlabel('adaptor lag (ms)')
set(gca,'TickDir','out');

% save figure
fignm = 'result1_CI';
saveas(gca, fullfile(outDir, fignm), 'epsc')

%% plot individual sigma against delta PSS

% scatter plot: x= sigma, y = |delta_PSS|, averaged aross sessions
% calculate mean and se across sessions, for each individual
mean_mu = mean(abs(delta_mu(sub_slc,:)), 2);
se_mu = std(abs(delta_mu(sub_slc,:)), [], 2)./sqrt(numel(all_adaptor));
mean_sigma = mean(sigma_ave(sub_slc,:), 2);
se_sigma = std(sigma_ave(sub_slc,:), [], 2)./sqrt(numel(all_adaptor));

% plotting
figure;hold on; box off;
set(gca, 'FontSize', 20, 'LineWidth', 2)
set(gcf, 'Position',[10 10 500 400])
% markers = ['o','s','d','^','v','>','<','p','*','h'];
markers = repmat('o',1,10);

for s = 1:numel(sub_slc)
    idx_sub = sub_slc(s);
    e = errorbar(mean_sigma(s), mean_mu(s), ...
        se_mu(s), se_mu(s), se_sigma(s), se_sigma(s), ...
        markers(idx_sub), 'LineWidth',1.5,'MarkerSize',7, 'Color','k');
    e.CapSize = 0;  
    e.MarkerFaceColor = 'k';
end

% add correlation 
[rho pval] = corr(mean_sigma, mean_mu);
corr_text = sprintf('r^2 = %.2f, p = %.3f', rho^2, pval);
corr_text = ['\it' corr_text];
text(140, 10, corr_text, ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom', ...
     'FontSize',20,...
     'FontName','Helvetica');

% look better
xlabel('noise (ms)'); ylabel('recalibration magnitude (ms)');
yl = ylim; yticks([0:30:150]);
xl = xlim; xticks([30:30:150])
set(gca,'TickDir','out');

% save figure
flnm = 'result2_corr_sigma_dPSS';
saveas(gca, fullfile(outDir, flnm), 'epsc')