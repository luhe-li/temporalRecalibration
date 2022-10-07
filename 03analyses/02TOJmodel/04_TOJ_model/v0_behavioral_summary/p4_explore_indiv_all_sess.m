% a04_p4_plot_indiv_all_sessions

% this script plots all the sessions for one subject to check the fitting
% and fitted parameters.

clear all; close all; clc; rng(1);

% set data path
currentDir            = pwd;
exptDir               = currentDir(1:regexp(pwd,'03expCode')-1);
outDir                = [pwd '/p4_figures'];
folder1 = [exptDir '03expCode/01pretest/data'];
folder2 = [exptDir '03expCode/04posttest/data'];

addpath(genpath([exptDir '03expCode/05functions']));
addpath(genpath([exptDir '03expCode/06helperFunctions']));
addpath(genpath(folder1));
addpath(genpath(folder2));

% best fitting parameters are orgnized by subject x session
load('best_para_2.mat');

% define subject and session to use
for isub                  = 2

    figure; hold on
    set(gca, 'LineWidth', 2, 'FontSize', 20)
    figure('units','normalized','outerposition',[0 0 1 1])
    sgtitle(['sub' num2str(isub)])

    for isess = 1:9

        %% extract fittied parameters
        bestP                 = best_para{isub, isess};
        i_estP_btst           = estP_btst{isub, isess};
        i_lb_68CI = lb_68CI{isub, isess};
        i_ub_68CI = ub_68CI{isub, isess};
        % numBtst = 1000;

        %% obtain raw data

        %%%%% pre-test
        load(['pretest_sub' num2str(isub) '_session' num2str(isess) '.mat'])
        pre_s_unique          = ExpInfo.SOA; % unique SOA levels, in s
        pre_ms_unique         = pre_s_unique * 1e3; % unique SOA levels, in ms
        pre_numTrials         = ExpInfo.nTrials; % num of trials per SOA
        % inititate
        pre_r_org             = NaN(length(pre_s_unique), pre_numTrials);
        pre_respCount         = NaN(3, length(pre_s_unique));
        for i                 = 1:length(pre_s_unique)
            iSOA                  = pre_s_unique(i);
            iResp                 = Response.order(ExpInfo.trialSOA == iSOA);
            pre_r_org(i,:)        = iResp; % this matrix has a size of length(s_unique) x numTrials
            for j                 = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
                pre_respCount(j,i)    = sum(iResp == j);
            end
        end
        pre_pResp             = pre_respCount/pre_numTrials;

        %%%%% post-test
        % load data and define key parameters
        load(['posttest_sub' num2str(isub) '_session' num2str(isess) '.mat'])
        post_s_unique         = ExpInfo.SOA; % unique SOA levels, in ms
        post_ms_unique        = post_s_unique * 1e3; % unique SOA levels, in s
        post_numTrials        = ExpInfo.nTrials; % num of trials per SOA
        % inititate
        post_r_org            = NaN(length(post_s_unique), post_numTrials);
        post_respCount        = NaN(3, length(post_s_unique));
        for i                 = 1:length(post_s_unique)
            iSOA                  = post_s_unique(i);
            iResp                 = Response.order(ExpInfo.trialSOA == iSOA);
            post_r_org(i,:)       = iResp; % this matrix has a size of length(s_unique) x numTrials
            for j                 = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
                post_respCount(j,i)   = sum(iResp == j);
            end
        end
        post_pResp            = post_respCount/post_numTrials;

        %% fit PMF

        %make a finer grid for the timing difference between the auditory and the
        %visual stimulus
        SOA_finer             = linspace(pre_ms_unique(1), pre_ms_unique(end), 1000);

        % define functions
        P_Afirst              = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
        P_Vfirst              = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
        P_simultaneous        = @(SOA, mu, sig, c, lambda) ...
            1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
            - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));
        psychFunctions = {P_Vfirst; P_simultaneous; P_Afirst};

        %plug in t_diff_finer and the best-fitting mu, sigma and lapse rate to
        %function P

        for i = 1:3
            pmf(i,:) = psychFunctions{i} (SOA_finer, bestP(1), bestP(3), bestP(5), bestP(7)); % pre, best fitting
            pmf(i+3,:) = psychFunctions{i} (SOA_finer, bestP(2), bestP(4), bestP(6), bestP(8)); % post, best fitting
            errorbar_lb(i,:) = psychFunctions{i} (SOA_finer,i_lb_68CI(1), i_lb_68CI(3), i_lb_68CI(5), i_lb_68CI(7)); % pre, lower
            errorbar_lb(i+3,:) = psychFunctions{i} (SOA_finer,i_lb_68CI(2), i_lb_68CI(4), i_lb_68CI(6), i_lb_68CI(8)); % post, lower
            errorbar_ub(i,:) = psychFunctions{i} (SOA_finer,i_ub_68CI(1), i_ub_68CI(3), i_ub_68CI(5), i_ub_68CI(7)); % pre, higher
            errorbar_ub(i+3,:) = psychFunctions{i} (SOA_finer,i_ub_68CI(2), i_ub_68CI(4), i_ub_68CI(6), i_ub_68CI(8)); % post, higher
        end

        %% plotting

        % set plotting parameters
        cmp =  [229, 158, 168; 203, 227, 172; 171,223,235;...
            216, 49, 91; 175, 213, 128; 88,193,238]./255;
        alpha_value = [repmat(0.6, 1, 3), repmat(0.6, 1, 3)];
        line_pattern = {'--','--','--','-','-','-'};

        % decide jitter direction by adaptor soa
        if ExpInfo.adaptor <= 0
            jitter = 10;
        else
            jitter = -10;
        end

        subplot(3,3,isess); hold on
        set(gca, 'LineWidth', 1.5, 'FontSize', 10)

        % plot pmf
        for i = [4:6, 1:3] % so that lighter lines are at the front
            patch([SOA_finer, fliplr(SOA_finer)], [errorbar_lb(i,:), fliplr(errorbar_ub(i,:))],cmp(i,:), 'EdgeColor','none','FaceAlpha', alpha_value(i));
            %     patch([SOA_finer, fliplr(SOA_finer)], [CI_lb{i}, fliplr(CI_ub{i})],cmp(i,:), 'EdgeColor','none','FaceAlpha', alpha_value(i));
            plot(SOA_finer, pmf(i,:), line_pattern{i}, 'Color', cmp(i,:), 'LineWidth', 1.5)
        end

        % plot raw data
        for i = 1:3
            plot(pre_ms_unique + jitter, pre_pResp(i,:), 'o','MarkerSize',7,'MarkerFaceColor', 'none','MarkerEdgeColor',cmp(i+3,:),'LineWidth',1.2);
            plot(post_ms_unique - jitter, post_pResp(i,:), 'o','MarkerSize',7,'MarkerFaceColor', cmp(i+3,:),'MarkerEdgeColor',cmp(i+3,:),'LineWidth',1.2);
        end

        % plot mu
        h1                    = xline(bestP(1),':','LineWidth',1,'Color', repmat(0.3, 1, 3));
        h4                    = xline(bestP(2),'LineWidth',1,'Color', repmat(0.3, 1, 3));

        % look better
        tick_soa = [-500, -300:100:300,500];
        tick_p = [0:0.25:1];
        xticks(tick_soa)
        xticklabels(strsplit(num2str(tick_soa)))
        yticks(tick_p)
        yticklabels(strsplit(num2str(tick_p)))
        xlabel(['adaptor SOA: ' num2str(ExpInfo.adaptor*1e3) ' ms'])
        ylabel('probability')
        title_para  = sprintf('%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f', round(bestP,2));
        title(title_para)

        % save figure
        flnm = ['sub' num2str(isub) '_all_sess'];
        saveas(gca, fullfile(outDir, flnm),'epsc')
    

    end
end