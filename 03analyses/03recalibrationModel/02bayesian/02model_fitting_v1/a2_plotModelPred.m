% This script takes the fitting results from a1_modelFitting.m and plots
% the model prediction to check whether the best-fitting parameters are
% reasonable

clear all; clc;

load("a1_modelFittingResults.mat")

% %% check distribution of fitted parameters
% figure; hold on
% set(gcf, 'Position', get(0, 'Screensize'));
% set(gca, 'LineWidth', 1, 'FontSize', 15)
% for i = 1:model.numPara
%     subplot(3,4,i)
%     histogram(model.estimatedP(:,i), 10)
%     title(model.paraID{i})
% end

currentDir = pwd;
outDir = [currentDir '/a2_figures'];
all_sub = 1:10;
all_ses = 1:9;

for sub = all_sub

    figure; set(gcf, 'Position', get(0, 'Screensize'));

    for  ses = all_ses

        %% extract the best-fitting parameters based on minNLL
        %find the index that corresponds to the minimum nLL
        [fits.min_val, fits.min_idx] = min(model(sub, ses).minNLL);
        %find the corresponding estimated parameters
        fits.p = model(sub, ses).estimatedP(fits.min_idx,:);
        p = fits.p;
        %         p_common(sub,ses) = p(7);

        %% simulate mu_shift to obtain mu_post
        mu_shift = NaN(1, model(sub, ses).expo_num_sim);
        for t   = 1:model(sub, ses).expo_num_sim
            mu_shift(t)   = simulate_mu_shift_MA(model(sub, ses).expo_num_trial, ...
                data(sub, ses).adaptor_soa, fits.p(1),...
                fits.p(7), fits.p(8), fits.p(9), fits.p(10), fits.p(11));
        end

        % mu_post = mu_pre + simulated mean(mu_shift)
        mu2 = fits.p(1) + mean(mu_shift);

        %% plot PMF predictions

        %m ake a finer grid for the timing difference between the auditory and the
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
            P_Vfirst(SOA_finer, mu2, p(5), p(6), p(4)), ...
            P_simultaneous(SOA_finer, mu2, p(5), p(6), p(4)), ...
            P_Afirst(SOA_finer, mu2, p(5), p(6), p(4))};

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
        h1 = xline(p(1),':','LineWidth',2,'Color', repmat(0.3, 1, 3));
        h2 = xline(mu2,'LineWidth',2,'Color', repmat(0.3, 1, 3));

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

    end

    % save figure
    flnm = ['sub' num2str(sub) 'ses' num2str(ses)];
    saveas(gca, fullfile(outDir, flnm),'epsc')

end