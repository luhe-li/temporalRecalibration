% this script is a individual version of p2
% fig1-3: delta_mu/delta_sigma/delta_criterion by adaptor SOA for each individual
clear all; close all; clc;

% set data path
currentDir            = pwd;
exptDir               = currentDir(1:regexp(pwd,'03expCode')-1);
outDir                = [currentDir '/p6_figures'];
addpath(genpath([exptDir '03ExpCode/05functions']));
addpath(genpath([exptDir '03ExpCode/06helperFunctions']));
addpath(genpath([exptDir '03ExpCode/01pretest/data']));
addpath(genpath([exptDir '03ExpCode/04posttest/data']));

%% extract all parameters

load('best_para_2.mat')

% define subjects and sessions to use
all_sub               = [1:10];
all_sess              = 1:9; % complete sessions
% avail_sessions = [7,5,repmat(9,1,8)]; % incomplete sessions

% initiate
best_para_            = cell(all_sub(end), all_sess(end));

% sort by adaptor soa order
for s                 = all_sub
    [B, I]                = sort(idx_adaptor_soa(s,:));
    [best_para_{s,:}]     = deal(best_para{s,I});
end

% initiate
% [delta_mu, delta_sigma, delta_c]       = deal(NaN(all_sub(end), all_sess(end)));

% initiate
delta_para       = cell(1,3);

for i                 = all_sub
    for j                 = all_sess
        % delta_mu = post_mu - pre_mu
        delta_para{1}(i,j)  = best_para_{i,j}(2) - best_para_{i,j}(1);
        % delta_sigma = post_sigma - pre_sigma
        delta_para{2}(i,j)  = best_para_{i,j}(4) - best_para_{i,j}(3);
        % delta_c = post_c - pre_c
        delta_para{3}(i,j)  = best_para_{i,j}(6) - best_para_{i,j}(5);
    end
end


%% plot mu
para = {'\Delta_{\mu}','\Delta_{\sigma}','\Delta_{criterion}'};

for p = 1:3
    figure; hold on
    set(gcf, 'Position', get(0, 'Screensize'));
    sgtitle(para{p})

    for i = all_sub
        plots{i} = subplot(3,4,i); hold on
        set(gca, 'LineWidth', 1, 'FontSize', 10)

        % plot
        plot(all_sess, delta_para{p}(i,:),'-o','LineWidth',1.5);
        yline(0,'--');
        yline(mean(delta_para{p}(i,:)),'-','LineWidth',1.5);% mean across sessions

        % look better
        title(['sub' num2str(i)])
        tick_adaptor_soa = [-700, -300:100:300, 700];
        xticks(all_sess)
        xticklabels(strsplit(num2str(tick_adaptor_soa)))
        xlabel('adaptor soa')

    end
    linkaxes([plots{all_sub}],'xy')

    % save figure
    flnm = ['para' num2str(p) '_all_sess'];
    saveas(gca, fullfile(outDir, flnm),'epsc')

end