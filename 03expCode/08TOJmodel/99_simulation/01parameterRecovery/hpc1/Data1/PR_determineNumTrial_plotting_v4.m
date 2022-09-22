
clear all; clc; close all;

folder = pwd;
filePattern = fullfile(folder, 'PR_determineNumTrials_04-May-2022 *.mat');
files = dir(filePattern);

% simulated number of trial for each SOA, indexed by t
simTrial = 10:5:30;
t_len = length(simTrial);

% set a counter by t to index dataset for plotting; counts from 1 to 500
% for each t
c = ones(1, 5);

% init ub and lb of mu_diff, size = ntrial x ndataset
LBs = NaN(t_len, 500);
UBs = NaN(t_len, 500);

% initiate figures
for t = 1:5
    f(t) = figure;
end

for i = 1:10 % for each mat files

    %%%%%% load each mat file %%%%%%
    filenm = files(i).name;
    load(filenm);

    %%%%%% organize simulation data %%%%%%
    % Data{1,4}.estP_btst: 1x5 cells, each is simulated for nTrials = 10:5:30
    % Data{1,4}.estP_btst{1,1}: 1x50 cells, each is bootstrapped parameters for
    % a single dataset

    % inside loops, plot mu difference for each data set
    D = Data{1,4};

    for t = 1:t_len % for each nTrial

        for j = 1:50 % for each dataset
            figure(f(t)); hold on

            % overwrite estP_btst for individual dataset
            i_estP_btst = D.estP_btst{1,t}{1,j};
            mu_diff = i_estP_btst(:,2) - i_estP_btst(:,1);

            % plot histogram of mu_diff
            [counts, edges] = histcounts(mu_diff,10);
            imagesc([c(t), c(t)],[min(edges), max(edges)],counts')

            % add 95% CI
            [lb, ub] = get95CI(mu_diff);
            plot([c(t), c(t)],[lb, ub],'k-','LineWidth',1.5);

            % save lb and ub for power calculation
            LBs(t, c(t)) = lb; UBs(t, c(t)) = ub;

            % update counter of dataset
            c(t) = c(t)+1;
        end
    end
end


%% calculate power
power = sum((LBs >= 0),2)./500;

%% look better
for t = 1:t_len
    figure(f(t))
    set(gca,'FontSize',15);
    set(f(t),'Position',[10 10 2000 600])
    ylabel(['\Delta_{\mu} bootstrapped for 1000 trials'])
    xlabel('datasets')
    xlim([0, 500+1])
    ylim([-300 300])
    yline(0,'-w','LineWidth',2);
    cb = colorbar;
    cb.Label.String = 'number of bootstrapped trials';
    title(['Simulate ' num2str(simTrial(t)) ' trials per SOA, power = ' num2str(power(t))])
end