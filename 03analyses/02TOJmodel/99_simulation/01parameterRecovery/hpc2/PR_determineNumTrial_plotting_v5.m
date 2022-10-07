
clear all; clc; close all;


folder = pwd;
filePattern = fullfile(folder, 'PR_determineNumTrials_test*.mat');
files = dir(filePattern);

for dd = 1:6
    filenm = files(dd).name;
    load(filenm);

    % inside loops, plot mu difference for each data set
    D = Data{1,4};

    % simulated number of trial for each SOA, indexed by t
    simTrial = Data{1,3}.nTrials;
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

    %%%%%% organize simulation data %%%%%%
    % Data{1,4}.estP_btst: 1x5 cells, each is simulated for nTrials = 10:5:30
    % Data{1,4}.estP_btst{1,1}: 1x50 cells, each is bootstrapped parameters for
    % a single dataset

    for t = 1:t_len % for each nTrial

        for j = 1:500 % for each dataset
            figure(f(t)); hold on

            % overwrite estP_btst for individual dataset
            i_estP_btst = D.estP_btst{1,j}{1,t};
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


    %% calculate power
    load('summary_power.mat')
    power = sum((LBs >= 0),2)./500;
    all_power(dd,:) = power';
    save('summary_power', 'all_power')

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
        fn = ['Test6 Simulate ' num2str(simTrial(t)) ' trials per SOA, power = ' num2str(power(t)*1000)];
        title(fn)
        saveas(gcf, fn, 'png')
    end

end

%% summary plot
clear all; clc; close all;
load('summary_power.mat')
  
testedTrials = [14 17 20 23 26];
figure; 
set(gcf,'Position',[0,0,900,400])
set(gca, 'LineWidth',1.5)
set(gca,'FontSize',15);hold on
plot(testedTrials, all_power,'-o','LineWidth',1.5)
xticks(testedTrials)
ylim([0.5 1])
yline(0.8,'--k','LineWidth',1.5)
xlabel('tested trials per SOA')
ylabel('power')

legend({'[-500,-200,-100,-50:16.66:50,100,200,500],len=13', ...
    '[-500,-300,-200,-100:33.33:100,200,300,500],len=13', ...
    '[-400:50:400],len=17', ...
    '[-500,-300,-200,-100,-50:16.66:50,100,200,300,500],len=15', ...
    '[-500,-350:50:350,500],len=17', ...
    '[-500,-300:50:300,500],len=15'},'Location','bestoutside')

saveas(gca,'summary_power','png')