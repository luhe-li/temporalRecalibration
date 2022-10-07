% parameter recovery: determine trial number

% This script generates several fake datasets of different number of
% trials of a ternary TOJ task for pre and post test, bootstrap for 1000
% trials to obtain the confidence interval of parameters, then repeat the
% whole process for 20 times (20 fake datasets)

clear all; close all; clc; rng(2);

%% generate fake data

%the levels of SOA
s_unique = [-600, -400:50:400, 600];
%the total number of levels
len_deltaT = length(s_unique);
%the number of trial for each audiovisual pair
nTrials = [21:24];
% initiate figure for each simulated numTrial
figure; hold on
% the number of simulated datasets
nDatasets = 20;

for t = 1:length(nTrials)
    f(t) = figure;
    g(t) = figure;
    for nData = 1:nDatasets
        simTrial = nTrials(t);
        %the number of total trials
        ntrials = len_deltaT*simTrial;
        % set true parameters
        muPrePost   = [70.02, 40.02];
        for s  = 1:2 % for pre and post test
            mu  = muPrePost(s);
            %The measurement distribution has a width of sigma_deltaT
            sigma_deltaT = 69.07;
            % lapse rate
            lapse = 0.06;
            % The absolute criterion for switching decision from ‘A-coincides-V’ to ‘A-precedes-V’ or ‘A-follows-V’
            c = 147.76;
            % put parameters together
            realP = [muPrePost , sigma_deltaT, lapse, c];
            %function for creating a measurement distribution
            m_deltaT_dist = @(x,mu,b,sig) normpdf(x, mu - b, sig);
            % function for creating a cumulative Gausian, the probability of reporting
            % 'A-precedes-V','A-follows-V','A-coincides-V' is:
            P_A_follows_V = 1 - normcdf(c, s_unique - mu, sigma_deltaT);
            P_A_precedes_V = normcdf(-c, s_unique - mu, sigma_deltaT);
            P_A_coincides_V = 1 - P_A_precedes_V - P_A_follows_V;
            %the psychometric function after taking lapses into account
            Ptilde = @(P, lambda) lambda/3 + (1-lambda).*P;
            % on \lambda of the trials they ignore the stimulus and gues, with equal
            % numbers of gueses for each response
            Ptilde_A_follows_V = Ptilde(P_A_follows_V, lapse);
            Ptilde_A_precedes_V = Ptilde(P_A_precedes_V, lapse);
            Ptilde_A_coincides_V = Ptilde(P_A_coincides_V, lapse);
            % simulate fake data
            [bool_Vfirst{t,s}, sim_prob_Vfirst{t,s}, bool_Afirst{t,s}, ...
                sim_prob_Afirst{t,s}, bool_simul{t,s}, sim_prob_simul{t,s}] ...
                =simTernaryTOJ(len_deltaT, simTrial, Ptilde_A_follows_V, Ptilde_A_precedes_V, Ptilde_A_coincides_V);
            % nT_Vfirst + nT_Afirst + nT_simul should add up to trial number
            % organize raw simulated response into a single matrix (Vfirst = 1,
            % simultaneous = 2, Afirst = 3)
            r_org{t,s} = bool_Vfirst{t,s} + bool_simul{t,s} * 2 + bool_Afirst{t,s}*3;
        end
    
        %% bootstrap fake data of pre and post tests
        P_Afirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
        P_Vfirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
        P_simultaneous = @(SOA, mu, sig, lambda, c) ...
            1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
            - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

        % define lower and upper bounds for each parameter (i.e.,
        % search space) as well as an initial point for MATLAB to start
        % searching. mu_pre, mu_post, sigma, lambda, criterion
        lb = [ -150, -150, 10, 1e-2, 50];
        ub = [150, 150, 200, 0.06, 250];
        init = rand(1,length(lb)).*(ub-lb) + lb;
        options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');

        numBtst = 1e3;
        [estP_btst{t}, minNLL_btst{t}, CI_btst_lb{t}, CI_btst_ub{t}] = ...
            BootstrapPrePosttestTOJ(s_unique,r_org{t,1},r_org{t,2}, simTrial, numBtst, P_Afirst, P_Vfirst, P_simultaneous, lb, ub, options);
        save('pr_pSummary','estP_btst', 'minNLL_btst', 'CI_btst_lb','CI_btst_ub','r_org')

        %% Plot the best-fitting parameters given each bootstrapped dataset
        %parameter names
        valName = {'\mu_{pre}','\mu_{post}', '\sigma', '\lambda', 'criterion'};
        numP = 5;
        cMAP = [200, 40, 40; 255, 128, 0; 13, 183, 200]./255;
        figure; set(gcf,'position',[0,0,400,800])
        fignm2 = ['simulated for ' num2str(simTrial) ' trials'];
        sgtitle(fignm2)
        hold on;
        for i = 1:numP
            subplot(numP,1,i)
            histogram(estP_btst{t}(:,i),'FaceColor', cMAP(2,:), ...
                'FaceAlpha', 0.5, 'EdgeColor',cMAP(2,:)); hold on
            plot([realP(i), realP(i)], [0, numBtst*0.3],'r--', 'lineWidth',3); hold on
            plot([CI_btst_lb{t}(:,i), CI_btst_ub{t}(:,i)],[5, 5],'-k','LineWidth', 2); hold on
            xlim([lb(i), ub(i)]); ylim([0, numBtst*0.3]);
            xticks(sort(unique([linspace(lb(i), ub(i),5), realP(i)])));
            title(valName{i}); box off
            set(gca,'FontSize',15);
        end
        saveas(gca,fignm2,'png')

        %% plot mu_change against each dataset (x-axis), each numTrial (figure)
        figure(f(t)); hold on
        mu_diff = estP_btst{t}(:,1) - estP_btst{t}(:,2);
        mu_diffS{t, nData} = mu_diff;
        [counts{nData} edges{nData}] = histcounts(mu_diff,10);
        imagesc([nData, nData],[min(edges{nData}), max(edges{nData})],counts{nData}')
        % look better
        set(gca,'FontSize',15);
        sumfig = ['simulated for ' num2str(simTrial) ' trials'];
        title(sumfig)
        ylabel(['\Delta_{\mu} bootstrapped  for ' num2str(numBtst) ' trials'])
        xlabel('data sets')
        xlim([0, nDatasets+1])
        ylim([-100 150])
        l1 = yline(0,'--','LineWidth',2);
        cb = colorbar;
        cb.Label.String = 'number of bootstrapped trials';

        %% plot mu_change CI against each dataset (x-axis), each numTrial (figure)
        % compute mean
        m = mean(mu_diff);
        % compute 95%CI
        n_sorted = sort(mu_diff);
        %compute how many entries the vector has
        lenN     = length(mu_diff);
        %lower bound
        CI_ub    = n_sorted(ceil(lenN*0.975));
        %upper bound
        CI_lb    = n_sorted(floor(lenN*0.025));
        % plot
        l2 = plot([nData, nData],[CI_lb, CI_ub],'k-','LineWidth',1.5);
%         legend([l1, l2],{'ground truth \Delta_{\mu}','95 CI'},'Location','northwest')
        legend(l2, '95 CI','Location','northwest')
        saveas(gca,sumfig,'png')
        
    end
end
