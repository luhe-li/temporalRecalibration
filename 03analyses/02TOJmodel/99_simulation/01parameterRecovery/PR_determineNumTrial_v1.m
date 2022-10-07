% parameter recovery: trialnumber
% This script generates fake data of the ternary TOJ task, then fit fake
% data with the same codes as analysis_TOJ.m, boostrap for 1000 trials to
% obtain the confidence interval of parameters

clear all; close all; clc; rng(2);

%% Setting up experimental info
%In each trial, participants are presented with an auditory and a visual
%stimulus with a temporal discrepancy between them. The discrepancy can
%have various levels, ranging from -350 to 350 ms with an increment of 50
%ms. Positive values represent the visual stimulus coming before the
%auditory stimulus; negative values represent the auditory stimulus coming
%first. After stimulus presentation, participants are asked to report
%whether they judge the temporal order, i.e., report which stimulus comes
%first (V, A or simutaneous). Each temporal discrepancy (a.k.a. stimulus onset
%asynchrony; SOA) is tested multiple times.

%let's first define some experimental info
%the levels of SOA
t_diff                             = -350:50:350;
%the total number of levels
len_deltaT                         = length(t_diff);
%the number of trial for each audiovisual pair
% nTrials= [2:2:20];
nTrials= [8:2:30];

for tt                          = 1:length(nTrials) % look through possible nTrials
    
    simTrial = nTrials(tt);
    
    %the number of total trials
    nTTrials                           = len_deltaT*simTrial;
    
    % initiate figures
    f1 = figure;
    f2 = figure;
    
   %% %%%%%%%%%%%%%%%% simulate fake data %%%%%%%%%%%%%%%%%% %%
    muPrePost                            = [70.02, 40.02];
    for ii = 1:2 % loop through pre_mu to post_mu
        
        mu = muPrePost(ii);
        %the minus sign was put to the model, P_A_precedes_V , P_A_follows_V
        %there is a minus sign here because the x-axis is t_A - t_V, so when the
        %physical temporal difference is 0ms, the perceived temporal difference
        %will become t_A - t_V + bias_t = -60ms (i.e., the auditory stimulus is
        %perceived as preceding the visual stimulus by 60ms).
        
        %The measurement distribution has a width of sigma_deltaT
        sigma_deltaT                       = 69.07;
        
        % lapse rate
        lapse = 0.06;
        
        % The absolute criterion for switching decision from ‘A-coincides-V’ to ‘A-precedes-V’ or ‘A-follows-V’
        c                    = 147.76;
        
        % put parameters together
        realP = [mu, sigma_deltaT, lapse, c];
        
        %function for creating a measurement distribution
        m_deltaT_dist                      = @(x,mu,b,sig) normpdf(x, mu - b, sig);
        
        % function for creating a cumulative Gaussian, the probability of reporting
        % 'A-precedes-V','A-follows-V','A-coincides-V' is:
        P_A_follows_V        = 1 - normcdf(c, t_diff - mu, sigma_deltaT);
        P_A_precedes_V       = normcdf(-c, t_diff - mu, sigma_deltaT);
        P_A_coincides_V      = 1 - P_A_precedes_V - P_A_follows_V;
        
        %% Add lapse
        %the psychometric function after taking lapses into account
        Ptilde               = @(P, lambda) lambda/3 + (1-lambda).*P;
        
        % on \lambda of the trials they ignore the stimulus and guess, with equal
        % numbers of guesses for each response
        Ptilde_A_follows_V   = Ptilde(P_A_follows_V, lapse);
        Ptilde_A_precedes_V  = Ptilde(P_A_precedes_V, lapse);
        Ptilde_A_coincides_V = Ptilde(P_A_coincides_V, lapse);
        
        %% simulate fake data
        [bool_Vfirst, sim_prob_Vfirst, bool_Afirst, sim_prob_Afirst, bool_simul, ...
            sim_prob_simul] = simTernaryTOJ(len_deltaT, simTrial, Ptilde_A_follows_V, Ptilde_A_precedes_V, Ptilde_A_coincides_V);
        nT_Vfirst = sum(bool_Vfirst'); % the number of A-first responses for each SOA
        nT_Afirst = sum(bool_Afirst'); % the number of V-first responses for each SOA
        nT_simul = sum(bool_simul');
        % nT_Vfirst + nT_Afirst + nT_simul should add up to trial number
        % organize raw simulated response into a single matrix (Vfirst = 1,
        % simultaneous = 2, Afirst = 3)
        r_org =bool_Vfirst + bool_simul * 2 + bool_Afirst*3;
        
        %% %%%%%%%%%%%%%%%% bootsrap fake data %%%%%%%%%%%%%%%%%% %%
        %% define the (unscaled/scaled) psychometric function and the cost function
        P_Afirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
        P_Vfirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
        P_simultaneous = @(SOA, mu, sig, lambda, c) ...
            1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
            - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));
       
          %Bootstrap uses random sampling with replacement (i.e., say that you have
        %[1,2,3,4,5] in a magic bag. You randomly draw one number from the bag
        %each time, and then put it back). This means as you repeatedly sample from
        %the bag, each number can be selected more/less than once. Bootstrap is
        %often used when you want to get a confidence interval on estimated
        %parameters.
        
        s_unique = t_diff;
        numBtst = 1e2;
        [estP_btst{tt,ii}, minNLL_btst{tt,ii}, CI_btst_lb{tt,ii}, CI_btst_ub{tt,ii}] = BootstrapTernaryTOJ(s_unique,...
            r_org, simTrial, numBtst, P_Afirst, P_Vfirst, P_simultaneous, lb, ub, options);
        save('pSummary2-20','estP_btst', 'minNLL_btst', 'CI_btst_lb','CI_btst_ub')
        
         %% optinal: fitting faka data, plot fake data and fitting curve  
%          %nLL cost function
%          nLogL = @(p) -nT_Afirst*log(P_Afirst(s_unique, p(1), p(2), p(3), p(4)))' ...
%              -nT_Vfirst*log(P_Vfirst(s_unique, p(1), p(2), p(3), p(4)))'...
%              -nT_simul * log(P_simultaneous(s_unique, p(1), p(2), p(3), p(4)))';
%          
%          %To find the combination of mu, sigma and lapse rate that minimizes the
%          %negative log likelihood, we will use fmincon.m in MATLAB. To use this
%          %function, we need to define lower and upper bounds for each parameter
%          %(i.e., search space) as well as an initial point for MATLAB to start
%          %searching.
%          lb      = [ -150, 10, 1e-2, 50];
%          ub      = [150, 200, 0.06, 250];
%          init    = rand(1,length(lb)).*(ub-lb) + lb;
%          %You can also define how many times you want MATLAB to search
%          options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');
%          
%          %fmincon returns best-fitting parameters that minimize the cost function as
%          %well as the corresponding value for the cost function (in this case, the
%          %negative log likelihood)
%          [estP, min_NLL] = fmincon(nLogL, init,[],[],[],[],lb,ub,[],options);
%          %display the best-fitting parameters
%          disp(estP);
%          disp(realP);
%          
%          cMAP                 = [200, 40, 40; 255, 128, 0; 13, 183, 200]./255;
%          
%          % figure
%          figure(f1)
%          subplot(1,2,ii)
%          hold on
%          scatter(s_unique, sim_prob_Afirst,'MarkerFaceColor', cMAP(3,:),...
%              'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.5)
%          scatter(s_unique, sim_prob_Vfirst,'MarkerFaceColor', cMAP(1,:),...
%              'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.5)
%          scatter(s_unique, sim_prob_simul,'MarkerFaceColor', cMAP(2,:),...
%              'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.5)
%          
%          t_diff_finer = linspace(s_unique(1), s_unique(end), 1000);
%          P_Vfirst_fit = P_Vfirst(t_diff_finer, estP(1), estP(2), estP(3), estP(4));
%          P_Afirst_fit = P_Afirst(t_diff_finer, estP(1), estP(2), estP(3), estP(4));
%          P_simultaneous_fit = P_simultaneous(t_diff_finer, estP(1), estP(2), estP(3), estP(4));
%          
%          plot(t_diff_finer, P_Afirst_fit, 'Color', cMAP(3,:), ...
%              'lineWidth', 3); hold on;
%          plot(t_diff_finer, P_Vfirst_fit, 'Color', cMAP(1,:), ...
%              'lineWidth', 3); hold on;
%          plot(t_diff_finer, P_simultaneous_fit, 'Color', cMAP(2,:), ...
%              'lineWidth', 3); hold on;
%          h1= xline(estP(1),'LineWidth',2);
%          h2 = xline(estP(1) - estP(4),'--','LineWidth',2);
%          h3 = xline(estP(1) + estP(4),'--','LineWidth',2);
%          %     legend({'$P(A-precedes-V)$','$P(A-follows-V)$',...
%          %         '$P(A-coincides-V)$'},'Location','bestoutside', 'Interpreter','Latex');
%          xlim([-400 400])
%          title(['model fitting of fake data, trial number = ' num2str(simTrial)])
%          
%          fignm1 = ['sim_test_nTrial' num2str(simTrial)];
%          saveas(gca,fignm1,'png')
%          
      
        %% Plot the best-fitting parameters given each bootstrapped dataset
        %parameter names
        valName = {'\mu', '\sigma', '\lambda', 'criterion'};
        numP = 4;
        figure(f2);        hold on
        
        for i = 1:numP
            subplot(numP,1,i)
            l{ii} = histogram(estP_btst{tt,ii}(:,i),'FaceColor', cMAP(ii+1,:), ...
                'FaceAlpha', 0.5, 'EdgeColor',cMAP(ii+1,:)); hold on
            l{3} = plot([realP(i), realP(i)], [0, numBtst*0.3],'r--', 'lineWidth',3); hold on
            plot([CI_btst_lb{tt,ii}(:,i), CI_btst_ub{tt,ii}(:,i)],[290,290],'--','Color',cMAP(ii+1,:),'LineWidth', 2); hold on
            xlim([lb(i), ub(i)]); ylim([0, numBtst*0.3]);
            xticks(sort(unique([linspace(lb(i), ub(i),5), realP(i)])));
            title(valName{i}); box off
            set(gca,'FontSize',15);
            if i == 1
                title('Estimated value given bootstrapped datasets')
            end
        end
        
        if ii == 2
            legend([l{1}, l{2}, l{3}],{'Pre-test','Post-test','True value'},...
                'Location', 'northwest');
            legend boxoff;
        end
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.7]);
        save
        
        fignm2 = ['btsp_parameters_nTrial' num2str(simTrial)];
        saveas(gca,fignm2,'png')
        
    end
end

%% %%%%%%%%%%%%%%%% summarize btst by trial %%%%%%%%%%%%%%%%%% %%
clear all; close all; clc;
load('pSummary2-20.mat')
% for each parameter
for pp = 1:4
    figure; set(gca,'FontSize',15); hold on
    ts = [10:10:100];
    valName = {'\mu_{pre} - \mu_{post}', '\sigma_{pre} - \sigma_{post}',...
        '\lambda_{pre} - \lambda_{post}', 'criterion_{pre} - criterion_{post}'};
    realPchange = [30, 0, 0, 0];
    for tt = 1:length(ts)
        t = ts(tt);
     
        diff = sort(estP_btst{tt,1}(:,pp)) - sort(estP_btst{tt,2}(:,pp));
        m = mean(diff);
        % compute 95%CI
        n_sorted = sort(diff);
        %compute how many entries the vector has
        lenN     = length(diff);
        %lower bound
        CI_ub    = n_sorted(ceil(lenN*0.975));
        %upper bound
        CI_lb    = n_sorted(floor(lenN*0.025));
        % plot
        eb = (CI_ub + CI_lb)./2;
        plot(t,m,'kx','LineWidth',1)
        errorbar(t,m,eb,'k','LineWidth',1)
    end
    yline(realPchange(pp))
    ylabel(valName{pp})
    xlim([ts(1)-5, ts(end)+5])
    xlabel('trial number')
    title('bootstrap for 1e3 trials')
     fignm3 = ['btsp_parameters_summary1_p' num2str(pp)];
        saveas(gca,fignm3,'png')
end
% 
% 
% clear all; close all; clc;
% load('pSummary2-20.mat')
% % for each parameter
% for pp = 1:4
%     figure; set(gca,'FontSize',15); hold on
%     ts = [2:2:20];
%     valName = {'\mu_{pre} - \mu_{post}', '\sigma_{pre} - \sigma_{post}',...
%         '\lambda_{pre} - \lambda_{post}', 'criterion_{pre} - criterion_{post}'};
%     realPchange = [30, 0, 0, 0];
%     for tt = 1:length(ts)
%         t = ts(tt);
%      
%         diff = sort(estP_btst{tt,1}(:,pp)) - sort(estP_btst{tt,2}(:,pp));
%         m = mean(diff);
%         % compute 95%CI
%         n_sorted = sort(diff);
%         %compute how many entries the vector has
%         lenN     = length(diff);
%         %lower bound
%         CI_ub    = n_sorted(ceil(lenN*0.975));
%         %upper bound
%         CI_lb    = n_sorted(floor(lenN*0.025));
%         % plot
%         eb = (CI_ub + CI_lb)./2;
%         plot(t,m,'kx','LineWidth',1)
%         errorbar(t,m,eb,'k','LineWidth',1)
%     end
%     yline(realPchange(pp))
%     ylabel(valName{pp})
%     xlim([ts(1)-2, ts(end)+2])
%     xlabel('trial number')
%     title('bootstrap for 1e3 trials')
%      fignm3 = ['btsp_parameters_summary1_p' num2str(pp)];
%         saveas(gca,fignm3,'png')
% end
