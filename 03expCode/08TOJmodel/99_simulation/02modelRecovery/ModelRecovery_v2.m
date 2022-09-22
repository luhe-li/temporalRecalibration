% model recovery

% 2022/03/07

% This script generates fake dataset of pre and post test of a ternary
% task, fit the fake datasets into a null-model (4 free parameters, mu is
% the same) and a real model (5 free parameters, mu is different from pre
% and post test). It then compare models by AIC and BIC.

% This process is repeated 100 times. We then summarize results in a
% confusion matrix quantifies the probability that each model is the best
% fit to data generated from the other models. In a perfect world, the
% confusion matrix will be the identity matrix, but in practice, this is
% not always the case.

clear all; close all; clc; rng(2);


%% initiate confusion matrix
CM = zeros(2);
counts = 100;

%% create movies
%create an empty struct
numFramesPerL = 20;
mov(1:counts*numFramesPerL)=struct('cdata',[],'colormap',[]);
v = VideoWriter('modelRecoveryNumTrial_15.avi');
open(v)

for count = 1:counts
    count

    figure(1); clf;
    FM = round(100*CM/sum(CM(1,:)))/100;
    t = imageTextMatrix(FM);
    set(t(FM'<0.3), 'color', 'w')
    hold on;
    [l1, l2] = addFacetLines(CM);
    set(t, 'fontsize', 22)
    title(['count = ' num2str(count)]);
    set(gca, 'xtick', [1:2], 'ytick', [1:2], 'fontsize', 28, ...
        'xaxislocation', 'top', 'tickdir', 'out')
    xlabel('fit model')
    ylabel('simulated model')

    pause(0.1);
    mov(((count-1)*numFramesPerL+1):(count*numFramesPerL)) = getframe(gcf);

    %% simulate fake data

    % define fake data parameters
    muSets = {[70.02, 40.02], [40.02, 40.02]};
    % initate output for model comparison
    numM     = 2; %number of models
    [min_NLL, L_test] = deal(NaN(2,numM));
    estP              = cell(2,numM);
    [AIC, BIC]        = deal(NaN(2, numM));

    for d = 1:2

        %the levels of SOA
        s_unique = [-600,-400:50:400,600];
        %the total number of levels
        len_deltaT = length(s_unique);
        %the number of trial for each audiovisual pair
%         nTrials = [5:5:20];
%         t = length(nTrials);
        simTrial = 15;
        %the number of total trials
        nTTrials = len_deltaT*simTrial;
        % set true parameters.
        % pre and post tests have two different mus
        muPrePost   = muSets{d};
        %The measurement distribution has a width of sigma_deltaT
        sigma_deltaT = 69.07;
        % lapse rate
        lapse = 0.06;
        % The absolute criterion for switching decision from ‘A-coincides-V’ to ‘A-precedes-V’ or ‘A-follows-V’
        c = 147.76;
        % put parameters together
        realP = [muPrePost , sigma_deltaT, lapse, c];

        for s  = 1:2 % for pre and post test
            mu  = muPrePost(s);
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
            [bool_Vfirst{s}, sim_prob_Vfirst{s}, bool_Afirst{s}, ...
                sim_prob_Afirst{s}, bool_simul{s}, sim_prob_simul{s}] ...
                =simTernaryTOJ(len_deltaT, simTrial, Ptilde_A_follows_V, Ptilde_A_precedes_V, Ptilde_A_coincides_V);
            % number of each response at each level
            nT_Afirst{s} = sum(bool_Afirst{s}'); % the number of V-first responses for each SOA
            nT_Vfirst{s} = sum(bool_Vfirst{s}'); % the number of A-first responses for each SOA
            nT_simul{s} = sum(bool_simul{s}');
            % nT_Vfirst + nT_Afirst + nT_simul should add up to trial number
            % organize raw simulated response into a single matrix (Vfirst = 1,
            % simultaneous = 2, Afirst = 3)
            r_org{s} = bool_Vfirst{s} + bool_simul{s} * 2 + bool_Afirst{s}*3;
        end
        %     fakeName = {'MR_fakeData_realModel', 'MR_fakeData_nullModel'};
        %     save(fakeName{d})

        %% fit real model data to both models
        % M1(real model). Assumes that mu changes from pre to post test. It has 5
        % free parameters: mu_pre, mu_post, sigma, lambda, criterion.

        % M2(null model). Assumes that mu is the same between pre and post test. It
        % has 4 free parameters: mu, sigma, lambda, criterion.

        % define models
        numM     = 2; %number of models
        numP     = [5, 4]; %number of free parameters for M1, M2 respectively
        P_Afirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
        P_Vfirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
        P_simultaneous = @(SOA, mu, sig, lambda, c) ...
            1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
            - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

        %define upper and lower bounds
        lb       = {[ -150, -150, 10, 1e-2, 50], [ -150, 10, 1e-2, 50]};
        ub       = {[150, 150, 200, 0.06, 250], [150, 200, 0.06, 250]};
        init_fun = @(a,b) rand(1,length(a)).*(b-a) + a;
        options  = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');

        %negative log likelihood
        nLogL{1} = @(p) -nT_Afirst{1} * log(P_Afirst(s_unique, p(1), p(3), p(4), p(5)))' ...
            -nT_Vfirst{1} *log(P_Vfirst(s_unique, p(1), p(3), p(4), p(5)))'...
            -nT_simul{1} * log(P_simultaneous(s_unique, p(1), p(3), p(4), p(5)))'...
            -nT_Afirst{2} * log(P_Afirst(s_unique, p(2), p(3), p(4), p(5)))' ...
            -nT_Vfirst{2} * log(P_Vfirst(s_unique, p(2), p(3), p(4), p(5)))'...
            -nT_simul{2}  * log(P_simultaneous(s_unique,  p(2), p(3), p(4), p(5)))';% M1
        nLogL{2} = @(p)  -nT_Afirst{1} * log(P_Afirst(s_unique, p(1), p(2), p(3), p(4)))' ...
            -nT_Vfirst{1} *log(P_Vfirst(s_unique, p(1), p(2), p(3), p(4)))'...
            -nT_simul{1} * log(P_simultaneous(s_unique, p(1), p(2), p(3), p(4)))'...
            -nT_Afirst{2} * log(P_Afirst(s_unique, p(1), p(2), p(3), p(4)))' ...
            -nT_Vfirst{2}*log(P_Vfirst(s_unique, p(1), p(2), p(3), p(4)))'...
            -nT_simul{2} * log(P_simultaneous(s_unique,  p(1), p(2), p(3), p(4)))'; % M2

        %loop through the two models
        for m = 1:numM
            %the initial point for matlab to start searching
            init = init_fun(lb{m}, ub{m});
            %use fmincon.m to fit
            [estP{d,m}, min_NLL(d,m)] = fmincon(nLogL{m}, init,[],[],[],[],...
                lb{m}, ub{m},[],options);
            %compute the AIC/BIC
            AIC(d, m) = 2*min_NLL(d, m) + 2*numP(m);
            BIC(d, m) = 2*min_NLL(d,m) + numP(m)*log(nTTrials);
        end

        [M iBEST] = min(AIC(d,:));
        BEST = AIC(d,:) == M;
        BEST = BEST / sum(BEST);
        CM(d,:) = CM(d,:) + BEST;

        % relative likelihood
        % pAIC = exp((min(AIC) - AIC)/2);

    end
end
writeVideo(v,mov)
close(v)
% set(gca,'FontSize',15);
% saveas(gca,'MR','png')
