% parameter recovery: check analysis_TOJ.m
% This script generates fake data of the ternary TOJ task, then fit fake
% data with the same codes as analysis_TOJ.m
clear all; close  all; clc

%% %%%%%%%%%%%%%%%% simulate fake data %%%%%%%%%%%%%%%%%% %%
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
simTrial                          = 1e3; 
%the number of total trials
nTTrials                           = len_deltaT*simTrial;

%% Model assumptions
%In each trial, the occurrence of an auditory and a visual stimulus leads
%to a measurement regarding the temporal difference between the two stimuli
%in an observer's head. The measurement is denoted as m_{delta_t}. This
%measurement is assumed to have a Gaussian distribution, centered at the
%actual physical temporal difference of that trial and a temporal bias,
%i.e., m_{delta_t} ~ N(t_A - t_V + bias_t, sigma_deltaT).

%The reason we include a bias term is because sounds are normally perceived 
%faster as visual stimuli by ~60ms. In other words, participants perceive 
%an auditory and a visual stimulus as simultaneous when the auditory 
%stimulus is delayed by 60ms.
mu_pre                            = 71.1; 

%the minus sign was put to the model, P_A_precedes_V , P_A_follows_V 
%there is a minus sign here because the x-axis is t_A - t_V, so when the
%physical temporal difference is 0ms, the perceived temporal difference
%will become t_A - t_V + bias_t = -60ms (i.e., the auditory stimulus is 
%perceived as preceding the visual stimulus by 60ms). 

%The measurement distribution has a width of sigma_deltaT
sigma_deltaT                       = 62.01;

% lapse rate
lapse = 0.1;

% The absolute criterion for switching decision from ‘A-coincides-V’ to ‘A-precedes-V’ or ‘A-follows-V’
c                    = 148.76;

% put parameters together
realP = [mu_pre, sigma_deltaT, lapse, c];

%function for creating a measurement distribution
m_deltaT_dist                      = @(x,mu,b,sig) normpdf(x, mu - b, sig);

% function for creating a cumulative Gaussian, the probability of reporting
% 'A-precedes-V','A-follows-V','A-coincides-V' is:
P_A_follows_V        = 1 - normcdf(c, t_diff - mu_pre, sigma_deltaT);
P_A_precedes_V       = normcdf(-c, t_diff - mu_pre, sigma_deltaT);
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
[bool_Vfirst, sim_prob_Vfirst, bool_Afirst, sim_prob_Afirst, bool_simul, sim_prob_simul] ...
    = simTernaryTOJ(len_deltaT, simTrial, Ptilde_A_follows_V, Ptilde_A_precedes_V, Ptilde_A_coincides_V);
nT_Vfirst = sum(bool_Vfirst');
nT_Afirst = sum(bool_Afirst');
nT_simul = sum(bool_simul');
fakeData.nT_Vfirst = sum(bool_Vfirst');
fakeData.nT_Afirst = sum(bool_Afirst');
fakeData.nT_simultaneous = sum(bool_simul');
fakeData.sim_prob_Vfirst = sim_prob_Vfirst;
fakeData.sim_prob_Afirst = sim_prob_Afirst;
fakeData.sim_prob_simul = sim_prob_simul;
fakeData.simTrial = simTrial;
fakeData.s_unique = t_diff;
save('fakeTernaryData','fakeData')

%% %%%%%%%%%%%%%%%% model fiting of fake data %%%%%%%%%%%%%%%%%% %%
%% define the (unscaled/scaled) psychometric function and the cost function
%P is a cumulative gaussian distribution. It takes three inputs
%1. x: the timing difference between the auditory and the visual stimulus
%2. mu: the center of the psychometric function (a.k.a., PSS,
%       perceived subjective simultaneity). Past studies have shown that
%       PSS is not equal to 0, as sounds are normally perceived faster than
%       visual stimuli by ~60ms. In other words, participants perceive
%       an auditory and a visual stimulus as simultaneous when the auditory
%       stimulus is delayed by 60ms.
%3. sig: the measurement noise regarding the temporal difference

%the scaled cumulative gaussian distribution by a lapse rate
P_Afirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous = @(SOA, mu, sig, lambda, c) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

%% fitting the scaled psychometric function by minimizing the cost function
% set parameters
nT_A1st = fakeData.nT_Afirst;% the number of A-first responses for each SOA
nT_V1st = fakeData.nT_Vfirst; % the number of V-first responses for each SOA
nT_simul = fakeData.nT_simultaneous;
s_unique = fakeData.s_unique;

%nLL cost function
nLogL = @(p) -nT_A1st*log(P_Afirst(s_unique, p(1), p(2), p(3), p(4)))' ...
    -nT_V1st*log(P_Vfirst(s_unique, p(1), p(2), p(3), p(4)))'...
    -nT_simul * log(P_simultaneous(s_unique, p(1), p(2), p(3), p(4)))';

%To find the combination of mu, sigma and lapse rate that minimizes the
%negative log likelihood, we will use fmincon.m in MATLAB. To use this
%function, we need to define lower and upper bounds for each parameter
%(i.e., search space) as well as an initial point for MATLAB to start
%searching.
lb      = [ -50, 10, 1e-2, 50];
ub      = [150, 200, 0.1, 150];
init    = rand(1,length(lb)).*(ub-lb) + lb;
%You can also define how many times you want MATLAB to search
options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');

%fmincon returns best-fitting parameters that minimize the cost function as
%well as the corresponding value for the cost function (in this case, the
%negative log likelihood)
[estP, min_NLL] = fmincon(nLogL, init,[],[],[],[],lb,ub,[],options);
%display the best-fitting parameters
disp(estP);
disp(realP);

%% plot fake data and fit data
cMAP                 = [200, 40, 40; 255, 128, 0; 13, 183, 200]./255;

figure
hold on
scatter(s_unique, fakeData.sim_prob_Afirst,'MarkerFaceColor', cMAP(3,:),...
    'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.5)
scatter(s_unique, fakeData.sim_prob_Vfirst,'MarkerFaceColor', cMAP(1,:),...
    'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.5)
scatter(s_unique, fakeData.sim_prob_simul,'MarkerFaceColor', cMAP(2,:),...
    'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.5)

t_diff_finer = linspace(s_unique(1), s_unique(end), 1000);
P_Vfirst_fit = P_Vfirst(t_diff_finer, estP(1), estP(2), estP(3), estP(4));
P_Afirst_fit = P_Afirst(t_diff_finer, estP(1), estP(2), estP(3), estP(4));
P_simultaneous_fit = P_simultaneous(t_diff_finer, estP(1), estP(2), estP(3), estP(4));

plot(t_diff_finer, P_Afirst_fit, 'Color', cMAP(3,:), ...
    'lineWidth', 3); hold on;
plot(t_diff_finer, P_Vfirst_fit, 'Color', cMAP(1,:), ...
    'lineWidth', 3); hold on;
plot(t_diff_finer, P_simultaneous_fit, 'Color', cMAP(2,:), ...
    'lineWidth', 3); hold on;
h1= xline(estP(1),'LineWidth',2);
h2 = xline(estP(1) - estP(4),'--','LineWidth',2);
h3 = xline(estP(1) + estP(4),'--','LineWidth',2);
legend({'$P(A-precedes-V)$','$P(A-follows-V)$',...
    '$P(A-coincides-V)$'},'Location','bestoutside', 'Interpreter','Latex');
xlim([-400 400])
title(['model fitting of fake data with simulated trials of ' num2str(fakeData.simTrial)])

fignm = ['sim_pretest'];
saveas(gca,fignm,'png')

