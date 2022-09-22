% a04_model_recovery_v1.m

% Parameter recovery for 5 models of the TOJ responses, also serves as
% debugging if you set ExpInfo.nTrials to be over 1e3. For each model, we
% used fixed parameters to simulate fake data, and use five models to
% fit. We return both the best fitting parameters for each model to see if
% ground-truth parameters are recovered, and return delta_AIC to see if
% true model is chosen as the best fitting model.

clear all; close all; clc; rng(1); %parpool(20);

%% initiate confusion matrix

nModel = 5;
CM = zeros(nModel);
counts = 100;

%% experimental parameters

ExpInfo.ifi              = 1000/60;
ExpInfo.SOA              = [-0.5000, -0.3000, -0.2500, -0.2000, -0.1500, -0.1000,...
    -0.0500, 0, 0.0500, 0.1000, 0.1500, 0.2000, 0.2500, 0.3000, 0.5000]*1000; % in ms
ExpInfo.lenS             = length(ExpInfo.SOA);
ExpInfo.nTrials          = 1000; % num of trials per SOA

%% ground truth parameters

% used the ground truth parameters from sub7 sess 7 SOA = -300 ms
TruePara                 = {[40, -20, 70, 50, 110, 80, 0.01, 0.06], ...%M1
    [40, -20, 70, 110, 80, 0.01, 0.06], ...%M2
    [40, -20, 70, 50, 110, 0.01, 0.06], ...%M3
    [40, -20, 70, 110, 0.01, 0.06], ...%M4
    [40, 70, 110, 0.01, 0.06]};%M5

%% functions

% define PMF
P_Afirst                 = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst                 = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous           = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

M1                       = {
    @(p) P_Vfirst(ExpInfo.SOA, p(1), p(3), p(5), p(7));...
    @(p) P_simultaneous(ExpInfo.SOA, p(1), p(3), p(5), p(7));...
    @(p) P_Afirst(ExpInfo.SOA, p(1), p(3), p(5), p(7));...
    @(p) P_Vfirst(ExpInfo.SOA, p(2), p(4), p(6), p(8));...
    @(p) P_simultaneous(ExpInfo.SOA, p(2), p(4), p(6), p(8));...
    @(p) P_Afirst(ExpInfo.SOA, p(2), p(4), p(6), p(8))};

M2                       = {
    @(p) P_Vfirst(ExpInfo.SOA, p(1), p(3), p(4), p(6));...
    @(p) P_simultaneous(ExpInfo.SOA, p(1), p(3), p(4), p(6));...
    @(p) P_Afirst(ExpInfo.SOA, p(1), p(3), p(4), p(6));...
    @(p) P_Vfirst(ExpInfo.SOA, p(2), p(3), p(5), p(7));...
    @(p) P_simultaneous(ExpInfo.SOA, p(2), p(3), p(5), p(7));...
    @(p) P_Afirst(ExpInfo.SOA, p(2), p(3), p(5), p(7))};

M3                       = {
    @(p) P_Vfirst(ExpInfo.SOA, p(1), p(3), p(5), p(6));...
    @(p) P_simultaneous(ExpInfo.SOA, p(1), p(3), p(5), p(6));...
    @(p) P_Afirst(ExpInfo.SOA, p(1), p(3), p(5), p(6));...
    @(p) P_Vfirst(ExpInfo.SOA, p(2), p(4), p(5), p(7));...
    @(p) P_simultaneous(ExpInfo.SOA, p(2), p(4), p(5), p(7));...
    @(p) P_Afirst(ExpInfo.SOA, p(2), p(4), p(5), p(7))};

M4                       = {
    @(p) P_Vfirst(ExpInfo.SOA, p(1), p(3), p(4), p(5));...
    @(p) P_simultaneous(ExpInfo.SOA, p(1), p(3), p(4), p(5));...
    @(p) P_Afirst(ExpInfo.SOA, p(1), p(3), p(4), p(5));...
    @(p) P_Vfirst(ExpInfo.SOA, p(2), p(3), p(4), p(6));...
    @(p) P_simultaneous(ExpInfo.SOA, p(2), p(3), p(4), p(6));...
    @(p) P_Afirst(ExpInfo.SOA, p(2), p(3), p(4), p(6))};

M5                       = {
    @(p) P_Vfirst(ExpInfo.SOA, p(1), p(2), p(3), p(4));...
    @(p) P_simultaneous(ExpInfo.SOA, p(1), p(2), p(3), p(4));...
    @(p) P_Afirst(ExpInfo.SOA, p(1), p(2), p(3), p(4));...
    @(p) P_Vfirst(ExpInfo.SOA, p(1), p(2), p(3), p(5));...
    @(p) P_simultaneous(ExpInfo.SOA, p(1), p(2), p(3), p(5));...
    @(p) P_Afirst(ExpInfo.SOA, p(1), p(2), p(3), p(5))};

models                   = {M1; M2; M3; M4; M5};

%% run model fitting for multiple times

% initialization
[deltaAIC, estP] = deal(cell(counts, nModel));

for count = 1:counts
    count

    %% simulate fake datasets for 5 models
    parfor iModel                   = 1:nModel

        % for the specific model and its corresponding ground-truth
        % parameter sets, generate the probability of reporting: Vfirst,
        % simultaneous, Afirst in pretest; Vfirst, simultaneous, Afirst in
        % posttest
        p_all_cond               = NaN(6, ExpInfo.lenS); % len of conditions x len of SOA levels
        Para                     = TruePara{iModel}; % true parameters for this model
        for iCondition           = 1:6
            p_all_cond(iCondition,:) = models{iModel}{iCondition}(Para);
        end

        % simulate data set
        sim_r_org                = cell(1,2);
        for s                    = 1:2 % loop through pre and post session
            p                        = p_all_cond(((s-1)*3+1):s*3, :); % pre and post probability of reporting V,simul, A
            sim_r_org{s}             = sampleMatrix(p, ExpInfo.nTrials);
        end

        %% fit fake data to 5 models 

        % obtain best-fitting parameters and AIC for each model
        [deltaAIC{count, iModel}, estP{count, iModel}]  = model_recovery_fit_fake_data_v2(sim_r_org, ExpInfo);

        % iBest is the index of model that best fits this specific fake data set
        [M iBEST] = min(deltaAIC{count, iModel});
        BEST = deltaAIC{count, iModel} == M;
        BEST = BEST / sum(BEST);
        CM(iModel,:) = CM(iModel,:) + BEST;

    end
end

%% summary plot of CM

figure(1); clf;
FM = round(100*CM/sum(CM(1,:)))/100;
t = imageTextMatrix(FM);
set(t(FM'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(CM);
set(t, 'fontsize', 22)
title(['count = ' num2str(count)]);
set(gca, 'xtick', [1:5], 'ytick', [1:5], 'fontsize', 28, ...
    'xaxislocation', 'top', 'tickdir', 'out')
xlabel('fit model')
ylabel('simulated model')