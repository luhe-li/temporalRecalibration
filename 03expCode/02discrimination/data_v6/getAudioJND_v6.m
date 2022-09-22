function JND = getAudioJND_v6(subjID)

% This function obtains auditory JND at 75% accuracy.

% In data-cleaning, it extracts two key outputs:
% r_org :  this matrix has a size of length(s_unique) x numTrials
% respCount: this matrix has the size of 2(type of response, 1 =
% standard-more, 2 = comparison more) x length(s_unique)

% In fitting, we used log_weibull to fit volume in log scale, fit the PMF use
% fmincon for 1e3 times and chose the best-fitting parameters based on the
% smallest nLL.

% Data were then bootsrapped for 1e3 times to calculate the CI of
% parameters.

% In plotting, we plotted the raw data, the fitted PMF with the
% best-fitting parameters, and errorbars(?), and saved the figure

% It saves the fitting results as filename per subject.

%% clean raw data
flnm = [num2str(subjID) '_A'];
load(['discrimination_v6_sub' num2str(subjID) '_A.mat']);

% % fake response
% Response.more = Shuffle([ones(1,70), 2*ones(1,70)]);
% ExpInfo.VolumeRange = logspace(log10(0.5), log10(1), 7);

% for the auditory level, use measured db instead
% levels: physical intensity (volume) to reorganize data
% s_unique: same as levels
levels = ExpInfo.VolumeRange; 
s_unique = ExpInfo.VolumeRange; % unique intensity levels

% replace right click index from 3 to 2
Response.more(Response.more == 3) = 2;

% convert response from "which order is more" to "standard(1) or comparison(2) is
% more"
trialType = {[1,2];[2,1]}; % 1 = standard first; 2 = comparison first
cbTrials = reshape([trialType{ExpInfo.counterbalance}],2, [])';
for t = 1:ExpInfo.numTotalTrials
    RespComp(t) = cbTrials(t, Response.more(t)); % 1 = chose standard more; 2 = chose comparison more
end

% organize response based on intensity level
for i = 1:length(levels)
    iLevel = levels(i);
    iResp = RespComp(AudInfo.shuffledIntensity == iLevel);
    r_org(i,:) = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j = unique(RespComp) % 1 = chose standard more; 2 = chose comparison more
        respCount(j,i) = sum(iResp == j);
    end
end
pResp = respCount/ExpInfo.numTrials;

% plot raw data
cMAP = [215,48,39; 252,141,89; 254,224,144; 24,243,248; 145,191,219; 69,117,180]./255;
f1 = figure; hold on
l1 = scatter(s_unique, pResp(2,:),40,'MarkerFaceColor', cMAP(6,:),...
    'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.5);

%% define psychometric function

% define PMF
PF_cumGaussian = @(x, lambda, mu, sig) lambda/2 + (1 - lambda).* normcdf(x, mu, sig);
% PF_log_weibull = @(x, gamma, lambda, alpha, beta) gamma + (1 - lambda - gamma).* (1 - exp(-10.^(beta.*(x-alpha))));
% PF_weibull = @(x, gamma, lambda, alpha, beta) gamma + (1 - lambda - gamma).* (1 - exp( -(x./alpha).^beta));

% check PMF with arbiturary parameters
% int_finer = linspace(s_unique(1), s_unique(end), 1000);
% P =[0.01, 62, 0.2];
% P_comp_more = PF_cumGaussian(int_finer, P(1), P(2), P(3));
% figure; hold on
% plot(int_finer, P_comp_more)
% ylim([0 1.0])
% ylabel('proportion of perceiving comparison stimulus longer')
% xlabel('comparison volumn')
% xticks(round(levels,2))

%% define cost function

% set parameters
nT_compMore = respCount(2,:); % number of comparison stimulus more responses for each level
nT_standardMore = respCount(1,:); % number of standard stimulus more responses for each level
numTrials = ExpInfo.numTrials;

% nLL cost function
nLL = @(p) -nT_compMore * log(PF_cumGaussian(s_unique, p(1), p(2), p(3)))'...
    -nT_standardMore * log(1 - PF_cumGaussian(s_unique, p(1), p(2), p(3)))';

% set lower and upper bounds
% lambda, mu, sigma
lb      = [0, s_unique(1), 0];
ub      = [0.06, s_unique(end), 5];

% choose random initial values for 1e3 times
for i = 1:1e3
    init(i,:)    = rand(1,length(lb)).*(ub-lb) + lb;
    %You can also define how many times you want MATLAB to search
    options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');

    %fmincon returns best-fitting parameters that minimize the cost function as
    %well as the corresponding value for the cost function (in this case, the
    %negative log likelihood)
    [estP(i,:), min_NLL(i)] = fmincon(nLL, init(i,:),[],[],[],[],lb,ub,[],options);
end

% use the best-fitting parameters with the smallest NLL among 1e3 fittings
[value idx] = min(min_NLL);
bestP = estP(idx,:);
fprintf('lambda = %4.2f mu = %4.2f sig = %4.2f\n', bestP)

%% bootstrap to obtain error bars

% definee bounds for each parameter: gamma, lambda, alpha, beta
b_lb = [0, 0, 0];
b_ub = [0.06, 1, 2];
b_init = rand(1,length(b_lb)).*(b_ub-b_lb) + b_lb;
b_options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');

% bootstrap response at each level to obtain p
numBtst = 1e3;
for i = 1:numBtst
    disp(i) %display the counter so you can see the progress
	%initialize resampled responses
    r_slc = NaN(length(s_unique), numTrials);
    %for each stimulus location, we resample the responses
    for j = 1:length(s_unique)
        %randomly select trial indices (indices are allowed to occur more
        %than once, since we resample with replacement).
        %Hint: you'll find function randi.m useful
        idx        = randi([1 numTrials],[1 numTrials]);
        %store the resampled responses
        r_slc(j,:) = r_org(j,idx);
    end
    %compute the total number of comparison-more responses given each
    %intensity level
    nT_compMore_slc(i,:) = sum(r_slc == 2,2)';
end

btst_p =  nT_compMore_slc/numTrials;
for x                   = 1:length(s_unique)
    [lb(x), ub(x)]      = get68CI(btst_p(:,x));
end

%% plot with best-fitting parameter
figure(f1); hold on
set(gca,'FontSize',15,'linewidth',2)
int_finer = linspace(s_unique(1), s_unique(end), 1000);
P = bestP;
P_comp_more = PF_cumGaussian(int_finer, P(1), P(2), P(3));
l2 = plot(int_finer, P_comp_more,  'Color', cMAP(6,:), 'LineWidth',3);
errorbar(s_unique, pResp(2,:), pResp(2,:)-lb, ub-pResp(2,:), '.','Color', cMAP(6,:),'MarkerSize',4, 'LineWidth',2);
%  look better

ylim([0 1.0])
xlim([min(s_unique), max(s_unique)])
ylabel('p(comparison louder)')
xlabel('comparison volume')
title(['sub' flnm])
xticks(round(s_unique,2))
legend([l1, l2],{'data','fit'},'Location','northwest')
saveas(f1, flnm, 'epsc')

%% obtain JND as (x.90 - x.10)./2
[value idx]  = min(abs(P_comp_more - 0.95));
x95 = int_finer(idx);

[value idx]  = min(abs(P_comp_more - 0.05));
x05 = int_finer(idx);

JND = (x95 - x05)/2;
