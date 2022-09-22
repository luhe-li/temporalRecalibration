% This function plots the raw data with model fitting with a specified
% model its corresponding best-fitting parameters

% inputs:
% sub       : a scalar, subject id (e.g., 8)
% sess      : a scalar, session number (e.g., 9)
% imodel    : a scalar, index of one specific model (1-6)
% bestP     : a vector (size: 1 x num_parameters), must correspond with the
%           specified model (M1-M6 have different number of parameters)

function model_comp_plotting_v5(sub, ses, iModel, bestP)

%% organize data

%%%%% pre-test
% load data and define key parameters
load(['pretest_sub' num2str(sub) '_session' num2str(ses) '.mat'])
pre_s_unique                            = ExpInfo.SOA; % unique SOA levels, in s
pre_ms_unique                           = pre_s_unique * 1e3; % unique SOA levels, in ms
pre_numTrials                           = ExpInfo.nTrials; % num of trials per SOA
% inititate
pre_r_org                               = NaN(length(pre_s_unique), pre_numTrials);
pre_respCount                           = NaN(3, length(pre_s_unique));
for i                                   = 1:length(pre_s_unique)
    iSOA                                    = pre_s_unique(i);
    iResp                                   = Response.order(ExpInfo.trialSOA == iSOA);
    pre_r_org(i,:)                          = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j                                   = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        pre_respCount(j,i)                      = sum(iResp == j);
    end
end
pre_pResp                               = pre_respCount/pre_numTrials;

%%%%% post-test
% load data and define key parameters
load(['posttest_sub' num2str(sub) '_session' num2str(ses) '.mat'])
post_s_unique                           = ExpInfo.SOA; % unique SOA levels, in ms
post_ms_unique                          = post_s_unique * 1e3; % unique SOA levels, in s
post_numTrials                          = ExpInfo.nTrials; % num of trials per SOA
% inititate
post_r_org                              = NaN(length(post_s_unique), post_numTrials);
post_respCount                          = NaN(3, length(post_s_unique));
for i                                   = 1:length(post_s_unique)
    iSOA                                    = post_s_unique(i);
    iResp                                   = Response.order(ExpInfo.trialSOA == iSOA);
    post_r_org(i,:)                         = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j                                   = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        post_respCount(j,i)                     = sum(iResp == j);
    end
end
post_pResp                              = post_respCount/post_numTrials;

%% define and select best model

%make a finer grid for the timing difference between the auditory and the
%visual stimulus
SOA_finer                               = linspace(pre_ms_unique(1), pre_ms_unique(end), 1000);

% define PMF
P_Afirst                                = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst                                = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous                          = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

M1                                      = {@(p) P_Vfirst(SOA_finer, p(1), p(3), p(5), p(7));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(5), p(7));...
    @(p) P_Afirst(SOA_finer, p(1), p(3), p(5), p(7));...
    @(p) P_Vfirst(SOA_finer, p(2), p(4), p(6), p(8));...
    @(p) P_simultaneous(SOA_finer, p(2), p(4), p(6), p(8));...
    @(p) P_Afirst(SOA_finer, p(2), p(4), p(6), p(8))};

M2                                      = {@(p) P_Vfirst(SOA_finer, p(1), p(3), p(4), p(6));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(4), p(6));...
    @(p) P_Afirst(SOA_finer, p(1), p(3), p(4), p(6));...
    @(p) P_Vfirst(SOA_finer, p(2), p(3), p(5), p(7));...
    @(p) P_simultaneous(SOA_finer, p(2), p(3), p(5), p(7));...
    @(p) P_Afirst(SOA_finer, p(2), p(3), p(5), p(7))};

M3                                      = {@(p) P_Vfirst(SOA_finer, p(1), p(3), p(5), p(6));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(5), p(6));...
    @(p) P_Afirst(SOA_finer, p(1), p(3), p(5), p(6));...
    @(p) P_Vfirst(SOA_finer, p(2), p(4), p(5), p(7));...
    @(p) P_simultaneous(SOA_finer, p(2), p(4), p(5), p(7));...
    @(p) P_Afirst(SOA_finer, p(2), p(4), p(5), p(7))};

M4                                      = {@(p) P_Vfirst(SOA_finer, p(1), p(3), p(4), p(5));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(4), p(5));...
    @(p) P_Afirst(SOA_finer, p(1), p(3), p(4), p(5));...
    @(p) P_Vfirst(SOA_finer, p(2), p(3), p(4), p(6));...
    @(p) P_simultaneous(SOA_finer, p(2), p(3), p(4), p(6));...
    @(p) P_Afirst(SOA_finer, p(2), p(3), p(4), p(6))};

M5                                      = {@(p) P_Vfirst(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_Afirst(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_Vfirst(SOA_finer, p(1), p(2), p(4), p(6));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(4), p(6));...
    @(p) P_Afirst(SOA_finer, p(1), p(2), p(4), p(6))};

M6                                      = {@(p) P_Vfirst(SOA_finer, p(1), p(2), p(3), p(4));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(3), p(4));...
    @(p) P_Afirst(SOA_finer, p(1), p(2), p(3), p(4));...
    @(p) P_Vfirst(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_Afirst(SOA_finer, p(1), p(2), p(3), p(5))};

models                                  = {M1; M2; M3; M4; M5; M6};

%% fit PMF

pmf                                     = arrayfun(@(i) models{iModel}{i}(bestP), 1:6,'UniformOutput',false);

%% plotting

% set plotting parameters
cmp                                     = [229, 158, 168; 203, 227, 172; 171,223,235;...
    216, 49, 91; 175, 213, 128; 88,193,238]./255;
alpha_value                             = [repmat(0.6, 1, 3), repmat(0.6, 1, 3)];
line_pattern                            = {'--','--','--','-','-','-'};

% decide jitter direction by adaptor soa
if ExpInfo.adaptor <= 0
    jitter                                  = 10;
else
    jitter                                  = -10;
end

% initiate figure
figure; hold on
set(gca, 'LineWidth', 2, 'FontSize', 20)
set(gcf, 'Position',[10 10 700 400])

% plot raw data
for i                                   = 1:3
    dots(i) = plot(pre_ms_unique + jitter, pre_pResp(i,:), 'o',...
        'MarkerSize',7,'MarkerFaceColor', 'none','MarkerEdgeColor',cmp(i+3,:),'LineWidth',1.2);
    dots(i+3) = plot(post_ms_unique - jitter, post_pResp(i,:), 'o',...
        'MarkerSize',7,'MarkerFaceColor', cmp(i+3,:),'MarkerEdgeColor',cmp(i+3,:),'LineWidth',1.2);
end

% plot pmf (without error bars)
for i                                   = [4:6, 1:3] % so that lighter lines are at the front
    lines(i)  = plot(SOA_finer, pmf{i}, line_pattern{i},'Color', cmp(i,:), 'LineWidth', 2);
end

% plot mu
h1 = xline(bestP(1),':','LineWidth',2,'Color', repmat(0.3, 1, 3));
h4 = xline(bestP(2),'LineWidth',2,'Color', repmat(0.3, 1, 3));

% % display best fitting parameters
% annotation( 'textbox', 'String', round(bestP,3), ...
%     'FontSize', 14, 'Units', 'normalized', 'EdgeColor', 'none', ...
%     'Position', [0.9,0.9,0.03,0.03], 'Interpreter','Latex')

% add legend
legendlabels = {'p_{post}(vision lead)', 'p_{post}(simultaneous)','p_{post}(audition lead)',...
    'p_{pre}(vision lead)', 'p_{pre}(simultaneous)', 'p_{pre}(audition lead)',...
    'PSS_{pre}','PSS_{post}'};
legend([lines, h1, h4],legendlabels,'Location','bestoutside')

% look better
tick_soa                                = [-500, -300:100:300,500];
tick_p                                  = [0:0.25:1];
xticks(tick_soa)
xticklabels(strsplit(num2str(tick_soa)))
yticks(tick_p)
yticklabels(strsplit(num2str(tick_p)))
xlabel('SOA (ms)')
ylabel('probability')
title(['M' num2str(iModel)])
% title(['sub' num2str(sub) ' adaptor SOA = ' num2str(ExpInfo.adaptor*1000) ' ms'])

end