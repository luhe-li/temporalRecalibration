subjID = 9;
sess = 7;

pre_filenm = ['pretest_sub' num2str(subjID) '_session' num2str(sess) '.mat'];
post_filenm = ['posttest_sub' num2str(subjID) '_session' num2str(sess) '.mat'];

%% plot the raw figure

load(pre_filenm);

pre_s_unique = ExpInfo.SOA; % unique SOA levels, in ms
pre_numTrials = ExpInfo.nTrials;
for i = 1:length(pre_s_unique)
    iSOA = pre_s_unique(i);
    iResp = Response.order(ExpInfo.trialSOA == iSOA);
    pre_r_org(i,:) = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        pre_respCount(j,i) = sum(iResp == j);
    end
end
pre_pResp = pre_respCount/pre_numTrials;

% convert s_unique from s to ms for the following fitting and plotting
pre_s_unique = pre_s_unique*1000;

% plot raw data
cMAP1 = [158,202,225; 161,217,155; 252,174,145]./255;
cMAP2 = [33,113,181; 35,139,69; 203,24,29]./255;
alpha = 0.5; % transparency
f1 = figure;
set(f1,'Position',[10 10 600 450])
set(gca, 'FontSize', 15)
hold on
for i = 1:3
    scatter(pre_s_unique+5, pre_pResp(i,:),60,'MarkerFaceColor', cMAP1(i,:),...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',alpha)
end
xlabel('SOA (ms)')
ylabel('probability')

% posttest: raw data

load(post_filenm);

post_s_unique = ExpInfo.SOA; % unique SOA levels, in ms
post_numTrials = ExpInfo.nTrials;
for i = 1:length(post_s_unique)
    iSOA = post_s_unique(i);
    iResp = Response.order(ExpInfo.trialSOA == iSOA);
    post_r_org(i,:) = iResp; % this matrix has a size of length(s_unique) x numTrials
    for j = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
        post_respCount(j,i) = sum(iResp == j);
    end
end
post_pResp = post_respCount/post_numTrials;

% convert s_unique from s to ms for the following fitting and plotting
post_s_unique = post_s_unique*1000;

% plot raw data
figure(f1)
hold on
for i = 1:3
    scatter(post_s_unique-5, post_pResp(i,:),60,'MarkerFaceColor', cMAP2(i,:),...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',alpha)
end
xlabel('SOA (ms)')
ylabel('probability')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% fittings %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define the scaled psychometric function

% define PMF
P_Afirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous = @(SOA, mu, sig, lambda, c) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

% % check PMF with arbiturary parameters
% figure
% t_diff_finer = linspace(s_unique(1), s_unique(end), 1000);
% P = [-50, 50, 0.06, 100];
% P_Afirst_finer = P_Afirst(t_diff_finer, P(1), P(2), P(3), P(4));
% P_Vfirst_finer = P_Vfirst(t_diff_finer, P(1), P(2), P(3), P(4));
% P_simultaneous_finer = P_simultaneous(t_diff_finer, P(1), P(2), P(3), P(4));
%
% figure
% plot(t_diff_finer, P_Afirst_finer, 'Color', cMAP(3,:), ...
%     'lineWidth', 3); hold on;
% plot(t_diff_finer, P_Vfirst_finer, 'Color', cMAP(1,:), ...
%     'lineWidth', 3); hold on;
% plot(t_diff_finer, P_simultaneous_finer, 'Color', cMAP(2,:), ...
%     'lineWidth', 3); hold on;
% legend({'$P(A-precedes-V)$','$P(A-follows-V)$',...
%     '$P(A-coincides-V)$'},'Location','bestoutside', 'Interpreter','Latex');

%% define cost function
% set parameters
pre_nT_A1st = pre_respCount(3,:);% the number of A-first responses for each SOA
pre_nT_V1st = pre_respCount(1,:); % the number of V-first responses for each SOA

post_nT_A1st = post_respCount(3,:);% the number of A-first responses for each SOA
post_nT_V1st = post_respCount(1,:); % the number of V-first responses for each SOA

%nLL cost function
nLogL = @(p) -pre_nT_A1st*log(P_Afirst(pre_s_unique, p(1), p(3), p(4), p(5)))' ...
    -pre_nT_V1st*log(P_Vfirst(pre_s_unique, p(1), p(3), p(4), p(5)))'...
    -(repmat(pre_numTrials,size(pre_nT_A1st)) - pre_nT_A1st - pre_nT_V1st)...
    * log(P_simultaneous(pre_s_unique, p(1), p(3), p(4), p(5)))'...
    -post_nT_A1st*log(P_Afirst(post_s_unique, p(2), p(3), p(4), p(5)))' ...
    -post_nT_V1st*log(P_Vfirst(post_s_unique, p(2), p(3), p(4), p(5)))'...
    -(repmat(post_numTrials,size(post_nT_A1st)) - post_nT_A1st - post_nT_V1st)...
    * log(P_simultaneous(post_s_unique, p(2), p(3), p(4), p(5)))';

% set lower and upper bounds
lb      = [-150, -150, 10, 1e-2, 50];
ub      = [150, 150, 350, 0.06, 350];

% choose random initial values for 1e3 times
for i = 1:1e3
    init    = rand(1,length(lb)).*(ub-lb) + lb;
    %You can also define how many times you want MATLAB to search
    options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');

    %fmincon returns best-fitting parameters that minimize the cost function as
    %well as the corresponding value for the cost function (in this case, the
    %negative log likelihood)
    [estP(i,:), min_NLL(i)] = fmincon(nLogL, init,[],[],[],[],lb,ub,[],options);
end

% use the best-fitting parameters with the smallest NLL among 1e3 fittings
[value idx] = min(min_NLL);
bestP = estP(idx,:);
fprintf('pre mu = %4.2f post mu = %4.2f sig = %4.2f lambda = %4.2f criterion = %4.2f\n', bestP)

%% bootstrap
numBtst = 1e3;
[estP_btst, minNLL, lb_68CI, ub_68CI] = bootstrapPost(pre_s_unique, post_s_unique, ...
    pre_r_org, post_r_org, pre_numTrials, post_numTrials,...
    numBtst, P_Afirst, P_Vfirst, P_simultaneous, lb, ub, options) ;


%% plot with the best-fitting parametrs

% jitter raw data point
if ExpInfo.adaptor <= 0
    jitter = 5;
else
    jitter = -5;
end

%make a finer grid for the timing difference between the auditory and the
%visual stimulus
SOA_finer = linspace(pre_s_unique(1), pre_s_unique(end), 1000);

% plot raw data
f2                      = figure;
xlabel('SOA (ms)')
ylabel('probability')
set(f2,'Position',[10 10 600 450])
set(gca, 'FontSize', 15)
title(['sub' num2str(ExpInfo.subjID) ' adaptor SOA = ' num2str(ExpInfo.adaptor*1000) ' ms'])
hold on

for i                   = 1:3
    scatter(pre_s_unique + jitter, pre_pResp(i,:),'MarkerFaceColor', cMAP1(i,:),...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.7)
    scatter(post_s_unique - jitter, post_pResp(i,:),'MarkerFaceColor', cMAP2(i,:),...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.7)
end

saveas(gca,'sample_fig1','epsc')

%%%%%%%%%%%%%%%%%%%% pre %%%%%%%%%%%%%%%%%%%%

%plug in t_diff_finer and the best-fitting mu, sigma and lapse rate to
%function P
Pre_Vfirst_fit = P_Vfirst(SOA_finer, bestP(1), bestP(3), bestP(4), bestP(5));
Pre_Afirst_fit = P_Afirst(SOA_finer, bestP(1), bestP(3), bestP(4), bestP(5));
Pre_simul_fit = P_simultaneous(SOA_finer, bestP(1), bestP(3), bestP(4), bestP(5));

% errorbar of PMF
for i = 1:numBtst
    Pre_Vfirst_btst(i,:)            = P_Vfirst(SOA_finer, estP_btst(i,1), estP_btst(i,3), estP_btst(i,4), estP_btst(i,5));
    Pre_Afirst_btst(i,:)         = P_Afirst(SOA_finer, estP_btst(i,1),  estP_btst(i,3), estP_btst(i,4), estP_btst(i,5));
    Pre_simul_btst(i,:)      = P_simultaneous(SOA_finer, estP_btst(i,1),  estP_btst(i,3), estP_btst(i,4), estP_btst(i,5));
end

Pre_V_bounds = [max(Pre_Vfirst_btst, [], 1) - Pre_Vfirst_fit; Pre_Vfirst_fit - min(Pre_Vfirst_btst, [], 1)];
Pre_simul_bounds = [max(Pre_simul_btst, [], 1) - Pre_simul_fit; Pre_simul_fit - min(Pre_simul_btst, [], 1)];
Pre_A_bounds = [max(Pre_Afirst_btst, [], 1) - Pre_Afirst_fit; Pre_Afirst_fit - min(Pre_Afirst_btst, [], 1)];

l3 = shadedErrorBar(SOA_finer, Pre_Afirst_fit, Pre_A_bounds, ...
    'lineProps',{'--','Color', cMAP1(3,:), 'lineWidth', 3});
saveas(gca,'sample_fig2','epsc')

l1 = shadedErrorBar(SOA_finer, Pre_Vfirst_fit, Pre_V_bounds, ...
    'lineProps',{'--','Color', cMAP1(1,:), 'lineWidth', 3});
saveas(gca,'sample_fig3','epsc')

l2 = shadedErrorBar(SOA_finer, Pre_simul_fit, Pre_simul_bounds, ...
    'lineProps',{'--','Color', cMAP1(2,:), 'lineWidth', 3});
saveas(gca,'sample_fig4','epsc')

%%%%%%%%%%%%%%%%%%%% post %%%%%%%%%%%%%%%%%%%%
%plug in t_diff_finer and the best-fitting mu, sigma and lapse rate to
%function P
Post_Vfirst_fit = P_Vfirst(SOA_finer, bestP(2), bestP(3), bestP(4), bestP(5));
Post_Afirst_fit = P_Afirst(SOA_finer, bestP(2), bestP(3), bestP(4), bestP(5));
Post_simul_fit = P_simultaneous(SOA_finer, bestP(2), bestP(3), bestP(4), bestP(5));

% errorbar of PMF
for i = 1:numBtst
    Post_Vfirst_btst(i,:)            = P_Vfirst(SOA_finer, estP_btst(i,2), estP_btst(i,3), estP_btst(i,4), estP_btst(i,5));
    Post_Afirst_btst(i,:)         = P_Afirst(SOA_finer, estP_btst(i,2),  estP_btst(i,3), estP_btst(i,4), estP_btst(i,5));
    Post_simul_btst(i,:)      = P_simultaneous(SOA_finer, estP_btst(i,2),  estP_btst(i,3), estP_btst(i,4), estP_btst(i,5));
end

Post_V_bounds = [max(Post_Vfirst_btst, [], 1) - Post_Vfirst_fit; Post_Vfirst_fit - min(Post_Vfirst_btst, [], 1)];
Post_simul_bounds = [max(Post_simul_btst, [], 1) - Post_simul_fit; Post_simul_fit - min(Post_simul_btst, [], 1)];
Post_A_bounds = [max(Post_Afirst_btst, [], 1) - Post_Afirst_fit; Post_Afirst_fit - min(Post_Afirst_btst, [], 1)];

l6 = shadedErrorBar(SOA_finer, Post_Afirst_fit, Post_A_bounds, ...
    'lineProps',{'-','Color', cMAP2(3,:), 'lineWidth', 3});

l4 = shadedErrorBar(SOA_finer, Post_Vfirst_fit, Post_V_bounds, ...
    'lineProps',{'-','Color', cMAP2(1,:), 'lineWidth', 3});

l5 = shadedErrorBar(SOA_finer, Post_simul_fit, Post_simul_bounds, ...
    'lineProps',{'-','Color', cMAP2(2,:), 'lineWidth', 3});


% plot parameters
%     h4 = xline(bestP(2),'-','LineWidth',2);
%     h5 = xline(bestP(2) - bestP(5),'-','LineWidth',2);
%     h6 = xline(bestP(2) + bestP(5),'-','LineWidth',2);

saveas(gca,'sample_fig5','epsc')

