function posttest_fitting_ind(pre_filenm, post_filenm, saveFigure, saveData)

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
hold on
for i = 1:3
    scatter(pre_s_unique+5, pre_pResp(i,:),60,'MarkerFaceColor', cMAP1(i,:),...
        'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',alpha)
end
xlabel('SOA (ms)')
ylabel('probability')

%% posttest: raw data

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

if saveFigure
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% plot btst mu difference
    figure; hold on
    mu_diff = estP_btst(:,2) - estP_btst(:,1);
    h = histogram(mu_diff, 100);
    h.FaceColor = cMAP2(1,:);
    % add 95% CI
    [lb_95, ub_95] = get95CI (mu_diff);
    yl = ylim;
    CIline = plot([lb_95, ub_95],[yl(2)-2, yl(2)-2],'-r','LineWidth', 2);
    % add x=0
    xline(0,'LineWidth', 2);
    % look better
    xlim([-150 150])
    xlabel('\mu_{post} - \mu_{pre}')
    ylabel('count')
    set(gca, 'FontSize', 15)
    legend(CIline, {'95% CI'})
    title(['sub' num2str(ExpInfo.subjID) ' adaptor SOA = ' num2str(ExpInfo.adaptor*1000) ' ms'])
    fignm = ['mu_diff_sub', num2str(ExpInfo.subjID) '_session' num2str(ExpInfo.session)];
    saveas(gca,fignm,'epsc')

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
    f2                      = figure; hold on
    for i                   = 1:3
        scatter(pre_s_unique + jitter, pre_pResp(i,:),'MarkerFaceColor', cMAP1(i,:),...
            'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.7)
        scatter(post_s_unique - jitter, post_pResp(i,:),'MarkerFaceColor', cMAP2(i,:),...
            'MarkerEdgeAlpha', 0, 'MarkerFaceAlpha',0.7)
    end

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

    l1 = shadedErrorBar(SOA_finer, Pre_Vfirst_fit, Pre_V_bounds, ...
        'lineProps',{'--','Color', cMAP1(1,:), 'lineWidth', 3});
    l2 = shadedErrorBar(SOA_finer, Pre_simul_fit, Pre_simul_bounds, ...
        'lineProps',{'--','Color', cMAP1(2,:), 'lineWidth', 3});
    l3 = shadedErrorBar(SOA_finer, Pre_Afirst_fit, Pre_A_bounds, ...
        'lineProps',{'--','Color', cMAP1(3,:), 'lineWidth', 3});

    % plot parameters
    h1                      = xline(bestP(1),'--','LineWidth',2);
    h2                      = xline(bestP(1) - bestP(5),'--','LineWidth',2);
    h3                      = xline(bestP(1) + bestP(5),'--','LineWidth',2);

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

    l4 = shadedErrorBar(SOA_finer, Post_Vfirst_fit, Post_V_bounds, ...
        'lineProps',{'-','Color', cMAP2(1,:), 'lineWidth', 3});
    l5 = shadedErrorBar(SOA_finer, Post_simul_fit, Post_simul_bounds, ...
        'lineProps',{'-','Color', cMAP2(2,:), 'lineWidth', 3});
    l6 = shadedErrorBar(SOA_finer, Post_Afirst_fit, Post_A_bounds, ...
        'lineProps',{'-','Color', cMAP2(3,:), 'lineWidth', 3});

    % plot parameters
    h4 = xline(bestP(2),'-','LineWidth',2);
    h5 = xline(bestP(2) - bestP(5),'-','LineWidth',2);
    h6 = xline(bestP(2) + bestP(5),'-','LineWidth',2);

    %%%%%%%%%%%%%%%%%%%% look better %%%%%%%%%%%%%%%%%%%%

    % display best-fitting parameters
    bestP = round(bestP,2);
    text = {['$\mu_{pre} = $' num2str(bestP(1))], ...
        ['$\mu_{post} = $' num2str(bestP(2))],...
        ['$\sigma = $' num2str(bestP(3))],...
        ['$\lambda = $' num2str(bestP(4))],...
        ['$criterion = $' num2str(bestP(5))]};
    annotation( 'textbox', 'String', text, ...
        'FontSize', 14, 'Units', 'normalized', 'EdgeColor', 'none', ...
        'Position', [0.72,0.5,0.4,0.2], 'Interpreter','Latex')

    % look better
    legend([l1.mainLine, l2.mainLine, l3.mainLine, l4.mainLine, l5.mainLine, l6.mainLine],...
        {'$P_{pre}(A-precedes-V)$','$P_{pre}(A-follows-V)$',...
        '$P_{pre}(A-coincides-V)$', '$P_{post}(A-precedes-V)$','$P_{post}(A-follows-V)$',...
        '$P_{post}(A-coincides-V)$'},'Location','bestoutside', 'Interpreter','Latex');

    xlabel('SOA (ms)')
    ylabel('probability')
    set(f2,'Position',[10 10 900 600])
    set(gca, 'FontSize', 15)
    title(['sub' num2str(ExpInfo.subjID) ' adaptor SOA = ' num2str(ExpInfo.adaptor*1000) ' ms'])
    fignm = ['posttest_sub', num2str(ExpInfo.subjID) '_session' num2str(ExpInfo.session)];

    saveas(gca,fignm,'epsc')
end

if saveData
    %% summarize and save data
    % 68% error bar of delta_mu
    mu_diff = estP_btst(:,2) - estP_btst(:,1);
    [lb_68, ub_68] = get68CI(mu_diff);

    % save data
    if exist('data_summary.mat') == 2
        load('data_summary.mat')
    end
    adaptor{ExpInfo.subjID, ExpInfo.session} = ExpInfo.adaptor*1000;
    delta_mu{ExpInfo.subjID, ExpInfo.session} = bestP(2) - bestP(1);
    btst_delta_mu{ExpInfo.subjID, ExpInfo.session} = mu_diff;
    delta_mu_lb68{ExpInfo.subjID, ExpInfo.session} = lb_68;
    delta_mu_ub68{ExpInfo.subjID, ExpInfo.session} = ub_68;
    save('data_summary','adaptor','delta_mu','btst_delta_mu','delta_mu_lb68','delta_mu_ub68')
end
end