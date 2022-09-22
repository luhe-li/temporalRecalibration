
%% plot SE on data

% calculate se of data
num_trials = pre_numTrials; % trial number per soa is the same for pre and post;
calSE = @(pv) sqrt(pv.*(1-pv)./num_trials);
pre_SE = arrayfun(@(i) calSE(pre_pResp(i,:)), 1:3, 'UniformOutput',false);
post_SE = arrayfun(@(i) calSE(post_pResp(i,:)), 1:3, 'UniformOutput',false);

% draw error bars on the data points
figure; hold on
for i = 1:3
    e = errorbar(pre_ms_unique + jitter + i*10, pre_pResp(i,:), pre_SE{i}, ...
        'o','MarkerSize',8,'MarkerFaceColor', cmp(i,:), 'MarkerEdgeColor','none');
    e.Color = cmp(i,:); e.LineWidth = 1; e.CapSize = 0;
    e = errorbar(post_ms_unique - jitter - i*10, post_pResp(i,:), post_SE{i}, ...
        'o','MarkerSize',8,'MarkerFaceColor', cmp(i+3,:), 'MarkerEdgeColor','none');
    e.Color = cmp(i+3,:); e.LineWidth = 1; e.CapSize = 0;
end

figure; hold on
% plot shaded regions around data points
for i = 1:3
    [l p] = boundedline(post_ms_unique - jitter, post_pResp(i,:), post_SE{i}, 'o','cmap', cmp(i,:),'lineWidth', 0.1,'alpha');
    l.MarkerSize = 6; 
    l.MarkerFaceColor = cmp(i,:);
    l.MarkerEdgeColor = 'none';
end

plot(SOA_finer, Pre_Vfirst_fit, 'Color', cmp(1,:),'LineWidth',1.5)









%% plot SE on psychometric curves
%plug in t_diff_finer and the best-fitting mu, sigma and lapse rate to
%function P
Pre_Vfirst_fit = P_Vfirst(SOA_finer, bestP(1), bestP(3), bestP(5), bestP(7));
Pre_Afirst_fit = P_Afirst(SOA_finer, bestP(1), bestP(3), bestP(5), bestP(7));
Pre_simul_fit = P_simultaneous(SOA_finer, bestP(1), bestP(3), bestP(5), bestP(7));

% errorbar of PMF
numBtst = 1000;
for i = 1:numBtst
    Pre_Vfirst_btst(i,:)            = P_Vfirst(SOA_finer, i_estP_btst(i,1), i_estP_btst(i,3), i_estP_btst(i,5), i_estP_btst(i,7));
    Pre_Afirst_btst(i,:)         = P_Afirst(SOA_finer, i_estP_btst(i,1), i_estP_btst(i,3), i_estP_btst(i,5), i_estP_btst(i,7));
    Pre_simul_btst(i,:)      = P_simultaneous(SOA_finer, i_estP_btst(i,1), i_estP_btst(i,3), i_estP_btst(i,5), i_estP_btst(i,7));
end

Pre_V_bounds = [max(Pre_Vfirst_btst, [], 1) - Pre_Vfirst_fit; Pre_Vfirst_fit - min(Pre_Vfirst_btst, [], 1)];
Pre_simul_bounds = [max(Pre_simul_btst, [], 1) - Pre_simul_fit; Pre_simul_fit - min(Pre_simul_btst, [], 1)];
Pre_A_bounds = [max(Pre_Afirst_btst, [], 1) - Pre_Afirst_fit; Pre_Afirst_fit - min(Pre_Afirst_btst, [], 1)];

l1 = shadedErrorBar(SOA_finer, Pre_Vfirst_fit, Pre_V_bounds, ...
    'lineProps',{'--','Color', cmp(4,:), 'lineWidth', 3});
l2 = shadedErrorBar(SOA_finer, Pre_simul_fit, Pre_simul_bounds, ...
    'lineProps',{'--','Color', cmp(5,:), 'lineWidth', 3});
l3 = shadedErrorBar(SOA_finer, Pre_Afirst_fit, Pre_A_bounds, ...
    'lineProps',{'--','Color', cmp(6,:), 'lineWidth', 3});

% plot parameters
h1                      = xline(bestP(1),'--','LineWidth',2);
h2                      = xline(bestP(1) - bestP(5),'--','LineWidth',2);
h3                      = xline(bestP(1) + bestP(5),'--','LineWidth',2);

%%%%%%%%%%%%%%%%%%%% post %%%%%%%%%%%%%%%%%%%%
%plug in t_diff_finer and the best-fitting mu, sigma and lapse rate to
%function P
Post_Vfirst_fit = P_Vfirst(SOA_finer, bestP(2), bestP(4), bestP(6), bestP(8));
Post_Afirst_fit = P_Afirst(SOA_finer, bestP(2), bestP(4), bestP(6), bestP(8));
Post_simul_fit = P_simultaneous(SOA_finer, bestP(2), bestP(4), bestP(6), bestP(8));

% errorbar of PMF
for i = 1:numBtst
    Post_Vfirst_btst(i,:)     = P_Vfirst(SOA_finer, i_estP_btst(i,2), i_estP_btst(i,4), i_estP_btst(i,6), i_estP_btst(i,8));
    Post_Afirst_btst(i,:)     = P_Afirst(SOA_finer, i_estP_btst(i,2), i_estP_btst(i,4), i_estP_btst(i,6), i_estP_btst(i,8));
    Post_simul_btst(i,:)      = P_simultaneous(SOA_finer, i_estP_btst(i,2), i_estP_btst(i,4), i_estP_btst(i,6), i_estP_btst(i,8));
end

Post_V_bounds = [max(Post_Vfirst_btst, [], 1) - Post_Vfirst_fit; Post_Vfirst_fit - min(Post_Vfirst_btst, [], 1)];
Post_simul_bounds = [max(Post_simul_btst, [], 1) - Post_simul_fit; Post_simul_fit - min(Post_simul_btst, [], 1)];
Post_A_bounds = [max(Post_Afirst_btst, [], 1) - Post_Afirst_fit; Post_Afirst_fit - min(Post_Afirst_btst, [], 1)];

l4 = shadedErrorBar(SOA_finer, Post_Vfirst_fit, Post_V_bounds, ...
    'lineProps',{'-','Color', cmp(1,:), 'lineWidth', 3});
l5 = shadedErrorBar(SOA_finer, Post_simul_fit, Post_simul_bounds, ...
    'lineProps',{'-','Color', cmp(2,:), 'lineWidth', 3});
l6 = shadedErrorBar(SOA_finer, Post_Afirst_fit, Post_A_bounds, ...
    'lineProps',{'-','Color', cmp(3,:), 'lineWidth', 3});

% plot parameters
h4 = xline(bestP(2),'-','LineWidth',2);
h5 = xline(bestP(2) - bestP(5),'-','LineWidth',2);
h6 = xline(bestP(2) + bestP(5),'-','LineWidth',2);

%% darker shade, lighter line
% set different colors for lines and shades (revers the two)
figure; hold on
for i = 1:6
ls(i) = shadedErrorBar(SOA_finer, pmf(i,:), bounds{i}, ...
    'lineProps',{line_pattern{i}, 'Color', cmp(i,:),'lineWidth', 2},'patchSaturation', alpha_value);
% change shade color and transparency here
    ls(i).patch.FaceColor = patch_cmp(i,:); ls(i).patch.FaceAlpha = 0.6;
end

