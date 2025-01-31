% Appendix 13: 
% Parameter recovery of causal-inference model with modality-specific
% precision

clear; close all;

%% manage path

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
dataDir = fullfile(projectDir,'fit_results');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load results

results_folder           = fullfile(dataDir,'recalibration_models','param_recovery');
files = dir(fullfile(results_folder, 'sample-*'));

for jj = 1:size(files)

    r = load(fullfile(results_folder, files(jj).name));

    try
        gt(jj,:) = r.summ.gt;
        est(jj,:) = r.summ.est;
    catch
        continue
    end

end

% info
model_str = r.model.currModelStr;
paraID = {'\beta','\tau_A','\tau_V','Criterion','\lambda','p_{common}','\alpha','\sigma_{C=1}','\sigma_{C=2}'};
num_para = numel(paraID);
lb = r.model.initVal.lb;
ub = r.model.initVal.ub;

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 0.5;
fontSZ = 7;
titleSZ = 9;
dotSZ = 10;

figure;
set(gcf, 'Position', [0, 0, 420, 300]);

nPlot = num_para;
nRow  = ceil(sqrt(nPlot));
nCol  = ceil(sqrt(nPlot));
if nPlot <= (nRow*nCol)-nCol, nRow = nRow-1; end

for jj = 1:num_para

    subplot(nRow,nCol,jj);
    set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out')
    set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
    set(gca, 'FontName', 'Helvetica');
    axis square;
    axis equal
    hold on
    scatter(gt(:,jj), est(:,jj),10,'MarkerEdgeColor','k','MarkerFaceColor','none');
    xlim([lb(jj) ub(jj)])
    ylim([lb(jj) ub(jj)])

    % identity line
    ax = gca;
    x_limits = ax.XLim;
    y_limits = ax.YLim;
    line_min = min([x_limits y_limits]);
    line_max = max([x_limits y_limits]);
    plot([line_min line_max], [line_min line_max], 'k--', 'LineWidth', lw);

    % calculate the Pearson correlation coefficient and p-value
    [R, P] = corrcoef(gt(:, jj), est(:, jj));

    % extract the correlation coefficient and p-value
    r = R(1,2);
    p_value = P(1,2);

    if p_value < 0.01
        title(sprintf('%s \n r = %.2f, p < 0.01', paraID{jj}, r),'FontSize',fontSZ,'FontWeight','normal');
    else
        title(sprintf('%s \n r = %.2f, p = %.3f', paraID{jj}, r, round(p_value, 3)),'FontSize',fontSZ,'FontWeight','normal');
    end

    if jj == 4
        ylabel('Model prediction','FontSize',titleSZ)
    elseif jj == 8
        xlabel('Ground-truth','FontSize',titleSZ)
    end

end

flnm = 'param_recovery';
saveas(gcf, fullfile(out_dir, flnm), 'pdf');
