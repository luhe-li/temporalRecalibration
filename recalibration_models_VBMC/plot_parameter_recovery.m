clear; close all;

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(cur_dir);
[git_dir, ~] = fileparts(project_dir);
dataDir = fullfile(git_dir,'temporalRecalibrationData');
addpath(genpath(fullfile(project_dir, 'utils')));
out_dir               = fullfile(cur_dir, mfilename);
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% load results

results_folder           = fullfile(dataDir,'recalibration_models_VBMC','param_recovery');
files = dir(fullfile(results_folder, 'sample-*'));

for jj = 1:size(files)

    r = load(fullfile(results_folder, files(jj).name));
    gt(jj,:) = r.summ.gt;
    est(jj,:) = r.summ.est;

end

% info
model_str = r.model.currModelStr;
paraID = r.model.initVal.paraID;
num_para = numel(paraID);
lb = r.model.initVal.lb;
ub = r.model.initVal.ub;

%% %%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%

nPlot = num_para;
nRow  = ceil(sqrt(nPlot));
nCol  = ceil(sqrt(nPlot));
if nPlot <= (nRow*nCol)-nCol, nRow = nRow-1; end

for jj = 1:num_para

    subplot(nRow,nCol,jj);
    axis square;
    axis equal
    hold on
    scatter(gt(:,jj), est(:,jj),50,'MarkerEdgeColor','k','MarkerFaceColor','none');
    xlim([lb(jj) ub(jj)])
    ylim([lb(jj) ub(jj)])

    % identity line
    ax = gca;
    x_limits = ax.XLim;
    y_limits = ax.YLim;
    line_min = min([x_limits y_limits]);
    line_max = max([x_limits y_limits]);
    plot([line_min line_max], [line_min line_max], 'k--', 'LineWidth', 1);

    % Calculate the Pearson correlation coefficient and p-value
    [R, P] = corrcoef(gt(:, jj), est(:, jj));

    % Extract the correlation coefficient and p-value
    r = R(1,2);
    p_value = P(1,2);

    % Label the r and p in the title
    title(sprintf('%s: r= %.2f, p=%.3f', paraID{jj}, r, p_value));

    if jj == 4
        ylabel('Prediction')
    elseif jj == 8
        xlabel('Ground-truth')
    end

end