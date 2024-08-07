clear; close all;

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(cur_dir);
[git_dir, ~] = fileparts(project_dir);
dataDir = fullfile(git_dir,'temporalRecalibrationData');
addpath(genpath(fullfile(project_dir, 'utils')));
addpath(genpath(fullfile(project_dir, 'vbmc')));
out_dir               = fullfile(cur_dir, mfilename);
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% load results

results_folder           = fullfile(dataDir,'recalibration_models_VBMC','param_recovery');
files = dir(fullfile(results_folder, 'sample-*'));

for pp = 1:size(files)

    r = load(fullfile(results_folder, files(pp).name));

    try
        gt(pp,:) = r.summ.gt;
        est(pp,:) = r.summ.est;
    catch
        continue
    end

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

for pp = 1:num_para

    subplot(nRow,nCol,pp);
    axis square;
    axis equal
    hold on
    scatter(gt(:,pp), est(:,pp),30,'MarkerEdgeColor','k','MarkerFaceColor','none');
    xlim([lb(pp) ub(pp)])
    ylim([lb(pp) ub(pp)])

    % identity line
    ax = gca;
    x_limits = ax.XLim;
    y_limits = ax.YLim;
    line_min = min([x_limits y_limits]);
    line_max = max([x_limits y_limits]);
    plot([line_min line_max], [line_min line_max], 'k--', 'LineWidth', 1);

    % Calculate the Pearson correlation coefficient and p-value
    [R, P] = corrcoef(gt(:, pp), est(:, pp));

    % Extract the correlation coefficient and p-value
    r = R(1,2);
    p_value = P(1,2);

    % Label the r and p in the title
    title(sprintf('%s: r= %.2f, p=%.3f', paraID{pp}, r, p_value));

    if pp == 4
        ylabel('Prediction')
    elseif pp == 8
        xlabel('Ground-truth')
    end

end

flnm = 'param_recovery';
saveas(gca,fullfile(out_dir, flnm),'png')

%% check which paramaters trade off with p_common 

idx_pcc = 6;
rest_p = [1:5,7:9];
figure
for jj = 1:num_para-1

    pp = rest_p(jj);
    subplot(nRow,nCol,jj);
    axis square;
    hold on
    scatter(est(:,idx_pcc), est(:,pp),30,'MarkerEdgeColor','k','MarkerFaceColor','none');
    xlim([lb(idx_pcc) ub(idx_pcc)])
    ylim([lb(pp) ub(pp)])

    % Calculate the Pearson correlation coefficient and p-value
    [R, P] = corrcoef(est(:,idx_pcc), est(:, pp));

    % Extract the correlation coefficient and p-value
    r = R(1,2);
    p_value = P(1,2);

    % Label the r and p in the title
    title(sprintf('%s: r= %.2f, p=%.3f', paraID{pp}, r, p_value));

    if pp == 4
        ylabel('Prediction')
    elseif pp == 8
        xlabel('Ground-truth')
    end

end

flnm = 'pcc_tradeoff';
saveas(gca,fullfile(out_dir, flnm),'png')