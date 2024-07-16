% for each model, display estimates, and plot posterior

clear; close all; clc;

%% model info

specifications = {'Heuristic, asymmetric', 'Heuristic, symmetric', 'Causal inference, asymmetric',  'Causal inference, symmetric','Atheoretical'}; % Column 2: specifications
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym','exp_shiftMu'}; % Column 3: folder names
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
[tempDir, ~] = fileparts(projectDir);
dataDir = fullfile(tempDir,'temporalRecalibrationData');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% atheoretical model for baseline
athe_path = fullfile(projectDir, 'atheoretical_models_VBMC','exp_shiftMu');

%% load recal models

model_slc = [1:4];
n_model = numel(model_slc);
sub_slc = [1:4,6:10];
saveEst = 1;
dispLatex = 1;

for mm = 1:n_model
    result_folder = fullfile(dataDir, 'recalibration_models_VBMC', folders{mm});
    R(mm, :) = load_subject_data(result_folder, sub_slc, 'sub-*');
    
    for ss = 1:numel(sub_slc)
        bestP{mm, ss} = R{mm, ss}.diag.post_mean;
        post_samples{mm, ss} = R{mm, ss}.diag.Xs;
    end
end

%% 1. display estimates

% summarize parameter estimates from one model
for mm = 1:n_model

    clearvars para est
    for ss                            = 1:numel(sub_slc)
        est(ss,:)                   = bestP{mm, ss};
    end

    est_mean                         = mean(est);
    est_std                          = std(est)./sqrt(numel(sub_slc));
    para.data                        = round([est; est_mean; est_std; R{mm, ss}.model.initVal.lb;R{mm, ss}.model.initVal.ub],3);
    para.tableRowLabels              = [sprintfc('S%i', 1:numel(sub_slc)),'Mean','S.D.','Lower bound','Upper bound'];
    para.tableColLabels              = R{mm, ss}.model.initVal.paraID;

    % save best parameter estimates
    if saveEst
        estTable = array2table(para.data, 'VariableNames', para.tableColLabels);
        estTable = addvars(estTable, para.tableRowLabels', 'Before', 1, 'NewVariableNames', 'RowLabels');
        writetable(estTable, fullfile(out_dir, sprintf('%s_bestParameterEstimates.csv', folders{mm})));
    end

    if dispLatex
%         para.dataFormat = {'%.2f', 4, '%.3f',3, '%.2f',2};
        latexTable(para);
    end
end

%% 2. check each model's posterior agaisnt prior to validate prior selection

for mm = 1:n_model

    f1 = figure;
    set(gcf, 'Position', get(0, 'Screensize'));

    for ss = 1:numel(sub_slc)

        % sample posterior
        Xs = post_samples{mm, ss};
        Vals = R{mm, ss}.model.initVal;
        n_para = numel(Vals.paraID);

        % plot posterior and check  tradeoff
        cornerplot(Xs, Vals.paraID);
        saveas(gca, fullfile(out_dir, sprintf('cornerplot_model-%i_sub-%i', mm, ss)),'png')

        %% plot individual parameters' posteriors and priors
        plot_lb = min(Xs, [], 1);
        plot_ub = max(Xs, [], 1);

        for pp = 1:n_para
            figure(f1)
            subplot(numel(sub_slc), n_para, (ss - 1) * n_para + pp)
            hold on

            % plot posterior
            % if posterior is too narrow, narrown the prior and fit again
            h = histogram(Xs(:,pp),'normalization', 'probability','NumBins',100,'FaceColor',[0.5, 0.5, 0.5],'EdgeColor','none');
            maxHistogramValue = max(h.Values);

            % plot prior
            x = linspace(Vals.lb(pp),Vals.ub(pp),1e3);
            y = msplinetrapezpdf(x', Vals.lb(pp), Vals.plb(pp),Vals.pub(pp),Vals.ub(pp));
            y = y / max(y) * maxHistogramValue; % scale the y-values of the line plot by the max value of the histogram
            plot(x,y,'r')

            if ss == 1
                title(Vals.paraID{pp});
                sgtitle(sprintf('Model: %s', R{mm, ss}.model.model_info.Specification{R{mm, ss}.model.i_model}))
            end

            if pp == 1
                ylabel(sprintf('S%i', ss));
            end
        end

    end

    saveas(gca, fullfile(out_dir, sprintf('checkPrior_model-%i', mm)),'png')

end