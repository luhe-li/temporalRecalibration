% for each model, display estimates, and plot posterior

clear; close all; clc;

%% model info

specifications = {'Full model','No \sigma_{C1}', 'No \sigma_{C2}','No \sigma_{C1}, \sigma_{C2}'}; % Column 2: specifications
folders = {'cauInf_asym', 'cauInf_asym_xSigmaC1', 'cauInf_asym_xSigmaC2', 'cauInf_asym_xSigmaC1C2'}; % Column 3: folder names
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});


%% set environment

useCluster = true;

% set cores
if ~exist('useCluster', 'var') || isempty(useCluster)
    useCluster                  = false;
end

% job/sample = 100, core/sim_m = 6, fit once
switch useCluster
    case true
        if ~exist('numCores', 'var') || isempty(numCores)
            numCores  = maxNumCompThreads;
        end
        hpc_job_number = str2double(getenv('SLURM_ARRAY_TASK_ID'));
        if isnan(hpc_job_number), error('Problem with array assigment'); end
        fprintf('Job number: %i \n', hpc_job_number);
        i_sample = hpc_job_number;

        % make sure Matlab does not exceed this
        fprintf('Number of cores: %i  \n', numCores);
        maxNumCompThreads(numCores);
        if isempty(gcp('nocreate'))
            parpool(numCores);
        end

    case false
        numCores = 6;
        i_sample = 1;
end

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

%% load recal models

model_slc = [1,2,4];
n_model = numel(model_slc);
sub_slc = [1,3,4,6:10];
saveEst = 1;
dispLatex = 1;

for mm = 1:n_model
    result_folder = fullfile(dataDir, 'recalibration_cauInf', folders{model_slc(mm)});
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
    para.tableRowLabels              = [sprintfc('S%i', 1:numel(sub_slc)),'Mean','S.E.','Lower bound','Upper bound'];
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
% 
% for mm = 1:n_model
% 
%     f1 = figure;
%     set(gcf, 'Position', get(0, 'Screensize'));
% 
%     for ss = 1:numel(sub_slc)
% 
%         % sample posterior
%         Xs = post_samples{mm, ss};
%         Vals = R{mm, ss}.model.initVal;
%         n_para = numel(Vals.paraID);
% 
%         % plot posterior and check  tradeoff
%         cornerplot(Xs, Vals.paraID);
%         saveas(gca, fullfile(out_dir, sprintf('cornerplot_model-%i_sub-%i', mm, ss)),'png')
% 
%         %% plot individual parameters' posteriors and priors
%         plot_lb = min(Xs, [], 1);
%         plot_ub = max(Xs, [], 1);
% 
%         for pp = 1:n_para
%             figure(f1)
%             subplot(numel(sub_slc), n_para, (ss - 1) * n_para + pp)
%             hold on
% 
%             % plot posterior
%             % if posterior is too narrow, narrown the prior and fit again
%             h = histogram(Xs(:,pp),'normalization', 'probability','NumBins',100,'FaceColor',[0.5, 0.5, 0.5],'EdgeColor','none');
%             maxHistogramValue = max(h.Values);
% 
%             % plot prior
%             x = linspace(Vals.lb(pp),Vals.ub(pp),1e3);
%             y = msplinetrapezpdf(x', Vals.lb(pp), Vals.plb(pp),Vals.pub(pp),Vals.ub(pp));
%             y = y / max(y) * maxHistogramValue; % scale the y-values of the line plot by the max value of the histogram
%             plot(x,y,'r')
% 
%             if ss == 1
%                 title(Vals.paraID{pp});
%                 sgtitle(sprintf('Model: %s', R{mm, ss}.model.model_info.Specification{R{mm, ss}.model.i_model}))
%             end
% 
%             if pp == 1
%                 ylabel(sprintf('S%i', ss));
%             end
%         end
% 
%     end
% 
%     saveas(gca, fullfile(out_dir, sprintf('checkPrior_model-%i', mm)),'png')
% 
% end