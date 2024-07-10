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
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% atheoretical model for baseline
athe_path = fullfile(projectDir, 'atheoretical_models_VBMC','exp_shiftMu');

%% load recal models

model_slc = [1:4];%[1,2,5];
n_model = numel(model_slc);
sub_slc = [3,4,9];%[1:4,6:10];
saveEst = false;
dispLatex = false;

for mm = 1:n_model

    recal_folder = fullfile(projectDir, 'recalibration_models_VBMC', folders{mm});
    files = dir(fullfile(recal_folder, 'sub-*'));

    for ss = 1:numel(sub_slc)

        i_sub = sub_slc(ss);
        i_data = load(fullfile(recal_folder, files(ss).name));
        DATA(mm, ss) = i_data;
        bestP{mm, ss} = i_data.diag.post_mean;
        post_samples{mm, ss} = i_data.diag.Xs;

    end
end

%% 1. display estimates

% summarize parameter estimates from one model
for mm = 1:n_model

    clearvars para
    for ss                            = 1:numel(sub_slc)
        para.data(ss,:)                   = bestP{mm, ss};
    end

    est_mean                         = mean(para.data);
    est_std                          = std(para.data)./sqrt(numel(sub_slc));
    para.data                        = [para.data; est_mean; est_std; DATA(mm, ss).model.initVal.lb;DATA(mm, ss).model.initVal.ub];
    para.tableRowLabels              = [sprintfc('S%i', sub_slc),'Mean','S.D.','Lower bound','Upper bound'];
    para.tableColLabels              = DATA(mm, ss).model.initVal.paraID;

    % save best parameter estimates
    if saveEst
        estTable = array2table(para.data, 'VariableNames', para.tableColLabels);
        estTable = addvars(estTable, para.tableRowLabels', 'Before', 1, 'NewVariableNames', 'RowLabels');
        writematrix([folders{mm} '_bestParameterEstimates.csv']);
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
        Vals = DATA(mm,ss).model.initVal;
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
                sgtitle(sprintf('Model: %s', DATA(mm,ss).model.model_info.Specification{DATA(mm,ss).model.i_model}))
            end

            if pp == 1
                ylabel(sprintf('S%i', ss));
            end
        end

    end

    saveas(gca, fullfile(out_dir, sprintf('checkPrior_model-%i', mm)),'png')

end