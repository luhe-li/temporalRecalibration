% Appendix 12-fig 1 
% Confusion matrix of model recovery of variants of
% causal-inference recalibration models with
% modality-specific precision

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

%% load fitting results

fit_model = [8,10]; % 8 = percept model, 10 = update model
sim_model = [8,10];
num_total_trials = 300*9;

for i = 1:numel(fit_model)

    for j = 1:numel(sim_model)

        fit_dir = fullfile(dataDir, 'recalibration_models','model_recovery_CIvariant',sprintf('M%i', fit_model(i)));
        flnm = sprintf('results_simM%i', sim_model(j));
        all_files = dir(fullfile(fit_dir, [flnm '*.mat']));

        nll = [];
        aic = [];
        bic = [];

        for k = 1:3 % num_batch

            load(fullfile(fit_dir,all_files(k).name))
            nll = [nll, modelRecov.minNLL];
            aic = [aic, 2*modelRecov.minNLL + 2*model.num_para];
            bic = [bic, 2*modelRecov.minNLL + model.num_para*log(num_total_trials)];

        end
       
        % sim model, fit model, repeat
        NLL(j, i, :) = nll;
        AIC(j, i, :) = aic;
        BIC(j, i, :) = bic;

    end

end

for i = 1:numel(sim_model)

    iAIC = squeeze(AIC(i,:,:));
    [~, bestM] = min(iAIC); % return index of model with smallest AIC
    CM(i,1) = sum(bestM == 1)./numel(bestM);
    CM(i,2) = sum(bestM == 2)./numel(bestM);

end

%% plot

figure; 
set(gcf, 'Position',[0,0,170,150]);

imagesc(CM); 
colormap('bone')
xticks(1:2)
yticks(1:2)
xticklabels({'Update','Percept'})
yticklabels({'Update','Percept'})
xlabel('Model used for fitting','FontWeight','bold');
ylabel('Data generating model','FontWeight','bold');

[num_rows, num_cols] = size(CM);
for row = 1:num_rows
    for col = 1:num_cols
        val = CM(row, col);
        
        textColor = 'w'; % default black
        if val > 0.5
            textColor = 'k'; % white for contrast
        end
        text(col, row, num2str(val, '%0.2f'), ...
            'HorizontalAlignment', 'center', ...
            'Color', textColor);
    end
end

saveas(gca, fullfile(out_dir, 'cm'),'pdf')