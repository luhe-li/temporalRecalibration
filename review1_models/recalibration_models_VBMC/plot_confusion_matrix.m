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

model_slc = [1,2];%[1,2,5];
n_model = numel(model_slc);
sub_slc = [1:4,6:10];

for mm = 1:n_model

    curr_folder = fullfile(pwd, folders{mm});
    if mm == 5;    curr_folder = athe_path;    end
    files = dir(fullfile(curr_folder, 'sub-*'));

    for ss = 1:numel(sub_slc)

        fprintf('Extract results of model-%i, sub-%i \n',mm, ss);

        i_sub = sub_slc(ss);
        i_data = load(fullfile(curr_folder, files(ss).name));
        DATA(mm, ss) = i_data;
        log_model_evi(mm, ss) = i_data.diag.bestELCBO;
        bestP{mm, ss} = i_data.diag.post_mean;
        pred{mm, ss} = i_data.pred;

    end
end

%% 1. plot model evidence

% max subtract other log model evidence
delta_LME = max(log_model_evi, [], 1) - log_model_evi; 

figure
h = heatmap(round(delta_LME, 1), 'XLabel','Participant', ...
    'Colormap', flipud(bone),...
    'ColorLimits', [0, 15], 'ColorbarVisible', 'on', 'GridVisible', 'off',...
    'FontSize', 8);
colorbar;
%     'ColorLimits', [0, 15], 'ColorbarVisible', 'on', 'GridVisible', 'off',...

h.YDisplayLabels = specifications;
h.XDisplayLabels = num2cell(1:numel(sub_slc));

% save figure
set(gca, 'FontSize', 8)
set(gcf, 'Position',[0 0 400 110])

flnm = 'ModelEvidence_recal_models';
% saveas(gca, fullfile(out_dir, flnm),'pdf')

%% 2. plot group average log model evidence

aa = 1;

%% 3. check each model's posterior agaisnt prior to validate prior selection

for mm = model_slc

    f1 = figure;
    set(gcf, 'Position', get(0, 'Screensize'));

    for ss = 1:numel(sub_slc)

        % sample posterior
        Xs = vbmc_rnd(bestELBO(mm, ss).vp,1e5);
        Vals = DATA(mm,ss).model.initVal;

        % plot posterior and check  tradeoff
        cornerplot(Xs, Vals.paraID);
        saveas(gca, fullfile(out_dir, sprintf('cornerplot_model-%i_sub-%i', mm, ss)),'pdf')

        %% plot individual parameters' posteriors and priors
        plot_lb = min(Xs, [], 1);
        plot_ub = max(Xs, [], 1);

        for pp = 1:n_para
            figure(f1)
            subplot(numel(sub_slc), n_para, (ss - 1) * n_para + pp)
            hold on

            % plot prior
            x = linspace(Vals.lb(pp),Vals.ub(pp),1e3);
            y = msplinetrapezpdf(x', Vals.lb(pp), Vals.plb(pp),Vals.pub(pp),Vals.ub(pp));
            plot(x,y,'r')

            % plot posterior
            % if posterior is too narrow, narrown the prior and fit again
            h = histogram(Xs(:,pp),'normalization', 'probability','NumBins',100,'FaceColor',[0.5, 0.5, 0.5],'EdgeColor','none');

            if ss == 1
                title(Vals.paraID{pp});
                sgtitle(sprintf('Model: %s', DATA(mm,ss).model.model_info.Specification{DATA(mm,ss).model.i_model}))
            end

            if pp == 1
                ylabel(sprintf('S%i', ss));
            end
        end
        saveas(gca, fullfile(out_dir, sprintf('checkPrior_model-%i_sub-%i', mm, ss)),'pdf')

    end

end
