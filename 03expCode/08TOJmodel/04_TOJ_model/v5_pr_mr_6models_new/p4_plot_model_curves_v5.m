% a04_p4_plot_model_curves_v5

% 1) loads best-fitting parameters from 6 models
% 2) load one session from one subject
% 3) compares psychometric function from different models on the same data
% to observe why the models differ, what effects changing one parameter
% impose on the fitting

clear all; close all; clc;

% set path
currentDir            = pwd;
exptDir               = currentDir(1:regexp(pwd,'03expCode')-1);
outDir                = [exptDir '02figures/fig3_model_comparison'];
addpath(genpath([exptDir '03ExpCode/05functions']));
addpath(genpath([exptDir '03ExpCode/06helperFunctions']));
addpath(genpath([exptDir '03ExpCode/01pretest/data']));
addpath(genpath([exptDir '03ExpCode/04posttest/data']));
addpath(genpath('model_comparison_function'))

% loads best-fitting parameters from 6 models
load('mc_results_v5_n10.mat')

% define subject and session to load, and models to compare
sub = 7;
ses = 7;
model2compare = 1:6;
num_model = model2compare(end);

for imodel = model2compare

    bestP = estP{sub}{ses, imodel};
    model_comp_plotting_v5(sub, ses, imodel, bestP);

    % save figure
    flnm = ['fig3_sub' num2str(sub) '_ses' num2str(ses) '_M' num2str(imodel)];
    saveas(gca, fullfile(outDir, flnm),'epsc')

end

%% scatter plot to compare fitting between models on one session, one subject

%% get data and fitting from models

% initiate
% data{imodel}(condition x test_soa)
% pmf_fit{imodel}{condition}{1 x test_soa}
[data, pmf_fit] = deal(cell(1,num_model));

for imodel = model2compare

    bestP = estP{sub}{ses, imodel};
    [data{imodel}, pmf_fit{imodel}] = model_comp_scatter_plot_fit_v5(sub, ses, imodel, bestP);

end

%% plotting

% color code = models
model2compare = 1:3;
n_color_level = numel(model2compare);

% plotting parameters
cmap = parula;
cidx = floor(linspace(1,230, n_color_level));
cond_title = {'p_{pre}("vision lead")', 'p_{pre}("simultaneous")', 'p_{pre}("audition lead")',...
    'p_{post}("vision lead")', 'p_{post}("simultaneous")','p_{post}("audition lead")'};

% initiate figure
figure; hold on
set(gcf,'position',[0 0 500 400]);
% title(cond_title{c});
xlabel('data'); ylabel('model fit');
line; axis square

for m = model2compare
    for c = 1:6
        i_data = data{m}(c,:);
        i_fit = pmf_fit{m}{c};
        scatter(i_data, i_fit, 50, 'o','MarkerFaceColor',cmap(cidx(m),:), 'MarkerEdgeColor','none','MarkerFaceAlpha', 0.5)
    end
    pause(1)
end

