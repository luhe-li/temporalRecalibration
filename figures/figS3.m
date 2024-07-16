% fig S3: plot confidence interval of asymmetry index

clear; clc; close all;

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
[tempDir, ~] = fileparts(projectDir);
dataDir = fullfile(tempDir,'temporalRecalibrationData');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% set up

sub_slc = [1:4,6:10];

%% load bootstrapped results, extract pss_shift and calculate asymmetry index

result_folder = fullfile(projectDir, 'atheoretical_models_VBMC', 'exp_shiftMu');
atheo_btst = load_subject_data(result_folder, exp_sub, 'diag_btst_sub-*');

for ss = 1:numel(sub_slc)

        % Get asymmetry index of each bootstrap trials
        for jj = 1:numel(atheo_btst{ss}.pred)
            btst_ai(ss, jj)= sum(atheo_btst{ss}.pred{jj}.pss_shift);
        end
        [CI_lb(ss), CI_ub(ss)]  = get95CI(btst_ai(ss,:));

end

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 0.5;
fontsz = 7;
titleFontSz = 10;

figure;
set(gcf, 'Position', [0,0,210,180]);
set(gca, 'LineWidth', lw, 'FontSize', fontsz,'TickDir', 'out');
hold on 

for ss = 1:numel(sub_slc)

    plot([ss,ss],  [CI_lb(ss), CI_ub(ss)],'k','LineWidth',2)
   
end

yline(0,'--')
xlim([0,10])
ylim([-700, 700])
yticks([-700, 0, 700])
yticklabels([-0.7, 0, 0.7])
xticks(1:9)
xticklabels(1:9)
xlabel('Participant')
ylabel('Asymmetry index')

flnm = 'asymmetry_index_ci';
saveas(gca, fullfile(out_dir, flnm),'pdf')
