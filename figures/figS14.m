% fig S14: outlier behavior

clear; clc; close all;

%% model info

specifications = {'Exponential likelihood, shift criterion', 'Exponential likelihood, shift bias', 'Gaussian likelihood, shift criterion',  'Gaussian likelihood, shift bias',};
folders = {'exp_shiftC', 'exp_shiftMu', 'gauss_shiftC', 'gauss_shiftMu'};
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
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load atheoretical results

save_fig = 1;
sub_slc = 5;

result_folder = fullfile(dataDir, 'atheoretical_models_VBMC', 'exp_shiftMu');
atheo = load_subject_data(result_folder, sub_slc, 'sub-*');
pred = atheo{1}.pred;

% reorganize data
for  ses  = 1:9
    tempD(ses)  = organizeData(sub_slc, ses);
end
[sorted_adaptor_soa, order] = sort([tempD.adaptor_soa]);
D.data = tempD(order);

%% %%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% figure set up

cmp1 = [229, 158, 168; 203, 227, 172; 171,223,235;]./255;
cmp2 = [216, 49, 91; 175, 213, 128; 88,193,238]./255;

lw = 0.5;
fontSZ = 7;
titleSZ = 9;
dotSZ = 10;

tick_y = 0:0.5:1;
tick_x = [-500, 0, 500];
test_soa = pred.test_soa; %ms
adaptor_soa = pred.adaptor_soa; %ms
num_ses = 9;

%% plot TOJ

for mm = 2%1:n_model

        figure; hold on
        set(gcf, 'Position', [0,0,800,150]);

        tl = tiledlayout(2,9);
        sgtitle('Temporal-order-judgmenst of excluded participant','FontSize',titleSZ,'FontWeight','bold')
       

        for adapter = 1:9

            % pre
            nexttile(adapter); hold on
            set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out','ColorOrder', cmp2)
            title(sprintf('Adapter SOA\n%.1f s', adaptor_soa(adapter)./1e3),'FontSize',fontSZ-1,'FontWeight','bold')

            % data
            scatter(D.data(adapter).pre_ms_unique, D.data(adapter).pre_pResp, dotSZ, 'filled');

            % prediction
            plot(pred.test_soa, pred.pre_pmf{adapter}, 'LineWidth',lw)

            % look better
            xlim([-550 550])
            xticks([])
            yticks([])

            if adapter == 1
                ylabel({'Pre-test', 'probability'}, 'FontSize',fontSZ,'FontWeight','bold')
                yticks(tick_y)
                yticklabels(strsplit(num2str(tick_y)))
            end

            % post
            nexttile(adapter+num_ses); hold on
            set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out','ColorOrder', cmp2)

            % data
            scatter(D.data(adapter).post_ms_unique, D.data(adapter).post_pResp, dotSZ, 'filled');

            % prediction
            plot(pred.test_soa, pred.post_pmf{adapter},'LineWidth',lw)
            if adapter == 5
                xlabel('Test SOA','FontSize',fontSZ,'FontWeight','bold')
            end

            % look better
            xlim([-550 550])
            xticks(tick_x)
            xticklabels(strsplit(num2str(tick_x./1e3)))
            yticks(tick_y)
            yticklabels(strsplit(num2str(tick_y)))

            if adapter == 1
                ylabel({'Post-test'; 'probability'},'FontSize',fontSZ,'FontWeight','bold')
                yticks(tick_y)
                yticklabels(strsplit(num2str(tick_y)))
            end

        end

        tl.TileSpacing = 'compact';

        if save_fig
            flnm  = sprintf('M-%s_S%s_TOJ', folders{mm}, sub_slc);
            saveas(gca, fullfile(out_dir, flnm),'pdf')
        end   

end
