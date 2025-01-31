% For each model, plot predictions of group/individual recalibration and TOJ responses

clear; close all; clc;

%% Model info

specifications = {'Exponential likelihood, shift criterion', 'Exponential likelihood, shift bias', 'Gaussian likelihood, shift criterion',  'Gaussian likelihood, shift bias'};
folders = {'exp_shiftC', 'exp_shiftMu', 'gauss_shiftC', 'gauss_shiftMu'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% Manage paths

restoredefaultpath;
current_dir = pwd;
[project_dir, ~] = fileparts(current_dir);
[git_dir, ~] = fileparts(project_dir);
addpath(genpath(fullfile(project_dir, 'data')));
addpath(genpath(fullfile(project_dir, 'utils')));
addpath(genpath(fullfile(git_dir, 'vbmc')));
addpath(genpath(fullfile(current_dir, curr_model_str)));
out_dir = fullfile(current_dir, curr_model_str);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% Load

save_fig = 1;
n_model = numel(folders);
sub_slc = [1:4, 6:10];

% TOJ prediction
for mm = 1:n_model

    result_folder = fullfile(project_dir, 'fit_results','atheoretical_models', folders{mm});
    files = dir(fullfile(result_folder, 'sub-*'));

    for ss = 1:numel(sub_slc)

        % Find the file that matches 'sub-XX'
        sub_str = sprintf('%02d', sub_slc(ss));
        file_name = '';
        for file = files'
            if contains(file.name, ['sub-', sub_str])
                file_name = file.name;
                break;
            end
        end

        % Load the data if the file was found
        if ~isempty(file_name)
            i_data = load(fullfile(result_folder, file_name));
            pred{mm, ss} = i_data.pred;
        else
            error('File for sub-%s not found.', sub_str);
        end

    end
    
end

% Reorganize data
for ss = 1:numel(sub_slc)
    sub = sub_slc(ss);
    for ses = 1:9
        temp_d(ses) = organize_data(sub, ses);
    end
    [sorted_adaptor_soa, order] = sort([temp_d.adaptor_soa]);
    D(ss).data = temp_d(order);
end

%% %%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure set up

cmp1 = [229, 158, 168; 203, 227, 172; 171, 223, 235] ./ 255;
cmp2 = [216, 49, 91; 175, 213, 128; 88, 193, 238] ./ 255;

lw = 0.5;
font_sz = 7;
title_sz = 9;
dot_sz = 10;

tick_y = 0:0.5:1;
tick_x = [-500, 0, 500];
test_soa = pred{1, 1}.test_soa; % ms
adaptor_soa = pred{1, 1}.adaptor_soa; % ms
num_ses = 9;

%% Plot TOJ

for mm = 1:n_model

    for ss = 1:numel(sub_slc)

        sub = sub_slc(ss);

        figure; hold on
        set(gcf, 'Position', [0, 0, 1000, 200]);

        tl = tiledlayout(2, 9);
        sgtitle(sprintf('%s, S%i', specifications{mm}, sub), 'FontSize', title_sz, 'FontWeight', 'bold')

        for adapter = 1:9

            % Pre
            nexttile(adapter); hold on 
            set(gca, 'FontSize', font_sz, 'LineWidth', lw, 'TickDir', 'out', 'ColorOrder', cmp2)
            title(sprintf('Adapter SOA = %.1f s', adaptor_soa(adapter) ./ 1e3), 'FontSize', font_sz, 'FontWeight', 'bold')

            % Data
            scatter(D(ss).data(adapter).pre_ms_unique, D(ss).data(adapter).pre_p_resp, dot_sz, 'filled');

            % Prediction
            plot(pred{mm, ss}.test_soa, pred{mm, ss}.pre_pmf{adapter}, 'LineWidth', lw)

            % Look better
            xlim([-550 550])
            xticks([])
            yticks([])

            if adapter == 1
                ylabel('Pretest probability', 'FontSize', font_sz, 'FontWeight', 'bold')
                yticks(tick_y)
                yticklabels(strsplit(num2str(tick_y)))
            end

            % Post
            nexttile(adapter + num_ses); hold on
            set(gca, 'FontSize', font_sz, 'LineWidth', lw, 'TickDir', 'out', 'ColorOrder', cmp2)

            % Data
            scatter(D(ss).data(adapter).post_ms_unique, D(ss).data(adapter).post_p_resp, dot_sz, 'filled');

            % Prediction
            plot(pred{mm, ss}.test_soa, pred{mm, ss}.post_pmf{adapter}, 'LineWidth', lw)
            xlabel('Test SOA', 'FontSize', font_sz, 'FontWeight', 'bold')

            % Look better
            xlim([-550 550])
            xticks(tick_x)
            xticklabels(strsplit(num2str(tick_x ./ 1e3)))
            yticks(tick_y)
            yticklabels(strsplit(num2str(tick_y)))

            if adapter == 1
                ylabel('Posttest probability', 'FontSize', font_sz, 'FontWeight', 'bold')
                yticks(tick_y)
                yticklabels(strsplit(num2str(tick_y)))
            end

        end

        tl.TileSpacing = 'compact';

        if save_fig
            flnm = sprintf('M-%s_S%i_TOJ', folders{mm}, sub);
            saveas(gca, fullfile(out_dir, flnm), 'png')
        end

    end

end
