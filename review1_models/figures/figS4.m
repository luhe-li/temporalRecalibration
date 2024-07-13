% fig S4: outlier behavior

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
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% load

save_fig = 1;
n_model = numel(folders);
sub_slc = 5;

% TOJ prediction
for mm = 2

    result_folder = fullfile(projectDir, 'atheoretical_models_VBMC', folders{mm});
    files = dir(fullfile(result_folder, 'sub-*'));

    for ss = 1:numel(sub_slc)

        i_sub = sub_slc(ss);
        i_data = load(fullfile(result_folder, files(i_sub).name));
        pred= i_data.pred;

    end
    
end

% reorganize data
for ss = 1:numel(sub_slc)
    sub = sub_slc(ss);
    for  ses  = 1:9
        tempD(ses)  = organizeData(sub, ses);
    end
    [sorted_adaptor_soa, order] = sort([tempD.adaptor_soa]);
    D(ss).data = tempD(order);
end

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

    for  ss  = 1:numel(sub_slc)

        sub              = sub_slc(ss);

        figure; hold on
        set(gcf, 'Position', [0,0,1000,200]);

        tl = tiledlayout(2,9);
        sgtitle(sprintf('%s, outlier', specifications{mm}),'FontSize',titleSZ,'FontWeight','bold')

        for adapter = 1:9

            % pre
            nexttile(adapter); hold on 
            set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out','ColorOrder', cmp2)
            title(sprintf('Adapter SOA = %.1f s', adaptor_soa(adapter)./1e3),'FontSize',fontSZ,'FontWeight','bold')

            % data
            scatter(D(ss).data(adapter).pre_ms_unique, D(ss).data(adapter).pre_pResp, dotSZ, 'filled');

            % prediction
            plot(pred.test_soa, pred.pre_pmf{adapter}, 'LineWidth',lw)

            % look better
            xlim([-550 550])
            xticks([])
            yticks([])

            if adapter == 1
                ylabel('Pretest probability', 'FontSize',fontSZ,'FontWeight','bold')
                yticks(tick_y)
                yticklabels(strsplit(num2str(tick_y)))
            end

            % post
            nexttile(adapter+num_ses); hold on
            set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out','ColorOrder', cmp2)

            % data
            scatter(D(ss).data(adapter).post_ms_unique, D(ss).data(adapter).post_pResp, dotSZ, 'filled');

            % prediction
            plot(pred.test_soa, pred.post_pmf{adapter},'LineWidth',lw)
            %             xlabel('Adapter SOA','FontSize',fontSZ,'FontWeight','bold')
            xlabel('Test SOA','FontSize',fontSZ,'FontWeight','bold')

            % look better
            xlim([-550 550])
            xticks(tick_x)
            xticklabels(strsplit(num2str(tick_x./1e3)))
            yticks(tick_y)
            yticklabels(strsplit(num2str(tick_y)))

            if adapter == 1
                ylabel('Posttest probability','FontSize',fontSZ,'FontWeight','bold')
                yticks(tick_y)
                yticklabels(strsplit(num2str(tick_y)))
            end

        end

        tl.TileSpacing = 'compact';
        %         tl.XLabel.String = sprintf('S%i: %s',sub,concatenated_str);
        %         tl.XLabel.FontSize = titleSZ;
        %         tl.XLabel.FontWeight = 'bold';
        %         tl.XLabel.String = 'Adapter SOA';
        %         tl.XLabel.FontSize = titleSZ;
        %         tl.XLabel.FontWeight = 'bold';

        if save_fig
            flnm  = sprintf('M-%s_S%i_TOJ', folders{mm}, sub);
            saveas(gca, fullfile(out_dir, flnm),'pdf')
        end

    end

end
