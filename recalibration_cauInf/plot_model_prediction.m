% for each model, plot predictions of group/individual recalibration and TOJ
% responses

clear; clc; close all;

%% model info

specifications = {'Full model','No \sigma_{C1}', 'No \sigma_{C2}','No \sigma_{C1}, \sigma_{C2}'}; % Column 2: specifications
folders = {'cauInf_asym', 'cauInf_asym_xSigmaC1', 'cauInf_asym_xSigmaC2', 'cauInf_asym_xSigmaC1C2'}; % Column 3: folder names
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

%% load recalibration model results

save_fig = 1;
plot_toj = 0;

model_slc = [1,2,4];
n_model = numel(model_slc);
sub_slc = [1,3,4,6:10];

for mm = 1:n_model
    result_folder = fullfile(dataDir, 'recalibration_cauInf', folders{model_slc(mm)});
  
    R(mm, :) = load_subject_data(result_folder, sub_slc, 'sub-*');
    
    for ss = 1:numel(sub_slc)
        pred{mm, ss} = R{mm, ss}.pred;
        % Extract summary data for plot
        pred_recal(mm, ss, :) = mean(pred{mm, ss}.pss_shift, 2);
    end
end

%% load atheoretical model results

result_folder = fullfile(dataDir, 'atheoretical_models_VBMC', 'exp_shiftMu');
atheo = load_subject_data(result_folder, sub_slc, 'sub-*');

% Initialize toj_pss based on the first loaded file
toj_pss = zeros(numel(sub_slc), size(atheo{1}.pred.pss_shift, 2));
for ss = 1:numel(sub_slc)
    toj_pss(ss, :) = atheo{ss}.pred.pss_shift;
end

% Calculate group mean and standard error
mean_toj_pss = mean(toj_pss, 1, 'omitnan');
se_toj_pss = std(toj_pss, [], 1, 'omitnan') ./ sqrt(numel(sub_slc));

% Reorganize data
D = struct([]);
for ss = 1:numel(sub_slc)
    for ses = 1:9
        tempD(ses) = organizeData(sub_slc(ss), ses);
    end
    [sorted_adaptor_soa, order] = sort([tempD.adaptor_soa]);
    D(ss).data = tempD(order);
end

%% load bootstrap results of atheoretical models

btst_data = load_subject_data(result_folder, sub_slc, 'diag_btst_sub-*');

btst_pss = [];
for ss = 1:numel(sub_slc)
    btst = btst_data{ss};
    for jj = 1:numel(btst.pred)
        btst_pss(jj, :) = btst.pred{jj}.pss_shift;
    end
    
    for tt = 1:size(btst_pss, 2)
        [lb(ss, tt), ub(ss, tt)] = get68CI(btst_pss(:, tt));
    end
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
test_soa = pred{1,1}.test_soa; %ms
adaptor_soa = pred{1,1}.adaptor_soa; %ms
num_ses = 9;

%% 1. plot group recalibration

figure;
set(gcf, 'Position',[0,0,420,250]);
set(gcf, 'DefaultAxesFontName', 'Helvetica');
set(gcf, 'DefaultTextFontName', 'Helvetica');
t = tiledlayout(2,3,'Padding', 'compact', 'TileSpacing', 'compact');

adaptor_soa = pred{1,1}.adaptor_soa; %ms

order = [6,2,4, 5,1,3];
yl = 100;
ytks = {[], [],[],[],[-yl, 0, yl], [-yl, 0, yl]};
ytklabels = {[], [], [], [], [-yl, 0, yl]./1e3, [-yl, 0, yl]./1e3};

for mm = 1:n_model

    nexttile; hold on

    set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out')
    set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
    set(gca, 'FontName', 'Helvetica');

    e = errorbar(adaptor_soa, mean_toj_pss, se_toj_pss, ...
        'o','LineWidth', 1.5);
    e.CapSize = 0; e.Color = 'k'; e.MarkerFaceColor = 'k';
    e.LineWidth = lw;
    e.MarkerSize   = 2;

    % color parameter for model fitting
    clt = repmat(0.7, 1, 3);

    mean_sim_recal = mean(squeeze(pred_recal(mm,:,:)), 1, 'omitnan');
    se_sim_recal = std(squeeze(pred_recal(mm,:,:)), [], 1, 'omitnan')./sqrt(numel(sub_slc));

    % calculate lower and upper bound using se
    lb_sim_recal  = mean_sim_recal - se_sim_recal;
    ub_sim_recal  = mean_sim_recal + se_sim_recal;

    % plot fitting range
    p = patch([adaptor_soa, fliplr(adaptor_soa)], [lb_sim_recal, fliplr(ub_sim_recal)], ...
        clt, 'FaceAlpha', 0.2, 'EdgeColor','none');

    ylim([-yl, yl])
    yticks(ytks{mm})
    yticklabels(ytklabels{mm})
    yline(0,'--','LineWidth',lw)
    xticks(adaptor_soa)
    xticklabels({'-0.7','-0.3','','','0','','','0.3','0.7'})
    xtickangle(0)
    xlim([min(adaptor_soa)-50, max(adaptor_soa)+50])

    xlabel(t, 'Adaptor SOA (s)','FontSize',titleSZ);
    ylabel(t,'Recalibration effect (s)','FontSize',titleSZ);

    title(specifications{model_slc(mm)},'FontSize',fontSZ,'FontWeight', 'normal');

end

saveas(gca, fullfile(out_dir,'Group_recal'),'png')

%% 2. plot individual recalibration

for mm = 1:n_model

    figure; hold on
    set(gcf, 'Position',[0, 0, 600, 400]);

    sgtitle(specifications{model_slc(mm)},'FontSize',titleSZ)

    for ss = 1:numel(sub_slc)

        sub = sub_slc(ss);

        %% plot
        subplot(3,3,ss); hold on
        set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out')
        set(gca, 'LineWidth', lw, 'FontSize', fontSZ,'TickDir', 'out')
        set(gca, 'FontName', 'Helvetica');
        title(['S' num2str(ss)],'FontSize',titleSZ)

        %% plot atheoretical prediction

        blackL  = plot(adaptor_soa, toj_pss(ss,:),'ko', 'MarkerFaceColor','k','MarkerSize',3);

        for jj = 1:9
            plot([adaptor_soa(jj),adaptor_soa(jj)],[lb(ss, jj), ub(ss, jj)],'k-','LineWidth',1)
        end

        %% plot model prediction

        redL = plot(adaptor_soa, squeeze(pred_recal(mm, ss,:)), '-o','LineWidth',lw, 'Color','r','MarkerSize',3);

        % look better
        yl = 250;
        ylim([-yl, yl])
%         yticks([-yl, 0, yl])
%         yticklabels([-yl, 0, yl]./1e3)
        yline(0,'--','LineWidth',1)
        xticks(adaptor_soa)
        xticklabels(adaptor_soa/1e3)
        xtickangle(60)
        xlim([min(adaptor_soa)-100, max(adaptor_soa)+100])

%         % Enable the grid
%         grid on;
%         ax = gca;
%         ax.Color = 0.95*ones(1,3); % Light gray background color
%         ax.GridColor = [1, 1, 1]; % White grid lines
%         ax.GridAlpha = 1.0; % Opaque grid lines
%         ax.XColor = 'k'; % X-axis color
%         ax.YColor = 'k'; % Y-axis color
%         ax.GridLineStyle = '-'; % Solid lines for the grid
        ax.YTick = [-0.2, -0.1, 0, 0.1, 0.2]*1e3;
        ax.YTickLabel = [-0.2, -0.1, 0, 0.1, 0.2];
       
        if ss == 4
            ylabel('Recalibration effect (s)','FontWeight','bold','FontSize',titleSZ)
        elseif ss == 8
            xlabel('Adapter SOA (s)','FontWeight','bold','FontSize',titleSZ)
        end


        if ss == numel(sub_slc)
           legend([blackL, redL],{'Data','Model prediction'})
           saveas(gca, fullfile(out_dir, sprintf('M-%s_indiv_recal', folders{mm})),'png')
        end
        
    end

end

%% 3. plot TOJ 

if plot_toj

    for mm = [1,3]

        for  ss  = 1:numel(sub_slc)

            sub              = sub_slc(ss);

            figure; hold on
            set(gcf, 'Position', [0,0,1000,200]);

            tl = tiledlayout(2,9);
            sgtitle(sprintf('%s, S%i', specifications{model_slc(mm)}, ss),'FontSize',titleSZ,'FontWeight','bold')

            for adapter = 1:9

                % pre
                nexttile(adapter); hold on
                set(gca, 'FontSize', fontSZ, 'LineWidth', lw, 'TickDir', 'out','ColorOrder', cmp2)
                title(sprintf('Adapter SOA = %.1f s', adaptor_soa(adapter)./1e3),'FontSize',fontSZ,'FontWeight','bold')

                % data
                scatter(D(ss).data(adapter).pre_ms_unique, D(ss).data(adapter).pre_pResp, dotSZ, 'filled');

                % prediction
                plot(pred{mm,ss}.test_soa, pred{mm,ss}.pre_pmf, 'LineWidth',lw)

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
                plot(pred{mm,ss}.test_soa, squeeze(pred{mm,ss}.post_pmf(adapter,:,:)),'LineWidth',lw)
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
                saveas(gca, fullfile(out_dir, flnm),'png')
            end

        end

    end

end