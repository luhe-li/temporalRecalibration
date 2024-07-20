clear; close all;

%% manage path

cur_dir               = pwd;
[project_dir, ~]      = fileparts(cur_dir);
[git_dir, ~] = fileparts(project_dir);
dataDir = fullfile(git_dir,'temporalRecalibrationData');
addpath(genpath(fullfile(project_dir, 'utils')));
addpath(genpath(fullfile(project_dir, 'vbmc')));
out_dir               = fullfile(cur_dir, mfilename);
if ~exist(out_dir,'dir') mkdir(out_dir); end

%% load results

results_folder = fullfile(dataDir,'recalibration_models_VBMC','sim_tradeoff_by_NLL');
files = dir(fullfile(results_folder, 'sim_num-*'));
pattern = 'sim_num\d{3}_i(\d{2})_j(\d{2})';

for ii = 1:size(files)
    flnm = files(ii).name;
    tokens = regexp(filename, pattern, 'tokens');
    ii = str2double(tokens{1}{1});
    jj = str2double(tokens{1}{2});

    r = load(fullfile(results_folder, files(ii).name));
    tempNLL = r.NLL;
    NLL(:,:,ii,jj) = tempNLL;

end

save('sim_tradeoff',"NLL");

%% plot

set(gcf, 'Position', get(0, 'Screensize'));
set(gca, 'LineWidth', 1.5, 'FontSize', 15, 'TickDir', 'out')

for pp = 1:size(NLL, 1)

    figure
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gca, 'LineWidth', 1.5, 'FontSize', 15, 'TickDir', 'out')

    for ss           = 1:size(NLL,2)

        subplot(3, 3, ss);

        d = squeeze(nll(pp, ss, :, :));
        imagesc(p2s, p1s, d, [min(d, [], "all"), min(d, [], "all")+4000]); hold on
        c                = colorbar;
        c.Label.String   = 'NLL';
        %         c.Label.FontSize = 10;

        xticks(p2s)
        yticks(p1s)
        xlabel(fits(sub).paraID(idx_p2))
        ylabel(fits(sub).paraID(idx_p1))
        title(sprintf('S%i', sub_slc(ss)))

    end

    if save_fig
        flnm = sprintf('tradeoff_comb %i', pp);
        saveas(gca, fullfile(outDir, flnm),'png')
    end
end