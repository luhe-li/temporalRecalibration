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

save_fig=1;

results_folder = fullfile(pwd,'sim_tradeoff_by_NLL');
files = dir(fullfile(results_folder, 'sim_num*'));
pattern = 'sim_num\d{3}_i(\d{2})_j(\d{2})';

for uu = 1:size(files)
    
    flnm = files(uu).name;
    tokens = regexp(flnm, pattern, 'tokens');
    ii = str2double(tokens{1}{1});
    jj = str2double(tokens{1}{2});

    r = load(fullfile(results_folder, files(uu).name));
    tempNLL = r.NLL;
    NLL(:,:,ii,jj) = tempNLL;

end

para_combi = r.para_combi;
paraID = r.paraID;
p1ss = r.p1ss;
p2ss = r.p2ss;

%% plot

for pp = 1:size(NLL, 1)

    figure
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gca, 'LineWidth', 1.5, 'FontSize', 15, 'TickDir', 'out')

    for ss           = 1:size(NLL,2)

        subplot(3, 3, ss);

        d = squeeze(NLL(pp, ss, :, :));
        imagesc(p2ss(pp,:), p1ss(pp,:), d); %[min(d, [], "all"), min(d, [], "all")]); 
        hold on
        c                = colorbar;
        c.Label.String   = 'NLL';
        c.Label.FontSize = 10;

        xticks(p2s)
        yticks(p1s)
        xlabel(paraID(para_combi(pp,2)))
        ylabel(paraID(para_combi(pp,1)))
        title(sprintf('S%i', (ss)))

    end

    if save_fig
        flnm = sprintf('tradeoff_comb %i', pp);
        saveas(gca, fullfile(out_dir, flnm),'png')
    end
end