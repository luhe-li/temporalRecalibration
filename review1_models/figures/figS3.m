% fig S3: plot confidence interval of asymmetry index

clear; clc; close all;

%% manage paths

restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'utils')));
out_dir = fullfile(currentDir, mfilename);
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% set up

sub_slc = [1:4,6:10];
adp_slc = 1:9;

%% extract delta_mu and calculate asymmetry index

for i = 1:numel(sub_slc)

    sub = sub_slc(i);
    flnm = sprintf('btstResults_sub%i_', sub);
    allFiles           = dir(fullfile(project_dir, TOJ_dir,'**', [flnm '*.mat']));
    load(allFiles.name);

    btst_ai(i,:) = sum(reshape([btst.pred.all_pss_shift],[9,1000]),1);
    [CI_lb(i), CI_ub(i)]  = get95CI(btst_ai(i,:));
end

%% plot

lw = 0.5;
fontsz = 7;
titleFontSz = 10;

figure;
set(gcf, 'Position', [0,0,210,180]);
set(gca, 'LineWidth', lw, 'FontSize', fontsz,'TickDir', 'out');
hold on 

for i = 1:numel(sub_slc)

plot([i,i],  [CI_lb(i), CI_ub(i)],'k','LineWidth',2)
   
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

%% function
function [CI_lb, CI_ub] = get95CI(v)
%first sort the vector in an ascending order
n_sorted = sort(v);
%compute how many entries the vector has
lenN     = length(v);
%lower bound
CI_ub    = n_sorted(ceil(lenN*0.975));
%upper bound
CI_lb    = n_sorted(floor(lenN*0.025));
end