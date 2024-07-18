% fig S13: data analysis of the oddball-detection task

clear; clc; close all;

%% Manage Paths

restoredefaultpath;
current_dir = pwd;
[project_dir, ~] = fileparts(current_dir);
addpath(genpath(fullfile(project_dir, 'data')));
out_dir = fullfile(current_dir, mfilename);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% Set Up

subject_list = [1:4, 6:10]; % Selected subjects
num_subjects = 10;
num_sessions = 9;

[hit_rate, false_alarm] = deal(NaN(2, num_subjects, 2)); % Preallocate Hit Rate and False Alarm arrays

%% Load and Process Data

for subj = subject_list
    % Preallocate response and oddball arrays
    [expo_resp_aud, expo_oddball_aud, expo_resp_vis, expo_oddball_vis] = deal([]);
    [post_resp_aud, post_oddball_aud, post_resp_vis, post_oddball_vis] = deal([]);
    
    for sess = 1:num_sessions
        % Load exposure session data
        load(sprintf('exposure_sub%i_session%i.mat', subj, sess));
        [resp_aud, oddball_aud, resp_vis, oddball_vis] = summarize_session_performance(1, ExpInfo.nTTrials, ExpInfo.idxOddballA, ExpInfo.idxOddballV, Response.oddball);
        expo_resp_aud = [expo_resp_aud, resp_aud];
        expo_oddball_aud = [expo_oddball_aud, oddball_aud];
        expo_resp_vis = [expo_resp_vis, resp_vis];
        expo_oddball_vis = [expo_oddball_vis, oddball_vis];
        
        % Load post-test session data
        load(sprintf('posttest_sub%i_session%i.mat', subj, sess));
        [resp_aud, oddball_aud, resp_vis, oddball_vis] = summarize_session_performance(1, ExpInfo.expoNTTrials, ExpInfo.idxOddballA, ExpInfo.idxOddballV, ExpoResponse.oddball);
        post_resp_aud = [post_resp_aud, resp_aud];
        post_oddball_aud = [post_oddball_aud, oddball_aud];
        post_resp_vis = [post_resp_vis, resp_vis];
        post_oddball_vis = [post_oddball_vis, oddball_vis];
    end

    % Exclude trials with both oddballs present
    expo_exclude_idx = expo_oddball_aud & expo_oddball_vis;
    post_exclude_idx = post_oddball_aud & post_oddball_vis;
    
    % Remove excluded trials
    [expo_resp_aud, expo_oddball_aud, expo_resp_vis, expo_oddball_vis] = remove_excluded_trials(expo_resp_aud, expo_oddball_aud, expo_resp_vis, expo_oddball_vis, expo_exclude_idx);
    [post_resp_aud, post_oddball_aud, post_resp_vis, post_oddball_vis] = remove_excluded_trials(post_resp_aud, post_oddball_aud, post_resp_vis, post_oddball_vis, post_exclude_idx);
    
    % Calculate Hit Rates (HR) and False Alarms (FA)
    hit_rate(:, subj, :) = [sum(expo_oddball_aud & expo_resp_aud) / sum(expo_oddball_aud), sum(expo_oddball_vis & expo_resp_vis) / sum(expo_oddball_vis); ...
                            sum(post_oddball_aud & post_resp_aud) / sum(post_oddball_aud), sum(post_oddball_vis & post_resp_vis) / sum(post_oddball_vis)];
    false_alarm(:, subj, :) = [sum(~expo_oddball_aud & expo_resp_aud) / sum(~expo_oddball_aud), sum(~expo_oddball_vis & expo_resp_vis) / sum(~expo_oddball_vis); ...
                               sum(~post_oddball_aud & post_resp_aud) / sum(~post_oddball_aud), sum(~post_oddball_vis & post_resp_vis) / sum(~post_oddball_vis)];
end

%% Calculate d-prime and Criterion

dprime = NaN(2, numel(subject_list), 2);
criteria = NaN(2, numel(subject_list), 2);

for task = 1:2
    for idx = 1:numel(subject_list)
        subj = subject_list(idx);
        for av = 1:2 % 1 = Auditory; 2 = Visual
            [dprime(task, idx, av), criteria(task, idx, av)] = calculate_dprime(hit_rate(task, subj, av), false_alarm(task, subj, av), [numel(expo_resp_aud), numel(post_resp_aud)]);
        end
    end
end

%% Display Stats

fprintf('[%s] Mean of auditory d-prime %.2f and visual d-prime %.2f\n', mfilename, mean(dprime(:, :, 1), 'all'), mean(dprime(:, :, 2), 'all'));
re_dprime = reshape(dprime, [2 * numel(subject_list), 2]);
fprintf('[%s] S.D. of auditory d-prime %.2f and visual d-prime %.2f\n', mfilename, std(re_dprime(:, 1)), std(re_dprime(:, 2)));

%% Plot

clt = [202, 0, 32; 5, 113, 176] ./ 255; % red and blue
lw = 1.5;
font_sz = 15;
title_sz = 20;
dot_sz = 120;
label = {'auditory'; 'visual'};

figure;
set(gcf, 'Position', [10 10 450 400]);
axis equal;
set(gca, 'LineWidth', lw, 'FontSize', font_sz, 'TickDir', 'out');
hold on;
scatter(dprime(1, :, 1), dprime(2, :, 1), dot_sz, 'filled', 'MarkerFaceColor', clt(1, :), 'MarkerEdgeColor', 'w', 'LineWidth', lw);
scatter(dprime(1, :, 2), dprime(2, :, 2), dot_sz, 'filled', 'MarkerFaceColor', clt(2, :), 'MarkerEdgeColor', 'w', 'LineWidth', lw);

lb = 1;
ub = 4.5;
x = linspace(lb, ub, 10);
plot(x, x, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

xlim([lb, ub]);
xticks(lb:ub);
ylim([lb, ub]);
yticks(lb:ub);

xlabel('d'', exposure phase')%, 'Interpreter', 'latex');
ylabel('d'', post-test phase')%, 'Interpreter', 'latex');

legend({'Audition', 'Vision'}, 'Location', 'southeast');

flnm = 'oddball_result';
saveas(gca, fullfile(out_dir, flnm), 'pdf');

%% Helper Functions

function [d, c] = calculate_dprime(hit_rate, false_alarm, n_trial)
    iHR = hit_rate;
    iFA = false_alarm;
    if hit_rate == 1
        iHR = 1 - 0.5 / n_trial; % or 1 - 0.5(hit_rate + false_alarm)
    elseif hit_rate == 0
        iHR = 0.5 / n_trial;
    elseif false_alarm == 1
        iFA = 1 - 0.5 / n_trial;
    elseif false_alarm == 0
        iFA = 0.5 / n_trial;
    end
    d = norminv(iHR, 0, 1) - norminv(iFA, 0, 1);
    c = norminv(1 - iFA, 0, 1);
end

function [resp_aud, oddball_aud, resp_vis, oddball_vis] = summarize_session_performance(start_trial, end_trial, idx_oddball_aud, idx_oddball_vis, response)
    oddball_aud = zeros(1, length(response));
    oddball_aud(idx_oddball_aud) = 1; % oddball presence
    resp_aud = ismember(response, [2, 3]); % oddball response
    oddball_aud = oddball_aud(start_trial:end_trial);
    resp_aud = resp_aud(start_trial:end_trial);
    
    oddball_vis = zeros(1, length(response));
    oddball_vis(idx_oddball_vis) = 1; % oddball presence
    resp_vis = ismember(response, [2, 1]); % oddball response
    oddball_vis = oddball_vis(start_trial:end_trial);
    resp_vis = resp_vis(start_trial:end_trial);
end

function [resp_aud, oddball_aud, resp_vis, oddball_vis] = remove_excluded_trials(resp_aud, oddball_aud, resp_vis, oddball_vis, exclude_idx)
    resp_aud(exclude_idx) = [];
    oddball_aud(exclude_idx) = [];
    resp_vis(exclude_idx) = [];
    oddball_vis(exclude_idx) = [];
end