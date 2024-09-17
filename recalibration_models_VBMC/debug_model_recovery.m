% use all the codes in model recovery to see if fixed-update model can be
% recovered by current fitting setup
clear; close all; clc;

model_slc = 5;

%% select models

rng('shuffle'); rng('Shuffle');
specifications = {'Heuristic, asymmetric', 'Heuristic, symmetric', 'Causal inference, asymmetric', 'Causal inference, symmetric','Fixed updated, asymmetric', 'Fixed updated, symmetric'};
folders = {'heu_asym', 'heu_sym', 'cauInf_asym', 'cauInf_sym','fixed_asym','fixed_sym'};
numbers = (1:numel(specifications))';
model_info = table(numbers, specifications', folders', 'VariableNames', {'Number', 'Specification', 'FolderName'});

%% manage paths

% restoredefaultpath;
currentDir= pwd;
[projectDir, ~]= fileparts(currentDir);
[tempDir, ~] = fileparts(projectDir);
dataDir = fullfile(tempDir,'temporalRecalibrationData');
addpath(genpath(fullfile(projectDir, 'data')));
addpath(genpath(fullfile(projectDir, 'vbmc')));
addpath(genpath(fullfile(projectDir, 'utils')));
outDir = fullfile(currentDir, mfilename);
if ~exist(outDir, 'dir'); mkdir(outDir); end
useCluster = false;
if useCluster == false; projectDir = dataDir; end

%% organize data

sub_slc = [1:4,6:10];
for sess = 1:9
    data(sess) = organizeData(1, sess);
end

%% define model

% set fixed & set-up parameters
model.num_ses = 9;
model.thres_R2 = 0.95;
model.expo_num_sim = 1e3; % number of simulation for exposure phase
model.expo_num_trial = 250; % number of *real* trials in exposure phase
model.num_runs = 10;
model.num_bin  = 100; % numer of bin to approximate tau_shift distribution
model.bound_full = 10*1e3; % in second, the bound for prior axis
model.bound_int = 1.4*1e3; % in second, where measurements are likely to reside
model.num_sample = 1e3; % number of samples for simulating psychometric function with causal inference, only used in pmf_exp_CI
model.test_soa = [-0.5, -0.3:0.05:0.3, 0.5]*1e3;
model.sim_adaptor_soa  = [-0.7, -0.3:0.1:0.3, 0.7]*1e3;
model.toj_axis_finer = 0; % simulate pmf with finer axis
model.adaptor_axis_finer = 0; % simulate with more adpators

%% sample ground-truth from best parameter estimates

i_sample = 1;

for sim_m = model_slc

    sim_str = folders{sim_m};
    sim_func = str2func(['nll_' sim_str]);
    addpath(genpath(fullfile(pwd, sim_str)));

    model.mode = 'initialize';
    Val = sim_func([], model, []);

    % load best estimates
    clearvars bestP
    result_folder = fullfile(projectDir, 'recalibration_models_VBMC', sim_str);
    R = load_subject_data(result_folder, sub_slc, 'sub-*');
    for ss = 1:numel(sub_slc)
        bestP(ss,:) = R{ss}.diag.post_mean;
    end

    % sample 100 ground truth
    mu_GT = mean(bestP, 1);
    sd_GT = std(bestP, [], 1);
    GT_samples = generate_samples(Val, mu_GT, sd_GT, 1);
    model.mode = 'predict';

    temp_sim_func = sim_func;
    pred =  temp_sim_func(GT_samples, model, []);
    sim_data = simulateData(pred, data);
    fake_data(sim_m, i_sample).data = sim_data;
    fake_data(sim_m, i_sample).gt_p = GT_samples;
    fake_data(sim_m, i_sample).mu_GT = mu_GT;
    fake_data(sim_m, i_sample).sd_GT = sd_GT;
    fake_data(sim_m, i_sample).pred = pred;

end

%% fitting

% set OPTIONS
options = vbmc('defaults');
options.TolStableCount = 10;
options.SpecifyTargetNoise = true;

sim_data = fake_data(model_slc, i_sample).data;
currModelStr = model_info.FolderName{model_slc};
currModel = str2func(['nll_' currModelStr]);

model.mode = 'initialize';
Val = currModel([], model, []);
model.initVal = Val;

% set priors
lpriorfun = @(x) msplinetrapezlogpdf(x, Val.lb, Val.plb, Val.pub, Val.ub);
% lpriorfun = @(x) msplinetrapezlogpdf(x, mu_GT - 3*sd_GT, mu_GT - sd_GT,  mu_GT + sd_GT, mu_GT + 3*sd_GT);

% set likelihood
model.mode = 'optimize';
llfun = @(x) currModel(x, model, sim_data);
% fun = @(x) llfun(x) + lpriorfun(x);
fun = @(x) lpostfun(x,llfun,lpriorfun); 

fprintf('[%s] Start sim model-%s, fit model-%s, recovery sample-%i \n', mfilename, folders{sim_m}, currModelStr, i_sample);

% vp: variational posterior
% elbo: Variational Evidence Lower d Bound
[temp_vp, temp_elbo, temp_elbo_sd, temp_exitflag,temp_output] = vbmc(fun, Val.init(randi(model.num_runs,1),:), Val.lb,...
    Val.ub, Val.plb, Val.pub, options);

% make predictions
Xs = vbmc_rnd(temp_vp,1e5);
post_mean = mean(Xs,1);
fit_est = post_mean;

%%
model.mode       = 'predict';
pred =  currModel(post_mean, model, sim_data);

%% test 1: compare estimates

disp(fake_data(model_slc, i_sample).gt_p')
disp(fit_est)

% Estimates are close to the ground truth.

%% test 2: see prediction

figure; hold on
plot(fake_data(model_slc, i_sample).pred.adaptor_soa, mean(fake_data(model_slc, i_sample).pred.pss_shift, 2),'-ok')
plot(pred.adaptor_soa, mean(pred.pss_shift,2),'--r')

% Test if LL of GT > LL of fitted estimates
GT = fake_data(model_slc, i_sample).gt_p';
EST = fit_est; 

ll_gt = llfun(GT);
ll_est = llfun(EST);
% this should be one, since estimation's negative log likelihood should be the smallest
-ll_gt>-ll_est 
p_gt = fun(GT);
p_est = fun(EST);
% this should be one, since estimation's log posterior should be the largest
p_gt < p_est

%% test 3: see how sigma and alpha influences NLL

model.mode       = 'predict';
adaptor_soa = pred.adaptor_soa;
sigmas = 40:10:80;
alphas = 0.0001:0.0001:0.0004;
figure;
set(gcf, 'Position', get(0, 'Screensize'));

for ii = 1:numel(sigmas)
for jj = 1:numel(alphas)
% 
    test_p = [41.1716, sigmas(ii), 41.0812, 240.9095, 0.0372, alphas(jj)];
    ll(ii,jj) = llfun(test_p);
    pred =  currModel(test_p, model, sim_data);
    psss(ii,jj,:) = mean(pred.pss_shift,2);

    subplot(numel(sigmas), numel(alphas), (ii-1)*numel(alphas) + jj)
    hold on
    plot(fake_data(model_slc, i_sample).pred.adaptor_soa, mean(fake_data(model_slc, i_sample).pred.pss_shift, 2),'-ok')
    plot(adaptor_soa, squeeze(psss(ii,jj,:)),'--r')
    ylim([-200 200])
    xlim([-700 700])
    title(sprintf('LL %.2f, sigma %.2f, alpha %.4f', ll(ii,jj), sigmas(ii), alphas(jj)))
end
end

saveas(gca, fullfile(outDir,'NLL test'),'png')
