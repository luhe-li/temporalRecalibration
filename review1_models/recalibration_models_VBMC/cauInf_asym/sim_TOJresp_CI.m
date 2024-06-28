function response = sim_TOJresp_CI(tau, sigma_a, sigma_v, criterion, lambda, p_common, model)

%% experimental info

n_test_soa = numel(model.test_soa);
n_trial = model.n_trial;

%% sample noisy measurements

% draw y from U[0, 1] for each test soa
y = rand(n_test_soa, n_trial);
bool_y = y < model.y_criterion;
soa_m = NaN(size(y));
m_test_soa = repmat(model.test_soa', 1, n_trial);
m_tau = repmat(tau, n_test_soa, n_trial);

% sample from iCDF of exponential measurement distribution
% For y smaller than threshold:
soa_m(bool_y)  = m_test_soa(bool_y) + m_tau(bool_y) + sigma_v .* log((sigma_a + sigma_v)./ sigma_v.* y(bool_y));
% for y larger than threshod:
soa_m(~bool_y) = m_test_soa(~bool_y) + m_tau(~bool_y) - sigma_a .* log((sigma_a + sigma_v)./ sigma_a .* (1 - y(~bool_y)));

%% perform causal inference on measurements

v_soa_m = reshape(soa_m, [1, numel(soa_m)]);
[shat, ~] = cau_inf_map(v_soa_m, p_common, model);

%% compare to criterion
response = ones(size(shat))*2; % report simultaneuous
response(shat >= criterion) = 1; % report V
response(shat < -criterion) = 3; % report A

%% add lapse

lapse_idx = rand(size(response)) < lambda;
response(lapse_idx) = randi([1, 3], size(response(lapse_idx)));

end

