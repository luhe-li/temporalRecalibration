function out = sim_cauInf_asym(param, T, model)

% response is a vector of [session[pretest[test_soa[trial]]]] [session[posttest[test_soa[trial]]]]
if ~exist("T",'var'); T=1:numel(model.vec_resp); end

if strcmp(model.mode, 'initialize')

    out.paraID   = {'\tau','\sigma_{A}','\sigma_{V}','c','\lambda','p_{common}','\alpha','\sigma_{C=1}','\sigma_{C=2}'};
    out.num_para = length(out.paraID);

    % hard bounds, the range for LB, UB, larger than soft bounds
    paraH.tau      = [-100,   100]; % ms
    paraH.sigma_a  = [  10,   120]; % ms
    paraH.sigma_v  = [  10,   200]; % ms
    paraH.criterion= [   1,   350]; % criterion, s
    paraH.lambda   = [1e-4,  0.06]; % percentage
    paraH.p_common = [1e-4, 1-1e-4]; % weight
    paraH.alpha    = [1e-4,  0.02]; % percentage
    paraH.sigma_C1 = [   1,   300]; % ms
    paraH.sigma_C2 = [ 100,   1e3]; % ms

    % soft bounds, the range for PLB, PUB
    paraS.tau      = [ -50,    50]; % ms
    paraS.sigma_a  = [  20,    50]; % ms
    paraS.sigma_v  = [  20,   120]; % ms
    paraS.criterion= [  30,   150]; % criterion, s
    paraS.lambda   = [0.01,  0.03]; % percentage
    paraS.p_common = [ 0.3,   0.7]; % weight
    paraS.alpha    = [1e-3,  2e-3]; % percentage
    paraS.sigma_C1 = [  10,   100]; % ms
    paraS.sigma_C2 = [ 500,   700]; % ms

    % reorganize parameter bounds to feed to bads
    fn = fieldnames(paraH);
    for k= 1:numel(fn)
        out.lb(:,k)  = paraH.(fn{k})(1);
        out.ub(:,k)  = paraH.(fn{k})(2);
        out.plb(:,k) = paraS.(fn{k})(1);
        out.pub(:,k) = paraS.(fn{k})(2);
    end

    % get grid initializations
    numSections = model.num_runs * 2;
    out.init = getInit(out.plb, out.pub, numSections, model.num_runs);

else

    %% free parameters

    tau = param(1);
    sigma_a = param(2);
    sigma_v = param(3);
    criterion = param(4);
    lambda = param(5);
    p_common = param(6);
    alpha = param(7);
    sigma_C1 = param(8);
    sigma_C2 = param(9);

    %% precompute

    % prior
    model.x_axis     = -model.bound_full:1:model.bound_full;
    model.x_axis_int = -model.bound_int:1:model.bound_int;
    model.l_window   = find(model.x_axis == model.x_axis_int(1));
    model.r_window   = find(model.x_axis == model.x_axis_int(end));
    model.prior_C1   = normpdf(model.x_axis_int, 0, sigma_C1);
    model.prior_C2   = normpdf(model.x_axis_int, 0, sigma_C2);

    % likelihood that centers around 0
    idx_peak   = ceil(length(model.x_axis)/2);
    lf = (1/(sigma_a + sigma_v)).* exp(1/sigma_a .* (model.x_axis(1:idx_peak)));
    rf = (1/(sigma_a + sigma_v)).* exp(-1/sigma_v .* (model.x_axis(idx_peak+1:end)));
    model.df_likelihood     = [lf, rf] + realmin;
    model.y_criterion = sigma_v/(sigma_a + sigma_v);

    %% simulate data for pre-test

    [pre_resp, post_resp] = deal(NaN(model.num_ses, numel(model.test_soa) * model.n_trial));

    % simulate recalibration of tau
    tau_shift = sim_recal_CI(tau, sigma_a, sigma_v, p_common, alpha, model);

    for ses = 1:model.num_ses

        pre_resp(ses,:) = sim_TOJresp_CI(tau, sigma_a, sigma_v, criterion, lambda, p_common, model);
        post_resp(ses,:) = sim_TOJresp_CI(tau + tau_shift(ses), sigma_a, sigma_v, criterion, lambda, p_common, model);

    end

    resp = [reshape(pre_resp, 1, numel(pre_resp)), reshape(post_resp, 1, numel(post_resp))]';
    out = resp(T);
end

end