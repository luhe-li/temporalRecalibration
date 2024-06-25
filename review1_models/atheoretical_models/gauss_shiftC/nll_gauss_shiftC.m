
function out         = nll_gauss_shiftC(mu, sigma, c_pre, dc_post1, dc_post2, dc_post3, dc_post4, dc_post5, dc_post6, dc_post7, dc_post8, dc_post9, lambda, ...
    model, data)

%--------------------------------------------------------------------------
% Inputs:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% out:

% if model.mode     = 'initialize'
% nLL               : total nLL of pre and posttest, all sessions

% if model.mode     = 'optimize'
% nLL               : total nLL of pre and posttest, all sessions

% if model.mode     = 'predict'
% out.pre_pmf       : psychometric function of pre-test
% out.all_pss_shift : simulated all mu_shift
% out.R2            : R^2 of using Gaussian to approxiamte mu_shift pdf
% out.post_pmf      : psychometric function of post-test;

%--------------------------------------------------------------------------

%% pre-compute

dc_post = [dc_post1, dc_post2, dc_post3, dc_post4, dc_post5, dc_post6, dc_post7, dc_post8, dc_post9];

%% loop for each session

if strcmp(model.mode, 'optimize')

    % reorder and fit by adaptor order
    adaptor_soas             = [data(1:model.num_ses).adaptor_soa];
    [~, order]               = sort(adaptor_soas); 

    % find lower and upper bounds for simultaneity window in pre-test
    lc_pre = mu - c_pre;
    uc_pre = mu + c_pre;

    [pre_afirst, pre_simul, pre_vfirst] = pmf_gauss(model.test_soa,...
        mu, sigma, lc_pre, uc_pre, lambda);

    nLL_ses = NaN(1, model.num_ses);
    for i = 1:model.num_ses

        adaptor = order(i);
        adaptor_soa = model.sim_adaptor_soa(i);

        pre_LL = data(adaptor).pre_nT_A1st*log(pre_afirst)'...
            + data(adaptor).pre_nT_V1st*log(pre_vfirst)'...
            + data(adaptor).pre_nT_simul*log(pre_simul)';

        if adaptor_soa > uc_pre
            uc_post = uc_pre + dc_post(i);
            lc_post = lc_pre;
        elseif adaptor_soa < lc_pre
            lc_post = lc_pre - dc_post(i);
            uc_post = uc_pre;
        else
            uc_post = uc_pre;
            lc_post = lc_pre;
        end

        [post_afirst, post_simul, post_vfirst] = pmf_gauss(model.test_soa,...
            mu, sigma, lc_post, uc_post, lambda);

        post_LL = data(adaptor).post_nT_A1st*log(post_afirst)'...
            + data(adaptor).post_nT_V1st*log(post_vfirst)'...
            + data(adaptor).post_nT_simul*log(post_simul)';

        nLL_ses(adaptor)         = - pre_LL - post_LL;

    end
    out                  = nansum(nLL_ses);

elseif strcmp(model.mode,'predict')

    out.test_soa         = model.test_soa;
    if model.test_axis_finer == 1
        out.test_soa = linspace(out.test_soa(1), out.test_soa(end), 1e3);
    end
    out.adaptor_soa      = model.sim_adaptor_soa;

    lc_pre = mu - c_pre;
    uc_pre = mu + c_pre;

    [pre_afirst, pre_simul, pre_vfirst] = pmf_gauss(model.test_soa,...
        mu, sigma, lc_pre, uc_pre, lambda);

    c_shift = NaN(1, model.num_ses);
    for adaptor = 1:9
        out.pre_pmf{adaptor}     = [pre_vfirst; pre_simul; pre_afirst];
        adaptor_soa = model.sim_adaptor_soa(adaptor);

        if adaptor_soa > uc_pre
            uc_post = uc_pre + dc_post(adaptor);
            lc_post = lc_pre;
            c_shift(adaptor_soa) = dc_post(adaptor);
        elseif adaptor_soa < lc_pre
            lc_post = lc_pre - dc_post(adaptor);
            uc_post = uc_pre;
            c_shift(adaptor_soa) = -dc_post(adaptor);
        else
            uc_post = uc_pre;
            lc_post = lc_pre;
            c_shift(adaptor_soa) = 0;
        end

        [post_afirst, post_simul, post_vfirst] = pmf_gauss(model.test_soa,...
            mu, sigma, lc_post, uc_post, lambda);

        out.post_pmf{adaptor}    = [post_vfirst; post_simul; post_afirst];

    end

    out.pss_shift = zeros(1:model.num_ses);
    out.c_shift = c_shift;

end

end
