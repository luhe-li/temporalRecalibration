
function out         = nll_gauss_shiftMu(mu, sigma, c, mu_post1, mu_post2, mu_post3, mu_post4, mu_post5, mu_post6, mu_post7, mu_post8, mu_post9, lambda, ...
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

mu_post = [mu_post1, mu_post2, mu_post3, mu_post4, mu_post5, mu_post6, mu_post7, mu_post8, mu_post9];

%% loop for each session

if strcmp(model.mode, 'optimize')

    % reorder and fit by adaptor order
    adaptor_soas             = [data(1:model.num_ses).adaptor_soa];
    [~, order]               = sort(adaptor_soas); 

    [pre_afirst, pre_simul, pre_vfirst] = pmf_gauss(model.test_soa,...
        mu, sigma, mu - c, mu + c, lambda);

    nLL_ses = NaN(1, model.num_ses);
    for i = 1:model.num_ses

        adaptor = order(i);

        pre_LL = data(adaptor).pre_nT_A1st*log(pre_afirst)'...
            + data(adaptor).pre_nT_V1st*log(pre_vfirst)'...
            + data(adaptor).pre_nT_simul*log(pre_simul)';

        [post_afirst, post_simul, post_vfirst] = pmf_gauss(model.test_soa,...
            mu_post(adaptor), sigma, mu_post(adaptor) - c, mu_post(adaptor) + c, lambda);

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

    [pre_afirst, pre_simul, pre_vfirst] = pmf_gauss(model.test_soa,...
        mu, sigma, mu - c, mu + c, lambda);

    for adaptor = 1:9
        
        out.pre_pmf{adaptor}     = [pre_vfirst; pre_simul; pre_afirst];

        [post_afirst, post_simul, post_vfirst] = pmf_gauss(model.test_soa,...
             mu_post(adaptor), sigma, mu_post(adaptor) - c, mu_post(adaptor) + c, lambda);
        out.post_pmf{adaptor}    = [post_vfirst; post_simul; post_afirst];

    end

    out.pss_shift = mu_post - mu;
    out.c_shift = zeros(1, model.num_ses);

end

end
