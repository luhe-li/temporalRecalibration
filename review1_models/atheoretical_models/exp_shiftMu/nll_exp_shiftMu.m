function out         = nll_exp_shiftMu(tau_pre, tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9,...
    sigma_a, sigma_v, criterion, lambda, ...
    model, data)

%% pre-compute

tau_post = [tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9];

%% loop for each session

if strcmp(model.mode, 'optimize')

    % fit by adaptor order
    adaptor_soas             = [data(1:model.num_ses).adaptor_soa];
    [~, order]               = sort(adaptor_soas);

    [pre_afirst, pre_simul, pre_vfirst] = pmf_exp(model.test_soa,...
        tau_pre, sigma_a, sigma_v, -criterion, criterion, lambda);

    nLL_ses = NaN(1, model.num_ses);
    for i = 1:model.num_ses

        adaptor = order(i);

        pre_LL = data(adaptor).pre_nT_A1st*log(pre_afirst)'...
            + data(adaptor).pre_nT_V1st*log(pre_vfirst)'...
            + data(adaptor).pre_nT_simul*log(pre_simul)';

        [post_afirst, post_simul, post_vfirst] = pmf_exp(model.test_soa,...
            tau_post(i), sigma_a, sigma_v, -criterion, criterion, lambda);

        post_LL = data(adaptor).post_nT_A1st*log(post_afirst)'...
            + data(adaptor).post_nT_V1st*log(post_vfirst)'...
            + data(adaptor).post_nT_simul*log(post_simul)';

        % sum the negative likelihood of pre and post test
        nLL_ses(adaptor)         = - pre_LL - post_LL;

    end
    out                  = nansum(nLL_ses);

elseif strcmp(model.mode,'predict')

    out.test_soa         = model.test_soa;
    if model.test_axis_finer == 1
        out.test_soa = linspace(out.test_soa(1), out.test_soa(end), 1e3);
    end
    out.adaptor_soa      = model.sim_adaptor_soa;

    % pre-test pmf
    [pre_afirst, pre_simul, pre_vfirst] = pmf_exp(out.test_soa,...
        tau_pre, sigma_a, sigma_v, -criterion, criterion, lambda);

    for adaptor = 1:9
        out.pre_pmf{adaptor}     = [pre_vfirst; pre_simul; pre_afirst];

        % post-test pmf
        [post_afirst, post_simul, post_vfirst] = pmf_exp(out.test_soa,...
            tau_post(adaptor), sigma_a, sigma_v, -criterion, criterion, lambda);
        out.post_pmf{adaptor}    = [post_vfirst; post_simul; post_afirst];

    end

    % pss shift = - tau shift
    out.all_pss_shift = - (tau_post - tau_pre);

end

end
