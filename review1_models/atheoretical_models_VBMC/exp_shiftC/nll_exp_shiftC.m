function out         = nll_exp_shiftC(freeParam, model, data)

if strcmp(model.mode, 'initialize')

    out.paraID     = {'\mu','\sigma_A','\sigma_V','c1_{pre}','dc_{post,1}','dc_{post,2}','dc_{post,3}','dc_{post,4}','dc_{post,5}','dc_{post,6}','dc_{post,7}','dc_{post,8}','dc_{post,9}','\lambda'};
    out.num_para   = length(out.paraID);

    % hard bounds, the range for LB, UB, larger than soft bounds
    paraH.mu         = [-300,   300]; % ms
    paraH.sigma_a    = [  10,   200]; % ms
    paraH.sigma_v    = [  10,   200]; % ms
    paraH.c1         = [   1,   350]; % ms
    paraH.dc_post1   = [0.01,   300]; % ms
    paraH.dc_post2   = [0.01,   300]; % ms
    paraH.dc_post3   = [0.01,   300]; % ms
    paraH.dc_post4   = [0.01,   300]; % ms
    paraH.dc_post5   = [0.01,   300]; % ms
    paraH.dc_post6   = [0.01,   300]; % ms
    paraH.dc_post7   = [0.01,   300]; % ms
    paraH.dc_post8   = [0.01,   300]; % ms
    paraH.dc_post9   = [0.01,   300]; % ms
    paraH.lambda     = [1e-4,  0.06]; % percentage

    % soft bounds, the range for PLB, PUB
    paraS.mu         = [-100,   100]; % ms
    paraS.sigma_a    = [  20,   120]; % ms
    paraS.sigma_v    = [  30,   150]; % ms
    paraS.c1         = [   1,   100]; % ms
    paraS.dc_post1   = [   1,   100]; % ms
    paraS.dc_post2   = [   1,   100]; % ms
    paraS.dc_post3   = [   1,   100]; % ms
    paraS.dc_post4   = [   1,   100]; % ms
    paraS.dc_post5   = [   1,   100]; % ms
    paraS.dc_post6   = [   1,   100]; % ms
    paraS.dc_post7   = [   1,   100]; % ms
    paraS.dc_post8   = [   1,   100]; % ms
    paraS.dc_post9   = [   1,   100]; % ms
    paraS.lambda     = [0.01,  0.03]; % percentage

    % reorganize parameter bounds to feed to bads
    fn    = fieldnames(paraH);
    for k = 1:numel(fn)
        out.lb(:,k)    = paraH.(fn{k})(1);
        out.ub(:,k)    = paraH.(fn{k})(2);
        out.plb(:,k)   = paraS.(fn{k})(1);
        out.pub(:,k)   = paraS.(fn{k})(2);
    end
    model.paraS      = paraS; model.paraH = paraH;

    % get grid initializations
    numSections = model.num_runs * 2;
    out.init = getInit(out.plb, out.pub, numSections, model.num_runs);

else

    %% assign free parameters

    tau = freeParam(1);
    sigma_a = freeParam(2);
    sigma_v = freeParam(3);
    c_pre = freeParam(4);
    dc_post1 = freeParam(5);
    dc_post2 = freeParam(6);
    dc_post3 = freeParam(7);
    dc_post4 = freeParam(8);
    dc_post5 = freeParam(9);
    dc_post6 = freeParam(10);
    dc_post7 = freeParam(11);
    dc_post8 = freeParam(12);
    dc_post9 = freeParam(13);
    lambda = freeParam(14);

    %% pre-compute

    dc_post = [dc_post1, dc_post2, dc_post3, dc_post4, dc_post5, dc_post6, dc_post7, dc_post8, dc_post9];

    %% loop for each session

    if strcmp(model.mode, 'optimize')

        % reorder and fit by adaptor order
        adaptor_soas             = [data(1:model.num_ses).adaptor_soa];
        [~, order]               = sort(adaptor_soas);

        % find lower and upper bounds for simultaneity window in pre-test
        lc_pre = tau - c_pre;
        uc_pre = tau + c_pre;

        [pre_afirst, pre_simul, pre_vfirst] = pmf_exp(model.test_soa,...
            tau, sigma_a, sigma_v, lc_pre, uc_pre, lambda);

        nLL_ses = NaN(1, model.num_ses);
        for i = 1:model.num_ses

            adaptor = order(i);
            adaptor_soa = model.sim_adaptor_soa(i);

            pre_LL = data(adaptor).pre_nT_A1st*log(pre_afirst)'...
                + data(adaptor).pre_nT_V1st*log(pre_vfirst)'...
                + data(adaptor).pre_nT_simul*log(pre_simul)';

            if adaptor_soa > uc_pre
                uc_post = uc_pre + dc_post(adaptor);
                lc_post = lc_pre;
            elseif adaptor_soa < lc_pre
                lc_post = lc_pre - dc_post(adaptor);
                uc_post = uc_pre;
            else
                uc_post = uc_pre;
                lc_post = lc_pre;
            end

            [post_afirst, post_simul, post_vfirst] = pmf_exp(model.test_soa,...
                tau, sigma_a, sigma_v, lc_post, uc_post, lambda);

            post_LL = data(adaptor).post_nT_A1st*log(post_afirst)'...
                + data(adaptor).post_nT_V1st*log(post_vfirst)'...
                + data(adaptor).post_nT_simul*log(post_simul)';

            nLL_ses(adaptor)         = pre_LL + post_LL;

        end
        out                  = nansum(nLL_ses);

    elseif strcmp(model.mode,'predict')

        out.test_soa         = model.test_soa;
        if model.test_axis_finer == 1
            out.test_soa = linspace(out.test_soa(1), out.test_soa(end), 1e3);
        end
        out.adaptor_soa      = model.sim_adaptor_soa;

        lc_pre = tau - c_pre;
        uc_pre = tau + c_pre;

        [pre_afirst, pre_simul, pre_vfirst] = pmf_exp(out.test_soa,...
            tau, sigma_a, sigma_v, lc_pre, uc_pre, lambda);

        c_shift = NaN(1, model.num_ses);
        for adaptor = 1:9
            out.pre_pmf{adaptor}     = [pre_vfirst; pre_simul; pre_afirst];
            adaptor_soa = model.sim_adaptor_soa(adaptor);

            if adaptor_soa > uc_pre
                uc_post = uc_pre + dc_post(adaptor);
                lc_post = lc_pre;
                c_shift(adaptor) = dc_post(adaptor);
            elseif adaptor_soa < lc_pre
                lc_post = lc_pre - dc_post(adaptor);
                uc_post = uc_pre;
                c_shift(adaptor) = -dc_post(adaptor);
            else
                uc_post = uc_pre;
                lc_post = lc_pre;
                c_shift(adaptor) = 0;
            end

            [post_afirst, post_simul, post_vfirst] = pmf_exp(out.test_soa,...
                tau, sigma_a, sigma_v, lc_post, uc_post, lambda);

            out.post_pmf{adaptor}    = [post_vfirst; post_simul; post_afirst];

        end

        out.pss_shift = zeros(1, model.num_ses);
        out.c_shift = c_shift;

    end

end

end