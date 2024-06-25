
function out         = nll_gauss_shiftMu(freeParam, model, data)

if strcmp(model.mode, 'initialize')

    out.paraID     = {'\mu_{pre}','\sigma','c','\mu_{post,1}','\mu_{post,2}','\mu_{post,3}','\mu_{post,4}','\mu_{post,5}','\mu_{post,6}','\mu_{post,7}','\mu_{post,8}','\mu_{post,9}','\lambda'};
    out.num_para   = length(out.paraID);

    % hard bounds, the range for LB, UB, larger than soft bounds
    paraH.mu         = [-300,   300]; % ms
    paraH.sigma      = [0.01,   300]; % ms
    paraH.c          = [0.01,   350]; % ms
    paraH.mu_post1   = [-300,   300]; % ms
    paraH.mu_post2   = [-300,   300]; % ms
    paraH.mu_post3   = [-300,   300]; % ms
    paraH.mu_post4   = [-300,   300]; % ms
    paraH.mu_post5   = [-300,   300]; % ms
    paraH.mu_post6   = [-300,   300]; % ms
    paraH.mu_post7   = [-300,   300]; % ms
    paraH.mu_post8   = [-300,   300]; % ms
    paraH.mu_post9   = [-300,   300]; % ms
    paraH.lambda     = [1e-4,  0.06]; % percentage

    % soft bounds, the range for PLB, PUB
    paraS.mu         = [-100,   100]; % ms
    paraS.sigma      = [  10,   100]; % ms
    paraS.c          = [  10,   100]; % ms
    paraS.mu_post1   = [-100,   100]; % ms
    paraS.mu_post2   = [-100,   100]; % ms
    paraS.mu_post3   = [-100,   100]; % ms
    paraS.mu_post4   = [-100,   100]; % ms
    paraS.mu_post5   = [-100,   100]; % ms
    paraS.mu_post6   = [-100,   100]; % ms
    paraS.mu_post7   = [-100,   100]; % ms
    paraS.mu_post8   = [-100,   100]; % ms
    paraS.mu_post9   = [-100,   100]; % ms
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
    out.init = getInit(out.lb, out.ub, numSections, model.num_runs);

else

    %% assign free parameters

    mu = freeParam(1);
    sigma = freeParam(2);
    c = freeParam(3);
    mu_post1 = freeParam(4);
    mu_post2 = freeParam(5);
    mu_post3 = freeParam(6);
    mu_post4 = freeParam(7);
    mu_post5 = freeParam(8);
    mu_post6 = freeParam(9);
    mu_post7 = freeParam(10);
    mu_post8 = freeParam(11);
    mu_post9 = freeParam(12);
    lambda = freeParam(13);

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
end