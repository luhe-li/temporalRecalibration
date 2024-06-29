function out         = nll_exp_shiftMu(freeParam, model, data)

if strcmp(model.mode, 'initialize')

    out.paraID     = {'\tau_{pre}','\tau_{post,1}', '\tau_{post,2}', '\tau_{post,3}', '\tau_{post,4}', '\tau_{post,5}', '\tau_{post,6}', '\tau_{post,7}', '\tau_{post,8}', '\tau__{post,9}','\sigma_{A}','\sigma_{V}','c','\lambda'};
    out.num_para   = length(out.paraID);

    % hard bounds, the range for LB, UB, larger than soft bounds
    paraH.tau_pre               = [-200,   200]; % ms
    paraH.tau1                  = [-200,   200]; % ms
    paraH.tau2                  = [-200,   200]; % ms
    paraH.tau3                  = [-200,   200]; % ms
    paraH.tau4                  = [-200,   200]; % ms
    paraH.tau5                  = [-200,   200]; % ms
    paraH.tau6                  = [-200,   200]; % ms
    paraH.tau7                  = [-200,   200]; % ms
    paraH.tau8                  = [-200,   200]; % ms
    paraH.tau9                  = [-200,   200]; % ms
    paraH.sigma_a               = [  10,   200]; % ms
    paraH.sigma_v               = [  10,   200]; % ms
    paraH.criterion             = [   1,   350]; % criterion, s
    paraH.lambda                = [1e-4,  0.06]; % percentage

    % soft bounds, the range for PLB, PUB
    paraS.tau_pre               = [ -40,    40]; % ms
    paraS.tau1                  = [-100,   100]; % ms
    paraS.tau2                  = [-100,   100]; % ms
    paraS.tau3                  = [-100,   100]; % ms
    paraS.tau4                  = [-100,   100]; % ms
    paraS.tau5                  = [-100,   100]; % ms
    paraS.tau6                  = [-100,   100]; % ms
    paraS.tau7                  = [-100,   100]; % ms
    paraS.tau8                  = [-100,   100]; % ms
    paraS.tau9                  = [-100,   100]; % ms
    paraS.sigma_a               = [  20,    50]; % ms
    paraS.sigma_v               = [  20,    50]; % ms
    paraS.criterion             = [  30,    80]; % criterion, s
    paraS.lambda                = [0.01,  0.03]; % percentage

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

    tau_pre = freeParam(1);
    tau1 = freeParam(2);
    tau2 = freeParam(3);
    tau3 = freeParam(4);
    tau4 = freeParam(5);
    tau5 = freeParam(6);
    tau6 = freeParam(7);
    tau7 = freeParam(8);
    tau8 = freeParam(9);
    tau9 = freeParam(10);
    sigma_a = freeParam(11);
    sigma_v = freeParam(12);
    criterion = freeParam(13);
    lambda = freeParam(14);

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
            nLL_ses(adaptor)         = pre_LL + post_LL;

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
        out.pss_shift = - (tau_post - tau_pre);
        out.c_shift = zeros(1, model.num_ses);

    end

end

end