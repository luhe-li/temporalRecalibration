
function out   = nll_cauInf_asym(freeParam, model, data)

if strcmp(model.mode, 'initialize')

    out.paraID   = {'\tau','\sigma_{A}','\sigma_{V}','c','\lambda','p_{common}','\alpha','\sigma_{C=1}','\sigma_{C=2}'};
    out.num_para = length(out.paraID);

    % hard bounds, the range for LB, UB, larger than soft bounds
    paraH.tau= [-100,   100]; % ms
    paraH.sigma_a  = [  10,   120]; % ms
    paraH.sigma_v  = [  10,   200]; % ms
    paraH.criterion= [   1,   350]; % criterion, s
    paraH.lambda   = [1e-4,  0.06]; % percentage
    paraH.p_common = [1e-4, 1-1e-4]; % weight
    paraH.alpha    = [1e-4,  0.02]; % percentage
    paraH.sigma_C1 = [   1,   300]; % ms
    paraH.sigma_C2 = [ 100,   1e3]; % ms

    % soft bounds, the range for PLB, PUB
    paraS.tau= [ -50,    50]; % ms
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
    out.init = getInit(out.lb, out.ub, numSections, model.num_runs);

else

    %% assign free parameters

    tau = freeParam(1);
    sigma_a = freeParam(2);
    sigma_v = freeParam(3);
    criterion = freeParam(4);
    lambda = freeParam(5);
    p_common = freeParam(6);
    alpha = freeParam(7);
    sigma_C1 = freeParam(8);
    sigma_C2 = freeParam(9);

    %% pre-compute

    checkPlot = 0;

    % prior
    fixP.x_axis     = -model.bound_full:1:model.bound_full;
    fixP.x_axis_int = -model.bound_int:1:model.bound_int;
    fixP.l_window   = find(fixP.x_axis == fixP.x_axis_int(1));
    fixP.r_window   = find(fixP.x_axis == fixP.x_axis_int(end));
    fixP.prior_C1   = normpdf(fixP.x_axis_int, 0, sigma_C1);
    fixP.prior_C2   = normpdf(fixP.x_axis_int, 0, sigma_C2);

    % likelihood that centers around 0
    idx_peak   = ceil(length(fixP.x_axis)/2);

    lf = (1/(sigma_a + sigma_v)).* exp(1/sigma_a .* (fixP.x_axis(1:idx_peak)));
    rf = (1/(sigma_a + sigma_v)).* exp(-1/sigma_v .* (fixP.x_axis(idx_peak+1:end)));
    fixP.df_likelihood     = [lf, rf] + realmin;
    %     fixP.df_likelihood_int = fixP.df_likelihood(fixP.l_window:fixP.r_window);

    % iCDF
    fixP.y_criterion = sigma_v/(sigma_a + sigma_v);
    fixP.num_sample = model.num_sample;

    %% loop for each session

    if strcmp(model.mode, 'optimize')

        %%  prediction of probability of three responses is shared across sessions
        [pre_afirst, pre_simul, pre_vfirst] = pmf_exp_CI(model.test_soa, fixP,...
            tau, sigma_a, sigma_v, criterion, lambda, p_common);

        if checkPlot
            figure; plot(model.test_soa, [pre_afirst; pre_simul; pre_vfirst],'-o');
        end

        %% simulate shift_mu of all adaptors

        % simulate using adaptor_soa in session order
        adaptor_soas = [data(1:model.num_ses).adaptor_soa];
        tau_shift = NaN(numel(adaptor_soas), model.expo_num_sim);

        for t    = 1:model.expo_num_sim

            % simulate by adaptor soas in each session, unsorted
            tau_shift(:, t)  = sim_recalibration(model.expo_num_trial, adaptor_soas, fixP,...
                tau, sigma_a, sigma_v, p_common, alpha);

        end

        if checkPlot
            [~, order] = sort(adaptor_soas);
            figure; plot(-mean(tau_shift(order,:),2))
            title('PSS shift')
        end

        nLL_ses = NaN(1, model.num_ses);
        for ses = 1:model.num_ses

            %% calculate pretest nLL

            pre_LL = data(ses).pre_nT_A1st*log(pre_afirst)'...
                + data(ses).pre_nT_V1st*log(pre_vfirst)'...
                + data(ses).pre_nT_simul*log(pre_simul)';

            %% extract the simulated mu_shift for this session

            tau_shift = tau_shift(ses, :);

            %% approximate the probability of shift_pss by Gaussian

            % find the max and min of pss_shift
            shift_min= min(tau_shift);
            shift_max= max(tau_shift);
            shift_range    = shift_max - shift_min;

            % define the lower and upper boundaries (i.e., min-(max-min), max+(max-min))
            shift_lb = shift_min - shift_range;
            shift_ub = shift_max + shift_range;
            delta_tau_shift = linspace(shift_lb, shift_ub, model.num_bin);

            % fit a Gaussian by using empirical mean and s.d.
            gauss_tau = mean(tau_shift);
            gauss_sigma    = sqrt(sum((tau_shift - gauss_tau).^2)./numel(tau_shift)); % denominator is N instead of N-1
            gauss_pdf= normpdf(delta_tau_shift, gauss_tau, gauss_sigma); % approximated gaussian pdf

            % compute R2 of using Gaussian to approxiamte pss_shift pdf
            binsize  = diff(delta_tau_shift(1:2));
            predicted_y    = gauss_pdf./sum(gauss_pdf); % normalize
            delta_shift_tau_edges = [delta_tau_shift, delta_tau_shift(end) + binsize] - binsize/2; % create edges around delta_pss_shift
            observed_y     = histcounts(tau_shift, delta_shift_tau_edges)./numel(tau_shift); % manually normalize counts to probability
            R  = corr(predicted_y(:), observed_y(:));
            R2 = R^2;

            %if R2 is big enough, we fit a Gaussian, else use ksdensity
            if R2 > model.thres_R2; pdf_delta = gauss_pdf;
            else; pdf_delta = ksdensity(tau_shift); end
            pdf_delta = pdf_delta./sum(pdf_delta);

            if checkPlot
                figure;
                plot(delta_tau_shift, pdf_delta)
            end

            %% calculate posttest nLL

            % compute the likelihood of approxiamated delta: P(resp|delta_pss_shift, M,
            % \theta)
            LL_delta = NaN(1, length(delta_tau_shift));

            for i    = 1:numel(delta_tau_shift)

                % delta_tau = tau_pre - tau_post, use tau_post to predict probability of three responses

                [post_afirst, post_simul, post_vfirst] = pmf_exp_CI(model.test_soa, fixP,...
                    tau + delta_tau_shift(i), sigma_a, sigma_v, criterion, lambda, p_common);

                if checkPlot
                    figure; plot(model.test_soa, [post_afirst; post_simul; post_vfirst] ,'-o');
                end

                LL_delta(i) = data(ses).post_nT_A1st*log(post_afirst)'...
                    + data(ses).post_nT_V1st*log(post_vfirst)'...
                    + data(ses).post_nT_simul*log(post_simul)';

            end

            % post LL is the log sum of (likelihood of approximated delta x
            % probability of approximated delta)
            % post_LL = log(sum(exp(LL_delta) .* pdf_delta));

            % re-written to avoid underflow of likelihood. Note that const can be
            % subtracted to the exponent and added later because it is NOT summed, and
            % log(exp(const)) = const
            const     = max(LL_delta + log(pdf_delta));
            post_LL   = log(sum(exp(LL_delta + log(pdf_delta) - const))) + const;

            % sum the negative likelihood of pre and post test
            nLL_ses(ses)   = - pre_LL - post_LL;

        end

        out= nansum(nLL_ses);

        if checkPlot
            [~, order] = sort(adaptor_soas);
            figure; hold on
            plot(1:9, [-post_LL_resp(order,:)])
            legend('a','v','simul')
            ylabel('nll')
            xlabel('adaptor')
            ylim([0, 350])
        end


    elseif strcmp(model.mode,'predict')

        %% pre-test TOJ

        if model.toj_axis_finer == 1
            out.test_soa = linspace(model.test_soa(1), model.test_soa(end), 1e3);
        else
            out.test_soa = model.test_soa;
        end

        if model.adaptor_axis_finer ==1
            out.num_adaptor = 30;
            out.adaptor_soa = linspace(model.sim_adaptor_soa(1), model.sim_adaptor_soa(end), out.num_adaptor );
        else
            out.adaptor_soa = model.sim_adaptor_soa;
            out.num_adaptor = numel(model.sim_adaptor_soa);
        end

        [pre_afirst, pre_simul, pre_vfirst] = pmf_exp_CI(out.test_soa, fixP,...
            tau, sigma_a, sigma_v, criterion, lambda, p_common);
        out.pre_pmf     = [pre_vfirst; pre_simul; pre_afirst];

        if checkPlot
            figure; plot(out.test_soa, out.pre_pmf)
        end
        %% simulate shift_mu of all adaptors

        % simulate in ordered adaptor_soa
        out.tau_shift  = NaN(out.num_adaptor , model.expo_num_sim);

        for t    = 1:model.expo_num_sim

            % simulate by sorted adaptor soas
            [out.tau_shift(:,t), out.post_C1(:,t), out.shat(:,t)] = sim_recalibration(model.expo_num_trial, out.adaptor_soa, fixP,...
                tau, sigma_a, sigma_v, p_common, alpha);

        end

        % pss shift = -tau shift
        out.pss_shift = -out.tau_shift;

        for jj = 1:out.num_adaptor

            %% extract the simulated mu_shift for this session

            tau_shift = out.tau_shift(jj, :);

            %% approximate the probability of shift_mu by Gaussian

            % find the max and min of mu_shift
            shift_min = min(tau_shift);
            shift_max = max(tau_shift);
            shift_range = shift_max - shift_min;

            % define the lower and upper boundaries (i.e., min-(max-min), max+(max-min))
            shift_lb = shift_min - shift_range;
            shift_ub = shift_max + shift_range;
            delta_tau_shift = linspace(shift_lb, shift_ub, model.num_bin);

            % fit a Gaussian by using empirical mean and s.d.
            gauss_tau = mean(tau_shift);
            gauss_sigma    = sqrt(sum((tau_shift - gauss_tau).^2)./numel(tau_shift)); % denominator is N instead of N-1
            gauss_pdf= normpdf(delta_tau_shift, gauss_tau, gauss_sigma); % approximated gaussian pdf

            % compute R2 of using Gaussian to approxiamte mu_shift pdf
            binsize  = diff(delta_tau_shift(1:2));
            predicted_y    = gauss_pdf./sum(gauss_pdf); % normalize
            delta_shift_mu_edges = [delta_tau_shift, delta_tau_shift(end) + binsize] - binsize/2; % create edges around delta_mu_shift
            observed_y     = histcounts(tau_shift, delta_shift_mu_edges)./numel(tau_shift); % manually normalize counts to probability
            R  = corr(predicted_y(:), observed_y(:));
            out.R2(jj)     = R^2;

            %% posttest TOJ

            [post_afirst, post_simul, post_vfirst] = pmf_exp_CI(out.test_soa, fixP,...
                tau + mean(tau_shift), sigma_a, sigma_v, criterion, lambda, p_common);

            out.post_tau(jj) = tau + mean(tau_shift);
            out.post_pmf(jj, :, :) = [post_vfirst; post_simul; post_afirst];

        end

    end

end

end
