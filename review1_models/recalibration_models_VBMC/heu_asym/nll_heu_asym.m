function out   = nll_heu_asym(freeParam, model, data)

if strcmp(model.mode, 'initialize')

    out.paraID   = {'\tau','\sigma_{A}','\sigma_{V}','c','\lambda','\alpha'};
    out.num_para = length(out.paraID);

    % hard bounds, the range for LB, UB, larger than soft bounds
    paraH.tau      = [-200,   200]; % ms
    paraH.sigma_a  = [  10,   200]; % ms
    paraH.sigma_v  = [  10,   200]; % ms
    paraH.criterion= [   1,   350]; % criterion, s
    paraH.lambda   = [1e-4,  0.06]; % percentage
    paraH.alpha    = [1e-4,  0.02]; % percentage

    % soft bounds, the range for PLB, PUB
    paraS.tau      = [ -40,    40]; % ms
    paraS.sigma_a  = [  20,    50]; % ms
    paraS.sigma_v  = [  20,    50]; % ms
    paraS.criterion= [  30,    80]; % criterion, s
    paraS.lambda   = [0.01,  0.03]; % percentage
    paraS.alpha    = [1e-3,  2e-3]; % percentage

    % reorganize parameter bounds 
    fn = fieldnames(paraH);
    for k= 1:numel(fn)
        out.lb(:,k)  = paraH.(fn{k})(1);
        out.ub(:,k)  = paraH.(fn{k})(2);
        out.plb(:,k) = paraS.(fn{k})(1);
        out.pub(:,k) = paraS.(fn{k})(2);
    end

    % get grid initializations (bounds not included)
    numSections = model.num_runs * 2;
    out.init = getInit(out.plb, out.pub, numSections, model.num_runs);

else

    %% assign free parameters

    tau = freeParam(1);
    sigma_a = freeParam(2);
    sigma_v = freeParam(3);
    criterion = freeParam(4);
    lambda = freeParam(5);
    alpha = freeParam(6);

    checkPlot = 0;

    %% loop for each session

    if strcmp(model.mode, 'optimize')

        %%  prediction of probability of three responses is shared across sessions

        [pre_afirst, pre_simul, pre_vfirst] = pmf_exp(model.test_soa,...
            tau, sigma_a, sigma_v, -criterion, criterion, lambda);

        if checkPlot
            figure; plot(model.test_soa, [pre_afirst; pre_simul; pre_vfirst],'-o');
        end

        %% simulate shift_mu of all adaptors

        % simulate using adaptor_soa in session order
        adaptor_soas = [data(1:model.num_ses).adaptor_soa];
        tau_shift = NaN(numel(adaptor_soas), model.expo_num_sim);

        for t    = 1:model.expo_num_sim

            % simulate by adaptor soas in each session, unsorted
            tau_shift(:, t)  = sim_recal_heu(model.expo_num_trial, adaptor_soas, ...
                tau, sigma_a, sigma_v, criterion, alpha);

        end

        if checkPlot
            [~, order] = sort(adaptor_soas);
            figure; plot(-mean(tau_shift(order,:),2))
            title('PSS shift')
        end

        LL_ses = NaN(1, model.num_ses);
        for ses = 1:model.num_ses

            %% calculate pretest nLL

            pre_LL = data(ses).pre_nT_A1st*log(pre_afirst)'...
                + data(ses).pre_nT_V1st*log(pre_vfirst)'...
                + data(ses).pre_nT_simul*log(pre_simul)';

            %% extract the simulated mu_shift for this session

            i_tau_shift = tau_shift(ses, :);

            %% approximate the probability of shift_pss by Gaussian

            % find the max and min of pss_shift
            shift_min= min(i_tau_shift);
            shift_max= max(i_tau_shift);
            shift_range    = shift_max - shift_min;

            % define the lower and upper boundaries (i.e., min-(max-min), max+(max-min))
            shift_lb = shift_min - shift_range;
            shift_ub = shift_max + shift_range;
            delta_tau_shift = linspace(shift_lb, shift_ub, model.num_bin);

            % fit a Gaussian by using empirical mean and s.d.
            gauss_tau = mean(i_tau_shift);
            gauss_sigma = sqrt(sum((i_tau_shift - gauss_tau).^2)./numel(i_tau_shift)); % denominator is N instead of N-1
            gauss_pdf= normpdf(delta_tau_shift, gauss_tau, gauss_sigma); % approximated gaussian pdf

            % compute R2 of using Gaussian to approxiamte pss_shift pdf
            binsize  = diff(delta_tau_shift(1:2));
            predicted_y    = gauss_pdf./sum(gauss_pdf); % normalize
            delta_shift_tau_edges = [delta_tau_shift, delta_tau_shift(end) + binsize] - binsize/2; % create edges around delta_pss_shift
            observed_y     = histcounts(i_tau_shift, delta_shift_tau_edges)./numel(i_tau_shift); % manually normalize counts to probability
            R  = corr(predicted_y(:), observed_y(:));
            R2 = R^2;

            %if R2 is big enough, we fit a Gaussian, else use ksdensity
            if R2 > model.thres_R2; pdf_delta = gauss_pdf; 
            else; pdf_delta = ksdensity(i_tau_shift); end
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
                [post_afirst, post_simul, post_vfirst] = pmf_exp(model.test_soa,...
                    tau + delta_tau_shift(i), sigma_a, sigma_v, -criterion, criterion, lambda);

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
            LL_ses(ses)   = pre_LL + post_LL;

        end

        out= nansum(LL_ses);

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
            out.test_soa = linspace(model.test_soa(1), model.test_soa(end), 500);
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

        [pre_afirst, pre_simul, pre_vfirst] = pmf_exp(out.test_soa,...
            tau, sigma_a, sigma_v, -criterion, criterion, lambda);
        out.pre_pmf     = [pre_vfirst; pre_simul; pre_afirst];

        if checkPlot
            figure; plot(out.test_soa, out.pre_pmf)
        end
        %% simulate shift_mu of all adaptors

        % simulate in ordered adaptor_soa
        out.tau_shift  = NaN(out.num_adaptor , model.expo_num_sim);

        for t    = 1:model.expo_num_sim

            % simulate by sorted adaptor soas
            out.tau_shift(:,t) = sim_recal_heu(model.expo_num_trial, out.adaptor_soa, ...
                tau, sigma_a, sigma_v, criterion, alpha);

        end

        % pss shift = -tau shift
        out.pss_shift = -out.tau_shift;

        for jj = 1:out.num_adaptor

            %% extract the simulated mu_shift for this session

            i_tau_shift = out.tau_shift(jj, :);

            %% approximate the probability of shift_mu by Gaussian

            % find the max and min of mu_shift
            shift_min = min(i_tau_shift);
            shift_max = max(i_tau_shift);
            shift_range = shift_max - shift_min;

            % define the lower and upper boundaries (i.e., min-(max-min), max+(max-min))
            shift_lb = shift_min - shift_range;
            shift_ub = shift_max + shift_range;
            delta_tau_shift = linspace(shift_lb, shift_ub, model.num_bin);

            % fit a Gaussian by using empirical mean and s.d.
            gauss_tau = mean(i_tau_shift);
            gauss_sigma    = sqrt(sum((i_tau_shift - gauss_tau).^2)./numel(i_tau_shift)); % denominator is N instead of N-1
            gauss_pdf= normpdf(delta_tau_shift, gauss_tau, gauss_sigma); % approximated gaussian pdf

            % compute R2 of using Gaussian to approxiamte mu_shift pdf
            binsize  = diff(delta_tau_shift(1:2));
            predicted_y    = gauss_pdf./sum(gauss_pdf); % normalize
            delta_shift_mu_edges = [delta_tau_shift, delta_tau_shift(end) + binsize] - binsize/2; % create edges around delta_mu_shift
            observed_y     = histcounts(i_tau_shift, delta_shift_mu_edges)./numel(i_tau_shift); % manually normalize counts to probability
            R  = corr(predicted_y(:), observed_y(:));
            out.R2(jj)     = R^2;

            %% posttest TOJ

            [post_afirst, post_simul, post_vfirst] = pmf_exp(out.test_soa,...
                tau + mean(i_tau_shift), sigma_a, sigma_v, -criterion, criterion, lambda);

            out.post_tau(jj) = tau + mean(i_tau_shift);
            out.post_pmf(jj, :, :) = [post_vfirst; post_simul; post_afirst];

        end

    end

end

end
