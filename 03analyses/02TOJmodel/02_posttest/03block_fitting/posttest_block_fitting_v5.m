% This script fits post- tests TOJ data by block, assuming the same pre-mu
% but 5 block post-mus.

clear all; close all; clc;
for subjID = 4:4
    figure; hold on
    set(gca, 'FontSize', 20, 'LineWidth', 1.5); box off
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.7, 1]);
    for sess = 1:9

        %% organize data

        %%%%% pre-test
        load(['pretest_sub' num2str(subjID) '_session' num2str(sess) '.mat'])
        pre_s_unique = ExpInfo.SOA; % unique SOA levels, in ms
        pre_numTrials = ExpInfo.nTrials; % num of trials per SOA
        % inititate
        pre_r_org = NaN(length(pre_s_unique), pre_numTrials);
        pre_respCount = NaN(3, length(pre_s_unique));
        for i = 1:length(pre_s_unique)
            iSOA = pre_s_unique(i);
            iResp = Response.order(ExpInfo.trialSOA == iSOA);
            pre_r_org(i,:) = iResp; % this matrix has a size of length(s_unique) x numTrials
            for j = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
                pre_respCount(j,i) = sum(iResp == j);
            end
        end
        pre_pResp = pre_respCount/pre_numTrials;

        %%%%% post-test
        % load data and define key parameters
        load(['posttest_sub' num2str(subjID) '_session' num2str(sess) '.mat'])
        post_s_unique = ExpInfo.SOA; % unique SOA levels, in ms
        post_numTrials = ExpInfo.nTrials; % num of trials per SOA

        % arbituarily bin post-test trials into 5 blocks
        nblock = 5; % how many blocks to fit
        tperblock = ExpInfo.nTTrials/nblock; % trials per block
        b_numTrials = post_numTrials/nblock; % trials per soa per block

        % initiate
        [post_r_org, post_respCount, post_pResp] = deal(cell(1,nblock));

        % organize post-test data by block
        for b = 1:5
            bstart = (b-1)*tperblock+1;
            bend = b*tperblock;
            bResp = Response.order(bstart:bend);
            bTrialSOA = ExpInfo.trialSOA(bstart:bend);

            % organize response
            for i = 1:length(post_s_unique)
                iSOA = post_s_unique(i);
                iResp = bResp(bTrialSOA == iSOA);
                post_r_org{b}(i,:) = iResp; % this matrix has a size of length(s_unique) x numTrials
                for j = unique(Response.order) % 1 = V first, 2 = simultaneous, 3 = A first
                    post_respCount{b}(j,i) = sum(iResp == j);
                end
            end
            % convert to percentage
            post_pResp{b} = post_respCount{b}/b_numTrials;
        end

        %% fitting
        % define PMF
        P_Afirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
        P_Vfirst = @(SOA, mu, sig, lambda, c) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
        P_simultaneous = @(SOA, mu, sig, lambda, c) ...
            1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
            - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

        % preprare behavioral responses
        pre_s_unique_s = pre_s_unique*1e3; post_s_unique_s = post_s_unique*1e3; % s->ms

        % cost function
        % response(3,:): the number of A-first responses for each SOA
        % response(1,:): the number of V-first responses for each SOA
        nLL = @(p) -pre_respCount(3,:)*log(P_Afirst(pre_s_unique_s, p(1), p(7), p(8), p(9)))' ...
            -pre_respCount(1,:)*log(P_Vfirst(pre_s_unique_s, p(1), p(7), p(8), p(9)))'...
            -(repmat(pre_numTrials,size(pre_respCount(3,:))) - pre_respCount(3,:) - pre_respCount(1,:))...
            * log(P_simultaneous(pre_s_unique_s, p(1), p(7), p(8), p(9)))'...% pre-test kokomade
            -post_respCount{1}(3,:)*log(P_Afirst(post_s_unique_s, p(2), p(7), p(8), p(9)))' ...
            -post_respCount{1}(1,:)*log(P_Vfirst(post_s_unique_s, p(2), p(7), p(8), p(9)))'...
            -(repmat(post_numTrials,size(post_respCount{1}(3,:))) - post_respCount{1}(3,:) - post_respCount{1}(1,:))...
            * log(P_simultaneous(post_s_unique_s, p(2), p(7), p(8), p(9)))'... % block1 post-test
            -post_respCount{2}(3,:)*log(P_Afirst(post_s_unique_s, p(3), p(7), p(8), p(9)))' ...
            -post_respCount{2}(1,:)*log(P_Vfirst(post_s_unique_s, p(3), p(7), p(8), p(9)))'...
            -(repmat(post_numTrials,size(post_respCount{1}(3,:))) - post_respCount{2}(3,:) - post_respCount{2}(1,:))...
            * log(P_simultaneous(post_s_unique_s, p(3), p(7), p(8), p(9)))'... % block2 post-test
            -post_respCount{3}(3,:)*log(P_Afirst(post_s_unique_s, p(4), p(7), p(8), p(9)))' ...
            -post_respCount{3}(1,:)*log(P_Vfirst(post_s_unique_s, p(4), p(7), p(8), p(9)))'...
            -(repmat(post_numTrials,size(post_respCount{1}(3,:))) - post_respCount{3}(3,:) - post_respCount{3}(1,:))...
            * log(P_simultaneous(post_s_unique_s, p(4), p(7), p(8), p(9)))'... % block3 post-test
            -post_respCount{4}(3,:)*log(P_Afirst(post_s_unique_s, p(5), p(7), p(8), p(9)))' ...
            -post_respCount{4}(1,:)*log(P_Vfirst(post_s_unique_s, p(5), p(7), p(8), p(9)))'...
            -(repmat(post_numTrials,size(post_respCount{1}(3,:))) - post_respCount{4}(3,:) - post_respCount{4}(1,:))...
            * log(P_simultaneous(post_s_unique_s, p(5), p(7), p(8), p(9)))'... % block4 post-test
            -post_respCount{5}(3,:)*log(P_Afirst(post_s_unique_s, p(6), p(7), p(8), p(9)))' ...
            -post_respCount{5}(1,:)*log(P_Vfirst(post_s_unique_s, p(6), p(7), p(8), p(9)))'...
            -(repmat(post_numTrials,size(post_respCount{1}(3,:))) - post_respCount{5}(3,:) - post_respCount{5}(1,:))...
            * log(P_simultaneous(post_s_unique_s, p(6), p(7), p(8), p(9)))'; % block5 post-test

        % set lower and upper bounds
        % pre_mu, b1_mu, b2_mu, b3_mu, b4_mu, b5_mu, sigma, lambda, criterion
        lb      = [-150, -150, -150, -150, -150, -150, 10, 1e-2, 50];
        ub      = [150, 150, 150, 150, 150, 150, 350, 0.06, 350];

        % choose random initial values for 1e3 times
        for i = 1:1
            init    = rand(1,length(lb)).*(ub-lb) + lb;
            %You can also define how many times you want MATLAB to search
            options = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');

            %fmincon returns best-fitting parameters that minimize the cost function as
            %well as the corresponding value for the cost function (in this case, the
            %negative log likelihood)
            [estP(i,:), min_NLL(i)] = fmincon(nLL, init,[],[],[],[],lb,ub,[],options);
        end

        % use the best-fitting parameters with the smallest NLL among 1e3 fittings
        [value idx] = min(min_NLL);
        bestP = estP(idx,:);
        fprintf('pre mu = %4.2f post mu 1= %4.2f post mu 2= %4.2f post mu 3= %4.2f post mu 4= %4.2f post mu 5= %4.2f sig = %4.2f lambda = %4.2f criterion = %4.2f\n', bestP)

        %% plot
        subplot(3,3,sess)
        h(sess) = plot(bestP(1:6),'.-');
        xticks([1:6])
        xlim([1 6])
        xticklabels({'pre mu','post mu 1', 'post mu 2', 'post mu 3', 'post mu 4', 'post mu 5'})
        title(['sub' num2str(subjID) ' adaptor SOA = ' num2str(ExpInfo.adaptor*1000) ' ms'])

    end
    linkaxes(h,'y')
end