function [AIC, NLL, estP] = fit_6M_v5(sim_r_org, ExpInfo)
    
%% model parameters outside the loop
% define PMF
P_Afirst                      = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst                      = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous                = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

numM                          = 6; %number of models
n_init                        = 10;
numP                          = [8,7,7,6,6,5]; %number of free parameters for M1, M2, M3, M4, M5 respectively

%define upper and lower bounds
lb                            = {[-250, -250, 10, 10, 10, 10, 1e-4, 1e-4], ...%M1
    [-250, -250, 10, 10, 10, 1e-4, 1e-4] ...%M2
    [-250, -250, 10, 10, 10, 1e-4, 1e-4], ...%M3
    [-250, -250, 10, 10, 1e-4, 1e-4], ...%M4
    [-250, 10, 10, 10, 1e-4, 1e-4], ...%M5
    [-250, 10, 10, 1e-4, 1e-4]};%M6
ub                            = {[250, 250, 350, 350, 350, 350, 0.06, 0.06], ...%M1
    [250, 250, 350, 350, 350, 0.06, 0.06], ...%M2
    [250, 250, 350, 350, 350, 0.06, 0.06], ...%M3
    [250, 250, 350, 350, 0.06, 0.06], ...%M4
    [250, 350, 350, 350, 0.06, 0.06], ...%M5
    [250, 350, 350, 0.06, 0.06]};%M6

init_fun                      = @(a,b) rand(1,length(a)).*(b-a) + a;
options                       = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');
estP                          = cell(1,numM); % 1 session x n_models
[NLL, AIC]                = deal(NaN(1, numM)); % 1 session x n_models

%% organize data outside the loop

% extract experimental parameters
pre_ms_unique                  = ExpInfo.SOA; % unique SOA levels, in s
post_ms_unique                 = ExpInfo.SOA; % unique SOA levels, in ms


% initiation
[respCount, pResp]            = deal(cell(1,2));

for sess                     = 1:2 % loop through pre and post
    % calculate response count and p response as fitting
    for i                         = 1:ExpInfo.lenS % for each SOA
        iResp                         = sim_r_org{sess}(i,:); % iResp is the responses in each SOA level
        for j                         = unique(iResp) % for each response type: % 1 = V first, 2 = simultaneous, 3 = A first
            respCount{sess}(j,i)             = sum(iResp == j); %
        end
    end
    pResp{sess}                      = respCount{sess}./ExpInfo.nTrials;
end

% plot simulated response to double-check
% figure;hold on; plot(pre_ms_unique, pResp{1},'o'); plot(post_ms_unique, pResp{2},'x');

pre_nT_V1st                   = respCount{1}(1,:);
pre_nT_simul                  = respCount{1}(2,:);
pre_nT_A1st                   = respCount{1}(3,:);
post_nT_V1st                  = respCount{2}(1,:);
post_nT_simul                 = respCount{2}(2,:);
post_nT_A1st                  = respCount{2}(3,:);

%% specify models
    M1                            = @(p) -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(3), p(5), p(7)))'...
        -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(3), p(5), p(7)))'...
        -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(3), p(5), p(7)))'...
        -post_nT_A1st*log(P_Afirst(post_ms_unique, p(2), p(4), p(6), p(8)))' ...
        -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(2), p(4), p(6), p(8)))'...
        -post_nT_simul*log(P_simultaneous(post_ms_unique, p(2), p(4), p(6), p(8)))';

    M2                            = @(p) -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(3), p(4), p(6)))' ...
        -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(3), p(4), p(6)))'...
        -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(3), p(4), p(6)))'...
        -post_nT_A1st*log(P_Afirst(post_ms_unique, p(2), p(3), p(5), p(7)))' ...
        -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(2), p(3), p(5), p(7)))'...
        -post_nT_simul*log(P_simultaneous(post_ms_unique, p(2), p(3), p(5), p(7)))';

    M3                            = @(p) -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(3), p(5), p(6)))' ...
        -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(3), p(5), p(6)))'...
        -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(3), p(5), p(6)))'...
        -post_nT_A1st*log(P_Afirst(post_ms_unique, p(2), p(4), p(5), p(7)))' ...
        -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(2), p(4), p(5), p(7)))'...
        -post_nT_simul*log(P_simultaneous(post_ms_unique, p(2), p(4), p(5), p(7)))';...

    M4                            = @(p) -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(3), p(4), p(5)))' ...
        -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(3), p(4), p(5)))'...
        -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(3), p(4), p(5)))'...
        -post_nT_A1st*log(P_Afirst(post_ms_unique, p(2), p(3), p(4), p(6)))' ...
        -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(2), p(3), p(4), p(6)))'...
        -post_nT_simul*log(P_simultaneous(post_ms_unique, p(2), p(3), p(4), p(6)))';

    M5                            = @(p) -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(2), p(3), p(5)))' ...
        -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(2), p(3), p(5)))'...
        -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(2), p(3), p(5)))'...
        -post_nT_A1st*log(P_Afirst(post_ms_unique, p(1), p(2), p(4), p(6)))' ...
        -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(1), p(2), p(4), p(6)))'...
        -post_nT_simul*log(P_simultaneous(post_ms_unique, p(1), p(2), p(4), p(6)))';

    M6                            = @(p) -pre_nT_A1st*log(P_Afirst(pre_ms_unique, p(1), p(2), p(3), p(4)))' ...
        -pre_nT_V1st*log(P_Vfirst(pre_ms_unique, p(1), p(2), p(3), p(4)))'...
        -pre_nT_simul*log(P_simultaneous(pre_ms_unique, p(1), p(2), p(3), p(4)))'...
        -post_nT_A1st*log(P_Afirst(post_ms_unique, p(1), p(2), p(3), p(5)))' ...
        -post_nT_V1st*log(P_Vfirst(post_ms_unique, p(1), p(2), p(3), p(5)))'...
        -post_nT_simul*log(P_simultaneous(post_ms_unique, p(1), p(2), p(3), p(5)))';

    nLLs                          = {M1; M2; M3; M4; M5; M6};

    %loop through models
    for m                         = 1:numM
        for ii                        = 1:n_init
            %the initial point for matlab to start searching
            init                          = init_fun(lb{m}, ub{m});
            %use fmincon.m to fit
            [ii_estP{ii}, ii_min_NLL(ii)] = fmincon(nLLs{m}, init,[],[],[],[],...
                lb{m}, ub{m},[],options);
        end
        % % use the best-fitting parameters with the smallest NLL among ninit fittings
        [NLL(m) idx]        = min(ii_min_NLL);
        estP{m}                 = ii_estP{idx};
        %compute the AIC/BIC
        AIC(m)                  = 2*NLL(m) + 2*numP(m);
        % BIC(m) = 2*min_NLL(m) + numP(m)*log(pre_numTrials);
    end

end