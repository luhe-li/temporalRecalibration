function nLL = exponentialModelNLL(p, s_unique, nT_A1st, nT_V1st, numTrials)

% load variables
tau = p(1); % systematic processing delay between vision and audition.  
% τ = τ_a − τ_v is a processing advantage such that τ < 0 indicates faster
% auditory processing and τ > 0 indicates faster visual processing.
criterion = p(2); % criterion of reporting simutalneous
mu_a = p(3); % mu is the rate parameter that defines exppdf. mu = 1/lambda.
mu_v = p(4); % mu is the rate parameter that defines exppdf. mu = 1/lambda.
epsilon = p(5); % lapse rate, the probability of mistakenly making other two responses, assumed to be the same for three responses
kappa = p(6); % second lapse rate, the probability of reporting the 2nd or the 3rd response, assumed to be 0.5

% probability of each response
p_afirst = cumulativeD(s_unique, tau, -criterion, 1/mu_a, 1/mu_v);
p_vfirst = 1 - cumulativeD(s_unique, tau, criterion, 1/mu_a, 1/mu_v);
p_simul = cumulativeD(s_unique, tau, criterion, 1/mu_a, 1/mu_v)...
    - cumulativeD(s_unique, tau, -criterion, 1/mu_a, 1/mu_v);

% add lapse to each response probability
% add_lapse = @(p1, p2, p3) (1-epislon) * p1 + epislon * kappa * p2 + epislon * kappa * p3;
p_afirst_lampse = add_lapse(p_afirst, p_vfirst, p_simul, epsilon, kappa);
p_vfirst_lampse = add_lapse(p_vfirst, p_afirst, p_simul, epsilon, kappa);
p_simul_lampse = add_lapse(p_simul, p_afirst, p_vfirst, epsilon, kappa);

% compute NLL
nLL = -nT_A1st*log(p_afirst_lampse)' ...
    -nT_V1st*log(p_vfirst_lampse)'...
    -(repmat(numTrials,size(nT_A1st)) - nT_A1st - nT_V1st)...
    * log(p_simul_lampse)';

end