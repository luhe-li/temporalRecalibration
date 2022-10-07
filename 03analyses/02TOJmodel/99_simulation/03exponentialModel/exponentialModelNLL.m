function nLL = exponentialModelNLL(p, s_unique, nT_A1st, nT_V1st, numTrials)

% load variables
global s_unique nT_A1st nT_V1st numTrials
tau = p(1);
criterion = p(2);
lambda_a = p(3);
lambda_v = p(4);
epsilon = p(5);
kappa = p(6);

% probability of each response
p_afirst = cumulativeD(s_unique, tau, -criterion, lambda_a, lambda_v);
p_vfirst = 1 - cumulativeD(s_unique, tau, criterion, lambda_a, lambda_v);
p_simul = cumulativeD(s_unique, tau, criterion, lambda_a, lambda_v)...
    - cumulativeD(s_unique, tau, -criterion, lambda_a, lambda_v);

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