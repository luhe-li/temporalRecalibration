# Fit the Bayesian causal inference model

## set PLB and PUB
 PLB and PUB set a plausible range in which you would expect to find almost all solutions.
 The plausible box defined by PLB and PUB naturally represents a good region where to randomly draw starting points for the optimization.

paraS.mu1 = [-0.25, 0.25]; % s
paraS.sigma1 = [0.01, 0.3]; %s
paraS.c1 = [0.01, 0.3]; % s
paraS.lambda = [1e-4, 0.06]; % percentage
paraS.sigma2 = [0.01, 0.3]; %s
paraS.c2 = [0.01, 0.3]; % s
paraS.p_common = [1e-2, 0.99]; % weight
paraS.sigma_soa  = [0.01, 0.35]; % s
paraS.sigma_c1 = [1e-4, 0.05]; % s
paraS.sigma_c2  = [0.5, 2]; % s
paraS.alpha = [1e-3, 0.03]; % percentage

## set LB and UB
LB and UB are the hard bounds of the optimization. They are chosen to set a larger range than PLB and PUB.

paraH.mu1 = [-0.3, 0.3]; % s
paraH.sigma1 = [0.01, 0.35]; %s
paraH.c1 = [0.01, 0.35]; % s
paraH.lambda = [1e-4, 0.06]; % percentage
paraH.sigma2 = [0.01, 0.35]; %s
paraH.c2 = [0.01, 0.35]; % s
paraH.p_common = [1e-2, 0.99]; % weight
paraH.sigma_soa  = [0.01, 0.35]; % s
paraH.sigma_c1 = [1e-4, 0.05]; % s
paraH.sigma_c2  = [0.2, 3]; % s
paraH.alpha = [1e-3, 0.05]; % percentage

## default/suggested values of free parameters (for debugging cal_NLL_CI)
mu1 = 0.03; % s
sigma1 = 0.07; %s
c1 = 0.09;
lambda = 0.03;
sigma2 = 0.09;
c2 = 0.12;
p_common = 0.3;
sigma_soa  = 0.16;
sigma_c1 = 0.01;
sigma_c2  = 1;
alpha  = 0.01;

nLL = cal_nLL_CI(mu1, sigma1, c1, lambda,... % pretest free para
    sigma2, c2, ... % posttest free para
    p_common, sigma_soa, sigma_c1, sigma_c2, alpha, ... % recal simul para
    model, data); % fixed para, data

## markdown configuration
https://sublimetext-markdown.github.io/MarkdownEditing/config/#configuration
