

import numpy as np
from scipy.stats import norm

# the probability of each response is cdf of measurement distribution
# evaluated at criterion
p_afirst_lapse = lambda_ / 3 + (1 - lambda_) * norm.cdf(lc, test_soa - mu, sigma) + np.finfo(float).tiny
p_vfirst_lapse = lambda_ / 3 + (1 - lambda_) * (1 - norm.cdf(uc, test_soa - mu, sigma)) + np.finfo(float).tiny
p_simul_lapse = np.ones_like(p_afirst_lapse) - p_afirst_lapse - p_vfirst_lapse

# check PMF if needed# Plot the results
plt.figure()
plt.plot(test_soa, p_afirst_lapse, label='A First')
plt.plot(test_soa, p_simul_lapse, label='Simultaneous')
plt.plot(test_soa, p_vfirst_lapse, label='V First')
plt.xlabel('Test SOA')
plt.ylabel('Probability')
plt.legend()
plt.title('Response Probabilities vs Test SOA')
plt.show()

