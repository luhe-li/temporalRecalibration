function sim_r_org = sampleMatrix(p, nT)

% This function sample nT trials from a probability distribution for
% multiple SOA conditions at the same time. It returns a matrix of
% simulated raw responses instead of a vector.

% inputs
% p         : each column is a prob distribution, and every column is a distribution
%           for each SOA; p should be a matrix = len of conditions x len of soa
% nT        : the number of trials to simulate

% outputs
% sim_r_org : a matrix of lenSOA x nT, where 1 = response of reporting
%           V first, 2 = response of reporting simultaneous, 3 = response of reporting
%           A first. The order is the same as experimental setting and
%           model fitting

[lenCond, lenSOA] = size(p);
m_rand = rand(lenSOA, nT); % 
cum_p = [zeros(1,lenSOA); cumsum(p)]; % added zero as the lowest boundary to create the range
sim_r_org = NaN(size(m_rand));
for i = 2:(lenCond+1)
    % when sample value falls into the range, sample is assigned as i-1 (because we added zero as the first index)
    sim_r_org(repmat(cum_p(i-1,:)', [1, nT]) <= m_rand & ...
        m_rand < repmat(cum_p(i,:)', [1, nT])) = i-1;
end