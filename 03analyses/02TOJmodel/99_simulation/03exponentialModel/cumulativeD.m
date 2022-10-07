function p = cumulativeD(SOA, tau, criterion, lambda_a, lambda_v)
    for i = 1:length(SOA)
        soa = SOA(i); % evaluate function at each level of SOA
        if criterion <= soa + tau
            p(i) =   lambda_a / (lambda_a + lambda_v) * exp(lambda_v * (criterion - soa - tau));
        else
            p(i) = 1 - lambda_v / (lambda_a + lambda_v) * exp( - lambda_a * (criterion - soa - tau));
        end
    end
end