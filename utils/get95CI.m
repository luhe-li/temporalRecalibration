function [CI_lb, CI_ub] = get95CI(v)

%This function takes a vector v, sorts the entries in an ascending order,
%and then returns two scalars that correspond the 2.5% and 97.5% perentiles

%first sort the vector in an ascending order
n_sorted = sort(v);
%compute how many entries the vector has
lenN     = length(v);
%lower bound
CI_ub    = n_sorted(ceil(lenN*0.975));
%upper bound
CI_lb    = n_sorted(floor(lenN*0.025));
end