%This function takes a vector v, sorts the entries in an ascending order,
%and then returns two scalars that correspond the 16% and 84% perentiles

function [CI_lb, CI_ub] = get68CI(v)
%first sort the vector in an ascending order
n_sorted = sort(v);
%compute how many entries the vector has
lenN     = length(v);
%upper bound 
CI_ub    = n_sorted(ceil(lenN*0.84));
%lower bound
CI_lb    = n_sorted(floor(lenN*0.16));