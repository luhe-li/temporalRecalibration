% generate randomized session order for exposure phase
clear all; clc; close all;
rng('shuffle');
adaptors = [-300:100:300]/1000; %secs

% because it is impossible to cover counterbalanced orders across limited
% number of subjects, we randomize the session order and check if there is
% the order effect on the results.

for s = 1:100 % suppose there are at most 20 subjects
    randomAdaptorOrder(s,:) = Shuffle(adaptors);
end

save('RandomAdaptorOrder')