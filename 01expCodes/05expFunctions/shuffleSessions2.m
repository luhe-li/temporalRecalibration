% shuffleSession2
clear all; clc; close all;
rng('shuffle');
adaptors = [-700, -300:100:300, 700]/1000; %secs

lenAdaptors = length(adaptors);
randomAdaptorOrder = NaN(100, lenAdaptors);

for s = 9:20 % suppose there are at most 20 subjects
    randomAdaptorOrder(s,:) = shuffleit(adaptors);
end

save('RandomAdaptorOrder')
