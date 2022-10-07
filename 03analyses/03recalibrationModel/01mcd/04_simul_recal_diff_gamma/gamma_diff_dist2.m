% define difference of gamma
clear all; close all; clc

alpha1s = [0.5, 1:3]; % test 4 alpha1 listed in the sample figure
alpha2 = 1;
beta1 = 1;
beta2 = 1;
figure; hold on
syms x 
% lo_DGG = c * exp(-beta1 * z) * sum(xs.^(alpha2 - 1) .* (xs + z).^(alpha1 - 1) .* exp(-(beta1 + beta2) .* xs));

for ii = 1:length(alpha1s)
    alpha1 = alpha1s(ii);
    c = (beta1^alpha1) * (beta2^alpha2)/(gamma(alpha1) * gamma(alpha2));

    % for z > 0
    up_zs = 0.001:0.001:20;
    for i = 1:length(up_zs)
        z = up_zs(i);
%         xs = z:0.001:20; % integrate from z to 'infinity', chosen as 20
        up_DGG = c * exp(beta2 * z) * x.^(alpha1 - 1) .* (x - z).^(alpha2 - 1) .* exp(-(beta1 + beta2) .* x);
        up(i) = int(up_DGG, x, [z 20]);
    end

%     % for z < 0
%     lo_zs = -20:0.001:-0.001;
%     for i = 1:length(lo_zs)
%         z = lo_zs(i);
%         xs = -z:0.001:20; % sum up x from -z to infinity
%         lo_DGG(i) = c * exp(-beta1 * z) * sum(xs.^(alpha2 - 1) .* (xs + z).^(alpha1 - 1) .* exp(-(beta1 + beta2) .* xs));
%     end
% 
%     plot(up_zs, up_DGG)
%     plot(lo_zs, lo_DGG)
plot(up_zs, up)
end
