% shift on mu, simply add it to the previous mu
x = [-10:0.1:10];
mu = 0;
sig = 1;
y = normpdf(x, mu, sig);
delta = 3;
y2 = normpdf(x, mu+delta, 1);
figure; hold on
plot(x,[y;y2])

% shift on x should add a minus sign
x = [-10:0.1:10];
mu = 0;
sig = 1;
y = normpdf(x, mu, sig);
delta = 3;
y2 = normpdf(x-delta, mu, 1);
figure; hold on
plot(x,[y;y2])

nsample = 100;
z = randn(1,nsample)*sig + mu;
z2 = randn(1,nsample)*sig + mu + delta;
figure; hold on
histogram(z)
histogram(z2)