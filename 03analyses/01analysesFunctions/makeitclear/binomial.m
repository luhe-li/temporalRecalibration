c = 1;
figure; hold on
sgtitle('The probability of getting exactly k successes in n independent Bernoulli trials')
for p = 0.1:0.1:1
n = 25;
x = 1:n;
y = binopdf(x,n,p);
subplot(2,5,c)
plot(x, y)
ylim([0, 1])
title(['p = ' num2str(p) ])
c = c + 1;
end

c = 1;
figure; hold on
sgtitle('The cumulative probability of getting exactly k successes in n independent Bernoulli trials')
for p = 0.1:0.1:1
n = 25;
x = 1:n;
y = binocdf(x,n,p);
subplot(2,5,c)
plot(x, y)
ylim([0, 1])
title(['p = ' num2str(p)])
c = c + 1;
end