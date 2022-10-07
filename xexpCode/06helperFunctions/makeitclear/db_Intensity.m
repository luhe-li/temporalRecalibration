% making sense of decibel and intensity

% k = I/I0; I0 is the baseline intensity
% plot db against k
k = 1:10;
db = 10*log10(k);
figure; hold on
plot(k, db)
xlim([0 10])
xlabel('Intensity')
ylabel('db')

% plot k against db
figure; hold on
db2 = 1:10;
k2 = 10.^(db2/10);
xlim([0 10])
xlabel('db')
ylabel('Intensity')
plot(db2, k2)

%% calculate intensity from evenly spaced db
dbs = linspace(0, 10*log10(2),7);
ints = 10.^(dbs/10);
I0 = 0.5;
ints = ints * I0;

% check if ratio is the same
ints(2:end)./ ints(1:end-1)

% This is the same as:
ints2 = logspace(log10(0.5), log10(1), 7)