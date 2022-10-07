%% Plot the best-fitting parameters given each bootstrapped dataset
simTrial = simTrial_all(1);
% initiate figure for each simulated numTrial
figure; hold on
%parameter names
valName = {'\mu_{pre}','\mu_{post}', '\sigma', '\lambda', 'criterion'};
numP = 5;
cMAP = [200, 40, 40; 255, 128, 0; 13, 183, 200]./255;
figure; set(gcf,'position',[0,0,400,800])
fignm2 = ['btst parameters ' num2str(simTrial) ' trials per SOA_' num2str(nData)];
sgtitle(fignm2)
hold on;
for i = 1:numP
    subplot(numP,1,i)
    histogram(estP_btst{t}(:,i),'FaceColor', cMAP(2,:), ...
        'FaceAlpha', 0.5, 'EdgeColor',cMAP(2,:)); hold on
    plot([realP(i), realP(i)], [0, numBtst*0.3],'r--', 'lineWidth',3); hold on
    plot([CI_btst_lb{t}(:,i), CI_btst_ub{t}(:,i)],[5, 5],'-k','LineWidth', 2); hold on
    xlim([lb(i), ub(i)]); ylim([0, numBtst*0.3]);
    xticks(sort(unique([linspace(lb(i), ub(i),5), realP(i)])));
    title(valName{i}); box off
    set(gca,'FontSize',15);
end
%saveas(gca,fignm2,'png')

%% plot mu_change against each dataset (x-axis), each numTrial (figure)
figure(f(t)); hold on
mu_diff = estP_btst{t}(:,2) - estP_btst{t}(:,1);
mu_diffS{t, nData} = mu_diff;
[counts{nData} edges{nData}] = histcounts(mu_diff,10);
imagesc([nData, nData],[min(edges{nData}), max(edges{nData})],counts{nData}')

% look better
set(gca,'FontSize',15);
ylabel(['\Delta_{\mu} bootstrapped  for ' num2str(numBtst) ' trials'])
xlabel('data sets')
xlim([0, nDatasets+1])
ylim([-100 150])
l1 = yline(0,'--','LineWidth',2);
cb = colorbar;
cb.Label.String = 'number of bootstrapped trials';

%% plot mu_change CI against each dataset (x-axis), each numTrial (figure)
% compute mean
m = mean(mu_diff);
% compute 95%CI
n_sorted = sort(mu_diff);
%compute how many entries the vector has
lenN     = length(mu_diff);
%lower bound
CI_ub    = n_sorted(ceil(lenN*0.975));
%upper bound
CI_lb    = n_sorted(floor(lenN*0.025));
% plot
l2 = plot([nData, nData],[CI_lb, CI_ub],'k-','LineWidth',1.5);
legend(l2, '95 CI','Location','northwest')
%         legend([l1, l2],{'ground truth \Delta_{\mu}','95 CI'},'Location','northwest')

% compute power
if CI_lb >= 0
    reject = reject + 1; 
end
power(t) = reject/nData;
%         save('mu_diffS', 'mu_diffS','power')

% title
sumfig = ['Simulate' num2str(simTrial) ' trials per SOA'];
title(['Simulate' num2str(simTrial) ' trials per SOA, power = ' num2str(power(t))])
saveas(gca,sumfig,'png')