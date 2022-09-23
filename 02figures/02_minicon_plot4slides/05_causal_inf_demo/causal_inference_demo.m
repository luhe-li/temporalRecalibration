
x = [-700, -300:100:300, 700]./1000;
y = [-0.3, -0.8, -0.6, -0.4, -0.1, 0.2, 0.5, 0.7, 0.2];
figure; hold on
set(gca, 'LineWidth', 2, 'FontSize', 20)
set(gcf, 'Position',[10 10 500 400])
% look better
xlabel('adaptor lag')
ylabel('recalibration effect')
set(gca,'TickDir','out');
xticks(0)
yticks(0)
ylim([-1.5 1.5])
xlim([-1 1])
yline(0)
plot(x,y,'ko','MarkerSize',10,'MarkerFaceColor','k')

saveas(gca,'f1','epsc')