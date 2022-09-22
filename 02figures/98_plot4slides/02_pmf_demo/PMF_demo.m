
cmp =  [229, 158, 168; 203, 227, 172; 171,223,235]./255;

% set parameters
t = -3;
sig   = 1;
x     = -6:0.01:6;
mu = 0;
c = 1;

% obtain lines
gauss = normpdf(x, mu, sig);
P_Afirst = normcdf(-c, SOA - mu, sig);
P_Vfirst = 1 - normcdf(c, SOA-mu, sig);
P_simul = 1 - P_Afirst - P_Vfirst;

figure; hold on

subplot 121
axis off; hold on
plot(x, gauss, 'Color','k','lineWidth',2);
xline(mu - c,'lineWidth',2)
xline(mu + c,'lineWidth',2)
% idx_start = find(abs(x - x_slc(i)) < 1e-3,1); xx = x(1:idx_start);
% patch([xx, fliplr(xx)], [Gauss(1:idx_start), ...
%     zeros(1,idx_start)], [13, 183, 200]./255, 'FaceAlpha',0.2); hold on;

subplot 122









%% ref code
x_slc = -4:0.2:4;
for i = 1:length(x_slc)
    subplot(1,2,1)
    plot(x, gauss, 'Color','k','lineWidth',2); hold on;
    idx_start = find(abs(x - x_slc(i)) < 1e-3,1); xx = x(1:idx_start);
    patch([xx, fliplr(xx)], [gauss(1:idx_start), ...
        zeros(1,idx_start)], [13, 183, 200]./255, 'FaceAlpha',0.2); hold on;
    plot([x_slc(i), x_slc(i)], [0,0.5],'Color', 'r','lineWidth',2); hold off
    xlabel('X=x');  ylabel('P(X = x)'); title('Gaussian distribution');
    set(gca,'Fontsize',15);

    subplot(1,2,2)
    plot(x(1:idx_start), cumGauss(1:idx_start), 'Color', [13, 183, 200]./255, ...
        'lineWidth', 3); hold off; ylim([0,1]); xlim([x(1), x(end)]); hold on
    plot([x_slc(i), x_slc(i)], [0, cumGauss(idx_start)],'r-','lineWidth',2); hold on
    plot([x_slc(1), x_slc(i)], [cumGauss(idx_start), cumGauss(idx_start)],'r--','lineWidth',2); hold off
    xlabel('X=x');  ylabel('P(X < x)'); title('Cumulative gaussian distribution');
    set(gca,'Fontsize',15);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.55, 0.35]);
    pause(0.2)
end