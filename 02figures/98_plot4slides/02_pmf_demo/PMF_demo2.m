
clear all; clc; close all;

% plot parameters
cmp =  [216, 49, 91; 175, 213, 128; 88,193,238]./255;
mksz = 10;

% set parameters
min_x = -1000;
max_x = 1000;
x_axis = [min_x:max_x];
soas = [-500:100:500];
sig = 200;
mu = 0;
c = 150;
lc = mu-c;
rc = mu+c;

figure; hold on
set(gcf, 'Position',[10 10 1000 400])
set(gcf, 'color', 'white')
flnm = 'PMF_demo.gif';

for i = 1:numel(soas)

    soa = soas(i);

    % obtain lines
    gauss = normpdf(x_axis, soa, sig);
    p_resp(1) = 1 - normcdf(c, soa-mu, sig); % p(v_lead)
    p_resp(3) = normcdf(-c, soa - mu, sig); % p(a_lead)
    p_resp(2) = 1 - p_resp(1) - p_resp(3); % p(simul)


    %%
    subplot 121;
    set(gca, 'LineWidth', 2, 'FontSize', 30); hold on

    % plot(x_axis, gauss)
    xline(lc,'LineWidth',2)
    xline(rc,'LineWidth',2)

    % delete previous patches
    if i >1
        delete([area1, area2, area3])
    end

    % plot area left to -c
    idx_lc = find(x_axis == lc);
    idx_rc = find(x_axis == rc);
    area1 = patch([min_x:lc, fliplr(min_x:lc)], ...
        [zeros(1, numel(min_x:lc)), fliplr(gauss(1: idx_lc))],cmp(3,:), ...
        'EdgeColor','none');

    % plot area between -c and c
    area2 = patch([lc:rc, fliplr(lc:rc)], ...
        [zeros(1, numel(lc:rc)), fliplr(gauss(idx_lc:idx_rc))],cmp(2,:), ...
        'EdgeColor','none');

    % plot area right to c
    area3 = patch([rc:max_x, fliplr(rc:max_x)], ...
        [zeros(1, numel(rc:max_x)), fliplr(gauss(idx_rc:end))], ...
        cmp(1,:),'EdgeColor','none');

    % look better
    ylim([0 max(gauss)*1.5])
    % tick options
    xticks(soa)
    xticklabels('lag')
    set(gca,'TickDir','out');
    ax = gca;
    ax.XAxis.LineWidth = 2;
    ax.XAxis.TickLength = [0.025 0.025];
    ax.YAxis.Visible = 'off';
    hold off

    %%
    subplot 122; hold on
    set(gca, 'LineWidth', 2, 'FontSize', 30)

    for ii = 1:3
        plot(soa, p_resp(ii),'o',...
            'MarkerSize',mksz,'MarkerFaceColor',cmp(ii,:),'MarkerEdgeColor','none')
    end

    % look better
    yticks([0 1])
    ylabel('probability')
    xlim([min(soas)-100,max(soas)+100])
    % tick options
    xticks(soa)
    xticklabels('lag')
    set(gca,'TickDir','out');
    ax = gca;
    ax.XAxis.TickLength = [0.025 0.025];

    % write gif
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
        imwrite(imind,cm,flnm,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,flnm,'gif','WriteMode','append');
    end

end