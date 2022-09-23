close all; clc; clear all;

pre_ms_unique = [-0.500000000000000	-0.300000000000000	-0.250000000000000	-0.200000000000000	-0.150000000000000	-0.100000000000000	-0.0500000000000000	0	0.0500000000000000	0.100000000000000	0.150000000000000	0.200000000000000	0.250000000000000	0.300000000000000	0.500000000000000]*1e3;

%make a finer grid for the timing difference between the auditory and the
%visual stimulus
SOA_finer = pre_ms_unique;

% define PMF
P_Afirst = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig);
P_Vfirst = @(SOA, mu, sig, c, lambda) lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig));
P_simultaneous = @(SOA, mu, sig, c, lambda) ...
    1 - (lambda/3 + (1-lambda).*normcdf(-c, SOA - mu, sig)) ...
    - (lambda/3 + (1-lambda).*(1 - normcdf(c, SOA-mu, sig)));

M1 = {@(p) P_Afirst(SOA_finer, p(1), p(3), p(5), p(7));...
    @(p) P_Vfirst(SOA_finer, p(1), p(3), p(5), p(7));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(5), p(7));...
    @(p) P_Afirst(SOA_finer, p(2), p(4), p(6), p(8));...
    @(p) P_Vfirst(SOA_finer, p(2), p(4), p(6), p(8));...
    @(p) P_simultaneous(SOA_finer, p(2), p(4), p(6), p(8))};

M2 = {@(p) P_Afirst(SOA_finer, p(1), p(3), p(4), p(6));...
    @(p) P_Vfirst(SOA_finer, p(1), p(3), p(4), p(6));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(4), p(6));...
    @(p) P_Afirst(SOA_finer, p(2), p(3), p(5), p(7));...
    @(p) P_Vfirst(SOA_finer, p(2), p(3), p(5), p(7));...
    @(p) P_simultaneous(SOA_finer, p(2), p(3), p(5), p(7))};

M3 = {@(p) P_Afirst(SOA_finer, p(1), p(3), p(5), p(6));...
    @(p) P_Vfirst(SOA_finer, p(1), p(3), p(5), p(6));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(5), p(6));...
    @(p) P_Afirst(SOA_finer, p(2), p(4), p(5), p(7));...
    @(p) P_Vfirst(SOA_finer, p(2), p(4), p(5), p(7));...
    @(p) P_simultaneous(SOA_finer, p(2), p(4), p(5), p(7))};

M4 = {@(p) P_Afirst(SOA_finer, p(1), p(3), p(4), p(5));...
    @(p) P_Vfirst(SOA_finer, p(1), p(3), p(4), p(5));...
    @(p) P_simultaneous(SOA_finer, p(1), p(3), p(4), p(5));...
    @(p) P_Afirst(SOA_finer, p(2), p(3), p(4), p(6));...
    @(p) P_Vfirst(SOA_finer, p(2), p(3), p(4), p(6));...
    @(p) P_simultaneous(SOA_finer, p(2), p(3), p(4), p(6))};

M5 = {@(p) P_Afirst(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_Vfirst(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_Afirst(SOA_finer, p(1), p(2), p(4), p(6));...
    @(p) P_Vfirst(SOA_finer, p(1), p(2), p(4), p(6));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(4), p(6))};

M6 = {@(p) P_Afirst(SOA_finer, p(1), p(2), p(3), p(4));...
    @(p) P_Vfirst(SOA_finer, p(1), p(2), p(3), p(4));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(3), p(4));...
    @(p) P_Afirst(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_Vfirst(SOA_finer, p(1), p(2), p(3), p(5));...
    @(p) P_simultaneous(SOA_finer, p(1), p(2), p(3), p(5))};

models = {M1; M2; M3; M4; M5; M6};

%% plug in best-fitting parameters and plot

% plot raw data
cMAP1 = [229, 158, 168; 203, 227, 172; 171,223,235]./255;
cMAP2 = [216, 49, 91; 175, 213, 128; 88,193,238]./255;

bestPs = {[42.84, -40, 66.131, 49.306, 109.76, 82.858, 0.01, 0.06];
    [42.84, -20.277, 66.131, 109.76, 82.858, 0.01, 0.06];
    [42.84, -20.277, 66.131, 49.306, 109.76, 0.01, 0.06];
    [42.84, -20.277, 66.131, 109.76, 0.01, 0.06];
    [42.84, 66.131, 109.76, 82.858, 0.01, 0.06];
    [42.84, 66.131, 109.76, 0.01, 0.06]};

% plot fitting curves

for m = 4:4
    bestP = bestPs{m};
    Pre_Afirst_fit = models{m}{1}(bestP);
    Pre_Vfirst_fit = models{m}{2}(bestP);
    Pre_simul_fit = models{m}{3}(bestP);
    Post_Afirst_fit = models{m}{4}(bestP);
    Post_Vfirst_fit = models{m}{5}(bestP);
    Post_simul_fit = models{m}{6}(bestP);

    figure; hold on
    set(gca, 'LineWidth', 2, 'FontSize', 20)
    set(gcf, 'Position',[10 10 500 400])
    % look better
    xlabel('test lag')
    ylabel('proportion of responses')
    set(gca,'TickDir','out');
    xlim([SOA_finer(1), SOA_finer(end)])
    xticks(0)

%     plot(SOA_finer, Pre_Vfirst_fit, 'o','w');
    saveas(gca,'f0','png')

    plot(SOA_finer, Pre_Vfirst_fit, 'o','Color', cMAP2(1,:), 'MarkerSize',12,'lineWidth', 3)
    saveas(gca,'f1','png')
    
    plot(SOA_finer, Pre_Afirst_fit, 'o','Color', cMAP2(3,:),'MarkerSize',12,'lineWidth', 3)
    saveas(gca,'f2','png')

    plot(SOA_finer, Pre_simul_fit, 'o','Color', cMAP2(2,:), 'MarkerSize',12,'lineWidth', 3)
    saveas(gca,'f3','png')

    plot(SOA_finer, Post_Vfirst_fit, 'o','MarkerFaceColor', cMAP2(1,:),'MarkerEdgeColor', cMAP2(1,:), 'MarkerSize',12,'lineWidth', 1)
    plot(SOA_finer, Post_Afirst_fit, 'o','MarkerFaceColor', cMAP2(3,:),'MarkerEdgeColor', cMAP2(3,:),'MarkerSize',12,'lineWidth', 1)
    plot(SOA_finer, Post_simul_fit, 'o','MarkerFaceColor', cMAP2(2,:),'MarkerEdgeColor', cMAP2(2,:),'MarkerSize',12,'lineWidth', 1)
    saveas(gca,'f4','png')
    %     % display best fitting parameters
    %     annotation( 'textbox', 'String', round(bestP,3), ...
    %         'FontSize', 14, 'Units', 'normalized', 'EdgeColor', 'none', ...
    %         'Position', [0.9,0.9,0.03,0.03], 'Interpreter','Latex')

    legend('"vision first"','"audition first"','"simultaneous"','"vision first"','"audition first"','"simultaneous"')


    saveas(gca,'f5','png')

end