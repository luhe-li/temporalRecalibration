
clear all; close all; clc;
cMAP = [200, 40, 40; 255, 128, 0; 13, 183, 200]./255;

% select mat file here
load('PR_determineNumTrials_04-May-2022 16:47:07.mat')
D = Data{1,4};
P = Data{1,2};

valName = {'\mu_{pre}','\mu_{post}', '\sigma', '\lambda', 'criterion'};
realP = [P.muPrePost, P.sigma_deltaT, P.lapse, P.c];
numP = length(realP);

lb = [ -150, -150,  10, 1e-2, 10];
ub = [  150,  150, 150, 0.06, 300];
numBtst = 1000;

% select dataset here
t = 4; % i_nTrial
j = 4; % i_dataset

i_estP_btst = D.estP_btst{1,t}{1,j};
i_CI_btst_lb = D.CI_btst_lb{1,t}{1,j};
i_CI_btst_ub = D.CI_btst_ub{1,t}{1,j};

figure; set(gcf,'position',[0,0,400,800]); hold on

for i = 1:numP
    subplot(numP,1,i)
    histogram(i_estP_btst(:,i),'FaceColor', cMAP(2,:), ...
        'FaceAlpha', 0.5, 'EdgeColor',cMAP(2,:)); hold on
    plot([realP(i), realP(i)], [0, numBtst*0.3],'r--', 'lineWidth',3); hold on
    plot([i_CI_btst_lb(:,i), i_CI_btst_ub(:,i)],[5, 5],'-k','LineWidth', 2); hold on
    xlim([lb(i), ub(i)]); ylim([0, numBtst*0.3]);
    xticks(sort(unique([linspace(lb(i), ub(i), 5), realP(i)])));
    title(valName{i}); box off
    set(gca,'FontSize',15);
end

% find out btstp for extreme case
% mu_diff = i_estP_btst(:,2) - i_estP_btst(:,1);
% figure
% histogram(mu_diff)
% [m, idx1] = min(mu_diff);
% [m, idx2] = max(mu_diff);
% i_estP_btst(idx1,:)
% i_estP_btst(idx2,:)