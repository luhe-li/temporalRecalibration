% This script takes the fitting results from a1_modelFitting.m and checks
% parameter trade-offs by plotting fitted parameters between each other.

load("a1_modelFittingResults.mat")

%% plot free parameters in the exposure phase
figure; hold on
set(gcf, 'Position', get(0, 'Screensize'));

combi = nchoosek(7:11,2);
n_combi = size(combi, 1);
n_col = max(factor(n_combi));
n_row = n_combi/n_col;

for i = 1:n_combi
    subplot(n_row, n_col, i); hold on
    set(gca, 'LineWidth', 1, 'FontSize', 15)
    axis square
    x = combi(i,1);
    y = combi(i,2);

    plot(model.estimatedP(:,x), model.estimatedP(:,y),'o','Color',[repmat(0.3, 1, 3)])
    [r, p] = corr(model.estimatedP(:,x), model.estimatedP(:,y));
    r2 = r^2;
    title(['r^2 = ' num2str(round(r2,2)) ', p = ' num2str(round(p,2))])
    xlabel(model.paraID{x})
    ylabel(model.paraID{y});

end