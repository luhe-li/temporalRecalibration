%%
% This script illustratea all possible combinations of RGB colors in Matlab
% IV: number of levels

%%
close all
num = 5; % number of levels within the range 0:1
sz = 30; % size of the color patch

rgb_vec = linspace(0,1,num);
rgb_m = combvec(combvec(rgb_vec,rgb_vec),rgb_vec);
num_comb = size(rgb_m,2);

for nn = 1:num_comb
    if mod(nn,num^2) == 1, figure('Position',[0, 300, 15*sz, 15*sz]), end
    ind = mod(nn,num^2);
    if ind == 0,ind = num^2;end
    
    subplot(num,num,ind)
    plot(1,1,'-o','MarkerFaceColor',rgb_m(:,nn),'MarkerEdgeColor','k','MarkerSize',20)
    axis off
    title(sprintf('[%.2f %.2f %.2f]',rgb_m(1,nn),rgb_m(2,nn),rgb_m(3,nn)))
    
    if mod(nn,num^2) == 0
        fileName = [pwd,'/',num2str(round(nn/num^2))];
        saveas(gcf,fileName, 'jpg')
    end
end

% close all
