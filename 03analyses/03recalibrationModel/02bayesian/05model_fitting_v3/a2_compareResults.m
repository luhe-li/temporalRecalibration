load("a1_modelFit_results_y_pre_mu.mat")

f1 = figure; hold on
histogram(model.minNLL)
f2 = figure; hold on
for p = 1:model.numPara
    subplot(3,4,p)
    histogram(model.estimatedP(:,p),10)
end

load("a1_modelFit_results_n_pre_mu.mat")
figure(f1)
histogram(model.minNLL)
figure(f2);
for p = 1:model.numPara
    subplot(3,4,p); hold on
    histogram(model.estimatedP(:,p),10)
end

