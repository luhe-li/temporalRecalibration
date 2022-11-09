load("a1_modelFit_results_y_pre_mu.mat")

f1 = figure; hold on
histogram(model.minNLL)
f2 = figure; hold on
for p = 1:model.numPara
    subplot(3,4,p)
    histogram(model.estimatedP(:,p),10)
end

%find the index that corresponds to the minimum nLL
[min_nll1, min_idx] = min(model(sub, ses).minNLL);
%find the corresponding estimated parameters
p1 = model(sub, ses).estimatedP(min_idx,:);
disp([min_nll1, p1])

load("a1_modelFit_results_n_pre_mu.mat")
figure(f1)
histogram(model.minNLL)
figure(f2);
for p = 1:model.numPara
    subplot(3,4,p); hold on
    histogram(model.estimatedP(:,p),10)
end

%find the index that corresponds to the minimum nLL
[min_nll2, min_idx] = min(model(sub, ses).minNLL);
%find the corresponding estimated parameters
p2 = model(sub, ses).estimatedP(min_idx,:);
disp([min_nll2, p2])