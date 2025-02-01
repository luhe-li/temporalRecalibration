useCluster = 0;
for i_model = 1:6
    for sub = [1:4, 6,8]
        fit_recal_model(i_model, useCluster, sub)
    end
end