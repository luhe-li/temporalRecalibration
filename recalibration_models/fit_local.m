useCluster = 0;
for i_model = 1:10
    for sub = [2:4, 6,8]

        fit_recal_model(i_model, useCluster, sub)

    end
end