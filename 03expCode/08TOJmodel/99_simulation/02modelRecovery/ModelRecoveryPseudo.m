for count = 1:100
    % use fixed model paramaters 
    % (Is it better to sample model parameters under the prior
    % distribution? If so, how to decide prior?)
    for sim = 1:2
        % simulate data of the real/null model
        for model = 1:2
            % fit real/null models
        end
        % select the best model with lowest AIC
        % update confusion matrix for each simulation
    end
end