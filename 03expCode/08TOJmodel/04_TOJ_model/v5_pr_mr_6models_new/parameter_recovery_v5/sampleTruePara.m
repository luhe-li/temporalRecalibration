% This function summarizes mean and SD for each parameter, each model, from
% available subjects and all sessions, and samples a 'true parameter' from
% the corresponding distribution

% outputs
% realP_mean    : mean of parameter, realP_mean{m}{p} is a scalar
% realP_SD      : SD of parameter, realP_SD{m}{p} is a scalar
% sampledP      : sampledP{m}(1,p) is a scalar

function [realP_mean, realP_SD, sampledP] = sampleTruePara

%% summarize fitted parameters

% load model fitted parameters
load('a04_v5_mc_results.mat')

% for each model, extract each parameter from ava_sub, 9 sessions
nModel = 6;
nPara = [8,7,7,6,6,5];
sub2use = [3,4,6,7,8,9,10];
nSub = numel(sub2use);
nSes = 9;
realP = cell(1, nModel);

% Only mu can be negative. Create a matrix for the
% index of mu in each model
mu_idx = [1,1;1,2;2,1;2,2;3,1;3,2;4,1;4,2;5,1;6,1];

for m = 1:nModel

    % number of parameters for this model
    m_nPara = nPara(m);

    % initiate parameter matrix (i_para x session x ava_sub)
    m_p = NaN(m_nPara, nSes, numel(sub2use));

    for p = 1:m_nPara
        for s = 1:nSub
            % for each subject, store a matrix (i_para x session)
            m_p(:,:,s) = reshape([estP{sub2use(s)}{:, m}], [], 9);
        end

        % reorganize matrix
        % for each parameter, get a matrix (session x sub)
        realP{m}{p} = squeeze(m_p(p, :, :));

    end
end

% calculate mean and sd for each model, each parameter
[realP_mean, realP_SD, sampledP] = deal(cell(1, nModel));

for m = 1:nModel

    % number of parameters for this model
    m_nPara = nPara(m);

    for p = 1:m_nPara
        if ismember([m,p], mu_idx,'rows') % if current parameter is mu, which can be negative, summarize them assuming they are from a normal dist
            realP_mean{m}{p} = mean(realP{m}{p},'all');
            realP_SD{m}{p} = std(realP{m}{p},[],'all');

        else % if current parameter is sigma/criterion/lambda, firstly take its log, assume the distribution is normal, take the mean and sd
            realP_mean{m}{p} = mean(log(realP{m}{p}),'all');
            realP_SD{m}{p} = std(log(realP{m}{p}),[],'all');

        end
    end

end


%% sample simulated parameters for each model from fitted parameters

for m = 1:nModel

    % number of parameters for this model
    m_nPara = nPara(m);

    for p = 1:m_nPara
        % extract mean and sd from fitted parameters
        temp_mean = realP_mean{m}{p};
        temp_sd = realP_SD{m}{p};

        % if current parameter is mu, which can be negative, sample from a normal distribution
        if ismember([m,p], mu_idx,'rows')
            temp_sampledP = randn * temp_sd + temp_mean;

            % if current parameter is sigma/criterion/lambda, sample log(p)
            % from a normal dist, and transfer back p = exp(log(p)), so
            % that p is always positive
        else
            temp_sampledP = exp(randn * temp_sd + temp_mean);

        end
        % store sampled sampled parameters
        sampledP{m}(1,p) = temp_sampledP;
    end
end

end