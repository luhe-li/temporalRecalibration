%This function (1) bootstraps the data (i.e., resample with replacement), 
%(2) fits a psychometric function to resampled datasets, and (3) computes 
%confidence intervals on the estimated parameters 

%input arguments:
%s_unique : unique levels of SOA
%r_org    : organized responses (1 = V first, 2 = simultaneous, 3 = A first)
%           this matrix has size of length(s_unique) x nT 
%nT       : number of trials per level
%numBtst  : how many times do we want to bootstrap
%P_afirst : the anonymous function for a scaled psychometric function
%P_vfirst : the anonymous function for a scaled psychometric function
%P_simul  : the anonymous function for a scaled psychometric function
%lb       : the lower bound for parameter estimation
%           this vector has size of 1 x 3 (mu, sigma, lapse rate)
%ub       : the upper bound for parameter estimation
%options  : optional settings for fmincon

%output arguments:
%estP_btst: the estimated parameters for each resampled dataset 
%lb_68CI  : the lower bounds for the 68% confidence intervals
%ub_68CI  : the upper bounds for the 68% confidence intervals

function [estP_btst, minNLL, lb_68CI, ub_68CI] = bootstrapPost(pre_s_unique, post_s_unique, ...
    pre_r_org, post_r_org, pre_numTrials, post_numTrials,...
    numBtst, P_Afirst, P_Vfirst, P_simul, lb, ub, options) 
%number of free parameters
numP               = 5; 
%initialize outputs
estP_btst          = NaN(numBtst,numP);
minNLL             = NaN(1, numBtst);
[lb_68CI, ub_68CI] = deal(NaN(1,numP));
%an anonymous function that randomly picks an initial point given a lower
%and a upper bound
randInit           = @(val_lb, val_ub) rand(1)*(val_ub - val_lb) + val_lb;

for i = 1:numBtst
%     disp(i) %display the counter so you can see the progress
    %----------------------------------------------------------------------
    %                           PART 1.1: BOOTSTRAP pre 
    %%----------------------------------------------------------------------
    %initialize resampled responses
    pre_r_slc = NaN(length(pre_s_unique), pre_numTrials);
    %for each stimulus location, we resample the responses
    for j = 1:length(pre_s_unique)
        %randomly select trial indices (indices are allowed to occur more
        %than once, since we resample with replacement).
        idx        = randi([1 pre_numTrials],[1 pre_numTrials]);
        %store the resampled responses
        pre_r_slc(j,:) = pre_r_org(j,idx);
    end
    %compute the total number of comparison-more responses given each
    %intensity level
    pre_nT_Vfirst_slc=  sum(pre_r_slc == 1,2)';
    pre_nT_simultaneous_slc = sum(pre_r_slc == 2,2)';
    pre_nT_Afirst_slc = sum(pre_r_slc == 3,2)';

    %----------------------------------------------------------------------
    %                           PART 1.2: BOOTSTRAP post 
    %%----------------------------------------------------------------------
    %initialize resampled responses
    post_r_slc = NaN(length(post_s_unique), post_numTrials);
    %for each stimulus location, we resample the responses
    for j = 1:length(post_s_unique)
        %randomly select trial indices (indices are allowed to occur more
        %than once, since we resample with replacement).
        idx        = randi([1 post_numTrials],[1 post_numTrials]);
        %store the resampled responses
        post_r_slc(j,:) = post_r_org(j,idx);
    end
    %compute the total number of comparison-more responses given each
    %intensity level
    post_nT_Vfirst_slc = sum(post_r_slc == 1,2)';
    post_nT_simultaneous_slc= sum(post_r_slc == 2,2)';
    post_nT_Afirst_slc = sum(post_r_slc == 3,2)';

    %----------------------------------------------------------------------
    %                       PART 2: MODEL FITTING 
    %----------------------------------------------------------------------
    %use the anonymous function defined above randInt to generate random
    %initializations for the three free parameters
    init = arrayfun(@(idx) randInit(lb(idx), ub(idx)), 1:numP);
    %define nLL function given each bootstrapped dataset
    nLL = @(p) -pre_nT_Afirst_slc*log(P_Afirst(pre_s_unique, p(1), p(3), p(4), p(5)))' ...
    -pre_nT_Vfirst_slc*log(P_Vfirst(pre_s_unique, p(1), p(3), p(4), p(5)))'...
    -pre_nT_simultaneous_slc * log(P_simul(pre_s_unique, p(1), p(3), p(4), p(5)))'...
    -post_nT_Afirst_slc*log(P_Afirst(post_s_unique, p(2), p(3), p(4), p(5)))' ...
    -post_nT_Vfirst_slc*log(P_Vfirst(post_s_unique, p(2), p(3), p(4), p(5)))'...
    -post_nT_simultaneous_slc * log(P_simul(post_s_unique, p(2), p(3), p(4), p(5)))';

    %fit a psychometric curve to each bootstrapped dataset
    [estP_btst(i,:), minNLL(i)] = fmincon(nLL,init,[],[],[],[],lb,ub,[],options);     
end

%--------------------------------------------------------------------------
%                    PART 3: FINDING CONFIDENCE INTERVALS
%--------------------------------------------------------------------------
for j = 1:numP
    [lb_68CI(j), ub_68CI(j)] = get68CI(estP_btst(:,j));
end     
    