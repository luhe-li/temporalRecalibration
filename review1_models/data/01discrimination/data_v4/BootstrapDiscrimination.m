%This function (1) bootstraps the data (i.e., resample with replacement), 
%(2) fits a psychometric function to resampled datasets, and (3) computes 
%confidence intervals on the estimated parameters 

%input arguments:
%s_unique : unique levels of SOA
%r_org    : organized responses (1: standard-first; 2: comparison-first)
%           this matrix has size of length(s_unique) x nT 
%nT       : number of trials per level
%numBtst  : how many times do we want to bootstrap
%Pcomp   : the anonymous function for a scaled psychometric function, for
%           comparison is perceived more intense
%lb       : the lower bound for parameter estimation
%           this vector has size of 1 x 3 (mu, sigma, lapse rate)
%ub       : the upper bound for parameter estimation
%options  : optional settings for fmincon

%output arguments:
%estP_btst: the estimated parameters for each resampled dataset 
%lb_95CI  : the lower bounds for the 95% confidence intervals
%ub_95CI  : the upper bounds for the 95% confidence intervals

function [estP_btst, minNLL, lb_95CI, ub_95CI, nT_compMore_slc] = BootstrapDiscrimination(s_unique, r_org, nT,...
    numBtst, Pcomp, lb, ub, options) 

%number of free parameters
numP               = length(lb); 
%initialize outputs
estP_btst          = NaN(numBtst,numP);
minNLL             = NaN(1, numBtst);
[lb_95CI, ub_95CI] = deal(NaN(1,numP));
%an anonymous function that randomly picks an initial point given a lower
%and a upper bound
randInit           = @(val_lb, val_ub) rand(1)*(val_ub - val_lb) + val_lb;

for i = 1:numBtst
    disp(i) %display the counter so you can see the progress
    %----------------------------------------------------------------------
    %                           PART 1: BOOTSTRAP
    %%----------------------------------------------------------------------
	%initialize resampled responses
    r_slc = NaN(length(s_unique), nT);
    %for each stimulus location, we resample the responses
    for j = 1:length(s_unique)
        %randomly select trial indices (indices are allowed to occur more
        %than once, since we resample with replacement).
        %Hint: you'll find function randi.m useful
        idx        = randi([1 nT],[1 nT]);
        %store the resampled responses
        r_slc(j,:) = r_org(j,idx);
    end
    %compute the total number of comparison-more responses given each
    %intensity level
    nT_compMore_slc(i,:) = sum(r_slc == 2,2)';

    %----------------------------------------------------------------------
    %                       PART 2: MODEL FITTING 
    %----------------------------------------------------------------------
    %use the anonymous function defined above randInt to generate random
    %initializations for the three free parameters
    init = arrayfun(@(idx) randInit(lb(idx), ub(idx)), 1:numP);
    %define nLL function given each bootstrapped dataset
    nLL = @(p) -nT_compMore_slc(i,:)*log(Pcomp(s_unique, p(1), p(2), p(3), p(4)))'...
            -(nT - nT_compMore_slc(i,:))*log(1-Pcomp(s_unique, p(1), p(2), p(3), p(4)))';
    %fit a psychometric curve to each bootstrapped dataset
    [estP_btst(i,:), minNLL(i)] = fmincon(nLL,init,[],[],[],[],lb,ub,[],options);     
end

%--------------------------------------------------------------------------
%                    PART 3: FINDING CONFIDENCE INTERVALS
%--------------------------------------------------------------------------
for j = 1:numP
    [lb_95CI(j), ub_95CI(j)] = get68CI(estP_btst(:,j));
end

%This function takes a vector v, sorts the entries in an ascending order,
%and then returns the entries that correspond the 2.5% and 97.5% perentiles
function [CI_lb, CI_ub] = get95CI(v)
%first sort the vector in an ascending order
n_sorted = sort(v);
%compute how many entries the vector has
lenN     = length(v);
%lower bound 
CI_ub    = n_sorted(ceil(lenN*0.975));
%upper bound
CI_lb    = n_sorted(floor(lenN*0.025));
