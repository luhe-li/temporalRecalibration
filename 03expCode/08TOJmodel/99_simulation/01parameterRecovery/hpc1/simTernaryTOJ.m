%This function simulates the probability of reporting 'V-first' responses
%in a temporal order judgment task (TOJ)
%It can also simulate the 'A-first' or 'simutaneous' response given the
%corect p_model

%input arguments:
%lenS     : the number of temporal discrepancies
%nT       : the number of trials for each temporal discrepancy
%p_model  : the probability of reporting 'V-first' defined by our model
%           given each of the temporal discrepancy
%           (this is a row vector: size = 1 x lenS)

%output arguments:
%bool_V1st : raw data of reporting 'V-first' for each trial
%sim_prob_V1st : simulated probability of reporting 'V-first' responses
%                given each of the temporal discrepancies
%                (this is a row vector: size = 1 x lenS)

function [bool_V1st, sim_prob_V1st,bool_A1st, sim_prob_A1st, bool_simul, sim_prob_simul]...
    = simTernaryTOJ(lenS, nT, p_model_V1st, p_model_A1st, p_model_simul)
    %----------------------------------------------------------------------
    % YOUR CODE STARTS HERE
    % Hint: you will probably find the function rand.m useful
    % if you set M = rand(1,10), and then compute how many elements in M
    % is greater than p=0.1, i.e., sum(M < 0.1)/10, what do you get?
    
    % if we increase trials to 1000, and print out numel(M < 0.1)/1000,
    % you would get a number closer to 0.1.
    
    % In other words, rand.m gives you noisy data when the number of trials
    % is small.
    
    M_rand        = rand(lenS, nT);
    bool_V1st     = M_rand < repmat(p_model_V1st',[1, nT]);
    sim_prob_V1st = sum(bool_V1st,2)/nT;
    bool_simul = M_rand > repmat((1-p_model_simul)',[1, nT]);
    sim_prob_simul = sum(bool_simul,2)/nT;
    bool_A1st = (M_rand >= repmat(p_model_V1st',[1, nT]))...
        &(M_rand <= repmat((1-p_model_simul)',[1, nT]));
    sim_prob_A1st = sum(bool_A1st,2)/nT;
    %----------------------------------------------------------------------
end