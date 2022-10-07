function [hitRate, FA]    = summarizePerformance(currentTrial, nTotalTrials, idxOddballA, idxOddballV, response)
    % sum up oddball presence and response
boolOddballA              = false(1, nTotalTrials);
boolOddballA(idxOddballA) = true; % oddball presence
boolRespA                 = ismember(response,[2,3]); % oddball response
boolOddballA              = boolOddballA(1:currentTrial); 
boolRespA                 = boolRespA(1:currentTrial); 
% repeat for V
boolOddballV              = false(1, nTotalTrials);
boolOddballV(idxOddballV) = true; % oddball presence
boolRespV                 = ismember(response,[2,1]);% oddball response
boolOddballV              = boolOddballV(1:currentTrial); 
boolRespV                 = boolRespV(1:currentTrial); 
% auditory summary
hitRate(1,1)              = sum(boolOddballA & boolRespA)/sum(boolOddballA);
FA(1,1)                   = sum(~boolOddballA  &  boolRespA)/sum(~boolOddballA); %sum((boolRespA - boolOddballA) == 1)/sum(~boolOddballA); 
% visual summary
hitRate(1,2)              = sum(boolOddballV & boolRespV)/sum(boolOddballV);
FA(1,2)                   = sum(~boolOddballV & boolRespV)/sum(~boolOddballV); %sum((boolRespV - boolOddballV) == 1)/sum(~boolOddballV);
end