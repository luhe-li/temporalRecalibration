function [hitRate, FA, btst_d, lb_68CI, ub_68CI]  = summarizeBlockPerformance2(startTrial, endTrial, idxOddballA, idxOddballV, response)
% sum up oddball presence and response
isOddballA              = zeros(1, length(response));
isOddballA(idxOddballA) = 1; % oddball presence
RespA                   = ismember(response,[2,3]); % oddball response
isOddballA              = isOddballA(startTrial:endTrial); 
RespA                   = RespA(startTrial:endTrial); 
% repeat for V
isOddballV              = zeros(1, length(response));
isOddballV(idxOddballV) = 1; % oddball presence
RespV                   = ismember(response,[2,1]);% oddball response
isOddballV              = isOddballV(startTrial:endTrial); 
RespV                   = RespV(startTrial:endTrial); 
% auditory summary
hitRate(1,1)            = sum(isOddballA & RespA)/sum(isOddballA);
FA(1,1)                 = sum(~isOddballA  &  RespA)/sum(~isOddballA); %sum((boolRespA - boolOddballA) == 1)/sum(~boolOddballA); 
% visual summary
hitRate(1,2)            = sum(isOddballV & RespV)/sum(isOddballV);
FA(1,2)                 = sum(~isOddballV & RespV)/sum(~isOddballV); %sum((boolRespV - boolOddballV) == 1)/sum(~boolOddballV);

%% bootstrap d'
numBtst = 1e3;
[btst_d(1,:), lb_68CI(1), ub_68CI(1)] = bootstrapDprime2(numBtst, RespA, isOddballA);
[btst_d(2,:), lb_68CI(2), ub_68CI(2)] = bootstrapDprime2(numBtst, RespV, isOddballV);

end