function [btst_d, lb_68CI, ub_68CI] = bootstrapDprime2(numBtst, resp, isOddball)

% output:
% btst_d: btst d'
% lb_68CI: lower bound of 68% CI of d'
% ub_68CI: higher bound of 68% CI of d'

nT = length(resp); 

for i = 1:numBtst
    % resample with replacement
    idx  = randi([1 nT],[1 nT]);
    %store the resampled responses
    btstResp = resp(idx);
    btstisOddball = isOddball(idx);
    % compute dprime
    HR = sum(btstisOddball & btstResp)/sum(btstisOddball);
    FAR = sum(~btstisOddball  &  btstResp)/sum(~btstisOddball);
    if HR == 1
        HR = 1 - (1/2*(HR+FAR));
    elseif HR == 0
        HR = 1/2*(HR+FAR);
    elseif FAR == 1
        FAR = 1 - (1/2*(HR+FAR));
    elseif FAR == 0
        FAR = 1/2*(HR+FAR);
    end
    btst_d(i) = norminv(HR,0,1) - norminv(FAR,0,1);
end

% obtain 68%CI
[lb_68CI, ub_68CI] = get68CI(btst_d);

end