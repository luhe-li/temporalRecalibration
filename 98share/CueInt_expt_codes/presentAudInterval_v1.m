function [order, RT, ExpInfo] = presentAudInterval_v1(t,...
    ExpInfo,ScreenInfo, VSinfo, AudInfo, pahandle, windowPtr)

%----------------------------------------------------------------------
%-------------------------- create stimuli ----------------------------
%----------------------------------------------------------------------

% jitter fixation and ITI duration, record it
jitter = @(lb, ub) (ub - lb)*rand + lb;
ExpInfo.jFixation(t) = jitter(ExpInfo.fixation(1),ExpInfo.fixation(2));
ExpInfo.jISI(t) = jitter(ExpInfo.ISI(1), ExpInfo.ISI(2));
ExpInfo.jITI(t) = jitter(ExpInfo.ITI(1), ExpInfo.ITI(2));

% generate sound samples
[audSamples, amp, ExpInfo] = generateSound(t, ExpInfo);

% create buffer
PsychPortAudio('FillBuffer', pahandle, [audSamples; zeros(size(audSamples))]);
% PsychPortAudio('FillBuffer', pahandle, audSamples);

%----------------------------------------------------------------------
%--------------------- display interval stimuli -----------------------
%----------------------------------------------------------------------

% allow top priority for better temporal accuracy
Priority(ScreenInfo.topPriorityLevel);

% fixation
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x1_lb,...
    ScreenInfo.y1_lb, ScreenInfo.x1_ub, ScreenInfo.y1_ub]);
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x2_lb,...
    ScreenInfo.y2_lb, ScreenInfo.x2_ub, ScreenInfo.y2_ub]);
Screen('Flip',windowPtr); 
WaitSecs(ExpInfo.jFixation(t));
Screen('Flip',windowPtr);

% start by playing 
PsychPortAudio('Start', pahandle);

% wait till the playback stops
PsychPortAudio('Stop', pahandle, 1); 

%----------------------------------------------------------------------
%---------------------------- record response -------------------------
%----------------------------------------------------------------------

% response cue
DrawFormattedText(windowPtr, 'First longer: 1\n Second longer: 2\n',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);

% wait for response
resp=1; tic;
while resp
    [~, ~, keyCode] = KbCheck();
    if keyCode(KbName('numLock'))
        ShowCursor;
        sca;
        error('Escape');
    elseif keyCode(KbName('1'))
        order = 1;
        resp=0;
    elseif keyCode(KbName('2'))
        order = 2;
        resp=0;
    end
end
RT   = toc;
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jITI(t)); %ITI

end