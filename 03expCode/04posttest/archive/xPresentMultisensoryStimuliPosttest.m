function [order, RT] = PresentMultisensoryStimuliPosttest(TrialNum,...
    ExpInfo,ScreenInfo,VSinfo, AudInfo,pahandle,windowPtr)

% -------------------------------------------------------------
%-----------Calculate the coordinates of the target stimuli------------
%----------------------------------------------------------------------
%Make visual stimuli
blob_coordinates = [ScreenInfo.xmid, ScreenInfo.liftingYaxis];
dotCloud = generateOneBlob(windowPtr,blob_coordinates,VSinfo,ScreenInfo);

%----------------------------------------------------------------------
%--------------------- re-exposure ----------------------
%----------------------------------------------------------------------
jitter = @(lb, ub) (ub - lb)*rand + lb;

% fixation
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x1_lb,...
    ScreenInfo.y1_lb, ScreenInfo.x1_ub, ScreenInfo.y1_ub]);
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x2_lb,...
    ScreenInfo.y2_lb, ScreenInfo.x2_ub, ScreenInfo.y2_ub]);
Screen('Flip',windowPtr); WaitSecs(jitter(ExpInfo.fixation(1),ExpInfo.fixation(2)));

% blank screen 1
Screen('Flip',windowPtr); WaitSecs(jitter(ExpInfo.blankDuration1(1), ExpInfo.blankDuration1(2)));

% top-up adaptor x nRep
otoi1 = abs(ExpInfo.adaptor) - ExpInfo.stimDura; % offset-to-onset interval
for tt = 1:ExpInfo.nRepAdaptor
    if ExpInfo.adaptor > 0 %present V first
        % visual stimulus
        Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        onset = Screen('Flip',windowPtr); % onset
        offset = Screen('Flip',windowPtr, onset + VSinfo.duration - ScreenInfo.ifi/2); % offset
        % ISI
        WaitSecs(otoi1 - ScreenInfo.ifi/2);
        %start playing the sound
        PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
        PsychPortAudio('Start', pahandle, 1, 0, 0);
        WaitSecs(AudInfo.stimDura);
        PsychPortAudio('Stop', pahandle);
        
    elseif ExpInfo.adaptor < 0 %present A first
        %start playing the sound
        PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
        PsychPortAudio('Start', pahandle, 1, 0, 0);
        WaitSecs(AudInfo.stimDura);
        PsychPortAudio('Stop', pahandle);
        % ISI
        WaitSecs(otoi1 - ScreenInfo.ifi/2);
        %present V second
        Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        onset = Screen('Flip',windowPtr); % onset
        offset = Screen('Flip',windowPtr, onset + VSinfo.duration - ScreenInfo.ifi/2); % offset
        
    else %present A and V at the same time
        %start playing the sound while displaying the visual stimulus
        PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
        PsychPortAudio('Start', pahandle, 1, 0, 0); % A onset
        Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        onset = Screen('Flip',windowPtr); % V onset
        offset = Screen('Flip',windowPtr, onset + VSinfo.duration - ScreenInfo.ifi/2); % V offset
        PsychPortAudio('Stop', pahandle); % A offset
    end
    
    % ISI between adaptors
    Screen('Flip',windowPtr); WaitSecs(jitter(ExpInfo.aISI(1), ExpInfo.aISI(2)));
    
end

% blank screen 2
Screen('Flip',windowPtr); WaitSecs(jitter(ExpInfo.blankDuration2(1), ExpInfo.blankDuration2(2)));

%----------------------------------------------------------------------
%--------------------- post-test ----------------------
%----------------------------------------------------------------------

% present test stimuli
trialSOA = ExpInfo.trialSOA(TrialNum);
otoi2 = abs(trialSOA) - VSinfo.duration; % offset-to-onset interval
if trialSOA > 0 %present V first
    % visual stimulus
    Screen('DrawTexture',windowPtr,dotCloud,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    onset = Screen('Flip',windowPtr); % onset
    offset = Screen('Flip',windowPtr, onset + VSinfo.duration - ScreenInfo.ifi/2); % offset
    % offset-to-onset interval
    WaitSecs(otoi2 - ScreenInfo.ifi/2);
    %start playing the sound
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    PsychPortAudio('Start', pahandle, 1, 0, 0);
    WaitSecs(AudInfo.stimDura);
    PsychPortAudio('Stop', pahandle);
    
elseif trialSOA < 0 %present A first
    %start playing the sound
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    PsychPortAudio('Start', pahandle, 1, 0, 0);
    WaitSecs(AudInfo.stimDura);
    PsychPortAudio('Stop', pahandle);
    % offset-to-onset interval
    WaitSecs(otoi2 - ScreenInfo.ifi/2);
    %present V second
    Screen('DrawTexture',windowPtr,dotCloud,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    onset = Screen('Flip',windowPtr); % onset
    offset = Screen('Flip',windowPtr, onset + VSinfo.duration - ScreenInfo.ifi/2); % offset
    
else %present A and V at the same time
    %start playing the sound while displaying the visual stimulus
    %     PsychPortAudio('FillBuffer', pahandle, AudInfo.GaussianWhiteNoise);
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    PsychPortAudio('Start', pahandle, 1, 0, 0); % A onset
    Screen('DrawTexture',windowPtr,dotCloud,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    onset = Screen('Flip',windowPtr); % V onset
    offset = Screen('Flip',windowPtr, onset + VSinfo.duration - ScreenInfo.ifi/2); % V offset
    PsychPortAudio('Stop', pahandle); % A offset
end

%----------------------------------------------------------------------
%--------------Record response ---------------
%----------------------------------------------------------------------
% left arrow: V-first
% down arrow: simultaneous
% right arrow: A-first
KbDevices = GetKeyboardIndices;
if ExpInfo.deviceIndices(3) == 0
    KbDeviceIndex = KbDevices(1); % internal
else
    KbDeviceIndex = KbDevices(end); % external
end

Screen('Flip',windowPtr);

resp=1; tic
while resp
    [~, ~, keyCode] = KbCheck(KbDeviceIndex);
    if keyCode(KbName('ESCAPE'))
        ShowCursor;
        sca;
        error('Escape');
    elseif keyCode(KbName('LeftArrow'))
        order = 1;
        resp=0;
    elseif keyCode(KbName('DownArrow'))
        order = 2;
        resp=0;
    elseif keyCode(KbName('RightArrow'))
        order = 3;
        resp=0;
    end
end
RT            = toc;
Screen('Flip',windowPtr); WaitSecs(jitter(ExpInfo.ITI(1), ExpInfo.ITI(2))); %ITI

end

