function [oddball, RT] = PresentMultisensoryStimuliExposure(TrialNum,...
    ExpInfo,ScreenInfo,VSinfo, AudInfo,pahandle,windowPtr)

% -------------------------------------------------------------
%-----------Calculate the coordinates of the target stimuli------------
%----------------------------------------------------------------------
%Make visual stimuli
blob_coordinates = [ScreenInfo.xmid, ScreenInfo.liftingYaxis];
dotCloud = generateOneBlob(windowPtr,blob_coordinates,VSinfo,ScreenInfo);

%----------------------------------------------------------------------
%---------------------display audiovisual stimuli----------------------
%----------------------------------------------------------------------

% fixation
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x1_lb,...
    ScreenInfo.y1_lb, ScreenInfo.x1_ub, ScreenInfo.y1_ub]);
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x2_lb,...
    ScreenInfo.y2_lb, ScreenInfo.x2_ub, ScreenInfo.y2_ub]);
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jFixation(TrialNum));

% blank screen 1
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jBlankDuration1(TrialNum));

%given the SOA, calculate how many frames we present the blank screen
% SOA_frames = abs(ExpInfo.trialSOA(TrialNum))/(1000/ScreenInfo.frameRate); %1000ms/60Hz
trialSOA = ExpInfo.adaptor; %s
otoi = abs(trialSOA) - VSinfo.duration; %s
if trialSOA > 0 %present V first
    % visual stimulus
    Screen('DrawTexture',windowPtr,dotCloud,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    onset = Screen('Flip',windowPtr); % onset
    offset = Screen('Flip',windowPtr, onset + VSinfo.duration - ScreenInfo.ifi/2); % offset
    
    % ISI = SOA - stimDura
    %     Screen('DrawTexture',windowPtr,VSinfo.blk_texture,[],...
    %             [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    %     Screen('Flip',windowPtr, offset + ISI  - ScreenInfo.ifi/2);
    WaitSecs(otoi - ScreenInfo.ifi/2);
    
    %start playing the sound
    %     PsychPortAudio('FillBuffer', pahandle, AudInfo.GaussianWhiteNoise);
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    PsychPortAudio('Start', pahandle, 1, 0, 0);
    WaitSecs(AudInfo.stimDura);
    PsychPortAudio('Stop', pahandle);
    
elseif trialSOA < 0 %present A first
    %start playing the sound
    %     PsychPortAudio('FillBuffer', pahandle, AudInfo.GaussianWhiteNoise);
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    
    PsychPortAudio('Start', pahandle, 1, 0, 0);
    WaitSecs(AudInfo.stimDura);
    PsychPortAudio('Stop', pahandle);
    % ISI
    WaitSecs(otoi - ScreenInfo.ifi/2);
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
%-------------- Post-cue ---------------
%----------------------------------------------------------------------
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jBlankDuration2(TrialNum));
DrawFormattedText(windowPtr, ExpInfo.postCue,...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);

%----------------------------------------------------------------------
%--------------Record response ---------------
%----------------------------------------------------------------------
click = 0; tic
while click == 0
    %click left button: not oddball; click right button: yes it's an
    %oddball
    [click,~,~,oddball] = GetClicks;
end
RT     = toc;
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jITI(TrialNum)); %ITI

end

