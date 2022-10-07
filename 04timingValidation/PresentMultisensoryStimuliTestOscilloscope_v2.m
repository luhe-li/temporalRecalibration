function PresentMultisensoryStimuliTestOscilloscope_v2(TrialNum,...
    ExpInfo,ScreenInfo,VSinfo, AudInfo,pahandle,windowPtr, KbInfo)

% This function flips frames to control time
%----------------------------------------------------------------------
%-----------Calculate the coordinates of the target stimuli------------
%----------------------------------------------------------------------
%Make visual stimuli
blob_coordinates = [ScreenInfo.xmid, ScreenInfo.liftingYaxis];
dotCloud = generateOneBlob(windowPtr,blob_coordinates,VSinfo,ScreenInfo);

%----------------------------------------------------------------------
%---------------------display audiovisual stimuli----------------------
%----------------------------------------------------------------------

% calculate frames
trialSOA = ExpInfo.trialSOA(TrialNum);
trialSOA_frame = abs(trialSOA) / ScreenInfo.frameDura;
stim_frame = ExpInfo.stimFrame;
otoi_frame = trialSOA_frame - stim_frame; % % offset-to-onset interval frame

% calculate time
stim_time = stim_frame*ScreenInfo.frameDura;
IFI = ScreenInfo.frameDura;
otoi_time = abs(trialSOA) - stim_time;

% present stimulus
if trialSOA > 0 %present V first
    % buffer
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);

    % play V
    t = Screen('Flip',windowPtr);
    Screen('Flip',windowPtr, t + stim_time - IFI/2); 

    % offset-to-onset interval
    WaitSecs(otoi_time)

    % play A
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    WaitSecs(AudInfo.stimDura);
    PsychPortAudio('Stop', pahandle);
    
elseif trialSOA < 0 %present A first
   % buffer
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);

    % play A
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    WaitSecs(AudInfo.stimDura);
    PsychPortAudio('Stop', pahandle);

    % offset-to-onset interval
    WaitSecs(otoi_time)

    % play V
    t = Screen('Flip',windowPtr);
    Screen('Flip',windowPtr, t + stim_time - IFI/2); 
    
else %present A and V at the same time
   % buffer
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);

    % play V&A
    t = Screen('Flip',windowPtr);
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    Screen('Flip',windowPtr, t + stim_time - IFI/2); 
    PsychPortAudio('Stop', pahandle);

end

% blank screen 2
Screen('Flip',windowPtr); WaitSecs(1-abs(trialSOA)); %ITI

end