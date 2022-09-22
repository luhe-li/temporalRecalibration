function PresentMultisensoryStimuliTestOscilloscope_v15(TrialNum,...
    ExpInfo,ScreenInfo,VSinfo, AudInfo,pahandle,windowPtr, KbInfo)

% This function uses 'when' to control stimulus onsets
% flip is always called first
%----------------------------------------------------------------------
%-----------Create stimuli------------
%----------------------------------------------------------------------
%Make visual stimuli
blob_coordinates = [ScreenInfo.xmid, ScreenInfo.liftingYaxis];
dotCloud = generateOneBlob(windowPtr,blob_coordinates,VSinfo,ScreenInfo);

% calculate frames
trialSOA = ExpInfo.trialSOA(TrialNum);
SOA_frame = abs(trialSOA) / ScreenInfo.ifi;
stim_frame = ExpInfo.stimFrame;

% calculate time (time is always positive)
SOA_time = abs(trialSOA);
stim_time = stim_frame*ScreenInfo.ifi;
IFI = ScreenInfo.ifi;
wait_frame = 2;

% create buffer
PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);

%----------------------------------------------------------------------
%---------------------display audiovisual stimuli----------------------
%----------------------------------------------------------------------

% allow top priority for better temporal accuracy
Priority(ScreenInfo.topPriorityLevel);

if trialSOA > 0 %present V first

    % play V immediately
    vbl2 = Screen('Flip', windowPtr);
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 

    tic;
    a = toc;
    while a < SOA_time - stim_time
        a = toc;
    end

    % Play A after SOA elapsed
    PsychPortAudio('Start', pahandle, 1, 0, 0);
    
elseif trialSOA < 0 %present A first

    % play A immediately
    PsychPortAudio('Start', pahandle, 1, 0);

    tic;
    a = toc;
    while a < SOA_time - stim_time
        a = toc;
    end
    

    % play V after SOA elapsed
    vbl2 = Screen('Flip', windowPtr);
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 

else %present A and V at the same time
    % pre-onset timestamp

    % play V
    vbl2 = Screen('Flip', windowPtr);
    PsychPortAudio('Start', pahandle, 1, 0);
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 
end

Priority(0);

% blank screen 2
Screen('Flip',windowPtr); WaitSecs(1-abs(trialSOA)); %ITI

end