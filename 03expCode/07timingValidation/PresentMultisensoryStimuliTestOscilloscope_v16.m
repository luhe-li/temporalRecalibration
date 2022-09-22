function PresentMultisensoryStimuliTestOscilloscope_v16(TrialNum,...
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

%----------------------------------------------------------------------
%---------------------display audiovisual stimuli----------------------
%----------------------------------------------------------------------

% allow top priority for better temporal accuracy
Priority(ScreenInfo.topPriorityLevel);

if trialSOA > 0 %present V first
    % stimulus onset
    vbl1 = Screen('Flip', windowPtr);

    % Play A after SOA elapsed
%     PsychPortAudio('Start', pahandle, 1, vbl1 + (SOA_frame) * IFI);
    PsychPortAudio('Start', pahandle, 1, vbl1 + (SOA_frame + 1) * IFI);

    % play V immediately
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    vbl2 = Screen('Flip', windowPtr, vbl1 + 0.5 * IFI);
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 

elseif trialSOA < 0 %present A first
    % stimulus onset
    vbl1 = Screen('Flip', windowPtr);

    % play A immediately
    PsychPortAudio('Start', pahandle, 1, vbl1 + IFI);

    % play V after SOA elapsed
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    vbl2 = Screen('Flip', windowPtr, vbl1 + (SOA_frame + 0.5) * IFI);
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 

else %present A and V at the same time
    % pre-onset timestamp
    vbl1 = Screen('Flip', windowPtr);

    % stimulus onset after two frames
%     tWhen = vbl1 + 1.5 * IFI;
%     tPredictedVisualOnset = PredictVisualOnsetForTime(windowPtr, tWhen);

    % play A
    PsychPortAudio('Start', pahandle, 1, vbl1 + IFI);

    % play V
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    vbl2 = Screen('Flip', windowPtr, vbl1 + .5 * IFI);
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 

end

Priority(0);

% blank screen 2
Screen('Flip',windowPtr); WaitSecs(0.3); %ITI

end