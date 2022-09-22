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
bias = ExpInfo.bias;

% create buffer
PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);

%----------------------------------------------------------------------
%---------------------display audiovisual stimuli----------------------
%----------------------------------------------------------------------

% allow top priority for better temporal accuracy
Priority(ScreenInfo.topPriorityLevel);

if trialSOA > 0 %present V first
    tic;
    % stimulus onset
    vbl1 = Screen('Flip', windowPtr);

    % Play A after SOA + 1 frame elapsed
    PsychPortAudio('Start', pahandle, 1, vbl1 + (SOA_frame + 1) * IFI + bias);

    % play V immediately after 1 frame
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    vbl2 = Screen('Flip', windowPtr, vbl1 + 0.5 * IFI);
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 
    
    % add waiting to make sure duration is around SOA+1 frame
    a = toc;
    while a < SOA_time + IFI
       a = toc; 
    end
    
elseif trialSOA < 0 %present A first
    % stimulus onset
    vbl1 = Screen('Flip', windowPtr);

    % play A immediately after 1 frame
    PsychPortAudio('Start', pahandle, 1, vbl1 + IFI + bias);

    % play V after SOA + 1 frame elapsed
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    vbl2 = Screen('Flip', windowPtr, vbl1 + (SOA_frame + 0.5) * IFI);
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 

else %present A and V at the same time
    % pre-onset timestamp
    vbl1 = Screen('Flip', windowPtr);

    % play A after 1 frame
    PsychPortAudio('Start', pahandle, 1, vbl1 + IFI + bias);

    % play V after 1 frame
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    vbl2 = Screen('Flip', windowPtr, vbl1 + 0.5 * IFI);
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 

end

Priority(0);

% blank screen 2
WaitSecs(0.3); %ITI

end