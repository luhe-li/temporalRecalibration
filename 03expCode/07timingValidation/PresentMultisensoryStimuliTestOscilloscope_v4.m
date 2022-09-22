function checkt = PresentMultisensoryStimuliTestOscilloscope_v4(TrialNum,...
    ExpInfo,ScreenInfo,VSinfo, AudInfo,pahandle,windowPtr, KbInfo)

% This function uses waitsecs to control time;
% It uses a white screen instead of blob when presenting visual stimulus
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

% calculate time (time is always positive)
SOA_time = abs(trialSOA);
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
    Screen('Flip',windowPtr, t + (stim_frame - 0.5) * IFI); 

    % offset-to-onset interval
    WaitSecs(otoi_frame * IFI)

    % play A
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    WaitSecs(AudInfo.stimDura);
    PsychPortAudio('Stop', pahandle);
    
    checkt = 0;
    
elseif trialSOA < 0 %present A first
   % buffer
     Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    
    % play A
%     t0 = GetSecs();
    tic;
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    WaitSecs(AudInfo.stimDura);
    PsychPortAudio('Stop', pahandle);

    a = toc;
    while a < SOA_time
        a = toc;
    end

% play V
    t1 = Screen('Flip', windowPtr);
%     t1 = Screen('Flip',windowPtr, t0 +(trialSOA_frame - 0.5) * IFI);
    Screen('Flip',windowPtr, t1 + (stim_frame - 0.5) * IFI); 
    checkt = a;
    
else %present A and V at the same time
   % buffer
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);

    % play V&A
    t = Screen('Flip',windowPtr);
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    Screen('Flip',windowPtr, t + (stim_frame - 0.5) * IFI); 
    PsychPortAudio('Stop', pahandle);

    checkt = 0;
end


% blank screen 2
Screen('Flip',windowPtr); WaitSecs(1-abs(trialSOA)); %ITI

end