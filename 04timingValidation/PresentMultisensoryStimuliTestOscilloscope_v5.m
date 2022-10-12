function checkt = PresentMultisensoryStimuliTestOscilloscope_v5(TrialNum,...
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
SOA_frame = abs(trialSOA) / ScreenInfo.ifi;
stim_frame = ExpInfo.stimFrame;
otoi_frame = SOA_frame - stim_frame; % % offset-to-onset interval frame

% calculate time (time is always positive)
SOA_time = abs(trialSOA);
stim_time = stim_frame*ScreenInfo.ifi;
IFI = ScreenInfo.ifi;

Priority(ScreenInfo.topPriorityLevel);
% present stimulus
if trialSOA > 0 %present V first
     % buffer
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);

    % play V
    vbl1 = Screen('Flip', win); % vbl1 marks AV stimulus onset
    Screen('Flip',windowPtr, vbl1 + (stim_frame - 0.5) * IFI); 

    % Play A
    PsychPortAudio('Start', pahandle, 0, vbl1 + (SOA_frame - 0.5) * IFI, 1);
    WaitSecs(stim_time); % or use tic toc while
    PsychPortAudio('Stop', pahandle);
    
elseif trialSOA < 0 %present A first
    % buffer A
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);

    % return timestamp of the black onset
    vbl1 = Screen('Flip', win); % vbl1 marks AV stimulus onset

    % play A
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    WaitSecs(stim_time)
    PsychPortAudio('Stop', pahandle);

   % buffer V
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);

    % play V 
    vbl2 = Screen('Flip', windowPtr, vbl1 + (SOA_frame - 0.5) * IFI);
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 
    % or by flip frames
    
else %present A and V at the same time
    % buffer A
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);

    % return timestamp of the black onset
    vbl1 = Screen('Flip', windowPtr);

   % buffer V
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);

    % Compute tWhen onset time for wanted visual onset at >= tWhen:
    tWhen = vbl1 + (stim_frame - 0.5) * IFI;
    tPredictedVisualOnset = PredictVisualOnsetForTime(windowPtr, tWhen);

    % A on
    PsychPortAudio('Start', pahandle, 0, tPredictedVisualOnset, 1);

    % V on
    vbl2 = Screen('Flip', windowPtr, tWhen);

    % V off
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 
    % or by flip frames

    % A off
    PsychPortAudio('Stop', pahandle);

end

Priority(0);

% blank screen 2
Screen('Flip',windowPtr); WaitSecs(1-abs(trialSOA)); %ITI

end