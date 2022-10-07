function PresentMultisensoryStimuliTestOscilloscope(TrialNum,...
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

% % blank screen 1
% Screen('Flip',windowPtr);
% Screen('DrawTexture',windowPtr,VSinfo.blk_texture,[],...
%     [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
% WaitSecs(0.5);

% present test stimuli
trialSOA = ExpInfo.trialSOA(TrialNum);
trialSOA_frame = abs(trialSOA) / ScreenInfo.frameDura;
stim_frame = ExpInfo.stimFrame;
otoi_frame = trialSOA_frame - stim_frame; % % offset-to-onset interval frame
if trialSOA > 0 %present V first
    % visual stimulus
    for j = 1:stim_frame
        Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        Screen('Flip',windowPtr); 
    end
    % offset-to-onset interval
    for j = 1:otoi_frame
        Screen('DrawTexture',windowPtr,VSinfo.blk_texture,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        Screen('Flip',windowPtr);
    end
    %start playing the sound
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    WaitSecs(AudInfo.stimDura);
    PsychPortAudio('Stop', pahandle);
    
elseif trialSOA < 0 %present A first
    %start playing the sound
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    WaitSecs(AudInfo.stimDura);
    PsychPortAudio('Stop', pahandle);
    % offset-to-onset interval
    for j = 1:otoi_frame
        Screen('DrawTexture',windowPtr,VSinfo.blk_texture,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        Screen('Flip',windowPtr);
    end
    %present V second
    for j = 1:stim_frame
        Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        Screen('Flip',windowPtr); 
    end
else %present A and V at the same time
    %start playing the sound while displaying the visual stimulus
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    PsychPortAudio('Start', pahandle, 0, 0, 1); % A onset
    for j = 1:stim_frame
        Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        Screen('Flip',windowPtr); 
    end
    PsychPortAudio('Stop', pahandle); % A offset
end

% blank screen 2
Screen('Flip',windowPtr); WaitSecs(1-trialSOA); %ITI

end