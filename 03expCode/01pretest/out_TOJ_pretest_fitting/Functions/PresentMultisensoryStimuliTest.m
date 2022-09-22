function [order, RT] = PresentMultisensoryStimuliTest(TrialNum,...
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

% fixation
Screen('FillRect', windowPtr,[255 0 0], [ScreenInfo.x1_lb,...
    ScreenInfo.y1_lb, ScreenInfo.x1_ub, ScreenInfo.y1_ub]);
Screen('FillRect', windowPtr,[255 0 0], [ScreenInfo.x2_lb,...
    ScreenInfo.y2_lb, ScreenInfo.x2_ub, ScreenInfo.y2_ub]);
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jFixation(TrialNum));

% blank screen 1
Screen('Flip',windowPtr);
Screen('DrawTexture',windowPtr,VSinfo.blk_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
WaitSecs(ExpInfo.jBlankDuration1(TrialNum));

% present test stimuli
trialSOA = ExpInfo.trialSOA(TrialNum);
trialSOA_frame = abs(trialSOA) / ScreenInfo.frameDura;
stim_frame = ExpInfo.stimFrame;
otoi_frame = trialSOA_frame - stim_frame; % offset-to-onset interval frame
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
    PsychPortAudio('Start', pahandle, 1, 0, 0); % A onset
    for j = 1:stim_frame
        Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
        Screen('Flip',windowPtr); 
    end
    PsychPortAudio('Stop', pahandle); % A offset
end

%----------------------------------------------------------------------
%--------------Record response ---------------
%----------------------------------------------------------------------
% blank screen 2
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jBlankDuration2(TrialNum));

% response cue
DrawFormattedText(windowPtr, 'V first: 1\n Simultaneous: 2\n A first: 3',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);

resp=1; tic
while resp
    [~, ~, keyCode] = KbCheck(KbInfo.KbDeviceIndex);
    if keyCode(KbName('ESCAPE'))
        ShowCursor;
        sca;
        error('Escape');
    elseif keyCode(KbName('1'))
        order = 1;
        resp=0;
    elseif keyCode(KbName('2'))
        order = 2;
        resp=0;
    elseif keyCode(KbName('3'))
        order = 3;
        resp=0;
    end
end
RT = toc;
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jITI(TrialNum)); %ITI

end