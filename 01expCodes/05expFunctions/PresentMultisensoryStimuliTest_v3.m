function [order, RT] = PresentMultisensoryStimuliTest_v3(TrialNum,...
    ExpInfo,ScreenInfo,VSinfo, AudInfo,pahandle,windowPtr)

% This function uses 'when' to control stimulus onsets
% PsychPortAudio is always called first
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

%----------------------------------------------------------------------
%--------------Record response ---------------
%----------------------------------------------------------------------
% blank screen 2
PsychPortAudio('Stop', pahandle, 0);
Screen('Flip',windowPtr); 
WaitSecs(ExpInfo.jBlankDuration2(TrialNum));

% response cue
DrawFormattedText(windowPtr, 'V first: 1\n Simultaneous: 2\n A first: 3',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);

resp=1; tic;
while resp
    [~, ~, keyCode] = KbCheck();
    if keyCode(KbName('numLock'))
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
RT            = toc;
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jITI(TrialNum)); %ITI

end