function [detection, RT] = PresentDetectionStimulus(TrialNum,...
    ExpInfo,ScreenInfo,VSinfo, AudInfo,pahandle,windowPtr, KbInfo)

% This function flips frames to control time

%----------------------------------------------------------------------
%-----------Pre-stimulus screen------------
%----------------------------------------------------------------------
% fixation
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x1_lb,...
    ScreenInfo.y1_lb, ScreenInfo.x1_ub, ScreenInfo.y1_ub]);
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x2_lb,...
    ScreenInfo.y2_lb, ScreenInfo.x2_ub, ScreenInfo.y2_ub]);
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jFixation(TrialNum));

% blank screen 1
Screen('Flip',windowPtr);
Screen('DrawTexture',windowPtr,VSinfo.blk_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
WaitSecs(ExpInfo.jBlankDuration1(TrialNum));

%----------------------------------------------------------------------
%-----------Display audio/visual stimulus------------
%----------------------------------------------------------------------

switch ExpInfo.session
    case 1 % auditory
        % start play
        PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
        PsychPortAudio('Start', pahandle, 1, 0, 0);
        WaitSecs(AudInfo.stimDura);
        PsychPortAudio('Stop', pahandle);
        verb = 'hear';
        
    case 2 % visual
        % make visual stimuli
        blob_coordinates = [ScreenInfo.xmid, ScreenInfo.liftingYaxis];
        dotCloud = generateOneBlob(windowPtr,blob_coordinates,VSinfo,ScreenInfo);
        for j = 1:ExpInfo.stimFrame
            Screen('DrawTexture',windowPtr,dotCloud,[],...
                [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
            Screen('Flip',windowPtr);
        end
        verb = 'see';
end

%----------------------------------------------------------------------
%--------------Record response ---------------
%----------------------------------------------------------------------
% blank screen 2
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jBlankDuration2(TrialNum));

% response cue
DrawFormattedText(windowPtr, ['Click left if you ' verb ' it \n Click right if you did not ' verb ' it '],...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);
click = 0; tic
while click == 0
    %click left button: not oddball; click right button: yes it's an
    %oddball
    [click,~,~,detection] = GetClicks;
end
RT            = toc;

Screen('Flip',windowPtr); WaitSecs(ExpInfo.jITI(TrialNum)); %ITI

end