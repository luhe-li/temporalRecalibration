function [oddball, RT] = PresentMultisensoryStimuliExposure_v3(TrialNum,...
    ExpInfo,ScreenInfo,VSinfo, AudInfo,pahandle,windowPtr, practice)

% This function uses 'when' to control stimulus onsets
% PA is always called first

%----------------------------------------------------------------------
%-----------Create stimuli------------
%----------------------------------------------------------------------
%Make visual stimuli
blob_coordinates = [ScreenInfo.xmid, ScreenInfo.liftingYaxis];
dotCloud = generateOneBlob(windowPtr,blob_coordinates,VSinfo,ScreenInfo);

% calculate frames
trialSOA = ExpInfo.adaptor;
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
% Screen('FillOval', windowPtr,[255 255 255]*0.2, ScreenInfo.fixDot);
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x1_lb,...
    ScreenInfo.y1_lb, ScreenInfo.x1_ub, ScreenInfo.y1_ub]);
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x2_lb,...
    ScreenInfo.y2_lb, ScreenInfo.x2_ub, ScreenInfo.y2_ub]);
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jExpoFixation(TrialNum));

% blank screen 1
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jExpoBlankDuration1(TrialNum));

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

% % % response cue
% DrawFormattedText(windowPtr, 'A: 1 \nV&A: 2 \nV: 3',...
%     'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
% Screen('Flip',windowPtr);

% record response within the response period
start=GetSecs();
t=GetSecs();
resp=1;
tic;
while t-start <=  ExpInfo.jExpoRespPeriod(TrialNum)
    t=GetSecs();
    if resp
        [~, ~, keyCode] = KbCheck();
        if keyCode(KbName('ESCAPE'))
            ShowCursor;
            sca;
            error('Escape');
        elseif keyCode(KbName('1'))
            RT = toc;
            oddball = 1; % visual oddball
            resp=0;
        elseif keyCode(KbName('2'))
            RT = toc;
            oddball = 2; % both
            resp=0;
        elseif keyCode(KbName('3')) 
            RT = toc;
            oddball = 3; % auditory
            resp=0;
        end
    end
end

% when no response is made
if ~exist('oddball', 'var')
   oddball = NaN;
   RT = NaN; 
end

% if it is practice trial, feedbacks about oddball are provided
if practice
    DrawFormattedText(windowPtr, ExpInfo.feedbackType{ExpInfo.feedback(TrialNum)},...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
    Screen('Flip',windowPtr);
    WaitSecs(1.0);
end

Screen('Flip',windowPtr);

end
