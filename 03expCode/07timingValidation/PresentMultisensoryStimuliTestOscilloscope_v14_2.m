function [t1, t2, t3, t4, t5] = PresentMultisensoryStimuliTestOscilloscope_v14_2(TrialNum,...
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
    % stimulus onset
    tic
    vbl1 = GetSecs();
    t1 = toc;
    
    % Play A after SOA elapsed
    tic;
    PsychPortAudio('Start', pahandle, 1, vbl1 + (SOA_frame) * IFI, 0);
   t2 = toc;
    % play V immediately
    tic;
    vbl2 = Screen('Flip', windowPtr, vbl1 + 0.5 * IFI);
    t3 = toc;
    tic;
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 
    t4=  toc;
    t5 = 0;
elseif trialSOA < 0 %present A first
    % stimulus onset
    tic;
    vbl1 = GetSecs();
    t1 = toc;
    % play A immediately
    tic;
    PsychPortAudio('Start', pahandle, 1, 0, 0);
    t2 = toc;
    % play V after SOA elapsed
    tic;
    vbl2 = Screen('Flip', windowPtr, vbl1 + (SOA_frame - 0.5) * IFI);
    t3 = toc;
    tic;
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 
    t4 = toc;
    t5 = 0
else %present A and V at the same time
    % pre-onset timestamp
    tic;
    vbl1 = GetSecs();
    t1 = toc;
    % stimulus onset after two frames
    tic;
    tWhen = vbl1 + 1.5 * IFI;
    tPredictedVisualOnset = PredictVisualOnsetForTime(windowPtr, tWhen);
    t2 = toc;
    % play A
    tic;
    PsychPortAudio('Start', pahandle, 1, tPredictedVisualOnset, 0);
    t3 = toc;
    % play V
    tic;
    vbl2 = Screen('Flip', windowPtr, tWhen);
    t4 = toc;
    tic;
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 
    t5 = toc;
end

Priority(0);

% blank screen 2
Screen('Flip',windowPtr); WaitSecs(1-abs(trialSOA)); %ITI

end