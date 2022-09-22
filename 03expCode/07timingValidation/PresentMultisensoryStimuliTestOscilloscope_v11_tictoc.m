function [t1, t2, t3, t4, t5, t6, t7, t8] = PresentMultisensoryStimuliTestOscilloscope_v11_tictoc(TrialNum,...
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
wait_frame = 3;

Priority(ScreenInfo.topPriorityLevel);
% present stimulus
if trialSOA > 0 %present V first
     % buffer
     tic;
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    t1 = toc;
    tic;
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    t2 = toc;

    % play V
    tic;
    vbl1 = Screen('Flip', windowPtr); % vbl1 marks AV stimulus onset
    t3 = toc;
    tic;
    Screen('Flip',windowPtr, vbl1 + (stim_frame - 0.5) * IFI); 
    t4 =toc;

    % Play A
    tic;
    PsychPortAudio('Start', pahandle, 0, vbl1 + (SOA_frame - 0.5) * IFI, 1);
    t5 = toc;
    tic;
    WaitSecs(stim_time); % or use tic toc while
    t6 = toc;
    tic;
    PsychPortAudio('Stop', pahandle);
    t7 = toc;
    t8 = 0;
    
elseif trialSOA < 0 %present A first
    % buffer A
    tic;
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    t1 = toc;

    % return timestamp of the black onset
    tic;
    vbl1 = Screen('Flip', windowPtr); % vbl1 marks AV stimulus onset
    t2 = toc;

    % play A
    tic;
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    t3 = toc;
    tic;
    WaitSecs(stim_time);
    t4 = toc;
    tic;
    PsychPortAudio('Stop', pahandle);
    t5 = toc;

   % buffer V
   tic;
    Screen('DrawTexture',windowPtr,dotCloud,[],...
            [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    t6 = toc;

    % play V 
    tic;
    vbl2 = Screen('Flip', windowPtr, vbl1 + (SOA_frame - 0.5) * IFI);
    t7 = toc;
    tic;
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 
    t8 = toc;
    
else %present A and V at the same time
    % return timestamp of the black onset
    tic;
    vbl1 = Screen('Flip', windowPtr);
t1 = toc;   

%     % buffer A
tic;
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    t2 = toc;
    tic;
    Screen('DrawTexture',windowPtr,dotCloud,[],...
                [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
    t3 = toc;
    
    % A on
    tic;
    PsychPortAudio('Start', pahandle, 0, vbl1 + (wait_frame - 0.5) * IFI, 1);
    t4 = toc;

    % V on
    tic;
    vbl2 = Screen('Flip', windowPtr, vbl1 + (wait_frame - 0.5) * IFI);
    t5 = toc;

    % V off
    tic;
    Screen('Flip',windowPtr, vbl2 + (stim_frame - 0.5) * IFI); 
    t6 = toc;

    % A off
    tic;
    PsychPortAudio('Stop', pahandle);
    t7 = toc;
    t8 = 0;

end

Priority(0);

% blank screen 2
Screen('Flip',windowPtr); WaitSecs(1-abs(trialSOA)); %ITI

end