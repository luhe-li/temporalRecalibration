function [resp, RT] = PresentDiscriminationStimulus_v3(TrialNum,...
    ExpInfo,ScreenInfo,VSinfo, AudInfo,pahandle,windowPtr)

% This function uses 'when' to control stimulus onsets
% PA is always called first

%----------------------------------------------------------------------
%-----------Create stimuli------------
%----------------------------------------------------------------------
% calculate frames
stim_frame = ExpInfo.stimFrame;

% calculate time (time is always positive)
IFI = ScreenInfo.ifi;

%----------------------------------------------------------------------
%-----------Pre-stimulus screen------------
%----------------------------------------------------------------------

% allow top priority for better temporal accuracy
Priority(ScreenInfo.topPriorityLevel);

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

switch ExpInfo.session
    case 1 % auditory
        if ExpInfo.counterbalance(TrialNum) == 1 % standard first
            %%%%%% first stimulus %%%%%
            % adjust intensity
            PsychPortAudio('Volume', pahandle, AudInfo.standard);
            % start play
            PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
            PsychPortAudio('Start', pahandle, 1, 0);
            
            %%%%% ISI %%%%%
            for j = 1:floor(ExpInfo.jISI(TrialNum)/IFI)
                Screen('DrawTexture',windowPtr,VSinfo.blk_texture,[],...
                    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
                Screen('Flip',windowPtr);
            end
            
            %%%%% second stimulus %%%%%
            % adjust intensity
            PsychPortAudio('Volume', pahandle, AudInfo.shuffledIntensity(TrialNum));
            % start play
            PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
            PsychPortAudio('Start', pahandle, 1, 0);
            
        elseif ExpInfo.counterbalance(TrialNum) == 2 % comparison first
            %%%%%% first stimulus %%%%%
            % adjust intensity
            PsychPortAudio('Volume', pahandle, AudInfo.shuffledIntensity(TrialNum));
            % start play
            PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
            PsychPortAudio('Start', pahandle, 1, 0);
            
            %%%%% ISI %%%%%
            for j = 1:floor(ExpInfo.jISI(TrialNum)/IFI)
                Screen('DrawTexture',windowPtr,VSinfo.blk_texture,[],...
                    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
                Screen('Flip',windowPtr);
            end
            
            %%%%% second stimulus %%%%%
            % adjust intensity
            PsychPortAudio('Volume', pahandle, AudInfo.standard);
            % start play
            PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
            PsychPortAudio('Start', pahandle, 1, 0);
        end
        
    case 2 % visual
        if ExpInfo.counterbalance(TrialNum) == 1 % standard first
            %%%%%% first stimulus %%%%%
            % adjust intensity
            VSinfo.intensity   = VSinfo.standard;
            VSinfo.Cloud       = 255.*VSinfo.intensity.*reshape(VSinfo.cloud,length(VSinfo.x),length(VSinfo.y));
            % make visual stimuli
            blob_coordinates = [ScreenInfo.xmid, ScreenInfo.liftingYaxis];
            dotCloud = generateOneBlob(windowPtr,blob_coordinates,VSinfo,ScreenInfo);
            % start play
            Screen('DrawTexture',windowPtr,dotCloud,[],...
                [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
            vbl1 = Screen('Flip',windowPtr);
            Screen('Flip',windowPtr, vbl1 + (stim_frame - 0.5) * IFI);
            
            %%%%% ISI %%%%%
            WaitSecs(ExpInfo.jISI(TrialNum));
            
            %%%%% second stimulus %%%%%
            % adjust intensity
            VSinfo.intensity   = VSinfo.shuffledIntensity(TrialNum);
            VSinfo.Cloud       = 255.*VSinfo.intensity.*reshape(VSinfo.cloud,length(VSinfo.x),length(VSinfo.y));
            % make visual stimuli
            blob_coordinates = [ScreenInfo.xmid, ScreenInfo.liftingYaxis];
            dotCloud = generateOneBlob(windowPtr,blob_coordinates,VSinfo,ScreenInfo);
            % start play
            Screen('DrawTexture',windowPtr,dotCloud,[],...
                [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
            vbl1 = Screen('Flip',windowPtr);
            Screen('Flip',windowPtr, vbl1 + (stim_frame - 0.5) * IFI);
            
        elseif ExpInfo.counterbalance(TrialNum) == 2 % comparison first
            %%%%%% first stimulus %%%%%
            % adjust intensity
            VSinfo.intensity   = VSinfo.shuffledIntensity(TrialNum);
            VSinfo.Cloud       = 255.*VSinfo.intensity.*reshape(VSinfo.cloud,length(VSinfo.x),length(VSinfo.y));
            % make visual stimuli
            blob_coordinates = [ScreenInfo.xmid, ScreenInfo.liftingYaxis];
            dotCloud = generateOneBlob(windowPtr,blob_coordinates,VSinfo,ScreenInfo);
            % start play
            Screen('DrawTexture',windowPtr,dotCloud,[],...
                [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
            vbl1 = Screen('Flip',windowPtr);
            Screen('Flip',windowPtr, vbl1 + (stim_frame - 0.5) * IFI);
            
            %%%%% ISI %%%%%
             WaitSecs(ExpInfo.jISI(TrialNum));
            
            %%%%% second stimulus %%%%%
            % adjust intensity
            VSinfo.intensity   = VSinfo.standard;
            VSinfo.Cloud       = 255.*VSinfo.intensity.*reshape(VSinfo.cloud,length(VSinfo.x),length(VSinfo.y));
            % make visual stimuli
            blob_coordinates = [ScreenInfo.xmid, ScreenInfo.liftingYaxis];
            dotCloud = generateOneBlob(windowPtr,blob_coordinates,VSinfo,ScreenInfo);
            % start play
            Screen('DrawTexture',windowPtr,dotCloud,[],...
                [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
            vbl1 = Screen('Flip',windowPtr);
            Screen('Flip',windowPtr, vbl1 + (stim_frame - 0.5) * IFI);
            
        end
end

%----------------------------------------------------------------------
%--------------Record response ---------------
%----------------------------------------------------------------------
% blank screen 2
Screen('Flip',windowPtr); WaitSecs(ExpInfo.jBlankDuration2(TrialNum));

% response cue
postcue = {'louder', 'brighter'};
DrawFormattedText(windowPtr, ['Which is ' postcue{ExpInfo.session} '?\n L: first\n R: second'],...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);
click = 0; tic
while click == 0
    %click left button: not oddball; click right button: yes it's an
    %oddball
    [click,~,~,resp] = GetClicks;
end
RT = toc;
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.jITI(TrialNum)); %ITI

end