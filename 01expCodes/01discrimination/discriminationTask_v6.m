%% discrimination task: TWO SIDED

% 2022/05/v6

% Equipment indices are set to work on Chartest with the projector, central
% speaker and number pad. Psychportaudio was calibrated by 32.53 ms.

% In v6, both auditory and visual intensity levels are sampled and fitted
% on a LINEAR SCALE.

% The discrimination task obtains contrast and loudness JND to set pamaters
% for exposure task. Enter session number to choose an auditory (1) or
% visual (2) task. In both tasks, constant stimuli methods were used, so
% that there wil be 7 intensity levels x 20 trials for the comparison
% stimulus. The standard stimulus has a constant intensity at the middle
% level. After stimulus presentation, participants respond which one is
% brighter or louder by mouse cliking.

% Response guide:
% 1: the first is louder/brighter
% 2: the second is louder/brighter

try
    %% enter subject's name
    clear all; close all; clc; rng('shuffle');
    
    %enter subject's name
    ExpInfo.subjID                       = [];
    while isempty(ExpInfo.subjID)        == 1
        try ExpInfo.subjID                   = input('Please enter participant ID#: ') ;
            ExpInfo.session                      = input('Please enter session, 1 = A, 2 = V#: ');
        catch
        end
    end
    
    % choose device indices
    audDevices                           = PsychPortAudio('GetDevices');
%     ExpInfo.deviceIndices = [0, 1]; % use internal display, speaker;
    ExpInfo.deviceIndices                = [0, 4]; % use external display, speaker;
    ExpInfo.sittingDistance              = 10; %in cm, the distance between the screen and participants
    sessioname                           = {'A','V'};
    out1FileName                         = ['discrimination_v6_sub', num2str(ExpInfo.subjID) '_' sessioname{ExpInfo.session} ]; %create file name
    
    % avoid rewriting data
    if exist([out1FileName '.mat'],'file')
        resp                                 = input('To replace the existing file, press y', 's');
        if ~strcmp(resp,'y')
            disp('Experiment stopped.')
            return
        end
    end

   %% define the experiment information
   % define bias
    ExpInfo.bias                         = 0.03253;
    
    % parameters of constant stimuli
    ExpInfo.intLevel                     = 7;
    ExpInfo.numTrials                    = 20;
    
    % intensity of standard stimulus
    AudInfo.standard                     = 0.5; % a ratio between 0 to 1
    VSinfo.standard                      = 0.2; % a ratio between 0 to 1 to be multipled by 255
    pblack                               = 1/8; % set contrast to 255*1/8 for the "black" background

    % intenisty range of comparison stimulus is half below and half above
    ExpInfo.VolumeRange = linspace(AudInfo.standard*0.5, AudInfo.standard*1.5, ExpInfo.intLevel);
%     ExpInfo.ContrastRange                = linspace(VSinfo.standard*0.5, VSinfo.standard*1.5, ExpInfo.intLevel);
    ExpInfo.ContrastRange                = linspace(0.01, 0.39, ExpInfo.intLevel);
    
    
    % numTotalTrials: total number of trials
    ExpInfo.numTotalTrials               = ExpInfo.intLevel * ExpInfo.numTrials;

    % shuffle order
    ExpInfo.trialOrder                   = shuffleit(repmat(1:ExpInfo.intLevel, [1,ExpInfo.numTrials]));
    AudInfo.shuffledIntensity            = ExpInfo.VolumeRange(ExpInfo.trialOrder);
    VSinfo.shuffledIntensity             = ExpInfo.ContrastRange(ExpInfo.trialOrder);
    
    % define counterbalance array
    ExpInfo.counterbalance               = ones(1,ExpInfo.numTotalTrials); % 1 = standard first
    idxComp                              = randsample(ExpInfo.numTotalTrials, ExpInfo.numTotalTrials/2); % randomly choose half of trials to present comparison stimulus first
    ExpInfo.counterbalance(idxComp)      = 2; % 2 = comparison first

    % nTBlock: number of trials in each block
    ExpInfo.numTrialsinBlock             = 70;
    
    % sound and blob presented the same duration
    ExpInfo.stimFrame                    = 2; % frame
    
    % ITI duration
    ExpInfo.fixation                     = [0.1, 0.2];
    ExpInfo.blankDuration1               = [0.4, 0.6]; % between fixation and the first stimulus
    ExpInfo.ISI                          = [0.6, 0.8];
    ExpInfo.blankDuration2               = [0.4, 0.6]; % between the second stimulus and post-cue
    ExpInfo.feedback                     = [0.4, 0.6];
    ExpInfo.ITI                          = [0.2, 0.4];
    
    %initialize a structure that stores all the responses and response time
    [Response.more, Response.RT]         = ...
        deal(NaN(1,ExpInfo.numTotalTrials ));
    
    %% screen setup
    PsychDefaultSetup(2);
    AssertOpenGL();
    GetSecs();
    WaitSecs(0.1);
    KbCheck();
    ListenChar(2); % silence the keyboard
    
    %Canvas size = 53.5" x 40"= 135.8cm x 101.6cm Screen size by the project =
    %1024 pixels x 768 pixels so each centimeter has 7.54 pixels horizontally
    %and 7.56 vertically
    Screen('Preference', 'VisualDebugLevel', 1);
    Screen('Preference', 'SkipSyncTests', 0);
%     Screen('Preference', 'SkipSyncTests', 1); % for testing
    [windowPtr,rect]                     = Screen('OpenWindow', ExpInfo.deviceIndices(1), [255, 255, 255] * pblack); % 1 = external display
%     [windowPtr,rect] = Screen('OpenWindow', ExpInfo.deviceIndices(1), [255, 255, 255] * pblack, [100 100 1000 600]); % for testing
    HideCursor(windowPtr);
    ScreenInfo.topPriorityLevel          = MaxPriority(windowPtr);
    [ScreenInfo.xaxis, ScreenInfo.yaxis] = Screen('WindowSize',windowPtr);
    Screen('TextSize', windowPtr, 25) ;
    Screen('TextFont',windowPtr,'Times');
    Screen('TextStyle',windowPtr,1);
    ScreenInfo.frameRate                 = Screen('FrameRate',0);
    ScreenInfo.ifi                       = Screen('GetFlipInterval', windowPtr);
    [center(1), center(2)]               = RectCenter(rect);
    ScreenInfo.xmid                      = center(1); % horizontal center
    ScreenInfo.ymid                      = center(2); % vertical center
    ScreenInfo.backgroundColor           = 0;
    ScreenInfo.numPixels_perCM           = 7.5;
    ScreenInfo.liftingYaxis              = 300;
    ScreenInfo.y_box_unity               = [-10, 22];
    
    %fixation cross
    ScreenInfo.x1_lb                     = ScreenInfo.xmid-7; ScreenInfo.x2_lb = ScreenInfo.xmid-1;
    ScreenInfo.x1_ub                     = ScreenInfo.xmid+7; ScreenInfo.x2_ub = ScreenInfo.xmid+1;
    ScreenInfo.y1_lb                     = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-1;
    ScreenInfo.y1_ub                     = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+1;
    ScreenInfo.y2_lb                     = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-7;
    ScreenInfo.y2_ub                     = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+7;
  
    %% open loudspeakers and create sound stimuli
    addpath(genpath(PsychtoolboxRoot))
    PsychDefaultSetup(2);
    % get correct sound card
    our_device                           = audDevices(ExpInfo.deviceIndices(2)).DeviceIndex;

    %sampling frequencies
    AudInfo.fs                           = audDevices(ExpInfo.deviceIndices(2)).DefaultSampleRate;
    audioSamples                         = linspace(1,AudInfo.fs,AudInfo.fs);
    
    % make a beep
    AudInfo.stimDura                     = ExpInfo.stimFrame * ScreenInfo.ifi; %s, the duration of auditory stimulus
    AudInfo.tf                           = 500;
    AudInfo.beepLengthSecs               = AudInfo.stimDura;
    beep                                 = MakeBeep(AudInfo.tf, AudInfo.beepLengthSecs, AudInfo.fs);
    AudInfo.Beep                         = [beep; zeros(size(beep))];
    pahandle                             = PsychPortAudio('Open', our_device, [], [], [], 2);%open device
    PsychPortAudio('Volume', pahandle, AudInfo.standard);
    
    % Initialize driver, request low-latency preinit:
    InitializePsychSound(1);
    
    % Perform one warmup trial, to get the sound hardware fully up and running,
    % performing whatever lazy initialization only happens at real first use.
    % This "useless" warmup will allow for lower latency for start of playback
    % during actual use of the audio driver in the real trials:
    PsychPortAudio('FillBuffer', pahandle, AudInfo.Beep);
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    PsychPortAudio('Stop', pahandle, 1);

    %% define the visual stimuli
    VSinfo.duration                      = ExpInfo.stimFrame * ScreenInfo.ifi;%s
    VSinfo.width                         = 401; %(pixel) Increasing this value will make the cloud more blurry
    VSinfo.boxSize                       = 201; %This is the box size for each cloud.
    %set the parameters for the visual stimuli
    VSinfo.blackScreen                   = 255 * pblack *ones(ScreenInfo.xaxis,ScreenInfo.yaxis);
    VSinfo.blankScreen                   = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);
    x                                    = 1:1:VSinfo.boxSize; y = x;
    VSinfo.x                             = x; VSinfo.y = y;
    [X,Y]                                = meshgrid(x,y);
    VSinfo.cloud                         = mvnpdf([X(:) Y(:)],[median(x) median(y)],...
    [VSinfo.width 0; 0 VSinfo.width]);
    pscale                               = (1-pblack)/max(VSinfo.cloud); % the max contrast of the blob adds the background contrast should <= 1
    VSinfo.cloud                         = VSinfo.cloud .* pscale;
    VSinfo.Cloud                         = 255.*VSinfo.standard.*reshape(VSinfo.cloud,length(x),length(y));
    VSinfo.blk_texture                   = Screen('MakeTexture', windowPtr, VSinfo.blackScreen,[],[],[],2);
    
    %% Run the experiment by calling the function InterleavedStaircase
    %  record tart time
    c                                    = clock;
    start                                = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));
    timestamp{1,1}                       = start;
    
    %start the experiment
    DrawFormattedText(windowPtr, 'Press any button to start the discrimination task.',...
        'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
    Screen('Flip',windowPtr);
    KbWait(); WaitSecs(1);
    Screen('Flip',windowPtr);
    
    jitter                               = @(lb, ub) (ub - lb)*rand + lb;
    for i                                = 1:ExpInfo.numTotalTrials
        % jitter irrevalant duration
        ExpInfo.jFixation(i)                 = jitter(ExpInfo.fixation(1),ExpInfo.fixation(2));
        ExpInfo.jBlankDuration1(i)           = jitter(ExpInfo.blankDuration1(1), ExpInfo.blankDuration1(2));
        ExpInfo.jISI(i)                      = jitter(ExpInfo.ISI(1), ExpInfo.ISI(2));
        ExpInfo.jBlankDuration2(i)           = jitter(ExpInfo.blankDuration2(1), ExpInfo.blankDuration2(2));
        ExpInfo.jfeedback(i)                 = jitter(ExpInfo.feedback(1), ExpInfo.feedback(2));
        ExpInfo.jITI(i)                      = jitter(ExpInfo.ITI(1), ExpInfo.ITI(2));
 
        %present stimuli
        [Response.more(i), Response.RT(i)]...
            =PresentDiscriminationStimulus_v3(i,ExpInfo,ScreenInfo,...
            VSinfo,AudInfo,pahandle,windowPtr);
        
        % save data just in case
        save(out1FileName,'Response', 'ExpInfo', 'ScreenInfo',...
            'VSinfo', 'AudInfo', 'pahandle', 'windowPtr');
        
        % break
        if rem(i,ExpInfo.numTrialsinBlock)   ==0 && i~=ExpInfo.numTotalTrials
            DrawFormattedText(windowPtr, 'You finished one block. Please take a break.\nPress 5 to continue.',...
                'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
            Screen('Flip',windowPtr);
            while 1
                [~, ~, keyCode]                      = KbCheck();
                if keyCode(KbName('5'))              ==1
                    break;
                end
            end
            
        elseif  i                            == ExpInfo.numTotalTrials
            DrawFormattedText(windowPtr, 'Thank you for participating.\nPress any key to exit.',...
                'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
            Screen('Flip',windowPtr);
            KbWait([],2);
        end 
    end
    
    %% Finish the experiment
    % end time
    c                                    = clock;
    finish                               = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));
    timestamp{2,1}                       = finish;
    % save data
    save(out1FileName,'Response', 'ExpInfo', 'ScreenInfo',...
        'VSinfo', 'AudInfo', 'pahandle', 'windowPtr','timestamp');
    ShowCursor();
    Screen('CloseAll');
    ListenChar(0);
    
catch e
    psychError                           = psychlasterror();
    save('error.mat','e','psychError')
    ShowCursor();
    Screen('CloseAll');
    ListenChar(0);
    psychrethrow(psychlasterror);
end