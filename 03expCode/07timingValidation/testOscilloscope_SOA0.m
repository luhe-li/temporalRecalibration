%% pre-exposure: ternary temporal order judgement task

% 2022/0310/v2
% Equipment indices are set to work on Chartest with the projector, central
% speaker and number pad.

% Participants will complete a 3-alternative-force-choice (3AFC) ternary
% task (Ulrich, 1987; García-Pérez & Alcalá-Quintana, 2012; Yarrow, Martin,
% Di Costa, Solomon, & Arnold, 2016). In each trial, they will be presented
% with audiovisual stimulus pair with various stimulus onset asynchronies
% (SOA), ranging from [-600, -400: 50: 400, 600] ms. A positive SOA refers
% to a visual lead and a negative SOA refers to an auditory lead. After the
% stimulus presentation, participants will report whether they perceive
% that the auditory stimulus precede the visual stimulus (‘A-precedes-V’
% response), both stimuli occur at the same time (‘A-coincides-V’
% response), or the visual stimulus precedes the auditory stimulus
% (‘A-follows-V’ response) by a button press. Each of the SOAs will be
% randomly interleaved within a method of a constant stimuli and presented
% 25 times, resulting in a total of 475 trials.

% Response guide:
% 1: V-first
% 2: simultaneous
% 3: A-first
% numLock: exit the experiment during response period

try
    %% enter subject's name
    clear all; close all; clc; rng('shuffle');
    
    % choose device indices
    audDevices = PsychPortAudio('GetDevices');
%     audDevicesEnd = size(audDevices, 2);
%     ExpInfo.deviceIndices = [0, 1, 0]; % use internal display, speaker, keyboard;
    ExpInfo.deviceIndices = [0, 4, 1]; % use external display, speaker, keyboard;
    ExpInfo.sittingDistance = 10; %in cm, the distance between the screen and participants

    %% define parameters
    % define SOA as the time difference between two stimulus onsets
    ExpInfo.SOA = sort(repmat([100:50:300]/1000,[1,10]));%sort(repmat(a, 1, 3)); %s
    ExpInfo.numSOA = length(ExpInfo.SOA);

    % trial number per SOA level
    ExpInfo.numTrials = 10;%1;
    
    % number of trials in each block to decide break
    ExpInfo.numTrialsinBlock = 95;
    
    % sound and blob presented the same duration
    % make sure the SOA increment is larger than stimulus duration
    ExpInfo.stimFrame = 2; % frame
    
    % define duration ranges by lower bound and higher bound
    ExpInfo.fixation = [0.1, 0.2];
    ExpInfo.blankDuration1 = [0.4, 0.6]; % between fixation and the first stimulus
    ExpInfo.blankDuration2 = [0.4, 0.6]; % between stimulus offset and response probe
    ExpInfo.ITI = [0.2, 0.4];

    %% define the experiment information
    ExpInfo.numTotalTrials    = ExpInfo.numTrials * ExpInfo.numSOA;
    ExpInfo.trialSOA = ExpInfo.SOA;
    
    %initialize a structure that stores all the responses and response time
    [Response.order, Response.RT] = ...
        deal(NaN(1,ExpInfo.numTotalTrials));
    
    %% screen setup
    %myKeyCheck();
    PsychDefaultSetup(2);
    AssertOpenGL();
    GetSecs();
    WaitSecs(0.1);
    KbCheck();
    ListenChar(2); % silence the keyboard
    HideCursor();
    
    %Canvas size = 53.5" x 40"= 135.8cm x 101.6cm Screen size by the project =
    %1024 pixels x 768 pixels so each centimeter has 7.54 pixels horizontally
    %and 7.56 vertically
    Screen('Preference', 'VisualDebugLevel', 1);
    Screen('Preference', 'SkipSyncTests', 0);
    [windowPtr,rect] = Screen('OpenWindow', ExpInfo.deviceIndices(1), [255, 255, 255]./4); % 1 = external display
%     [windowPtr,rect] = Screen('OpenWindow', ExpInfo.deviceIndices(1), [255, 255, 255]./4, [100 100 1000 600]); % for testing
    [ScreenInfo.xaxis, ScreenInfo.yaxis] = Screen('WindowSize',windowPtr);
    Screen('TextSize', windowPtr, 25) ;
    Screen('TextFont',windowPtr,'Times');
    Screen('TextStyle',windowPtr,1);
    ScreenInfo.frameRate=Screen('FrameRate',0);
%     ScreenInfo.frameDura=1/ScreenInfo.frameRate; % in s
    ScreenInfo.ifi = Screen('GetFlipInterval', windowPtr);
    ScreenInfo.topPriorityLevel = MaxPriority(windowPtr);
    [center(1), center(2)]     = RectCenter(rect);
    ScreenInfo.xmid            = center(1); % horizontal center
    ScreenInfo.ymid            = center(2); % vertical center
    ScreenInfo.backgroundColor = 0;
    ScreenInfo.numPixels_perCM = 7.5;
    ScreenInfo.liftingYaxis    = 300;
  
    %fixation locations
    ScreenInfo.x1_lb = ScreenInfo.xmid-7; ScreenInfo.x2_lb = ScreenInfo.xmid-1;
    ScreenInfo.x1_ub = ScreenInfo.xmid+7; ScreenInfo.x2_ub = ScreenInfo.xmid+1;
    ScreenInfo.y1_lb = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-1;
    ScreenInfo.y1_ub = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+1;
    ScreenInfo.y2_lb = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-7;
    ScreenInfo.y2_ub = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+7;
    
    %% open loudspeakers and create sound stimuli
    addpath(genpath(PsychtoolboxRoot))
    PsychDefaultSetup(2);
    % get correct sound card
    our_device = audDevices(ExpInfo.deviceIndices(2)).DeviceIndex;

    %sampling frequencies
    AudInfo.fs   = audDevices(ExpInfo.deviceIndices(2)).DefaultSampleRate;
    audioSamples = linspace(1,AudInfo.fs,AudInfo.fs);
    
    % make a beep
    AudInfo.testIntensity = 0.5;
    AudInfo.stimDura  = ExpInfo.stimFrame * ScreenInfo.ifi; %s, the duration of auditory stimulus
    AudInfo.tf = 500;
    AudInfo.beepLengthSecs = AudInfo.stimDura;
    AudInfo.Beep = [MakeBeep(AudInfo.tf, AudInfo.beepLengthSecs, AudInfo.fs);...
        zeros(1, round(AudInfo.beepLengthSecs * AudInfo.fs))];
%             zeros(1, floor(AudInfo.beepLengthSecs * AudInfo.fs))];
    pahandle = PsychPortAudio('Open', our_device, [], [], [], 2);%open device
    PsychPortAudio('Volume', pahandle, AudInfo.testIntensity);
    
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
    VSinfo.duration        = ExpInfo.stimFrame * ScreenInfo.ifi;%s
    VSinfo.width            = 401; %(pixel) Increasing this value will make the cloud more blurry
    VSinfo.boxSize          = 201; %This is the box size for each cloud.
    VSinfo.testIntensity        = 10; %This determines the contrast of the clouds. Lowering this value will make
    %them have lower contrast
    %set the parameters for the visual stimuli
    VSinfo.blackScreen = 255/4*ones(ScreenInfo.xaxis,ScreenInfo.yaxis);
    VSinfo.whiteScreen = 255*ones(ScreenInfo.xaxis,ScreenInfo.yaxis);
    VSinfo.blankScreen = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);
    x                  = 1:1:VSinfo.boxSize; y = x;
    [X,Y]              = meshgrid(x,y);
    cloud              = 1e2.*mvnpdf([X(:) Y(:)],[median(x) median(y)],...
        [VSinfo.width 0; 0 VSinfo.width]);
    VSinfo.Cloud       = 255.*VSinfo.testIntensity.*reshape(cloud,length(x),length(y));
    VSinfo.blk_texture = Screen('MakeTexture', windowPtr, VSinfo.blackScreen,[],[],[],2);
    VSinfo.white_texture = Screen('MakeTexture', windowPtr, VSinfo.whiteScreen,[],[],[],2);
    
    %% keyboard information
    KbInfo.KbDevices = GetKeyboardIndices;
    if ExpInfo.deviceIndices(3) == 0
        KbInfo.KbDeviceIndex = KbInfo.KbDevices(1); % internal
    else
        KbInfo.KbDeviceIndex = KbInfo.KbDevices(end - 2); % external
    end
    
    %% Run the experiment by calling the function InterleavedStaircase
    %  record tart time
    c = clock;
    start = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));
    timestamp{1,1} = start;
    
    %start the experiment
    DrawFormattedText(windowPtr, 'Press any button to start the temporal order judgement task.',...
        'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
    Screen('Flip',windowPtr);
    KbWait(-3); WaitSecs(1);
    Screen('Flip',windowPtr);
    
    jitter = @(lb, ub) (ub - lb)*rand + lb;
    for i = 1:ExpInfo.numTotalTrials
        % jitter irrevalant duration
        ExpInfo.jFixation(i) = jitter(ExpInfo.fixation(1),ExpInfo.fixation(2));
        ExpInfo.jBlankDuration1(i) = jitter(ExpInfo.blankDuration1(1), ExpInfo.blankDuration1(2));
        ExpInfo.jBlankDuration2(i) = jitter(ExpInfo.blankDuration2(1), ExpInfo.blankDuration2(2));
        ExpInfo.jITI(i) = jitter(ExpInfo.ITI(1), ExpInfo.ITI(2));
        
        %present multisensory stimuli
      PresentMultisensoryStimuliTestOscilloscope_v16(i,ExpInfo,ScreenInfo,...
            VSinfo, AudInfo,pahandle,windowPtr, KbInfo);
        
    end
    
    %% Finish the experiment
    % end time
    c = clock;
    finish = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));
    timestamp{2,1} = finish;
      ShowCursor();
    Screen('CloseAll');
    PsychPortAudio('Close');
    ListenChar(0);
    
catch e
    psychError = psychlasterror();
    save('error.mat','e','psychError')
    ShowCursor();
    Screen('CloseAll');
    PsychPortAudio('Close');
    ListenChar(0);
    psychrethrow(psychlasterror);
end