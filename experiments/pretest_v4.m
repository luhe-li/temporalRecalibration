%% pre-exposure: ternary temporal order judgement task
% 2022/04/v4
% Equipment indices are set to work on Chartest with the projector, central
% speaker and number pad.

% For calibration, Psychportaudio was delayed by 32.53 ms.

% Participants will complete a 3-alternative-force-choice (3AFC) ternary
% task (Ulrich, 1987; García-Pérez & Alcalá-Quintana, 2012; Yarrow, Martin,
% Di Costa, Solomon, & Arnold, 2016). In each trial, they will be presented
% with audiovisual stimulus pair with various stimulus onset asynchronies
% (SOA), ranging from [-500 -300 -250 -200 -150 -100 -50 0 50 100 150 200 250 300 500]...
% ms. A positive SOA refers to a visual lead
% and a negative SOA refers to an auditory lead. After the stimulus
% presentation, participants will report whether they perceive that the
% auditory stimulus precede the visual stimulus (‘A-precedes-V’ response),
% both stimuli occur at the same time (‘A-coincides-V’ response), or the
% visual stimulus precedes the auditory stimulus (‘A-follows-V’ response)
% by a button press. Each of the SOAs will be randomly interleaved within a
% method of a constant stimuli and presented 25 times, resulting in a total
% of 375 trials.

% Response guide:
% 1: V-first
% 2: simultaneous
% 3: A-first
% numLock: exit the experiment during response period

try
    %% enter subject's name
    clear all; close all; clc; rng('Shuffle');
    
    %enter subject's name
    ExpInfo.subjID                       = [];
    while isempty(ExpInfo.subjID)        == 1
        try ExpInfo.subjID                   = input('Please enter participant ID#: ') ; %'s'
            ExpInfo.session                      = input('Please enter session#: ');
        catch
        end
    end
    
    % choose device indices
    audDevices                           = PsychPortAudio('GetDevices');
%     ExpInfo.deviceIndices = [0, 1, 0]; % use internal display, speaker, keyboard;
    ExpInfo.deviceIndices                = [0, 4, 1]; % use external display, speaker, keyboard;
    ExpInfo.sittingDistance              = 10; %in cm, the distance between the screen and participants
    out1FileName                         = ['pretest_sub', num2str(ExpInfo.subjID) '_session' num2str(ExpInfo.session)]; %create file name
    
    % avoid rewriting data
    if exist([out1FileName '.mat'],'file')
        resp                                 = input('To replace the existing file, press y', 's');
        if ~strcmp(resp,'y')
            disp('Experiment stopped.')
            return
        end
    end
    
    %% define parameters
    % define bias
    ExpInfo.bias                         = 0.03253;

    % define SOA as the time difference between two stimulus onsets
%     ExpInfo.SOA                          = [-400, -300:50:300, 400] /1000; %s
    ExpInfo.SOA                          = [-500 -300 -250 -200 -150 -100 -50 0 50 100 150 200 250 300 500]/1000; %s
    ExpInfo.numSOA                       = length(ExpInfo.SOA);

    % trial number per SOA level
    ExpInfo.nTrials                      = 20;
    
    % number of trials in each block to decide break
    ExpInfo.numTrialsinBlock             = 60;
    
    % sound and blob presented the same duration
    % make sure the SOA increment is larger than stimulus duration
    ExpInfo.stimFrame                    = 2; % frame
    
    % intensity of standard stimulus
    AudInfo.standard                     = 0.5; % a ratio between 0 to 1
    VSinfo.standard                      = 0.2; % a ratio between 0 to 1 to be multipled by 255
    pblack                               = 1/8; % set contrast to 255*1/8 for the "black" background

    % define duration ranges by lower bound and higher bound
    ExpInfo.fixation                     = [0.1, 0.2];
    ExpInfo.blankDuration1               = [0.4, 0.6]; % between fixation and the first stimulus
    ExpInfo.blankDuration2               = [0.4, 0.6]; % between stimulus offset and response probe
    ExpInfo.ITI                          = [0.2, 0.4];
    
    %% define the experiment information
    ExpInfo.numTotalTrials               = ExpInfo.nTrials * ExpInfo.numSOA;
    for nt                               = 1:ExpInfo.nTrials
        ExpInfo.trialOrder(:,nt)             = randperm(ExpInfo.numSOA); % SOA idx for each trial, call SOA value by ExpInfo.numSOA(ExpInfo.trialOrder(trial))
    end
    ExpInfo.trialOrder                   = reshape(ExpInfo.trialOrder, [], 1)';
    ExpInfo.trialSOA                     = ExpInfo.SOA(ExpInfo.trialOrder);
    
    %initialize a structure that stores all the responses and response time
    [Response.order, Response.RT]        = ...
        deal(NaN(1,ExpInfo.numTotalTrials));
    
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
  
    %fixation locations
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
    InitializePsychSound
    our_device                           = audDevices(ExpInfo.deviceIndices(2)).DeviceIndex;
    
    % sampling frequencies
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
    
    % initialize driver, request low-latency preinit:
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
    
    % start the experiment
    DrawFormattedText(windowPtr, 'Press any button to start the temporal order judgement task.',...
        'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
    Screen('Flip',windowPtr);
    KbWait(); WaitSecs(1);
    Screen('Flip',windowPtr);
    
    jitter                               = @(lb, ub) (ub - lb)*rand + lb;
    for i                                = 1:ExpInfo.numTotalTrials
        % jitter irrevalant duration
        ExpInfo.jFixation(i)                 = jitter(ExpInfo.fixation(1),ExpInfo.fixation(2));
        ExpInfo.jBlankDuration1(i)           = jitter(ExpInfo.blankDuration1(1), ExpInfo.blankDuration1(2));
        ExpInfo.jBlankDuration2(i)           = jitter(ExpInfo.blankDuration2(1), ExpInfo.blankDuration2(2));
        ExpInfo.jITI(i)                      = jitter(ExpInfo.ITI(1), ExpInfo.ITI(2));
        
        %present multisensory stimuli
        [Response.order(i), Response.RT(i)]...
        =PresentMultisensoryStimuliTest_v3(i,ExpInfo,ScreenInfo,...
            VSinfo, AudInfo,pahandle,windowPtr);
        
        % save data just in case
        save(out1FileName,'Response', 'ExpInfo', 'ScreenInfo',...
            'VSinfo', 'AudInfo', 'pahandle', 'windowPtr');
        
        % break
        if rem(i,ExpInfo.numTrialsinBlock)   == 0 && i~=ExpInfo.numTotalTrials
            DrawFormattedText(windowPtr, 'You finished one block. Please take a break.\nPress 5 to continue.',...
                'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
            Screen('Flip',windowPtr);
            while 1
                [~, ~, keyCode]                      = KbCheck();
                if keyCode(KbName('5'))              == 1
                    break;
                end
            end
        elseif  i                            == ExpInfo.numTotalTrials
            DrawFormattedText(windowPtr, 'Thank you for participating!\nPress any key to exit.',...
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