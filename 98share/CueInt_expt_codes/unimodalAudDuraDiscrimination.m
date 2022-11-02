% unimodal duration discrimination task

% response guide:
% 1: the first interval was longer
% 2: the second interval was longer

try
    %% enter subject's name

    clear all; close all; clc; rng('Shuffle');

    %enter subject's name
    ExpInfo.subjID = [];
    while isempty(ExpInfo.subjID)  == 1
        try ExpInfo.subjID  = input('Please enter participant ID#: ') ; %'s'
            ExpInfo.session  = input('Please enter session#: ');
            ExpInfo.mode  = input('Experiment mode: 1; Debug mode: 2#: ');
            ExpInfo.practice  = input('Experiment mode: 1; Practice mode: 2#: ');
        catch
        end
    end

    % switch practice mode
    switch ExpInfo.practice
        case 1 % experiment mode
            ExpInfo.nTrials = 20; % number of trial per condition level
            out1FileName  = ['pilot1_sub', num2str(ExpInfo.subjID) '_ses' num2str(ExpInfo.session)]; %create file name
        case 2 % practice mode
            ExpInfo.nTrials = 1; % number of trial per condition level
            out1FileName  = ['practice_pilot1_sub', num2str(ExpInfo.subjID) '_ses' num2str(ExpInfo.session)]; %create file name
    end

    % avoid rewriting data
    if exist([out1FileName '.mat'],'file')
        resp                                 = input('To replace the existing file, press y', 's');
        if ~strcmp(resp,'y')
            disp('Experiment stopped.')
            return
        end
    end

    %% switch between debug mode

    switch ExpInfo.mode
        case 1 % experiment mode (in the experiment room)
            ExpInfo.deviceIndices = [0, 4, 1]; % use external display, speaker, keyboard;
            windowSize = [];
            Screen('Preference', 'SkipSyncTests', 0);
            HideCursor();
        case 2 % debug mode
            ExpInfo.deviceIndices = [0, 2, 0]; % use internal display, speaker, keyboard;
            windowSize = [100 100 1000 600]; % open a smaller window
            Screen('Preference', 'SkipSyncTests', 1);
    end

    %% define parameters

    % define duration by frame
    %     ExpInfo.ifi = 1/60; % in sec, use ScreenInfo.ifi later instead
    %     ExpInfo.duraFrame = 30; % 0.5 sec, standard stimulus duration
    %     ExpInfo.compRatioFrame = [-9, -3, -1, 0, 1, 3, 9];
    %     ExpInfo.compFrame = ExpInfo.compRatioFrame + 30; % comparison stimulus frame
    ExpInfo.duraStandard = 0.5; % in sec, standard stimulus duration
    ExpInfo.compRatio = [-0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5] + 1; % percentage
    ExpInfo.duraComp = ExpInfo.compRatio .* ExpInfo.duraStandard; % comparison stimulus duration

    % define trial numbers
    ExpInfo.numComp = length(ExpInfo.duraComp); % number of comparison conditions
    ExpInfo.numTrialsinBlock = 60; % number of trials in each block to decide break

    % define auditory stimulus amplitude
    % manipulation of temporal reliability by sampling amplitudes from uniform dist
    ExpInfo.bin = 1/100; % each bin is 1/100 sec
    ExpInfo.samp_rate = 48000; % overwritten by real sampling rate later
    ExpInfo.min_noise = 1; % lower bound of noise
    ExpInfo.min_signal = 3; % lower bound of signal
    ExpInfo.sigmas = [0, 2, 3]; % range of uniform dist, larger sigma/range/uncertainty
    ExpInfo.numSigmas = length(ExpInfo.sigmas); % number of reliability conditions

    % define auditory stimulus intensity (scaling of all amplitudes)
    AudInfo.standard = 0.02; % a ratio between 0 to 1

    % define visual stimulus intensity
    VSinfo.off = 0.2; % a ratio between 0 to 1 to be multipled by 255
    VSinfo.on = 0.5; % a ratio between 0 to 1 to be multipled by 255
    VSinfo.pblack = 1/8; % set contrast to 255*1/8 for the "black" background

    % define equipment temporal bias between visual and audio
    % representation
    ExpInfo.bias  = 0.03253;
    ExpInfo.sittingDistance  = 10; %in cm, the distance between the screen and participants

    % define irrelavant duration ranges by lower bound and higher bound
    ExpInfo.trialDura = 3.5; % s
    ExpInfo.fixation  = [0.1, 0.3];
    ExpInfo.ISI = [0.5, 1];
    ExpInfo.ITI  = [0.2, 0.4];

    %% define the experiment information

    ExpInfo.numTotalTrials               = ExpInfo.nTrials * ExpInfo.numComp * ExpInfo.numSigmas;
    ExpInfo.cond = combvec([1:ExpInfo.numComp], [1:ExpInfo.numSigmas]);
    for i = 1:ExpInfo.nTrials % repeat all conditions for each nTrial
        tempOrder = randperm(length(ExpInfo.cond));
        ExpInfo.trialCompDura(:,i) = ExpInfo.duraComp(ExpInfo.cond(1, tempOrder));
        ExpInfo.trialSigma(:,i) = ExpInfo.sigmas(ExpInfo.cond(2, tempOrder));
        ExpInfo.trialOrder(:,i) = tempOrder;
    end
    ExpInfo.trialOrder = reshape(ExpInfo.trialOrder, [], 1)';
    ExpInfo.trialSigma = reshape(ExpInfo.trialSigma, [], 1)';
    ExpInfo.trialCompDura = reshape(ExpInfo.trialCompDura, [], 1)';

    % define counterbalance array
    ExpInfo.counterbalance               = ones(1,ExpInfo.numTotalTrials); % 1 = standard first
    idxComp                              = randsample(ExpInfo.numTotalTrials, ExpInfo.numTotalTrials/2); % randomly choose half of trials to present comparison stimulus first
    ExpInfo.counterbalance(idxComp)      = 2; % 2 = comparison first

    %initialize a structure that stores all the responses and response time
    [Response.order, Response.RT] = deal(NaN(1,ExpInfo.numTotalTrials));
    ExpInfo.realDura = NaN(ExpInfo.numTotalTrials, 5); % record real durations within each trial, columns are: pre-stim, stim1, ISI, stim2, post-stim

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
    [windowPtr,rect] = Screen('OpenWindow', ExpInfo.deviceIndices(1), [255, 255, 255] * VSinfo.pblack, windowSize);
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

    % initialize driver, request low-latency preinit
    audDevices  = PsychPortAudio('GetDevices');
    InitializePsychSound(1);
    addpath(genpath(PsychtoolboxRoot))

    % get the correct sound card
    our_device = audDevices(ExpInfo.deviceIndices(2)).DeviceIndex;

    % make a warmup beep
    ExpInfo.sampRate = audDevices(ExpInfo.deviceIndices(2)).DefaultSampleRate;
    tf = 500;
    beepLengthSecs = 0.01;
    beep    = MakeBeep(tf, beepLengthSecs, ExpInfo.sampRate);
    beeps     = [beep; zeros(size(beep))];
    pahandle  = PsychPortAudio('Open', our_device, [], [], ExpInfo.sampRate, 2);%open device
    %     pahandle = PsychPortAudio('Open', our_device, 1, 1, ExpInfo.sampRate, 2);

    % Perform one warmup trial, to get the sound hardware fully up and running,
    % performing whatever lazy initialization only happens at real first use.
    % This "useless" warmup will allow for lower latency for start of playback
    % during actual use of the audio driver in the real trials:
    PsychPortAudio('Volume', pahandle, AudInfo.standard);
    PsychPortAudio('FillBuffer', pahandle, beeps);
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    PsychPortAudio('Stop', pahandle, 1);

    %% define the visual stimuli
    VSinfo.width                         = 401; %(pixel) Increasing this value will make the cloud more blurry
    VSinfo.boxSize                       = 201; %This is the box size for each cloud.
    %set the parameters for the visual stimuli
    VSinfo.blackScreen                   = 255 * VSinfo.pblack *ones(ScreenInfo.xaxis,ScreenInfo.yaxis);
    VSinfo.blankScreen                   = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);
    x                                    = 1:1:VSinfo.boxSize; y = x;
    VSinfo.x                             = x; VSinfo.y = y;
    [X,Y]                                = meshgrid(x,y);
    VSinfo.cloud                         = mvnpdf([X(:) Y(:)],[median(x) median(y)],...
        [VSinfo.width 0; 0 VSinfo.width]);
    VSinfo.pscale                               = (1-VSinfo.pblack)/max(VSinfo.cloud); % the max contrast of the blob adds the background contrast should <= 1
    VSinfo.cloud                         = VSinfo.cloud .* VSinfo.pscale;
    VSinfo.Cloud                         = 255.*VSinfo.off.*reshape(VSinfo.cloud,length(x),length(y));
    VSinfo.blk_texture                   = Screen('MakeTexture', windowPtr, VSinfo.blackScreen,[],[],[],2);

    %% Run the experiment by calling the function InterleavedStaircase
    %  record tart time
    c        = clock;
    start        = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));
    timestamp{1,1}        = start;

    % start the experiment
    DrawFormattedText(windowPtr, 'Press any button to start the duration judgement task.',...
        'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
    Screen('Flip',windowPtr);
    KbWait(); WaitSecs(1);
    Screen('Flip',windowPtr);

    for i = 1:ExpInfo.numTotalTrials

        %present multisensory stimuli
        [Response.order(i), Response.RT(i), ExpInfo]...
            =presentAudInterval_v1(i,...
            ExpInfo, ScreenInfo, VSinfo, AudInfo, pahandle, windowPtr);

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