%% exposure phase

% Exposure phase: participants will be exposed to a train of 250
% audiovisual stimulus pairs with a consistent temporal discrepancy between
% them (±750, ±625, ±500, ±375, ±250, ±125 ms; fixed within a session).

% During the exposure, participants were required to do a nonrelevant task
% to keep attention to both modalities. Across the session, an oddball
% stimulus with higher intensity would be presented at 10% chance. Half of
% the oddball stimuli were the visual stimulus with higher contrast, and
% the other half were auditory stimuli with higher volume. In every trial,
% after the stimulus presentation, a post-cue ("V” or “A”) would appear to
% ask whether the stimulus was different than the average of stimuli
% they've met in the relevant modality.

% There are 12 trials in this practice phase. Adaptor duration is the same
% as the session adaptor. This practice will show 3 common AV pairs, 2
% visual oddballs, 2 common AV pairs, 2 auditory oddballs, and 3 commom AV
% pairs. There will be notification whether the next trial will be an
% oddball or not before every stimulus presentation, which will be not
% shown in the formal experiment.

% Response guide:
% Respond by clicking the mouse
% left: not oddball
% right: oddball

try
    %% enter subject's name
    clear all; close all; clc; rng('shuffle');
    
    %enter subject's name
    ExpInfo.subjID = [];
    while isempty(ExpInfo.subjID) == 1
        try ExpInfo.subjID = input('Please enter participant ID#: ') ;
            ExpInfo.session = input('Please enter exposure session#: ');
        catch
        end
    end
    
    % choose device indices
    audDevices = PsychPortAudio('GetDevices');
    audDevicesEnd = size(audDevices, 2);
    ExpInfo.deviceIndices = [0, audDevicesEnd, 0]; % use internal display, speaker, keyboard;
    %     ExpInfo.deviceIndices = [1, 1, 1]; % use external display, speaker, keyboard;
    
    ExpInfo.sittingDistance = 10; %in cm, the distance between the screen and participants
    out1FileName = ['exposure_sub', num2str(ExpInfo.subjID) '_practice' num2str(ExpInfo.session)]; %create file name
    
    % avoid rewriting data
    if exist([out1FileName '.mat'],'file')
        resp=input('To replace the existing file, press y', 's');
        if ~strcmp(resp,'y')
            disp('Experiment stopped.')
            return
        end
    end
    
    %% define parameters
    % this section is for easy changing unfixed parameters in testing and
    % will be removed once they are fixed
    
    % define SOA/adaptor as the time difference between two stimulus
    % onsets;
    % make sure adaptors are divisible by frameDuration; if not, frame
    % number will be floored;
    % load counterbalance chart and select the adaptor by session
    load('RandomAdaptorOrder.mat')
    ExpInfo.adaptor = randomAdaptorOrder(ExpInfo.session); %secs
    
    % nTTrials: total number of trials
    nTTrials = 12;
    
    % nTBlock: number of trials in each block
    
    
    % sound and blob presented the same duration
    % make sure the SOA increment is larger than stimulus duration
    ExpInfo.stimFrame = 2; % frame
    
    % ITI duration
    ExpInfo.fixation = [0.1, 0.2];
    ExpInfo.blankDuration1 = [0.4, 0.6]; % between fixation and the first stimulus
    ExpInfo.blankDuration2 =  [0.4, 0.6]; % between the second stimulus and post-cue
    ExpInfo.ITI = [0.2, 0.4];
    
    %% screen setup
    %myKeyCheck();
    AssertOpenGL();
    GetSecs();
    WaitSecs(0.1);
    KbCheck();
    ListenChar(2); % silence the keyboard
    
    %Canvas size = 53.5" x 40"= 135.8cm x 101.6cm Screen size by the project =
    %1024 pixels x 768 pixels so each centimeter has 7.54 pixels horizontally
    %and 7.56 vertically
    Screen('Preference', 'VisualDebugLevel', 1);
    Screen('Preference', 'SkipSyncTests', 1);
%         [windowPtr,rect] = Screen('OpenWindow', ExpInfo.deviceIndices(1), [0,0,0]); % 1 = external display
    [windowPtr,rect] = Screen('OpenWindow', ExpInfo.deviceIndices(1), [0,0,0],[100 100 1000 600]); % for testing
    [ScreenInfo.xaxis, ScreenInfo.yaxis] = Screen('WindowSize',windowPtr);
    Screen('TextSize', windowPtr, 35) ;
    Screen('TextFont',windowPtr,'Times');
    Screen('TextStyle',windowPtr,1);
    %     ScreenInfo.frameRate=Screen('FrameRate',0);
    ScreenInfo.frameRate=60;
    ScreenInfo.frameDura=1/ScreenInfo.frameRate; % secs
    ScreenInfo.ifi = round(1000/(ScreenInfo.frameRate))/1000;  % inter-flip-interval in sec
    
    [center(1), center(2)]     = RectCenter(rect);
    ScreenInfo.xmid            = center(1); % horizontal center
    ScreenInfo.ymid            = center(2); % vertical center
    ScreenInfo.backgroundColor = 0;
    ScreenInfo.numPixels_perCM = 7.5;
    %     ScreenInfo.liftingYaxis    = 300;
    ScreenInfo.liftingYaxis    =  ScreenInfo.ymid;
    ScreenInfo.cursorColor     = [0,0,255]; %A: blue, V:red
    ScreenInfo.dispModality    = 'A'; %always localize the auditory component
    ScreenInfo.x_box_unity     = [-95, -32; 35, 98];
    ScreenInfo.y_box_unity     = [-10, 22];
    
    %fixation locations
    ScreenInfo.x1_lb = ScreenInfo.xmid-7; ScreenInfo.x2_lb = ScreenInfo.xmid-1;
    ScreenInfo.x1_ub = ScreenInfo.xmid+7; ScreenInfo.x2_ub = ScreenInfo.xmid+1;
    ScreenInfo.y1_lb = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-1;
    ScreenInfo.y1_ub = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+1;
    ScreenInfo.y2_lb = ScreenInfo.yaxis-ScreenInfo.liftingYaxis-7;
    ScreenInfo.y2_ub = ScreenInfo.yaxis-ScreenInfo.liftingYaxis+7;
    
    % fixation dot
    r1 = 4;
    ScreenInfo.fixDot = [ScreenInfo.xmid-r1 ScreenInfo.ymid-r1...
        ScreenInfo.xmid+r1 ScreenInfo.ymid+r1];
    
    %% open loudspeakers and create sound stimuli
    addpath(genpath(PsychtoolboxRoot))
    PsychDefaultSetup(2);
    % get correct sound card
    InitializePsychSound
    devices = PsychPortAudio('GetDevices');
    our_device = audDevices(ExpInfo.deviceIndices(2)).DeviceIndex;
    
    %sampling frequencies
    AudInfo.fs   = audDevices(ExpInfo.deviceIndices(2)).DefaultSampleRate;
    audioSamples = linspace(1,AudInfo.fs,AudInfo.fs);
    
    % make a beep
    AudInfo.stimDura  = ExpInfo.stimFrame * ScreenInfo.frameDura; %s, the duration of auditory stimulus
    AudInfo.tf = 500;
    AudInfo.beepLengthSecs = 1;
    AudInfo.Beep = [zeros(1,AudInfo.beepLengthSecs*AudInfo.fs);...
        MakeBeep(AudInfo.tf , AudInfo.beepLengthSecs, AudInfo.fs)];
    pahandle                    = PsychPortAudio('Open', our_device, [], [], [], 2);%open device
    
    %% define the visual stimuli
    VSinfo.duration        = ExpInfo.stimFrame * ScreenInfo.frameDura;%s
    VSinfo.width            = 401; %(pixel) Increasing this value will make the cloud more blurry
    VSinfo.boxSize          = 201; %This is the box size for each cloud.
    VSinfo.intensity        = 10; %This determines the contrast of the clouds. Lowering this value will make
    %them have lower contrast
    %set the parameters for the visual stimuli
    VSinfo.blackScreen = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);
    VSinfo.blankScreen = zeros(ScreenInfo.xaxis,ScreenInfo.yaxis);
    x                  = 1:1:VSinfo.boxSize; y = x;
    [X,Y]              = meshgrid(x,y);
    cloud              = 1e2.*mvnpdf([X(:) Y(:)],[median(x) median(y)],...
        [VSinfo.width 0; 0 VSinfo.width]);
    VSinfo.Cloud       = 255.*VSinfo.intensity.*reshape(cloud,length(x),length(y));
    VSinfo.blk_texture = Screen('MakeTexture', windowPtr, VSinfo.blackScreen,[],[],[],2);
    
    %% keyboard information
    KbInfo.KbDevices = GetKeyboardIndices;
    if ExpInfo.deviceIndices(3) == 0
        KbInfo.KbDeviceIndex = KbInfo.KbDevices(1); % internal
    else
        KbInfo.KbDeviceIndex = KbInfo.KbDevices(end); % external
    end

    %% define the experiment information
    ExpInfo.nTTrials         = nTTrials;% total trial number per session
    
    ExpInfo.respModality = repmat('V', [1,nTTrials]);
    ExpInfo.respAidx = randsample(nTTrials, nTTrials/2)'; % randomly choose half of trials to report auditory intensity
    ExpInfo.respModality(ExpInfo.respAidx) = 'A';
    ExpInfo.respVidx = find(ExpInfo.respModality == 'V');% the index of the other half trials to report visual intensity
    ExpInfo.respModality(4) = 'V';ExpInfo.respModality(5) = 'A';
    ExpInfo.respModality(8) = 'A';ExpInfo.respModality(9) = 'V';
    
    % test intensity manipulation with oddball in 5:10 trials
    ExpInfo.VIntensity = repmat(10,[1,nTTrials]);
    ExpInfo.VIntensity(4) = 20; ExpInfo.VIntensity(5) = 20;
    ExpInfo.AIntensity = repmat(0.3,[1,nTTrials]); % original volumn is 0.3
    ExpInfo.AIntensity(8) = 1; ExpInfo.AIntensity(9) = 1;
    
    %initialize a structure that stores all the responses and response time
    [Response.oddball, Response.RT] = ...
        deal(NaN(1,ExpInfo.nTTrials));
    
    %% Run the experiment by calling the function InterleavedStaircase
    %  record tart time
    c = clock;
    start = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));
    timestamp{1,1} = start;
    
    %start the experiment
    DrawFormattedText(windowPtr, 'Press any button to start the intensity oddball task.',...
        'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
    Screen('Flip',windowPtr);
    KbWait(-3); WaitSecs(1);
    Screen('Flip',windowPtr);
    
       jitter = @(lb, ub) (ub - lb)*rand + lb;
    for i = 1:ExpInfo.nTTrials
          % jitter irrevalant duration
        ExpInfo.jExpoFixation(i) = jitter(ExpInfo.fixation(1),ExpInfo.fixation(2));
        ExpInfo.jExpoBlankDuration1(i) = jitter(ExpInfo.blankDuration1(1), ExpInfo.blankDuration1(2));
        ExpInfo.jExpoBlankDuration2(i) = jitter(ExpInfo.blankDuration2(1), ExpInfo.blankDuration2(2));
        ExpInfo.jExpoITI(i) = jitter(ExpInfo.ITI(1), ExpInfo.ITI(2));
        
        % post-cue
        ExpInfo.postCue = ExpInfo.respModality(i);
        
        % set stimulus intensity depending on whether it's an oddball
        % visual stim
        VSinfo.intensity   = ExpInfo.VIntensity(i);
        VSinfo.Cloud       = 255.*VSinfo.intensity.*reshape(cloud,length(x),length(y));
        
        % auditory stim
        PsychPortAudio('Volume', pahandle, ExpInfo.AIntensity(i));
        
        %present multisensory stimuli
        [Response.oddball(i), Response.RT(i)]...
            = PresentMultisensoryStimuliExposure(i,ExpInfo,ScreenInfo,...
            VSinfo, AudInfo,pahandle,windowPtr);
        
        % save data just in case
        save(out1FileName,'Response', 'ExpInfo', 'ScreenInfo',...
            'VSinfo', 'AudInfo', 'pahandle', 'windowPtr');
        
        % display feedback
         if i == 4 || i == 5
            DrawFormattedText(windowPtr, 'It was a visual oddball.',...
                'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
            Screen('Flip',windowPtr); WaitSecs(1);
        elseif i == 8 || i == 9
            DrawFormattedText(windowPtr, 'It was an anditory oddball.',...
                'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
            Screen('Flip',windowPtr); WaitSecs(1);
        else 
            DrawFormattedText(windowPtr, 'It was NOT an oddball.',...
                'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
            Screen('Flip',windowPtr); WaitSecs(1);
        end
    end
    
    %% Finish the experiment
    % end time
    c = clock;
    finish = sprintf('%04d/%02d/%02d_%02d:%02d:%02d',c(1),c(2),c(3),c(4),c(5),ceil(c(6)));
    timestamp{2,1} = finish;
    % save data
    save(out1FileName,'Response', 'ExpInfo', 'ScreenInfo',...
        'VSinfo', 'AudInfo', 'pahandle', 'windowPtr','timestamp');
    ShowCursor();
    Screen('CloseAll');
    ListenChar(0);
    
catch e
    psychError = psychlasterror();
    save('error.mat','e','psychError')
    ShowCursor();
    Screen('CloseAll');
    ListenChar(0);
    psychrethrow(psychlasterror);
end