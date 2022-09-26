%% Continuous Random Dots Motion task
% Code written by Simon P. Kelly and data collected Hannah Craddock in 2016
% for her Master's thesis, 2016:
%   'Electrophysiological and behavioural indices of decision criterion
%   adjustments across context of weak and strong evidence.'
% Current code is additionally annotated by A.C. Geuzebroek.
%
% Participants were asked to continous monitor a cloud of randomly moving
% dots for intermittered target defined as upwards coherent moving dots.
% Afterwhich they were asked to response as fast and accurate as possible.
% Difficulty context was manipulated with:
%   1) Hard  (25% motion coherence)
%   2) Easy  (70% motion coherence)
%   3) Mixed (25% and 70% motion coherence with equal probabiltiy).
%
% NOTE this code will include an steady state visually evoked potentials,
% e.g. SSVEP of CONTRAST not of MOTION COHERENCE. As this was not task
% specific it might is not very informative for this specific experiment. 
%
% The High-density EEG data are recorded in sampling rate of 512Hz using a
% 128-channel EEG BioSemi in University College Dublin.
%
% IMPORTANT for this code to run you will need to do the following:
%   1)	ensure that the physical width and viewing distance are set
%       appropiately.
%   2)  Check if the triggers for the EEG are sent properly. You will have
%       to set the port to the one you are using. 
% DEPENDABLES:
%   
%
% UNFORTUANTLY DOT POSITION WAS NOT SAVED ORGINIALLY!!!

try
    %%  -------------------------------------------------------------------
    %                           BASIC SET-UP
    %   -------------------------------------------------------------------
    
    % PsychDebugWindowConfiguration     % debugging mode, trans. screen
    clear all
    commandwindow  % puts focus on command window, so if any keyboard buttons are pressed, they'll come up in the comand window rather than in your program!
            
    % Opens a graphics window on the main monitor
    whichScreen = 1;
    
    monitorwidth_cm = 40;   % physical width of monitor in cm
    dist_cm         = 57;	% viewing distance in cm
    
    [scresw, scresh] = Screen('WindowSize', whichScreen); % Get screen resolution
    center  = [scresw scresh]/2;    % useful to have the pixel coordinates of the very center of the screen (usually where you have someone fixate)
    fixRect = [center-2 center+2];  % fixation dot
    
    hz = Screen('FrameRate', whichScreen, 1);

    cm2px = scresw/monitorwidth_cm;  % multiplication factor to convert cm to pixels
    deg2px = dist_cm*cm2px*pi/180;	 % multiplication factor to convert degrees to pixels (uses aproximation tanT ~= T).
    
    
    %%  -------------------------------------------------------------------
    %                   PRESETTING CONDITION PARAMETERS
    %   -------------------------------------------------------------------
    % (note many of the time values in ms won't work out exactly - depends
    % on refresh rate & flicker - it'll be the closest it can be, should be
    % checked in EEG)
    
    par.runID       ='999'; % participant name - NOTE never use initials.
    par.recordEEG   = 0;  	% 0,1 record EEG. 
    par.recordEL    = 0;    % 0,1 to use eyelink
    
    dlg_title = 'Dot Motion Task';
    while 1   
        prompt = {'Enter SUBJECT/RUN/TASK IDENTIFIER:','EEG? (1=yes, 0=no)', 'Condition (Easy, Hard, Mixed)'};
        def     = {par.runID, num2str(par.recordEEG), 'Easy'};
        answer  = inputdlg(prompt, dlg_title, 1, def);
     
        par.runID     = answer{1};
        par.recordEEG = str2num(answer{2});
        par.Condition = answer{3};
        
        if exist([par.runID '.mat'],'file')
            dlg_title = [par.runID '.mat EXISTS ALREADY - CHOOSE ANOTHER, OR DELETE THAT ONE IF IT IS RUBBISH'];
        else
            break;
        end
        clear prompt def answer
    end
    
    %  ------------   FURTHER SCEEN PARAMETERS  ---------------------------
    par.videoFrate  = 60; % Monitor refresh rate
    if hz ~= par.videoFrate
       cleanup; error('Screen is not set at the right refresh rate, this will have consequences to the dot speed and the SSVEP!');
    end
    
    par.BGcolor     = 0;
    
    
    %  ------------   EXPERIMENTAL PARAMETERS   ---------------------------
    % temporal structure
    par.numtargets       = 24;
    par.targetDur        = 1;         % Target duration in sec.
    par.secs_btw_targs   = [2 4 6 8]; % In sec.
    par.MinTimeBtwClicks = 0.5;       % In sec., When polling for mouse clicks, don't want to be counting the same click any more than once.
    par.leadintime       = 1000;
    
    % Fixation
    par.FixWinSize = 3;    % Degrees, radius of fixation (circular).

    %  ------------------   DOT PARAMETERS  -------------------------------
    par.numPatches = 1;     % number of dot patches. 
    
    % Now enter parameters for coherence transitions for each condition:
    par.correctTransitionTimes = {[0] [0]}; % one entry per condition (not counting motion direction as a condition)
    par.condProportn = 1; % proportions - one number for each ROW of correctCoh and counterCoh (coherence waveform condition)

    switch  lower(par.Condition)
        case 'easy'
            par.correctCoh = {70};    % all the same - no actual transitions at the transition times
        case 'hard'
            par.correctCoh = {25};    % all the same - no actual transitions at the transition times
        case 'mixed'
            par.correctCoh   = {25, 70};    % all the same - no actual transitions at the transition times
            par.condProportn = [1, 1]; % proportions - one number for each ROW of correctCoh and counterCoh (coherence waveform condition)
    end
    
    par.counterTransitionTimes = {[] []};
    par.counterCoh             = {[] []};    % no counter motion AND later counterpulse
    
    % physical dots 
    par.dotspeed = 6;           % in degrees per second
    par.dotsize  = 6;           % in pixels
    par.numdots  = 150;         % this will be the number in a square; the ones outside the circle will be taken away
    par.motionDir    = 90;         % Degrees relative to positive x-axis (0 = rightward)
    par.dotpatchsize = 8;       % Degrees of visual angle - diameter of dot patch
    par.patchloc     = [0 0];   % Degrees of visual angle, relative to center
    par.rgbGRAY      = [1 1 1]*255;  % white
    
    par.FlickFdots = 15;    % SSVEP frequency (NOTE contrast)

    %  ------------------   SET UP DEVICES  -------------------------------
    HideCursor;
    window = Screen('OpenWindow', whichScreen, par.BGcolor);
    
    if par.recordEL
        par.TgWinSize = 3;    % RADIUS of fixation (circular) window in degrees
        ELsetupCalib
        Eyelink('Command', 'clear_screen 0')
        Eyelink('command', 'draw_box %d %d %d %d 15', center(1)-deg2px*par.FixWinSize, center(2)-deg2px*par.FixWinSize, center(1)+deg2px*par.FixWinSize, center(2)+deg2px*par.FixWinSize);
    end
    
    if par.recordEEG   % EEG port set-up for triggers
        port  = hex2dec('C010'); % NOTE this requires to be changed depending on the system used. 
        ioObj = io32;
        Stat  = io32(ioObj);
    end
    
    % trigger codes - can only use these 15: [1 4 5 8 9 12 13 16 17 20 21 24 25 28 29]  (don't ask!)
    par.CD_RESP         = 1;    % start response period
    par.CD_DOTS_ON      = 5;	% start ITI
    par.CD_COHMOTION_ON = 9;    % start target onset
    par.CD_BUTTONS      = [12 13];   % left and right mouse [NOT ONLY left relevant here]
    
    
    %%  -------------------------------------------------------------------
    %                   SET-UP dot motion
    %   -------------------------------------------------------------------
    % coherence, motion direction waveforms
    nMD = length(par.motionDir);    % number of motion directions
    nDP = size(par.patchloc,1);     % number of dot patches
    
    % flicker parameters:
    par.FlickTdots = round(par.videoFrate./par.FlickFdots);	% Period of flicker in number of video frames
    nrefON  = floor(par.FlickTdots/2);	% number of refreshes where the dots are ON
    nrefOFF = par.FlickTdots-nrefON;    % number of refreshes where the dots are OFF
    
    % For-loop to pre-set the timecourse of the continous randomly moving
    % dots with the intermitted targets! This is done, per ITI, per
    % coherence level, per patch (here one patch)
    numcond = 0;
    clear conditiondescrip coh trigger triggerFr correctDir patchWithMo
    
    for currITI = 1:length(par.secs_btw_targs)
        for currCoh = 1:length(par.correctCoh)
            for currPatch = 1:par.numPatches
                numFr1 = round((par.targetDur+par.secs_btw_targs(currITI)).*par.FlickFdots(currPatch));
                correctTransitionFr = round((par.secs_btw_targs(currITI) + par.correctTransitionTimes{currCoh}/1000)*par.FlickFdots(currPatch))+1;
                counterTransitionFr = round((par.secs_btw_targs(currITI) + par.counterTransitionTimes{currCoh}/1000)*par.FlickFdots(currPatch))+1;
                
                for currDir = 1:nMD
                    numcond = numcond+1;
                    coh{numcond}(1:nMD, 1:numFr1) = 0;
                    for t = 1:length(correctTransitionFr)
                        coh{numcond}(currDir,correctTransitionFr(t):end) = par.correctCoh{currCoh}(t);
                    end
                    
                    for t = 1:length(counterTransitionFr)-1
                        si = linspace(par.counterCoh{currCoh}(t),par.counterCoh{currCoh}(t+1),counterTransitionFr(t+1)-counterTransitionFr(t)+1);
                        coh{numcond}(3-currDir,counterTransitionFr(t):counterTransitionFr(t+1)-1) = si(2:end);
                    end
                    
                    trigger{numcond} = [par.CD_DOTS_ON 100+par.correctCoh{currCoh}]; %100+numcond]; % %length(par.correctCoh)*(m-1)];    % what triggers do we send?
                    
                    triggerFr{numcond} = [1 correctTransitionFr(1)];    % in what frame do we send a trigger
                    
                    correctDir(numcond)     = currDir;
                    patchWithMo(numcond)    = currPatch;
                    dotcolor{numcond}       = par.rgbGRAY'*ones(1,size(coh{numcond},2));
                    numperminblock(numcond) = par.condProportn(currCoh);             
                    
                    par.conditiondescrip{numcond} = ['Trigger ' num2str(numcond) ': correct coh ' num2str(par.correctCoh{currCoh}) ' counter coh ' num2str(par.counterCoh{currCoh}) ', motion dir ' num2str(par.motionDir(currDir)) ', ITI ' num2str(par.secs_btw_targs(currITI)) ', patch ' num2str(currPatch)];

                end
                
            end
        end
    end
    disp(['Number of conditions: ' num2str(numcond)])
    
    %%  -------------------------------------------------------------------
    %                   Trial condition randomization
    %   -------------------------------------------------------------------
    minblock = [];
    for currTrial = 1:numcond
        minblock = [minblock ones(1, numperminblock(currTrial))*currTrial];
    end
    
    temp = repmat(minblock,[1,ceil(par.numtargets/size(minblock,2))]);
    temp = temp(:,randperm(size(temp,2))); 	% randomize the columns
    temp(:,par.numtargets+1:end)=[];        % if not evenly distributed, shave off the end
    trialCond = temp;
    
    %%  -------------------------------------------------------------------
    %                   MAIN PRESENTATION LOOP
    %   -------------------------------------------------------------------
  
    %  ------------------   INSTRUCTIONS     ------------------------------
    leftmargin = 0.1;
    Screen('DrawText', window, 'Click left button on mouse for Upward motion.', leftmargin*scresw, 0.2*scresh, 255);
    Screen('DrawText', window, 'Press any button to begin task.', leftmargin*scresw, 0.9*scresh, 255);
    Screen('Flip', window);
    
    % Preallocated things that will be save on a trial by trial basis
    [PTBtrigT, PTBtrig, ClickT, Click, RespLR, RespT, RespTcopy] = deal([]);
    nPTBtrig = 0;
    numResp  = 0;
    
    %  ------------------   START     ------------------------------
    % wait for the participant to click a bottum
    [~, ~, ~, whichButton] = GetClicks(whichScreen,0);
   
    if par.recordEEG, sendtrigger64(par.CD_RESP,port,ioObj), end
    if par.recordEL, Eyelink('Message', 'TASK_START'); end
   
    ClickT(1) = GetSecs;
    Click(1)  = whichButton(1);    % The first response will be the one that sets the task going, after subject reads instructions
    
    Screen('DrawText', window, 'Loading...', 0.35*scresw, 0.5*scresh, 255);
    Screen('Flip', window);
    
    %  ------------------   CREATE DOTS     ------------------------------
    clear PT
    % First make ALL dot stimuli and store:
    dots = cell(2,par.numtargets);   % This will contain dot locations relative to center of screen, in DEGREES
    for currTrial = 1:par.numtargets
        % First, the patch with the motion (the only patch if numPatches=1)
        pwm = patchWithMo(trialCond(currTrial));
        dots{pwm,currTrial}=[]; % every dots cell has three dimensions: dot number x coordinate (x,y) x frame
        numFr = size(coh{trialCond(currTrial)},2);
        % First generate dots at random locations on each frame
        for i=1:numFr
            for d=1:par.numdots
                dots{pwm,currTrial}(d,:,i) = [(rand-0.5)*par.dotpatchsize (rand-0.5)*par.dotpatchsize];
            end
        end
        % then add the coherence by selecting dots to move in certain direction relative to
        % previous frame. A different random set is selected each frame.
        for i=2:numFr
            r = randperm(par.numdots);
            for currDir=1:nMD
                ncd = round(par.numdots*coh{trialCond(currTrial)}(currDir,i)/100);
                randsel = r(1:ncd);
                % for the selected dots, move them in a particular direction
                dots{pwm,currTrial}(randsel,1,i) = dots{pwm,currTrial}(randsel,1,i-1)+cos(par.motionDir(currDir)*pi/180)*par.dotspeed/par.FlickFdots(pwm);         % x-coordinate
                dots{pwm,currTrial}(randsel,2,i) = dots{pwm,currTrial}(randsel,2,i-1)-sin(par.motionDir(currDir)*pi/180)*par.dotspeed/par.FlickFdots(pwm);         % y-coordinate
                r(1:ncd)=[];
            end
            % if it's gone off to the left, wrap it around to the far right
            dots{pwm,currTrial}(find(dots{pwm,currTrial}(:,1,i)<par.dotpatchsize/2),1,i) = dots{pwm,currTrial}(find(dots{pwm,currTrial}(:,1,i)<par.dotpatchsize/2),1,i)+par.dotpatchsize;
            % if it's gone off to the right, wrap it around to the far left
            dots{pwm,currTrial}(find(dots{pwm,currTrial}(:,1,i)>par.dotpatchsize/2),1,i) = dots{pwm,currTrial}(find(dots{pwm,currTrial}(:,1,i)>par.dotpatchsize/2),1,i)-par.dotpatchsize;
        end
        % Finally, go through the dots and get rid of the dots falling outside the
        % circle - put them off the screen.
        for i=1:numFr
            for d=1:par.numdots
                if sqrt(sum(dots{pwm,currTrial}(d,:,i).^2)) > par.dotpatchsize/2
                    dots{pwm,currTrial}(d,:,i) = 2*center/deg2px + 0.01;
                end
            end
        end
        
        PT{currTrial}=[];
        % Make the pulse train:
        PT1 = [];
        for i=1:numFr
            PT1 = [PT1 ; i*ones(nrefON(pwm),1);zeros(nrefOFF(pwm),1)];
        end
        PT{currTrial}(:,pwm) = PT1;
    end
    
    
    %  ------------------   START TRIALS    ------------------------------
    
    % Lead-in:
    Screen('FillRect',window, 255, fixRect);
    Screen('Flip', window);
    WaitSecs(par.leadintime/1000);
    
    % Dot presentation
    ButtonDown = 0;
    for currTrial = 1:par.numtargets
        pwm = patchWithMo(trialCond(currTrial));    % patch with motion (only one if one patch)

        % Dot motion
        trigs_sent = 0;
        for i=1:size(PT{currTrial},1)   % PT is pulse train. 'i' here goes through every refresh of the monitor
            
            for currPatch = 1:par.numPatches
                if PT{currTrial}(i,currPatch)
                    Screen('DrawDots', window, dots{currPatch,currTrial}(:,:,PT{currTrial}(i,currPatch))'*deg2px, par.dotsize,...
                        dotcolor{trialCond(currTrial)}(:,min([PT{currTrial}(i,currPatch) size(dotcolor{trialCond(currTrial)},2)])),...
                        round(center+par.patchloc(currPatch,:).*[1 -1]*deg2px));
                end
            end
            Screen('FillRect',window, 255, fixRect); % put up fixation on top
            
            % Now send triggers when they need to be sent
            trg = find(triggerFr{trialCond(currTrial)} == PT{currTrial}(i,pwm));
            if ~isempty(trg)
                if trg > trigs_sent
                    
                    if par.recordEEG, sendtrigger64(trigger{trialCond(currTrial)}(trg),port,ioObj); end
                    if par.recordEL, Eyelink('Message', ['TRIAL' num2str(currTrial) '_' num2str(trigger{trialCond(currTrial)}(trg))]); end
                    
                    nPTBtrig = nPTBtrig+1;
                    [VBLTimestamp, PTBtrigT(nPTBtrig)] = Screen('Flip', window);
                    PTBtrig(nPTBtrig) = trigger{trialCond(currTrial)}(trg);
                    trigs_sent = trigs_sent+1;
                else
                    Screen('Flip', window);
                end
                checkButton
            else
                Screen('Flip', window);
                checkButton
            end
        end
    end
    
    % And then finish out with half of the first trial (only incoherent bit),
    % so it doesn't end on a target
    currTrial = 1;
    pwm = patchWithMo(trialCond(currTrial));

    % Dot motion
    trigs_sent = 0;
    for i = 1:round(size(PT{currTrial},1)/3)
        
        for currPatch=1:par.numPatches
            if PT{currTrial}(i,currPatch)
                Screen('DrawDots', window, dots{currPatch,currTrial}(:,:,PT{currTrial}(i,currPatch))'*deg2px, par.dotsize, dotcolor{trialCond(currTrial)}(:,min([PT{currTrial}(i,currPatch) size(dotcolor{trialCond(currTrial)},2)])), round(center+par.patchloc(currPatch,:).*[1 -1]*deg2px));
            end
        end
        Screen('FillRect',window, 255, fixRect);
        
        trg = find(triggerFr{trialCond(currTrial)}==PT{currTrial}(i,pwm));
        if ~isempty(trg)
            if trg>trigs_sent
                if par.recordEEG, sendtrigger64(trigger{trialCond(currTrial)}(trg),port,ioObj); end
                if par.recordEL, Eyelink('Message', ['TRIAL' num2str(currTrial) '_' num2str(trigger{trialCond(currTrial)}(trg))]); end
                
                nPTBtrig = nPTBtrig+1;
                [VBLTimestamp, PTBtrigT(nPTBtrig)] = Screen('Flip', window);
                PTBtrig(nPTBtrig) = trigger{trialCond(currTrial)}(trg);
                trigs_sent = trigs_sent+1;
            else
                Screen('Flip', window);
            end
            checkButton
        else
            Screen('Flip', window);
            checkButton
        end
    end
    
    % Feed back
    response_deadline = par.targetDur+0.5;  % in seconds, OFFLINE response deadline for counting performance
    cohmo_trigs = find(PTBtrig>100);
    if length(cohmo_trigs)~=length(trialCond), error('trigger number mismatch!!'); end
    
    numHits = 0; clear RTs acc
    for currTrial=1:length(cohmo_trigs)
        RTs(currTrial) = nan;
        acc(currTrial) = 0;
        stimtime = PTBtrigT(cohmo_trigs(currTrial));
        nextresp = find(RespT>stimtime & RespT<stimtime+response_deadline,1);
        if ~isempty(nextresp)
            RTs(currTrial) = RespT(nextresp) - stimtime;
            numHits = numHits + 1;
            if RespLR(nextresp)==correctDir(trialCond(currTrial))
                acc(currTrial)=1;
            end
        end
    end
    
    disp(['Accuracy: ' num2str(round(mean(acc)*100)) '%'])
    disp([num2str(length(find(isnan(RTs)))) ' Misses'])
    
    %False Alarms
    numFA = length(RespT)- numHits;
    
    txtstart = 0.1*scresw;
    hitRate  = round(mean(acc)*100);
    missRate = round(100*length(find(isnan(RTs)))/par.numtargets);
    
    Screen('DrawText', window, ['On this block you had: '], txtstart, 0.25*scresh, 255);
    Screen('DrawText', window, ['A Hit Rate of ' num2str(round(mean(acc)*100)) '%.'], txtstart, 0.35*scresh, 255);
    Screen('DrawText', window, ['A Miss Rate of ' num2str( round(100*length(find(isnan(RTs)))/par.numtargets) ) '%.'], txtstart, 0.45*scresh, 255);
    Screen('DrawText', window, [num2str(numFA) ' False Alarms.'], txtstart, 0.55*scresh, 255);
    
    Screen('DrawText', window, 'Click to Exit', txtstart, 0.70*scresh, 255);
    Screen('Flip', window);
    GetClicks(whichScreen,0);
    
    save([par.runID],'ClickT','Click','nPTBtrig','PTBtrigT','PTBtrig',...
        'RespT','RespLR','trialCond','coh','dotcolor','trigger','triggerFr',...
        'correctDir','patchWithMo','par', 'hitRate', 'missRate', 'numFA', 'PT') % NOTE THAT FOR ORIGNAL EXPERIMENT PT WAS UNFORTUANTLY NOT SAVED!!!
    save parDC par
    sca;
    ShowCursor;
    
catch ME
    sca
    ShowCursor;
    rethrow(ME)
end

