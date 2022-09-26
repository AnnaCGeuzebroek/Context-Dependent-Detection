% Continuous Dots paradigm, with 3 interleaved frames of dot motion

% reset(RandStream.getDefaultStream,sum(100*clock)); % to avoid randomisation problem
% reset(RandStream.getGlobalStream,sum(100*clock)); % to avoid randomisation problem

%% ----------- INITIALIZE PARAMETERS    -----------------------------------
rng('shuffle');  % to avoid randomisation problem
PsychDebugWindowConfiguration % use for debugging (give transperant screen)

clear
close all

SITE = 'T';     % T = TCD, C = City college!
commandwindow;

if SITE=='C' || SITE=='E'
    TheUsualParamsCRT_Dell_lores      % this script defines some useful parameters of the monitor, trigger port, etc common to all experiments
elseif SITE == 'T'
    TheUsualParamsCRT_TCD
end
load parDC

%%%%%%%%%% IMPORTANT SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note many of the time values in ms won't work out exactly - depends on
% refresh rate & flicker - it'll be the closest it can be, should be
% checked in EEG
par.videoFrate = 100;   % Monitor refresh rate
par.FlickFdots = 100;	% Flicker frequency in Hz for the dots, set to 60 so that there's no flicker.
par.skipframes = 3;     % How many frames to skip in interleaved dot motion.

par.numtargets  = 42; %210
par.BGcolor     = 0;
par.cohLevels   = 25; % coherence levels (%)
par.rampTime    = 0.2; % GL: in sec
par.targetDur   = 1.88;    % coherent motion duration in seconds
par.dotspeed    = 6;       % coherent motion in degrees per second
par.dotspeed    = par.dotspeed*par.skipframes;
par.dotsize     = 4;   % in pixels
par.numdots     = 75;  % this will be the number in a square; the ones outside the circle will be taken away
% par.patchloc = [0 0];
par.patchloc    = [0 0]; % patch location coordinates [x y] in degrees relative to center of screen
par.motionDir   = 90;    % in degrees relative to positive x-axis (0=rightward), GL: 0 = right, 180 = left. 90 270
par.dotpatchsize = 8;   % degrees of visual angle - diameter of dot patch

par.leadintime      = 1000;
par.secs_btw_targs  = [2 4 6]; % ITI in sec
par.gap_start_secs  = 0.175;   % gap timing
par.gap_end_secs    = 0.375;

par.rgbGRAY = [221 221 221];

% par.EEGtrigCodes = [25 35 50 70; 26 36 51 71];
par.MinTimeBtwClicks = 0.5; % In sec. When polling for mouse clicks, don't want to be counting the same click any more than once.

par.useEL     = 0; %Use Eyelink? (1=yes, 0=no)
par.recordEEG = 0;
par.runID     = 'Test';
par.prac      = 0;

dlg_title = 'Dot Motion Task';
while 1
    prompt  = {'Enter SUBJECT/RUN/TASK IDENTIFIER:','EEG? (1=yes, 0=no)','Practice? (1=yes,0=no)'};
    def     = {par.runID,num2str(par.recordEEG), num2str(par.prac)};
    answer  = inputdlg(prompt,dlg_title,1,def);
    par.runID     = answer{1};
    par.recordEEG = str2double(answer{2});
    par.prac      = str2double(answer{3});
    
    if exist([par.runID '.mat'],'file')
        dlg_title = [par.runID '.mat EXISTS ALREADY - CHOOSE ANOTHER, OR DELETE THAT ONE IF IT IS RUBBISH'];
    else
        break;
    end
end

if par.prac == 1
    par.gaps      = 1;
    par.cohLevels = 50; % coherence levels (%)
else
    par.gaps = [1,2];
end

% Set up for triggers
if par.recordEEG
    if SITE == 'C'
        % USB port (posing as Serial Port) for triggers
        [port, errmsg] = IOPort('OpenSerialPort', 'COM3','BaudRate=115200');
        IOPort('Write', port, uint8([setpulsedur 2 0 0 0]))   % pulse width given by last 4 numbers (each a byte, little-endian)
    elseif SITE == 'E'
        port = hex2dec('1030');
        lptwrite(port,0);
    elseif SITE == 'T'
        % Parallel Port for triggers - set to zero to start
        port = hex2dec('2010');
        lptwrite(port,0);
    end
end

if par.useEL, ELCalibrateDialog, end

% if abs(hz-par.videoFrate)>1
%     error(['The monitor is NOT SET to the desired frame rate of ' num2str(par.videoFrate) ' Hz. Change it.'])
% end

window = Screen('OpenWindow', whichScreen, par.BGcolor);

%% ------------------- CODES AND TRIAL SEQUENCE ---------------------------
% trigger codes - can only use these 15: [1 4 5 8 9 12 13 16 17 20 21 24 25 28 29]  (don't ask!)
par.CD_RESP     = 1;
par.CD_DOTS_ON  = 5;
par.CD_BUTTONS  = [12 13];   % left and right mouse
par.CD_COHMOTION_ON = [8 9];

%%%%%  coherence, motion direction waveforms
nMD = length(par.motionDir);

% par.FlickTdots = round(par.videoFrate./par.FlickFdots); % In number of video frames
% nrefON = floor(par.FlickTdots/2);                       % number of refreshes where the dots are ON
% nrefOFF = par.FlickTdots-nrefON;                        % number of refreshes where the dots are OFF

par.FlickTdots = round(par.videoFrate/par.FlickFdots); % In number of video frames
nrefON  = par.FlickTdots;                              % number of refreshes where the dots are ON
nrefOFF = 0;                                           % number of refreshes where the dots are OFF

numcond = 0;
clear conditiondescrip coh trigger triggerFr correctDir

for g = 1:length(par.gaps)
    for c = 1:length(par.cohLevels)
        for o = 1:length(par.secs_btw_targs)
            
            numFr1 = round((par.targetDur+par.secs_btw_targs(o)).*par.FlickFdots);
            cohMotionOnFr1 = round(par.secs_btw_targs(o)*par.FlickFdots)+1;
            numCohMotionFr = numFr1-cohMotionOnFr1;
            par.numRampFr = ceil(par.rampTime/(1/par.videoFrate)); % approx 100ms at 60Hz is 6 frames.
            fullCohMotionOnFr1 = cohMotionOnFr1+par.numRampFr;
            
            t_frames = -cohMotionOnFr1+1:1:par.numRampFr;
            unitstep = t_frames>=0;
            ramp = (t_frames.*unitstep)/par.numRampFr;
            quad = (t_frames.^2.*unitstep)/(par.numRampFr^2);
            quad_down = -quad;
            quad_down(unitstep) = quad_down(unitstep)+2;
            ramp_down = -ramp;
            ramp_down(unitstep) = ramp_down(unitstep)+2;
            
            % figure, plot(t_frames,[unitstep;ramp;ramp_down]), set(gca,'xlim',[-50,50],'ylim',[-0.1,2.5])
            
            unitstep = ceil([unitstep,ones(1,numFr1-fullCohMotionOnFr1)]*par.cohLevels(c));
            ramp = ceil([ramp,ones(1,numFr1-fullCohMotionOnFr1)]*par.cohLevels(c));
            quad = ceil([quad,ones(1,numFr1-fullCohMotionOnFr1)]*par.cohLevels(c));
            quad_down = ceil([quad_down,ones(1,numFr1-fullCohMotionOnFr1)]*par.cohLevels(c));
            ramp_down = ceil([ramp_down,ones(1,numFr1-fullCohMotionOnFr1)]*par.cohLevels(c));
            t_frames = -cohMotionOnFr1+1:1:numCohMotionFr;
            time = t_frames*1000/par.videoFrate;
            
            %{
    %             coh_plots = [unitstep;ramp;ramp_down];
    %             figure, clear h
    %             for kk = 1:3
    %                 h(kk) = plot(time,coh_plots(kk,:),'LineWidth',2);  hold on
    %             end
    %             set(gca,'xlim',[-700,700],'ylim',[-2,50],'FontSize',12);
    %             legend(h,{'Step','Ramp Up','Ramp Down'},'Location','NorthWest')
    %             xlabel('Time from coherent motion onset (ms)')
    %             ylabel('Coherence (%)')

    %             figure, plot(t_frames,[unitstep;quad;quad_down]), set(gca,'xlim',[-50,50],'ylim',[-2,100])
            %}
            
            par.rampers = [unitstep;ramp;ramp_down];
            for ramper = 1%:3 % This just makes all trials a step function
                for m=1:nMD
                    numcond = numcond+1;
                    
                    coh{numcond}(1:nMD,1:numFr1) = 0;
                    coh{numcond}(m,1:numFr1) = par.rampers(ramper,:);
                    
                    if par.gaps(g)==2
                        coh{numcond}((par.gap_start_secs+par.secs_btw_targs(o))*par.videoFrate:(par.gap_end_secs+par.secs_btw_targs(o))*par.videoFrate)=0;
                    end
                    
                    if SITE=='T'|SITE=='E'
                        trigger{numcond} = [par.CD_DOTS_ON par.CD_COHMOTION_ON(g)];
                    elseif SITE=='C'
                        trigger{numcond} = [par.CD_DOTS_ON par.CD_COHMOTION_ON];
                    end
                    triggerFr{numcond} = [1 cohMotionOnFr1];
                    
                    correctDir(numcond) = m;
                    rampcond(numcond) = ramper;
                    conditiondescrip{numcond} = ['Trigger ' num2str(numcond) ': coherence ' num2str(par.cohLevels(c)) ', motion dir ' num2str(par.motionDir(m)) ', ITI ' num2str(par.secs_btw_targs(o)) ', ramp ' num2str(ramper)];
                    
                end
            end
        end
    end
end

% No colored cues implemented for continuous dots...
for n=1:numcond
    dotcolor{n} = par.rgbGRAY'*ones(1,size(coh{n},2));
    numperminblock(n) = 1;
end

disp(['Number of conditions: ' num2str(numcond)])

% Trial condition randomization:
minblock = [];
for n=1:numcond
    minblock = [minblock ones(1,numperminblock(n))*n];
end
temp = repmat(minblock,[1,ceil(par.numtargets/size(minblock,2))]);

% Now RANDOMIZE the sequence of trials
temp = temp(:,randperm(size(temp,2)));      % jumble the columns in temp
% if numtrials is not evenly divisible by the minimum block length, shave off some trials from the end
temp(:,par.numtargets+1:end)=[];
trialCond = temp;

if par.useEL
    ELsetupCalib
end
Screen('FillRect',window, par.BGcolor); % screen blank
% window = Screen('OpenWindow', whichScreen, par.BGcolor);

%% ------------------------- START TASK -----------------------------------
% Instructions:
leftmargin = 0.1;
Screen('DrawText', window, 'Click left button with right hand for UPWARD motion.', leftmargin*scres(1), 0.2*scres(2), 255);
Screen('DrawText', window, 'Press any button to begin task.', leftmargin*scres(1), 0.9*scres(2), 255);
Screen('Flip', window);
HideCursor;
% Things that we'll save on a trial by trial basis
clear PTBtrigT PTBtrig ClickT Click RespLR
RespT=[];
nPTBtrig=0;
numResp=1;

% Waits for the user to press a button.
[clicks,x,y,whichButton] = GetClicks(whichScreen,0);
if par.recordEEG, sendtrigger(par.CD_RESP,port,SITE,0), end
if par.useEL
    Eyelink('Message', ['TASK_START']);
end
ClickT(1) = GetSecs;
Click(1)  = whichButton(1);    % The first response will be the one that sets the task going, after subject reads instructions

%%%%%%%%%%%%%%%%%%%% START TRIALS

% initial lead-in:
% Screen('FillRect',window, 255, fixRect);
% Screen('Flip', window);

Screen('DrawText', window, 'Loading...', 0.35*scres(1), 0.5*scres(2), 255);
Screen('Flip', window);
WaitSecs(par.leadintime/1000);

clear PT
% First make ALL dot stimuli and store:
dots = cell(2,par.numtargets);   % This will contain dot locations relative to center of screen, in DEGREES
for n=1:par.numtargets
    
    dots{n}=[];
    numFr = round(size(coh{trialCond(n)},2)*par.FlickTdots/par.FlickTdots);
    % First generate dots at random locations on each frame
    for i=1:numFr
        for d=1:par.numdots
            dots{n}(d,:,i) = [(rand-0.5)*par.dotpatchsize (rand-0.5)*par.dotpatchsize];
        end
    end
    % if this is the patch with dots, make coherent motion:
    
    % then add the coherence by selecting dots to move in certain direction relative to
    % previous frame. A different random set is selected each frame.
    for i=par.skipframes+1:numFr
        r = randperm(par.numdots);
        for m=1:nMD
            ncd = round(par.numdots*coh{trialCond(n)}(m,i)/100);
            randsel = r(1:ncd);
            % for the selected dots, move them in a particular direction
            dots{n}(randsel,1,i) = dots{n}(randsel,1,i-par.skipframes)+cos(par.motionDir(m)*pi/180)*par.dotspeed/par.FlickFdots;         % x-coordinate
            dots{n}(randsel,2,i) = dots{n}(randsel,2,i-par.skipframes)-sin(par.motionDir(m)*pi/180)*par.dotspeed/par.FlickFdots;         % y-coordinate
            r(1:ncd)=[];
        end
        % if it's gone off to the left, wrap it around to the far right
        dots{n}(find(dots{n}(:,1,i)<par.dotpatchsize/2),1,i) = dots{n}(find(dots{n}(:,1,i)<par.dotpatchsize/2),1,i)+par.dotpatchsize;
        % if it's gone off to the right, wrap it around to the far left
        dots{n}(find(dots{n}(:,1,i)>par.dotpatchsize/2),1,i) = dots{n}(find(dots{n}(:,1,i)>par.dotpatchsize/2),1,i)-par.dotpatchsize;
        % if it's gone off to the left, wrap it around to the far right
        dots{n}(find(dots{n}(:,2,i)<par.dotpatchsize/2),2,i) = dots{n}(find(dots{n}(:,2,i)<par.dotpatchsize/2),2,i)+par.dotpatchsize;
        % if it's gone off to the right, wrap it around to the far left
        dots{n}(find(dots{n}(:,2,i)>par.dotpatchsize/2),2,i) = dots{n}(find(dots{n}(:,2,i)>par.dotpatchsize/2),2,i)-par.dotpatchsize;
    end
    % Finally, go through the dots and get rid of the dots falling outside the
    % circle - put them off the screen.
    for i=1:numFr
        for d=1:par.numdots
            if sqrt(sum(dots{n}(d,:,i).^2)) > par.dotpatchsize/2
                dots{n}(d,:,i) = 2*center/deg2px + 0.01;
            end
        end
    end
    PT1 = [];
    for i=1:numFr
        PT1 = [PT1 ; i*ones(nrefON,1);zeros(nrefOFF,1)];
    end
    PT{n} = PT1;
end

% initial lead-in:
Screen('FillRect',window, 255, fixRect);
Screen('Flip', window);
WaitSecs(3);

% START STIMULATION
portUP=0; lastTTL=0; ButtonDown=0;

tic
for n=1%:par.numtargets
    disp(['Condition ',num2str(trialCond(n))])
    % DOT MOTION
    trigs_sent = 0;
    for i=1:size(PT{n},1)
        if par.recordEEG, if SITE=='T'|SITE=='E', if portUP & GetSecs-lastTTL>0.01, lptwrite(port,0); portUP=0; end, end, end
        
        if PT{n}(i)
            Screen('DrawDots', window, dots{n}(:,:,PT{n}(i))'*deg2px, par.dotsize, dotcolor{trialCond(n)}(:,min([PT{n}(i) size(dotcolor{trialCond(n)},2)])), round(center+par.patchloc.*[1 -1]*deg2px));
        end
        
        Screen('FillRect',window, 255, fixRect);
        
        trg = find(triggerFr{trialCond(n)}==PT{n}(i));
        if ~isempty(trg)
            if trg>trigs_sent
                if par.recordEEG, sendtrigger(trigger{trialCond(n)}(trg),port,SITE,1); portUP=1; lastTTL=GetSecs; end
                if par.useEL, Eyelink('Message', ['TRIAL' num2str(n) '_' num2str(trigger{trialCond(n)}(trg))]); end
                nPTBtrig = nPTBtrig+1;
                [VBLTimestamp PTBtrigT(nPTBtrig)] = Screen('Flip', window);
                PTBtrig(nPTBtrig) = trigger{trialCond(n)}(trg);
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
n=1;
% DOT MOTION
trigs_sent = 0;
for i=1:round(size(PT{n},1))
    if par.recordEEG 
        if SITE=='T' || SITE=='E'
            if portUP && GetSecs-lastTTL > 0.01
                lptwrite(port,0); portUP=0; 
            end
        end
    end
    
    if PT{n}(i)
        Screen('DrawDots', window, dots{n}(:,:,PT{n}(i))'*deg2px, par.dotsize, dotcolor{trialCond(n)}(:,min([PT{n}(i) size(dotcolor{trialCond(n)},2)])), round(center+par.patchloc.*[1 -1]*deg2px));
    end
    
    Screen('FillRect',window, 255, fixRect);
    
    trg = find(triggerFr{trialCond(n)}==PT{n}(i));
    %{
    %     if ~isempty(trg)
    %         if trg>trigs_sent
    %             if par.recordEEG, sendtrigger(trigger{trialCond(n)}(trg),port,SITE,1); portUP=1; lastTTL=GetSecs; end
    %             if par.useEL, Eyelink('Message', ['TRIAL' num2str(n) '_' num2str(trigger{trialCond(n)}(trg))]); end
    %             nPTBtrig = nPTBtrig+1;
    %             [VBLTimestamp PTBtrigT(nPTBtrig)] = Screen('Flip', window);
    %             PTBtrig(nPTBtrig) = trigger{trialCond(n)}(trg);
    %             trigs_sent = trigs_sent+1;
    %         else
    %             Screen('Flip', window);
    %         end
    %         checkButton
    %     else
    %}
    Screen('Flip', window);
    %     checkButton
    %     end
end

toc

if par.useEL
    Eyelink('CloseFile');
    Eyelink('ReceiveFile',[answer{1},'.edf']);
end

% FEEDBACK %DN:added Feedback here, took code from Redmond's 'reverse'scrip
response_deadline = par.targetDur+0.2;  % in seconds, OFFLINE response deadline for counting performance
if SITE=='T'
    cohmo_trigs=find(PTBtrig==par.CD_COHMOTION_ON(1)|PTBtrig==par.CD_COHMOTION_ON(2));
elseif SITE=='C'
    cohmo_trigs = find(PTBtrig==par.CD_COHMOTION_ON);
end
if length(cohmo_trigs) ~=length(trialCond), error('trigger number mismatch!!'); end %DN: Made length(cohmo_trigs)-1 instead of justlength(cohmo_trigs) to get rid of that last half trial maybe delete this change?
clear RTs acc condo
for n=1:length(cohmo_trigs) %DN: Made length(cohmo_trigs)-1 instead of justlength(cohmo_trigs) to get rid of that last trial maybe delete this change?
    condo(n) = PTBtrig(cohmo_trigs(n));
    RTs(n) = nan;
    acc(n) = 0;
    stimtime = PTBtrigT(cohmo_trigs(n));
    nextresp = find(RespT>stimtime & RespT<stimtime+response_deadline,1);
    if ~isempty(nextresp)
        RTs(n) = RespT(nextresp) - stimtime;
        if RespLR(nextresp)==correctDir(trialCond(n))
            acc(n)=1;
        end
    end
end


% figure;
% hist(RTs*1000,[0:100:2500]), hold on
% title(['MEDIAN RT: ' num2str(nanmedian(RTs)*1000)])

ramp_conds(1,:) = [101,104,107];
ramp_conds(2,:) = [102,105,108];
ramp_conds(3,:) = [103,106,109];

for i = 1:3
    disp(['Ramp cond ',num2str(i),': ',num2str(nanmedian(RTs(find(ismember(condo,ramp_conds(i,:))))*1000))])
end

disp(['Accuracy: ' num2str(round(mean(acc)*100)) '%'])
disp([num2str(length(find(isnan(RTs)))) ' Misses'])
txtstart = 0.1*scres(1);
Screen('DrawText', window, ['On this block you got ' num2str(round(mean(acc)*100)) '% correct.'], txtstart, 0.25*scres(2), 255);
Screen('DrawText', window, ['You MISSED ' num2str(length(find(isnan(RTs)))) ' trials out of the ' num2str(length(RTs)) '.'], txtstart, 0.35*scres(2), 255);
Screen('DrawText', window, ['Click to Exit'], txtstart, 0.55*scres(2), 255);
Screen('Flip', window);
GetClicks(whichScreen,0);

save([par.runID],'ClickT','Click','nPTBtrig','PTBtrigT','PTBtrig','RespT','RespLR','trialCond','coh','dotcolor','trigger','triggerFr','correctDir','par','conditiondescrip','RTs','acc','dots')

save parDC2 par

sca

return
%% Calc overall RTs etc
RT_block=[]; condo_block=[]; acc_block=[];
for i = 1:8
    load(['CJ_',num2str(i),'.mat'])
    clear RTs acc condo
    for n=1:(length(cohmo_trigs)-1) %DN: Made length(cohmo_trigs)-1 instead of justlength(cohmo_trigs) to get rid of that last trial maybe delete this change?
        condo(n) = PTBtrig(cohmo_trigs(n));
        RTs(n) = nan;
        acc(n) = 0;
        stimtime = PTBtrigT(cohmo_trigs(n));
        nextresp = find(RespT>stimtime & RespT<stimtime+response_deadline,1);
        if ~isempty(nextresp)
            RTs(n) = RespT(nextresp) - stimtime;
            if RespLR(nextresp)==correctDir(trialCond(n))
                acc(n)=1;
            end
        end
    end
    RT_block = [RT_block,RTs*1000];
    condo_block = [condo_block,condo];
    acc_block = [acc_block,acc];
end
ramp_conds(1,:) = [101,104,107];
ramp_conds(2,:) = [102,105,108];
ramp_conds(3,:) = [103,106,109];

figure, clear h
for i = 1:3
    [N(i,:),edges] = histcounts(RT_block(find(ismember(condo_block,ramp_conds(i,:)))),10);
    h(i) = plot(edges(1:end-1),N(i,:)), hold on
    disp(['Ramp cond ',num2str(i),': ',num2str(nanmedian(RT_block(find(ismember(condo_block,ramp_conds(i,:))))))])
end
legend(h,{'Step','Ramp Up','Ramp Down'})

