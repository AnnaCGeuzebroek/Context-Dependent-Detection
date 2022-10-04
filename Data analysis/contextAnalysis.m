%% Data preprocessing and EEG waveform analysis.
% Data collected Hannah Craddock for her Master's thesis, 2016:
%   'Electrophysiological and behavioural indices of decision criterion
%   adjustments across context of weak and strong evidence.'
%
% Current code is written by A.C. Geuzebroek and part of:
%   Geuzebroek AC, Craddock H, O’Connell RG, & Kelly SP (2022).
%   Balancing true and false detection of intermittent sensory targets
%   by adjusting the inputs to the evidence accumulation process 
%   (https://biorxiv.org/cgi/content/short/2022.09.01.505650v1)
%
% In short, participants were asked to continous monitor a cloud of randomly moving
% dots for intermittered target defined as upwards coherent moving dots.
% Afterwhich they were asked to response as fast and accurate as possible.
% Difficulty context was manipulated with:
%   1) Hard  (25% motion coherence)
%   2) Easy  (70% motion coherence)
%   3) Mixed (25% and 70% motion coherence with equal probabiltiy).
% 
% Code requires you to download the preprocessed data at:
% https://osf.io/yjvku/?view_only=7ed5aee5d09a4d5ca13de1ba169b0588
%
% Add this data in the follow way:
%   InputFolder     = BASEFOLDER\Data\Processed\
%       Participant folder = BASEFOLDER\Data\Processed\PPNAME\ (note never use initials!)
%           EEG folder           = BASEFOLDER\Data\Processed\EEG data\PPNAME\
%
%
% Depending on:
%   1) EEGLAB (including Biosig extention).
%   2) CSD toolbox and lay-out.
%   3) findNoisyChannels
%   4) Brewermap (to get colors for plot. Can be easily replaced with just choosing colours)
%   5) panels
%   6) costum-made code dataAnalysis(object to access all the functions)

clear
clc
clear global

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------     Input parameters    -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set all paths.
if strcmp(computer, 'PCWIN64')
    %error('Please set your data folder')
    inputFolder = 'C:\Users\eng\Documents\2. Postdoc_Kelly\Project 5 - Neural correlates of static and dynamic decision-bound adjustments\Project 5.1 - Signatures of bound adjustment\Data';
    % inputFolder = 'YOURFOLDER\Data';    % add base folder of data. 
    addpath(genpath('C:\Users\eng\Downloads\Context-Dependent-Detection-main\dataAnalysis'))
    
    %error('Please add the path to EEGLAB!')

    if ~exist('eeglab', 'file')
        run('YOURPATHTOEEGLAB\EEGLAB\eeglab'); close all;
    end
end

% In some behavioral task it might be intresting to look at time on task
% (Within-block effects) or time doing the experiments
blockRandomized = 1; % yes: easy/hard/mixed
numRTBins  = 2; TimeOnTask = 0; TimeOnExperiment = 0;

% Set conditions, here this is:
%   Context (Fixed vs. Mixed)
%   Coherence level (25 vs. 70)

% IMPORTANT condNames should refer to the way that it was saved in the
% experiment. This can be checked by loading the trial information file and
% see how it will be set in the workspace. As the presentation code was written
% with other analysis code in might, it requires a quick for-loop to actually add 
% condition and coherence (see line TODO).

conditions{1}  = [1 2]; % easy difficult mixed
condNames{1}   = 'condition';
conditions{2}  = [25 70];
condNames{2}   = 'coherence';

% preset the signal processing object.
parameters  = dataAnalysis.getInstance();

parameters.numTrials  = 24; % number of trials per block.
parameters.numBlocks  = 12; % number of trials per block --> USUALLY SET BY PER INDIVUDAL in obj.getInformation

% set a couple of standard parameters within this object.
parameters.system = 'bioSemi';   % Object will choice the approtiated loading function of EEGLAB
parameters.analysisEEG      = 1; % Analyse EEG (obviously we want to look at this)

parameters.analysisEyelink  = 0;
parameters.analysisEOG      = 1; % Use additionally EOG electrodes and Frontal electrodes
                                 % to look for blinks and exclude trials with blinks (We do this to prevent
                                 % any possible effects of the blink on the sensory perception)
parameters.analysisEMG      = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------     GET EXPERIMENT INFORMATION    ------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function looping through the data files of each participants. Identifying
% the behavioural information. including numBlocks, numBlocks per
% condition, numStimulus in each conditions and their specific condition (
% here initialy all stimulus in each block are 1 condition.) However, it is
% possible that the conditions information are in the EEG files. So here
% simply put all the par and other parameters of each participant away.

prepData = 0;
parameters.getInformation(prepData, inputFolder, conditions, condNames, blockRandomized);
clear cond* *Folder block*

% As the inter-trial-interval length is an interesting and crucial part of the
% this experiment, we straigth away add ITI as an variable. 
parameters.conditions{3}   = [2 4 6 8];
parameters.condNames{3}    = 'ITI';

% NOTE below was needed to created the right trialMatrix construction, usually,
% this can be done simply in epoching through setting the right condNames before
% calling parameters.getInformation. If designing or writing your own presenation
% code this can be easily prevented by having a save parameters refering to the 
% conditions structure. 
%{
for indPP = 1:length(parameters.ppNames)
    for indBlock = 1:parameters.numBlocks(indPP)
        
        tmpITI = zeros(size(parameters.experiment{indPP}{indBlock}.trialCond));

        if length(parameters.experiment{indPP}{indBlock}.par.correctCoh) == 1
            parameters.experiment{indPP}{indBlock}.coherence  = parameters.experiment{indPP}{indBlock}.par.correctCoh{:};
            parameters.experiment{indPP}{indBlock}.condition  = 1;
            
            for indITI = 1:length(parameters.conditions{3})
                tmpITI(parameters.experiment{indPP}{indBlock}.trialCond == indITI) = parameters.conditions{3}(indITI);
            end
        else
            tmpCoh = zeros(size(parameters.experiment{indPP}{indBlock}.trialCond));
            tmpCoh(ismember(parameters.experiment{indPP}{indBlock}.trialCond,1:2:8)) = 70;
            tmpCoh(ismember(parameters.experiment{indPP}{indBlock}.trialCond,2:2:8)) = 25;
            
            parameters.experiment{indPP}{indBlock}.coherence  = tmpCoh;
            parameters.experiment{indPP}{indBlock}.condition  = 2;
            indexITI = 1:2:length(parameters.conditions{3})*2;
            for indITI = 1:length(indexITI)
                tmpITI(parameters.experiment{indPP}{indBlock}.trialCond == indexITI(indITI))  = parameters.conditions{3}(indITI);
                tmpITI(parameters.experiment{indPP}{indBlock}.trialCond == indexITI(indITI)+1) = parameters.conditions{3}(indITI);
            end
        end
        parameters.experiment{indPP}{indBlock}.ITI = tmpITI;
    end
end
%}

% Set everything for the stimulus presentation. 
% Assuming that parameters are the same for all participants and 
% within blocks.
parameters.DetectOrDisc = 0;         % 0 as for Detection task

parameters.stim.refreshRate  = 60;   % or parameters.experiment{1}{1}.par.videoFrate;	      % in Hz.,  Screens refreshRate
parameters.stim.duration     = 1; % or parameters.experiment{1}{1}.par.targetDur.*1000;    % in ms., duration of the stimulus (NOTE CHANGE TO SEC.)
parameters.stim.freqSSVEP    = 15;   % in Hz (SSVEP), is however not task related 
                                     % and therefore only used to determine window size e.g.
                                     % have whole cycles presents.
parameters.stim.epoch        = [-1.8 2.5];  % in ms., duration of the stimulus (NOTE CHANGE TO SEC.)

parameters.stim.timing     = 'past';      
parameters.stim.RTdeadLine = [-1.8 1.65]; % As it is past, for new codes this can be set to 0, now as the data is already preprocessed
parameters.stim.RTCutOff   = 0.25;        % only include reaction times after the gap offset.

parameters.stim.FA         = 1;
parameters.stim.FACutOff   = [];

parameters.stim.lengthITI  = [2 4 6 8];
parameters.stim.namesITI   = 'ITI';

%% %%%%%%%%%%%%%%%%%%%% Pre-set figure parameters    %%%%%%%%%%%%%%%%%%%%%%
% get your color scheme set here
ColoursCond12  = [128 205 193; 1 133 113; 223 194 125; 166 97 26]/255;  % Context by Motion coherence
ColoursCond3 =  brewermap(4,'Greys'); ColoursCond3 = ColoursCond3-0.16; % ITI

parameters.figLayOut.letterSize  = 9;
parameters.figLayOut.letterType  = 'Arial';
parameters.figLayOut.lineWidth   = 1.2;
parameters.figLayOut.lineType    = {'-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-'};
parameters.figLayOut.legNames{1} = {'Constant', 'Mixed'};
parameters.figLayOut.legTitle{1} = 'Context';

parameters.figLayOut.legNames{2} = {'Weak', 'Strong'}; % shift as they are 'order' by strength
parameters.figLayOut.legTitle{2} = 'Evidence strength';

parameters.figLayOut.legNames{3} = {'2', '4', '6', '8'};
parameters.figLayOut.legTitle{3} = 'ITI duration';

parameters.figLayOut.legends = {'Weak Context', 'Strong Context', 'Mixed (Weak)', 'Mixed (Strong)'};

parameters.figLayOut.plotCI      = 0.05; % get shaded area of not.
parameters.figLayOut.removeInter = 1;
parameters.figLayOut.saveDim     = [5 11];
parameters.figLayOut.colours     = ColoursCond12;

% these range are used for the initial plots to explore the full range.
% We later zoomed in to the new ranges to get a better view. 
parameters.figLayOut.targetLim   = [-0.2 0:0.5:1];
parameters.figLayOut.responseLim = [-0.4 0 0.2];
parameters.figLayOut.plotRT      = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------    EEG PROCESSING PARAMETERS    -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All parameters needed to set-up the pre-processing. We use a most
% simplest pre-processing protocol with the following steps:
%   1) detrending for all (1) or only for the channels that it would
%   benefit (2).
%   2) Low-pass filtered using a filter designed by Simon. (see obj.setFilters)
%   3) High-pass filtered using filtfilt.
%   4) Additionally, the findNoisyChannels within AutoMagic.

%---------------     set pre-processing parameters   ----------------------
parameters.eeg.NumberOfChannels = 128;
parameters.eeg.SampleRate   	= 512;	% EEG sample rate

% load in a structure 'chanlocs' containing the x,y,z locations of each of the 128 scalp channels in the cap.
load chanlocsBioSemi128;
elab = strings(parameters.eeg.NumberOfChannels,1);
for indChan = 1:parameters.eeg.NumberOfChannels
    elab(indChan) = chanlocs(indChan).labels;
end
elab(end+1)='FCz'; clear ind*

parameters.eeg.ChannelsName   	= elab; 
parameters.eeg.chanlocs         = chanlocs;	% electrode locations
parameters.eeg.transChanlocs	= load('coordsCSDBioSemi128.mat');	% electrode locations
clear elab chanlocs

parameters.eeg.ampSaturation	= 3276;	% What's the maximum limiting amplitude of the amplifiers, to detect saturation?

parameters.eeg.postDCC  = 700; % how many samples after DCC is EEG rubbish?
parameters.eeg.applyCSD = 1;

parameters.eeg.applyFindDCcorrect   = 0;

parameters.eeg.badChannel   = 1;	% set on one if you want to apply bad channel selection and interpolation before epoching.
parameters.eeg.applyLPF 	= 1;    % attention here that the lowpass filter is different in New York as in Ireland.
parameters.eeg.cuttoffLPF	= 38;	% 38 cut-off threshold for low-pass filter
                                    % low-pass filter will be different as NY has 60 hz electricity rate.
parameters.eeg.LHamming     = 77;   % 77
parameters.eeg.applyHPF     = 0;    % look at raw data if is nescessary.
parameters.eeg.simonHPF     = 0;    % RODO CHECK THIS OUT??
parameters.eeg.HPFcutoff    = 0;
parameters.eeg.applyDetrend	= 1;    % option 1(all)/ 2(only that will benefit from it)

parameters.eeg.artifactThres = 80;
parameters.eeg.artifactEOG   = 200;
parameters.eeg.intArtifact   = 1;

parameters.eeg.channelVEOG   = 1:2;
parameters.eeg.channelHEOG   = 3:4;

%--------------- set triggers and EPOCHS definition -----------------------
% The trigger codes we used for the different events were as follows:
parameters.triggers.start       = 1;            % Start of the trial
parameters.triggers.stimulusOFF = 5;            % Target off
parameters.triggers.stimulusON  = [170 125];    % 170 for 70% and 125 for 25%
parameters.triggers.response    = [12 13];      % left vs. right

% EPOCH DEFINITION
parameters.eeg.epochLock = parameters.triggers.stimulusON;
parameters.eeg.baseline  = 0 + [-1 0]*1/parameters.stim.freqSSVEP;  

% Initially we set a rather larger epoch based on the stimulus definitition
% 'parameters.stim.epoch'. This is used to cut out a larger epoch to apply
% all the ERP pre-processing on. Afterwards smaller epochs are created to
% get the target- and response-locked epochs. This is usefull especially
% for the ERD (event related desynchronization) and the SSVEP later
% onwards.

% these range are used for the initial plots to explore the full range.
% We later zoomed in to the new ranges to get a better view. 

parameters.eeg.targetEpoch   = -ceil(0.5/(1/parameters.eeg.SampleRate))*(1/parameters.eeg.SampleRate):1/parameters.eeg.SampleRate:ceil(parameters.stim.RTdeadLine(2)/(1/parameters.eeg.SampleRate))*(1/parameters.eeg.SampleRate);
parameters.eeg.responseEpoch = -ceil(0.8/(1/parameters.eeg.SampleRate))*(1/parameters.eeg.SampleRate):1/parameters.eeg.SampleRate:ceil(0.4/(1/parameters.eeg.SampleRate))*(1/parameters.eeg.SampleRate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------    Perform eeg prep  -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE usually this needs to be run to get the preprocessed data in order
% to actually created the epoched data. 
% parameters.applyPreprocessing;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------- Reset behaviour data    -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE below was needed to created the right trialMatrix construction, this
% has already been preformed in the epoched data
%{
for indPP = 1:length(parameters.ppNames)
    for indBlock = 1:parameters.numBlocks(indPP)
        Responses = find(ismember(parameters.event{indPP}{indBlock}.Number, parameters.triggers.response));
        StimOn    = find(ismember(parameters.event{indPP}{indBlock}.Number, parameters.triggers.stimulusON));
        if length(StimOn) > parameters.numTrials
            parameters.event{indPP}{indBlock}.Number(StimOn(end):end) = [];
            parameters.event{indPP}{indBlock}.Times(StimOn(end):end) = [];
            
            StimOn(end) = [];
            Responses(Responses > StimOn(end)) = [];
        end
        
        parameters.experiment{indPP}{indBlock}.RespT   = parameters.event{indPP}{indBlock}.Times(Responses)*(1/parameters.eeg.SampleRate);
        parameters.experiment{indPP}{indBlock}.TargOnT = parameters.event{indPP}{indBlock}.Times(StimOn)*(1/parameters.eeg.SampleRate);
        parameters.experiment{indPP}{indBlock}.RespLR = ones(size(parameters.experiment{indPP}{indBlock}.RespT));
    end
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------   Epoching data    -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters.epochingData;
% parameters.applyCSD; % NOTE already applied in the epoched EEG data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------   Behavioural data    -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we create equal sized Reaction Time (RT) bins per participant per conditions
% to check the alignment of the CPP with the actual median RT as the CPP
% peak should be highly linked with the reaction time.
parameters.binRTs(1, [1/3 2/3], [1 2]); % as there are sometimes to little trials per all conditions, we need to resort to a smaller number. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------     Waveform plotting       -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now transfer parameters to EEG object. These can be run without each
% other.
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------      get CPP        -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveCPPhere = fullfile(parameters.figFolder, 'groupAverage\HPF0_CSD1\CPP\');

% define topoplot information
trangeTopo  =  0 + [-2 -1]*1/parameters.stim.freqSSVEP; % note if you add a second row, these will be used for FA plots. Only nesc. when you are using target-locked topos
TargetOrResponse = 2;       % as written 1) target-locked 2) response-locked

% define rest
baselineCorrect  = 1;
methodUsed       = 2;       % Method used can be 1) choose spec. electrodes, 2) choose cluster of electrodes and get the 3 best based on SNR, 3) apply lucalize
plotThis         = 3;       % 0) only topoplots
                            % 1) individual plots, 
                            % 2) average Hits and seperated average Misses and FA
                            % 3) ,, include all average target-locked Misses
                            % 4) plot target and response seperately.
                            % 5) ,, include all average target-locked
                            % Misses  (NOTE this doesn't change the
                            % response-locked ones).
                            % 6) get first derivative
                            
grouping         = 0;       % if you want to include seperated plots for a specific condition in plotCondition (see below) for example 1 == plotting first condition in plotCondition
plotCondition    = [1 2];   % order in which it is plotted, look through trialmatrix and getConditions to better understand everything. 

cppElec          = 'A5 A19 A32 A18 A20 A31 A17 A21 A30';

% 1)    Plot target and response locked CPP (target-locked including all misses)
%       as a function of target context as well as target evidence strength. 
%       (SEE FIGURE 2D AND TOPOGRAPHY)
[Quality.CPP, fig]  = parameters.plotERP(methodUsed, 'CPP', cppElec, baselineCorrect, ...
   plotThis, plotCondition, grouping, trangeTopo, TargetOrResponse);

% some small addition to finish the figures for publication
figure(fig.Handle{1});
fig.Info{1}(1,2).select()
legend('off')
ylim([-3 35])
yticks(0:10:30)

getAxes = gca;
shaded = shadedErrorBar(trangeTopo, repmat(nanmean(getAxes.YLim),2,1), repmat(getAxes.YLim(2) - nanmean(getAxes.YLim),2,1));
shaded.mainLine.LineStyle = 'none';
shaded.edge(1).LineStyle  = 'none'; shaded.edge(2).LineStyle  = 'none';
shaded.patch.FaceAlpha    = 0.1;
fig.Info{1}(1,1).select()
yticks(0:10:30)

plotSave(gca, 'finalFigCPP.png', saveCPPhere, parameters.figLayOut.saveDim)

% 2)    Plot to empiral target-locked and response-locked to
%       compare with the sim. CPP (SEE FIGURE 3D and FIGURE 3 suppl. 1A )
parameters.figLayOut.saveDim     = [4 6];
parameters.figLayOut.targetLim   = [0 0.5:0.5:1];
parameters.figLayOut.plotCI      = 0;

plotThis = 5; % Plot first derivative

[~, fig]  = parameters.plotERP(methodUsed, 'CPP',  cppElec, baselineCorrect,...
    plotThis, plotCondition, grouping, trangeTopo, TargetOrResponse);

figure(fig.Handle{1})
ylim([-1 30])
yticks(0:10:30)
legend('off')
ylabel('')
plotSave(gca, 'finalFigCPPModel.png', saveCPPhere, parameters.figLayOut.saveDim);

% 3)    Plot to first derivative empiral target-locked and 
%       response-locked to compare with the sim. CPP (SEE FIGURE 3D and FIGURE 3 suppl. 1A )
plotThis  = 6; 
[~, fig]  = parameters.plotERP(methodUsed, 'CPP',  cppElec, baselineCorrect,...
    plotThis, plotCondition, grouping, trangeTopo, TargetOrResponse);

figure(fig.Handle{1})
ylim([-0.2 0.3])
yticks(-0.2:0.2:.2)
legend('off')
plotSave(gca, 'finalFigderCPPModel.png', saveCPPhere, [parameters.figLayOut.saveDim(1) 6.27]);

% 4)    Plot difference in topoplot to check selectivity of effects.
%       (NOT IN PAPER)
%{
parameters.plotDifferenceTopo('CPP', cppElec, [2 1], trangeTopo,TargetOrResponse); % TODO check the slope stuff here
%}

% 5)    Plot staticis of pre-response CPP.
%       (SEE FIGURE 2E)
statsRange(1,:)  = trangeTopo;
statsRange(2,:)  = [-.25 -1*1/parameters.stim.freqSSVEP];
PeakMeanOrMax    = 1;
negOrPos         = 2;
RTbins           = [];

parameters.figLayOut.saveDim   =  [4.173 2.5];
tblCPP = parameters.getWaveformParameters(methodUsed, 'CPP', RTbins, baselineCorrect,...
    1, plotCondition, grouping, statsRange, TargetOrResponse, negOrPos, PeakMeanOrMax, [0 30]);

% preform linear mixed model on the pre-response CPP amplitude.
tblCPP.Context = categorical(tblCPP.Context);
tblCPP.Evidencestrength = categorical(tblCPP.Evidencestrength);
lmeCPP = fitglme(tblCPP,  'Peak ~  1 +  Context*Evidencestrength*RT + (1|ppNames)',...
    'FitMethod','Laplace'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------      get N2      ----------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for supplementary figures
saveN2here = fullfile(parameters.figFolder, 'groupAverage\HPF0_CSD1\N2\');

trangeTopo       = [0.24 0.34]; % target-locked now --> based on interaction effect found.
TargetOrResponse = 1;           % as written 1) target-locked 2) response-locked

baselineCorrect  = 1;
methodUsed       = 1;   % Method used can be 1) choose spec. electrodes, 2) choose cluster of electrodes and get the 3 best based on SNR, 3) apply lucalize
plotThis         = 5;   % TODO check
N2Elec           = 'D30 D31 B12 B11';

% 1)    Plot target and response locked N2 (target-locked including all misses)
%       as a function of target context as well as target evidence strength.
%       (SEE FIGURE 2 - supplement 1A and topography);
[QualityN2, fig]  = parameters.plotERP(methodUsed, 'N2', N2Elec, baselineCorrect,...
   plotThis, plotCondition, grouping, trangeTopo, TargetOrResponse);

figure(fig.Handle{1});
legend('off')
ylim([-25 5])
getAxes = gca;
ylabel('')

% check for effects of context and motion coherence
shaded = shadedErrorBar(trangeTopo,...
    repmat(nanmean(getAxes.YLim),2,1), repmat(getAxes.YLim(2) - nanmean(getAxes.YLim),2,1));
shaded.mainLine.LineStyle = 'none';
shaded.edge(1).LineStyle = 'none'; shaded.edge(2).LineStyle = 'none';
shaded.patch.FaceAlpha = 0.05;
legend('off')
      
plotSave(gca, 'N2_shaded.png', saveN2here, [4.1 6.1])

% 2)    Plot target and response locked N2 (target-locked including all misses)
%       as a function of target context as well as target evidence strength.
%       (SEE FIGURE 2 - supplement 1C);
parameters.plotDifferenceTopo('N2',  N2Elec, [2 1], trangeTopo, TargetOrResponse);

% 3)    Plot target and response locked N2 (target-locked including all misses)
%       as a function of target context as well as target evidence strength.
%       (SEE FIGURE 2 - supplement 1B);
TargetOrResponse = 1;
statsRange(1,:)  = trangeTopo;
statsRange(2,:)  = [0 0.3];
PeakMeanOrMax    = 1;
negOrPos         = 1;

tblN2 = parameters.getWaveformParameters(methodUsed, 'N2', [1/3 2/3], baselineCorrect,...
    1, plotCondition, grouping, statsRange, TargetOrResponse, negOrPos, PeakMeanOrMax, [-25 5]);

xlim([0.3 1.05])
xticks(0.4:0.2:1.1);
xlabel(' ')

plotSave(gca, 'finalN2peak.png', saveN2here, [4.1 3.8]);

% get stats
tblN2.Context = categorical(tblN2.Context);
tblN2.Evidencestrength = categorical(tblN2.Evidencestrength);
lme.N2 = fitglme(tblN2,  'Peak ~  1 + Context*Evidencestrength*RT + (1|ppNames)',...
    'FitMethod','Laplace');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------   EVENT RELATED BETA DESYNCHRONIZATION --------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveBetahere = fullfile(parameters.figFolder, 'groupAverage\HPF0_CSD1\Beta ERD\');

trangeTopo  = 0 + [-1 0]*1/parameters.stim.freqSSVEP;
TargetOrResponse = 2;

baselineCorrect  = 2;    % baseline-correct towards response period!
parameters.eeg.baseline =  0 + [-2 -1]*1/parameters.stim.freqSSVEP; % Motor Execution threshold

parameters.figLayOut.targetLim   = [-0.2 0.5:0.5:1];
parameters.figLayOut.plotCI      = 0.05;
parameters.figLayOut.saveDim     = [5 11];

methodUsed       = 1;
plotThis         = 2; % 3) plot everything together

% use specific code as the people werent instructed to use the right hand. 
% handmess(parameters, 'Beta ERD', [], trangeTopo, TargetOrResponse, 1, [13:30]);
possElec{1} = 'D20 D19 D18 D27 D28 D17'; % right handed response
possElec{2} = 'B21 B22 B23 B19 B18 B17'; % left handed response

betaElec{1} = possElec{1};
betaElec{2} = possElec{1};
betaElec{3} = possElec{1};
betaElec{4} = possElec{1};
betaElec{5} = possElec{2};
betaElec{6} = possElec{1};
betaElec{7} = possElec{1};
betaElec{8} = possElec{1};
betaElec{9} = possElec{1};
betaElec{10} = possElec{1};
betaElec{11} = possElec{1};
betaElec{12} = possElec{2};
betaElec{13} = possElec{1};
betaElec{14} = possElec{1};

% 1)    Plot target and response locked Beta (target-locked including all misses)
%       as a function of target context as well as target evidence strength. 
%       First, baseline-correct to pre-response Beta threshold (SEE FIGURE 2A).
% first for topoplot, with usual baseline, just to see the pre-response
% BETA topographic distribution... This will cost a bit more time this
% way
parameters.plotERP(methodUsed, 'Beta ERD', betaElec,...
    1, 0, plotCondition, grouping, trangeTopo,...
    TargetOrResponse, 2, 15:30); % NOTE, 15 hz is excluded for the SSVEP

[QualityBeta, fig] = parameters.plotERP(methodUsed, 'Beta ERD', betaElec,...
    baselineCorrect, plotThis, plotCondition, grouping, trangeTopo,...
    TargetOrResponse, 2, 15:30); % NOTE, 15 hz is excluded for the SSVEP

figure(fig.Handle{1});
fig.Info{1}(1,2).select()

legend('off')
ylim([-0.1 0.4])
yticks(0:0.2:0.4)
line([0 0], [-0.1 0.8], 'Color', 'k', 'LineWidth', 2)

fig.Info{1}(1,1).select()
yticks(0:0.2:0.8)
line([0 0], [-0.1 0.8], 'Color', 'k', 'LineWidth', 2)

% add shaded grey bar to show topoplot range in figure.
getAxes = gca;
shaded = shadedErrorBar(trangeTopo, repmat(nanmean(getAxes.YLim),2,1), repmat(getAxes.YLim(2) - nanmean(getAxes.YLim),2,1));
shaded.mainLine.LineStyle = 'none';
shaded.edge(1).LineStyle  = 'none'; shaded.edge(2).LineStyle  = 'none';
shaded.patch.FaceAlpha    = 0.1;

plotSave(gca, 'finalFigBeta.png', saveBetahere, parameters.figLayOut.saveDim)

% 2)    Plot difference, showing specificity of context effects 
%       (SEE topographic plots in FIGURE 2A).
parameters.plotDifferenceTopo('Beta ERD', [], plotCondition, trangeTopo,...
    TargetOrResponse, 1, 15:30, 0, 0);

% 3)    Get all stats of excursion
%       (SEE FIG 2B and FIG 2A topography)
baselineCorrect  = 2;
TargetOrResponse = 1;  % take baseline, to get differences rather than topoplot. 
negOrPos         = 1;
PeakMeanOrMax    = 1;
trangeTopo(2,:)  =  [-0.5 0];  % for slope, not used here though. 
parameters.figLayOut.saveDim   =  [4.173 2.5];

tbl.Excursion =  parameters.getWaveformParameters(methodUsed, 'Beta ERD', [], baselineCorrect,...
    1, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos,PeakMeanOrMax);

% preform linear mixed model on the pre-response CPP amplitude.
tbl.Excursion.Context = categorical(tbl.Excursion.Context);
tbl.Excursion.Evidencestrength = categorical(tbl.Excursion.Evidencestrength);

% check for RT interaction effects:
lme.RTinteraction = fitglme(tbl.Excursion,  'Peak ~  1 + Context*Evidencestrength*RT + (1|ppNames)',...
    'FitMethod','Laplace');

lme.Excursion = fitglme(tbl.Excursion,  'Peak ~  1 + Context*Evidencestrength + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.FixedExcursion = tbl.Excursion(double(tbl.Excursion.Context) == 1,:);
lme.FixedExcursion  = fitglme(tbl.FixedExcursion ,  'Peak ~  1 + Evidencestrength + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.MixedExcursion = tbl.Excursion(double(tbl.Excursion.Context) == 2,:);
lme.MixedExcursion = fitglme(tbl.MixedExcursion,  'Peak ~  1 + Evidencestrength + (1|ppNames)',...
    'FitMethod','Laplace');   

tbl.WeakExcursion = tbl.Excursion(double(tbl.Excursion.Evidencestrength) == 1,:);
lme.WeakExcursion  = fitglme(tbl.WeakExcursion ,  'Peak ~  1 + Context + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.StrongExcursion = tbl.Excursion(double(tbl.Excursion.Evidencestrength) == 2,:);
lme.StrongExcursion = fitglme(tbl.StrongExcursion,  'Peak ~  1 + Context + (1|ppNames)',...
    'FitMethod','Laplace'); 

diary(fullfile(parameters.logFolder, 'WaveformsStats', 'lmeERD.txt'))
clc
fprintf('Waveform stats checking for excursion effects.\nInitial assumption that they are context and RT-invariant by baseline correction to motor execution threshold.\n');
fprintf('--------------------------------------------------------------------------------------------------------------\n');
lme.RTinteraction
fprintf('While there is an interaction effect (e.g. smaller excursion relates to faster RTs,\nthere are no interaction effects\n');
fprintf('Remove RT effects\n');
fprintf('--------------------------------------------------------------------------------------------------------------\n');
lme.Excursion
fprintf('Interaction effect of Context and Motion Coherence\n');
fprintf('--------------------------------------------------------------------------------------------------------------\n');
fprintf('Posthocs of Context effects show that there only is a difference between Weak and Strong in fixed condition=\n');
Fixed = lme.FixedExcursion
fprintf('And not in the mixed condition\n');
Mixed = lme.MixedExcursion
fprintf('--------------------------------------------------------------------------------------------------------------\n');
fprintf('--------------------------------------------------------------------------------------------------------------\n');
fprintf('Posthocs of Motion Coherence effects show that there only is a difference between Strong in fixed and mixed condition=\n');
Fixed = lme.StrongExcursion
fprintf('but not for Weak\n');
Mixed = lme.WeakExcursion
fprintf('--------------------------------------------------------------------------------------------------------------\n');
diary off

% 3)    Plot target and response locked Beta (target-locked including all misses)
%       as a function of target context as well as target evidence strength. 
%       No baseline-correction (SEE FIGURE 4A).
baselineCorrect = 0;
parameters.figLayOut.saveDim     = [5 11];
[Quality, fig] = parameters.plotERP(methodUsed, 'Beta ERD', betaElec,...
    baselineCorrect, plotThis, plotCondition, grouping, trangeTopo,...
    TargetOrResponse, 2, 15:30);

figure(fig.Handle{1});
fig.Info{1}(1,2).select()

legend('off')
ylim([3.5 5])
yticks([3.5:0.5:5])
line([0 0], [3.5 5], 'Color', 'k', 'LineWidth', 2)

figure(fig.Handle{2});
fig.Info{2}(1,2).select()

legend('off')
ylim([3.5 5])
yticks([3.5:0.5:5])
line([0 0], [3.5 5], 'Color', 'k', 'LineWidth', 2)

% add stats time area
getAxes = gca;
shaded = shadedErrorBar(0 + [-2 -1]*1000/parameters.stim.freqSSVEP, repmat(nanmean(getAxes.YLim),2,1), repmat(getAxes.YLim(2) - nanmean(getAxes.YLim),2,1));
shaded.mainLine.LineStyle = 'none';
shaded.edge(1).LineStyle  = 'none'; shaded.edge(2).LineStyle  = 'none';
shaded.patch.FaceAlpha    = 0.1;

fig.Info{1}(1,1).select()
ylim([3.5 5])
yticks([3.5:0.5:5])
line([0 0], [3.5 5], 'Color', 'k', 'LineWidth', 2); 
getAxes = gca;
shaded = shadedErrorBar( 0 + [-2 0]*1000/parameters.stim.freqSSVEP, repmat(nanmean(getAxes.YLim),2,1), repmat(getAxes.YLim(2) - nanmean(getAxes.YLim),2,1));
shaded.mainLine.LineStyle = 'none';
shaded.edge(1).LineStyle  = 'none'; shaded.edge(2).LineStyle  = 'none';
shaded.patch.FaceAlpha    = 0.1;

plotSave(gca, 'finalFigBeta_noBaseline.png', saveBetahere, parameters.figLayOut.saveDim)



% ITI period
%{
keyboard

baselineCorrect = 3; %ColoursCond12;%
parameters.figLayOut.colours  = ColoursCond12([1 2 4],:); %ColoursCond12([1 1 1 1 2 2 2 2 4 4 4 4],:); %
parameters.figLayOut.lineType = repmat({'--'},4,1) ;%
parameters.figLayOut.saveDim 
for indPP = 1:length(parameters.ppNames)
    parameters.behaviour{indPP}.trialMatrix(:,5) = nan;
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,3) == 8,5) = 1;
end
parameters.conditions{5} = 1;
parameters.figLayOut.legNames{5} = '8 sec';
parameters.figLayOut.legTitle{5} = 'ITI';

% ITI plotting.
plotCondition = [4 5]; grouping =0 ;
parameters.figLayOut.plotCI = 0;
fig = parameters.plotITI(methodUsed,  'Beta ERD', Elec, baselineCorrect,...
    plotCondition, grouping, -67 + [-1 0]*1000/parameters.stim.freqSSVEP);

getAxis = gca;
%getAxis.YAxis.Visible = 'off'; 
xlim([-8000 500])
ylim([0 0.5])
yticks([0:0.25:0.5])

xlabel (' ')
xticks([-8000:2000:0])
xticklabels([0:2:8])
%line([-8000 -8000], [180 0], 'Color', 'black', 'LineWidth', 2)
line( [-8000 500],[0 0], 'Color', 'black', 'LineWidth', 2)
legend('off')

% legend(fig.legendThis, {'Weak', 'Strong',  'Mixed'}, 'Location', 'northwest')
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);

plotSave(gca,  ['Beta ERD_ITI.png'], hereFigFolder, [4 6]);% 

fig.legendThis.delete

ylim([50 135])
getAxis.YAxis.Visible = 'off'; 
getAxis.XAxis.Visible = 'off'; 
xticks([-8000:2000:0])
xlim([-8000 1000])
xlabel('off')

plotSave(gca,  ['Beta ERD_ITI_model.png'], hereFigFolder, [1.4 2.1].*1.5);% 

keyboard
%}


keyboard
% 2) baseline effects
%
PeakMeanOrMax    = 1;
plotCondition    = [1 2];
negOrPos = 1;
baselineCorrect  = 0;
TargetOrResponse = 1;
trangeTopo(1,:)  =  0 + [-4 -1]*1000/parameters.stim.freqSSVEP; %

% parameters.getWaveformParameters(methodUsed, 'Beta ERD', [1/3 2/3], baselineCorrect,...
%     1, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos,PeakMeanOrMax);
% 
% ylim([3.4 5])
% xlim([300 1050])
% xticks(400:200:1100);
% xticklabels((400:200:1100)./1000)
% set(gca,'FontSize', parameters.figLayOut.letterSize);
% set(gca,'FontName', parameters.figLayOut.letterType);
% % set(gca, 'YDir','reverse')
% 
% plotSave(gca, ['amplitudeBaseline.png'], hereFigFolder, [3.7 4.5]);%[4.3 obj.figLayOut.saveDim(2)/3]);

parameters.plotDifferenceTopo('Beta ERD_baseline', [], [2 1], trangeTopo(1,:),...
    TargetOrResponse, 0, 1, 15:30, 0, 0);


parameters.plotDifferenceTopo('Beta ERD_baseline', [], [5 6], trangeTopo(1,:),...
    TargetOrResponse, 0, 1, 15:30, 0,0);

tbl.BetaBaseline = parameters.getWaveformParameters(methodUsed, 'Beta ERD', [], baselineCorrect,...
    0, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos,PeakMeanOrMax);


lme.RTinteraction = fitglme(tbl.BetaBaseline,  'Peak ~  1 + Context*Evidencestrength*RT + (1|ppNames)',...
    'FitMethod','Laplace');

% check context
tbl.Weak = tbl.BetaBaseline(double(tbl.BetaBaseline.Evidencestrength) == 1,:);
lme.Weak = fitglme(tbl.Weak,  'Peak ~  1 + Context*RT + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.Strong = tbl.BetaBaseline(double(tbl.BetaBaseline.Evidencestrength) == 2,:);
lme.Strong = fitglme(tbl.Strong,  'Peak ~  1 + Context*RT + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.Strong = tbl.BetaBaseline(double(tbl.BetaBaseline.Evidencestrength) == 2 & double(tbl.BetaBaseline.Context) == 1,:);
lme.Strong1 = fitglme(tbl.Strong,  'Peak ~  1 + RT + (1|ppNames)',...
    'FitMethod','Laplace');
tbl.Strong= tbl.BetaBaseline(double(tbl.BetaBaseline.Evidencestrength) == 2 & double(tbl.BetaBaseline.Context) == 2,:);
lme.Strong2 = fitglme(tbl.Strong,  'Peak ~  1 + RT + (1|ppNames)',...
    'FitMethod','Laplace');

diary(fullfile(parameters.logFolder, 'WaveformsStats', 'lmeERD_baseline.txt'))
clc
fprintf('Waveform stats checking for baseline effects.\nBaseline activity shows a signifant interaction effect bewteen RT*context*evidencestrength\n');
fprintf('--------------------------------------------------------------------------------------------------------------\n');
lme.RTinteraction
fprintf('While lawfully lower beta synchronization causes faster reaction times in all conditions but the Strong Context\n');
fprintf('--------------------------------------------------------------------------------------------------------------\n');
fprintf('--------------------------------- WEAK -------------------------------------------\n');
lme.Weak
fprintf('--------------------------------- STRONG -------------------------------------------\n');
lme.Strong
fprintf('--------------------------------- interaction effect -------------------------------------------\n');
lme.Strong1
lme.Strong2
fprintf('--------------------------------------------------------------------------------------------------------------\n');
diary off
%}


% 2) Pre response effects
%
baselineCorrect  = 0;
TargetOrResponse = 2;
PeakMeanOrMax    = 2;

trangeTopo(1,:)  =  0 + [-2 0]*1000/parameters.stim.freqSSVEP; %
trangeTopo(2,:)  =  0 + [-500 0]*1000/parameters.stim.freqSSVEP; %

% parameters.getWaveformParameters(methodUsed, 'Beta ERD', [1/3 2/3], baselineCorrect,...
%     1, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos,PeakMeanOrMax);
% 
% ylim([3.4 5])
% xlim([300 1050])
% xticks(400:200:1100);
% xticklabels((400:200:1100)./1000)
% 
% set(gca,'FontSize', parameters.figLayOut.letterSize);
% set(gca,'FontName', parameters.figLayOut.letterType);
% % set(gca, 'YDir','reverse')
% 
% plotSave(gca, ['amplitudePreResponse.png'], hereFigFolder, [3.7 4.5]);%[4.3 obj.figLayOut.saveDim(2)/3]);

tbl.BetaPreresponse = parameters.getWaveformParameters(methodUsed, 'Beta ERD', [], baselineCorrect,...
    0, plotCondition, grouping, trangeTopo, TargetOrResponse, negOrPos,PeakMeanOrMax);

parameters.plotDifferenceTopo('Beta ERD_Preresponse', [], [2 1],...
    0 + [-2 2]*1000/parameters.stim.freqSSVEP, TargetOrResponse,...
    0, 1, 15:30,PeakMeanOrMax,0);

parameters.plotDifferenceTopo('Beta ERD_Preresponse', [], [6 5],...
    0 + [-2 2]*1000/parameters.stim.freqSSVEP, TargetOrResponse,...
    0, 1, 15:30,PeakMeanOrMax,0);
% 

lme.RTinteraction = fitglme(tbl.BetaPreresponse,  'Peak ~  1 + Context*Evidencestrength*RT + (1|ppNames)',...
    'FitMethod','Laplace');

% check context
tbl.Weak = tbl.BetaPreresponse(double(tbl.BetaPreresponse.Evidencestrength) == 1,:);
lme.Weak = fitglme(tbl.Weak,  'Peak ~  1 + Context*RT + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.Strong = tbl.BetaPreresponse(double(tbl.BetaPreresponse.Evidencestrength) == 2,:);
lme.Strong = fitglme(tbl.Strong,  'Peak ~  1 + Context*RT + (1|ppNames)',...
    'FitMethod','Laplace');

tbl.Strong  = tbl.BetaPreresponse(double(tbl.BetaPreresponse.Evidencestrength) == 2 & double(tbl.BetaPreresponse.Context) == 1,:);
lme.Strong1 = fitglme(tbl.Strong,  'Peak ~  1 + RT + (1|ppNames)',...
    'FitMethod','Laplace');
tbl.Strong  = tbl.BetaPreresponse(double(tbl.BetaPreresponse.Evidencestrength) == 2 & double(tbl.BetaPreresponse.Context) == 2,:);
lme.Strong2 = fitglme(tbl.Strong,  'Peak ~  1 + RT + (1|ppNames)',...
    'FitMethod','Laplace');
% 
% diary(fullfile(parameters.logFolder, 'WaveformsStats', 'lmeERD_preResponse.txt'))
% clc
% fprintf('Waveform stats checking for baseline effects.\nBaseline activity shows a signifant interaction effect bewteen RT*context*evidencestrength\n');
% fprintf('--------------------------------------------------------------------------------------------------------------\n');
% lme.RTinteraction
% fprintf('While lawfully lower beta synchronization causes faster reaction times in all conditions but the Strong Context\n');
% fprintf('--------------------------------------------------------------------------------------------------------------\n');
% fprintf('--------------------------------- WEAK -------------------------------------------\n');
% lme.Weak
% fprintf('--------------------------------- STRONG -------------------------------------------\n');
% lme.Strong
% fprintf('--------------------------------- interaction effect -------------------------------------------\n');
% lme.Strong1
% lme.Strong2
% fprintf('--------------------------------------------------------------------------------------------------------------\n');
% diary off
% %}
%{
all = tbl.BetaPreresponse;
all.Baseline = tbl.BetaBaseline.Peak;

tbl.Weak  = all(double(tbl.BetaPreresponse.Evidencestrength) == 1 & double(tbl.BetaPreresponse.Context) == 1,:);
lme.Weak1 = fitglme(tbl.Weak,  'Peak ~  1 + RT*Baseline + (Baseline-1|ppNames)',...
    'FitMethod','Laplace');

tbl.Weak  = all(double(tbl.BetaPreresponse.Evidencestrength) == 1 & double(tbl.BetaPreresponse.Context) == 2,:);
lme.Weak2 = fitglme(tbl.Weak,  'Peak ~  1 + RT*Baseline + (Baseline-1|ppNames)',...
    'FitMethod','Laplace');

tbl.Strong  = all(double(tbl.BetaPreresponse.Evidencestrength) == 2 & double(tbl.BetaPreresponse.Context) == 1,:);
lme.Strong1 = fitglme(tbl.Strong,  'Peak ~  1 + RT*Baseline + (Baseline-1|ppNames)',...
    'FitMethod','Laplace');
tbl.Strong  = all(double(tbl.BetaPreresponse.Evidencestrength) == 2 & double(tbl.BetaPreresponse.Context) == 2,:);
lme.Strong2 = fitglme(tbl.Strong,  'Peak ~  1 + RT*Baseline + (Baseline-1|ppNames)',...
    'FitMethod','Laplace');

diary(fullfile(parameters.logFolder, 'WaveformsStats', 'lmeERD_preResponseconvBaseline.txt'))
clc
fprintf('--------------------------------- WEAK fixed -------------------------------------------\n');
lme.Weak1 

fprintf('--------------------------------- WEAK fixed -------------------------------------------\n');
lme.Weak2

fprintf('--------------------------------- Strong fixed -------------------------------------------\n');
lme.Strong1 

fprintf('--------------------------------- Strong fixed -------------------------------------------\n');
lme.Strong2
fprintf('--------------------------------------------------------------------------------------------------------------\n');
diary off
%}

