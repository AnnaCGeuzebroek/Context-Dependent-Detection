%% Behaviorual and neurally-informed modelling
% Data collected Hannah Craddock for her Master's thesis, 2016:
%   'Electrophysiological and behavioural indices of decision criterion
%   adjustments across context of weak and strong evidence.'
%
% Current code is written by A.C. Geuzebroek and part of:
%   Geuzebroek AC, Craddock H, O’Connell RG, & Kelly SP (2022).
%   Balancing true and false detection of intermittent sensory targets
%   by adjusting the inputs to the evidence accumulation process (https://biorxiv.org/cgi/content/short/2022.09.01.505650v1)
%
% In short, participants were asked to continous monitor a cloud of randomly moving
% dots for intermittered target defined as upwards coherent moving dots.
% Afterwhich they were asked to response as fast and accurate as possible.
% Difficulty context was manipulated with:
%   1) Hard  (25% motion coherence)
%   2) Easy  (70% motion coherence)
%   3) Mixed (25% and 70% motion coherence with equal probabiltiy).
% 
% Code requires the data to be structured as follows:
%   InputFolder     = BASEFOLDER\Data\Raw\
%       Participant folder = BASEFOLDER\Data\Raw\PPNAME\ (note never use initials!)
%           EEG folder           = BASEFOLDER\Data\Raw\PPNAME\EEG data\
%           Eyelink folder       = BASEFOLDER\Data\Raw\PPNAME\Eyelink data\
%           Trial files folder   = BASEFOLDER\Data\Raw\PPNAME\Trial files\

% costum-made code:
%   1) dataAnalysis (object to access all the functions, OR get data
%   structure set-up as see below) --> here the data is cut to get the
%   conditions/epochs and will allow and get the EEG data needed for
%   constraining the NI models. 
%   2) dataModelling 
%
%
% Depending on:
%   1) EEGLAB (including Biosig extention).
%   2) CSD toolbox and lay-out.
%   3) findNoisyChannels
%   4) Brewermap (to get colors for plot. Can be easily replaced with just choosing colours)
%   5) panels
%   6) fminsearchbnd

clear
clc
clear global

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------     Input parameters    -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set all paths.
if strcmp(computer, 'PCWIN64')
    error('Please set your data folder')
    inputFolder = 'YOURFOLDER\Data';    % add base folder of data. 
    addpath(genpath('C:\Users\eng\Documents\MATLAB\dataAnalysis'))
    
    error('Please add the path to EEGLAB!')

    if ~exist('eeglab', 'file')
        run('YOURPATHTOEEGLAB\EEGLAB\eeglab'); close all;
    end
end

% some behavioral task that can be check. looking at time on task
% (Within-block effects) or time doing the experiments
blockRandomized = 1; % yes, as you have young and old groups.
numRTBins  = 2; 
TimeOnTask = 0;        
TimeOnExperiment = 0;

% Set conditions, here this is left or right tilted as well as inter-trial
% duration. For initial purposes we do not have to look at these, except
% that we bin the reaction times for each conditions to account of possible
% effects.

% IMPORTANT condNames should refer to the way that it was saved in the
% experiment. This can be checked by loading the trial information file and
% see how it will be set in the workspace. Same is happening with
% conditions, e.g. for the left vs. right tilded contrast change, it is
% saved as trialLR as condNames and 1 for left and 2 for right in the so-called
% trialMatrix (matrix that will track the experimental conditions).
conditions{1}  = [1 2]; % easy difficult mixed
condNames{1}   = 'condition';
conditions{2}  = [25 70];
condNames{2}   = 'coherence';

% get your color scheme set here
tmp =  flipud(brewermap(4,'RdBu')); ColoursCond12  = [128 205 193; 1 133 113; 223 194 125; 166 97 26]/255;% [brewermap(2, 'Greens'); tmp([2 1],:)];
ColoursCond3 =  brewermap(4,'Greys'); ColoursCond3 = ColoursCond3-0.16;

% preset the signal processing object.
parameters  = dataAnalysis.getInstance();

% set a couple of standard parameters within this object.
parameters.system = 'bioSemi';   % Object will choice the approtiated loading function of EEGLAB
parameters.analysisEEG      = 1; % Analyse EEG (obviously we want to look at this)
parameters.analysisEyelink  = 0;
parameters.analysisEOG      = 1; % Use additionally EOG electrodes and Frontal electrodes
% to look for blinks and
% exclude trials with blinks (We do this to prevent
% any possible effects of the blink on the sensory perception)
parameters.analysisEMG = 0;
parameters.numTrials   = 24; % number of trials per block.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------     GET EXPERIMENT INFORMATION    ------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function looping through the data files of each participants. Identifying
% the behavioural information. including numBlocks, numBlocks per
% condition, numStimulus in each conditions and their specific condition (
% here initialy all stimulus in each block are 1 condition.) However, it is
% possible that the conditions information are in the EEG files. So here
% simply put all the par and other parameters of each participant away.
parameters.getInformation(inputFolder, conditions, condNames, blockRandomized);
clear cond* *Folder


% recode the experiment.trialCond from 1:4 as 1:2 is no gap, 4:6 is gap. As
% well as 1,2,3 being [2 4 6] ITI respectively.
parameters.conditions{3}   = [2 4 6 8];
parameters.condNames{3}    = 'ITI';

parameters.stim.lengthITI  = [2 4 6 8];
parameters.stim.namesITI   = 'ITI';

for indPP = 1:length(parameters.ppNames)
    for indBlock = 1:parameters.numBlocks(indPP)
        
        tmpITI = zeros(size(parameters.experiment{indPP}{indBlock}.trialCond));
        
        if length(parameters.experiment{indPP}{indBlock}.par.correctCoh) == 1
            parameters.experiment{indPP}{indBlock}.coherence  = parameters.experiment{indPP}{indBlock}.par.correctCoh{:};
            parameters.experiment{indPP}{indBlock}.condition  = 1;%find(parameters.conditions{2} == parameters.experiment{indPP}{indBlock}.par.correctCoh{:});
            
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

% assuming that parameters are the same for all participants and with all
% blocks
parameters.DetectOrDisc = 0;

parameters.stim.refreshRate  = parameters.experiment{1}{1}.par.videoFrate;	 % Screens refreshRate
parameters.stim.duration     = parameters.experiment{1}{1}.par.targetDur.*1000;    % duration of the stimulus
parameters.stim.freqSSVEP    = 15;
parameters.stim.epoch        = [-1800 2500];

parameters.stim.RTdeadLine = [-1.8 1.65];
parameters.stim.RTCutOff   = 0.25; % only include reaction times after the gap offset.
parameters.stim.FACutOff   = []; % only include reaction times after the gap offset.

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
% load in a structure 'chanlocs' containing the x,y,z locations of each of the 128 scalp channels in the cap.
load chanlocsBioSemi128;
elab = strings(parameters.eeg.NumberOfChannels,1);
for indChan = 1:parameters.eeg.NumberOfChannels
    elab(indChan) = chanlocs(indChan).labels;
end
elab(end+1)='FCz'; clear ind*

parameters.eeg.ChannelsName   	= elab;
parameters.eeg.SampleRate   	= 512;	% EEG sample rate

parameters.eeg.chanlocs         = chanlocs;	% electrode locations
parameters.eeg.transChanlocs	= load('CSD_coords.mat');	% electrode locations

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
parameters.eeg.simonHPF     = 0;
parameters.eeg.HPFcutoff    = 0;
parameters.eeg.applyDetrend	= 1;    % option 1(all)/ 2(only that will benefit from it)

parameters.eeg.artifactThres = 80;
parameters.eeg.artifactEOG   = 200;
parameters.eeg.intArtifact   = 1;

parameters.eeg.channelVEOG   = 1:2;
parameters.eeg.channelHEOG   = 3:4;

parameters.eeg.timing = 'past';

%--------------- set triggers and EPOCHS definition -----------------------
% The trigger codes we used for the different events were as follows:
parameters.triggers.start       = 1;            % Start of the trial
parameters.triggers.stimulusOFF = 5;            % Target off
parameters.triggers.stimulusON  = [170 125];    % 170 for 70% and 125 for 25%
parameters.triggers.response    = [12 13];      % left vs. right

% EPOCH DEFINITION
parameters.eeg.epochLock = parameters.triggers.stimulusON;
parameters.eeg.baseline  = 0 + [-1 1]*1000/parameters.stim.freqSSVEP;

% Initially we set a rather larger epoch based on the stimulus definitition
% 'parameters.stim.epoch'. This is used to cut out a larger epoch to apply
% all the ERP pre-processing on. Afterwards smaller epochs are created to
% get the target- and response-locked epochs. This is usefull especially
% for the ERD (event related desynchronization) and the SSVEP later
% onwards.

parameters.eeg.targetEpoch   = -ceil(500/(1000/parameters.eeg.SampleRate))*(1000/parameters.eeg.SampleRate):1000/parameters.eeg.SampleRate:ceil(parameters.stim.RTdeadLine(2)*1000/(1000/parameters.eeg.SampleRate))*(1000/parameters.eeg.SampleRate);
parameters.eeg.responseEpoch = -ceil(800/(1000/parameters.eeg.SampleRate))*(1000/parameters.eeg.SampleRate):1000/parameters.eeg.SampleRate:ceil(400/(1000/parameters.eeg.SampleRate))*(1000/parameters.eeg.SampleRate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------    Perform eeg prep  -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters.applyPreprocessing;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------- Reset behaviour data    -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------   Epoching data    -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
parameters.epochingData;
% applyPreProcess = toc
% tic
% parameters.applyCSD;
% applyCSD = toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------   Behavioural data    -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we create equal sized Reaction Time (RT) bins per participant per conditions
% to check the alignment of the CPP with the actual median RT as the CPP
% peak should be highly linked with the reaction time.
parameters.binRTs(1, [.1 .3 .5 .7 .9], [1 2]); % as there are sometimes to little trials per all conditions, we need to resort to a smaller number.


% add the context
for indPP = 1:length(parameters.ppNames)
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,1) == 2, 4) = 3;
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,1) == 2 & parameters.behaviour{indPP}.trialMatrix(:,2) == 25, 5) = 3;
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,1) == 2 & parameters.behaviour{indPP}.trialMatrix(:,2) == 70, 5) = 4;

    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,1) == 1 & parameters.behaviour{indPP}.trialMatrix(:,2) == 70, 4) = 2;
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,1) == 1 & parameters.behaviour{indPP}.trialMatrix(:,2) == 25, 4) = 1;
    
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,1) == 1 & parameters.behaviour{indPP}.trialMatrix(:,2) == 70, 5) = 2;
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,1) == 1 & parameters.behaviour{indPP}.trialMatrix(:,2) == 25, 5) = 1;
end

parameters.conditions{4} = [1 2 3];
parameters.conditions{5} = [1 2 3 4];

%% %%%%%%%%%%%%%%%%%%%% Pre-set figure parameters    %%%%%%%%%%%%%%%%%%%%%%

parameters.figLayOut.letterSize  = 9;
parameters.figLayOut.letterType  = 'Arial';
parameters.figLayOut.lineWidth   = 1.2;
parameters.figLayOut.lineType    = {'-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-' '-'};
parameters.figLayOut.legNames{1} = {'Constant', 'Mixed'};
parameters.figLayOut.legNames{2} = {'Weak', 'Strong'}; % shift as they are 'order' by strength
parameters.figLayOut.legNames{3} = {'2', '4', '6', '8'};
parameters.figLayOut.legTitle{1} = 'Context';
parameters.figLayOut.legTitle{2} = 'Evidence strength';
parameters.figLayOut.legTitle{3} = 'ITI duration';

parameters.figLayOut.legNames{4}  = {'Weak', 'Strong', 'Mixed'};
parameters.figLayOut.legTitle{4} = 'Context';

parameters.figLayOut.legNames{5}  = {'Weak', 'Strong', 'Mixed (Weak)', 'Mixed (Strong)'};
parameters.figLayOut.legTitle{5} = 'Evidence strength';

parameters.figLayOut.plotCI      = 0.3; % get shaded area of not.
parameters.figLayOut.removeInter = 1;
parameters.figLayOut.saveDim     = [5 11];
parameters.figLayOut.colours     = ColoursCond12;
parameters.figLayOut.legends     = {'Weak', 'Strong', 'Mixed (Weak)', 'Mixed (Strong)'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------ Calculate DME ---------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unfortuantly, the dot coordinates aren't saved in Hannah presenation code
% and therefore we are using on/off timeseries as input to the modelfit.
% e.g. here instead of DME we are reshuffling the timeseries to be 0 when
% off and 1 when coherence motion. And add these in the behaviour file.

for indPP = 1:length(parameters.ppNames)
    fprintf('Calculate Momentary Evidence for %i out of %i.\n', indPP, length(parameters.ppNames))
    
    currOutput = fullfile(parameters.outputFolder, 'MomEvidence data', parameters.ppNames{indPP});
    if ~exist(currOutput, 'dir'), mkdir(currOutput); end
    
    
    
    if ~exist(fullfile(currOutput, 'MomEvidence.mat'), 'file')
        acc = 0;
        MomEvidence = [];
        
        for indBlock = 1:parameters.numBlocks(indPP)
            load(fullfile(parameters.inputFolder, parameters.ppNames{indPP}, 'Trial files', parameters.dataFiles{indPP}{indBlock}), 'par', 'coh');
            
            % get timecourse and resample --> weirdly save with a sample
            % rate of 15...
            trialCond = parameters.experiment{indPP}{indBlock}.ITI.*60;
            %trialCond(1) = trialCond(1);
            stim2frames= 60/(parameters.stim.duration/1000);
            maxLength = max(trialCond) + stim2frames;
            
            for indTrial = 1:length(trialCond)
                acc = acc + 1;
                
                MomEvidence(acc,:) = nan(1,maxLength);
                MomEvidence(acc,1:(stim2frames+trialCond(indTrial))) = [zeros(1, trialCond(indTrial) ,1)...
                    ones(1,stim2frames).* find(parameters.conditions{2} == parameters.behaviour{indPP}.trialMatrix(acc,2))];
            end
        end
        
        save(fullfile(currOutput, 'MomEvidence.mat'), 'MomEvidence');
        
        parameters.behaviour{indPP}.MomEvidence = MomEvidence;
    else
        load(fullfile(currOutput, 'MomEvidence.mat'), 'MomEvidence')
        parameters.behaviour{indPP}.MomEvidence = MomEvidence;
    end
end
clear MomEvidence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------ Behavioural data analysis ---------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TimeOnTask       = 0;
TimeOnExperiment = 0;
plotExperiment   = 0; % 1 is counts, 2 is percentages. This will be depending on which are used for the model.
fprintf('We started of with running the whole model to get an general understanding of the data.\n')
fprintf('Please go through the .log file that is created to see how the different conditions effect behaviour.\n')
% [rm, results] = parameters.doRANOVA(TimeOnTask, TimeOnExperiment, plotExperiment, [3 2 1]); % [x,'possible grouping', seperated plots]

fprintf('When it is run and inspected, we will rerun the code with the addition of plotting the behaviour.\n')
fprintf('Starting with the effect of context and evidence strength.\n')


parameters.figLayOut.colours = ColoursCond12;
% plotExperiment = 1;

% [rm, results] = parameters.doRANOVA(TimeOnTask, TimeOnExperiment, plotExperiment, [2 1]);
% parameters.plotRThistogram([3 1 2], 3);
% parameters.figLayOut.colours     = ColoursCond12([2 4],:);
% [rm, results] = parameters.doRANOVA(plotExperiment, [1 2]);

% parameters.figLayOut.colours     = ColoursCond12([2 4],:);
% [~,~, tbl] = parameters.doRANOVA(0, [4]);
%{
lme1 = fitglme(tbl, 'FAcount ~ 1 + Context', 'Distribution','Poisson',...
	'Link','log','FitMethod','Laplace', ...
	'DummyVarCoding','effects');
lme2 = fitglme(tbl, 'FAcount ~ 1 + Context + (1 | ppNames)', 'Distribution','Poisson',...
	'Link','log','FitMethod','Laplace', ...
	 'DummyVarCoding','effects');
compare(lme1, lme2, 'CheckNesting', true)
%}
% strongFA = tbl.FAcount(double(tbl.Context) == 2);
% [~, p_strongFA, diff] = permtest(strongFA, 0,10000);
% 
% lme_effectWeak = fitglme(tbl, 'FAcount ~ 1 + Context + (Context | ppNames)', 'Distribution','Poisson',...
%     'Link','log','FitMethod','Laplace');
% 
% tbl.Context = reordercats(tbl.Context, [2 1 3]);
% lme_effectStrong = fitglme(tbl, 'FAcount ~ 1 + Context + (Context | ppNames)', 'Distribution','Poisson',...
%     'Link','log','FitMethod','Laplace');
%compare(lme2, lme3, 'CheckNesting', true)


% [rm,~, tblAccu] = parameters.doRANOVA(plotExperiment, [1 2]);


%{
lme1 = fitlme(tblAccu, 'Accuracy ~ 1 + Context + Evidencestrength'%, 'Distribution','Binomial',...
    'FitMethod','Laplace', 'DummyVarCoding','effects');
lme2 = fitglme(tblAccu, 'Accuracy ~ 1 + Context * Evidencestrength', 'Distribution','Binomial',...
    'FitMethod','Laplace', 'DummyVarCoding','effects');
compare(lme1, lme2, 'CheckNesting', true)

lme3 = fitglme(tblAccu, 'Accuracy ~ 1 + Context + Evidencestrength + (1 | ppNames)', 'Distribution','Binomial',...
    'FitMethod','Laplace');
compare(lme1, lme3, 'CheckNesting', true)
%}
% lme_accuracy = fitglme(tblAccu, 'Accuracy ~ 1 + Context * Evidencestrength + (1 | ppNames) + (Evidencestrength | Context)',...
%     'Distribution','Binomial',...
%     'FitMethod','Laplace');
%}            

% plot ITI of context effect
%{

% first, we look at the possible effects of motion direction by looking
% % at possible interaction or main effects.
parameters.figLayOut.saveDim = [5 11];

figInfo = parameters.plotRThistogram([1 2]);
legend({'Weak Context', 'Strong Context', 'Mixed (Weak)', 'Mixed (Strong)'},...
    'Location', 'northeast');
plotSave(gcf, 'RTPlotForModel.png', figureFolder, [5 11]);

parameters.figLayOut.saveDim = [5 4];
% parameters.figLayOut.colours = [brewermap(2, 'Greens'); tmp([2 1],:)];
parameters.figLayOut.colours = parameters.figLayOut.colours([2 4],:);
multcompare(rm.FA,  'Evidencestrength', 'by', 'Context')

% parameters.figLayOut.colours = brewermap(4, 'Greys')-0.3;
% parameters.figLayOut.saveDim = [5 11];

% figInfo = parameters.plotRThistogram([4 3]);

parameters.figLayOut.saveDim = [5 4];
output = parameters.doRANOVA( 1, [3 5]); % [x,'possible grouping', seperated plots]
fprintf('Uncommand and use previous section to inspected direction effects.\n')

parameters.figLayOut.saveDim = [5 4];
parameters.figLayOut.colours = brewermap(4, 'Greys')-0.3;
output = parameters.doRANOVA(1, [3 2]); % [x,'possible grouping', seperated plots]
fprintf('Uncommand and use previous section to inspected direction effects.\n')
%}


%% Replot False alarms
% add third condition, compare the three context for FA plotting/
figureFolder = 'C:\Users\eng\Documents\2. Postdoc_Kelly\Project 5 - Neural correlates of static and dynamic decision-bound adjustments\Project 5.1 - Signatures of bound adjustment\Data\Figures\groupAverage\Behaviour\';

% behaviouralModel(parameters, timeSeries, [4 2], {'Weak context', 'Strong context', 'Weak mixed', 'Strong mixed'})
% false alarms are quite crucial for Hannah data, as we are looking into
% how a stable context is changing the 'bound'. Mixed condition is giving
% us about the same number of False alarms as the hard conditions.
%{
keyboard

% counter = 0;
for indPP = 1:length(parameters.ppNames)
    for indCond = 1:3
        posComb = unique(parameters.behaviour{indPP}.trialMatrix(:,4), 'rows');
        currCond  = all(parameters.behaviour{indPP}.trialMatrix(:, 4) == indCond,2);
        
        % then add a final bin that counts the false alarms:
        tmpFART  = parameters.behaviour{indPP}.indFalseAlarm(currCond,:);
        currFA   = parameters.behaviour{indPP}.FalseAlarm(currCond,:);
        
        qn(indCond,indPP) = sum(sum(currFA(:)))./ (sum(parameters.behaviour{indPP}.ITI(currCond))./2);
    end
end

FATable = array2table(qn');
withinTable = table(categorical([1:3]'), 'VariableNames', {'Context'});
rm.FA = fitrm(FATable, 'Var1-Var3 ~ 1', 'WithinDesign', withinTable);

spherTest.FA = mauchly(rm.FA);
results.FA   = ranova(rm.FA, 'WithinModel', 'Context');
multcompare(rm.FA, 'Context')
plotData = margmean(rm.FA, {'Context'});

FAcount = plotData.Mean; % and qn are the average number of trials in each bin
FAci   =  1.96*plotData.StdErr;
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
for indCond = 1:3
    bb = bar(indCond, FAcount(indCond));
    set(bb,'FaceColor', parameters.figLayOut.colours(indCond,:))
    errorbar(indCond, FAcount(indCond), FAci(indCond),'k','LineWidth',parameters.figLayOut.lineWidth)
end
fig2 = gca;
xtickangle( 45 )
% ylabel('False alarms');
% ylabel(' ');
% yticklabels('Proportion FA'');
xticks([1:indCond])
xticklabels({'Weak', 'Strong', 'Mixed'});
xlim([0 4])
ylim([0  0.2])
yticks(0:0.1:0.2)
% ylabel({'Proportion' 'False Alarm'})
% fig2.YAxis.Visible = 'off';   % remove y-axis

% xlabel('Evidence strength')
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'FA.png', figureFolder, [3.5 3.5]);

fprintf('Behaviour data has already been run.\n')
%}

%% FA and ITI length
%{
[within, between, modelFun, allTrialMatrix] = getConditions(parameters, [3 4]);

acc = 0;
for indITI = 2:2:8
    for indCond = 1:3
        acc = acc+1;
        forHistogram{acc} = [];  interFAInterval{acc} = []; numResponse{acc} = 0;
        for indPP = 1:length(parameters.ppNames)
            
            currCond  = parameters.behaviour{indPP}.trialMatrix(:, 4) == indCond & parameters.behaviour{indPP}.trialMatrix(:, 3) == indITI;
            possibilities = (sum(parameters.behaviour{indPP}.ITI(currCond))./2);
            % then add a final bin that counts the false alarms:
            tmpFART  = parameters.behaviour{indPP}.indFalseAlarm(currCond,:);
            currFA   = parameters.behaviour{indPP}.FalseAlarm(currCond,:);
            tmpFART(currFA == 0) = NaN;
            qn(acc,indPP) = sum(currFA(:))./possibilities;
            
            findForMisOrHit = find(currCond);
            findForMisOrHit = findForMisOrHit - 1; % find previous hit or mis
            findForMisOrHit(findForMisOrHit < 1) = [];
            
            qnHits(acc,indPP)   = sum(sum(currFA(~parameters.behaviour{indPP}.Misses(findForMisOrHit), :)))./possibilities;
            qnMisses(acc,indPP) = sum(sum(currFA(parameters.behaviour{indPP}.Misses(findForMisOrHit) == 1, :)))./possibilities;

            tmpInterval = diff([parameters.behaviour{indPP}.trialMatrix(currCond, 3)*-1 (tmpFART)],[], 2);
            interFAInterval{acc} = [interFAInterval{acc}; tmpInterval(~isnan(tmpInterval))];
            forHistogram{acc}    = [forHistogram{acc}; tmpFART(~isnan(tmpFART))];
            numResponse{acc}     = numResponse{acc} + (sum(parameters.behaviour{indPP}.ITI(parameters.behaviour{indPP}.trialMatrix(:, 4) == indCond))./2);
        end
    end
end
%
FATable = array2table(qn'.*100);
withinTable = table(repmat(categorical([1:3]'),4,1), repmat([2:2:8]', 3,1),'VariableNames', {'Context', 'ITIDuration'});
rm.FA = fitrm(FATable, 'Var1-Var12 ~ 1', 'WithinDesign', within.Table);

spherTest.FA = mauchly(rm.FA);
results.FA   = ranova(rm.FA, 'WithinModel', modelFun);
multcompare(rm.FA, 'ITIduration', 'by', 'Context')
plotData = margmean(rm.FA, {'Context',  'ITIduration'});

FAcount = plotData.Mean; % and qn are the average number of trials in each bin
FAci   =  1.96*plotData.StdErr;
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
for indCond = 1:3
    errorbar(2:2:8, FAcount(double(plotData.Context) == indCond)', FAci(double(plotData.Context) == indCond)',...
        'Color', parameters.figLayOut.colours(indCond,:),'LineWidth',parameters.figLayOut.lineWidth)
end
fig2 = gca;
ylabel('False alarms (%)');
xticks([2:2:8])
xticklabels({'2', '4', '6', '8'});
xlabel('ITI duration');

xlim([1.5 8.5])
ylim([0  0.2].*100)
yticks((0:0.05:0.2)*100)
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'FAITIall.png', figureFolder, [5 4]);
fprintf('Behaviour data has already been run.\n')

% hits
FATable = array2table(qnHits'.*100);
rm.FA = fitrm(FATable, 'Var1-Var12 ~ 1', 'WithinDesign', within.Table);

spherTest.FA = mauchly(rm.FA);
results.FA   = ranova(rm.FA, 'WithinModel', modelFun);
multcompare(rm.FA,  'Context', 'by', 'ITIduration')
plotData = margmean(rm.FA, {'Context',  'ITIduration'});

FAcount = plotData.Mean; % and qn are the average number of trials in each bin
FAci   =  1.96*plotData.StdErr;
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
for indCond = 1:3
    errorbar(2:2:8, FAcount(double(plotData.Context) == indCond)', FAci(double(plotData.Context) == indCond)',...
        'Color', parameters.figLayOut.colours(indCond,:),'LineWidth',parameters.figLayOut.lineWidth)
end
fig2 = gca;
ylabel('False alarms (%)');
xticks([2:2:8])
xticklabels({'2', '4', '6', '8'});
xlabel('ITI duration');

xlim([1.5 8.5])
ylim([0  0.2].*100)
yticks((0:0.05:0.2)*100)
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'FAITIHits.png', figureFolder, [5 4]);
fprintf('Behaviour data has already been run.\n')

% miss
FATable = array2table(qnMisses'.*100);
rm.FA = fitrm(FATable, 'Var1-Var12 ~ 1', 'WithinDesign', within.Table);

spherTest.FA = mauchly(rm.FA);
results.FA   = ranova(rm.FA, 'WithinModel',  modelFun);
multcompare(rm.FA,  'Context', 'by', 'ITIduration')
plotData = margmean(rm.FA, {'Context',  'ITIduration'});

FAcount = plotData.Mean; % and qn are the average number of trials in each bin
FAci   =  1.96*plotData.StdErr;
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
for indCond = 1:3
    errorbar(2:2:8, FAcount(double(plotData.Context) == indCond)', FAci(double(plotData.Context) == indCond)',...
        'Color', parameters.figLayOut.colours(indCond,:),'LineWidth',parameters.figLayOut.lineWidth)
end
fig2 = gca;
ylabel('False alarms (%)');
xlim([1.5 8.5])
xticks([2:2:8])
xticklabels({'2', '4', '6', '8'});
xlabel('ITI duration');

ylim([0  0.2].*100)
yticks((0:0.05:0.2)*100)
% xlabel('Evidence strength')
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
% leg = legend({'Weak', 'Strong', 'Mixed'});
% leg.Title.String = 'Difficulty Context';
plotSave(gca, 'FAITIMisse.png', figureFolder, [5 4]);
fprintf('Behaviour data has already been run.\n')
%}
%% False alarms occurent through the whole ITI
%{

parameters.figLayOut.colours = repmat(ColoursCond12([1 2 4 ],:),4,1); %
parameters.figLayOut.lineType = {'-', '--', '-.', ':'};

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
acc = 0;
timeSeries{1} = fliplr(-500:1000:1500).*-1;
timeSeries{2} = fliplr(-500:1000:3500).*-1;
timeSeries{3} = fliplr(-500:1000:5500).*-1;
timeSeries{4} = fliplr(-500:1000:7500).*-1;

for indCond = 1:3:length(forHistogram)
    acc = acc+1;
    x = hist(forHistogram{indCond}.*1000,timeSeries{acc});
    xgrid = linspace(timeSeries{acc}(1),timeSeries{acc}(end),length(x))';
    
    % cumHazard
       line(timeSeries{4}(1:length(xgrid)-1), (x(1:end-1)/numResponse{indCond}),...
        'Color', parameters.figLayOut.colours(indCond,:), 'LineWidth', 2,...
        'LineStyle', parameters.figLayOut.lineType{acc}); 
end

acc = 0;
for indCond = 2:3:length(forHistogram)
    acc = acc+1;
    x = hist(forHistogram{indCond}.*1000,timeSeries{acc});
    xgrid = linspace(timeSeries{acc}(1),timeSeries{acc}(end),length(x))';
    
    % cumHazard
       line(timeSeries{4}(1:length(xgrid)-1), (x(1:end-1)/numResponse{indCond}),...
        'Color', parameters.figLayOut.colours(indCond,:), 'LineWidth', 2,...
        'LineStyle', parameters.figLayOut.lineType{acc}); 
end

acc = 0;
for indCond = 3:3:length(forHistogram)
    acc = acc+1;
    x = hist(forHistogram{indCond}.*1000,timeSeries{acc});
    xgrid = linspace(timeSeries{acc}(1),timeSeries{acc}(end),length(x))';
    
    % cumHazard
       line(timeSeries{4}(1:length(xgrid)-1), (x(1:end-1)/numResponse{indCond}),...
        'Color', parameters.figLayOut.colours(indCond,:), 'LineWidth', 2,...
        'LineStyle', parameters.figLayOut.lineType{acc}); 
end

ylabel(' ');
xlim([-8000 350])
xticks([-8000:2000:500])
xticklabels([0:2000:8000]./1000);
xlabel(' ');

ylim([0 0.01]); yticks([0:0.005:0.01]); 
line([0 0], [0 0.01], 'Color', 'black', 'LineWidth', 1.5)
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'FAITIHist1.png', figureFolder, [4 6.5154]);
%}
%{
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
acc = 0;

for indCond = 2:3:length(forHistogram)
    acc = acc+1;
    x = hist(forHistogram{indCond}.*1000,timeSeries{acc});
    xgrid = linspace(timeSeries{acc}(1),timeSeries{acc}(end),length(x))';
    
   line(timeSeries{4}(1:length(xgrid)), (x/numResponse{indCond}*100),...
        'Color', parameters.figLayOut.colours(indCond,:), 'LineWidth', 2,...
        'LineStyle', parameters.figLayOut.lineType{acc});  

end

ylabel('False alarms');
xlim([-8000 350])
xticks([-8000:2000:500])
xticklabels([0:2000:8000]./1000);
xlabel(' ');

ylim([0 1.5]); yticks(0:0.5:1.5);

set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'FAITIHist2.png', figureFolder, [5 4]);

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
acc = 0;
for indCond = 3:3:length(forHistogram)
    acc = acc+1;
    x = hist(forHistogram{indCond}.*1000,timeSeries{acc});
    xgrid = linspace(timeSeries{acc}(1),timeSeries{acc}(end),length(x))';
    
    line(timeSeries{4}(1:length(xgrid)), (x/numResponse{indCond}*100),...
        'Color', parameters.figLayOut.colours(indCond,:), 'LineWidth', 2,...
        'LineStyle', parameters.figLayOut.lineType{acc}); 

end
ylabel('False alarms');
xlim([-8000 350])
xticks([-8000:2000:500])
xticklabels([0:2000:8000]./1000);
xlabel(' ');

ylim([0 1.5]); yticks(0:0.5:1.5);
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'FAITIHist3.png', figureFolder, [5 4]);

fprintf('Behaviour data has already been run.\n')
%{
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
acc = 0;
for indCond = 3:3:length(forHistogram)
    acc = acc + 1;
      x = hist(forHistogram{indCond}.*1000,timeSeries);%,'Normalization','pdf');
      xgrid  = linspace(timeSeries(1),0,length(x))';
      
     onlyLegend(acc) = line(xgrid,x/numResponse{indCond}.*100, 'Color','k', 'LineWidth', 2,...
         'LineStyle', parameters.figLayOut.lineType{acc});
end
ylabel('False alarms (%)');
xlim([-8000 0])
xticks([-8000:2000:0])
xticklabels([-8000:2000:0]./1000);
xlabel(' ');

ylim([0 50]); yticks(0:10:50);
% xlabel('Evidence strength')
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
leg = legend(parameters.figLayOut.legNames{3});
leg.Title.String = 'ITI duration (s)';
plotSave(gca, 'legend.png', figureFolder, [5 4]);
%}
%}
%% Inter false alarm time
%{
timeSeries{1} = (0:1000:2000);
timeSeries{2} = (0:1000:4000);
timeSeries{3} = (0:1000:6000);
timeSeries{4} = (0:1000:8000);

parameters.figLayOut.colours = repmat(parameters.figLayOut.colours(1:3,:),4,1);
parameters.figLayOut.lineType = {'-', '--', '-.', ':'};

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
acc = 0;
for indCond = 1:3:length(interFAInterval)

    acc = acc+1;
    x = hist(interFAInterval{indCond}.*1000,timeSeries{acc});
    
    xgrid  = linspace(timeSeries{acc}(1),timeSeries{acc}(end),length(x))';
    
%     ecdf(gca,xgrid','Censoring',x);%,'function','survivor');
     line(xgrid,(x/numResponse{indCond}.*100), 'Color', parameters.figLayOut.colours(indCond,:), 'LineWidth', 2,...
         'LineStyle', parameters.figLayOut.lineType{acc});
end

ylabel('False alarms');
xlim([0 8000])
xticklabels([0:1000:8000]./1000);
xlabel(' ');

ylim([0 1]); yticks(0:0.2:1);
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'interFAITIHist1.png', figureFolder, [5 4]);


figure('units','normalized','outerposition',[0 0 1 1]); hold on;
acc = 0;

for indCond = 2:3:length(interFAInterval)
    acc = acc + 1;
   

    x = hist(interFAInterval{indCond}.*1000,timeSeries{acc});
    
    xgrid  = linspace(timeSeries{acc}(1),timeSeries{acc}(end),length(x))';
    
%     ecdf(gca,xgrid','Censoring',x);
     onlyLegend(acc) = line(xgrid,(x/numResponse{indCond}.*100), 'Color', parameters.figLayOut.colours(indCond,:), 'LineWidth', 2,...
         'LineStyle', parameters.figLayOut.lineType{acc});
end
ylabel('False alarms');
xlim([0 8000])
xticklabels([0:2000:8000]./1000);
xlabel(' ');
xlabel(' ');

ylim([0 1]); yticks(0:0.2:1);
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'interFAITIHist2.png', figureFolder, [5 4]);

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
acc = 0;
for indCond = 3:3:length(interFAInterval)
    acc = acc+1;
  
    x = hist(interFAInterval{indCond}.*1000,timeSeries{acc});
    xgrid  = linspace(timeSeries{acc}(1),timeSeries{acc}(end),length(x))';
    
%     ecdf(gca,xgrid','Censoring',x,'function','survivor');
     onlyLegend(acc) = line(xgrid,(x/numResponse{indCond}.*100), 'Color', parameters.figLayOut.colours(indCond,:), 'LineWidth', 2,...
         'LineStyle', parameters.figLayOut.lineType{acc});
end

ylabel('False alarms');
xlim([0 8000])
xticklabels([0:2000:8000]./1000);
xlabel(' ');
xlabel(' ');

ylim([0 1]); yticks(0:0.2:1);
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'interFAITIHist3.png', figureFolder, [5 4]);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------ Behavioural modelling  ------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------   preset everything for modelling -------------------
% 
parameters.stim.epoch = parameters.eeg.epoch;
parameters.stim.epochPlot = parameters.eeg.epochPlot;
parameters.stim.targetEpoch = parameters.eeg.targetEpoch;
parameters.stim.responseEpoch = parameters.eeg.responseEpoch;
parameters.stim.baselineCorrect = 0 + [-1 1]*1000/parameters.stim.freqSSVEP; % note if you add a second row, these will be used for FA plots. Only nesc. when you are using target-locked topos

parModelling = dataModelling.getInstance(parameters); clear parameters

parModelling.modelBehaviour.ChiOrG = 2;
parModelling.figLayOut.colours  = ColoursCond12;
parModelling.figLayOut.lineWidth   = 1.2;
parModelling.figLayOut.targetLim   = [-200 0:500:1000];
parModelling.figLayOut.responseLim = [-400 0 200];

parModelling.modelBehaviour.cohs = parModelling.conditions{2}; % for transferability we have to indix the coherences level as drift specificly would be depending on this.

parModelling.plotRTquantiles([4 2]); % 1 is taking the count, 2 is taking the proportions
legend({'Weak Context', 'Strong Context', 'Mixed (Weak)', 'Mixed (Strong)'},...
    'Location', 'southeast');
plotSave(gcf, 'RTPlotForModel.png', fullfile(parModelling.figFolder,  'groupAverage/Modelling/'), [5 11]);

% validation compared to Simon code
%{
% datsumFA = parameters.plotRTquantiles(4, ChiOrG); % 1 is taking the count, 2 is taking the proportions

% parameters.modelBehaviour.datsum.pij(end-1:end,end) = NaN;
% parameters.modelBehaviour.datsum.pij(end-1,end) = datsumFA.pij(end,end);
% parameters.modelBehaviour.datsum.nShortestITI = datsumFA.nShortestITI;
%}


% These will allow bound to be a free parameter (which might be strategically adjusted) and within-trial noise
% set as the scaling parameter, set at s=0.1 for all conditions - this is effectively the noise in the sensory evidence

% There are many ways to set up the drift rate parameters and I think all of them will have to be run and the best
% chosen for each model structure based on AIC/BIC
% 1) one drift rate parameter for W, with 70/25* that value assumed for S, same in M. These I'm calling '1d'
% 2a) one drift rate parameter for W, with 70/25* that value assumed for S, and a drift rate boost parameter k for the M condition. These I'm calling '1d_1k'
% 2b) two drift rate parameters, one for W one for S, assumed same in M. Called '2d'
% 3a) 3 free d.r. params: one for W, one for S, and a scaled *k version of those in M. Called '2d_1k'
% 3b) 3 free d.r. params: one for W, one for S, and one for loh coh in M, with M hicoh scaled *70/25. Called '3d'
% 4) four free parameters for the 4 different conditions - no constraints at all (could be too flexible). called '4d'
% I might be missing some other possible ways??

% fminsearch (SIMPLEX) Settings:
parModelling.options = optimset('fminsearch');
parModelling.options.MaxIter = 500; % maximum iterations
parModelling.options.TolFun = 1;     % tolerance - if changes in the objective function are consistently less than this on lots of consecutive evaluations, then it's probably as good as you're going to get, so fminsearch will quit and move on. Don't know why set at 5, presumably trial and error, but this is worth trying to reduce in the refinement step to make sure it doens't stop until it really has a good fit
parModelling.options.TolX = 1e-2;

StopEarly = @(x,optimvalues,state) optimvalues.fval>4000;
parModelling.options = optimset(parModelling.options, 'Display', 'iter', 'OutputFcn', StopEarly); % @stopEarly give us an oppertunity to prematurely stop if something has more than 100 interation > 1000


parModelling.modelBehaviour.noiseSTD  = 0.1;
parModelling.modelBehaviour.learnBoost = 3;
parModelling.modelBehaviour.simulateMoreX = 1; % it can be good to simulate more trials than the subjects actually did because then at least the simulated data are reliable.

repeats = 1000;
parModelling.figLayOut.saveDim     = [5 11];

%% 1) 3 bounds, 1 drift, 1 leak and 1 tnd
% First a model ONLY allowing Bound-adjustment
% same non-decision time across conds (should we try freeing that?) and the same leak
%

ModelName = '3Bound1Drift1Leak1tnd';
modelPar = table(char(repmat('bound',3,1), 'drift', 'leak', 'tnd'),...
    [repmat(0.1,3,1); 0;  0.001; 0.05],...
    [repmat(2,3,1); 0.4; 1; 0.35],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

[Params_3Bound, err] = parModelling.applyModelling(ModelName,  [], modelPar, 1000);
% [pred, AIC, BIC]= parModelling.plotModellingResults(Params_3Bound(1,:), ModelName, modelPar.Names);


% refine
parModelling.options.TolFun = 1e-4;
parModelling.options.MaxIter = 1000; % maximum iterations
parModelling.modelBehaviour.simulateMoreX = 5;
% [Params_3Bound, Err3Bound] = parModelling.applyModelling(ModelName, Params_3Bound(1:50,:), modelPar);

% [pred, output] =  parModelling.plotModellingResults(Params_3Bound(1,:), ModelName, modelPar.Names);

% adaptedFAPlot(parModelling, pred, ModelName);
% parModelling.plotParameters([Params_3Bound(1:60,[6 1:3])],ModelName,...
%     {'tnd', 'Bound weak', 'Bound strong', 'Bound mixed'},  [modelPar.Lower([6 1:3],:) modelPar.Upper([6 1:3],:)]);


%% 2) 1 bounds, 1 drift, 3 leak and 1 tnd
% Next model: same but now bound is fixed and LEAK is adjustable across the 3 regimes instead.
parModelling.options.MaxIter = 500; % maximum iterations
parModelling.options.TolFun = 1;    % tolerance - if changes in the objective function are consistently less than this on lots of consecutive evaluations, then it's probably as good as you're going to get, so fminsearch will quit and move on. Don't know why set at 5, presumably trial and error, but this is worth trying to reduce in the refinement step to make sure it doens't stop until it really has a good fit
parModelling.modelBehaviour.simulateMoreX = 1;

ModelName = '1Bound1Drift3Leak1tnd';
modelPar = table(char('bound', 'drift', repmat('leak',3,1), 'tnd'),...
    [0.1; 0;  repmat(0.001,3,1); 0.05],...
    [2; 0.4; ones(3,1); 0.35],...
    'VariableNames', {'Names', 'Lower', 'Upper'});
Params_3Leak = parModelling.applyModelling(ModelName, [], modelPar, 1000);

% refine
parModelling.options.TolFun = 1e-2;
parModelling.modelBehaviour.simulateMoreX = 5;

% [Params_3Leak, Err3Leak] = parModelling.applyModelling(ModelName, Params_3Leak(1:50,:), modelPar);
% keyboard
% [pred, output] = parModelling.plotModellingResults(Params_3Leak(1,:), ModelName, modelPar.Names);
% adaptedFAPlot(parModelling, pred, ModelName);
% parModelling.plotParameters([Params_3Leak(:,[6 3:5])], [ModelName],...
%     {'tnd', 'Leak weak', 'Leak strong', 'Leak mixed'},  [modelPar.Lower([6 3:5],:) modelPar.Upper([6 3:5],:)]);

%% The crude fit above will do for starting any of these models with
% various drift rate set-ups - saves time

close all
clc
repeats = 200;
parModelling.options.TolFun = 1e-2;
parModelling.modelBehaviour.simulateMoreX = 5;

%% 3) refine the bound models
% 3.1) introducing a fit 2 drifts due to partice in mixed condition   
%{
ModelName = '3Bound2Drift1Leak1tnd';
fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)
modelPar = table(char(repmat('bound',3,1), 'leak', 'tnd', repmat('drift',2,1)),...
    [repmat(0.1,3,1); 0.001; 0.05;0; 0],...
    [repmat(2,3,1); 1; 0.35;0.4; 0.4],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

[Params_3Bound2Drift, Err] = parModelling.applyModelling(ModelName, Params_3Bound(1:20,[1:3 5:end]),modelPar, repeats, [Params_3Bound(1:20,4) Params_3Bound(1:20,4)*70/25], 6:7);

pred = parModelling.plotModellingResults(Params_3Bound2Drift(1,:), ModelName, modelPar.Names);
adaptedFAPlot(parModelling, pred, ModelName);
%}

% 3.2) individually determine the drifts, assuming same for contast and
% mixed conditions full fit 
%{
close all
clc
try
    ModelName = '3Bound2Drift1Leak1tnd1Boost';
    fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)
    fittedDrift = Params_3Bound2Drift(1:20,6:7);
    modelPar = table(char(repmat('bound',3,1), 'leak', 'tnd', repmat('drift',2,1), 'boost'), ...
        [repmat(0.1,3,1); 0.001; 0.05; 0; 0; 0.8],...
        [repmat(2,3,1); 1; 0.35; 0.4; 0.4; 1.5],...
        'VariableNames', {'Names', 'Lower', 'Upper'});
    
    [Params_3Bound2Drift1Boost, Err] = parModelling.applyModelling(ModelName,Params_3Bound2Drift(1:20,:), modelPar, repeats, [], 8);
    
    
%     pred = parModelling.plotModellingResults(Params_3Bound2Drift1Boost(1,:), ModelName, modelPar.Names);
%     adaptedFAPlot(parModelling, pred, ModelName);
end
%}

% 3.3) individually determine all drifts to account for the boost (as there are different
% effect of easy vs. hard)
%

ModelName = '3Bound4Drift1Leak1tnd1';

modelPar = table(char(repmat('bound',3,1), 'leak', 'tnd', repmat('drift',4,1)),...
    [repmat(0.1,3,1);   0.001;  0.05;   0;      0;      0;      0],...
    [repmat(2,3,1);     1;      0.35;   0.4;    0.4;    0.4;    0.4],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

[Params_3Bound4Drift, Err] = parModelling.applyModelling(ModelName, Params_3Bound(1:20,[1:3 5:end]),modelPar, repeats,...
    [Params_3Bound(1:20,4) Params_3Bound(1:20,4)*70/25 Params_3Bound(1:20,4) Params_3Bound(1:20,4)*70/25], 6:9);

% smaller: [4 6]
% smaller: [obj.figLayOut.saveDim(1) obj.figLayOut.saveDim(2)./2]

% [pred, outputParameters.Bound] = parModelling.plotModellingResults(Params_3Bound4Drift(1,:), ModelName, modelPar.Names);
% adaptedFAPlot(parModelling, pred, ModelName);


% quickPlotPredictedBounds(parModelling, Params_3Bound4Drift(1,1:3))
% createExamplePlot(parModelling,ModelName, Params_3Bound4Drift(1,:), modelPar.Names, [3 6])
%}
% 3.4) individually determine all drifts to account for the boost (as there are different
% effect of easy vs. hard)
%{
ModelName = '3Bound2Drift2Boost1Leak1tnd1';

modelPar = table(char(repmat('bound',3,1), 'leak', 'tnd',  repmat('boost',2,1), repmat('drift',2,1)),...
    [repmat(0.1,3,1);   0.001;  0.05;   1;	1;	0;      0],...
    [repmat(2,3,1);     1;      0.35;   1.6;  1.6;  0.4;    0.4;   ],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

[Params_3Bound4Drift2Boost, Err] = parModelling.applyModelling(ModelName, Params_3Bound(1:20,[1:3 5:end]),modelPar, repeats,...
    [Params_3Bound(1:20,4) Params_3Bound(1:20,4)*70/25], [6 7; 8 9]);

% [pred, GoF] = parModelling.plotModellingResults(Params_3Bound4Drift2Boost(1,:), ModelName, modelPar.Names);
% adaptedFAPlot(parModelling, pred, ModelName);
% quickPlotPredictedBounds(parModelling, Params_3Bound4Drift(1,1:3))

% createExamplePlot(parModelling, Params_3Bound4Drift(1,:), modelPar.Names)
clear Params_3Bound*

clear Params_3Bound*
%}
%% 4) refine the Leak models
%{
close all
clc
% 4.1) introducing a fit 2 drifts due to partice in mixed condition
ModelName = '1Bound2Drift3Leak1tnd';
fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)
modelPar = table(char('bound', repmat('leak',3,1), 'tnd', repmat('drift',2,1)),...
    [0.1; repmat(0.001,3,1); 0.05; 0; 0],...
    [2; ones(3,1); 0.35; 0.4; 0.4],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

[Params_3Leak2Drift, Err] = parModelling.applyModelling(ModelName, Params_3Leak(1:20,[1 3:end]), modelPar, repeats,...
    [Params_3Leak(1:20, 2), Params_3Leak(1:20, 2)*70/25], 6:7);

% pred = parModelling.plotModellingResults(Params_3Leak2Drift(1,:), ModelName, modelPar.Names);
% adaptedFAPlot(parModelling, pred, ModelName);


% 4.2) individually determine the drifts, assuming same for contast and
% mixed conditions
%{
try
    ModelName = '1Bound2Drift3Leak1tnd1Boost';
    fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)
    fittedDrift = Params_3Leak2Drift(1:20,6:7);

    modelPar = table(char('bound', repmat('leak',3,1), 'tnd', repmat('drift',2,1), 'boost'),...
        [0.1; repmat(0.001,3,1); 0.05; nanmean(fittedDrift(:,1))  - std(fittedDrift(:,1)).*4; nanmean(fittedDrift(:,2)) - std(fittedDrift(:,2)).*4; 0.8],...
        [2; repmat(1,3,1); 0.35; nanmean(fittedDrift(:,1)) + std(fittedDrift(:,1)).*4; nanmean(fittedDrift(:,2)) + std(fittedDrift(:,2)).*4; 1.5],...
        'VariableNames', {'Names', 'Lower', 'Upper'});
    
    [Params_3Leak2Drift1boost,err] = parModelling.applyModelling(ModelName, Params_3Leak2Drift(1:20,:), modelPar, repeats, 8);
    
    
    
    pred = parModelling.plotModellingResults(Params_3Leak2Drift1boost(1,:), ModelName, modelPar.Names);
    adaptedFAPlot(parModelling, pred, ModelName);
end
close all
clc
%}

% 3.3) individually determine all drifts to account for the boost (as there are different
% effect of easy vs. hard)

ModelName = '1Bound4Drift3Leak1tnd1';
fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)
fittedDrift = Params_3Leak2Drift(1:20,2);

modelPar = table(char(repmat('leak',3,1), 'bound', 'tnd', repmat('drift',4,1)), ...
    [repmat(0.001,3,1); 0; 0.05; 0; 0; 0; 0],...
    [ones(3,1); 2; 0.35; 0.4; 0.4; 0.4; 0.4],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

[Params_3Leak4Drift, Err] = parModelling.applyModelling(ModelName, Params_3Leak(1:20,[1 3:end]), modelPar, repeats,...
    [Params_3Leak(1:20, 2), Params_3Leak(1:20, 2)*70/25 Params_3Leak(1:20, 2), Params_3Leak(1:20, 2)*70/25], 6:9);

% pred = parModelling.plotModellingResults(Params_3Leak4Drift1(1,:), ModelName, modelPar.Names);
% adaptedFAPlot(parModelling, pred, ModelName);

clear Params_3Leak*
%}

%% 5) Now neurally-constrained model with urgency signal determined by beta in long ITI:
% 5.1) model with reflective bound and fitting 3 criteria
% Finally, it might make sense to go simpler with one drift rate parameter 
% first, to verify model technically works before going to more drift rates...
clc
% close 
repeats = 1000;
StopEarly = @(x,optimvalues,state) optimvalues.fval>5000;
parModelling.options = optimset(parModelling.options, 'Display', 'iter', 'OutputFcn', StopEarly); % @stopEarly give us an oppertunity to prematurely stop if something has more than 100 interation > 1000

load('urg.mat')
applyReflectiveBound = 1;

ModelName = 'NI_3Criteria1Drift1Tnd1stdNoise';
fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)

modelPar = table(char(repmat('criteria',3,1), 'drift', 'tnd', 'trialnoise'),...
    [zeros(3,1); 0; 0.05; 0.001],...
    [repmat(0.2,3,1); 0.4; 0.35; 2],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

parModelling.modelBehaviour.noiseSTD  = 1;
parModelling.options.TolFun = 1;
parModelling.options.MaxIter = 500; % maximum iterations
parModelling.modelBehaviour.simulateMoreX = 1;
[Params_NI3Criteria, err] = applyNIModelling(parModelling, ModelName, urg, applyReflectiveBound, [], modelPar, repeats);

% refine NI 3 criteria model
parModelling.options.TolFun = 1e-4;
parModelling.options.MaxIter = 1000; % maximum iterations
parModelling.modelBehaviour.simulateMoreX = 5;
[Params_NI3Criteria, err]= applyNIModelling(parModelling, ModelName, urg, applyReflectiveBound, Params_NI3Criteria(1:25,:), modelPar);

% noise = parModelling.modelBehaviour.noiseSTD*randn(length(parModelling.modelBehaviour.datsum.TOW), max(parModelling.modelBehaviour.datsum.maxln), parModelling.modelBehaviour.simulateMoreX); % note the 0.1* because assuming s = 0.1. Not done for NI models below
% [err, simdat, ~, DV] = NIDecisionModels(parModelling, Params_NI3Criteria(1,:),  modelPar.Names, urg, noise, 1, applyReflectiveBound);

% pred = parModelling.plotModellingResults(Params_NI3Criteria(1,:), ModelName, modelPar.Names, err, simdat, DV);
% adaptedFAPlot(parModelling, pred, ModelName);


numRepeats = 200;

% 5.2) model with reflective bound and fitting 3 criteria and 2 drifts
% Finally, it might make sense to go simpler with one drift rate parameter first, to verify model technically works
% before going to more drift rates...
%{
clc
ModelName = 'NI_3Criteria2Drift1Tnd1stdNoise';

fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)
try    
    fittedDrift = Params_NI3Criteria(1:5,2);

    modelPar = table(char(repmat('criteria',3,1), 'tnd', 'trialnoise', repmat('drift',2,1)),...
        [repmat(0,3,1); 0.05; 0.001; nanmean(fittedDrift)  - std(fittedDrift).*4; nanmean(fittedDrift).*(70/25)  - std(fittedDrift).*4],...
        [repmat(0.2,3,1); 0.35; 2; nanmean(fittedDrift)  + std(fittedDrift).*4; nanmean(fittedDrift).*(70/25)  + std(fittedDrift).*4],...
        'VariableNames', {'Names', 'Lower', 'Upper'});

    parModelling.options.TolFun = 1e-4;
    parModelling.options.MaxIter = 1000; % maximum iterations
    parModelling.modelBehaviour.simulateMoreX = 5;   
    [Params_NI3Criteria2drift,err] = applyNIModelling(parModelling, ModelName, Params_NI3Criteria(1:5, [1:3 5:end]), modelPar, 200, 6:7, urg,1);
   
    noise = parModelling.modelBehaviour.noiseSTD*randn(length(parModelling.modelBehaviour.datsum.TOW), max(parModelling.modelBehaviour.datsum.maxln), parModelling.modelBehaviour.simulateMoreX); % note the 0.1* because assuming s = 0.1. Not done for NI models below
    [err, simdat, ~, DV] = NIDecisionModels(parModelling, Params_NI3Criteria2drift(1,:),  modelPar.Names, urg, noise,1,1);
    
    pred = parModelling.plotModellingResults(Params_NI3Criteria2drift(1,:), ModelName, modelPar.Names, err, simdat, DV);
    adaptedFAPlot(parModelling, pred, ModelName);

catch
    fprintf('didnt work  %s \n', ModelName)
end
%}

% 5.3) model with reflective bound and fitting 3 criteria and 2 drifts and
% boost
% Finally, it might make sense to go simpler with one drift rate parameter first, to verify model technically works
% before going to more drift rates...
%{
clc
ModelName = 'NI_3Criteria2Drift1Tnd1stdNoiseBoost';

fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)
try
    fittedDrift = Params_NI3Criteria2drift(1:20,6:7);

    modelPar = table(char(repmat('criteria',3,1), 'tnd', 'trialnoise', repmat('drift',2,1), 'boost'),...
        [repmat(0,3,1); 0.05; 0.001;  nanmean(fittedDrift(:,1))  - std(fittedDrift(:,1)).*4; nanmean(fittedDrift(:,2))  - std(fittedDrift(:,2)).*4;0.8],...
        [repmat(0.2,3,1); 0.35; 2;  nanmean(fittedDrift(:,1))  + std(fittedDrift(:,1)).*4; nanmean(fittedDrift(:,2))  + std(fittedDrift(:,2)).*4; 1.5],...
        'VariableNames', {'Names', 'Lower', 'Upper'});
    
    parModelling.options.TolFun = 1e-4;
    parModelling.options.MaxIter = 1000; % maximum iterations
    parModelling.modelBehaviour.simulateMoreX = 5;   
    [Params_NI3Criteria2driftboost,err] = applyNIModelling(parModelling, ModelName, Params_NI3Criteria2drift, modelPar, 200, 8, urg,1);

    noise = parModelling.modelBehaviour.noiseSTD*randn(length(parModelling.modelBehaviour.datsum.TOW), max(parModelling.modelBehaviour.datsum.maxln), parModelling.modelBehaviour.simulateMoreX); % note the 0.1* because assuming s = 0.1. Not done for NI models below
    [err, simdat, ~, DV] = NIDecisionModels(parModelling, Params_NI3Criteria2driftboost(1,:),  modelPar.Names, urg, noise,1,1);
    
    pred = parModelling.plotModellingResults(Params_NI3Criteria2driftboost(1,:), ModelName, modelPar.Names, err, simdat, DV);
    adaptedFAPlot(parModelling, pred, ModelName);

catch
    fprintf('didnt work  %s \n', ModelName)
end
%}

% 5.4) model with reflecting bound with 3 criteria and 4 drifts
%
ModelName = 'NI_3Criteria4Drift1Tnd1stdNoise';
fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)

modelPar = table(char(repmat('criteria',3,1), 'tnd', 'trialnoise', repmat('drift',4,1)),...
    [zeros(3,1); 0.05; 0.001; 0; 0; 0; 0],...
    [repmat(0.2,3,1); 0.35; 2; 0.4; 0.4; 0.4; 0.4],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

[Params_NI3Criteria4drift, err] = applyNIModelling(parModelling, ModelName, urg, applyReflectiveBound,  Params_NI3Criteria(:, [1:3 5:end]), modelPar, numRepeats,...
   [Params_NI3Criteria(:,4) Params_NI3Criteria(:,4)*70/25 Params_NI3Criteria(:,4) Params_NI3Criteria(:,4)*70/25], 6:9);


noise = parModelling.modelBehaviour.noiseSTD*randn(length(parModelling.modelBehaviour.datsum.TOW), max(parModelling.modelBehaviour.datsum.maxln), parModelling.modelBehaviour.simulateMoreX); % note the 0.1* because assuming s = 0.1. Not done for NI models below
[err, simdat, ~, DV] = NIDecisionModels(parModelling, Params_NI3Criteria4drift(1,:),  modelPar.Names, urg, noise,1,1);
 

keyboard
parModelling.figLayOut.saveDim = [4 6];
[pred, outputParameters.Criteria]  = parModelling.plotModellingResults(Params_NI3Criteria4drift(1,:), ModelName, modelPar.Names, err, simdat, DV);
% adaptedFAPlot(parModelling, pred, ModelName);

% createExamplePlot(parModelling,ModelName,Params_NI3Criteria4drift(1,:), modelPar.Names, [1.4 2.5])


%%  5) Now neurally-constrained model with urgency signal determined by beta in long ITI:
% Finally, it might make sense to go simpler with one drift rate parameter 
% first, to verify model technically works before going to more drift
% rates...
% There's lots more neurally-informed models to try after this - other drift rate setups as above, but also adding back
% in the motor fluctuation parameter sz somehow (e.g. what if we just set it so that motor fluct have the same overall
% s.d. as the within trial sensory evidence noise s, and so we would still only have one free parameter) and it will be
% important to try the leak adjustment model with the neurally-constrained urgency.

% 6.1) model fitting 3 leak
% Finally, it might make sense to go simpler with one drift rate parameter first, to verify model technically works
% before going to more drift rates...
%
clc

applyReflectiveBound = 0;
numRepeats = 1000;
ModelName = 'NI_3Leak1Drift1Tnd1stdNoise';

StopEarly = @(x,optimvalues,state) optimvalues.fval>4000;
parModelling.options = optimset(parModelling.options, 'Display', 'iter', 'OutputFcn', StopEarly); % @stopEarly give us an oppertunity to prematurely stop if something has more than 100 interation > 1000

fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)

modelPar = table(char(repmat('Leak',3,1), 'drift', 'tnd', 'trialnoise'),...
    [repmat(0.001,3,1); 0; 0.05; 0.001],...
    [ones(3,1); 0.4; 0.35; 0.4],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

parModelling.modelBehaviour.noiseSTD  = 1;
parModelling.options.TolFun = 1;
parModelling.options.MaxIter = 500; % maximum iterations
parModelling.modelBehaviour.simulateMoreX = 1;
[Params_NI3Leak, err] = applyNIModelling(parModelling, ModelName, urg, applyReflectiveBound, [], modelPar, numRepeats);

parModelling.options.TolFun = 1e-4;
parModelling.options.MaxIter = 1000; % maximum iterations
parModelling.modelBehaviour.simulateMoreX = 5;
[Params_NI3Leak, err]= applyNIModelling(parModelling, ModelName, urg, applyReflectiveBound, Params_NI3Leak(1:25,:), modelPar);

% noise = parModelling.modelBehaviour.noiseSTD*randn(length(parModelling.modelBehaviour.datsum.TOW), max(parModelling.modelBehaviour.datsum.maxln), parModelling.modelBehaviour.simulateMoreX); % note the 0.1* because assuming s = 0.1. Not done for NI models below
% [err, simdat, ~, DV] = NIDecisionModels(parModelling, Params_NI3Leak(1,:),  modelPar.Names, urg, noise,1);
% [pred, A, B] = parModelling.plotModellingResults(Params_NI3Leak(1,:), ModelName, modelPar.Names, err, simdat, DV);
% adaptedFAPlot(parModelling, pred, ModelName);

% parModelling.plotParameters([Params_NI3Leak(1:100,[6 1:3])],ModelName,...
%     {'trial noise', 'Leak weak', 'Leak strong', 'Leak mixed'},  [modelPar.Lower([4 1:3],:) modelPar.Upper([4 1:3],:)]);


% 6.4 model with reflecting bound with 3 leak and 4 drifts
numRepeats = 200;
ModelName = 'NI_3Leak4Drift1Tnd1stdNoise';
fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)

modelPar = table(char(repmat('Leak',3,1), 'tnd', 'trialnoise', repmat('drift',4,1)),...
    [repmat(0.001,3,1); 0.05; 0.001; 0; 0; 0; 0],...
    [ones(3,1); 0.35; 2; 0.4; 0.4; 0.4; 0.4],...
    'VariableNames', {'Names', 'Lower', 'Upper'});
[Params_NI3Leak4drift, err]= applyNIModelling(parModelling, ModelName,urg, applyReflectiveBound,  Params_NI3Leak(:, [1:3 5:end]), modelPar, numRepeats,...
    [Params_NI3Leak(:,4) Params_NI3Leak(:,4)*70/25 Params_NI3Leak(:,4) Params_NI3Leak(:,4)*70/25], 6:9);

 
noise = parModelling.modelBehaviour.noiseSTD*randn(length(parModelling.modelBehaviour.datsum.TOW), max(parModelling.modelBehaviour.datsum.maxln), parModelling.modelBehaviour.simulateMoreX); % note the 0.1* because assuming s = 0.1. Not done for NI models below
[err, simdat, ~, DV] = NIDecisionModels(parModelling, Params_NI3Leak4drift(1,:),  modelPar.Names, urg, noise,1);

[pred, outputParameters.Leak] = parModelling.plotModellingResults(Params_NI3Leak4drift(1,:), ModelName, modelPar.Names, err, simdat, DV);
adaptedFAPlot(parModelling, pred, ModelName);
% createExamplePlot(parModelling,ModelName,Params_NI3Criteria4drift(1,:), modelPar.Names, [1.4 2.5])

% parModelling.plotParameters([Params_NI3Leak4drift(1:100,[2 1:3])],ModelName,...
%     {'tnd', 'Leak weak', 'Leak strong', 'Leak mixed'},  [modelPar.Lower([2 1:3],:) modelPar.Upper([2 1:3],:)]);
