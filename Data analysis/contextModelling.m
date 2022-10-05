%% Behaviorual and neurally-informed modelling
% Data collected Hannah Craddock for her Master's thesis, 2016:
%   'Electrophysiological and behavioural indices of decision criterion
%   adjustments across context of weak and strong evidence.'
%
% Current code is written by A.C. Geuzebroek and part of:
%   Geuzebroek AC, Craddock H, O’Connell RG, & Kelly SP (2022).
%   Balancing true and false detection of intermittent sensory targets
%   by adjusting the inputs to the evidence accumulation process 
%   https://biorxiv.org/cgi/content/short/2022.09.01.505650v1)
%
% In short, participants were asked to continous monitor a cloud of randomly moving
% dots for intermittered target defined as upwards coherent moving dots.
% Afterwhich they were asked to response as fast and accurate as possible.
% Difficulty context was manipulated with:
%   1) Hard  (25% motion coherence)
%   2) Easy  (70% motion coherence)
%   3) Mixed (25% and 70% motion coherence with equal probabiltiy).
%
% Add this data in the follow way:
%   InputFolder     = BASEFOLDER\Data\Processed\
%       Participant folder = BASEFOLDER\Data\Processed\PPNAME\ (note never use initials!)
%           EEG folder           = BASEFOLDER\Data\Processed\EEG data\PPNAME\
%
% costum-made code:
%   1) dataAnalysis
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
    %error('Please set your data folder')
    inputFolder = 'C:\Users\eng\Documents\2. Postdoc_Kelly\Project 5 - Neural correlates of static and dynamic decision-bound adjustments\Project 5.1 - Signatures of bound adjustment\Data';
    % inputFolder = 'YOURFOLDER\Data';    % add base folder of data. 
    addpath(genpath('C:\Users\eng\Downloads\Context-Dependent-Detection-main\Neurally-Informed-Modelling-main\'))
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
% As we use preprocessed data here, we can put everything on 0!
parameters.system = 'bioSemi';   % Object will choice the approtiated loading function of EEGLAB
parameters.analysisEEG      = 0; % Analyse EEG (obviously we want to look at this)

parameters.analysisEyelink  = 0;
parameters.analysisEOG      = 0; % Use additionally EOG electrodes and Frontal electrodes
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
parameters.getInformation(prepData,inputFolder, conditions, condNames, blockRandomized);
clear cond* *Folder prepData block*


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
% all parameters from EEG to use for NI modelling, or actually mostly to
% plot the simulated decision-variable.

%---------------     set some EEG parameters  ----------------------
parameters.eeg.SampleRate   = 512;
parameters.eeg.HPFcutoff    = 0;
parameters.eeg.baseline     = 0 + [-1 0]*1/parameters.stim.freqSSVEP;  

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------------   Behavioural data    -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we create equal sized Reaction Time (RT) bins per participant per conditions
% to check the alignment of the CPP with the actual median RT as the CPP
% peak should be highly linked with the reaction time.
parameters.binRTs(1, [.1 .3 .5 .7 .9], [1 2]); % as there are sometimes to little trials per all conditions, we need to resort to a smaller number.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------ Calculate DME ---------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unfortuantly, the dot coordinates aren't saved in Hannah presenation code
% and therefore we are using on/off timeseries as input to the modelfit.
% e.g. here instead of DME we create a timeseries to be 0 when
% off and 1 for low coherence motion and 2 for high motion coherence
% to create the  determinatistic part and later onwards noise will be fitted 
% to mimick the stocastic part.

for indPP = 1:length(parameters.ppNames)
    fprintf('Calculate Momentary Evidence for %i out of %i.\n', indPP, length(parameters.ppNames))
    
    currOutput = fullfile(parameters.outputFolder, 'MomEvidence data', parameters.ppNames{indPP});
    if ~exist(currOutput, 'dir'), mkdir(currOutput); end
    
    if ~exist(fullfile(currOutput, 'MomEvidence.mat'), 'file')
        acc = 0;
        MomEvidence = [];
        
        for indBlock = 1:parameters.numBlocks(indPP)
            load(fullfile(parameters.inputFolder, parameters.ppNames{indPP},...
                'Trial files', parameters.dataFiles{indPP}{indBlock}), 'par', 'coh');
            
            % get timecourse and resample --> weirdly save with a sample
            % rate of 15...
            currITI   = parameters.experiment{indPP}{indBlock}.ITI.*60;
            stim2frames = parameters.stim.refreshRate/parameters.stim.duration;
            maxLength   = max(currITI) + stim2frames;
            
            for indTrial = 1:length(currITI)
                acc = acc + 1;
                
                MomEvidence(acc,:) = nan(1,maxLength);
                MomEvidence(acc,1:(stim2frames+currITI(indTrial))) = [zeros(1, currITI(indTrial) ,1)...
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
plotExperiment   = 0; % 1 is counts, 2 is percentages. This will be depending on which are used for the model.

% We started of with running the whole model to get an general
% understanding of the data. Stats are reported in 'Behavioural signatures
% of contextual adjustments'. See the .log file that is created to see for
% effects of RT tested with ranova as normal distribution is assumed. 

[rm, ~, tbl] = parameters.doRANOVA(plotExperiment, [1 2]);

% HIT RATE has a binominal distribution so we use a generally linear model
% to test for significance. 
lme_accuracy = fitglme(tbl, 'Accuracy ~ 1 + Context * Evidencestrength + (1 | ppNames)',...
    'Distribution','Binomial',...
    'FitMethod','Laplace');

% FALSE ALARM count has a poission distribution so we use a generally linear model
% to test for significance. 
lme_FAcount= fitglme(tbl, 'FAcount ~ 1 + Context * Evidencestrength + (1 | ppNames)',...
    'Distribution','Poisson',...
    'FitMethod','Laplace');

% Significant interaction between context and evidence strength. Where
% there basically is only an effect of context (easy fixed < mixed < hard fixed)
tblFixed = tbl(double(tbl.Context) == 1,:);
lme_FAFixed= fitglme(tblFixed, 'FAcount ~ 1 + Evidencestrength + (1 | ppNames)',...
    'Distribution','Poisson',...
    'FitMethod','Laplace');

tblMixed = tbl(double(tbl.Context) == 2,:);
lme_FAMixed= fitglme(tblMixed, 'FAcount ~ 1 + Evidencestrength + (1 | ppNames)',...
    'Distribution','Poisson',...
    'FitMethod','Laplace');


%% Replot False alarms
% As the code plots everything depending on the condition combintation and
% as FA are independent on the upcoming motion coherence (presumanly). 
figureFolder = fullfile(parameters.figFolder, 'groupAverage', 'Behaviour');

% False alarms are quite crucial for Hannah data, as we are looking into
% how a stable context is changing the 'bound'. Mixed condition is giving
% us about the same number of False alarms as the hard conditions.

% Illustrating how to easily add conditions --> addition of Context
for indPP = 1:length(parameters.ppNames)
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,1) == 2, 4) = 3;
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,1) == 1 & parameters.behaviour{indPP}.trialMatrix(:,2) == 70, 4) = 2;
    parameters.behaviour{indPP}.trialMatrix(parameters.behaviour{indPP}.trialMatrix(:,1) == 1 & parameters.behaviour{indPP}.trialMatrix(:,2) == 25, 4) = 1;
end

parameters.conditions{4} = [1 2 3];
parameters.figLayOut.legTitle{4} = 'Difficulty Context';

% Uncomment to just plot False alarm rate, but for paper we use the model
% adjusted FA plotting. 
%{
for indPP = 1:length(parameters.ppNames)
    for indCond = 1:3
        posComb   = unique(parameters.behaviour{indPP}.trialMatrix(:,4), 'rows');
        currCond  = all(parameters.behaviour{indPP}.trialMatrix(:, 4) == indCond,2);
        
        % then add a final bin that counts the false alarms:
        tmpFART  = parameters.behaviour{indPP}.indFalseAlarm(currCond,:);
        currFA   = parameters.behaviour{indPP}.FalseAlarm(currCond,:);
        
        qn(indCond,indPP) = sum(sum(currFA(:)))./(sum(parameters.behaviour{indPP}.ITI(currCond))./2); % number of False alarm per 2 sec.
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

xtickangle(45)
xticks(1:indCond)
xticklabels({'Weak', 'Strong', 'Mixed'});
xlim([0 4])
ylim([0 0.2])
yticks(0:0.1:0.2)

set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'FA.png', figureFolder, [3.5 3.5]);

fprintf('Behaviour data has already been run.\n')
%}

% ----------------------  FA and ITI length -------------------------------
% Uncomment to plot as a function of ITI and context, but not used for the paper

[within, between, modelFun, allTrialMatrix] = parameters.getConditions([3 4]);

acc = 0;
for indITI = 2:2:8
    for indCond = 1:3
        
        acc = acc+1;
        forHistogram{acc} = [];  numResponse{acc} = 0;
        for indPP = 1:length(parameters.ppNames)
            
            currCond  = parameters.behaviour{indPP}.trialMatrix(:, 4) == indCond & parameters.behaviour{indPP}.trialMatrix(:, 3) == indITI;
           
            possibilities = (sum(parameters.behaviour{indPP}.ITI(currCond))./2); % per 2 sec.
            % then add a final bin that counts the false alarms:
            currFA   = parameters.behaviour{indPP}.FalseAlarm(currCond,:);
            tmpFART  = parameters.behaviour{indPP}.indFalseAlarm(currCond,:);
            
            tmpFART(currFA == 0) = NaN;
            
            qn(acc,indPP) = sum(currFA(:))./possibilities;
            
            forHistogram{acc}    = [forHistogram{acc}; tmpFART(~isnan(tmpFART))];
            numResponse{acc}     = numResponse{acc} + (sum(parameters.behaviour{indPP}.ITI(parameters.behaviour{indPP}.trialMatrix(:, 4) == indCond))./2);
        end
    end
end
%{
FATable = array2table(qn'.*100);
withinTable = table(repmat(categorical([1:3]'),4,1), repmat([2:2:8]', 3,1),'VariableNames', {'Context', 'ITIDuration'});
rm.FA = fitrm(FATable, 'Var1-Var12 ~ 1', 'WithinDesign', within.Table);

spherTest.FA = mauchly(rm.FA);
results.FA   = ranova(rm.FA, 'WithinModel', modelFun);
multcompare(rm.FA, 'ITIduration', 'by', 'DifficultyContext')
plotData = margmean(rm.FA, {'DifficultyContext',  'ITIduration'});

FAcount = plotData.Mean; % and qn are the average number of trials in each bin
FAci   =  1.96*plotData.StdErr;
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
for indCond = 1:3
    errorbar(2:2:8, FAcount(double(plotData.DifficultyContext) == indCond)', FAci(double(plotData.DifficultyContext) == indCond)',...
        'Color', parameters.figLayOut.colours(indCond,:),'LineWidth',parameters.figLayOut.lineWidth)
end

ylabel('False alarms');
xticks([2:2:8])
xticklabels({'2', '4', '6', '8'});
xlabel('ITI duration');

xlim([1.5 8.5])
ylim([0 18])
yticks(0:2:18)
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'FAITI.png', figureFolder, [5 4]);
fprintf('Behaviour data has already been run.\n')
%}

% ---------------------- FA through the whole ITI -------------------------
% Plotting FA distribution throughout the ITI for FIGURE 3A.

parameters.figLayOut.colours = repmat(ColoursCond12([1 2 4],:),4,1); 
parameters.figLayOut.lineType = {'-', '--', '-.', ':'};

figure('units','normalized','outerposition',[0 0 1 1]); hold on;
timeSeries{1} = fliplr(-0.5:1:1.5).*-1;
timeSeries{2} = fliplr(-0.5:1:3.5).*-1;
timeSeries{3} = fliplr(-0.5:1:5.5).*-1;
timeSeries{4} = fliplr(-0.5:1:7.5).*-1;

% plot hard context
acc = 0;
for indCond = 1:3:length(forHistogram)
    acc = acc+1;
    x = hist(forHistogram{indCond},timeSeries{acc});
    xgrid = linspace(timeSeries{acc}(1),timeSeries{acc}(end),length(x))';
    
    line(timeSeries{4}(1:length(xgrid)-1), (x(1:end-1)/numResponse{indCond}),...
        'Color', parameters.figLayOut.colours(indCond,:), 'LineWidth', 2,...
        'LineStyle', parameters.figLayOut.lineType{acc});
end

% plot easy context
acc = 0;
for indCond = 2:3:length(forHistogram)
    acc = acc+1;
    x = hist(forHistogram{indCond},timeSeries{acc});
    xgrid = linspace(timeSeries{acc}(1),timeSeries{acc}(end),length(x))';
    
    line(timeSeries{4}(1:length(xgrid)-1), (x(1:end-1)/numResponse{indCond}),...
        'Color', parameters.figLayOut.colours(indCond,:), 'LineWidth', 2,...
        'LineStyle', parameters.figLayOut.lineType{acc});
end

acc = 0;
for indCond = 3:3:length(forHistogram)
    acc = acc+1;
    x = hist(forHistogram{indCond},timeSeries{acc});
    xgrid = linspace(timeSeries{acc}(1),timeSeries{acc}(end),length(x))';
    
    line(timeSeries{4}(1:length(xgrid)-1), (x(1:end-1)/numResponse{indCond}),...
        'Color', parameters.figLayOut.colours(indCond,:), 'LineWidth', 2,...
        'LineStyle', parameters.figLayOut.lineType{acc});
end

ylabel(' ');
xlim([-8 0.350])
xticks(-8:2:0.5)
xticklabels(0:2:8);
xlabel(' ');

ylim([0 0.01]); yticks(0:0.005:0.01); 
line([0 0], [0 0.01], 'Color', 'black', 'LineWidth', 1.5)
set(gca,'FontSize', parameters.figLayOut.letterSize);
set(gca,'FontName', parameters.figLayOut.letterType);
plotSave(gca, 'FAITIHistogram.png', figureFolder, [4 6.5154]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------     SET-UP MODELLING  -------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now we move towards the modelling code. This clearly is a big TODO in the
% sense of making the code independent. Mostly to also allow purely
% behavioural modelling. 
parameters.stim.epoch           = parameters.eeg.epoch;
parameters.stim.epochPlot       = parameters.eeg.epochPlot;
parameters.stim.targetEpoch     = parameters.eeg.targetEpoch;
parameters.stim.responseEpoch   = parameters.eeg.responseEpoch;

parameters.stim.baselineCorrect = 0 + [-1 1]*1/parameters.stim.freqSSVEP; % note if you add a second row, these will be used for FA plots. Only nesc. when you are using target-locked topos

parModelling = dataModelling.getInstance(parameters); 
clear parameters

% ------------------   preset everything for modelling -------------------
% set some figure stuff
parModelling.figLayOut.colours     = ColoursCond12;
parModelling.figLayOut.lineWidth   = 1.2;
parModelling.figLayOut.lineType    = {'-' '-' '-' '-'};
parModelling.figLayOut.targetLim   = [-0.2 0:0.5:1];
parModelling.figLayOut.responseLim = [-0.4 0 0.2];
parModelling.figLayOut.saveDim     = [4 6];

% for transferability we have to indix the coherences level as drift
% specificily depend on this.
parModelling.modelBehaviour.ChiOrG = 2;    % 1 is taking the count, 2 is taking the proportions
parModelling.modelBehaviour.cohs = parModelling.conditions{2}; 
parModelling.modelBehaviour.qps  = [.1 .3 .5 .7 .9];   

% get RT quantiles plots and put it on the 
parModelling.plotRTquantiles([4 2]); 

close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------     Behavioural modelling ---------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of testing the hypothesis that this behavioural profile is due to
% strategical context-dependent adjustment of the BOUND. Meaning that the
% initial models have a free bound parameters per condition. The
% within-trial noise e.g. stocastic part in the model is set as a scaling
% parameter set as at s = 0.1 for all condition. This is effectiviely the
% noise in the sensory evidence, e.g. momentary evidence. 
% fminsearch (SIMPLEX) Settings:
parModelling.options = optimset('fminsearch');

% Initial fits are down relatively freely with little iterations and high
% tolarence
numRepeats = 1000;
parModelling.options.MaxIter = 500; % maximum iterations
parModelling.options.TolFun  = 1;	% tolerance - if changes in the objective function are consistently less than this on lots of consecutive evaluations, then it's probably as good as you're going to get, so fminsearch will quit and move on. Don't know why set at 5, presumably trial and error, but this is worth trying to reduce in the refinement step to make sure it doens't stop until it really has a good fit
parModelling.options.TolX    = 1e-2;

StopEarly = @(x,optimvalues,state) optimvalues.fval>4000;
parModelling.options = optimset(parModelling.options, 'Display', 'iter', 'OutputFcn', StopEarly); % @stopEarly give us an oppertunity to prematurely stop if something has more than 100 interation > 1000

parModelling.modelBehaviour.learnBoost = 3;     % which 'context' are we expecting to be boosted.

parModelling.modelBehaviour.noiseSTD   = 0.1;
parModelling.modelBehaviour.simulateMoreX = 1; % it can be good to simulate more trials than the subjects actually did because then at least the simulated data are reliable.
parModelling.modelBehaviour.reflectingBound = 0;

%% 1) There are many ways to set up the drift rate parameters and for the paper
% we runned the following:
%   1.1) One set  Drift parameter for Weak, scaled for Strong as 70/25* in both Fixed and Mixed.
%   1.2) Two free Drift parameters, one for Weak one for Strong, assumed same in Mixed.

% Models accounting for the increase in reaction time in the Mixed
% condition.
% Initially assuming that this is a training effect (as mixed was always last)
%   1.3) Three free Drift parameters; one for Weak and one for Strong in the Fixed, 
%      and a scaled *k version of those in Mixed (assuming training effect).
% Additionally, just checking if there are interaction effects.  
%   1.4) Four free Drift parameters for the 4 different conditions - no
%      constraints at all.


%% 1.1) 3Bound1Drift1Leak1ndt, e.g. 3 bounds, 1 drift, 1 leak and 1 ndt
% First a leaky-accumulation model ONLY allowing Bound-adjustment, with the 
% same non-decision time across conditions, same leak and same drift.
% This is going to be used to refine the models with the flexible drift
% rates as described above. 
ModelName = 'new3Bound1Drift1Leak1ndt';
modelPar = table(char(repmat('bound',3,1), 'drift', 'leak', 'ndt'),...
    [repmat(0.1,3,1);	0;     .001;    .05],...
    [repmat(2,3,1);     .4;    1;       .35],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

Params_3Bound = parModelling.applyModelling(ModelName, numRepeats, modelPar);

% refine this model with more strict SIMPLEX parameters
parModelling.options.MaxIter = 1000; % maximum iterations
parModelling.options.TolFun  = 1e-4;
parModelling.modelBehaviour.simulateMoreX = 5; % model several 'participants'
Params_3Bound = parModelling.applyModelling(ModelName, numRepeats, modelPar, Params_3Bound(1:50,:)); % refine first 50 parameters.

% Uncomment below to see the fit of the simplest model. 
% [pred, output] = parModelling.plotModellingResults(Params_3Bound(1,:), ModelName, modelPar.Names);
% adaptedFAPlot(parModelling, ModelName, pred);

% Additionally, you can plot the parameters. This might help determine if
% some parameters are edging into the one of the limits determined in
% modelPar
%
% parModelling.plotParameters(Params_3Bound(1:50,[5 1:3]), ModelName,...
%     {'Leak', 'Bound weak', 'Bound strong', 'Bound mixed'},...
%     [modelPar.Lower([5 1:3],:) modelPar.Upper([5 1:3],:)]);

%% The crude fit above will do for starting any of these models 
% These variables will now be used to determine if the drift rate are
% better of being fitted as well. 

close all
clc
numRepeats = 200;
parModelling.options.TolFun = 1e-2;
parModelling.modelBehaviour.simulateMoreX = 5;

%%  refine the bound models
% 1.2) Two free Drift parameters, one for Weak one for Strong, assumed same
% in Mixed. These models are NOW uncommented but can be looked into to
% review.
%{

ModelName = '3Bound2Drift1Leak1ndt';
fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)
modelPar = table(char(repmat('bound',3,1), 'leak', 'ndt', repmat('drift',2,1)),...
    [repmat(0.1,3,1);	.001;   .05;   0;	0],...
    [repmat(2,3,1);     1;      .35;   .4;	.4],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

Params_3Bound2Drift = parModelling.applyModelling(ModelName, numRepeats, modelPar, Params_3Bound(1:20,[1:3 5:end]),...
    [Params_3Bound(1:20,4) Params_3Bound(1:20,4)*70/25], 6:7);

pred = parModelling.plotModellingResults(Params_3Bound2Drift(1,:), ModelName, modelPar.Names);
adaptedFAPlot(parModelling, pred, ModelName);
%}

% 1.3) Three free Drift parameters; one for Weak and one for Strong in the Fixed, 
% and a scaled *k version of those in Mixed (assuming training effect).
%{
ModelName = '3Bound2Drift1Leak1ndt1Boost';
fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)
modelPar = table(char(repmat('bound',3,1), 'leak', 'ndt', repmat('drift',2,1), 'boost'), ...
    [repmat(0.1,3,1); 	.001; 	.05;    0;      0;      .8],...
    [repmat(2,3,1);     1;      .35; 	.4; 	.4;     1.5],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

Params_3Bound2Drift1Boost = parModelling.applyModelling(ModelName, numRepeats, modelPar, Params_3Bound2Drift(1:20,:),  [],  8);

pred = parModelling.plotModellingResults(Params_3Bound2Drift1Boost(1,:), ModelName, modelPar.Names);
adaptedFAPlot(parModelling, pred, ModelName);
%}

% 1.4) Four free Drift parameters for the 4 different conditions - no
% constraints at all.
ModelName = '3Bound4Drift1Leak1ndt1';

modelPar = table(char(repmat('bound',3,1), 'leak', 'ndt', repmat('drift',4,1)),...
    [repmat(0.1,3,1);   .001;  .05;   0;      0;      0;      0],...
    [repmat(2,3,1);     1;     .35;   .4;    .4;    .4;    .4],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

[Params_3Bound4Drift, Err] = parModelling.applyModelling(ModelName, numRepeats, modelPar, Params_3Bound(1:20,[1:3 5:end]),...
    [Params_3Bound(1:20,4) Params_3Bound(1:20,4)*70/25 Params_3Bound(1:20,4) Params_3Bound(1:20,4)*70/25], 6:9);

% for FIGURE 1B (DV sim. are not used but plotted and saved)
[pred, outputParameters.Bound] = parModelling.plotModellingResults(Params_3Bound4Drift(1,:), ModelName, modelPar.Names);
adaptedFAPlot(parModelling, pred, ModelName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------     neurally-constrained model ---------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models with urgency signal determined by beta in long ITI (see makeUrgencyFn.m)
% Two ways are implemented first 2) standard leaky-accumulator with the
% addition of the urgency signal as well as a 3) non-leaky accumulator with a
% reflective bound and a criteria adjustment
%
% There's lots more neurally-informed models to try after this - other drift rate setups as above, but also adding back
% in the motor fluctuation parameter sz somehow (e.g. what if we just set it so that motor fluct have the same overall
% s.d. as the within trial sensory evidence noise s, and so we would still only have one free parameter) and it will be
% important to try the leak adjustment model with the neurally-constrained urgency.

clc
numRepeats  = 1000;
StopEarly = @(x,optimvalues,state) optimvalues.fval>5000;
parModelling.options = optimset(parModelling.options, 'Display', 'iter', 'OutputFcn', StopEarly); % @stopEarly give us an oppertunity to prematurely stop if something has more than 100 interation > 1000

load('Urgency.mat') % Use to contrain e.g. motor preperation
parModelling.modelBehaviour.Urgency = urg;

%%  2) Now neurally-constrained model with urgency signal determined by beta in long ITI:
% First, leaky accumulator with different drift rate parameters. 

parModelling.modelBehaviour.reflectingBound = 0;

% 2.1) 3Leak1Drift1ndt1stdNoise, e.g. 1 leak, 1 drift, 1 leak, 1 ndt and a
% noise parameter.
% NOTE here the noise is fitted as we have normalized the bound to be set
% at 1.
% first estimated the parameters
ModelName   = 'NI_3Leak1Drift1ndt1stdNoise';

StopEarly = @(x,optimvalues,state) optimvalues.fval>4000;
parModelling.options = optimset(parModelling.options, 'Display', 'iter', 'OutputFcn', StopEarly); % @stopEarly give us an oppertunity to prematurely stop if something has more than 100 interation > 1000

modelPar = table(char(repmat('Leak',3,1), 'drift', 'ndt', 'trialnoise'),...
    [repmat(.001,3,1);	0;      .05; 	.001],...
    [ones(3,1);         .4; 	.35; 	.4],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

parModelling.options.TolFun     = 1;
parModelling.options.MaxIter    = 500; % maximum iterations
parModelling.modelBehaviour.simulateMoreX = 1;
parModelling.modelBehaviour.noiseSTD  = 1; % set on 1 as it will be multiplied with the trialnoise. 

Params_NI3Leak = applyNIModelling(parModelling, ModelName, numRepeats, modelPar);

% refine fitting
parModelling.options.TolFun  = 1e-4;
parModelling.options.MaxIter = 1000; % maximum iterations
parModelling.modelBehaviour.simulateMoreX = 5;
Params_NI3Leak = applyNIModelling(parModelling, ModelName, numRepeats, modelPar, Params_NI3Leak(1:25,:));

% possibility to plot these initial fit with the scaled drift rate.
%{
noise = parModelling.modelBehaviour.noiseSTD*randn(length(parModelling.modelBehaviour.datsum.TOW), max(parModelling.modelBehaviour.datsum.maxln), parModelling.modelBehaviour.simulateMoreX); % note the 0.1* because assuming s = 0.1. Not done for NI models below
[err, simdat, ~, DV] = NIDecisionModels(parModelling, Params_NI3Leak(1,:),  modelPar.Names, urg, noise,1);
pred = parModelling.plotModellingResults(Params_NI3Leak(1,:), ModelName, modelPar.Names, err, simdat, DV);
adaptedFAPlot(parModelling, pred, ModelName);
%}

% 2.1) NI_3Leak4Drift1ndt1stdNoise, e.g. 3 leak, 4 drift, 1 leak, 1 ndt and a
% noise parameter.
numRepeats = 200;

ModelName = 'NI_3Leak4Drift1ndt1stdNoise';
modelPar = table(char(repmat('Leak',3,1), 'ndt', 'trialnoise', repmat('drift',4,1)),...
    [repmat(.001,3,1);	.05;   .001;	0;      0;      0;      0],...
    [ones(3,1);         .35;   2;       .4;     .4; 	.4; 	.4],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

Params_NI3Leak4drift = applyNIModelling(parModelling, ModelName, numRepeats, modelPar, Params_NI3Leak(:, [1:3 5:end]),...
    [Params_NI3Leak(:,4) Params_NI3Leak(:,4)*70/25 Params_NI3Leak(:,4) Params_NI3Leak(:,4)*70/25], 6:9);

 
noise = parModelling.modelBehaviour.noiseSTD*randn(length(parModelling.modelBehaviour.datsum.TOW), max(parModelling.modelBehaviour.datsum.maxln), parModelling.modelBehaviour.simulateMoreX); % note the 0.1* because assuming s = 0.1. Not done for NI models below
[err, simdat, ~, DV] = NIDecisionModels(parModelling, Params_NI3Leak4drift(1,:),  modelPar.Names, 1);

[pred, outputParameters.Leak] = parModelling.plotModellingResults(Params_NI3Leak4drift(1,:), ModelName, modelPar.Names, err, simdat, DV);
adaptedFAPlot(parModelling, pred, ModelName);

%%  3) Now neurally-constrained model with urgency signal determined by beta in long ITI:
% finally, non-leaky accumulator with criteria adjustment and a refelctive
% bound (to account for 'leak').

clc
close all
numRepeats   = 1000;
parModelling.modelBehaviour.reflectingBound = 1;

% 3.1) model with reflective bound and fitting 3 criteria
ModelName = 'NI_3Criteria1Drift1ndt1stdNoise';

modelPar = table(char(repmat('criteria',3,1), 'drift', 'ndt', 'trialnoise'),...
    [zeros(3,1);        0;      .05;	.001],...
    [repmat(0.2,3,1);   .4; 	.35;    2],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

parModelling.options.TolFun     = 1;
parModelling.options.MaxIter    = 500; % maximum iterations
parModelling.modelBehaviour.noiseSTD     = 1;
parModelling.modelBehaviour.simulateMoreX = 1;
Params_NI3Criteria = applyNIModelling(parModelling, ModelName, numRepeats, modelPar);

% refine NI 3 criteria model
parModelling.options.TolFun = 1e-4;
parModelling.options.MaxIter = 1000; % maximum iterations
parModelling.modelBehaviour.simulateMoreX = 5;
Params_NI3Criteria = applyNIModelling(parModelling, ModelName, numRepeats, modelPar, Params_NI3Criteria(1:25,:));

%{
noise = parModelling.modelBehaviour.noiseSTD*randn(length(parModelling.modelBehaviour.datsum.TOW), max(parModelling.modelBehaviour.datsum.maxln), parModelling.modelBehaviour.simulateMoreX); % note the 0.1* because assuming s = 0.1. Not done for NI models below
[err, simdat, ~, DV] = NIDecisionModels(parModelling, Params_NI3Criteria(1,:),  modelPar.Names, urg, noise, 1, applyReflectiveBound);

pred = parModelling.plotModellingResults(Params_NI3Criteria(1,:), ModelName, modelPar.Names, err, simdat, DV);
adaptedFAPlot(parModelling, pred, ModelName);
%}

numRepeats = 200;

% 3.2) model with reflecting bound with 3 criteria and 4 drifts
ModelName = 'NI_3Criteria4Drift1ndt1stdNoise';
fprintf('NOW FITTING THE FOLLOWING REFINEMENT: %s \n', ModelName)

modelPar = table(char(repmat('criteria',3,1), 'ndt', 'trialnoise', repmat('drift',4,1)),...
    [zeros(3,1);        .05; 	.001;   0;      0;      0;      0],...
    [repmat(.2,3,1);    .35;    2;      .4;     .4;     .4;     .4],...
    'VariableNames', {'Names', 'Lower', 'Upper'});

Params_NI3Criteria4drift = applyNIModelling(parModelling, ModelName, numRepeats, modelPar,  Params_NI3Criteria(:, [1:3 5:end]),...
   [Params_NI3Criteria(:,4) Params_NI3Criteria(:,4)*70/25 Params_NI3Criteria(:,4) Params_NI3Criteria(:,4)*70/25], 6:9);


noise = parModelling.modelBehaviour.noiseSTD*randn(length(parModelling.modelBehaviour.datsum.TOW), max(parModelling.modelBehaviour.datsum.maxln), parModelling.modelBehaviour.simulateMoreX); % note the 0.1* because assuming s = 0.1. Not done for NI models below
[err, simdat, ~, DV] = NIDecisionModels(parModelling, Params_NI3Criteria4drift(1,:), modelPar.Names, noise, 1);
[pred, outputParameters.Criteria]  = parModelling.plotModellingResults(Params_NI3Criteria4drift(1,:), ModelName, modelPar.Names, err, simdat, DV);
adaptedFAPlot(parModelling, pred, ModelName);
