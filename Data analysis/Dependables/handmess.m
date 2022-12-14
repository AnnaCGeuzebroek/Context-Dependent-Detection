function [Quality, figInfo, tbl] = handmess(obj, plotName, Elec, trangeTopo, TargetOrResponse, AverageOrSlope, MorletOrsFFT, freqRange)
%% [Quality, figInfo, tbl] = plotERP(obj, plotName, Elec, plotThis, plotComb, grouping, trangeTopo, TargetOrResponse, MorletOrsFFT, AverageOrSlope)
% Similar to topoplotting in plotERP. It was not clearly instructed with what hand 
% participants had to respond and fortuantly after checking people seem to have used
% there dominated hand (two left-handed). To determine this, we used Beta lateralization 
% (e.g. 16-30 Hz) as well as lateralization of the preresponse slope of the LRP.
%
% Input parameters:
%   MethodUsed          =  gives the possibility to just use
%                          1) the pre-set electrodes
%                          2) get the 3 best electrodes with highest SNR
%                          3) Lucalize a cluster of electrodes
%   plotName            =  given name to the plots (e.g. CPP)
%   Elec                =  selected electrodes (either number
%                          or names), when empty the electrodes
%                          will be asked to be selected during
%                          the topoplot phase.
%   trangeTopo          =  Range in ms.
%   TargetOrResponse    =  1 - Target, 2 - Response
%   MorletOrsFFT        =  0 - just ERP, 1 - Morlet, 2- sFFT
%   AverageOrSlope      =  1 - Average, 2 - Slope

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- PRE-SET PARAMETERS  -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for parameters.
if ~exist('TargetOrResponse', 'var'); TargetOrResponse = 2; end  % standard on response-locked
if ~exist('AverageOrSlope', 'var');   AverageOrSlope = 1; end    % standard on average
if ~exist('MorletOrsFFT', 'var');     MorletOrsFFT = 0; end      % standard on average

trangeBaseline = obj.eeg.epochPlot > obj.eeg.baseline(1) & obj.eeg.epochPlot < obj.eeg.baseline(2);

% preset folders
currOutput = fullfile(obj.outputFolder, 'EEG data', 'groupAverage', plotName);
if ~exist(currOutput, 'dir'); mkdir(currOutput); end

currTopo = fullfile(currOutput, [plotName 'Topo_HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD) '.mat']);

% preset figure save folder
indivFolder = fullfile(obj.figFolder,  'individualPlots', ['HPF' num2str(obj.eeg.HPFcutoff) '_CSD' num2str(obj.eeg.applyCSD) '/HandMess/']);
if ~exist(indivFolder, 'dir'); mkdir(indivFolder); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- CREATE AVERAGE TOPOPLOT  --------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to solve the 'hand mess', we first create an individual topoplot.
% This allows us to visualize and properly chooise the cluster
% of electrodes the check.

if TargetOrResponse == 1
    trangeTopo = obj.eeg.targetEpoch > trangeTopo(1) & obj.eeg.targetEpoch < trangeTopo(2);
elseif TargetOrResponse == 2
    trangeTopo = obj.eeg.responseEpoch > trangeTopo(1) & obj.eeg.responseEpoch < trangeTopo(2);
end

fprintf('Extract data for avarage topoplot to determine best %s channels.\n', plotName)

% 1) extract all topoplot per participant
ERPTopo = nan(1,obj.eeg.NumberOfChannels, max(obj.numBlocks)*obj.numTrials , length(obj.ppNames));

for indPP = [11 14]; %1:length(obj.ppNames)
    clear ERP tmp*
    fprintf(['Now processing participant ' obj.ppNames{indPP} ' to get topoplot\n'])
    currInput  = fullfile(obj.outputFolder,'EEG data', obj.ppNames{indPP});
    
    % load the EEG data.
    if obj.eeg.applyCSD
        load(fullfile(currInput,[obj.ppNames{indPP} '_epochedEEG_HPF' num2str(obj.eeg.HPFcutoff) obj.eeg.timing '.mat']), 'csdERP');
        ERP = csdERP;
    else
        load(fullfile(currInput,[obj.ppNames{indPP} '_epochedEEG_HPF' num2str(obj.eeg.HPFcutoff) obj.eeg.timing '.mat']), 'ERP');
    end
    
    % Baseline-correct the data for the target and
    % response. Again not for the ERPwhole
    ERP = ERP - repmat(nanmean(ERP(:, trangeBaseline, :),2), 1, size(ERP,2), 1);
     
    if MorletOrsFFT ~= 0
        clear down*
        parfor indChan = 1:obj.eeg.NumberOfChannels-1
            downERD(indChan,:,:,:) = obj.shortfft(ERP(indChan,:,:), freqRange);
        end
        [downERD(obj.eeg.NumberOfChannels,:,:,:), timeSeries] = obj.shortfft(ERP(obj.eeg.NumberOfChannels,:,:), freqRange);
        downERD = squeeze(nanmean(downERD,2));
        
        ERP = nan(size(downERD,1), size(obj.eeg.epochPlot,2), size(downERD,3));
        
        for indChan = 1:size(downERD,1)
            for indEpoch = 1:size(downERD,3)
                ERP(indChan,timeSeries(1):timeSeries(end)+ diff(timeSeries(1:2))-1, indEpoch) = interp(downERD(indChan,:, indEpoch), diff(timeSeries(1:2)));
            end
        end
        ERP = ERP - repmat(nanmean(ERP(:, trangeBaseline, :),2), 1, size(ERP,2), 1);
    end
    
    % calculated ERP topography
    if TargetOrResponse == 1
        [~, tmpPlot] = sortERPs(obj, ERP, indPP, 1, 1);
    elseif TargetOrResponse == 2
        [~, ~, tmpPlot] = sortERPs(obj, ERP, indPP, 2, 1);
    end
    
    if AverageOrSlope == 1
        ERPTopo(1,:,:,indPP) = squeeze(nanmean(tmpPlot(:,trangeTopo,:,1),2));
    elseif AverageOrSlope == 2
        ERPTopo(1,:,:,indPP) = squeeze(nanmean(diff(tmpPlot(:,trangeTopo,:,1),[],2)./diff(obj.eeg.targetEpoch(trangeTopo),[],2),2));
    end
    
    orderBlocks = [];
    for indBlock = 1:length(obj.order{indPP})
       orderBlocks = [orderBlocks; repmat(obj.order{indPP}(indBlock), length(obj.experiment{1}{1}.TargOnT),1)];
    end

    figure('Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    for indBlock = 1:length(obj.order{indPP})
        subplot(4, ceil(length(obj.order{indPP})/4), indBlock)
        topoplot(  nanmean(ERPTopo(1,:,orderBlocks == indBlock, indPP),3),obj.eeg.chanlocs); % , 'electrodes', 'labels'        
    end
    plotSave(gca,  [obj.ppNames{indPP} '_' plotName '.png'], indivFolder, [25 20]);
end
