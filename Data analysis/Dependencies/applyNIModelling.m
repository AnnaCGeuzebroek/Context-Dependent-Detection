function [bestpmFits, EE] = applyNIModelling(obj, modelName, numRepeats,  modelPar, initialGuess, jitter, freeParams)
%% applyNIModelling
% as if for now, this is a seperated code, but eventually this will be
% intergrated in applymodelling in the dataModelling object.
% Needs a urgency signal and reflectingbound in obj.modelbehaviour.
% set the intial parameters, for now it is just an extensive
% search within the bounds given in modelPar. After which it
% will run through a parfor to optimize the intial search
% parameters.
% Use as:
%   modelName   =   Just in order to save stuff (such as
%                   parameters fits).
%   numRepeats  =   Number use to sample the parameters space.
%                   otherwise just use the intialGuess to
%                   refine.
%   intialGuess =   when left empty 'generateRandParameters' is
%                   used to search whole parameters space given
%                   by the limits in modelPar.
%                   Otherwise, it will take the intially
%                   crudely fitted parameters to start of with.
% There are also a couple of ways to use intial guesses to
% jitter around or created completely new parameters.
%	jitter      =	Sometimes a new estimated can be informed
%                   by previously fitted parameters. This
%                   can use intialGuess.
%	freeParams  =   Create completely new parameters with
%                   modelPar which will be added to

saveFolder = fullfile(obj.outputFolder, 'Modelling');
if ~exist(saveFolder, 'dir'); mkdir(saveFolder); end

if ~exist('initialGuess', 'var');   initialGuess   = []; end
if ~exist('freeParams', 'var');     freeParams = []; end
if ~exist('jitter', 'var');         jitter = []; end

if isempty(initialGuess) || ~isempty(freeParams)
    % for initial modelling, we start with generate a random matrix with
    % starting parameter vectors the specific range given in
    % modelPar.
    if ~isempty(freeParams)
        
        index  = 1:size(modelPar,1);
        
        % preset newParams, with 'little' var around
        % intialGuess that are 'fixed'.
        newParams = nan(numRepeats, size(modelPar,1));
        newParams(:, ~ismember(index, freeParams)) = initialGuess(randi(size(initialGuess,1),numRepeats,1),:) + (randn(numRepeats, size(initialGuess,2)).*mean(initialGuess).*0.1);
        
        for indRow = 1:size(freeParams,1)
            if isempty(jitter) ||(indRow == 1 && size(freeParams,1) > 1)
                newParams(:,ismember(index, freeParams(indRow,:)))  = generateRandParameters(obj, modelPar(freeParams(indRow,:),:), numRepeats);
            else
                newParams(:,ismember(index, freeParams(indRow,:)))  = jitter(randi(size(jitter,1),numRepeats,1),:) +  (randn(numRepeats, size(jitter,2)).*mean(jitter).*0.15);
            end
        end
        initialGuess = newParams;
    else
        initialGuess = generateRandParameters(obj, modelPar, numRepeats);
    end
    
    saveFile = ['paramfits_' modelName '_1'];
else
    % when model is already run, we are going to reuse the
    % refitted initialGuesss and refine them. We here save the
    % files as refine.
    saveFile = ['paramfits_' modelName '_2'];
end

% include check if it already exist.
if ~exist(fullfile(saveFolder, [saveFile '.mat']), 'file')
    EE  = nan(size(initialGuess,1),1);
    pmfits = nan(size(initialGuess));
    
    numInitialGuess = length(initialGuess);
    modelNames = modelPar.Names;
    lowerBound = modelPar.Lower';
    upperBound = modelPar.Upper';
    
    noise = obj.modelBehaviour.noiseSTD*randn(length(obj.modelBehaviour.datsum.TOW), max(obj.modelBehaviour.datsum.maxln), obj.modelBehaviour.simulateMoreX); % note the 0.1* because assuming s = 0.1. Not done for NI models below
    parfor indPM = 1:numInitialGuess
        fprintf('start pt %i:', indPM)
        
        [pm, Err] = fminsearchbnd(@(pm)NIDecisionModels(obj, pm, modelNames, noise, 0), initialGuess(indPM,:), lowerBound, upperBound,  obj.options);
        
        % save them:
        EE(indPM) = Err;
        pmfits(indPM,:) = pm';
        
        % only plot if error semi-decent
        if Err <= 100
            fprintf('start pt %i: E = %0.2f, pm = %s\n', indPM, EE(indPM), num2str(pmfits(indPM,:)))
        end
    end
    
    save(fullfile(saveFolder, saveFile), 'EE', 'pmfits', 'initialGuess')
else
    load(fullfile(saveFolder, saveFile), 'EE', 'pmfits')
end

% Get best parameter fits as output
[~, BestEE] = sort(EE);
EE = EE(BestEE);
bestpmFits  = pmfits(BestEE, :);
end