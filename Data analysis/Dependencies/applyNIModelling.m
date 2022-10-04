function [bestpmFits, EE] = applyNIModelling(obj, modelName, urgency, reflectingBound, initialGuess, modelPar, numRepeats, jitterPar, freeParameters)
saveFolder = fullfile(obj.outputFolder, 'Modelling');
if ~exist(saveFolder, 'dir'); mkdir(saveFolder); end
if ~exist('reflectingBound','var'); reflectingBound = 0; end

if isempty(initialGuess) || exist('freeParameters', 'var')
    % for initial modelling, we start with generate a random matrix with
    % starting parameter vectors the specific range given in
    % modelPar.
    if exist('freeParameters', 'var')
        index  = 1:size(modelPar,1);
        newParameters = nan(numRepeats, size(modelPar,1));
        newParameters(:,~ismember( index, freeParameters))	= initialGuess(randi(size(initialGuess,1),numRepeats,1),:) + (randn(numRepeats, size(initialGuess,2)).*mean(initialGuess).*0.1);
        
        if ~isempty(jitterPar)
            newParameters(:,ismember(index, freeParameters))   = jitterPar(randi(size(jitterPar,1),numRepeats,1),:) +  (randn(numRepeats, size(jitterPar,2)).*mean(jitterPar).*0.15);
        else
            newParameters(:,ismember(index, freeParameters))  = generateRandParameters(obj, modelPar, numRepeats);
        end
        initialGuess = newParameters;
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
        
        [pm, Err] = fminsearchbnd(@(pm)NIDecisionModels(obj, pm, modelNames, urgency, noise, 0, reflectingBound), initialGuess(indPM,:), lowerBound, upperBound,  obj.options);
        
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