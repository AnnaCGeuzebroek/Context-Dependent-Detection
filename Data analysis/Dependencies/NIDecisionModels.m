%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------  supporting functions     ----------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------  Leaky accumulation  ---------------------------
function [err, simdat, pred, DecValue] = NIDecisionModels(obj, pm, parNames, noise, getsimDV)
%% diffusionModel
% Implements the sequential sampling modelling based on
% accumulation where several parameters can be set to be free.
% Function can be used to optimize within SIMPLEX function (see
% obj.applyModelling), as well as just on each own, allowing
% you to extract a prediction of a certain parameters
% combination and the simulated decision value.  Model based on
% Ossmy, et al 2013 (https://pubmed.ncbi.nlm.nih.gov/23684972/)
% Parameters can be set free:
%   bound   =   Threshold determining the amount of accumulate
%               before commiting to a decision.
%   drift   =   average amount of evidence accumulated per unit
%               time, this can be an index of task difficulty
%               or on participants ability. When having several
%               difficutly levels, 1) this can either be set as
%               one parameters which will be scaled accordingly
%               2) one for each difficulty level or 3)
%               indiviudally fitted drift per difficulty level
%               per condition.
%   leak    =   How much previous accumulated evidence (t-1) is
%               weighted to current incoming evidence (t). Set
%               between 0 = perfect accumulation and 1 = full
%               leak e.g. decide based on current evidence
%               only.
%   ndt     =   non-decision time accounting for motor response
%               delays
%   boost   =   possibilty boosting the drift rate, if we
%               expect learning effects --> index to which
%               condition in the first column of the
%               trialmatrix.
%   criteria =  reference on momentary evidence, determine what is see as
%               signal or noise (similar to signal detection theory)
%
% inputParameters include
%   1) pm       =   which are either estimated in fminsearchbnd using
%                   going through the fit parameters space or a
%                   single set of parameters to get the output.
%   2) parNames =   Names corresponding to the pm, should be
%                   ones as described above and for each pm
%                   one!
%   3) noise    =   set noise levels, e.g. the stocastic part.
%   4) getsimDV =   0 - don't, 1 - get simulated DV.
%
if ~exist('getsimDV','var'); getsimDV = 0; end


% extract parameters from modelBehaviour.
urgency         = obj.modelBehaviour.Urgency;
reflectingBound = obj.modelBehaviour.reflectingBound;

TOW    = obj.modelBehaviour.datsum.TOW;
bt     = obj.modelBehaviour.datsum.bt;
trialMatrix = obj.modelBehaviour.datsum.trialMatrix;

% get timing parameters
dt = 1/obj.stim.refreshRate;      % sample period
maxRT = obj.stim.RTdeadLine(end); % the upper limit for RTs to be classed as hits
minRT = obj.stim.RTCutOff;        % minium allowable RT to call a response a 'hit', in msec. There was only one RT in the S blocks across all subjects that was shorter than this, and next shortest was nearly 100 ms later


% don't start accumulating until this number of samples after response. This is simply a reasonable guess.
% HOWEVER, it might be that targets appear together and this needs to be adjusted!!
% TODO make this independent on refresh rate and rather on
% timing!!
postRespPause = 1*obj.stim.refreshRate;

% Without loss of generality, and to facilitate DV simulation and comparison with NI models later, will assign 67 ms of
% the nondecision time to the motor end, and the free part will be for 'PRE' decisional NDT:

MT = 4*dt;
continAccDur = obj.stim.refreshRate/10;    % for how many samples (100 ms) should the CPP keep accumulating after reaching commitment? From R-locked CPPs it looks like around 100 ms
waitAfterTarget = obj.stim.refreshRate/2;  % start the urgency signal half a second into the ITI

% preset and extract parameters.
bound = []; drift = []; leak = []; ndt = []; boost = []; criteria = []; trialnoise = [];
for indParName = 1:size(parNames,1)
    eval(sprintf('%s(end+1) = pm(:, %i);', lower(parNames(indParName,:)), indParName))
end

% Some parameters are not going to be used, they do need to be
% defined, therefore to save time in the for loop, we set them
% here..
if isempty(bound),      bound = 1; end
if isempty(leak),       leak = 0; end
if isempty(criteria),   criteria = 0; end
if isempty(boost),      boost = 1; end
if isempty(trialnoise), trialnoise = 1; end

% adjust drift accordingly.
if size(drift,2) == 1
    % if only one drift rate is fitted we are going to assume that the
    % drift rate scale depending on the coherence levels.
    
    % however, this is how it is implemented now, making it more
    % transfable (for example with more coherence levels):
    drift = drift*(obj.modelBehaviour.cohs./min(obj.modelBehaviour.cohs)); % scale up the drift rate parameters for other conditions
elseif size(drift,2) > length(obj.modelBehaviour.cohs)
    drift = reshape(drift, length(obj.modelBehaviour.cohs), []);
elseif size(boost,2) > 1
    drift = repmat(drift, 2,1).* boost';
end


% preallocated parameters.
simdat  = [];

acc = 1;
hit = 1; miss = 2; false_alarm = 3;     % codes for response types


for indNoise = 1:size(noise,3) % for robustness, it can help to simulate more trials than there are in the real data, This is in line with 'several' participants
    
    % loop through block trials TODO now only continuous, might
    % need to add a trial-loop instead of block for discrete
    % experiments.
    for indBlock = 1:length(TOW)
        
        % get the appropiated fit parametrs per conditions.
        % When there is only 1 it will be fixed accross
        % conditions. If the user set several it will find the
        % appropiated current pm for that condition. Currently,
        % this is just per block (e.g. bt(indBlock, 1)), TODO
        % switch to the condtion using getCondition ect.
        
        if length(bound) > 1;      currBound = bound(bt(indBlock,1));         else; currBound = bound; end
        if length(leak) > 1;       currLeak = leak(bt(indBlock,1));           else; currLeak = leak; end
        if length(criteria) > 1;   currCriteria = criteria(bt(indBlock,1));	  else; currCriteria = criteria; end
        if length(trialnoise) > 1; currTrialnoise = trialnoise(bt(indBlock)); else; currTrialnoise = trialnoise; end
        
        if obj.modelBehaviour.learnBoost == bt(indBlock,1)
            currBoost = boost;
        else
            currBoost = 1;
        end
        
        if size(drift,1) > 1
            if obj.modelBehaviour.learnBoost == bt(indBlock,1)
                currDrift = drift(2,:);
            else
                currDrift = drift(1,:);
            end
        else
            currDrift = drift.*currBoost;
        end
        
        % make sensory evidence waveform - the target-on pulse train plus
        % Gaussian noise   changed this to make it more transferable,
        % e.g. get the timeline which should be indexing the
        % difficulty levels to get the right drift rate
        
        timeLine = TOW{indBlock}';
        timeLine(timeLine ~= 0) = currDrift(timeLine(timeLine ~= 0));
        
        if obj.DetectOrDisc
            direction = TOW{indBlock}(:,1)';
            timeLine(direction < 1) = timeLine(direction < 1)*-1;
        end
        
        sensEv = timeLine + currTrialnoise.*noise(indBlock,1:length(TOW{indBlock}), indNoise);
        
        
        % initalization of parameters.
        DV(1:round(ndt/dt)) = 0;    % initialize DV to a single 0 for the first time points of the block, up until the pre-decision nondecision time
        lastresp = -postRespPause;  % we only accumulate if it's a certain time since last response. Initialise to this value so accumulation begins right away at start of block
        
        RespT   = []; % keep track of all responses times
        targT   = []; % keep track of all targets onset times
        if obj.DetectOrDisc
            Resp = [];  % keep track of all actual response
        end
        
        DVendval = nan; % **************
        
        % now simulate the block by looping through all sample points:
        for indTime = round(ndt/dt)+1:length(sensEv)
            
            if  indTime <= lastresp+continAccDur || indTime > lastresp + postRespPause %  indTime <= lastresp+continAccDur || don't accumulate unless it has been a sufficient time since last response
                % this gives us the 'reflecting bound'
                if reflectingBound
                    DV(indTime) = max((1-currLeak)*DV(indTime-1) + sensEv(indTime-round(ndt/dt)) - currCriteria, 0);% There is the main model equation! (see e.g. Ossmy et al 2013)
                else
                    DV(indTime) = (1-currLeak)*DV(indTime-1) + sensEv(indTime-round(ndt/dt)) - currCriteria; % There is the main model equation! (see e.g. Ossmy et al 2013)
                end
            else
                % As we are using the empiraical CPP as an
                % indeicator of the evidence accumulation, we
                % want to see how the decision variable, e.g.
                % DV, encodes and compares with this. In the
                % future we can also use DV and CPP to contrain
                % our models. However, this requires us to not
                % have a 'absolute off' from the DV.
                
                % Continue accumulation for
                if isnan(DVendval) % I'm using this as a way to know when it is time to linearly ramp the DV down to zero
                    DVendval = DV(indTime-1);
                end
                
                DV(indTime) = max(0, DVendval*(1-(indTime-(lastresp+continAccDur))/25)); % linear decrease to a floor of zero
                if DV(indTime) == 0, DVendval = nan; end % when the CPP has reached back down to zero, turn off the linear ramp-down
            end
            
            % detect target transitions
            if TOW{indBlock}(indTime) > TOW{indBlock}(indTime-1)
                targT = [targT indTime*dt];
            end
            
            % detect responses
            if indTime > waitAfterTarget
                if indTime-waitAfterTarget <= size(urgency,1)
                    urgencyDV(indTime) = DV(indTime) + urgency(indTime-waitAfterTarget, bt(indBlock)); % total decision signal including urgency
                else  % when there's a miss, hard to know what to do with urgency - most neutral thing might be to keep it steady where it is
                    urgencyDV(indTime) = DV(indTime) + urgency(end, bt(indBlock));
                end
            else
                urgencyDV(indTime) = DV(indTime);
            end
            
            % Detect responses: This is when the DV crosses the bound.
            if indTime > lastresp+postRespPause
                if obj.DetectOrDisc
                    if abs(urgencyDV(indTime)) > currBound     % TODO for disc task < -currBound beside RespT check if correct or mistake.
                        RespT    = [RespT indTime*dt + MT];    % log the response after adding the non-decision time
                        lastresp = indTime;                 % and now this is the last response that happened, at sample n
                        
                        if sign(DV(indTime)) == sign(TOW{indBlock}(indTime))
                            Resp = [Resp 1];
                        else
                            Resp = [Resp 0];
                        end
                    end
                else
                    if urgencyDV(indTime) > currBound % TODO for disc task < -currBound beside RespT check if correct or mistake.
                        RespT = [RespT indTime*dt+MT]; % log the response after adding the non-decision time
                        lastresp = indTime; % and now this is the last response that happened, at sample n
                        
                        if ~isempty(targT)
                            waitAfterTarget = round(targT(end)/dt)+maxRT*obj.stim.refreshRate; % start counting the urg again at 0.5 sec after target offset
                        end
                    end
                end
            end
            
            % now make the output matrix that has all the RTs w.r.t. target onset
            % and false alarms. This has to be equivalent to how the real data were
            % analysed! Deals in seconds
            for indTarget = 1:length(targT)
                nextrespind = find(RespT > targT(indTarget)+minRT & RespT < targT(indTarget)+maxRT, 1); % find the index of the next response which is within the allowable @hit@ window
                if ~isempty(nextrespind) % TODO give additional Hit -- correct or wrong.
                    simdat  = [simdat; trialMatrix{indBlock}(indTarget,:) hit RespT(nextrespind)-targT(indTarget)]; % 1 = hit
                else % if there WAS no next response, set the response parameters for this trial as 'not a number'
                    simdat  = [simdat; trialMatrix{indBlock}(indTarget,:) miss nan]; % 2 = miss
                end
                
                timeFrame = [floor((obj.stim.targetEpoch(1))/dt):ceil((obj.stim.targetEpoch(end))/dt)];
                if getsimDV == 1
                    if (targT(indTarget)/dt) + ceil((obj.stim.targetEpoch(end))/dt) < length(DV)
                        DecValue.Target(:,acc) = DV(round(targT(indTarget)/dt) + timeFrame);
                        if ~isempty(nextrespind)
                            DecValue.Response(1:length((floor((obj.stim.responseEpoch(1))/dt):ceil((obj.stim.responseEpoch(end))/dt))),acc) = DV(round(RespT(nextrespind)/dt) + [floor((obj.stim.responseEpoch(1))/dt):ceil((obj.stim.responseEpoch(end))/dt)]);
                        else
                            DecValue.Response(1:length((floor((obj.stim.responseEpoch(1))/dt):ceil((obj.stim.responseEpoch(end))/dt))),acc) = nan;
                        end
                    else
                        DecValue.Target(:,acc) = nan;
                        if ~isempty(nextrespind)
                            DecValue.Response(:,acc) = nan;
                        end
                    end
                    acc = acc + 1;
                end
            end
            % TODO check what to do with FA in disc task
            ITIstartT = [0 targT+1];
            for indITI = 1:length(ITIstartT)
                nexttargind = find(targT > ITIstartT(indITI),1); % index of next target
                % from this establish the end of the window starting from the current ITI start where we will check for FAs
                if ~isempty(nexttargind)
                    endtime = targT(nexttargind) + minRT;
                else
                    endtime = length(DV)*dt-.125; % like in real data, if we didn't find a next target then this must be the end of the block, so check up as far as we ould possibly extract an erpr
                end
                % now find responses in this ITI window:
                nextrespind = find(RespT > ITIstartT(indITI)+maxRT-(obj.stim.duration) & RespT < endtime); % find indices of responses the ITI window, ruling out any at the very start that are within the hit window from the previous target. Target duration is 1 sec, so fs in sample points
                
                for m = 1:length(nextrespind)
                    simdat  = [simdat; trialMatrix{indBlock}(indITI,:) false_alarm RespT(nextrespind(m))-ITIstartT(indITI)];
                    
                    DecValue.Target(:,acc) = nan;
                    % TODO here we can start seeing psychophysical
                    % kernel if we have the actual traces
                    DecValue.Response(:,acc) = nan;
                    
                    acc = acc + 1;
                end
            end
            
        end
    end
    % It's possible that with certain parameters the DV NEVER crosses the bound, so need a quick fix so no error happens:
    % this happened once by total fluke - there was no 'E' below, so there must have been a false alarm before the first target (cond=nan) and nothing else!
    
    if isempty(simdat)
        tmpCond = unique(reshape([trialMatrix{:}]', size(trialMatrix{1},2), [])', 'rows');
        simdat = [tmpCond repmat(3, size(tmpCond,1),1) repmat(10, size(tmpCond,1),1)];
    end
    
    if obj.modelBehaviour.ChiOrG == 1
        [err,pred] = Chisquared(obj, simdat);
    elseif obj.modelBehaviour.ChiOrG == 2
        [err,pred] = Gsquared(obj,simdat);
    end
end