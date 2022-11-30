%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------  supporting functions     ----------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------  Leaky accumulation  ---------------------------
function [err, simdat, pred, DecValue] = NIDecisionModels(obj, pm, parNames, noise, getsimDV)
% Simulates beh matrix for continuous leaky accumulator with certain given parameters
% It has 3 bound parameters (one for W one for S one for M) but single leak parameter and single non-decision time (ndt)
% assumed to be before decision process (think of it as a delay from changes on the screen to being counted in the accumulator)
% and 1 free drift rate parameter for lo coh scaled up for high coh
% datsum = obj.modelBehaviour.datsum;
% Simulates beh matrix for continuous leaky accumulator with certain given parameters

% It has 3 bound parameters (one for W one for S one for M) but single leak parameter and single non-decision time (ndt)
% assumed to be before decision process (think of it as a delay from changes on the screen to being counted in the accumulator)
% and 1 free drift rate parameter for lo coh scaled up for high coh

if ~exist('getCPP','var'); getCPP = 0; end

% pm are parameters
TOW    = obj.modelBehaviour.datsum.TOW;
bt     = obj.modelBehaviour.datsum.bt;
trialMatrix = obj.modelBehaviour.datsum.trialMatrix;
reflectingBound = obj.modelBehaviour.reflectingBound;
urgency         = obj.modelBehaviour.Urgency;

% get timing parameters
dt = obj.stim.duration/obj.stim.refreshRate;      % sample period
maxRT = obj.stim.RTdeadLine(end); % the upper limit for RTs to be classed as hits
minRT = obj.stim.RTCutOff;        % minium allowable RT to call a response a 'hit', in msec. There was only one RT in the S blocks across all subjects that was shorter than this, and next shortest was nearly 100 ms later


% don't start accumulating until this number of samples after response. This is simply a reasonable guess.
% HOWEVER, it might be that targets appear together and this needs to be adjusted!!
% TODO make this independent on refresh rate and rather on
% timing!!
postRespPause = 1*obj.stim.refreshRate;

waitAfterTarget = obj.stim.refreshRate/2;  % start the urgency signal half a second into the ITI

% Without loss of generality, and to facilitate DV simulation and comparison with NI models later, will assign 67 ms of
% the nondecision time to the motor end, and the free part will be for 'PRE' decisional NDT:

MT = 4*dt;
continAccDur = obj.stim.refreshRate/10; % for how many samples (100 ms) should the CPP keep accumulating after reaching commitment? From R-locked CPPs it looks like around 100 ms

targetPlot = floor((obj.stim.targetEpoch(1))/dt):ceil((obj.stim.targetEpoch(end))/dt);
responsePlot = floor((obj.stim.responseEpoch(1))/dt):ceil((obj.stim.responseEpoch(end))/dt);


% preset and extract parameters.
bound = []; drift = []; leak = []; ndt = []; boost = []; criteria = []; trialnoise = [];
for indParName = 1:size(parNames,1)
    eval(sprintf('%s(end+1) = pm(:, %i);', lower(parNames(indParName,:)), indParName))
end

% Some parameters are not going to be used, they do need to be
% defined, therefore to save time in the for loop, we set them
% here.

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
simdat = [];

acc = 1;
hit = 1; miss = 2; false_alarm = 3; % codes for response types

for indNoise = 1:size(noise,3) % for robustness, it can help to simulate more trials than there are in the real data...
    
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
        if length(trialnoise) > 1; currTrialNoise = trialnoise(bt(indBlock)); else; currTrialNoise = trialnoise; end
        
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
        timeLine = TOW{indBlock}(:,1)';
        timeLine(timeLine ~= 0) = currDrift(abs(timeLine(timeLine ~= 0)));
        
        sensEv = timeLine + currTrialNoise.*noise(indBlock,1:length(TOW{indBlock}), indNoise);
        
        % initalization of parameters.
        DV(1:round(ndt/dt)) = 0;   % initialize DV to a single 0 for the first time points of the block, up until the pre-decision nondecision time
        lastresp = -postRespPause; % we only accumulate if it's a certain time since last response. Initialise to this value so accumulation begins right away at start of block
        
        RespT   = []; % keep track of all responses times
        targT   = []; % keep track of all targets onset times
        
        DVendval = nan;
        
        % now simulate the block by looping through all sample points:
        for indTime = round(ndt/dt)+1:length(sensEv)
            
            if indTime <= lastresp+continAccDur || indTime > lastresp + postRespPause % don't accumulate unless it has been a sufficient time since last response
                
                % this gives us the 'reflecting bound'
                if reflectingBound
                    DV(indTime) = max((1-currLeak)*DV(indTime-1) + sensEv(indTime-round(ndt/dt)) - currCriteria, 0);% There is the main model equation! (see e.g. Ossmy et al 2013)
                    % TODO REFLECTING BOUND FOR
                    % DISCRIMINATION!!!
                else
                    DV(indTime) = (1-currLeak)*DV(indTime-1) + sensEv(indTime-round(ndt/dt)) - currCriteria;   % There is the main model equation! (see e.g. Ossmy et al 2013)
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
            
            % Detect target transitions
            if TOW{indBlock}(indTime) > TOW{indBlock}(indTime-1)
                targT = [targT indTime*dt];
            end
            
            % detect responses
            if indTime > waitAfterTarget
                if indTime-waitAfterTarget <= size(urgency,1)
                    urgencyDV(indTime) = DV(indTime) + urgency(indTime-waitAfterTarget,bt(indBlock)); % total decision signal including urgency
                else  % when there's a miss, hard to know what to do with urgency - most neutral thing might be to keep it steady where it is
                    urgencyDV(indTime) = DV(indTime) + urgency(end, bt(indBlock));
                end
            else
                urgencyDV(indTime) = DV(indTime);
            end
            
            % Now check for bound crossing. Only do this if it is not within postRespPause since last bound crossing
            if indTime > lastresp+postRespPause
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
            nextrespind = find(RespT > targT(indTarget) + minRT & RespT < targT(indTarget) + maxRT, 1); % find the index of the next response which is within the allowable @hit@ window
            if ~isempty(nextrespind)
                simdat  = [simdat; trialMatrix{indBlock}(indTarget,:) hit RespT(nextrespind)-targT(indTarget)]; % 1 = hit
            else % if there WAS no next response, set the response parameters for this trial as 'not a number'
                simdat  = [simdat; trialMatrix{indBlock}(indTarget,:) miss nan];
            end
            
            if getsimDV == 1
                if (targT(indTarget)/dt) + ceil((obj.stim.targetEpoch(end))/dt) < length(DV)
                    DecValue.Target(1:length(targetPlot),acc) = DV(round(targT(indTarget)/dt) + targetPlot);
                    
                    if ~isempty(nextrespind)
                        DecValue.Response(1:length(responsePlot),acc) = DV(round(RespT(nextrespind)/dt) + responsePlot);
                    else
                        DecValue.Response(1:length(responsePlot),acc) = nan;
                    end
                    
                else
                    DecValue.Target(1:length(targetPlot),acc) = nan;
                    if ~isempty(nextrespind)
                        DecValue.Response(1:length(responsePlot),acc) = nan;
                    end
                end
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
                simdat  = [simdat; trialMatrix{indBlock}(indITI,:) false_alarm RespT(nextrespind(m))-endtime];
                
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
